from __future__ import annotations

import collections
import itertools as it
import operator
import warnings
from functools import partial

import numpy as np
import pandas as pd

from dask.base import tokenize
from dask.core import flatten
from dask.dataframe._compat import (
    PANDAS_GE_220,
    PANDAS_GE_300,
    check_groupby_axis_deprecation,
    check_observed_deprecation,
)
from dask.dataframe.core import _convert_to_numeric
from dask.dataframe.methods import concat
from dask.dataframe.utils import (
    get_numeric_only_kwargs,
    is_dataframe_like,
    is_series_like,
)
from dask.typing import no_default
from dask.utils import M, funcname, itemgetter

NUMERIC_ONLY_NOT_IMPLEMENTED = [
    "mean",
    "std",
    "var",
]

GROUP_KEYS_DEFAULT: bool | None = True


def _determine_levels(by):
    """Determine the correct levels argument to groupby."""
    if isinstance(by, (tuple, list)) and len(by) > 1:
        return list(range(len(by)))
    else:
        return 0


def _is_aligned(df, by):
    """Check if ``df`` and ``by`` have aligned indices"""
    if is_series_like(by) or is_dataframe_like(by):
        return df.index.equals(by.index)
    elif isinstance(by, (list, tuple)):
        return all(_is_aligned(df, i) for i in by)
    else:
        return True


def _groupby_raise_unaligned(df, convert_by_to_list=True, **kwargs):
    """Groupby, but raise if df and `by` key are unaligned.

    Pandas supports grouping by a column that doesn't align with the input
    frame/series/index. However, the reindexing does not seem to be
    threadsafe, and can result in incorrect results. Since grouping by an
    unaligned key is generally a bad idea, we just error loudly in dask.

    For more information see pandas GH issue #15244 and Dask GH issue #1876."""
    by = kwargs.get("by", None)
    if by is not None and not _is_aligned(df, by):
        msg = (
            "Grouping by an unaligned column is unsafe and unsupported.\n"
            "This can be caused by filtering only one of the object or\n"
            "grouping key. For example, the following works in pandas,\n"
            "but not in dask:\n"
            "\n"
            "df[df.foo < 0].groupby(df.bar)\n"
            "\n"
            "This can be avoided by either filtering beforehand, or\n"
            "passing in the name of the column instead:\n"
            "\n"
            "df2 = df[df.foo < 0]\n"
            "df2.groupby(df2.bar)\n"
            "# or\n"
            "df[df.foo < 0].groupby('bar')\n"
            "\n"
            "For more information see dask GH issue #1876."
        )
        raise ValueError(msg)
    elif by is not None and len(by) and convert_by_to_list:
        # since we're coming through apply, `by` will be a tuple.
        # Pandas treats tuples as a single key, and lists as multiple keys
        # We want multiple keys
        if isinstance(by, str):
            by = [by]
        kwargs.update(by=list(by))
    with check_observed_deprecation():
        return df.groupby(**kwargs)


def _groupby_slice_apply(
    df,
    grouper,
    key,
    func,
    *args,
    group_keys=GROUP_KEYS_DEFAULT,
    dropna=None,
    observed=None,
    **kwargs,
):
    # No need to use raise if unaligned here - this is only called after
    # shuffling, which makes everything aligned already
    dropna = {"dropna": dropna} if dropna is not None else {}
    observed = {"observed": observed} if observed is not None else {}
    g = df.groupby(grouper, group_keys=group_keys, **observed, **dropna)
    if key:
        g = g[key]
    return g.apply(func, *args, **kwargs)


def _groupby_slice_transform(
    df,
    grouper,
    key,
    func,
    *args,
    group_keys=GROUP_KEYS_DEFAULT,
    dropna=None,
    observed=None,
    **kwargs,
):
    # No need to use raise if unaligned here - this is only called after
    # shuffling, which makes everything aligned already
    dropna = {"dropna": dropna} if dropna is not None else {}
    observed = {"observed": observed} if observed is not None else {}
    g = df.groupby(grouper, group_keys=group_keys, **observed, **dropna)
    if key:
        g = g[key]

    # Cannot call transform on an empty dataframe
    if len(df) == 0:
        return g.apply(func, *args, **kwargs)

    return g.transform(func, *args, **kwargs)


def _groupby_slice_shift(
    df,
    grouper,
    key,
    shuffled,
    group_keys=GROUP_KEYS_DEFAULT,
    dropna=None,
    observed=None,
    **kwargs,
):
    # No need to use raise if unaligned here - this is only called after
    # shuffling, which makes everything aligned already
    dropna = {"dropna": dropna} if dropna is not None else {}
    observed = {"observed": observed} if observed is not None else {}
    if shuffled:
        df = df.sort_index()
    g = df.groupby(grouper, group_keys=group_keys, **observed, **dropna)
    if key:
        g = g[key]
    with check_groupby_axis_deprecation():
        result = g.shift(**kwargs)
    return result


def _groupby_get_group(df, by_key, get_key, columns):
    # SeriesGroupBy may pass df which includes group key
    grouped = _groupby_raise_unaligned(df, by=by_key, convert_by_to_list=False)

    try:
        if is_dataframe_like(df):
            grouped = grouped[columns]
        return grouped.get_group(get_key)
    except KeyError:
        # to create empty DataFrame/Series, which has the same
        # dtype as the original
        if is_dataframe_like(df):
            # may be SeriesGroupBy
            df = df[columns]
        return df.iloc[0:0]


###############################################################
# Aggregation
###############################################################


class Aggregation:
    """User defined groupby-aggregation.

    This class allows users to define their own custom aggregation in terms of
    operations on Pandas dataframes in a map-reduce style. You need to specify
    what operation to do on each chunk of data, how to combine those chunks of
    data together, and then how to finalize the result.

    See :ref:`dataframe.groupby.aggregate` for more.

    Parameters
    ----------
    name : str
        the name of the aggregation. It should be unique, since intermediate
        result will be identified by this name.
    chunk : callable
        a function that will be called with the grouped column of each
        partition, takes a Pandas SeriesGroupBy in input.
        It can either return a single series or a tuple of series.
        The index has to be equal to the groups.
    agg : callable
        a function that will be called to aggregate the results of each chunk.
        Again the argument(s) will be a Pandas SeriesGroupBy. If ``chunk``
        returned a tuple, ``agg`` will be called with all of them as
        individual positional arguments.
    finalize : callable
        an optional finalizer that will be called with the results from the
        aggregation.

    Examples
    --------
    We could implement ``sum`` as follows:

    >>> custom_sum = dd.Aggregation(
    ...     name='custom_sum',
    ...     chunk=lambda s: s.sum(),
    ...     agg=lambda s0: s0.sum()
    ... )  # doctest: +SKIP
    >>> df.groupby('g').agg(custom_sum)  # doctest: +SKIP

    We can implement ``mean`` as follows:

    >>> custom_mean = dd.Aggregation(
    ...     name='custom_mean',
    ...     chunk=lambda s: (s.count(), s.sum()),
    ...     agg=lambda count, sum: (count.sum(), sum.sum()),
    ...     finalize=lambda count, sum: sum / count,
    ... )  # doctest: +SKIP
    >>> df.groupby('g').agg(custom_mean)  # doctest: +SKIP

    Though of course, both of these are built-in and so you don't need to
    implement them yourself.
    """

    def __init__(self, name, chunk, agg, finalize=None):
        self.chunk = chunk
        self.agg = agg
        self.finalize = finalize
        self.__name__ = name


def _groupby_aggregate(
    df, aggfunc=None, levels=None, dropna=None, sort=False, observed=None, **kwargs
):
    dropna = {"dropna": dropna} if dropna is not None else {}
    observed = {"observed": observed} if observed is not None else {}

    with check_observed_deprecation():
        grouped = df.groupby(level=levels, sort=sort, **observed, **dropna)
    return aggfunc(grouped, **kwargs)


def _groupby_aggregate_spec(
    df, spec, levels=None, dropna=None, sort=False, observed=None, **kwargs
):
    """
    A simpler version of _groupby_aggregate that just calls ``aggregate`` using
    the user-provided spec.
    """
    dropna = {"dropna": dropna} if dropna is not None else {}
    observed = {"observed": observed} if observed is not None else {}
    return df.groupby(level=levels, sort=sort, **observed, **dropna).aggregate(
        spec, **kwargs
    )


def _non_agg_chunk(df, *by, key, dropna=None, observed=None, **kwargs):
    """
    A non-aggregation agg function. This simulates the behavior of an initial
    partitionwise aggregation, but doesn't actually aggregate or throw away
    any data.
    """
    if is_series_like(df):
        # Handle a series-like groupby. `by` could be columns that are not the series,
        # but are like-indexed, so we handle that case by temporarily converting to
        # a dataframe, then setting the index.
        result = df.to_frame().set_index(by[0] if len(by) == 1 else list(by))[df.name]
    else:
        # Handle a frame-like groupby.
        result = df.set_index(list(by))
        if isinstance(key, (tuple, list, set, pd.Index)):
            key = list(key)
        result = result[key]

    # If observed is False, we have to check for categorical indices and possibly enrich
    # them with unobserved values. This function is intended as an initial partitionwise
    # aggregation, so you might expect that we could enrich the frame with unobserved
    # values at the end. However, if you have multiple output partitions, that results
    # in duplicated unobserved values in each partition. So we have to do this step
    # at the start before any shuffling occurs so that we can consolidate all of the
    # unobserved values in a single partition.
    if observed is False:
        has_categoricals = False

        # Search for categorical indices and get new index objects that have all the
        # categories in them.
        if isinstance(result.index, pd.CategoricalIndex):
            has_categoricals = True
            full_index = result.index.categories.copy().rename(result.index.name)
        elif isinstance(result.index, pd.MultiIndex) and any(
            isinstance(level, pd.CategoricalIndex) for level in result.index.levels
        ):
            has_categoricals = True
            full_index = pd.MultiIndex.from_product(
                result.index.levels, names=result.index.names
            )
        if has_categoricals:
            # If we found any categoricals, append unobserved values to the end of the
            # frame.
            new_cats = full_index[~full_index.isin(result.index)]
            empty_data = {
                c: pd.Series(index=new_cats, dtype=result[c].dtype)
                for c in result.columns
            }
            empty = pd.DataFrame(empty_data)
            result = pd.concat([result, empty])

    return result


def _apply_chunk(df, *by, dropna=None, observed=None, **kwargs):
    func = kwargs.pop("chunk")
    columns = kwargs.pop("columns")
    dropna = {"dropna": dropna} if dropna is not None else {}
    observed = {"observed": observed} if observed is not None else {}

    g = _groupby_raise_unaligned(df, by=by, **observed, **dropna)
    if is_series_like(df) or columns is None:
        return func(g, **kwargs)
    else:
        if isinstance(columns, (tuple, list, set, pd.Index)):
            columns = list(columns)
        return func(g[columns], **kwargs)


def _var_chunk(df, *by, numeric_only=no_default, observed=False, dropna=True):
    numeric_only_kwargs = get_numeric_only_kwargs(numeric_only)
    if is_series_like(df):
        df = df.to_frame()

    df = df.copy()

    g = _groupby_raise_unaligned(df, by=by, observed=observed, dropna=dropna)
    x = g.sum(**numeric_only_kwargs)

    n = g[x.columns].count().rename(columns=lambda c: (c, "-count"))

    cols = x.columns
    df[cols] = df[cols] ** 2

    g2 = _groupby_raise_unaligned(df, by=by, observed=observed, dropna=dropna)
    x2 = g2.sum(**numeric_only_kwargs).rename(columns=lambda c: (c, "-x2"))

    return concat([x, x2, n], axis=1)


def _var_agg(
    g, levels, ddof, sort=False, numeric_only=no_default, observed=False, dropna=True
):
    numeric_only_kwargs = get_numeric_only_kwargs(numeric_only)
    g = g.groupby(level=levels, sort=sort, observed=observed, dropna=dropna).sum(
        **numeric_only_kwargs
    )
    nc = len(g.columns)
    x = g[g.columns[: nc // 3]]
    # chunks columns are tuples (value, name), so we just keep the value part
    x2 = g[g.columns[nc // 3 : 2 * nc // 3]].rename(columns=lambda c: c[0])
    n = g[g.columns[-nc // 3 :]].rename(columns=lambda c: c[0])

    # TODO: replace with _finalize_var?
    result = x2 - x**2 / n
    div = n - ddof
    div[div < 0] = 0
    result /= div
    result[(n - ddof) == 0] = np.nan
    assert is_dataframe_like(result)
    result[result < 0] = 0  # avoid rounding errors that take us to zero
    return result


def _cov_finalizer(df, cols, std=False):
    num_elements = len(list(it.product(cols, repeat=2)))
    num_cols = len(cols)
    vals = list(range(num_elements))
    col_idx_mapping = dict(zip(cols, range(num_cols)))
    for i, j in it.combinations_with_replacement(df[cols].columns, 2):
        x = col_idx_mapping[i]
        y = col_idx_mapping[j]
        idx = x + num_cols * y
        mul_col = f"{i}{j}"
        ni = df["%s-count" % i]
        nj = df["%s-count" % j]

        n = np.sqrt(ni * nj)
        div = n - 1
        div[div < 0] = 0
        val = (df[mul_col] - df[i] * df[j] / n).values[0] / div.values[0]
        if std:
            ii = f"{i}{i}"
            jj = f"{j}{j}"
            std_val_i = (df[ii] - (df[i] ** 2) / ni).values[0] / div.values[0]
            std_val_j = (df[jj] - (df[j] ** 2) / nj).values[0] / div.values[0]
            sqrt_val = np.sqrt(std_val_i * std_val_j)
            if sqrt_val == 0:
                val = np.nan
            else:
                val = val / sqrt_val

        vals[idx] = val
        if i != j:
            idx = num_cols * x + y
            vals[idx] = val

    level_1 = cols
    index = pd.MultiIndex.from_product([level_1, level_1])
    return pd.Series(vals, index=index)


def _mul_cols(df, cols):
    """Internal function to be used with apply to multiply
    each column in a dataframe by every other column

    a b c -> a*a, a*b, b*b, b*c, c*c
    """
    _df = df.__class__()
    for i, j in it.combinations_with_replacement(cols, 2):
        col = f"{i}{j}"
        _df[col] = df[i] * df[j]

    # Fix index in a groupby().apply() context
    # https://github.com/dask/dask/issues/8137
    # https://github.com/pandas-dev/pandas/issues/43568
    # Make sure index dtype is int (even if _df is empty)
    # https://github.com/dask/dask/pull/9701
    _df.index = np.zeros(len(_df), dtype=int)
    return _df


def _cov_chunk(df, *by, numeric_only=no_default):
    """Covariance Chunk Logic

    Parameters
    ----------
    df : Pandas.DataFrame
    std : bool, optional
        When std=True we are calculating with Correlation

    Returns
    -------
    tuple
        Processed X, Multiplied Cols,
    """
    numeric_only_kwargs = get_numeric_only_kwargs(numeric_only)
    if is_series_like(df):
        df = df.to_frame()
    df = df.copy()
    if numeric_only is False:
        dt_df = df.select_dtypes(include=["datetime", "timedelta"])
        for col in dt_df.columns:
            df[col] = _convert_to_numeric(dt_df[col], True)

    # mapping columns to str(numerical) values allows us to easily handle
    # arbitrary column names (numbers, string, empty strings)
    col_mapping = collections.OrderedDict()
    for i, c in enumerate(df.columns):
        col_mapping[c] = str(i)
    df = df.rename(columns=col_mapping)
    cols = df._get_numeric_data().columns

    # when grouping by external series don't exclude columns
    is_mask = any(is_series_like(s) for s in by)
    if not is_mask:
        by = [col_mapping[k] for k in by]
        cols = cols.difference(pd.Index(by))

    g = _groupby_raise_unaligned(df, by=by)
    x = g.sum(**numeric_only_kwargs)

    include_groups = {"include_groups": False} if PANDAS_GE_220 else {}
    mul = g.apply(_mul_cols, cols=cols, **include_groups).reset_index(
        level=-1, drop=True
    )

    n = g[x.columns].count().rename(columns=lambda c: f"{c}-count")
    return (x, mul, n, col_mapping)


def _cov_agg(_t, levels, ddof, std=False, sort=False):
    sums = []
    muls = []
    counts = []

    # sometime we get a series back from concat combiner
    t = list(_t)

    cols = t[0][0].columns
    for x, mul, n, col_mapping in t:
        sums.append(x)
        muls.append(mul)
        counts.append(n)
        col_mapping = col_mapping

    total_sums = concat(sums).groupby(level=levels, sort=sort).sum()
    total_muls = concat(muls).groupby(level=levels, sort=sort).sum()
    total_counts = concat(counts).groupby(level=levels).sum()
    result = (
        concat([total_sums, total_muls, total_counts], axis=1)
        .groupby(level=levels)
        .apply(_cov_finalizer, cols=cols, std=std)
    )

    inv_col_mapping = {v: k for k, v in col_mapping.items()}
    idx_vals = result.index.names
    idx_mapping = list()

    # when index is None we probably have selected a particular column
    # df.groupby('a')[['b']].cov()
    if len(idx_vals) == 1 and all(n is None for n in idx_vals):
        idx_vals = list(inv_col_mapping.keys() - set(total_sums.columns))

    for val in idx_vals:
        idx_name = inv_col_mapping.get(val, val)
        idx_mapping.append(idx_name)

        if len(result.columns.levels[0]) < len(col_mapping):
            # removing index from col_mapping (produces incorrect multiindexes)
            try:
                col_mapping.pop(idx_name)
            except KeyError:
                # when slicing the col_map will not have the index
                pass

    keys = list(col_mapping.keys())
    for level in range(len(result.columns.levels)):
        result.columns = result.columns.set_levels(keys, level=level)

    result.index.set_names(idx_mapping, inplace=True)

    # stacking can lead to a sorted index
    if PANDAS_GE_300:
        s_result = result.stack()
    else:
        s_result = result.stack(dropna=False)
    assert is_dataframe_like(s_result)
    return s_result


###############################################################
# nunique
###############################################################


def _nunique_df_chunk(df, *by, **kwargs):
    name = kwargs.pop("name")
    try:
        # This is a lot faster but kind of a pain to implement when by
        # has a boolean series in it.
        return df.drop_duplicates(subset=list(by) + [name]).set_index(list(by))
    except Exception:
        pass

    g = _groupby_raise_unaligned(df, by=by, group_keys=True)
    if len(df) > 0:
        grouped = g[name].unique().explode().to_frame()
    else:
        # Manually create empty version, since groupby-apply for empty frame
        # results in df with no columns
        grouped = g[[name]].nunique()
        grouped = grouped.astype(df.dtypes[grouped.columns].to_dict())

    return grouped


def _nunique_df_combine(df, levels, sort=False):
    result = (
        df.groupby(level=levels, sort=sort, observed=True)[df.columns[0]]
        .unique()
        .explode()
        .to_frame()
    )
    return result


###############################################################
# Aggregate support
#
# Aggregate is implemented as:
#
# 1. group-by-aggregate all partitions into intermediate values
# 2. collect all partitions into a single partition
# 3. group-by-aggregate the result into intermediate values
# 4. transform all intermediate values into the result
#
# In Step 1 and 3 the dataframe is grouped on the same columns.
#
###############################################################
def _make_agg_id(func, column):
    return f"{func!s}-{column!s}-{tokenize(func, column)}"


def _normalize_spec(spec, non_group_columns):
    """
    Return a list of ``(result_column, func, input_column)`` tuples.

    Spec can be

    - a function
    - a list of functions
    - a dictionary that maps input-columns to functions
    - a dictionary that maps input-columns to a lists of functions
    - a dictionary that maps input-columns to a dictionaries that map
      output-columns to functions.

    The non-group columns are a list of all column names that are not used in
    the groupby operation.

    Usually, the result columns are mutli-level names, returned as tuples.
    If only a single function is supplied or dictionary mapping columns
    to single functions, simple names are returned as strings (see the first
    two examples below).

    Examples
    --------
    >>> _normalize_spec('mean', ['a', 'b', 'c'])
    [('a', 'mean', 'a'), ('b', 'mean', 'b'), ('c', 'mean', 'c')]

    >>> spec = collections.OrderedDict([('a', 'mean'), ('b', 'count')])
    >>> _normalize_spec(spec, ['a', 'b', 'c'])
    [('a', 'mean', 'a'), ('b', 'count', 'b')]

    >>> _normalize_spec(['var', 'mean'], ['a', 'b', 'c'])
    ... # doctest: +NORMALIZE_WHITESPACE
    [(('a', 'var'), 'var', 'a'), (('a', 'mean'), 'mean', 'a'), \
     (('b', 'var'), 'var', 'b'), (('b', 'mean'), 'mean', 'b'), \
     (('c', 'var'), 'var', 'c'), (('c', 'mean'), 'mean', 'c')]

    >>> spec = collections.OrderedDict([('a', 'mean'), ('b', ['sum', 'count'])])
    >>> _normalize_spec(spec, ['a', 'b', 'c'])
    ... # doctest: +NORMALIZE_WHITESPACE
    [(('a', 'mean'), 'mean', 'a'), (('b', 'sum'), 'sum', 'b'), \
      (('b', 'count'), 'count', 'b')]

    >>> spec = collections.OrderedDict()
    >>> spec['a'] = ['mean', 'size']
    >>> spec['b'] = collections.OrderedDict([('e', 'count'), ('f', 'var')])
    >>> _normalize_spec(spec, ['a', 'b', 'c'])
    ... # doctest: +NORMALIZE_WHITESPACE
    [(('a', 'mean'), 'mean', 'a'), (('a', 'size'), 'size', 'a'), \
     (('b', 'e'), 'count', 'b'), (('b', 'f'), 'var', 'b')]
    """
    if not isinstance(spec, dict):
        spec = collections.OrderedDict(zip(non_group_columns, it.repeat(spec)))

    res = []

    if isinstance(spec, dict):
        for input_column, subspec in spec.items():
            if isinstance(subspec, dict):
                res.extend(
                    ((input_column, result_column), func, input_column)
                    for result_column, func in subspec.items()
                )

            else:
                if not isinstance(subspec, list):
                    subspec = [subspec]

                res.extend(
                    ((input_column, funcname(func)), func, input_column)
                    for func in subspec
                )

    else:
        raise ValueError(f"unsupported agg spec of type {type(spec)}")

    compounds = (list, tuple, dict)
    use_flat_columns = not any(
        isinstance(subspec, compounds) for subspec in spec.values()
    )

    if use_flat_columns:
        res = [(input_col, func, input_col) for (_, func, input_col) in res]

    return res


def _build_agg_args(spec):
    """
    Create transformation functions for a normalized aggregate spec.

    Parameters
    ----------
    spec: a list of (result-column, aggregation-function, input-column) triples.
        To work with all argument forms understood by pandas use
        ``_normalize_spec`` to normalize the argument before passing it on to
        ``_build_agg_args``.

    Returns
    -------
    chunk_funcs: a list of (intermediate-column, function, keyword) triples
        that are applied on grouped chunks of the initial dataframe.

    agg_funcs: a list of (intermediate-column, functions, keyword) triples that
        are applied on the grouped concatenation of the preprocessed chunks.

    finalizers: a list of (result-column, function, keyword) triples that are
        applied after the ``agg_funcs``. They are used to create final results
        from intermediate representations.
    """
    known_np_funcs = {
        np.min: "min",
        np.max: "max",
        np.median: "median",
        np.std: "std",
        np.var: "var",
    }

    # check that there are no name conflicts for a single input column
    by_name = {}
    for _, func, input_column in spec:
        key = funcname(known_np_funcs.get(func, func)), input_column
        by_name.setdefault(key, []).append((func, input_column))

    for funcs in by_name.values():
        if len(funcs) != 1:
            raise ValueError(f"conflicting aggregation functions: {funcs}")

    chunks = {}
    aggs = {}
    finalizers = []

    # a partial may contain some arguments, pass them down
    # https://github.com/dask/dask/issues/9615
    for result_column, func, input_column in spec:
        func_args = ()
        func_kwargs = {}
        if isinstance(func, partial):
            func_args, func_kwargs = func.args, func.keywords

        if not isinstance(func, Aggregation):
            func = funcname(known_np_funcs.get(func, func))

        impls = _build_agg_args_single(
            result_column, func, func_args, func_kwargs, input_column
        )

        # overwrite existing result-columns, generate intermediates only once
        for spec in impls["chunk_funcs"]:
            chunks[spec[0]] = spec
        for spec in impls["aggregate_funcs"]:
            aggs[spec[0]] = spec

        finalizers.append(impls["finalizer"])

    chunks = sorted(chunks.values())
    aggs = sorted(aggs.values())

    return chunks, aggs, finalizers


def _build_agg_args_single(result_column, func, func_args, func_kwargs, input_column):
    simple_impl = {
        "sum": (M.sum, M.sum),
        "min": (M.min, M.min),
        "max": (M.max, M.max),
        "count": (M.count, M.sum),
        "size": (M.size, M.sum),
        "first": (M.first, M.first),
        "last": (M.last, M.last),
        "prod": (M.prod, M.prod),
        "median": (
            None,
            M.median,
        ),  # No chunk func for median, we can only take it when aggregating
    }

    if func in simple_impl.keys():
        return _build_agg_args_simple(
            result_column, func, input_column, simple_impl[func]
        )

    elif func == "var":
        return _build_agg_args_var(
            result_column, func, func_args, func_kwargs, input_column
        )

    elif func == "std":
        return _build_agg_args_std(
            result_column, func, func_args, func_kwargs, input_column
        )

    elif func == "mean":
        return _build_agg_args_mean(result_column, func, input_column)

    elif func == "list":
        return _build_agg_args_list(result_column, func, input_column)

    elif isinstance(func, Aggregation):
        return _build_agg_args_custom(result_column, func, input_column)

    else:
        raise ValueError(f"unknown aggregate {func}")


def _build_agg_args_simple(result_column, func, input_column, impl_pair):
    intermediate = _make_agg_id(func, input_column)
    chunk_impl, agg_impl = impl_pair

    return dict(
        chunk_funcs=[
            (
                intermediate,
                _apply_func_to_column,
                dict(column=input_column, func=chunk_impl),
            )
        ],
        aggregate_funcs=[
            (
                intermediate,
                _apply_func_to_column,
                dict(column=intermediate, func=agg_impl),
            )
        ],
        finalizer=(result_column, itemgetter(intermediate), dict()),
    )


def _build_agg_args_var(result_column, func, func_args, func_kwargs, input_column):
    int_sum = _make_agg_id("sum", input_column)
    int_sum2 = _make_agg_id("sum2", input_column)
    int_count = _make_agg_id("count", input_column)

    # we don't expect positional args here
    if func_args:
        raise TypeError(
            f"aggregate function '{func}' doesn't support positional arguments, but got {func_args}"
        )

    # and we only expect ddof=N in kwargs
    expected_kwargs = {"ddof"}
    unexpected_kwargs = func_kwargs.keys() - expected_kwargs
    if unexpected_kwargs:
        raise TypeError(
            f"aggregate function '{func}' supports {expected_kwargs} keyword arguments, but got {unexpected_kwargs}"
        )

    return dict(
        chunk_funcs=[
            (int_sum, _apply_func_to_column, dict(column=input_column, func=M.sum)),
            (int_count, _apply_func_to_column, dict(column=input_column, func=M.count)),
            (int_sum2, _compute_sum_of_squares, dict(column=input_column)),
        ],
        aggregate_funcs=[
            (col, _apply_func_to_column, dict(column=col, func=M.sum))
            for col in (int_sum, int_count, int_sum2)
        ],
        finalizer=(
            result_column,
            _finalize_var,
            dict(
                sum_column=int_sum,
                count_column=int_count,
                sum2_column=int_sum2,
                **func_kwargs,
            ),
        ),
    )


def _build_agg_args_std(result_column, func, func_args, func_kwargs, input_column):
    impls = _build_agg_args_var(
        result_column, func, func_args, func_kwargs, input_column
    )

    result_column, _, kwargs = impls["finalizer"]
    impls["finalizer"] = (result_column, _finalize_std, kwargs)

    return impls


def _build_agg_args_mean(result_column, func, input_column):
    int_sum = _make_agg_id("sum", input_column)
    int_count = _make_agg_id("count", input_column)

    return dict(
        chunk_funcs=[
            (int_sum, _apply_func_to_column, dict(column=input_column, func=M.sum)),
            (int_count, _apply_func_to_column, dict(column=input_column, func=M.count)),
        ],
        aggregate_funcs=[
            (col, _apply_func_to_column, dict(column=col, func=M.sum))
            for col in (int_sum, int_count)
        ],
        finalizer=(
            result_column,
            _finalize_mean,
            dict(sum_column=int_sum, count_column=int_count),
        ),
    )


def _build_agg_args_list(result_column, func, input_column):
    intermediate = _make_agg_id("list", input_column)

    return dict(
        chunk_funcs=[
            (
                intermediate,
                _apply_func_to_column,
                dict(column=input_column, func=lambda s: s.apply(list)),
            )
        ],
        aggregate_funcs=[
            (
                intermediate,
                _apply_func_to_column,
                dict(
                    column=intermediate,
                    func=lambda s0: s0.apply(
                        lambda chunks: list(it.chain.from_iterable(chunks))
                    ),
                ),
            )
        ],
        finalizer=(result_column, itemgetter(intermediate), dict()),
    )


def _build_agg_args_custom(result_column, func, input_column):
    col = _make_agg_id(funcname(func), input_column)

    if func.finalize is None:
        finalizer = (result_column, operator.itemgetter(col), dict())

    else:
        finalizer = (
            result_column,
            _apply_func_to_columns,
            dict(func=func.finalize, prefix=col),
        )

    return dict(
        chunk_funcs=[
            (col, _apply_func_to_column, dict(func=func.chunk, column=input_column))
        ],
        aggregate_funcs=[
            (col, _apply_func_to_columns, dict(func=func.agg, prefix=col))
        ],
        finalizer=finalizer,
    )


def _groupby_apply_funcs(df, *by, **kwargs):
    """
    Group a dataframe and apply multiple aggregation functions.

    Parameters
    ----------
    df: pandas.DataFrame
        The dataframe to work on.
    by: list of groupers
        If given, they are added to the keyword arguments as the ``by``
        argument.
    funcs: list of result-colum, function, keywordargument triples
        The list of functions that are applied on the grouped data frame.
        Has to be passed as a keyword argument.
    kwargs:
        All keyword arguments, but ``funcs``, are passed verbatim to the groupby
        operation of the dataframe

    Returns
    -------
    aggregated:
        the aggregated dataframe.
    """
    if len(by):
        # since we're coming through apply, `by` will be a tuple.
        # Pandas treats tuples as a single key, and lists as multiple keys
        # We want multiple keys
        kwargs.update(by=list(by))

    funcs = kwargs.pop("funcs")
    grouped = _groupby_raise_unaligned(df, **kwargs)

    result = collections.OrderedDict()
    for result_column, func, func_kwargs in funcs:
        r = func(grouped, **func_kwargs)

        if isinstance(r, tuple):
            for idx, s in enumerate(r):
                result[f"{result_column}-{idx}"] = s

        else:
            result[result_column] = r

    if is_dataframe_like(df):
        return df.__class__(result)
    else:
        # Get the DataFrame type of this Series object
        return df.head(0).to_frame().__class__(result)


def _compute_sum_of_squares(grouped, column):
    # Note: CuDF cannot use `groupby.apply`.
    # Need to unpack groupby to compute sum of squares
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            "DataFrameGroupBy.grouper is deprecated and will be removed in a future version of pandas.",
            FutureWarning,
        )
        # TODO: Avoid usage of grouper
        if hasattr(grouped, "grouper"):
            keys = grouped.grouper
        elif hasattr(grouped, "_grouper"):
            keys = grouped._grouper
        else:
            # Handle CuDF groupby object (different from pandas)
            keys = grouped.grouping.keys
    df = grouped.obj[column].pow(2) if column else grouped.obj.pow(2)
    return df.groupby(keys).sum()


def _agg_finalize(
    df,
    aggregate_funcs,
    finalize_funcs,
    level,
    sort=False,
    arg=None,
    columns=None,
    is_series=False,
    **kwargs,
):
    # finish the final aggregation level
    df = _groupby_apply_funcs(
        df, funcs=aggregate_funcs, level=level, sort=sort, **kwargs
    )

    # and finalize the result
    result = collections.OrderedDict()
    for result_column, func, finalize_kwargs in finalize_funcs:
        result[result_column] = func(df, **finalize_kwargs)

    result = df.__class__(result)
    if columns is not None:
        try:
            result = result[columns]
        except KeyError:
            pass
    if (
        is_series
        and arg is not None
        and not isinstance(arg, (list, dict))
        and result.ndim == 2
    ):
        result = result[result.columns[0]]
    return result


def _apply_func_to_column(df_like, column, func):
    if column is None:
        return func(df_like)

    return func(df_like[column])


def _apply_func_to_columns(df_like, prefix, func):
    if is_dataframe_like(df_like):
        columns = df_like.columns
    else:
        # handle GroupBy objects
        columns = df_like.obj.columns

    columns = sorted(col for col in columns if col.startswith(prefix))

    columns = [df_like[col] for col in columns]
    return func(*columns)


def _finalize_mean(df, sum_column, count_column):
    result = df[sum_column] / df[count_column]
    return _adjust_for_arrow_na(result, df[count_column])


def _adjust_for_arrow_na(result, df, check_for_isna=False):
    if isinstance(result.dtype, pd.ArrowDtype):
        # Our mean computation results in np.nan here but pandas doesn't
        if check_for_isna:
            result[df.isna()] = pd.NA
        else:
            result[df == 0] = pd.NA
    return result


def _finalize_var(
    df, count_column, sum_column, sum2_column, adjust_arrow=True, **kwargs
):
    # arguments are being checked when building the finalizer. As of this moment,
    # we're only using ddof, and raising an error on other keyword args.
    ddof = kwargs.get("ddof", 1)
    n = df[count_column]
    x = df[sum_column]
    x2 = df[sum2_column]

    result = x2 - x**2 / n
    div = n - ddof
    div[div < 0] = 0
    result /= div
    result[(n - ddof) == 0] = np.nan
    if adjust_arrow:
        return _adjust_for_arrow_na(result, div)
    else:
        return result


def _finalize_std(df, count_column, sum_column, sum2_column, **kwargs):
    result = _finalize_var(
        df, count_column, sum_column, sum2_column, adjust_arrow=False, **kwargs
    )
    res = np.sqrt(result)
    if res.dtype != result.dtype:
        res = res.astype(result.dtype)
    return _adjust_for_arrow_na(res, result, check_for_isna=True)


def _cum_agg_aligned(part, cum_last, index, columns, func, initial):
    align = cum_last.reindex(part.set_index(index).index, fill_value=initial)
    align.index = part.index
    return func(part[columns], align)


def _cum_agg_filled(a, b, func, initial):
    union = a.index.union(b.index)
    return func(
        a.reindex(union, fill_value=initial),
        b.reindex(union, fill_value=initial),
        fill_value=initial,
    )


def _cumcount_aggregate(a, b, fill_value=None):
    return a.add(b, fill_value=fill_value) + 1


def _aggregate_docstring(based_on=None):
    # Insert common groupby-aggregation docstring.
    # Use `based_on` parameter to add note about the
    # Pandas method the Dask version is based on.
    based_on_str = "\n" if based_on is None else f"\nBased on {based_on}\n"

    def wrapper(func):
        func.__doc__ = f"""Aggregate using one or more specified operations
        {based_on_str}
        Parameters
        ----------
        arg : callable, str, list or dict, optional
            Aggregation spec. Accepted combinations are:

            - callable function
            - string function name
            - list of functions and/or function names, e.g. ``[np.sum, 'mean']``
            - dict of column names -> function, function name or list of such.
            - None only if named aggregation syntax is used
        split_every : int, optional
            Number of intermediate partitions that may be aggregated at once.
            This defaults to 8. If your intermediate partitions are likely to
            be small (either due to a small number of groups or a small initial
            partition size), consider increasing this number for better performance.
        split_out : int, optional
            Number of output partitions. Default is 1.
        shuffle : bool or str, optional
            Whether a shuffle-based algorithm should be used. A specific
            algorithm name may also be specified (e.g. ``"tasks"`` or ``"p2p"``).
            The shuffle-based algorithm is likely to be more efficient than
            ``shuffle=False`` when ``split_out>1`` and the number of unique
            groups is large (high cardinality). Default is ``False`` when
            ``split_out = 1``. When ``split_out > 1``, it chooses the algorithm
            set by the ``shuffle`` option in the dask config system, or ``"tasks"``
            if nothing is set.
        kwargs: tuple or pd.NamedAgg, optional
            Used for named aggregations where the keywords are the output column
            names and the values are tuples where the first element is the input
            column name and the second element is the aggregation function.
            ``pandas.NamedAgg`` can also be used as the value. To use the named
            aggregation syntax, arg must be set to None.
        """
        return func

    return wrapper


def _unique_aggregate(series_gb, name=None):
    data = {k: v.explode().unique() for k, v in series_gb}
    ret = type(series_gb.obj)(data, name=name)
    ret.index.names = series_gb.obj.index.names
    ret.index = ret.index.astype(series_gb.obj.index.dtype, copy=False)
    return ret


def _value_counts(x, **kwargs):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", "`groups` by one element list returns", FutureWarning
        )
        if not x.groups or all(
            pd.isna(key) for key in flatten(x.groups.keys(), container=tuple)
        ):
            return pd.Series(dtype=int)
        else:
            return x.value_counts(**kwargs)


def _value_counts_aggregate(series_gb):
    data = {k: v.groupby(level=-1).sum() for k, v in series_gb}
    if not data:
        data = [pd.Series(index=series_gb.obj.index[:0], dtype="float64")]
    res = pd.concat(data, names=series_gb.obj.index.names)
    typed_levels = {
        i: res.index.levels[i].astype(series_gb.obj.index.levels[i].dtype)
        for i in range(len(res.index.levels))
    }
    res.index = res.index.set_levels(
        typed_levels.values(), level=typed_levels.keys(), verify_integrity=False
    )
    return res


def _tail_chunk(series_gb, **kwargs):
    keys, groups = zip(*series_gb) if len(series_gb) else ((True,), (series_gb,))
    return pd.concat([group.tail(**kwargs) for group in groups], keys=keys)


def _tail_aggregate(series_gb, **kwargs):
    levels = kwargs.pop("index_levels")
    return series_gb.tail(**kwargs).droplevel(list(range(levels)))


def _head_chunk(series_gb, **kwargs):
    keys, groups = zip(*series_gb) if len(series_gb) else ((True,), (series_gb,))
    return pd.concat([group.head(**kwargs) for group in groups], keys=keys)


def _head_aggregate(series_gb, **kwargs):
    levels = kwargs.pop("index_levels")
    return series_gb.head(**kwargs).droplevel(list(range(levels)))
