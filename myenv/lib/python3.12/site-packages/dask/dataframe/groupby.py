from __future__ import annotations

import collections
import itertools as it
import operator
import uuid
import warnings
from functools import partial, wraps
from numbers import Integral

import numpy as np
import pandas as pd

from dask.base import is_dask_collection, tokenize
from dask.core import flatten
from dask.dataframe._compat import (
    PANDAS_GE_140,
    PANDAS_GE_150,
    PANDAS_GE_200,
    PANDAS_GE_210,
    PANDAS_GE_220,
    PANDAS_GE_300,
    check_groupby_axis_deprecation,
    check_numeric_only_deprecation,
    check_observed_deprecation,
)
from dask.dataframe.core import (
    GROUP_KEYS_DEFAULT,
    DataFrame,
    Series,
    _convert_to_numeric,
    _determine_split_out_shuffle,
    _extract_meta,
    _Frame,
    aca,
    map_partitions,
    new_dd_object,
    split_out_on_index,
)
from dask.dataframe.dispatch import grouper_dispatch
from dask.dataframe.methods import concat, drop_columns
from dask.dataframe.utils import (
    get_numeric_only_kwargs,
    insert_meta_param_description,
    is_dataframe_like,
    is_index_like,
    is_series_like,
    make_meta,
    raise_on_meta_error,
)
from dask.highlevelgraph import HighLevelGraph
from dask.typing import no_default
from dask.utils import (
    M,
    _deprecated,
    _deprecated_kwarg,
    derived_from,
    funcname,
    itemgetter,
)

if PANDAS_GE_140:
    from pandas.core.apply import reconstruct_func, validate_func_kwargs

# #############################################
#
# GroupBy implementation notes
#
# Dask groupby supports reductions, i.e., mean, sum and alike, and apply. The
# former do not shuffle the data and are efficiently implemented as tree
# reductions. The latter is implemented by shuffling the underlying partitions
# such that all items of a group can be found in the same partition.
#
# The argument to ``.groupby`` (``by``), can be a ``str``, ``dd.DataFrame``,
# ``dd.Series``, or a list thereof. In operations on the grouped object, the
# divisions of the the grouped object and the items of ``by`` have to align.
# Currently, there is no support to shuffle the ``by`` values as part of the
# groupby operation. Therefore, the alignment has to be guaranteed by the
# caller.
#
# To operate on matching partitions, most groupby operations exploit the
# corresponding support in ``apply_concat_apply``. Specifically, this function
# operates on matching partitions of frame-like objects passed as varargs.
#
# After the initial chunk step, ``by``` is implicitly passed along to
# subsequent operations as the index of the partitions. Groupby operations on
# the individual partitions can then access ``by`` via the ``levels``
# parameter of the ``groupby`` function. The correct argument is determined by
# the ``_determine_levels`` function.
#
# To minimize overhead, any ``by`` that is a series contained within the
# dataframe is passed as a column key. This transformation is implemented as
# ``_normalize_by``.
#
# #############################################

NUMERIC_ONLY_NOT_IMPLEMENTED = [
    "mean",
    "std",
    "var",
]


def _determine_levels(by):
    """Determine the correct levels argument to groupby."""
    if isinstance(by, (tuple, list)) and len(by) > 1:
        return list(range(len(by)))
    else:
        return 0


def _normalize_by(df, by):
    """Replace series with column names wherever possible."""
    if not isinstance(df, DataFrame):
        return by

    elif isinstance(by, list):
        return [_normalize_by(df, col) for col in by]

    elif is_series_like(by) and by.name in df.columns and by._name == df[by.name]._name:
        return by.name

    elif (
        isinstance(by, DataFrame)
        and set(by.columns).issubset(df.columns)
        and by._name == df[by.columns]._name
    ):
        return list(by.columns)

    else:
        return by


def _maybe_slice(grouped, columns):
    """
    Slice columns if grouped is pd.DataFrameGroupBy
    """
    # FIXME: update with better groupby object detection (i.e.: ngroups, get_group)
    if "groupby" in type(grouped).__name__.lower():
        if columns is not None:
            if isinstance(columns, (tuple, list, set, pd.Index)):
                columns = list(columns)
            return grouped[columns]
    return grouped


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


def numeric_only_deprecate_default(func):
    """Decorator for methods that should warn when numeric_only is default"""

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if isinstance(self, DataFrameGroupBy):
            numeric_only = kwargs.get("numeric_only", no_default)
            # Prior to `pandas=1.5`, `numeric_only` support wasn't uniformly supported
            # in pandas. We don't support `numeric_only=False` in this case.
            if not PANDAS_GE_150 and numeric_only is False:
                raise NotImplementedError(
                    "'numeric_only=False' is not implemented in Dask."
                )
            if PANDAS_GE_150 and not PANDAS_GE_200 and not self._all_numeric():
                if numeric_only is no_default:
                    warnings.warn(
                        "The default value of numeric_only will be changed to False in "
                        "the future when using dask with pandas 2.0",
                        FutureWarning,
                    )
                elif numeric_only is False and funcname(func) in ("sum", "prod"):
                    warnings.warn(
                        "Dropping invalid columns is deprecated. In a future version, a TypeError will be raised. "
                        f"Before calling .{funcname(func)}, select only columns which should be valid for the function",
                        FutureWarning,
                    )

        return func(self, *args, **kwargs)

    return wrapper


def numeric_only_not_implemented(func):
    """Decorator for methods that can't handle numeric_only=False"""

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if isinstance(self, DataFrameGroupBy):
            maybe_raise = not (
                func.__name__ == "agg"
                and len(args) > 0
                and args[0] not in NUMERIC_ONLY_NOT_IMPLEMENTED
            )
            if maybe_raise:
                numeric_only = kwargs.get("numeric_only", no_default)
                # Prior to `pandas=1.5`, `numeric_only` support wasn't uniformly supported
                # in pandas. We don't support `numeric_only=False` in this case.
                if not PANDAS_GE_150 and numeric_only is False:
                    raise NotImplementedError(
                        "'numeric_only=False' is not implemented in Dask."
                    )
                if not self._all_numeric():
                    if numeric_only is False or (
                        PANDAS_GE_200 and numeric_only is no_default
                    ):
                        raise NotImplementedError(
                            "'numeric_only=False' is not implemented in Dask."
                        )
                    if (
                        PANDAS_GE_150
                        and not PANDAS_GE_200
                        and numeric_only is no_default
                    ):
                        warnings.warn(
                            "The default value of numeric_only will be changed to False "
                            "in the future when using dask with pandas 2.0",
                            FutureWarning,
                        )
        return func(self, *args, **kwargs)

    return wrapper


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
        partition. It can either return a single series or a tuple of series.
        The index has to be equal to the groups.
    agg : callable
        a function that will be called to aggregate the results of each chunk.
        Again the argument(s) will be grouped series. If ``chunk`` returned a
        tuple, ``agg`` will be called with all of them as individual positional
        arguments.
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

    # we emit a warning earlier in stack about default numeric_only being deprecated,
    # so there's no need to propagate the warning that pandas emits as well
    with check_numeric_only_deprecation():
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
    with check_numeric_only_deprecation():
        x = g.sum(**numeric_only_kwargs)

    n = g[x.columns].count().rename(columns=lambda c: (c, "-count"))

    cols = x.columns
    df[cols] = df[cols] ** 2

    g2 = _groupby_raise_unaligned(df, by=by, observed=observed, dropna=dropna)
    with check_numeric_only_deprecation():
        x2 = g2.sum(**numeric_only_kwargs).rename(columns=lambda c: (c, "-x2"))

    return concat([x, x2, n], axis=1)


def _var_combine(g, levels, sort=False):
    return g.groupby(level=levels, sort=sort).sum()


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


def _cov_combine(g, levels):
    return g


def _cov_finalizer(df, cols, std=False):
    vals = []
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
    group_keys = {}
    if PANDAS_GE_150:
        group_keys["group_keys"] = True

    g = _groupby_raise_unaligned(df, by=by, **group_keys)
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


def _nunique_df_aggregate(df, levels, name, sort=False):
    return df.groupby(level=levels, sort=sort, observed=True)[name].nunique()


def _nunique_series_chunk(df, *by, **_ignored_):
    # convert series to data frame, then hand over to dataframe code path
    assert is_series_like(df)

    df = df.to_frame()
    kwargs = dict(name=df.columns[0], levels=_determine_levels(by))
    return _nunique_df_chunk(df, *by, **kwargs)


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
    return df[sum_column] / df[count_column]


def _finalize_var(df, count_column, sum_column, sum2_column, **kwargs):
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

    return result


def _finalize_std(df, count_column, sum_column, sum2_column, **kwargs):
    result = _finalize_var(df, count_column, sum_column, sum2_column, **kwargs)
    return np.sqrt(result)


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


def _drop_apply(group, *, by, what, **kwargs):
    # apply keeps the grouped-by columns, so drop them to stay consistent with pandas
    # groupby-fillna in pandas<2.2
    return getattr(group.drop(columns=by), what)(**kwargs)


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


class _GroupBy:
    """Superclass for DataFrameGroupBy and SeriesGroupBy

    Parameters
    ----------

    obj: DataFrame or Series
        DataFrame or Series to be grouped
    by: str, list or Series
        The key for grouping
    slice: str, list
        The slice keys applied to GroupBy result
    group_keys: bool | None
        Passed to pandas.DataFrame.groupby()
    dropna: bool
        Whether to drop null values from groupby index
    sort: bool
        Passed along to aggregation methods. If allowed,
        the output aggregation will have sorted keys.
    observed: bool, default False
        This only applies if any of the groupers are Categoricals.
        If True: only show observed values for categorical groupers.
        If False: show all values for categorical groupers.
    """

    def __init__(
        self,
        df,
        by=None,
        slice=None,
        group_keys=GROUP_KEYS_DEFAULT,
        dropna=None,
        sort=True,
        observed=None,
    ):
        by_ = by if isinstance(by, (tuple, list)) else [by]
        if any(isinstance(key, pd.Grouper) for key in by_):
            raise NotImplementedError("pd.Grouper is currently not supported by Dask.")

        # slicing key applied to _GroupBy instance
        self._slice = slice

        # Check if we can project columns
        projection = None
        if (
            np.isscalar(self._slice)
            or isinstance(self._slice, (str, list, tuple))
            or (
                (is_index_like(self._slice) or is_series_like(self._slice))
                and not is_dask_collection(self._slice)
            )
        ):
            projection = set(by_).union(
                {self._slice}
                if (np.isscalar(self._slice) or isinstance(self._slice, str))
                else self._slice
            )
            projection = [c for c in df.columns if c in projection]

        assert isinstance(df, (DataFrame, Series))
        self.group_keys = group_keys
        self.obj = df[projection] if projection else df
        # grouping key passed via groupby method
        self.by = _normalize_by(df, by)
        self.sort = sort

        partitions_aligned = all(
            item.npartitions == df.npartitions if isinstance(item, Series) else True
            for item in (self.by if isinstance(self.by, (tuple, list)) else [self.by])
        )

        if not partitions_aligned:
            raise NotImplementedError(
                "The grouped object and 'by' of the groupby must have the same divisions."
            )

        if isinstance(self.by, list):
            by_meta = [
                item._meta if isinstance(item, Series) else item for item in self.by
            ]

        elif isinstance(self.by, Series):
            by_meta = self.by._meta
        else:
            by_meta = self.by

        self.dropna = {}
        if dropna is not None:
            self.dropna["dropna"] = dropna

        # Hold off on setting observed by default: https://github.com/dask/dask/issues/6951
        self.observed = {}
        if observed is not None:
            self.observed["observed"] = observed

        # raises a warning about observed=False with pandas>=2.1.
        # We want to raise here, and not later down the stack.
        self._meta = self.obj._meta.groupby(
            by_meta, group_keys=group_keys, **self.observed, **self.dropna
        )

    @property
    @_deprecated()
    def index(self):
        return self.by

    @index.setter
    def index(self, value):
        self.by = value

    @property
    def _groupby_kwargs(self):
        return {
            "by": self.by,
            "group_keys": self.group_keys,
            **self.dropna,
            "sort": self.sort,
            **self.observed,
        }

    def __iter__(self):
        raise NotImplementedError(
            "Iteration of DataFrameGroupBy objects requires computing the groups which "
            "may be slow. You probably want to use 'apply' to execute a function for "
            "all the columns. To access individual groups, use 'get_group'. To list "
            "all the group names, use 'df[<group column>].unique().compute()'."
        )

    @property
    def _meta_nonempty(self):
        """
        Return a pd.DataFrameGroupBy / pd.SeriesGroupBy which contains sample data.
        """
        sample = self.obj._meta_nonempty

        if isinstance(self.by, list):
            by_meta = [
                item._meta_nonempty if isinstance(item, Series) else item
                for item in self.by
            ]

        elif isinstance(self.by, Series):
            by_meta = self.by._meta_nonempty

        else:
            by_meta = self.by

        with check_observed_deprecation():
            grouped = sample.groupby(
                by_meta,
                group_keys=self.group_keys,
                **self.observed,
                **self.dropna,
            )
        return _maybe_slice(grouped, self._slice)

    def _single_agg(
        self,
        token,
        func,
        aggfunc=None,
        meta=None,
        split_every=None,
        split_out=1,
        shuffle_method=None,
        chunk_kwargs=None,
        aggregate_kwargs=None,
        columns=None,
    ):
        """
        Aggregation with a single function/aggfunc rather than a compound spec
        like in GroupBy.aggregate
        """
        shuffle_method = _determine_split_out_shuffle(shuffle_method, split_out)

        if aggfunc is None:
            aggfunc = func

        if chunk_kwargs is None:
            chunk_kwargs = {}

        if aggregate_kwargs is None:
            aggregate_kwargs = {}

        if meta is None:
            with check_numeric_only_deprecation():
                meta = func(self._meta_nonempty, **chunk_kwargs)

        if columns is None:
            columns = meta.name if is_series_like(meta) else meta.columns

        args = [self.obj] + (self.by if isinstance(self.by, list) else [self.by])

        token = self._token_prefix + token
        levels = _determine_levels(self.by)

        if shuffle_method:
            return _shuffle_aggregate(
                args,
                chunk=_apply_chunk,
                chunk_kwargs={
                    "chunk": func,
                    "columns": columns,
                    **self.observed,
                    **self.dropna,
                    **chunk_kwargs,
                },
                aggregate=_groupby_aggregate,
                aggregate_kwargs={
                    "aggfunc": aggfunc,
                    "levels": levels,
                    **self.observed,
                    **self.dropna,
                    **aggregate_kwargs,
                },
                token=token,
                split_every=split_every,
                split_out=split_out,
                shuffle_method=shuffle_method,
                sort=self.sort,
            )

        return aca(
            args,
            chunk=_apply_chunk,
            chunk_kwargs=dict(
                chunk=func,
                columns=columns,
                **self.observed,
                **chunk_kwargs,
                **self.dropna,
            ),
            aggregate=_groupby_aggregate,
            meta=meta,
            token=token,
            split_every=split_every,
            aggregate_kwargs=dict(
                aggfunc=aggfunc,
                levels=levels,
                **self.observed,
                **aggregate_kwargs,
                **self.dropna,
            ),
            split_out=split_out,
            split_out_setup=split_out_on_index,
            sort=self.sort,
        )

    def _cum_agg(self, token, chunk, aggregate, initial, numeric_only=no_default):
        """Wrapper for cumulative groupby operation"""
        numeric_only_kwargs = get_numeric_only_kwargs(numeric_only)
        meta = chunk(self._meta, **numeric_only_kwargs)
        columns = meta.name if is_series_like(meta) else meta.columns
        by_cols = self.by if isinstance(self.by, list) else [self.by]

        # rename "by" columns internally
        # to fix cumulative operations on the same "by" columns
        # ref: https://github.com/dask/dask/issues/9313
        if columns is not None:
            # Handle series case (only a single column, but can't
            # enlist above because that fails in pandas when
            # constructing meta later)
            grouping_columns = [columns] if is_series_like(meta) else columns
            to_rename = set(grouping_columns) & set(by_cols)
            by = []
            for col in by_cols:
                if col in to_rename:
                    suffix = str(uuid.uuid4())
                    self.obj = self.obj.assign(**{col + suffix: self.obj[col]})
                    by.append(col + suffix)
                else:
                    by.append(col)
        else:
            by = by_cols

        name = self._token_prefix + token
        name_part = name + "-map"
        name_last = name + "-take-last"
        name_cum = name + "-cum-last"

        # cumulate each partitions
        cumpart_raw = map_partitions(
            _apply_chunk,
            self.obj,
            *by,
            chunk=chunk,
            columns=columns,
            token=name_part,
            meta=meta,
            **self.dropna,
        )

        cumpart_raw_frame = (
            cumpart_raw.to_frame() if is_series_like(meta) else cumpart_raw
        )

        cumpart_ext = cumpart_raw_frame.assign(
            **{
                i: self.obj[i]
                if np.isscalar(i) and i in getattr(self.obj, "columns", [])
                else self.obj.index
                for i in by
            }
        )

        # Use pd.Grouper objects to specify that we are grouping by columns.
        # Otherwise, pandas will throw an ambiguity warning if the
        # DataFrame's index (self.obj.index) was included in the grouping
        # specification (self.by). See pandas #14432
        grouper = grouper_dispatch(self._meta.obj)
        by_groupers = [grouper(key=ind) for ind in by]
        cumlast = map_partitions(
            _apply_chunk,
            cumpart_ext,
            *by_groupers,
            columns=0 if columns is None else columns,
            chunk=M.last,
            meta=meta,
            token=name_last,
            **self.dropna,
        )

        # aggregate cumulated partitions and its previous last element
        _hash = tokenize(self, token, chunk, aggregate, initial)
        name += "-" + _hash
        name_cum += "-" + _hash
        dask = {}
        dask[(name, 0)] = (cumpart_raw._name, 0)

        for i in range(1, self.obj.npartitions):
            # store each cumulative step to graph to reduce computation
            if i == 1:
                dask[(name_cum, i)] = (cumlast._name, i - 1)
            else:
                # aggregate with previous cumulation results
                dask[(name_cum, i)] = (
                    _cum_agg_filled,
                    (name_cum, i - 1),
                    (cumlast._name, i - 1),
                    aggregate,
                    initial,
                )
            dask[(name, i)] = (
                _cum_agg_aligned,
                (cumpart_ext._name, i),
                (name_cum, i),
                by,
                0 if columns is None else columns,
                aggregate,
                initial,
            )

        dependencies = [cumpart_raw]
        if self.obj.npartitions > 1:
            dependencies += [cumpart_ext, cumlast]

        graph = HighLevelGraph.from_collections(name, dask, dependencies=dependencies)
        return new_dd_object(
            graph, name, chunk(self._meta, **numeric_only_kwargs), self.obj.divisions
        )

    def compute(self, **kwargs):
        raise NotImplementedError(
            "DataFrameGroupBy does not allow compute method."
            "Please chain it with an aggregation method (like ``.mean()``) or get a "
            "specific group using ``.get_group()`` before calling ``compute()``"
        )

    def _shuffle(self, meta):
        df = self.obj

        if isinstance(self.obj, Series):
            # Temporarily convert series to dataframe for shuffle
            df = df.to_frame("__series__")
            convert_back_to_series = True
        else:
            convert_back_to_series = False

        if isinstance(self.by, DataFrame):  # add by columns to dataframe
            df2 = df.assign(**{"_by_" + c: self.by[c] for c in self.by.columns})
        elif isinstance(self.by, Series):
            df2 = df.assign(_by=self.by)
        else:
            df2 = df

        df3 = df2.shuffle(on=self.by)  # shuffle dataframe and index

        if isinstance(self.by, DataFrame):
            # extract by from dataframe
            cols = ["_by_" + c for c in self.by.columns]
            by2 = df3[cols]
            if is_dataframe_like(meta):
                df4 = df3.map_partitions(drop_columns, cols, meta.columns.dtype)
            else:
                df4 = df3.drop(cols, axis=1)
        elif isinstance(self.by, Series):
            by2 = df3["_by"]
            by2.name = self.by.name
            if is_dataframe_like(meta):
                df4 = df3.map_partitions(drop_columns, "_by", meta.columns.dtype)
            else:
                df4 = df3.drop("_by", axis=1)
        else:
            df4 = df3
            by2 = self.by

        if convert_back_to_series:
            df4 = df4["__series__"].rename(self.obj.name)

        return df4, by2

    @_deprecated_kwarg("axis")
    @derived_from(pd.core.groupby.GroupBy)
    def cumsum(self, axis=no_default, numeric_only=no_default):
        axis = self._normalize_axis(axis, "cumsum")
        if axis:
            if axis in (1, "columns") and isinstance(self, SeriesGroupBy):
                raise ValueError(f"No axis named {axis} for object type Series")
            return self.obj.cumsum(axis=axis)
        else:
            return self._cum_agg(
                "cumsum",
                chunk=M.cumsum,
                aggregate=M.add,
                initial=0,
                numeric_only=numeric_only,
            )

    @_deprecated_kwarg("axis")
    @derived_from(pd.core.groupby.GroupBy)
    def cumprod(self, axis=no_default, numeric_only=no_default):
        axis = self._normalize_axis(axis, "cumprod")
        if axis:
            if axis in (1, "columns") and isinstance(self, SeriesGroupBy):
                raise ValueError(f"No axis named {axis} for object type Series")
            return self.obj.cumprod(axis=axis)
        else:
            return self._cum_agg(
                "cumprod",
                chunk=M.cumprod,
                aggregate=M.mul,
                initial=1,
                numeric_only=numeric_only,
            )

    @_deprecated_kwarg("axis")
    @derived_from(pd.core.groupby.GroupBy)
    def cumcount(self, axis=no_default):
        return self._cum_agg(
            "cumcount", chunk=M.cumcount, aggregate=_cumcount_aggregate, initial=-1
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.core.groupby.GroupBy)
    @numeric_only_deprecate_default
    def sum(
        self,
        split_every=None,
        split_out=1,
        shuffle_method=None,
        min_count=None,
        numeric_only=no_default,
    ):
        numeric_kwargs = get_numeric_only_kwargs(numeric_only)
        result = self._single_agg(
            func=M.sum,
            token="sum",
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            chunk_kwargs=numeric_kwargs,
            aggregate_kwargs=numeric_kwargs,
        )
        if min_count:
            return result.where(self.count() >= min_count, other=np.nan)
        else:
            return result

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.core.groupby.GroupBy)
    @numeric_only_deprecate_default
    def prod(
        self,
        split_every=None,
        split_out=1,
        shuffle_method=None,
        min_count=None,
        numeric_only=no_default,
    ):
        numeric_kwargs = get_numeric_only_kwargs(numeric_only)
        result = self._single_agg(
            func=M.prod,
            token="prod",
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            chunk_kwargs=numeric_kwargs,
            aggregate_kwargs=numeric_kwargs,
        )
        if min_count:
            return result.where(self.count() >= min_count, other=np.nan)
        else:
            return result

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.core.groupby.GroupBy)
    def min(
        self,
        split_every=None,
        split_out=1,
        shuffle_method=None,
        numeric_only=no_default,
    ):
        numeric_kwargs = get_numeric_only_kwargs(numeric_only)
        return self._single_agg(
            func=M.min,
            token="min",
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            chunk_kwargs=numeric_kwargs,
            aggregate_kwargs=numeric_kwargs,
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.core.groupby.GroupBy)
    def max(
        self,
        split_every=None,
        split_out=1,
        shuffle_method=None,
        numeric_only=no_default,
    ):
        numeric_kwargs = get_numeric_only_kwargs(numeric_only)
        return self._single_agg(
            func=M.max,
            token="max",
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            chunk_kwargs=numeric_kwargs,
            aggregate_kwargs=numeric_kwargs,
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.DataFrame)
    @numeric_only_deprecate_default
    def idxmin(
        self,
        split_every=None,
        split_out=1,
        shuffle_method=None,
        axis=no_default,
        skipna=True,
        numeric_only=no_default,
    ):
        if axis != no_default:
            warnings.warn(
                "`axis` parameter is deprecated and will be removed in a future version.",
                FutureWarning,
            )
        if axis in (1, "columns"):
            raise NotImplementedError(
                f"The axis={axis} keyword is not implemented for groupby.idxmin"
            )
        self._normalize_axis(axis, "idxmin")
        chunk_kwargs = dict(skipna=skipna)
        numeric_kwargs = get_numeric_only_kwargs(numeric_only)
        chunk_kwargs.update(numeric_kwargs)
        return self._single_agg(
            func=M.idxmin,
            token="idxmin",
            aggfunc=M.first,
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            chunk_kwargs=chunk_kwargs,
            aggregate_kwargs=numeric_kwargs,
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.DataFrame)
    @numeric_only_deprecate_default
    def idxmax(
        self,
        split_every=None,
        split_out=1,
        shuffle_method=None,
        axis=no_default,
        skipna=True,
        numeric_only=no_default,
    ):
        if axis != no_default:
            warnings.warn(
                "`axis` parameter is deprecated and will be removed in a future version.",
                FutureWarning,
            )

        if axis in (1, "columns"):
            raise NotImplementedError(
                f"The axis={axis} keyword is not implemented for groupby.idxmax"
            )
        self._normalize_axis(axis, "idxmax")
        chunk_kwargs = dict(skipna=skipna)
        numeric_kwargs = get_numeric_only_kwargs(numeric_only)
        chunk_kwargs.update(numeric_kwargs)
        return self._single_agg(
            func=M.idxmax,
            token="idxmax",
            aggfunc=M.first,
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            chunk_kwargs=chunk_kwargs,
            aggregate_kwargs=numeric_kwargs,
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.core.groupby.GroupBy)
    def count(self, split_every=None, split_out=1, shuffle_method=None):
        return self._single_agg(
            func=M.count,
            token="count",
            aggfunc=M.sum,
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.core.groupby.GroupBy)
    @numeric_only_not_implemented
    def mean(
        self,
        split_every=None,
        split_out=1,
        shuffle_method=None,
        numeric_only=no_default,
    ):
        # We sometimes emit this warning ourselves. We ignore it here so users only see it once.
        with check_numeric_only_deprecation():
            s = self.sum(
                split_every=split_every,
                split_out=split_out,
                shuffle_method=shuffle_method,
                numeric_only=numeric_only,
            )
        c = self.count(
            split_every=split_every, split_out=split_out, shuffle_method=shuffle_method
        )
        if is_dataframe_like(s):
            c = c[s.columns]
        return s / c

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.core.groupby.GroupBy)
    def median(
        self,
        split_every=None,
        split_out=1,
        shuffle_method=None,
        numeric_only=no_default,
    ):
        if shuffle_method is False:
            raise ValueError(
                "In order to aggregate with 'median', you must use shuffling-based "
                "aggregation (e.g., shuffle='tasks')"
            )

        shuffle_method = shuffle_method or _determine_split_out_shuffle(True, split_out)
        numeric_only_kwargs = get_numeric_only_kwargs(numeric_only)

        with check_numeric_only_deprecation(name="median"):
            meta = self._meta_nonempty.median(**numeric_only_kwargs)
        columns = meta.name if is_series_like(meta) else meta.columns
        by = self.by if isinstance(self.by, list) else [self.by]
        return _shuffle_aggregate(
            [self.obj] + by,
            token="non-agg",
            chunk=_non_agg_chunk,
            chunk_kwargs={
                "key": columns,
                **self.observed,
                **self.dropna,
            },
            aggregate=_groupby_aggregate,
            aggregate_kwargs={
                "aggfunc": _median_aggregate,
                "levels": _determine_levels(self.by),
                **self.observed,
                **self.dropna,
                **numeric_only_kwargs,
            },
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            sort=self.sort,
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.core.groupby.GroupBy)
    def size(self, split_every=None, split_out=1, shuffle_method=None):
        return self._single_agg(
            token="size",
            func=M.size,
            aggfunc=M.sum,
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
        )

    @derived_from(pd.core.groupby.GroupBy)
    @numeric_only_not_implemented
    def var(self, ddof=1, split_every=None, split_out=1, numeric_only=no_default):
        if not PANDAS_GE_150 and numeric_only is not no_default:
            raise TypeError("numeric_only not supported for pandas < 1.5")

        levels = _determine_levels(self.by)
        result = aca(
            [self.obj, self.by]
            if not isinstance(self.by, list)
            else [self.obj] + self.by,
            chunk=_var_chunk,
            aggregate=_var_agg,
            combine=_var_combine,
            token=self._token_prefix + "var",
            aggregate_kwargs={
                "ddof": ddof,
                "levels": levels,
                "numeric_only": numeric_only,
                **self.observed,
                **self.dropna,
            },
            chunk_kwargs={"numeric_only": numeric_only, **self.observed, **self.dropna},
            combine_kwargs={"levels": levels},
            split_every=split_every,
            split_out=split_out,
            split_out_setup=split_out_on_index,
            sort=self.sort,
        )

        if isinstance(self.obj, Series):
            result = result[result.columns[0]]
        if self._slice:
            result = result[self._slice]

        return result

    @derived_from(pd.core.groupby.GroupBy)
    @numeric_only_not_implemented
    def std(self, ddof=1, split_every=None, split_out=1, numeric_only=no_default):
        if not PANDAS_GE_150 and numeric_only is not no_default:
            raise TypeError("numeric_only not supported for pandas < 1.5")
        # We sometimes emit this warning ourselves. We ignore it here so users only see it once.
        with check_numeric_only_deprecation():
            v = self.var(
                ddof,
                split_every=split_every,
                split_out=split_out,
                numeric_only=numeric_only,
            )
        result = map_partitions(np.sqrt, v, meta=v)
        return result

    @derived_from(pd.DataFrame)
    def corr(self, ddof=1, split_every=None, split_out=1, numeric_only=no_default):
        """Groupby correlation:
        corr(X, Y) = cov(X, Y) / (std_x * std_y)
        """
        if not PANDAS_GE_150 and numeric_only is not no_default:
            raise TypeError("numeric_only not supported for pandas < 1.5")
        return self.cov(
            split_every=split_every,
            split_out=split_out,
            std=True,
            numeric_only=numeric_only,
        )

    @derived_from(pd.DataFrame)
    def cov(
        self, ddof=1, split_every=None, split_out=1, std=False, numeric_only=no_default
    ):
        """Groupby covariance is accomplished by

        1. Computing intermediate values for sum, count, and the product of
           all columns: a b c -> a*a, a*b, b*b, b*c, c*c.

        2. The values are then aggregated and the final covariance value is calculated:
           cov(X, Y) = X*Y - Xbar * Ybar

        When `std` is True calculate Correlation
        """
        if not PANDAS_GE_150 and numeric_only is not no_default:
            raise TypeError("numeric_only not supported for pandas < 1.5")
        numeric_only_kwargs = get_numeric_only_kwargs(numeric_only)

        levels = _determine_levels(self.by)

        is_mask = any(is_series_like(s) for s in self.by)
        if self._slice:
            if is_mask:
                self.obj = self.obj[self._slice]
            else:
                sliced_plus = list(self._slice) + list(self.by)
                self.obj = self.obj[sliced_plus]

        result = aca(
            [self.obj, self.by]
            if not isinstance(self.by, list)
            else [self.obj] + self.by,
            chunk=_cov_chunk,
            aggregate=_cov_agg,
            combine=_cov_combine,
            token=self._token_prefix + "cov",
            aggregate_kwargs={"ddof": ddof, "levels": levels, "std": std},
            combine_kwargs={"levels": levels},
            chunk_kwargs=numeric_only_kwargs,
            split_every=split_every,
            split_out=split_out,
            split_out_setup=split_out_on_index,
            sort=self.sort,
        )

        if isinstance(self.obj, Series):
            result = result[result.columns[0]]
        if self._slice:
            result = result[self._slice]
        return result

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.core.groupby.GroupBy)
    def first(
        self,
        split_every=None,
        split_out=1,
        shuffle_method=None,
        numeric_only=no_default,
    ):
        numeric_kwargs = get_numeric_only_kwargs(numeric_only)
        return self._single_agg(
            func=M.first,
            token="first",
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            chunk_kwargs=numeric_kwargs,
            aggregate_kwargs=numeric_kwargs,
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.core.groupby.GroupBy)
    def last(
        self,
        split_every=None,
        split_out=1,
        shuffle_method=None,
        numeric_only=no_default,
    ):
        numeric_kwargs = get_numeric_only_kwargs(numeric_only)
        return self._single_agg(
            token="last",
            func=M.last,
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            chunk_kwargs=numeric_kwargs,
            aggregate_kwargs=numeric_kwargs,
        )

    @derived_from(
        pd.core.groupby.GroupBy,
        inconsistencies="If the group is not present, Dask will return an empty Series/DataFrame.",
    )
    def get_group(self, key):
        token = self._token_prefix + "get_group"

        meta = self._meta.obj
        if is_dataframe_like(meta) and self._slice is not None:
            meta = meta[self._slice]
        columns = meta.columns if is_dataframe_like(meta) else meta.name

        return map_partitions(
            _groupby_get_group,
            self.obj,
            self.by,
            key,
            columns,
            meta=meta,
            token=token,
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @_aggregate_docstring()
    def aggregate(
        self, arg=None, split_every=None, split_out=1, shuffle_method=None, **kwargs
    ):
        if split_out is None:
            warnings.warn(
                "split_out=None is deprecated, please use a positive integer, "
                "or allow the default of 1",
                category=FutureWarning,
            )
            split_out = 1
        shuffle_method = _determine_split_out_shuffle(shuffle_method, split_out)

        relabeling = None
        columns = None
        order = None
        column_projection = None
        if PANDAS_GE_140:
            if isinstance(self, DataFrameGroupBy):
                if arg is None:
                    relabeling, arg, columns, order = reconstruct_func(arg, **kwargs)

            elif isinstance(self, SeriesGroupBy):
                relabeling = arg is None
                if relabeling:
                    columns, arg = validate_func_kwargs(kwargs)

        if isinstance(self.obj, DataFrame):
            if isinstance(self.by, tuple) or np.isscalar(self.by):
                group_columns = {self.by}

            elif isinstance(self.by, list):
                group_columns = {
                    i for i in self.by if isinstance(i, tuple) or np.isscalar(i)
                }

            else:
                group_columns = set()

            if self._slice:
                # pandas doesn't exclude the grouping column in a SeriesGroupBy
                # like df.groupby('a')['a'].agg(...)
                non_group_columns = self._slice
                if not isinstance(non_group_columns, list):
                    non_group_columns = [non_group_columns]
            else:
                # NOTE: this step relies on the by normalization to replace
                #       series with their name.
                non_group_columns = [
                    col for col in self.obj.columns if col not in group_columns
                ]

            spec = _normalize_spec(arg, non_group_columns)

            # Check if the aggregation involves implicit column projection
            if isinstance(arg, dict):
                column_projection = group_columns.union(arg.keys()).intersection(
                    self.obj.columns
                )

        elif isinstance(self.obj, Series):
            if isinstance(arg, (list, tuple, dict)):
                # implementation detail: if self.obj is a series, a pseudo column
                # None is used to denote the series itself. This pseudo column is
                # removed from the result columns before passing the spec along.
                spec = _normalize_spec({None: arg}, [])
                spec = [
                    (result_column, func, input_column)
                    for ((_, result_column), func, input_column) in spec
                ]

            else:
                spec = _normalize_spec({None: arg}, [])
                spec = [
                    (self.obj.name, func, input_column)
                    for (_, func, input_column) in spec
                ]

        else:
            raise ValueError(f"aggregate on unknown object {self.obj}")

        chunk_funcs, aggregate_funcs, finalizers = _build_agg_args(spec)

        if isinstance(self.by, (tuple, list)) and len(self.by) > 1:
            levels = list(range(len(self.by)))
        else:
            levels = 0

        # Add an explicit `getitem` operation if the groupby
        # aggregation involves implicit column projection.
        # This makes it possible for the column-projection
        # to be pushed into the IO layer
        _obj = self.obj[list(column_projection)] if column_projection else self.obj

        if not isinstance(self.by, list):
            chunk_args = [_obj, self.by]

        else:
            chunk_args = [_obj] + self.by

        # If any of the agg funcs contain a "median", we *must* use the shuffle
        # implementation.
        has_median = any(s[1] in ("median", np.median) for s in spec)
        if has_median and not shuffle_method:
            raise ValueError(
                "In order to aggregate with 'median', you must use shuffling-based "
                "aggregation (e.g., shuffle='tasks')"
            )

        if shuffle_method:
            # Shuffle-based aggregation
            #
            # This algorithm is more scalable than a tree reduction
            # for larger values of split_out. However, the shuffle
            # step requires that the result of `chunk` produces a
            # proper DataFrame type

            # If we have a median in the spec, we cannot do an initial
            # aggregation.
            if has_median:
                result = _shuffle_aggregate(
                    chunk_args,
                    chunk=_non_agg_chunk,
                    chunk_kwargs={
                        "key": [
                            c for c in _obj.columns.tolist() if c not in group_columns
                        ],
                        **self.observed,
                        **self.dropna,
                    },
                    aggregate=_groupby_aggregate_spec,
                    aggregate_kwargs={
                        "spec": arg,
                        "levels": _determine_levels(self.by),
                        **self.observed,
                        **self.dropna,
                    },
                    token="aggregate",
                    split_every=split_every,
                    split_out=split_out,
                    shuffle_method=shuffle_method,
                    sort=self.sort,
                )
            else:
                result = _shuffle_aggregate(
                    chunk_args,
                    chunk=_groupby_apply_funcs,
                    chunk_kwargs={
                        "funcs": chunk_funcs,
                        "sort": self.sort,
                        **self.observed,
                        **self.dropna,
                    },
                    aggregate=_agg_finalize,
                    aggregate_kwargs=dict(
                        aggregate_funcs=aggregate_funcs,
                        finalize_funcs=finalizers,
                        level=levels,
                        **self.observed,
                        **self.dropna,
                    ),
                    token="aggregate",
                    split_every=split_every,
                    split_out=split_out,
                    shuffle_method=shuffle_method,
                    sort=self.sort,
                )
        else:
            # Check sort behavior
            if self.sort and split_out > 1:
                raise NotImplementedError(
                    "Cannot guarantee sorted keys for `split_out>1` and `shuffle=False`"
                    " Try using `shuffle=True` if you are grouping on a single column."
                    " Otherwise, try using split_out=1, or grouping with sort=False."
                )

            result = aca(
                chunk_args,
                chunk=_groupby_apply_funcs,
                chunk_kwargs=dict(
                    funcs=chunk_funcs,
                    sort=False,
                    **self.observed,
                    **self.dropna,
                ),
                combine=_groupby_apply_funcs,
                combine_kwargs=dict(
                    funcs=aggregate_funcs,
                    level=levels,
                    sort=False,
                    **self.observed,
                    **self.dropna,
                ),
                aggregate=_agg_finalize,
                aggregate_kwargs=dict(
                    aggregate_funcs=aggregate_funcs,
                    finalize_funcs=finalizers,
                    level=levels,
                    **self.observed,
                    **self.dropna,
                ),
                token="aggregate",
                split_every=split_every,
                split_out=split_out,
                split_out_setup=split_out_on_index,
                sort=self.sort,
            )

        if relabeling and result is not None:
            if order is not None:
                result = result.iloc[:, order]
            result.columns = columns

        return result

    @insert_meta_param_description(pad=12)
    def apply(self, func, *args, **kwargs):
        """Parallel version of pandas GroupBy.apply

        This mimics the pandas version except for the following:

        1.  If the grouper does not align with the index then this causes a full
            shuffle.  The order of rows within each group may not be preserved.
        2.  Dask's GroupBy.apply is not appropriate for aggregations. For custom
            aggregations, use :class:`dask.dataframe.groupby.Aggregation`.

        .. warning::

           Pandas' groupby-apply can be used to to apply arbitrary functions,
           including aggregations that result in one row per group. Dask's
           groupby-apply will apply ``func`` once on each group, doing a shuffle
           if needed, such that each group is contained in one partition.
           When ``func`` is a reduction, e.g., you'll end up with one row
           per group. To apply a custom aggregation with Dask,
           use :class:`dask.dataframe.groupby.Aggregation`.

        Parameters
        ----------
        func: function
            Function to apply
        args, kwargs : Scalar, Delayed or object
            Arguments and keywords to pass to the function.
        $META

        Returns
        -------
        applied : Series or DataFrame depending on columns keyword
        """
        meta = kwargs.get("meta", no_default)

        if meta is no_default:
            with raise_on_meta_error(f"groupby.apply({funcname(func)})", udf=True):
                meta_args, meta_kwargs = _extract_meta((args, kwargs), nonempty=True)
                meta = self._meta_nonempty.apply(func, *meta_args, **meta_kwargs)

            msg = (
                "`meta` is not specified, inferred from partial data. "
                "Please provide `meta` if the result is unexpected.\n"
                "  Before: .apply(func)\n"
                "  After:  .apply(func, meta={'x': 'f8', 'y': 'f8'}) for dataframe result\n"
                "  or:     .apply(func, meta=('x', 'f8'))            for series result"
            )
            warnings.warn(msg, stacklevel=2)

        meta = make_meta(meta, parent_meta=self._meta.obj)

        # Validate self.by
        if isinstance(self.by, list) and any(
            isinstance(item, Series) for item in self.by
        ):
            raise NotImplementedError(
                "groupby-apply with a multiple Series is currently not supported"
            )

        df = self.obj
        should_shuffle = not (df.known_divisions and df._contains_index_name(self.by))

        if should_shuffle:
            df2, by = self._shuffle(meta)
        else:
            df2 = df
            by = self.by

        # Perform embarrassingly parallel groupby-apply
        kwargs["meta"] = meta
        df3 = map_partitions(
            _groupby_slice_apply,
            df2,
            by,
            self._slice,
            func,
            *args,
            token=funcname(func),
            group_keys=self.group_keys,
            **self.observed,
            **self.dropna,
            **kwargs,
        )

        return df3

    @insert_meta_param_description(pad=12)
    def transform(self, func, *args, **kwargs):
        """Parallel version of pandas GroupBy.transform

        This mimics the pandas version except for the following:

        1.  If the grouper does not align with the index then this causes a full
            shuffle.  The order of rows within each group may not be preserved.
        2.  Dask's GroupBy.transform is not appropriate for aggregations. For custom
            aggregations, use :class:`dask.dataframe.groupby.Aggregation`.

        .. warning::

           Pandas' groupby-transform can be used to apply arbitrary functions,
           including aggregations that result in one row per group. Dask's
           groupby-transform will apply ``func`` once on each group, doing a shuffle
           if needed, such that each group is contained in one partition.
           When ``func`` is a reduction, e.g., you'll end up with one row
           per group. To apply a custom aggregation with Dask,
           use :class:`dask.dataframe.groupby.Aggregation`.

        Parameters
        ----------
        func: function
            Function to apply
        args, kwargs : Scalar, Delayed or object
            Arguments and keywords to pass to the function.
        $META

        Returns
        -------
        applied : Series or DataFrame depending on columns keyword
        """
        meta = kwargs.get("meta", no_default)

        if meta is no_default:
            with raise_on_meta_error(f"groupby.transform({funcname(func)})", udf=True):
                meta_args, meta_kwargs = _extract_meta((args, kwargs), nonempty=True)
                meta = self._meta_nonempty.transform(func, *meta_args, **meta_kwargs)

            msg = (
                "`meta` is not specified, inferred from partial data. "
                "Please provide `meta` if the result is unexpected.\n"
                "  Before: .transform(func)\n"
                "  After:  .transform(func, meta={'x': 'f8', 'y': 'f8'}) for dataframe result\n"
                "  or:     .transform(func, meta=('x', 'f8'))            for series result"
            )
            warnings.warn(msg, stacklevel=2)

        meta = make_meta(meta, parent_meta=self._meta.obj)

        # Validate self.by
        if isinstance(self.by, list) and any(
            isinstance(item, Series) for item in self.by
        ):
            raise NotImplementedError(
                "groupby-transform with a multiple Series is currently not supported"
            )

        df = self.obj
        should_shuffle = not (df.known_divisions and df._contains_index_name(self.by))

        if should_shuffle:
            df2, by = self._shuffle(meta)
        else:
            df2 = df
            by = self.by

        # Perform embarrassingly parallel groupby-transform
        kwargs["meta"] = meta
        df3 = map_partitions(
            _groupby_slice_transform,
            df2,
            by,
            self._slice,
            func,
            *args,
            token=funcname(func),
            group_keys=self.group_keys,
            **self.observed,
            **self.dropna,
            **kwargs,
        )

        if isinstance(self, DataFrameGroupBy):
            index_name = df3.index.name
            df3 = df3.reset_index().set_index(index_name or "index")
            df3.index = df3.index.rename(index_name)

        return df3

    @insert_meta_param_description(pad=12)
    def shift(
        self,
        periods=1,
        freq=no_default,
        axis=no_default,
        fill_value=no_default,
        meta=no_default,
    ):
        """Parallel version of pandas GroupBy.shift

        This mimics the pandas version except for the following:

        If the grouper does not align with the index then this causes a full
        shuffle.  The order of rows within each group may not be preserved.

        Parameters
        ----------
        periods : Delayed, Scalar or int, default 1
            Number of periods to shift.
        freq : Delayed, Scalar or str, optional
            Frequency string.
        axis : axis to shift, default 0
            Shift direction.
        fill_value : Scalar, Delayed or object, optional
            The scalar value to use for newly introduced missing values.
        $META

        Returns
        -------
        shifted : Series or DataFrame shifted within each group.

        Examples
        --------
        >>> import dask
        >>> ddf = dask.datasets.timeseries(freq="1h")
        >>> result = ddf.groupby("name").shift(1, meta={"id": int, "x": float, "y": float})
        """
        if axis != no_default:
            warnings.warn(
                "`axis` parameter is deprecated and will be removed in a future version.",
                FutureWarning,
            )
        axis = self._normalize_axis(axis, "shift")
        kwargs = {"periods": periods, "axis": axis}
        if freq is not no_default:
            kwargs.update({"freq": freq})
        if fill_value is not no_default:
            kwargs.update({"fill_value": fill_value})
        if meta is no_default:
            with raise_on_meta_error("groupby.shift()", udf=False):
                meta_kwargs = _extract_meta(
                    kwargs,
                    nonempty=True,
                )
                with check_groupby_axis_deprecation():
                    meta = self._meta_nonempty.shift(**meta_kwargs)

            msg = (
                "`meta` is not specified, inferred from partial data. "
                "Please provide `meta` if the result is unexpected.\n"
                "  Before: .shift(1)\n"
                "  After:  .shift(1, meta={'x': 'f8', 'y': 'f8'}) for dataframe result\n"
                "  or:     .shift(1, meta=('x', 'f8'))            for series result"
            )
            warnings.warn(msg, stacklevel=2)

        meta = make_meta(meta, parent_meta=self._meta.obj)

        # Validate self.by
        if isinstance(self.by, list) and any(
            isinstance(item, Series) for item in self.by
        ):
            raise NotImplementedError(
                "groupby-shift with a multiple Series is currently not supported"
            )
        df = self.obj
        should_shuffle = not (df.known_divisions and df._contains_index_name(self.by))

        if should_shuffle:
            df2, by = self._shuffle(meta)
        else:
            df2 = df
            by = self.by

        # Perform embarrassingly parallel groupby-shift
        result = map_partitions(
            _groupby_slice_shift,
            df2,
            by,
            self._slice,
            should_shuffle,
            token="groupby-shift",
            group_keys=self.group_keys,
            meta=meta,
            **self.observed,
            **self.dropna,
            **kwargs,
        )
        return result

    def rolling(self, window, min_periods=None, center=False, win_type=None, axis=0):
        """Provides rolling transformations.

        .. note::

            Since MultiIndexes are not well supported in Dask, this method returns a
            dataframe with the same index as the original data. The groupby column is
            not added as the first level of the index like pandas does.

            This method works differently from other groupby methods. It does a groupby
            on each partition (plus some overlap). This means that the output has the
            same shape and number of partitions as the original.

        Parameters
        ----------
        window : str, offset
           Size of the moving window. This is the number of observations used
           for calculating the statistic. Data must have a ``DatetimeIndex``
        min_periods : int, default None
            Minimum number of observations in window required to have a value
            (otherwise result is NA).
        center : boolean, default False
            Set the labels at the center of the window.
        win_type : string, default None
            Provide a window type. The recognized window types are identical
            to pandas.
        axis : int, default 0

        Returns
        -------
        a Rolling object on which to call a method to compute a statistic

        Examples
        --------
        >>> import dask
        >>> ddf = dask.datasets.timeseries(freq="1h")
        >>> result = ddf.groupby("name").x.rolling('1D').max()
        """
        from dask.dataframe.rolling import RollingGroupby

        if isinstance(window, Integral):
            raise ValueError(
                "Only time indexes are supported for rolling groupbys in dask dataframe. "
                "``window`` must be a ``freq`` (e.g. '1H')."
            )

        if min_periods is not None:
            if not isinstance(min_periods, Integral):
                raise ValueError("min_periods must be an integer")
            if min_periods < 0:
                raise ValueError("min_periods must be >= 0")

        return RollingGroupby(
            self,
            window=window,
            min_periods=min_periods,
            center=center,
            win_type=win_type,
            axis=axis,
        )

    def _normalize_axis(self, axis, method: str):
        if PANDAS_GE_210 and axis is not no_default:
            if axis in (0, "index"):
                warnings.warn(
                    f"The 'axis' keyword in {type(self).__name__}.{method} is deprecated and will "
                    "be removed in a future version. Call without passing 'axis' instead.",
                    FutureWarning,
                )
            else:
                warnings.warn(
                    f"{type(self).__name__}.{method} with axis={axis} is deprecated and will be removed "
                    "in a future version. Operate on the un-grouped DataFrame instead",
                    FutureWarning,
                )
        if axis is no_default:
            axis = 0

        if axis in ("index", 1):
            warnings.warn(
                "Using axis=1 in GroupBy does not require grouping and will be removed "
                "entirely in a future version.",
                FutureWarning,
            )

        return axis

    @_deprecated(message="Please use `ffill`/`bfill` or `fillna` without a GroupBy.")
    def fillna(self, value=None, method=None, limit=None, axis=no_default):
        """Fill NA/NaN values using the specified method.

        Parameters
        ----------
        value : scalar, default None
            Value to use to fill holes (e.g. 0).
        method : {'bfill', 'ffill', None}, default None
            Method to use for filling holes in reindexed Series. ffill: propagate last
            valid observation forward to next valid. bfill: use next valid observation
            to fill gap.
        axis : {0 or 'index', 1 or 'columns'}
            Axis along which to fill missing values.
        limit : int, default None
            If method is specified, this is the maximum number of consecutive NaN values
            to forward/backward fill. In other words, if there is a gap with more than
            this number of consecutive NaNs, it will only be partially filled. If method
            is not specified, this is the maximum number of entries along the entire
            axis where NaNs will be filled. Must be greater than 0 if not None.

        Returns
        -------
        Series or DataFrame
            Object with missing values filled

        See also
        --------
        pandas.core.groupby.DataFrameGroupBy.fillna
        """
        axis = self._normalize_axis(axis, "fillna")
        if not np.isscalar(value) and value is not None:
            raise NotImplementedError(
                "groupby-fillna with value=dict/Series/DataFrame is not supported"
            )

        kwargs = dict(value=value, method=method, limit=limit, axis=axis)
        if PANDAS_GE_220:
            func = M.fillna
            kwargs.update(include_groups=False)
        else:
            func = _drop_apply
            kwargs.update(by=self.by, what="fillna")

        meta = self._meta_nonempty.apply(func, **kwargs)
        result = self.apply(func, meta=meta, **kwargs)

        if PANDAS_GE_150 and self.group_keys:
            return result.map_partitions(M.droplevel, self.by)

        return result

    @derived_from(pd.core.groupby.GroupBy)
    def ffill(self, limit=None):
        kwargs = dict(limit=limit)
        if PANDAS_GE_220:
            func = M.ffill
            kwargs.update(include_groups=False)
        else:
            func = _drop_apply
            kwargs.update(by=self.by, what="ffill")

        meta = self._meta_nonempty.apply(func, **kwargs)
        result = self.apply(func, meta=meta, **kwargs)

        if PANDAS_GE_150 and self.group_keys:
            return result.map_partitions(M.droplevel, self.by)
        return result

    @derived_from(pd.core.groupby.GroupBy)
    def bfill(self, limit=None):
        kwargs = dict(limit=limit)
        if PANDAS_GE_220:
            func = M.bfill
            kwargs.update(include_groups=False)
        else:
            func = _drop_apply
            kwargs.update(by=self.by, what="bfill")

        meta = self._meta_nonempty.apply(func, **kwargs)
        result = self.apply(func, meta=meta, **kwargs)

        if PANDAS_GE_150 and self.group_keys:
            return result.map_partitions(M.droplevel, self.by)
        return result


class DataFrameGroupBy(_GroupBy):
    _token_prefix = "dataframe-groupby-"

    def __getitem__(self, key):
        with check_observed_deprecation():
            if isinstance(key, list):
                g = DataFrameGroupBy(
                    self.obj,
                    by=self.by,
                    slice=key,
                    sort=self.sort,
                    **self.dropna,
                    **self.observed,
                )
            else:
                g = SeriesGroupBy(
                    self.obj,
                    by=self.by,
                    slice=key,
                    sort=self.sort,
                    **self.dropna,
                    **self.observed,
                )

        # Need a list otherwise pandas will warn/error
        if isinstance(key, tuple):
            key = list(key)
        g._meta = g._meta[key]
        return g

    def __dir__(self):
        return sorted(
            set(
                dir(type(self))
                + list(self.__dict__)
                + list(filter(M.isidentifier, self.obj.columns))
            )
        )

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError as e:
            raise AttributeError(e) from e

    def _all_numeric(self):
        """Are all columns that we're not grouping on numeric?"""
        numerics = self.obj._meta._get_numeric_data()
        # This computes a groupby but only on the empty meta
        post_group_columns = self._meta.count().columns
        return len(set(post_group_columns) - set(numerics.columns)) == 0

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @_aggregate_docstring(based_on="pd.core.groupby.DataFrameGroupBy.aggregate")
    def aggregate(
        self, arg=None, split_every=None, split_out=1, shuffle_method=None, **kwargs
    ):
        if arg == "size":
            return self.size()

        return super().aggregate(
            arg=arg,
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            **kwargs,
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @_aggregate_docstring(based_on="pd.core.groupby.DataFrameGroupBy.agg")
    @numeric_only_not_implemented
    def agg(
        self, arg=None, split_every=None, split_out=1, shuffle_method=None, **kwargs
    ):
        return self.aggregate(
            arg=arg,
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            **kwargs,
        )


class SeriesGroupBy(_GroupBy):
    _token_prefix = "series-groupby-"

    def __init__(self, df, by=None, slice=None, observed=None, **kwargs):
        # for any non series object, raise pandas-compat error message
        # Hold off on setting observed by default: https://github.com/dask/dask/issues/6951
        observed = {"observed": observed} if observed is not None else {}

        if isinstance(df, Series):
            if isinstance(by, Series):
                pass
            elif isinstance(by, list):
                if len(by) == 0:
                    raise ValueError("No group keys passed!")

                non_series_items = [item for item in by if not isinstance(item, Series)]
                # raise error from pandas, if applicable

                df._meta.groupby(non_series_items, **observed)
            else:
                # raise error from pandas, if applicable
                df._meta.groupby(by, **observed)

        super().__init__(df, by=by, slice=slice, **observed, **kwargs)

    @derived_from(pd.core.groupby.SeriesGroupBy)
    def nunique(self, split_every=None, split_out=1):
        """
        Examples
        --------
        >>> import pandas as pd
        >>> import dask.dataframe as dd
        >>> d = {'col1': [1, 2, 3, 4], 'col2': [5, 6, 7, 8]}
        >>> df = pd.DataFrame(data=d)
        >>> ddf = dd.from_pandas(df, 2)
        >>> ddf.groupby(['col1']).col2.nunique().compute()
        """
        name = self._meta.obj.name
        levels = _determine_levels(self.by)

        if isinstance(self.obj, DataFrame):
            chunk = _nunique_df_chunk

        else:
            chunk = _nunique_series_chunk

        return aca(
            [self.obj, self.by]
            if not isinstance(self.by, list)
            else [self.obj] + self.by,
            chunk=chunk,
            aggregate=_nunique_df_aggregate,
            combine=_nunique_df_combine,
            token="series-groupby-nunique",
            chunk_kwargs={"levels": levels, "name": name},
            aggregate_kwargs={"levels": levels, "name": name},
            combine_kwargs={"levels": levels},
            split_every=split_every,
            split_out=split_out,
            split_out_setup=split_out_on_index,
            sort=self.sort,
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @_aggregate_docstring(based_on="pd.core.groupby.SeriesGroupBy.aggregate")
    def aggregate(
        self, arg=None, split_every=None, split_out=1, shuffle_method=None, **kwargs
    ):
        result = super().aggregate(
            arg=arg,
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            **kwargs,
        )
        if self._slice:
            try:
                result = result[self._slice]
            except KeyError:
                pass

        if (
            arg is not None
            and not isinstance(arg, (list, dict))
            and isinstance(result, DataFrame)
        ):
            result = result[result.columns[0]]

        return result

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @_aggregate_docstring(based_on="pd.core.groupby.SeriesGroupBy.agg")
    def agg(
        self, arg=None, split_every=None, split_out=1, shuffle_method=None, **kwargs
    ):
        return self.aggregate(
            arg=arg,
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            **kwargs,
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.core.groupby.SeriesGroupBy)
    def value_counts(self, split_every=None, split_out=1, shuffle_method=None):
        return self._single_agg(
            func=_value_counts,
            token="value_counts",
            aggfunc=_value_counts_aggregate,
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            # in pandas 2.0, Series returned from value_counts have a name
            # different from original object, but here, column name should
            # still reflect the original object name
            columns=self._meta.apply(pd.Series).name,
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.core.groupby.SeriesGroupBy)
    def unique(self, split_every=None, split_out=1, shuffle_method=None):
        name = self._meta.obj.name
        return self._single_agg(
            func=M.unique,
            token="unique",
            aggfunc=_unique_aggregate,
            aggregate_kwargs={"name": name},
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.core.groupby.SeriesGroupBy)
    def tail(self, n=5, split_every=None, split_out=1, shuffle_method=None):
        index_levels = len(self.by) if isinstance(self.by, list) else 1
        return self._single_agg(
            func=_tail_chunk,
            token="tail",
            aggfunc=_tail_aggregate,
            meta=M.tail(self._meta_nonempty),
            chunk_kwargs={"n": n},
            aggregate_kwargs={"n": n, "index_levels": index_levels},
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.core.groupby.SeriesGroupBy)
    def head(self, n=5, split_every=None, split_out=1, shuffle_method=None):
        index_levels = len(self.by) if isinstance(self.by, list) else 1
        return self._single_agg(
            func=_head_chunk,
            token="head",
            aggfunc=_head_aggregate,
            meta=M.head(self._meta_nonempty),
            chunk_kwargs={"n": n},
            aggregate_kwargs={"n": n, "index_levels": index_levels},
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
        )


def _unique_aggregate(series_gb, name=None):
    data = {k: v.explode().unique() for k, v in series_gb}
    ret = type(series_gb.obj)(data, name=name)
    ret.index.names = series_gb.obj.index.names
    ret.index = ret.index.astype(series_gb.obj.index.dtype, copy=False)
    return ret


def _value_counts(x, **kwargs):
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


def _median_aggregate(series_gb, **kwargs):
    with check_numeric_only_deprecation():
        return series_gb.median(**kwargs)


def _shuffle_aggregate(
    args,
    chunk=None,
    aggregate=None,
    token=None,
    chunk_kwargs=None,
    aggregate_kwargs=None,
    split_every=None,
    split_out=1,
    sort=True,
    ignore_index=False,
    shuffle_method="tasks",
):
    """Shuffle-based groupby aggregation

    This algorithm may be more efficient than ACA for large ``split_out``
    values (required for high-cardinality groupby indices), but it also
    requires the output of ``chunk`` to be a proper DataFrame object.

    Parameters
    ----------
    args :
        Positional arguments for the `chunk` function. All `dask.dataframe`
        objects should be partitioned and indexed equivalently.
    chunk : function [block-per-arg] -> block
        Function to operate on each block of data
    aggregate : function concatenated-block -> block
        Function to operate on the concatenated result of chunk
    token : str, optional
        The name to use for the output keys.
    chunk_kwargs : dict, optional
        Keywords for the chunk function only.
    aggregate_kwargs : dict, optional
        Keywords for the aggregate function only.
    split_every : int, optional
        Number of partitions to aggregate into a shuffle partition.
        Defaults to eight, meaning that the initial partitions are repartitioned
        into groups of eight before the shuffle. Shuffling scales with the number
        of partitions, so it may be helpful to increase this number as a performance
        optimization, but only when the aggregated partition can comfortably
        fit in worker memory.
    split_out : int, optional
        Number of output partitions.
    ignore_index : bool, default False
        Whether the index can be ignored during the shuffle.
    sort : bool
        If allowed, sort the keys of the output aggregation.
    shuffle_method : str, default "tasks"
        Shuffle method to be used by ``DataFrame.shuffle``.
    """

    if chunk_kwargs is None:
        chunk_kwargs = dict()
    if aggregate_kwargs is None:
        aggregate_kwargs = dict()

    if not isinstance(args, (tuple, list)):
        args = [args]

    dfs = [arg for arg in args if isinstance(arg, _Frame)]

    npartitions = {arg.npartitions for arg in dfs}
    if len(npartitions) > 1:
        raise ValueError("All arguments must have same number of partitions")
    npartitions = npartitions.pop()

    if split_every is None:
        split_every = 8
    elif split_every is False:
        split_every = npartitions
    elif split_every < 1 or not isinstance(split_every, Integral):
        raise ValueError("split_every must be an integer >= 1")

    # Shuffle-based groupby aggregation. N.B. we have to use `_meta_nonempty`
    # as some chunk aggregation functions depend on the index having values to
    # determine `group_keys` behavior.
    chunk_name = f"{token or funcname(chunk)}-chunk"
    chunked = map_partitions(
        chunk,
        *args,
        meta=chunk(
            *[arg._meta_nonempty if isinstance(arg, _Frame) else arg for arg in args],
            **chunk_kwargs,
        ),
        token=chunk_name,
        **chunk_kwargs,
    )
    if is_series_like(chunked):
        # Temporarily convert series to dataframe for shuffle
        series_name = chunked._meta.name
        chunked = chunked.to_frame("__series__")
        convert_back_to_series = True
    else:
        series_name = None
        convert_back_to_series = False

    shuffle_npartitions = max(
        chunked.npartitions // split_every,
        split_out,
    )

    # Handle sort kwarg
    if sort is not None:
        aggregate_kwargs = aggregate_kwargs or {}
        aggregate_kwargs["sort"] = sort

    if sort is None and split_out > 1:
        idx = set(chunked._meta.columns) - set(chunked._meta.reset_index().columns)
        if len(idx) > 1:
            warnings.warn(
                "In the future, `sort` for groupby operations will default to `True`"
                " to match the behavior of pandas. However, `sort=True` can have "
                " significant performance implications when `split_out>1`. To avoid "
                " global data shuffling, set `sort=False`.",
                FutureWarning,
            )

    # Perform global sort or shuffle
    if sort and split_out > 1:
        cols = set(chunked.columns)
        chunked = chunked.reset_index()
        index_cols = sorted(set(chunked.columns) - cols)
        if len(index_cols) > 1:
            # Cannot use `set_index` for multi-column sort
            result = chunked.sort_values(
                index_cols,
                npartitions=shuffle_npartitions,
                shuffle_method=shuffle_method,
            ).map_partitions(
                M.set_index,
                index_cols,
                meta=chunked._meta.set_index(list(index_cols)),
                enforce_metadata=False,
            )
        else:
            result = chunked.set_index(
                index_cols,
                npartitions=shuffle_npartitions,
                shuffle_method=shuffle_method,
            )
    else:
        result = chunked.shuffle(
            chunked.index,
            ignore_index=ignore_index,
            npartitions=shuffle_npartitions,
            shuffle_method=shuffle_method,
        )

    # Aggregate
    result = result.map_partitions(aggregate, **aggregate_kwargs)

    if convert_back_to_series:
        result = result["__series__"].rename(series_name)

    if split_out < shuffle_npartitions:
        return result.repartition(npartitions=split_out)
    return result
