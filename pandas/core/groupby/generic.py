"""
Define the SeriesGroupBy and DataFrameGroupBy
classes that hold the groupby interfaces (and some implementations).

These are user facing as the result of the ``df.groupby(...)`` operations,
which here returns a DataFrameGroupBy object.
"""
from collections import abc, namedtuple
import copy
from functools import partial
from textwrap import dedent
import typing
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    FrozenSet,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Type,
    Union,
    cast,
)
import warnings

import numpy as np

from pandas._libs import lib
from pandas._typing import FrameOrSeries, FrameOrSeriesUnion
from pandas.util._decorators import Appender, Substitution, doc

from pandas.core.dtypes.cast import (
    maybe_cast_result,
    maybe_cast_result_dtype,
    maybe_convert_objects,
    maybe_downcast_numeric,
    maybe_downcast_to_dtype,
)
from pandas.core.dtypes.common import (
    ensure_int64,
    ensure_platform_int,
    is_bool,
    is_integer_dtype,
    is_interval_dtype,
    is_numeric_dtype,
    is_object_dtype,
    is_scalar,
    needs_i8_conversion,
)
from pandas.core.dtypes.missing import isna, notna

from pandas.core.aggregation import (
    maybe_mangle_lambdas,
    reconstruct_func,
    validate_func_kwargs,
)
import pandas.core.algorithms as algorithms
from pandas.core.base import DataError, SpecificationError
import pandas.core.common as com
from pandas.core.construction import create_series_with_explicit_dtype
from pandas.core.frame import DataFrame
from pandas.core.generic import ABCDataFrame, ABCSeries, NDFrame
from pandas.core.groupby import base
from pandas.core.groupby.groupby import (
    GroupBy,
    _agg_template,
    _apply_docs,
    _transform_template,
    get_groupby,
)
from pandas.core.indexes.api import Index, MultiIndex, all_indexes_same
import pandas.core.indexes.base as ibase
from pandas.core.internals import BlockManager, make_block
from pandas.core.series import Series
from pandas.core.util.numba_ import (
    NUMBA_FUNC_CACHE,
    generate_numba_func,
    maybe_use_numba,
    split_for_numba,
)

from pandas.plotting import boxplot_frame_groupby

if TYPE_CHECKING:
    from pandas.core.internals import Block


NamedAgg = namedtuple("NamedAgg", ["column", "aggfunc"])
# TODO(typing) the return value on this callable should be any *scalar*.
AggScalar = Union[str, Callable[..., Any]]
# TODO: validate types on ScalarResult and move to _typing
# Blocked from using by https://github.com/python/mypy/issues/1484
# See note at _mangle_lambda_list
ScalarResult = typing.TypeVar("ScalarResult")


def generate_property(name: str, klass: Type[FrameOrSeries]):
    """
    Create a property for a GroupBy subclass to dispatch to DataFrame/Series.

    Parameters
    ----------
    name : str
    klass : {DataFrame, Series}

    Returns
    -------
    property
    """

    def prop(self):
        return self._make_wrapper(name)

    parent_method = getattr(klass, name)
    prop.__doc__ = parent_method.__doc__ or ""
    prop.__name__ = name
    return property(prop)


def pin_allowlisted_properties(klass: Type[FrameOrSeries], allowlist: FrozenSet[str]):
    """
    Create GroupBy member defs for DataFrame/Series names in a allowlist.

    Parameters
    ----------
    klass : DataFrame or Series class
        class where members are defined.
    allowlist : frozenset[str]
        Set of names of klass methods to be constructed

    Returns
    -------
    class decorator

    Notes
    -----
    Since we don't want to override methods explicitly defined in the
    base class, any such name is skipped.
    """

    def pinner(cls):
        for name in allowlist:
            if hasattr(cls, name):
                # don't override anything that was explicitly defined
                #  in the base class
                continue

            prop = generate_property(name, klass)
            setattr(cls, name, prop)

        return cls

    return pinner


@pin_allowlisted_properties(Series, base.series_apply_allowlist)
class SeriesGroupBy(GroupBy[Series]):
    _apply_allowlist = base.series_apply_allowlist

    def _iterate_slices(self) -> Iterable[Series]:
        yield self._selected_obj

    @property
    def _selection_name(self):
        """
        since we are a series, we by definition only have
        a single name, but may be the result of a selection or
        the name of our object
        """
        if self._selection is None:
            return self.obj.name
        else:
            return self._selection

    _agg_examples_doc = dedent(
        """
    Examples
    --------
    >>> s = pd.Series([1, 2, 3, 4])

    >>> s
    0    1
    1    2
    2    3
    3    4
    dtype: int64

    >>> s.groupby([1, 1, 2, 2]).min()
    1    1
    2    3
    dtype: int64

    >>> s.groupby([1, 1, 2, 2]).agg('min')
    1    1
    2    3
    dtype: int64

    >>> s.groupby([1, 1, 2, 2]).agg(['min', 'max'])
       min  max
    1    1    2
    2    3    4

    The output column names can be controlled by passing
    the desired column names and aggregations as keyword arguments.

    >>> s.groupby([1, 1, 2, 2]).agg(
    ...     minimum='min',
    ...     maximum='max',
    ... )
       minimum  maximum
    1        1        2
    2        3        4"""
    )

    @Appender(
        _apply_docs["template"].format(
            input="series", examples=_apply_docs["series_examples"]
        )
    )
    def apply(self, func, *args, **kwargs):
        return super().apply(func, *args, **kwargs)

    @doc(
        _agg_template, examples=_agg_examples_doc, klass="Series",
    )
    def aggregate(self, func=None, *args, engine=None, engine_kwargs=None, **kwargs):

        relabeling = func is None
        columns = None
        if relabeling:
            columns, func = validate_func_kwargs(kwargs)
            kwargs = {}

        if isinstance(func, str):
            return getattr(self, func)(*args, **kwargs)

        elif isinstance(func, abc.Iterable):
            # Catch instances of lists / tuples
            # but not the class list / tuple itself.
            func = maybe_mangle_lambdas(func)
            ret = self._aggregate_multiple_funcs(func)
            if relabeling:
                ret.columns = columns
        else:
            cyfunc = self._get_cython_func(func)
            if cyfunc and not args and not kwargs:
                return getattr(self, cyfunc)()

            if self.grouper.nkeys > 1:
                return self._python_agg_general(
                    func, *args, engine=engine, engine_kwargs=engine_kwargs, **kwargs
                )

            try:
                return self._python_agg_general(
                    func, *args, engine=engine, engine_kwargs=engine_kwargs, **kwargs
                )
            except (ValueError, KeyError):
                # Do not catch Numba errors here, we want to raise and not fall back.
                # TODO: KeyError is raised in _python_agg_general,
                #  see see test_groupby.test_basic
                result = self._aggregate_named(func, *args, **kwargs)

            index = Index(sorted(result), name=self.grouper.names[0])
            ret = create_series_with_explicit_dtype(
                result, index=index, dtype_if_empty=object
            )

        if not self.as_index:  # pragma: no cover
            print("Warning, ignoring as_index=True")

        if isinstance(ret, dict):
            from pandas import concat

            ret = concat(ret.values(), axis=1, keys=[key.label for key in ret.keys()])
        return ret

    agg = aggregate

    def _aggregate_multiple_funcs(self, arg):
        if isinstance(arg, dict):

            # show the deprecation, but only if we
            # have not shown a higher level one
            # GH 15931
            if isinstance(self._selected_obj, Series):
                raise SpecificationError("nested renamer is not supported")

            columns = list(arg.keys())
            arg = arg.items()
        elif any(isinstance(x, (tuple, list)) for x in arg):
            arg = [(x, x) if not isinstance(x, (tuple, list)) else x for x in arg]

            # indicated column order
            columns = next(zip(*arg))
        else:
            # list of functions / function names
            columns = []
            for f in arg:
                columns.append(com.get_callable_name(f) or f)

            arg = zip(columns, arg)

        results: Dict[base.OutputKey, Union[Series, DataFrame]] = {}
        for idx, (name, func) in enumerate(arg):
            obj = self

            # reset the cache so that we
            # only include the named selection
            if name in self._selected_obj:
                obj = copy.copy(obj)
                obj._reset_cache()
                obj._selection = name
            results[base.OutputKey(label=name, position=idx)] = obj.aggregate(func)

        if any(isinstance(x, DataFrame) for x in results.values()):
            # let higher level handle
            return results

        output = self._wrap_aggregated_output(results)
        return self.obj._constructor_expanddim(output, columns=columns)

    def _wrap_series_output(
        self, output: Mapping[base.OutputKey, Union[Series, np.ndarray]], index: Index,
    ) -> Union[Series, DataFrame]:
        """
        Wraps the output of a SeriesGroupBy operation into the expected result.

        Parameters
        ----------
        output : Mapping[base.OutputKey, Union[Series, np.ndarray]]
            Data to wrap.
        index : pd.Index
            Index to apply to the output.

        Returns
        -------
        Series or DataFrame

        Notes
        -----
        In the vast majority of cases output and columns will only contain one
        element. The exception is operations that expand dimensions, like ohlc.
        """
        indexed_output = {key.position: val for key, val in output.items()}
        columns = Index(key.label for key in output)

        result: Union[Series, DataFrame]
        if len(output) > 1:
            result = self.obj._constructor_expanddim(indexed_output, index=index)
            result.columns = columns
        elif not columns.empty:
            result = self.obj._constructor(
                indexed_output[0], index=index, name=columns[0]
            )
        else:
            result = self.obj._constructor_expanddim()

        return result

    def _wrap_aggregated_output(
        self, output: Mapping[base.OutputKey, Union[Series, np.ndarray]]
    ) -> Union[Series, DataFrame]:
        """
        Wraps the output of a SeriesGroupBy aggregation into the expected result.

        Parameters
        ----------
        output : Mapping[base.OutputKey, Union[Series, np.ndarray]]
            Data to wrap.

        Returns
        -------
        Series or DataFrame

        Notes
        -----
        In the vast majority of cases output will only contain one element.
        The exception is operations that expand dimensions, like ohlc.
        """
        result = self._wrap_series_output(
            output=output, index=self.grouper.result_index
        )
        return self._reindex_output(result)

    def _wrap_transformed_output(
        self, output: Mapping[base.OutputKey, Union[Series, np.ndarray]]
    ) -> Series:
        """
        Wraps the output of a SeriesGroupBy aggregation into the expected result.

        Parameters
        ----------
        output : dict[base.OutputKey, Union[Series, np.ndarray]]
            Dict with a sole key of 0 and a value of the result values.

        Returns
        -------
        Series

        Notes
        -----
        output should always contain one element. It is specified as a dict
        for consistency with DataFrame methods and _wrap_aggregated_output.
        """
        assert len(output) == 1
        result = self._wrap_series_output(output=output, index=self.obj.index)

        # No transformations increase the ndim of the result
        assert isinstance(result, Series)
        return result

    def _wrap_applied_output(
        self, keys: Index, values: Optional[List[Any]], not_indexed_same: bool = False
    ) -> FrameOrSeriesUnion:
        """
        Wrap the output of SeriesGroupBy.apply into the expected result.

        Parameters
        ----------
        keys : Index
            Keys of groups that Series was grouped by.
        values : Optional[List[Any]]
            Applied output for each group.
        not_indexed_same : bool, default False
            Whether the applied outputs are not indexed the same as the group axes.

        Returns
        -------
        DataFrame or Series
        """
        if len(keys) == 0:
            # GH #6265
            return self.obj._constructor(
                [], name=self._selection_name, index=keys, dtype=np.float64
            )
        assert values is not None

        def _get_index() -> Index:
            if self.grouper.nkeys > 1:
                index = MultiIndex.from_tuples(keys, names=self.grouper.names)
            else:
                index = Index(keys, name=self.grouper.names[0])
            return index

        if isinstance(values[0], dict):
            # GH #823 #24880
            index = _get_index()
            result: FrameOrSeriesUnion = self._reindex_output(
                self.obj._constructor_expanddim(values, index=index)
            )
            # if self.observed is False,
            # keep all-NaN rows created while re-indexing
            result = result.stack(dropna=self.observed)
            result.name = self._selection_name
            return result
        elif isinstance(values[0], (Series, DataFrame)):
            return self._concat_objects(keys, values, not_indexed_same=not_indexed_same)
        else:
            # GH #6265 #24880
            result = self.obj._constructor(
                data=values, index=_get_index(), name=self._selection_name
            )
            return self._reindex_output(result)

    def _aggregate_named(self, func, *args, **kwargs):
        result = {}

        for name, group in self:
            group.name = name
            output = func(group, *args, **kwargs)
            if isinstance(output, (Series, Index, np.ndarray)):
                raise ValueError("Must produce aggregated value")
            result[name] = output

        return result

    @Substitution(klass="Series")
    @Appender(_transform_template)
    def transform(self, func, *args, engine=None, engine_kwargs=None, **kwargs):
        func = self._get_cython_func(func) or func

        if not isinstance(func, str):
            return self._transform_general(
                func, *args, engine=engine, engine_kwargs=engine_kwargs, **kwargs
            )

        elif func not in base.transform_kernel_allowlist:
            msg = f"'{func}' is not a valid function name for transform(name)"
            raise ValueError(msg)
        elif func in base.cythonized_kernels:
            # cythonized transform or canned "agg+broadcast"
            return getattr(self, func)(*args, **kwargs)
        elif func in base.transformation_kernels:
            return getattr(self, func)(*args, **kwargs)

        # If func is a reduction, we need to broadcast the
        # result to the whole group. Compute func result
        # and deal with possible broadcasting below.
        # Temporarily set observed for dealing with categoricals.
        with com.temp_setattr(self, "observed", True):
            result = getattr(self, func)(*args, **kwargs)
        return self._transform_fast(result)

    def _transform_general(
        self, func, *args, engine="cython", engine_kwargs=None, **kwargs
    ):
        """
        Transform with a non-str `func`.
        """

        if maybe_use_numba(engine):
            numba_func, cache_key = generate_numba_func(
                func, engine_kwargs, kwargs, "groupby_transform"
            )

        klass = type(self._selected_obj)

        results = []
        for name, group in self:
            object.__setattr__(group, "name", name)
            if maybe_use_numba(engine):
                values, index = split_for_numba(group)
                res = numba_func(values, index, *args)
                if cache_key not in NUMBA_FUNC_CACHE:
                    NUMBA_FUNC_CACHE[cache_key] = numba_func
            else:
                res = func(group, *args, **kwargs)

            if isinstance(res, (ABCDataFrame, ABCSeries)):
                res = res._values

            indexer = self._get_index(name)
            ser = klass(res, indexer)
            results.append(ser)

        # check for empty "results" to avoid concat ValueError
        if results:
            from pandas.core.reshape.concat import concat

            result = concat(results).sort_index()
        else:
            result = self.obj._constructor(dtype=np.float64)

        # we will only try to coerce the result type if
        # we have a numeric dtype, as these are *always* user-defined funcs
        # the cython take a different path (and casting)
        dtype = self._selected_obj.dtype
        if is_numeric_dtype(dtype):
            result = maybe_downcast_to_dtype(result, dtype)

        result.name = self._selected_obj.name
        result.index = self._selected_obj.index
        return result

    def _transform_fast(self, result) -> Series:
        """
        fast version of transform, only applicable to
        builtin/cythonizable functions
        """
        ids, _, ngroup = self.grouper.group_info
        result = result.reindex(self.grouper.result_index, copy=False)
        out = algorithms.take_1d(result._values, ids)
        return self.obj._constructor(out, index=self.obj.index, name=self.obj.name)

    def filter(self, func, dropna=True, *args, **kwargs):
        """
        Return a copy of a Series excluding elements from groups that
        do not satisfy the boolean criterion specified by func.

        Parameters
        ----------
        func : function
            To apply to each group. Should return True or False.
        dropna : Drop groups that do not pass the filter. True by default;
            if False, groups that evaluate False are filled with NaNs.

        Examples
        --------
        >>> df = pd.DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
        ...                           'foo', 'bar'],
        ...                    'B' : [1, 2, 3, 4, 5, 6],
        ...                    'C' : [2.0, 5., 8., 1., 2., 9.]})
        >>> grouped = df.groupby('A')
        >>> df.groupby('A').B.filter(lambda x: x.mean() > 3.)
        1    2
        3    4
        5    6
        Name: B, dtype: int64

        Returns
        -------
        filtered : Series
        """
        if isinstance(func, str):
            wrapper = lambda x: getattr(x, func)(*args, **kwargs)
        else:
            wrapper = lambda x: func(x, *args, **kwargs)

        # Interpret np.nan as False.
        def true_and_notna(x, *args, **kwargs) -> bool:
            b = wrapper(x, *args, **kwargs)
            return b and notna(b)

        try:
            indices = [
                self._get_index(name) for name, group in self if true_and_notna(group)
            ]
        except (ValueError, TypeError) as err:
            raise TypeError("the filter must return a boolean result") from err

        filtered = self._apply_filter(indices, dropna)
        return filtered

    def nunique(self, dropna: bool = True) -> Series:
        """
        Return number of unique elements in the group.

        Returns
        -------
        Series
            Number of unique values within each group.
        """
        ids, _, _ = self.grouper.group_info

        val = self.obj._values

        codes, _ = algorithms.factorize(val, sort=False)
        sorter = np.lexsort((codes, ids))
        codes = codes[sorter]
        ids = ids[sorter]

        # group boundaries are where group ids change
        # unique observations are where sorted values change
        idx = np.r_[0, 1 + np.nonzero(ids[1:] != ids[:-1])[0]]
        inc = np.r_[1, codes[1:] != codes[:-1]]

        # 1st item of each group is a new unique observation
        mask = codes == -1
        if dropna:
            inc[idx] = 1
            inc[mask] = 0
        else:
            inc[mask & np.r_[False, mask[:-1]]] = 0
            inc[idx] = 1

        out = np.add.reduceat(inc, idx).astype("int64", copy=False)
        if len(ids):
            # NaN/NaT group exists if the head of ids is -1,
            # so remove it from res and exclude its index from idx
            if ids[0] == -1:
                res = out[1:]
                idx = idx[np.flatnonzero(idx)]
            else:
                res = out
        else:
            res = out[1:]
        ri = self.grouper.result_index

        # we might have duplications among the bins
        if len(res) != len(ri):
            res, out = np.zeros(len(ri), dtype=out.dtype), res
            res[ids[idx]] = out

        result = self.obj._constructor(res, index=ri, name=self._selection_name)
        return self._reindex_output(result, fill_value=0)

    @doc(Series.describe)
    def describe(self, **kwargs):
        result = self.apply(lambda x: x.describe(**kwargs))
        if self.axis == 1:
            return result.T
        return result.unstack()

    def value_counts(
        self, normalize=False, sort=True, ascending=False, bins=None, dropna=True
    ):

        from pandas.core.reshape.tile import cut
        from pandas.core.reshape.merge import _get_join_indexers

        if bins is not None and not np.iterable(bins):
            # scalar bins cannot be done at top level
            # in a backward compatible way
            return self.apply(
                Series.value_counts,
                normalize=normalize,
                sort=sort,
                ascending=ascending,
                bins=bins,
            )

        ids, _, _ = self.grouper.group_info
        val = self.obj._values

        # groupby removes null keys from groupings
        mask = ids != -1
        ids, val = ids[mask], val[mask]

        if bins is None:
            lab, lev = algorithms.factorize(val, sort=True)
            llab = lambda lab, inc: lab[inc]
        else:

            # lab is a Categorical with categories an IntervalIndex
            lab = cut(Series(val), bins, include_lowest=True)
            lev = lab.cat.categories
            lab = lev.take(lab.cat.codes)
            llab = lambda lab, inc: lab[inc]._multiindex.codes[-1]

        if is_interval_dtype(lab.dtype):
            # TODO: should we do this inside II?
            sorter = np.lexsort((lab.left, lab.right, ids))
        else:
            sorter = np.lexsort((lab, ids))

        ids, lab = ids[sorter], lab[sorter]

        # group boundaries are where group ids change
        idx = np.r_[0, 1 + np.nonzero(ids[1:] != ids[:-1])[0]]

        # new values are where sorted labels change
        lchanges = llab(lab, slice(1, None)) != llab(lab, slice(None, -1))
        inc = np.r_[True, lchanges]
        inc[idx] = True  # group boundaries are also new values
        out = np.diff(np.nonzero(np.r_[inc, True])[0])  # value counts

        # num. of times each group should be repeated
        rep = partial(np.repeat, repeats=np.add.reduceat(inc, idx))

        # multi-index components
        codes = self.grouper.reconstructed_codes
        codes = [rep(level_codes) for level_codes in codes] + [llab(lab, inc)]
        levels = [ping.group_index for ping in self.grouper.groupings] + [lev]
        names = self.grouper.names + [self._selection_name]

        if dropna:
            mask = codes[-1] != -1
            if mask.all():
                dropna = False
            else:
                out, codes = out[mask], [level_codes[mask] for level_codes in codes]

        if normalize:
            out = out.astype("float")
            d = np.diff(np.r_[idx, len(ids)])
            if dropna:
                m = ids[lab == -1]
                np.add.at(d, m, -1)
                acc = rep(d)[mask]
            else:
                acc = rep(d)
            out /= acc

        if sort and bins is None:
            cat = ids[inc][mask] if dropna else ids[inc]
            sorter = np.lexsort((out if ascending else -out, cat))
            out, codes[-1] = out[sorter], codes[-1][sorter]

        if bins is None:
            mi = MultiIndex(
                levels=levels, codes=codes, names=names, verify_integrity=False
            )

            if is_integer_dtype(out):
                out = ensure_int64(out)
            return self.obj._constructor(out, index=mi, name=self._selection_name)

        # for compat. with libgroupby.value_counts need to ensure every
        # bin is present at every index level, null filled with zeros
        diff = np.zeros(len(out), dtype="bool")
        for level_codes in codes[:-1]:
            diff |= np.r_[True, level_codes[1:] != level_codes[:-1]]

        ncat, nbin = diff.sum(), len(levels[-1])

        left = [np.repeat(np.arange(ncat), nbin), np.tile(np.arange(nbin), ncat)]

        right = [diff.cumsum() - 1, codes[-1]]

        _, idx = _get_join_indexers(left, right, sort=False, how="left")
        out = np.where(idx != -1, out[idx], 0)

        if sort:
            sorter = np.lexsort((out if ascending else -out, left[0]))
            out, left[-1] = out[sorter], left[-1][sorter]

        # build the multi-index w/ full levels
        def build_codes(lev_codes: np.ndarray) -> np.ndarray:
            return np.repeat(lev_codes[diff], nbin)

        codes = [build_codes(lev_codes) for lev_codes in codes[:-1]]
        codes.append(left[-1])

        mi = MultiIndex(levels=levels, codes=codes, names=names, verify_integrity=False)

        if is_integer_dtype(out):
            out = ensure_int64(out)
        return self.obj._constructor(out, index=mi, name=self._selection_name)

    def count(self) -> Series:
        """
        Compute count of group, excluding missing values.

        Returns
        -------
        Series
            Count of values within each group.
        """
        ids, _, ngroups = self.grouper.group_info
        val = self.obj._values

        mask = (ids != -1) & ~isna(val)
        ids = ensure_platform_int(ids)
        minlength = ngroups or 0
        out = np.bincount(ids[mask], minlength=minlength)

        result = self.obj._constructor(
            out,
            index=self.grouper.result_index,
            name=self._selection_name,
            dtype="int64",
        )
        return self._reindex_output(result, fill_value=0)

    def _apply_to_column_groupbys(self, func):
        """ return a pass thru """
        return func(self)

    def pct_change(self, periods=1, fill_method="pad", limit=None, freq=None):
        """Calculate pct_change of each value to previous entry in group"""
        # TODO: Remove this conditional when #23918 is fixed
        if freq:
            return self.apply(
                lambda x: x.pct_change(
                    periods=periods, fill_method=fill_method, limit=limit, freq=freq
                )
            )
        if fill_method is None:  # GH30463
            fill_method = "pad"
            limit = 0
        filled = getattr(self, fill_method)(limit=limit)
        fill_grp = filled.groupby(self.grouper.codes)
        shifted = fill_grp.shift(periods=periods, freq=freq)

        return (filled / shifted) - 1


@pin_allowlisted_properties(DataFrame, base.dataframe_apply_allowlist)
class DataFrameGroupBy(GroupBy[DataFrame]):

    _apply_allowlist = base.dataframe_apply_allowlist

    _agg_examples_doc = dedent(
        """
    Examples
    --------
    >>> df = pd.DataFrame(
    ...     {
    ...         "A": [1, 1, 2, 2],
    ...         "B": [1, 2, 3, 4],
    ...         "C": [0.362838, 0.227877, 1.267767, -0.562860],
    ...     }
    ... )

    >>> df
       A  B         C
    0  1  1  0.362838
    1  1  2  0.227877
    2  2  3  1.267767
    3  2  4 -0.562860

    The aggregation is for each column.

    >>> df.groupby('A').agg('min')
       B         C
    A
    1  1  0.227877
    2  3 -0.562860

    Multiple aggregations

    >>> df.groupby('A').agg(['min', 'max'])
        B             C
      min max       min       max
    A
    1   1   2  0.227877  0.362838
    2   3   4 -0.562860  1.267767

    Select a column for aggregation

    >>> df.groupby('A').B.agg(['min', 'max'])
       min  max
    A
    1    1    2
    2    3    4

    Different aggregations per column

    >>> df.groupby('A').agg({'B': ['min', 'max'], 'C': 'sum'})
        B             C
      min max       sum
    A
    1   1   2  0.590715
    2   3   4  0.704907

    To control the output names with different aggregations per column,
    pandas supports "named aggregation"

    >>> df.groupby("A").agg(
    ...     b_min=pd.NamedAgg(column="B", aggfunc="min"),
    ...     c_sum=pd.NamedAgg(column="C", aggfunc="sum"))
       b_min     c_sum
    A
    1      1  0.590715
    2      3  0.704907

    - The keywords are the *output* column names
    - The values are tuples whose first element is the column to select
      and the second element is the aggregation to apply to that column.
      Pandas provides the ``pandas.NamedAgg`` namedtuple with the fields
      ``['column', 'aggfunc']`` to make it clearer what the arguments are.
      As usual, the aggregation can be a callable or a string alias.

    See :ref:`groupby.aggregate.named` for more."""
    )

    @doc(
        _agg_template, examples=_agg_examples_doc, klass="DataFrame",
    )
    def aggregate(self, func=None, *args, engine=None, engine_kwargs=None, **kwargs):

        relabeling, func, columns, order = reconstruct_func(func, **kwargs)

        if maybe_use_numba(engine):
            return self._python_agg_general(
                func, *args, engine=engine, engine_kwargs=engine_kwargs, **kwargs
            )

        result, how = self._aggregate(func, *args, **kwargs)
        if how is None:
            return result

        if result is None:

            # grouper specific aggregations
            if self.grouper.nkeys > 1:
                return self._python_agg_general(func, *args, **kwargs)
            elif args or kwargs:
                result = self._aggregate_frame(func, *args, **kwargs)

            elif self.axis == 1:
                # _aggregate_multiple_funcs does not allow self.axis == 1
                result = self._aggregate_frame(func)

            else:

                # try to treat as if we are passing a list
                try:
                    result = self._aggregate_multiple_funcs([func], _axis=self.axis)

                    # select everything except for the last level, which is the one
                    # containing the name of the function(s), see GH 32040
                    result.columns = result.columns.rename(
                        [self._selected_obj.columns.name] * result.columns.nlevels
                    ).droplevel(-1)

                except ValueError as err:
                    if "no results" not in str(err):
                        # raised directly by _aggregate_multiple_funcs
                        raise
                    result = self._aggregate_frame(func)
                except AttributeError:
                    # catch exception from line 969
                    # (Series does not have attribute "columns"), see GH 35246
                    result = self._aggregate_frame(func)

        if relabeling:

            # used reordered index of columns
            result = result.iloc[:, order]
            result.columns = columns

        if not self.as_index:
            self._insert_inaxis_grouper_inplace(result)
            result.index = np.arange(len(result))

        return result._convert(datetime=True)

    agg = aggregate

    def _iterate_slices(self) -> Iterable[Series]:
        obj = self._selected_obj
        if self.axis == 1:
            obj = obj.T

        if isinstance(obj, Series) and obj.name not in self.exclusions:
            # Occurs when doing DataFrameGroupBy(...)["X"]
            yield obj
        else:
            for label, values in obj.items():
                if label in self.exclusions:
                    continue

                yield values

    def _cython_agg_general(
        self, how: str, alt=None, numeric_only: bool = True, min_count: int = -1
    ) -> DataFrame:
        agg_blocks, agg_items = self._cython_agg_blocks(
            how, alt=alt, numeric_only=numeric_only, min_count=min_count
        )
        return self._wrap_agged_blocks(agg_blocks, items=agg_items)

    def _cython_agg_blocks(
        self, how: str, alt=None, numeric_only: bool = True, min_count: int = -1
    ) -> "Tuple[List[Block], Index]":
        # TODO: the actual managing of mgr_locs is a PITA
        # here, it should happen via BlockManager.combine

        data: BlockManager = self._get_data_to_aggregate()

        if numeric_only:
            data = data.get_numeric_data(copy=False)

        agg_blocks: List[Block] = []
        new_items: List[np.ndarray] = []
        deleted_items: List[np.ndarray] = []
        # Some object-dtype blocks might be split into List[Block[T], Block[U]]
        split_items: List[np.ndarray] = []
        split_frames: List[DataFrame] = []

        no_result = object()
        for block in data.blocks:
            # Avoid inheriting result from earlier in the loop
            result = no_result
            locs = block.mgr_locs.as_array
            try:
                result, _ = self.grouper.aggregate(
                    block.values, how, axis=1, min_count=min_count
                )
            except NotImplementedError:
                # generally if we have numeric_only=False
                # and non-applicable functions
                # try to python agg

                if alt is None:
                    # we cannot perform the operation
                    # in an alternate way, exclude the block
                    assert how == "ohlc"
                    deleted_items.append(locs)
                    continue

                # call our grouper again with only this block
                obj = self.obj[data.items[locs]]
                if obj.shape[1] == 1:
                    # Avoid call to self.values that can occur in DataFrame
                    #  reductions; see GH#28949
                    obj = obj.iloc[:, 0]

                # Create SeriesGroupBy with observed=True so that it does
                # not try to add missing categories if grouping over multiple
                # Categoricals. This will done by later self._reindex_output()
                # Doing it here creates an error. See GH#34951
                s = get_groupby(obj, self.grouper, observed=True)
                try:
                    result = s.aggregate(lambda x: alt(x, axis=self.axis))
                except TypeError:
                    # we may have an exception in trying to aggregate
                    # continue and exclude the block
                    deleted_items.append(locs)
                    continue
                else:
                    result = cast(DataFrame, result)
                    # unwrap DataFrame to get array
                    if len(result._mgr.blocks) != 1:
                        # We've split an object block! Everything we've assumed
                        # about a single block input returning a single block output
                        # is a lie. To keep the code-path for the typical non-split case
                        # clean, we choose to clean up this mess later on.
                        split_items.append(locs)
                        split_frames.append(result)
                        continue

                    assert len(result._mgr.blocks) == 1
                    result = result._mgr.blocks[0].values
                    if isinstance(result, np.ndarray) and result.ndim == 1:
                        result = result.reshape(1, -1)

            assert not isinstance(result, DataFrame)

            if result is not no_result:
                # see if we can cast the block to the desired dtype
                # this may not be the original dtype
                dtype = maybe_cast_result_dtype(block.dtype, how)
                result = maybe_downcast_numeric(result, dtype)

                if block.is_extension and isinstance(result, np.ndarray):
                    # e.g. block.values was an IntegerArray
                    # (1, N) case can occur if block.values was Categorical
                    #  and result is ndarray[object]
                    # TODO(EA2D): special casing not needed with 2D EAs
                    assert result.ndim == 1 or result.shape[0] == 1
                    try:
                        # Cast back if feasible
                        result = type(block.values)._from_sequence(
                            result.ravel(), dtype=block.values.dtype
                        )
                    except (ValueError, TypeError):
                        # reshape to be valid for non-Extension Block
                        result = result.reshape(1, -1)

                agg_block: Block = block.make_block(result)

            new_items.append(locs)
            agg_blocks.append(agg_block)

        if not (agg_blocks or split_frames):
            raise DataError("No numeric types to aggregate")

        if split_items:
            # Clean up the mess left over from split blocks.
            for locs, result in zip(split_items, split_frames):
                assert len(locs) == result.shape[1]
                for i, loc in enumerate(locs):
                    new_items.append(np.array([loc], dtype=locs.dtype))
                    agg_blocks.append(result.iloc[:, [i]]._mgr.blocks[0])

        # reset the locs in the blocks to correspond to our
        # current ordering
        indexer = np.concatenate(new_items)
        agg_items = data.items.take(np.sort(indexer))

        if deleted_items:

            # we need to adjust the indexer to account for the
            # items we have removed
            # really should be done in internals :<

            deleted = np.concatenate(deleted_items)
            ai = np.arange(len(data))
            mask = np.zeros(len(data))
            mask[deleted] = 1
            indexer = (ai - mask.cumsum())[indexer]

        offset = 0
        for blk in agg_blocks:
            loc = len(blk.mgr_locs)
            blk.mgr_locs = indexer[offset : (offset + loc)]
            offset += loc

        return agg_blocks, agg_items

    def _aggregate_frame(self, func, *args, **kwargs) -> DataFrame:
        if self.grouper.nkeys != 1:
            raise AssertionError("Number of keys must be 1")

        axis = self.axis
        obj = self._obj_with_exclusions

        result: Dict[Union[int, str], Union[NDFrame, np.ndarray]] = {}
        if axis != obj._info_axis_number:
            for name, data in self:
                fres = func(data, *args, **kwargs)
                result[name] = fres
        else:
            for name in self.indices:
                data = self.get_group(name, obj=obj)
                fres = func(data, *args, **kwargs)
                result[name] = fres

        return self._wrap_frame_output(result, obj)

    def _aggregate_item_by_item(self, func, *args, **kwargs) -> DataFrame:
        # only for axis==0

        obj = self._obj_with_exclusions
        result: Dict[Union[int, str], NDFrame] = {}
        cannot_agg = []
        for item in obj:
            data = obj[item]
            colg = SeriesGroupBy(data, selection=item, grouper=self.grouper)

            cast = self._transform_should_cast(func)
            try:
                result[item] = colg.aggregate(func, *args, **kwargs)

            except ValueError as err:
                if "Must produce aggregated value" in str(err):
                    # raised in _aggregate_named, handle at higher level
                    #  see test_apply_with_mutated_index
                    raise
                # otherwise we get here from an AttributeError in _make_wrapper
                cannot_agg.append(item)
                continue

            else:
                if cast:
                    result[item] = maybe_cast_result(result[item], data)

        result_columns = obj.columns
        if cannot_agg:
            result_columns = result_columns.drop(cannot_agg)

        return self.obj._constructor(result, columns=result_columns)

    def _wrap_applied_output(self, keys, values, not_indexed_same=False):
        if len(keys) == 0:
            return self.obj._constructor(index=keys)

        key_names = self.grouper.names

        # GH12824
        first_not_none = next(com.not_none(*values), None)

        if first_not_none is None:
            # GH9684. If all values are None, then this will throw an error.
            # We'd prefer it return an empty dataframe.
            return self.obj._constructor()
        elif isinstance(first_not_none, DataFrame):
            return self._concat_objects(keys, values, not_indexed_same=not_indexed_same)
        else:
            if len(self.grouper.groupings) > 1:
                key_index = self.grouper.result_index

            else:
                ping = self.grouper.groupings[0]
                if len(keys) == ping.ngroups:
                    key_index = ping.group_index
                    key_index.name = key_names[0]

                    key_lookup = Index(keys)
                    indexer = key_lookup.get_indexer(key_index)

                    # reorder the values
                    values = [values[i] for i in indexer]

                    # update due to the potential reorder
                    first_not_none = next(com.not_none(*values), None)
                else:

                    key_index = Index(keys, name=key_names[0])

                # don't use the key indexer
                if not self.as_index:
                    key_index = None

            # make Nones an empty object
            if first_not_none is None:
                return self.obj._constructor()
            elif isinstance(first_not_none, NDFrame):

                # this is to silence a DeprecationWarning
                # TODO: Remove when default dtype of empty Series is object
                kwargs = first_not_none._construct_axes_dict()
                if isinstance(first_not_none, Series):
                    backup = create_series_with_explicit_dtype(
                        **kwargs, dtype_if_empty=object
                    )
                else:
                    backup = first_not_none._constructor(**kwargs)

                values = [x if (x is not None) else backup for x in values]

            v = values[0]

            if isinstance(v, (np.ndarray, Index, Series)) or not self.as_index:
                if isinstance(v, Series):
                    applied_index = self._selected_obj._get_axis(self.axis)
                    all_indexed_same = all_indexes_same([x.index for x in values])
                    singular_series = len(values) == 1 and applied_index.nlevels == 1

                    # GH3596
                    # provide a reduction (Frame -> Series) if groups are
                    # unique
                    if self.squeeze:
                        # assign the name to this series
                        if singular_series:
                            values[0].name = keys[0]

                            # GH2893
                            # we have series in the values array, we want to
                            # produce a series:
                            # if any of the sub-series are not indexed the same
                            # OR we don't have a multi-index and we have only a
                            # single values
                            return self._concat_objects(
                                keys, values, not_indexed_same=not_indexed_same
                            )

                        # still a series
                        # path added as of GH 5545
                        elif all_indexed_same:
                            from pandas.core.reshape.concat import concat

                            return concat(values)

                    if not all_indexed_same:
                        # GH 8467
                        return self._concat_objects(keys, values, not_indexed_same=True)

                if self.axis == 0 and isinstance(v, ABCSeries):
                    # GH6124 if the list of Series have a consistent name,
                    # then propagate that name to the result.
                    index = v.index.copy()
                    if index.name is None:
                        # Only propagate the series name to the result
                        # if all series have a consistent name.  If the
                        # series do not have a consistent name, do
                        # nothing.
                        names = {v.name for v in values}
                        if len(names) == 1:
                            index.name = list(names)[0]

                    # normally use vstack as its faster than concat
                    # and if we have mi-columns
                    if (
                        isinstance(v.index, MultiIndex)
                        or key_index is None
                        or isinstance(key_index, MultiIndex)
                    ):
                        stacked_values = np.vstack([np.asarray(v) for v in values])
                        result = self.obj._constructor(
                            stacked_values, index=key_index, columns=index
                        )
                    else:
                        # GH5788 instead of stacking; concat gets the
                        # dtypes correct
                        from pandas.core.reshape.concat import concat

                        result = concat(
                            values,
                            keys=key_index,
                            names=key_index.names,
                            axis=self.axis,
                        ).unstack()
                        result.columns = index
                elif isinstance(v, ABCSeries):
                    stacked_values = np.vstack([np.asarray(v) for v in values])
                    result = self.obj._constructor(
                        stacked_values.T, index=v.index, columns=key_index
                    )
                elif not self.as_index:
                    # We add grouping column below, so create a frame here
                    result = DataFrame(
                        values, index=key_index, columns=[self._selection]
                    )
                else:
                    # GH#1738: values is list of arrays of unequal lengths
                    #  fall through to the outer else clause
                    # TODO: sure this is right?  we used to do this
                    #  after raising AttributeError above
                    return self.obj._constructor_sliced(
                        values, index=key_index, name=self._selection_name
                    )

                # if we have date/time like in the original, then coerce dates
                # as we are stacking can easily have object dtypes here
                so = self._selected_obj
                if so.ndim == 2 and so.dtypes.apply(needs_i8_conversion).any():
                    result = _recast_datetimelike_result(result)
                else:
                    result = result._convert(datetime=True)

                if not self.as_index:
                    self._insert_inaxis_grouper_inplace(result)

                return self._reindex_output(result)

            # values are not series or array-like but scalars
            else:
                # self._selection_name not passed through to Series as the
                # result should not take the name of original selection
                # of columns
                return self.obj._constructor_sliced(values, index=key_index)

    def _transform_general(
        self, func, *args, engine="cython", engine_kwargs=None, **kwargs
    ):
        from pandas.core.reshape.concat import concat

        applied = []
        obj = self._obj_with_exclusions
        gen = self.grouper.get_iterator(obj, axis=self.axis)
        if maybe_use_numba(engine):
            numba_func, cache_key = generate_numba_func(
                func, engine_kwargs, kwargs, "groupby_transform"
            )
        else:
            fast_path, slow_path = self._define_paths(func, *args, **kwargs)

        for name, group in gen:
            object.__setattr__(group, "name", name)

            if maybe_use_numba(engine):
                values, index = split_for_numba(group)
                res = numba_func(values, index, *args)
                if cache_key not in NUMBA_FUNC_CACHE:
                    NUMBA_FUNC_CACHE[cache_key] = numba_func
                # Return the result as a DataFrame for concatenation later
                res = self.obj._constructor(
                    res, index=group.index, columns=group.columns
                )
            else:
                # Try slow path and fast path.
                try:
                    path, res = self._choose_path(fast_path, slow_path, group)
                except TypeError:
                    return self._transform_item_by_item(obj, fast_path)
                except ValueError as err:
                    msg = "transform must return a scalar value for each group"
                    raise ValueError(msg) from err

            if isinstance(res, Series):

                # we need to broadcast across the
                # other dimension; this will preserve dtypes
                # GH14457
                if not np.prod(group.shape):
                    continue
                elif res.index.is_(obj.index):
                    r = concat([res] * len(group.columns), axis=1)
                    r.columns = group.columns
                    r.index = group.index
                else:
                    r = self.obj._constructor(
                        np.concatenate([res.values] * len(group.index)).reshape(
                            group.shape
                        ),
                        columns=group.columns,
                        index=group.index,
                    )

                applied.append(r)
            else:
                applied.append(res)

        concat_index = obj.columns if self.axis == 0 else obj.index
        other_axis = 1 if self.axis == 0 else 0  # switches between 0 & 1
        concatenated = concat(applied, axis=self.axis, verify_integrity=False)
        concatenated = concatenated.reindex(concat_index, axis=other_axis, copy=False)
        return self._set_result_index_ordered(concatenated)

    @Substitution(klass="DataFrame")
    @Appender(_transform_template)
    def transform(self, func, *args, engine=None, engine_kwargs=None, **kwargs):

        # optimized transforms
        func = self._get_cython_func(func) or func

        if not isinstance(func, str):
            return self._transform_general(
                func, *args, engine=engine, engine_kwargs=engine_kwargs, **kwargs
            )

        elif func not in base.transform_kernel_allowlist:
            msg = f"'{func}' is not a valid function name for transform(name)"
            raise ValueError(msg)
        elif func in base.cythonized_kernels:
            # cythonized transformation or canned "reduction+broadcast"
            return getattr(self, func)(*args, **kwargs)
        elif func in base.transformation_kernels:
            return getattr(self, func)(*args, **kwargs)

        # GH 30918
        # Use _transform_fast only when we know func is an aggregation
        if func in base.reduction_kernels:
            # If func is a reduction, we need to broadcast the
            # result to the whole group. Compute func result
            # and deal with possible broadcasting below.
            # Temporarily set observed for dealing with categoricals.
            with com.temp_setattr(self, "observed", True):
                result = getattr(self, func)(*args, **kwargs)

            if isinstance(result, DataFrame) and result.columns.equals(
                self._obj_with_exclusions.columns
            ):
                return self._transform_fast(result)

        return self._transform_general(
            func, engine=engine, engine_kwargs=engine_kwargs, *args, **kwargs
        )

    def _transform_fast(self, result: DataFrame) -> DataFrame:
        """
        Fast transform path for aggregations
        """
        obj = self._obj_with_exclusions

        # for each col, reshape to to size of original frame
        # by take operation
        ids, _, ngroup = self.grouper.group_info
        result = result.reindex(self.grouper.result_index, copy=False)
        output = []
        for i, _ in enumerate(result.columns):
            output.append(algorithms.take_1d(result.iloc[:, i].values, ids))

        return self.obj._constructor._from_arrays(
            output, columns=result.columns, index=obj.index
        )

    def _define_paths(self, func, *args, **kwargs):
        if isinstance(func, str):
            fast_path = lambda group: getattr(group, func)(*args, **kwargs)
            slow_path = lambda group: group.apply(
                lambda x: getattr(x, func)(*args, **kwargs), axis=self.axis
            )
        else:
            fast_path = lambda group: func(group, *args, **kwargs)
            slow_path = lambda group: group.apply(
                lambda x: func(x, *args, **kwargs), axis=self.axis
            )
        return fast_path, slow_path

    def _choose_path(self, fast_path: Callable, slow_path: Callable, group: DataFrame):
        path = slow_path
        res = slow_path(group)

        # if we make it here, test if we can use the fast path
        try:
            res_fast = fast_path(group)
        except AssertionError:
            raise
        except Exception:
            # GH#29631 For user-defined function, we cant predict what may be
            #  raised; see test_transform.test_transform_fastpath_raises
            return path, res

        # verify fast path does not change columns (and names), otherwise
        # its results cannot be joined with those of the slow path
        if not isinstance(res_fast, DataFrame):
            return path, res

        if not res_fast.columns.equals(group.columns):
            return path, res

        if res_fast.equals(res):
            path = fast_path

        return path, res

    def _transform_item_by_item(self, obj: DataFrame, wrapper) -> DataFrame:
        # iterate through columns
        output = {}
        inds = []
        for i, col in enumerate(obj):
            try:
                output[col] = self[col].transform(wrapper)
            except TypeError:
                # e.g. trying to call nanmean with string values
                pass
            else:
                inds.append(i)

        if len(output) == 0:
            raise TypeError("Transform function invalid for data types")

        columns = obj.columns
        if len(output) < len(obj.columns):
            columns = columns.take(inds)

        return self.obj._constructor(output, index=obj.index, columns=columns)

    def filter(self, func, dropna=True, *args, **kwargs):
        """
        Return a copy of a DataFrame excluding filtered elements.

        Elements from groups are filtered if they do not satisfy the
        boolean criterion specified by func.

        Parameters
        ----------
        func : function
            Function to apply to each subframe. Should return True or False.
        dropna : Drop groups that do not pass the filter. True by default;
            If False, groups that evaluate False are filled with NaNs.

        Returns
        -------
        filtered : DataFrame

        Notes
        -----
        Each subframe is endowed the attribute 'name' in case you need to know
        which group you are working on.

        Examples
        --------
        >>> df = pd.DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
        ...                           'foo', 'bar'],
        ...                    'B' : [1, 2, 3, 4, 5, 6],
        ...                    'C' : [2.0, 5., 8., 1., 2., 9.]})
        >>> grouped = df.groupby('A')
        >>> grouped.filter(lambda x: x['B'].mean() > 3.)
             A  B    C
        1  bar  2  5.0
        3  bar  4  1.0
        5  bar  6  9.0
        """
        indices = []

        obj = self._selected_obj
        gen = self.grouper.get_iterator(obj, axis=self.axis)

        for name, group in gen:
            object.__setattr__(group, "name", name)

            res = func(group, *args, **kwargs)

            try:
                res = res.squeeze()
            except AttributeError:  # allow e.g., scalars and frames to pass
                pass

            # interpret the result of the filter
            if is_bool(res) or (is_scalar(res) and isna(res)):
                if res and notna(res):
                    indices.append(self._get_index(name))
            else:
                # non scalars aren't allowed
                raise TypeError(
                    f"filter function returned a {type(res).__name__}, "
                    "but expected a scalar bool"
                )

        return self._apply_filter(indices, dropna)

    def __getitem__(self, key):
        # per GH 23566
        if isinstance(key, tuple) and len(key) > 1:
            # if len == 1, then it becomes a SeriesGroupBy and this is actually
            # valid syntax, so don't raise warning
            warnings.warn(
                "Indexing with multiple keys (implicitly converted to a tuple "
                "of keys) will be deprecated, use a list instead.",
                FutureWarning,
                stacklevel=2,
            )
        return super().__getitem__(key)

    def _gotitem(self, key, ndim: int, subset=None):
        """
        sub-classes to define
        return a sliced object

        Parameters
        ----------
        key : string / list of selections
        ndim : 1,2
            requested ndim of result
        subset : object, default None
            subset to act on
        """
        if ndim == 2:
            if subset is None:
                subset = self.obj
            return DataFrameGroupBy(
                subset,
                self.grouper,
                selection=key,
                grouper=self.grouper,
                exclusions=self.exclusions,
                as_index=self.as_index,
                observed=self.observed,
            )
        elif ndim == 1:
            if subset is None:
                subset = self.obj[key]
            return SeriesGroupBy(
                subset, selection=key, grouper=self.grouper, observed=self.observed
            )

        raise AssertionError("invalid ndim for _gotitem")

    def _wrap_frame_output(self, result, obj) -> DataFrame:
        result_index = self.grouper.levels[0]

        if self.axis == 0:
            return self.obj._constructor(
                result, index=obj.columns, columns=result_index
            ).T
        else:
            return self.obj._constructor(result, index=obj.index, columns=result_index)

    def _get_data_to_aggregate(self) -> BlockManager:
        obj = self._obj_with_exclusions
        if self.axis == 1:
            return obj.T._mgr
        else:
            return obj._mgr

    def _insert_inaxis_grouper_inplace(self, result):
        # zip in reverse so we can always insert at loc 0
        izip = zip(
            *map(
                reversed,
                (
                    self.grouper.names,
                    self.grouper.get_group_levels(),
                    [grp.in_axis for grp in self.grouper.groupings],
                ),
            )
        )
        columns = result.columns
        for name, lev, in_axis in izip:
            # GH #28549
            # When using .apply(-), name will be in columns already
            if in_axis and name not in columns:
                result.insert(0, name, lev)

    def _wrap_aggregated_output(
        self, output: Mapping[base.OutputKey, Union[Series, np.ndarray]]
    ) -> DataFrame:
        """
        Wraps the output of DataFrameGroupBy aggregations into the expected result.

        Parameters
        ----------
        output : Mapping[base.OutputKey, Union[Series, np.ndarray]]
           Data to wrap.

        Returns
        -------
        DataFrame
        """
        indexed_output = {key.position: val for key, val in output.items()}
        name = self._obj_with_exclusions._get_axis(1 - self.axis).name
        columns = Index([key.label for key in output], name=name)

        result = self.obj._constructor(indexed_output)
        result.columns = columns

        if not self.as_index:
            self._insert_inaxis_grouper_inplace(result)
            result = result._consolidate()
        else:
            index = self.grouper.result_index
            result.index = index

        if self.axis == 1:
            result = result.T

        return self._reindex_output(result)

    def _wrap_transformed_output(
        self, output: Mapping[base.OutputKey, Union[Series, np.ndarray]]
    ) -> DataFrame:
        """
        Wraps the output of DataFrameGroupBy transformations into the expected result.

        Parameters
        ----------
        output : Mapping[base.OutputKey, Union[Series, np.ndarray]]
            Data to wrap.

        Returns
        -------
        DataFrame
        """
        indexed_output = {key.position: val for key, val in output.items()}
        columns = Index(key.label for key in output)

        result = self.obj._constructor(indexed_output)
        result.columns = columns
        result.index = self.obj.index

        return result

    def _wrap_agged_blocks(self, blocks: "Sequence[Block]", items: Index) -> DataFrame:
        if not self.as_index:
            index = np.arange(blocks[0].values.shape[-1])
            mgr = BlockManager(blocks, axes=[items, index])
            result = self.obj._constructor(mgr)

            self._insert_inaxis_grouper_inplace(result)
            result = result._consolidate()
        else:
            index = self.grouper.result_index
            mgr = BlockManager(blocks, axes=[items, index])
            result = self.obj._constructor(mgr)

        if self.axis == 1:
            result = result.T

        return self._reindex_output(result)._convert(datetime=True)

    def _iterate_column_groupbys(self):
        for i, colname in enumerate(self._selected_obj.columns):
            yield colname, SeriesGroupBy(
                self._selected_obj.iloc[:, i],
                selection=colname,
                grouper=self.grouper,
                exclusions=self.exclusions,
            )

    def _apply_to_column_groupbys(self, func):
        from pandas.core.reshape.concat import concat

        return concat(
            (func(col_groupby) for _, col_groupby in self._iterate_column_groupbys()),
            keys=self._selected_obj.columns,
            axis=1,
        )

    def count(self):
        """
        Compute count of group, excluding missing values.

        Returns
        -------
        DataFrame
            Count of values within each group.
        """
        data = self._get_data_to_aggregate()
        ids, _, ngroups = self.grouper.group_info
        mask = ids != -1

        # TODO(2DEA): reshape would not be necessary with 2D EAs
        vals = ((mask & ~isna(blk.values).reshape(blk.shape)) for blk in data.blocks)
        locs = (blk.mgr_locs for blk in data.blocks)

        counted = (
            lib.count_level_2d(x, labels=ids, max_bin=ngroups, axis=1) for x in vals
        )
        blocks = [make_block(val, placement=loc) for val, loc in zip(counted, locs)]

        return self._wrap_agged_blocks(blocks, items=data.items)

    def nunique(self, dropna: bool = True):
        """
        Return DataFrame with counts of unique elements in each position.

        Parameters
        ----------
        dropna : bool, default True
            Don't include NaN in the counts.

        Returns
        -------
        nunique: DataFrame

        Examples
        --------
        >>> df = pd.DataFrame({'id': ['spam', 'egg', 'egg', 'spam',
        ...                           'ham', 'ham'],
        ...                    'value1': [1, 5, 5, 2, 5, 5],
        ...                    'value2': list('abbaxy')})
        >>> df
             id  value1 value2
        0  spam       1      a
        1   egg       5      b
        2   egg       5      b
        3  spam       2      a
        4   ham       5      x
        5   ham       5      y

        >>> df.groupby('id').nunique()
              value1  value2
        id
        egg        1       1
        ham        1       2
        spam       2       1

        Check for rows with the same id but conflicting values:

        >>> df.groupby('id').filter(lambda g: (g.nunique() > 1).any())
             id  value1 value2
        0  spam       1      a
        3  spam       2      a
        4   ham       5      x
        5   ham       5      y
        """
        from pandas.core.reshape.concat import concat

        # TODO: this is duplicative of how GroupBy naturally works
        # Try to consolidate with normal wrapping functions

        obj = self._obj_with_exclusions
        axis_number = obj._get_axis_number(self.axis)
        other_axis = int(not axis_number)
        if axis_number == 0:
            iter_func = obj.items
        else:
            iter_func = obj.iterrows

        results = concat(
            [
                SeriesGroupBy(content, selection=label, grouper=self.grouper).nunique(
                    dropna
                )
                for label, content in iter_func()
            ],
            axis=1,
        )

        if axis_number == 1:
            results = results.T

        results._get_axis(other_axis).names = obj._get_axis(other_axis).names

        if not self.as_index:
            results.index = ibase.default_index(len(results))
            self._insert_inaxis_grouper_inplace(results)
        return results

    boxplot = boxplot_frame_groupby


def _recast_datetimelike_result(result: DataFrame) -> DataFrame:
    """
    If we have date/time like in the original, then coerce dates
    as we are stacking can easily have object dtypes here.

    Parameters
    ----------
    result : DataFrame

    Returns
    -------
    DataFrame

    Notes
    -----
    - Assumes Groupby._selected_obj has ndim==2 and at least one
    datetimelike column
    """
    result = result.copy()

    obj_cols = [
        idx
        for idx in range(len(result.columns))
        if is_object_dtype(result.dtypes.iloc[idx])
    ]

    # See GH#26285
    for n in obj_cols:
        converted = maybe_convert_objects(
            result.iloc[:, n].values, convert_numeric=False
        )

        result.iloc[:, n] = converted
    return result
