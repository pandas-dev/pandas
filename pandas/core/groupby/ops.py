"""
Provide classes to perform the groupby aggregate operations.

These are not exposed to the user and provide implementations of the grouping
operations, primarily in cython. These classes (BaseGrouper and BinGrouper)
are contained *in* the SeriesGroupBy and DataFrameGroupBy objects.
"""

import collections
from typing import (
    Dict,
    Generic,
    Hashable,
    Iterator,
    List,
    Optional,
    Sequence,
    Tuple,
    Type,
)

import numpy as np

from pandas._libs import NaT, iNaT, lib
import pandas._libs.groupby as libgroupby
import pandas._libs.reduction as libreduction
from pandas._typing import ArrayLike, F, FrameOrSeries, Label, Shape, final
from pandas.errors import AbstractMethodError
from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.cast import (
    maybe_cast_result,
    maybe_cast_result_dtype,
    maybe_downcast_to_dtype,
)
from pandas.core.dtypes.common import (
    ensure_float,
    ensure_float64,
    ensure_int64,
    ensure_int_or_float,
    ensure_platform_int,
    is_bool_dtype,
    is_categorical_dtype,
    is_complex_dtype,
    is_datetime64_any_dtype,
    is_datetime64tz_dtype,
    is_extension_array_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_numeric_dtype,
    is_period_dtype,
    is_sparse,
    is_timedelta64_dtype,
    needs_i8_conversion,
)
from pandas.core.dtypes.missing import isna, maybe_fill

import pandas.core.algorithms as algorithms
from pandas.core.base import SelectionMixin
import pandas.core.common as com
from pandas.core.frame import DataFrame
from pandas.core.generic import NDFrame
from pandas.core.groupby import base, grouper
from pandas.core.indexes.api import Index, MultiIndex, ensure_index
from pandas.core.series import Series
from pandas.core.sorting import (
    compress_group_index,
    decons_obs_group_ids,
    get_flattened_list,
    get_group_index,
    get_group_index_sorter,
    get_indexer_dict,
)


class BaseGrouper:
    """
    This is an internal Grouper class, which actually holds
    the generated groups

    Parameters
    ----------
    axis : Index
    groupings : Sequence[Grouping]
        all the grouping instances to handle in this grouper
        for example for grouper list to groupby, need to pass the list
    sort : bool, default True
        whether this grouper will give sorted result or not
    group_keys : bool, default True
    mutated : bool, default False
    indexer : intp array, optional
        the indexer created by Grouper
        some groupers (TimeGrouper) will sort its axis and its
        group_info is also sorted, so need the indexer to reorder

    """

    def __init__(
        self,
        axis: Index,
        groupings: Sequence["grouper.Grouping"],
        sort: bool = True,
        group_keys: bool = True,
        mutated: bool = False,
        indexer: Optional[np.ndarray] = None,
        dropna: bool = True,
    ):
        assert isinstance(axis, Index), axis

        self._filter_empty_groups = self.compressed = len(groupings) != 1
        self.axis = axis
        self._groupings: List[grouper.Grouping] = list(groupings)
        self.sort = sort
        self.group_keys = group_keys
        self.mutated = mutated
        self.indexer = indexer
        self.dropna = dropna

    @property
    def groupings(self) -> List["grouper.Grouping"]:
        return self._groupings

    @property
    def shape(self) -> Shape:
        return tuple(ping.ngroups for ping in self.groupings)

    def __iter__(self):
        return iter(self.indices)

    @property
    def nkeys(self) -> int:
        return len(self.groupings)

    def get_iterator(
        self, data: FrameOrSeries, axis: int = 0
    ) -> Iterator[Tuple[Label, FrameOrSeries]]:
        """
        Groupby iterator

        Returns
        -------
        Generator yielding sequence of (name, subsetted object)
        for each group
        """
        splitter = self._get_splitter(data, axis=axis)
        keys = self._get_group_keys()
        for key, (i, group) in zip(keys, splitter):
            yield key, group.__finalize__(data, method="groupby")

    @final
    def _get_splitter(self, data: FrameOrSeries, axis: int = 0) -> "DataSplitter":
        """
        Returns
        -------
        Generator yielding subsetted objects

        __finalize__ has not been called for the subsetted objects returned.
        """
        comp_ids, _, ngroups = self.group_info
        return get_splitter(data, comp_ids, ngroups, axis=axis)

    def _get_grouper(self):
        """
        We are a grouper as part of another's groupings.

        We have a specific method of grouping, so cannot
        convert to a Index for our grouper.
        """
        return self.groupings[0].grouper

    @final
    def _get_group_keys(self):
        if len(self.groupings) == 1:
            return self.levels[0]
        else:
            comp_ids, _, ngroups = self.group_info

            # provide "flattened" iterator for multi-group setting
            return get_flattened_list(comp_ids, ngroups, self.levels, self.codes)

    @final
    def apply(self, f: F, data: FrameOrSeries, axis: int = 0):
        mutated = self.mutated
        splitter = self._get_splitter(data, axis=axis)
        group_keys = self._get_group_keys()
        result_values = None

        sdata: FrameOrSeries = splitter._get_sorted_data()
        if sdata.ndim == 2 and np.any(sdata.dtypes.apply(is_extension_array_dtype)):
            # calling splitter.fast_apply will raise TypeError via apply_frame_axis0
            #  if we pass EA instead of ndarray
            #  TODO: can we have a workaround for EAs backed by ndarray?
            pass

        elif (
            com.get_callable_name(f) not in base.plotting_methods
            and isinstance(splitter, FrameSplitter)
            and axis == 0
            # fast_apply/libreduction doesn't allow non-numpy backed indexes
            and not sdata.index._has_complex_internals
        ):
            try:
                result_values, mutated = splitter.fast_apply(f, sdata, group_keys)

            except libreduction.InvalidApply as err:
                # This Exception is raised if `f` triggers an exception
                # but it is preferable to raise the exception in Python.
                if "Let this error raise above us" not in str(err):
                    # TODO: can we infer anything about whether this is
                    #  worth-retrying in pure-python?
                    raise

            else:
                # If the fast apply path could be used we can return here.
                # Otherwise we need to fall back to the slow implementation.
                if len(result_values) == len(group_keys):
                    return group_keys, result_values, mutated

        for key, (i, group) in zip(group_keys, splitter):
            object.__setattr__(group, "name", key)

            # result_values is None if fast apply path wasn't taken
            # or fast apply aborted with an unexpected exception.
            # In either case, initialize the result list and perform
            # the slow iteration.
            if result_values is None:
                result_values = []

            # If result_values is not None we're in the case that the
            # fast apply loop was broken prematurely but we have
            # already the result for the first group which we can reuse.
            elif i == 0:
                continue

            # group might be modified
            group_axes = group.axes
            res = f(group)
            if not _is_indexed_like(res, group_axes, axis):
                mutated = True
            result_values.append(res)

        return group_keys, result_values, mutated

    @cache_readonly
    def indices(self):
        """ dict {group name -> group indices} """
        codes_list = [ping.codes for ping in self.groupings]
        keys = [ping.group_index for ping in self.groupings]
        return get_indexer_dict(codes_list, keys)

    @property
    def codes(self) -> List[np.ndarray]:
        return [ping.codes for ping in self.groupings]

    @property
    def levels(self) -> List[Index]:
        return [ping.group_index for ping in self.groupings]

    @property
    def names(self) -> List[Label]:
        return [ping.name for ping in self.groupings]

    @final
    def size(self) -> Series:
        """
        Compute group sizes.
        """
        ids, _, ngroup = self.group_info
        ids = ensure_platform_int(ids)
        if ngroup:
            out = np.bincount(ids[ids != -1], minlength=ngroup)
        else:
            out = []
        return Series(out, index=self.result_index, dtype="int64")

    @cache_readonly
    def groups(self) -> Dict[Hashable, np.ndarray]:
        """ dict {group name -> group labels} """
        if len(self.groupings) == 1:
            return self.groupings[0].groups
        else:
            to_groupby = zip(*(ping.grouper for ping in self.groupings))
            to_groupby = Index(to_groupby)
            return self.axis.groupby(to_groupby)

    @final
    @cache_readonly
    def is_monotonic(self) -> bool:
        # return if my group orderings are monotonic
        return Index(self.group_info[0]).is_monotonic

    @cache_readonly
    def group_info(self):
        comp_ids, obs_group_ids = self._get_compressed_codes()

        ngroups = len(obs_group_ids)
        comp_ids = ensure_int64(comp_ids)
        return comp_ids, obs_group_ids, ngroups

    @final
    @cache_readonly
    def codes_info(self) -> np.ndarray:
        # return the codes of items in original grouped axis
        codes, _, _ = self.group_info
        if self.indexer is not None:
            sorter = np.lexsort((codes, self.indexer))
            codes = codes[sorter]
        return codes

    @final
    def _get_compressed_codes(self) -> Tuple[np.ndarray, np.ndarray]:
        all_codes = self.codes
        if len(all_codes) > 1:
            group_index = get_group_index(all_codes, self.shape, sort=True, xnull=True)
            return compress_group_index(group_index, sort=self.sort)

        ping = self.groupings[0]
        return ping.codes, np.arange(len(ping.group_index))

    @final
    @cache_readonly
    def ngroups(self) -> int:
        return len(self.result_index)

    @property
    def reconstructed_codes(self) -> List[np.ndarray]:
        codes = self.codes
        comp_ids, obs_ids, _ = self.group_info
        return decons_obs_group_ids(comp_ids, obs_ids, self.shape, codes, xnull=True)

    @cache_readonly
    def result_index(self) -> Index:
        if not self.compressed and len(self.groupings) == 1:
            return self.groupings[0].result_index.rename(self.names[0])

        codes = self.reconstructed_codes
        levels = [ping.result_index for ping in self.groupings]
        return MultiIndex(
            levels=levels, codes=codes, verify_integrity=False, names=self.names
        )

    @final
    def get_group_levels(self) -> List[Index]:
        if not self.compressed and len(self.groupings) == 1:
            return [self.groupings[0].result_index]

        name_list = []
        for ping, codes in zip(self.groupings, self.reconstructed_codes):
            codes = ensure_platform_int(codes)
            levels = ping.result_index.take(codes)

            name_list.append(levels)

        return name_list

    # ------------------------------------------------------------
    # Aggregation functions

    _cython_functions = {
        "aggregate": {
            "add": "group_add",
            "prod": "group_prod",
            "min": "group_min",
            "max": "group_max",
            "mean": "group_mean",
            "median": "group_median",
            "var": "group_var",
            "first": "group_nth",
            "last": "group_last",
            "ohlc": "group_ohlc",
        },
        "transform": {
            "cumprod": "group_cumprod",
            "cumsum": "group_cumsum",
            "cummin": "group_cummin",
            "cummax": "group_cummax",
            "rank": "group_rank",
        },
    }

    _cython_arity = {"ohlc": 4}  # OHLC

    @final
    def _is_builtin_func(self, arg):
        """
        if we define a builtin function for this argument, return it,
        otherwise return the arg
        """
        return SelectionMixin._builtin_table.get(arg, arg)

    @final
    def _get_cython_function(
        self, kind: str, how: str, values: np.ndarray, is_numeric: bool
    ):

        dtype_str = values.dtype.name
        ftype = self._cython_functions[kind][how]

        # see if there is a fused-type version of function
        # only valid for numeric
        f = getattr(libgroupby, ftype, None)
        if f is not None and is_numeric:
            return f

        # otherwise find dtype-specific version, falling back to object
        for dt in [dtype_str, "object"]:
            f2 = getattr(libgroupby, f"{ftype}_{dt}", None)
            if f2 is not None:
                return f2

        if hasattr(f, "__signatures__"):
            # inspect what fused types are implemented
            if dtype_str == "object" and "object" not in f.__signatures__:
                # disallow this function so we get a NotImplementedError below
                #  instead of a TypeError at runtime
                f = None

        func = f

        if func is None:
            raise NotImplementedError(
                f"function is not implemented for this dtype: "
                f"[how->{how},dtype->{dtype_str}]"
            )

        return func

    @final
    def _get_cython_func_and_vals(
        self, kind: str, how: str, values: np.ndarray, is_numeric: bool
    ):
        """
        Find the appropriate cython function, casting if necessary.

        Parameters
        ----------
        kind : str
        how : str
        values : np.ndarray
        is_numeric : bool

        Returns
        -------
        func : callable
        values : np.ndarray
        """
        try:
            func = self._get_cython_function(kind, how, values, is_numeric)
        except NotImplementedError:
            if is_numeric:
                try:
                    values = ensure_float64(values)
                except TypeError:
                    if lib.infer_dtype(values, skipna=False) == "complex":
                        values = values.astype(complex)
                    else:
                        raise
                func = self._get_cython_function(kind, how, values, is_numeric)
            else:
                raise
        return func, values

    @final
    def _disallow_invalid_ops(self, values: ArrayLike, how: str):
        """
        Check if we can do this operation with our cython functions.

        Raises
        ------
        NotImplementedError
            This is either not a valid function for this dtype, or
            valid but not implemented in cython.
        """
        dtype = values.dtype

        if is_categorical_dtype(dtype) or is_sparse(dtype):
            # categoricals are only 1d, so we
            #  are not setup for dim transforming
            raise NotImplementedError(f"{dtype} dtype not supported")
        elif is_datetime64_any_dtype(dtype):
            # we raise NotImplemented if this is an invalid operation
            #  entirely, e.g. adding datetimes
            if how in ["add", "prod", "cumsum", "cumprod"]:
                raise NotImplementedError(
                    f"datetime64 type does not support {how} operations"
                )
        elif is_timedelta64_dtype(dtype):
            if how in ["prod", "cumprod"]:
                raise NotImplementedError(
                    f"timedelta64 type does not support {how} operations"
                )

    @final
    def _ea_wrap_cython_operation(
        self, kind: str, values, how: str, axis: int, min_count: int = -1, **kwargs
    ) -> Tuple[np.ndarray, Optional[List[str]]]:
        """
        If we have an ExtensionArray, unwrap, call _cython_operation, and
        re-wrap if appropriate.
        """
        # TODO: general case implementation overrideable by EAs.
        orig_values = values

        if is_datetime64tz_dtype(values.dtype) or is_period_dtype(values.dtype):
            # All of the functions implemented here are ordinal, so we can
            #  operate on the tz-naive equivalents
            values = values.view("M8[ns]")
            res_values = self._cython_operation(
                kind, values, how, axis, min_count, **kwargs
            )
            if how in ["rank"]:
                # preserve float64 dtype
                return res_values

            res_values = res_values.astype("i8", copy=False)
            result = type(orig_values)._simple_new(res_values, dtype=orig_values.dtype)
            return result

        elif is_integer_dtype(values.dtype) or is_bool_dtype(values.dtype):
            # IntegerArray or BooleanArray
            values = ensure_int_or_float(values)
            res_values = self._cython_operation(
                kind, values, how, axis, min_count, **kwargs
            )
            dtype = maybe_cast_result_dtype(orig_values.dtype, how)
            if is_extension_array_dtype(dtype):
                cls = dtype.construct_array_type()
                return cls._from_sequence(res_values, dtype=dtype)
            return res_values

        elif is_float_dtype(values.dtype):
            # FloatingArray
            values = values.to_numpy(values.dtype.numpy_dtype, na_value=np.nan)
            res_values = self._cython_operation(
                kind, values, how, axis, min_count, **kwargs
            )
            result = type(orig_values)._from_sequence(res_values)
            return result

        raise NotImplementedError(values.dtype)

    @final
    def _cython_operation(
        self, kind: str, values, how: str, axis: int, min_count: int = -1, **kwargs
    ) -> np.ndarray:
        """
        Returns the values of a cython operation.
        """
        orig_values = values
        assert kind in ["transform", "aggregate"]

        if values.ndim > 2:
            raise NotImplementedError("number of dimensions is currently limited to 2")
        elif values.ndim == 2:
            # Note: it is *not* the case that axis is always 0 for 1-dim values,
            #  as we can have 1D ExtensionArrays that we need to treat as 2D
            assert axis == 1, axis

        # can we do this operation with our cython functions
        # if not raise NotImplementedError
        self._disallow_invalid_ops(values, how)

        if is_extension_array_dtype(values.dtype):
            return self._ea_wrap_cython_operation(
                kind, values, how, axis, min_count, **kwargs
            )

        is_datetimelike = needs_i8_conversion(values.dtype)
        is_numeric = is_numeric_dtype(values.dtype)

        if is_datetimelike:
            values = values.view("int64")
            is_numeric = True
        elif is_bool_dtype(values.dtype):
            values = ensure_int_or_float(values)
        elif is_integer_dtype(values):
            # we use iNaT for the missing value on ints
            # so pre-convert to guard this condition
            if (values == iNaT).any():
                values = ensure_float64(values)
            else:
                values = ensure_int_or_float(values)
        elif is_numeric and not is_complex_dtype(values):
            values = ensure_float64(ensure_float(values))
        else:
            values = values.astype(object)

        arity = self._cython_arity.get(how, 1)

        vdim = values.ndim
        swapped = False
        if vdim == 1:
            values = values[:, None]
            out_shape = (self.ngroups, arity)
        else:
            if axis > 0:
                swapped = True
                assert axis == 1, axis
                values = values.T
            if arity > 1:
                raise NotImplementedError(
                    "arity of more than 1 is not supported for the 'how' argument"
                )
            out_shape = (self.ngroups,) + values.shape[1:]

        func, values = self._get_cython_func_and_vals(kind, how, values, is_numeric)

        if how == "rank":
            out_dtype = "float"
        else:
            if is_numeric:
                out_dtype = f"{values.dtype.kind}{values.dtype.itemsize}"
            else:
                out_dtype = "object"

        codes, _, _ = self.group_info

        if kind == "aggregate":
            result = maybe_fill(np.empty(out_shape, dtype=out_dtype), fill_value=np.nan)
            counts = np.zeros(self.ngroups, dtype=np.int64)
            result = self._aggregate(result, counts, values, codes, func, min_count)
        elif kind == "transform":
            result = maybe_fill(
                np.empty_like(values, dtype=out_dtype), fill_value=np.nan
            )

            # TODO: min_count
            result = self._transform(
                result, values, codes, func, is_datetimelike, **kwargs
            )

        if is_integer_dtype(result) and not is_datetimelike:
            mask = result == iNaT
            if mask.any():
                result = result.astype("float64")
                result[mask] = np.nan

        if kind == "aggregate" and self._filter_empty_groups and not counts.all():
            assert result.ndim != 2
            result = result[counts > 0]

        if vdim == 1 and arity == 1:
            result = result[:, 0]

        if swapped:
            result = result.swapaxes(0, axis)

        if how not in base.cython_cast_blocklist:
            # e.g. if we are int64 and need to restore to datetime64/timedelta64
            # "rank" is the only member of cython_cast_blocklist we get here
            dtype = maybe_cast_result_dtype(orig_values.dtype, how)
            result = maybe_downcast_to_dtype(result, dtype)

        return result

    @final
    def _aggregate(
        self, result, counts, values, comp_ids, agg_func, min_count: int = -1
    ):
        if agg_func is libgroupby.group_nth:
            # different signature from the others
            agg_func(result, counts, values, comp_ids, min_count, rank=1)
        else:
            agg_func(result, counts, values, comp_ids, min_count)

        return result

    @final
    def _transform(
        self, result, values, comp_ids, transform_func, is_datetimelike: bool, **kwargs
    ):

        comp_ids, _, ngroups = self.group_info
        transform_func(result, values, comp_ids, ngroups, is_datetimelike, **kwargs)

        return result

    def agg_series(self, obj: Series, func: F):
        # Caller is responsible for checking ngroups != 0
        assert self.ngroups != 0

        if len(obj) == 0:
            # SeriesGrouper would raise if we were to call _aggregate_series_fast
            return self._aggregate_series_pure_python(obj, func)

        elif is_extension_array_dtype(obj.dtype):
            # _aggregate_series_fast would raise TypeError when
            #  calling libreduction.Slider
            # In the datetime64tz case it would incorrectly cast to tz-naive
            # TODO: can we get a performant workaround for EAs backed by ndarray?
            return self._aggregate_series_pure_python(obj, func)

        elif obj.index._has_complex_internals:
            # Preempt TypeError in _aggregate_series_fast
            return self._aggregate_series_pure_python(obj, func)

        try:
            return self._aggregate_series_fast(obj, func)
        except ValueError as err:
            if "Must produce aggregated value" in str(err):
                # raised in libreduction
                pass
            else:
                raise
        return self._aggregate_series_pure_python(obj, func)

    @final
    def _aggregate_series_fast(self, obj: Series, func: F):
        # At this point we have already checked that
        #  - obj.index is not a MultiIndex
        #  - obj is backed by an ndarray, not ExtensionArray
        #  - len(obj) > 0
        #  - ngroups != 0
        func = self._is_builtin_func(func)

        group_index, _, ngroups = self.group_info

        # avoids object / Series creation overhead
        dummy = obj.iloc[:0]
        indexer = get_group_index_sorter(group_index, ngroups)
        obj = obj.take(indexer)
        group_index = algorithms.take_nd(group_index, indexer, allow_fill=False)
        grouper = libreduction.SeriesGrouper(obj, func, group_index, ngroups, dummy)
        result, counts = grouper.get_result()
        return result, counts

    @final
    def _aggregate_series_pure_python(self, obj: Series, func: F):
        group_index, _, ngroups = self.group_info

        counts = np.zeros(ngroups, dtype=int)
        result = np.empty(ngroups, dtype="O")
        initialized = False

        splitter = get_splitter(obj, group_index, ngroups, axis=0)

        for label, group in splitter:

            # Each step of this loop corresponds to
            #  libreduction._BaseGrouper._apply_to_group
            res = func(group)
            res = libreduction.extract_result(res)

            if not initialized:
                # We only do this validation on the first iteration
                libreduction.check_result_array(res, 0)
                initialized = True

            counts[label] = group.shape[0]
            result[label] = res

        result = lib.maybe_convert_objects(result, try_float=0)
        result = maybe_cast_result(result, obj, numeric_only=True)

        return result, counts


class BinGrouper(BaseGrouper):
    """
    This is an internal Grouper class

    Parameters
    ----------
    bins : the split index of binlabels to group the item of axis
    binlabels : the label list
    filter_empty : boolean, default False
    mutated : boolean, default False
    indexer : a intp array

    Examples
    --------
    bins: [2, 4, 6, 8, 10]
    binlabels: DatetimeIndex(['2005-01-01', '2005-01-03',
        '2005-01-05', '2005-01-07', '2005-01-09'],
        dtype='datetime64[ns]', freq='2D')

    the group_info, which contains the label of each item in grouped
    axis, the index of label in label list, group number, is

    (array([0, 0, 1, 1, 2, 2, 3, 3, 4, 4]), array([0, 1, 2, 3, 4]), 5)

    means that, the grouped axis has 10 items, can be grouped into 5
    labels, the first and second items belong to the first label, the
    third and forth items belong to the second label, and so on

    """

    def __init__(
        self,
        bins,
        binlabels,
        filter_empty: bool = False,
        mutated: bool = False,
        indexer=None,
    ):
        self.bins = ensure_int64(bins)
        self.binlabels = ensure_index(binlabels)
        self._filter_empty_groups = filter_empty
        self.mutated = mutated
        self.indexer = indexer

        # These lengths must match, otherwise we could call agg_series
        #  with empty self.bins, which would raise in libreduction.
        assert len(self.binlabels) == len(self.bins)

    @cache_readonly
    def groups(self):
        """ dict {group name -> group labels} """
        # this is mainly for compat
        # GH 3881
        result = {
            key: value
            for key, value in zip(self.binlabels, self.bins)
            if key is not NaT
        }
        return result

    @property
    def nkeys(self) -> int:
        return 1

    def _get_grouper(self):
        """
        We are a grouper as part of another's groupings.

        We have a specific method of grouping, so cannot
        convert to a Index for our grouper.
        """
        return self

    def get_iterator(self, data: FrameOrSeries, axis: int = 0):
        """
        Groupby iterator

        Returns
        -------
        Generator yielding sequence of (name, subsetted object)
        for each group
        """
        if axis == 0:
            slicer = lambda start, edge: data.iloc[start:edge]
        else:
            slicer = lambda start, edge: data.iloc[:, start:edge]

        length = len(data.axes[axis])

        start = 0
        for edge, label in zip(self.bins, self.binlabels):
            if label is not NaT:
                yield label, slicer(start, edge)
            start = edge

        if start < length:
            yield self.binlabels[-1], slicer(start, None)

    @cache_readonly
    def indices(self):
        indices = collections.defaultdict(list)

        i = 0
        for label, bin in zip(self.binlabels, self.bins):
            if i < bin:
                if label is not NaT:
                    indices[label] = list(range(i, bin))
                i = bin
        return indices

    @cache_readonly
    def group_info(self):
        ngroups = self.ngroups
        obs_group_ids = np.arange(ngroups)
        rep = np.diff(np.r_[0, self.bins])

        rep = ensure_platform_int(rep)
        if ngroups == len(self.bins):
            comp_ids = np.repeat(np.arange(ngroups), rep)
        else:
            comp_ids = np.repeat(np.r_[-1, np.arange(ngroups)], rep)

        return (
            comp_ids.astype("int64", copy=False),
            obs_group_ids.astype("int64", copy=False),
            ngroups,
        )

    @cache_readonly
    def reconstructed_codes(self) -> List[np.ndarray]:
        # get unique result indices, and prepend 0 as groupby starts from the first
        return [np.r_[0, np.flatnonzero(self.bins[1:] != self.bins[:-1]) + 1]]

    @cache_readonly
    def result_index(self):
        if len(self.binlabels) != 0 and isna(self.binlabels[0]):
            return self.binlabels[1:]

        return self.binlabels

    @property
    def levels(self) -> List[Index]:
        return [self.binlabels]

    @property
    def names(self) -> List[Label]:
        return [self.binlabels.name]

    @property
    def groupings(self) -> "List[grouper.Grouping]":
        return [
            grouper.Grouping(lvl, lvl, in_axis=False, level=None, name=name)
            for lvl, name in zip(self.levels, self.names)
        ]

    def agg_series(self, obj: Series, func: F):
        # Caller is responsible for checking ngroups != 0
        assert self.ngroups != 0
        assert len(self.bins) > 0  # otherwise we'd get IndexError in get_result

        if is_extension_array_dtype(obj.dtype):
            # preempt SeriesBinGrouper from raising TypeError
            return self._aggregate_series_pure_python(obj, func)

        dummy = obj[:0]
        grouper = libreduction.SeriesBinGrouper(obj, func, self.bins, dummy)
        return grouper.get_result()


def _is_indexed_like(obj, axes, axis: int) -> bool:
    if isinstance(obj, Series):
        if len(axes) > 1:
            return False
        return obj.axes[axis].equals(axes[axis])
    elif isinstance(obj, DataFrame):
        return obj.axes[axis].equals(axes[axis])

    return False


# ----------------------------------------------------------------------
# Splitting / application


class DataSplitter(Generic[FrameOrSeries]):
    def __init__(self, data: FrameOrSeries, labels, ngroups: int, axis: int = 0):
        self.data = data
        self.labels = ensure_int64(labels)
        self.ngroups = ngroups

        self.axis = axis
        assert isinstance(axis, int), axis

    @cache_readonly
    def slabels(self):
        # Sorted labels
        return algorithms.take_nd(self.labels, self.sort_idx, allow_fill=False)

    @cache_readonly
    def sort_idx(self):
        # Counting sort indexer
        return get_group_index_sorter(self.labels, self.ngroups)

    def __iter__(self):
        sdata = self._get_sorted_data()

        if self.ngroups == 0:
            # we are inside a generator, rather than raise StopIteration
            # we merely return signal the end
            return

        starts, ends = lib.generate_slices(self.slabels, self.ngroups)

        for i, (start, end) in enumerate(zip(starts, ends)):
            yield i, self._chop(sdata, slice(start, end))

    def _get_sorted_data(self) -> FrameOrSeries:
        return self.data.take(self.sort_idx, axis=self.axis)

    def _chop(self, sdata, slice_obj: slice) -> NDFrame:
        raise AbstractMethodError(self)


class SeriesSplitter(DataSplitter):
    def _chop(self, sdata: Series, slice_obj: slice) -> Series:
        # fastpath equivalent to `sdata.iloc[slice_obj]`
        mgr = sdata._mgr.get_slice(slice_obj)
        # __finalize__ not called here, must be applied by caller if applicable
        return sdata._constructor(mgr, name=sdata.name, fastpath=True)


class FrameSplitter(DataSplitter):
    def fast_apply(self, f: F, sdata: FrameOrSeries, names):
        # must return keys::list, values::list, mutated::bool
        starts, ends = lib.generate_slices(self.slabels, self.ngroups)
        return libreduction.apply_frame_axis0(sdata, f, names, starts, ends)

    def _chop(self, sdata: DataFrame, slice_obj: slice) -> DataFrame:
        # Fastpath equivalent to:
        # if self.axis == 0:
        #     return sdata.iloc[slice_obj]
        # else:
        #     return sdata.iloc[:, slice_obj]
        mgr = sdata._mgr.get_slice(slice_obj, axis=1 - self.axis)
        # __finalize__ not called here, must be applied by caller if applicable
        return sdata._constructor(mgr)


def get_splitter(
    data: FrameOrSeries, labels: np.ndarray, ngroups: int, axis: int = 0
) -> DataSplitter:
    if isinstance(data, Series):
        klass: Type[DataSplitter] = SeriesSplitter
    else:
        # i.e. DataFrame
        klass = FrameSplitter

    return klass(data, labels, ngroups, axis)
