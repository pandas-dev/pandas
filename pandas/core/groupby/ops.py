"""
Provide classes to perform the groupby aggregate operations.

These are not exposed to the user and provide implementations of the grouping
operations, primarily in cython. These classes (BaseGrouper and BinGrouper)
are contained *in* the SeriesGroupBy and DataFrameGroupBy objects.
"""

import copy
import collections
import numpy as np

from pandas._libs import lib, reduction, NaT, iNaT, groupby as libgroupby
from pandas.util._decorators import cache_readonly

from pandas.compat import zip, range, lzip

from pandas.core.base import SelectionMixin
from pandas.core.dtypes.missing import isna, _maybe_fill
from pandas.core.index import (
    Index, MultiIndex, ensure_index)
from pandas.core.dtypes.common import (
    ensure_float64,
    ensure_platform_int,
    ensure_int64,
    ensure_object,
    needs_i8_conversion,
    is_integer_dtype,
    is_complex_dtype,
    is_bool_dtype,
    is_numeric_dtype,
    is_timedelta64_dtype,
    is_datetime64_any_dtype,
    is_categorical_dtype)
from pandas.core.series import Series
from pandas.core.frame import DataFrame
from pandas.core.generic import NDFrame
import pandas.core.common as com
from pandas.core.groupby import base
from pandas.core.sorting import (get_group_index_sorter, get_group_index,
                                 compress_group_index, get_flattened_iterator,
                                 decons_obs_group_ids, get_indexer_dict)
import pandas.core.algorithms as algorithms


def generate_bins_generic(values, binner, closed):
    """
    Generate bin edge offsets and bin labels for one array using another array
    which has bin edge values. Both arrays must be sorted.

    Parameters
    ----------
    values : array of values
    binner : a comparable array of values representing bins into which to bin
        the first array. Note, 'values' end-points must fall within 'binner'
        end-points.
    closed : which end of bin is closed; left (default), right

    Returns
    -------
    bins : array of offsets (into 'values' argument) of bins.
        Zero and last edge are excluded in result, so for instance the first
        bin is values[0:bin[0]] and the last is values[bin[-1]:]
    """
    lenidx = len(values)
    lenbin = len(binner)

    if lenidx <= 0 or lenbin <= 0:
        raise ValueError("Invalid length for values or for binner")

    # check binner fits data
    if values[0] < binner[0]:
        raise ValueError("Values falls before first bin")

    if values[lenidx - 1] > binner[lenbin - 1]:
        raise ValueError("Values falls after last bin")

    bins = np.empty(lenbin - 1, dtype=np.int64)

    j = 0  # index into values
    bc = 0  # bin count

    # linear scan, presume nothing about values/binner except that it fits ok
    for i in range(0, lenbin - 1):
        r_bin = binner[i + 1]

        # count values in current bin, advance to next bin
        while j < lenidx and (values[j] < r_bin or
                              (closed == 'right' and values[j] == r_bin)):
            j += 1

        bins[bc] = j
        bc += 1

    return bins


class BaseGrouper(object):
    """
    This is an internal Grouper class, which actually holds
    the generated groups

    Parameters
    ----------
    axis : int
        the axis to group
    groupings : array of grouping
        all the grouping instances to handle in this grouper
        for example for grouper list to groupby, need to pass the list
    sort : boolean, default True
        whether this grouper will give sorted result or not
    group_keys : boolean, default True
    mutated : boolean, default False
    indexer : intp array, optional
        the indexer created by Grouper
        some groupers (TimeGrouper) will sort its axis and its
        group_info is also sorted, so need the indexer to reorder

    """

    def __init__(self, axis, groupings, sort=True, group_keys=True,
                 mutated=False, indexer=None):
        self._filter_empty_groups = self.compressed = len(groupings) != 1
        self.axis = axis
        self.groupings = groupings
        self.sort = sort
        self.group_keys = group_keys
        self.mutated = mutated
        self.indexer = indexer

    @property
    def shape(self):
        return tuple(ping.ngroups for ping in self.groupings)

    def __iter__(self):
        return iter(self.indices)

    @property
    def nkeys(self):
        return len(self.groupings)

    def get_iterator(self, data, axis=0):
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
            yield key, group

    def _get_splitter(self, data, axis=0):
        comp_ids, _, ngroups = self.group_info
        return get_splitter(data, comp_ids, ngroups, axis=axis)

    def _get_group_keys(self):
        if len(self.groupings) == 1:
            return self.levels[0]
        else:
            comp_ids, _, ngroups = self.group_info

            # provide "flattened" iterator for multi-group setting
            return get_flattened_iterator(comp_ids,
                                          ngroups,
                                          self.levels,
                                          self.labels)

    def apply(self, f, data, axis=0):
        mutated = self.mutated
        splitter = self._get_splitter(data, axis=axis)
        group_keys = self._get_group_keys()

        # oh boy
        f_name = com.get_callable_name(f)
        if (f_name not in base.plotting_methods and
                hasattr(splitter, 'fast_apply') and axis == 0):
            try:
                values, mutated = splitter.fast_apply(f, group_keys)
                return group_keys, values, mutated
            except reduction.InvalidApply:
                # we detect a mutation of some kind
                # so take slow path
                pass
            except Exception:
                # raise this error to the caller
                pass

        result_values = []
        for key, (i, group) in zip(group_keys, splitter):
            object.__setattr__(group, 'name', key)

            # group might be modified
            group_axes = _get_axes(group)
            res = f(group)
            if not _is_indexed_like(res, group_axes):
                mutated = True
            result_values.append(res)

        return group_keys, result_values, mutated

    @cache_readonly
    def indices(self):
        """ dict {group name -> group indices} """
        if len(self.groupings) == 1:
            return self.groupings[0].indices
        else:
            label_list = [ping.labels for ping in self.groupings]
            keys = [com.values_from_object(ping.group_index)
                    for ping in self.groupings]
            return get_indexer_dict(label_list, keys)

    @property
    def labels(self):
        return [ping.labels for ping in self.groupings]

    @property
    def levels(self):
        return [ping.group_index for ping in self.groupings]

    @property
    def names(self):
        return [ping.name for ping in self.groupings]

    def size(self):
        """
        Compute group sizes

        """
        ids, _, ngroup = self.group_info
        ids = ensure_platform_int(ids)
        if ngroup:
            out = np.bincount(ids[ids != -1], minlength=ngroup)
        else:
            out = ids
        return Series(out,
                      index=self.result_index,
                      dtype='int64')

    @cache_readonly
    def groups(self):
        """ dict {group name -> group labels} """
        if len(self.groupings) == 1:
            return self.groupings[0].groups
        else:
            to_groupby = lzip(*(ping.grouper for ping in self.groupings))
            to_groupby = Index(to_groupby)
            return self.axis.groupby(to_groupby)

    @cache_readonly
    def is_monotonic(self):
        # return if my group orderings are monotonic
        return Index(self.group_info[0]).is_monotonic

    @cache_readonly
    def group_info(self):
        comp_ids, obs_group_ids = self._get_compressed_labels()

        ngroups = len(obs_group_ids)
        comp_ids = ensure_int64(comp_ids)
        return comp_ids, obs_group_ids, ngroups

    @cache_readonly
    def label_info(self):
        # return the labels of items in original grouped axis
        labels, _, _ = self.group_info
        if self.indexer is not None:
            sorter = np.lexsort((labels, self.indexer))
            labels = labels[sorter]
        return labels

    def _get_compressed_labels(self):
        all_labels = [ping.labels for ping in self.groupings]
        if len(all_labels) > 1:
            group_index = get_group_index(all_labels, self.shape,
                                          sort=True, xnull=True)
            return compress_group_index(group_index, sort=self.sort)

        ping = self.groupings[0]
        return ping.labels, np.arange(len(ping.group_index))

    @cache_readonly
    def ngroups(self):
        return len(self.result_index)

    @property
    def recons_labels(self):
        comp_ids, obs_ids, _ = self.group_info
        labels = (ping.labels for ping in self.groupings)
        return decons_obs_group_ids(
            comp_ids, obs_ids, self.shape, labels, xnull=True)

    @cache_readonly
    def result_index(self):
        if not self.compressed and len(self.groupings) == 1:
            return self.groupings[0].result_index.rename(self.names[0])

        labels = self.recons_labels
        levels = [ping.result_index for ping in self.groupings]
        result = MultiIndex(levels=levels,
                            labels=labels,
                            verify_integrity=False,
                            names=self.names)
        return result

    def get_group_levels(self):
        if not self.compressed and len(self.groupings) == 1:
            return [self.groupings[0].result_index]

        name_list = []
        for ping, labels in zip(self.groupings, self.recons_labels):
            labels = ensure_platform_int(labels)
            levels = ping.result_index.take(labels)

            name_list.append(levels)

        return name_list

    # ------------------------------------------------------------
    # Aggregation functions

    _cython_functions = {
        'aggregate': {
            'add': 'group_add',
            'prod': 'group_prod',
            'min': 'group_min',
            'max': 'group_max',
            'mean': 'group_mean',
            'median': {
                'name': 'group_median'
            },
            'var': 'group_var',
            'first': {
                'name': 'group_nth',
                'f': lambda func, a, b, c, d, e: func(a, b, c, d, 1, -1)
            },
            'last': 'group_last',
            'ohlc': 'group_ohlc',
        },

        'transform': {
            'cumprod': 'group_cumprod',
            'cumsum': 'group_cumsum',
            'cummin': 'group_cummin',
            'cummax': 'group_cummax',
            'rank': {
                'name': 'group_rank',
                'f': lambda func, a, b, c, d, **kwargs: func(
                    a, b, c, d,
                    kwargs.get('ties_method', 'average'),
                    kwargs.get('ascending', True),
                    kwargs.get('pct', False),
                    kwargs.get('na_option', 'keep')
                )
            }
        }
    }

    _cython_arity = {
        'ohlc': 4,  # OHLC
    }

    _name_functions = {
        'ohlc': lambda *args: ['open', 'high', 'low', 'close']
    }

    def _is_builtin_func(self, arg):
        """
        if we define an builtin function for this argument, return it,
        otherwise return the arg
        """
        return SelectionMixin._builtin_table.get(arg, arg)

    def _get_cython_function(self, kind, how, values, is_numeric):

        dtype_str = values.dtype.name

        def get_func(fname):
            # see if there is a fused-type version of function
            # only valid for numeric
            f = getattr(libgroupby, fname, None)
            if f is not None and is_numeric:
                return f

            # otherwise find dtype-specific version, falling back to object
            for dt in [dtype_str, 'object']:
                f = getattr(libgroupby, "%s_%s" % (fname, dtype_str), None)
                if f is not None:
                    return f

        ftype = self._cython_functions[kind][how]

        if isinstance(ftype, dict):
            func = afunc = get_func(ftype['name'])

            # a sub-function
            f = ftype.get('f')
            if f is not None:

                def wrapper(*args, **kwargs):
                    return f(afunc, *args, **kwargs)

                # need to curry our sub-function
                func = wrapper

        else:
            func = get_func(ftype)

        if func is None:
            raise NotImplementedError("function is not implemented for this"
                                      "dtype: [how->%s,dtype->%s]" %
                                      (how, dtype_str))
        return func

    def _cython_operation(self, kind, values, how, axis, min_count=-1,
                          **kwargs):
        assert kind in ['transform', 'aggregate']

        # can we do this operation with our cython functions
        # if not raise NotImplementedError

        # we raise NotImplemented if this is an invalid operation
        # entirely, e.g. adding datetimes

        # categoricals are only 1d, so we
        # are not setup for dim transforming
        if is_categorical_dtype(values):
            raise NotImplementedError(
                "categoricals are not support in cython ops ATM")
        elif is_datetime64_any_dtype(values):
            if how in ['add', 'prod', 'cumsum', 'cumprod']:
                raise NotImplementedError(
                    "datetime64 type does not support {} "
                    "operations".format(how))
        elif is_timedelta64_dtype(values):
            if how in ['prod', 'cumprod']:
                raise NotImplementedError(
                    "timedelta64 type does not support {} "
                    "operations".format(how))

        arity = self._cython_arity.get(how, 1)

        vdim = values.ndim
        swapped = False
        if vdim == 1:
            values = values[:, None]
            out_shape = (self.ngroups, arity)
        else:
            if axis > 0:
                swapped = True
                values = values.swapaxes(0, axis)
            if arity > 1:
                raise NotImplementedError("arity of more than 1 is not "
                                          "supported for the 'how' argument")
            out_shape = (self.ngroups,) + values.shape[1:]

        is_datetimelike = needs_i8_conversion(values.dtype)
        is_numeric = is_numeric_dtype(values.dtype)

        if is_datetimelike:
            values = values.view('int64')
            is_numeric = True
        elif is_bool_dtype(values.dtype):
            values = ensure_float64(values)
        elif is_integer_dtype(values):
            # we use iNaT for the missing value on ints
            # so pre-convert to guard this condition
            if (values == iNaT).any():
                values = ensure_float64(values)
            else:
                values = values.astype('int64', copy=False)
        elif is_numeric and not is_complex_dtype(values):
            values = ensure_float64(values)
        else:
            values = values.astype(object)

        try:
            func = self._get_cython_function(
                kind, how, values, is_numeric)
        except NotImplementedError:
            if is_numeric:
                values = ensure_float64(values)
                func = self._get_cython_function(
                    kind, how, values, is_numeric)
            else:
                raise

        if how == 'rank':
            out_dtype = 'float'
        else:
            if is_numeric:
                out_dtype = '%s%d' % (values.dtype.kind, values.dtype.itemsize)
            else:
                out_dtype = 'object'

        labels, _, _ = self.group_info

        if kind == 'aggregate':
            result = _maybe_fill(np.empty(out_shape, dtype=out_dtype),
                                 fill_value=np.nan)
            counts = np.zeros(self.ngroups, dtype=np.int64)
            result = self._aggregate(
                result, counts, values, labels, func, is_numeric,
                is_datetimelike, min_count)
        elif kind == 'transform':
            result = _maybe_fill(np.empty_like(values, dtype=out_dtype),
                                 fill_value=np.nan)

            # TODO: min_count
            result = self._transform(
                result, values, labels, func, is_numeric, is_datetimelike,
                **kwargs)

        if is_integer_dtype(result) and not is_datetimelike:
            mask = result == iNaT
            if mask.any():
                result = result.astype('float64')
                result[mask] = np.nan

        if kind == 'aggregate' and \
           self._filter_empty_groups and not counts.all():
            if result.ndim == 2:
                try:
                    result = lib.row_bool_subset(
                        result, (counts > 0).view(np.uint8))
                except ValueError:
                    result = lib.row_bool_subset_object(
                        ensure_object(result),
                        (counts > 0).view(np.uint8))
            else:
                result = result[counts > 0]

        if vdim == 1 and arity == 1:
            result = result[:, 0]

        if how in self._name_functions:
            # TODO
            names = self._name_functions[how]()
        else:
            names = None

        if swapped:
            result = result.swapaxes(0, axis)

        return result, names

    def aggregate(self, values, how, axis=0, min_count=-1):
        return self._cython_operation('aggregate', values, how, axis,
                                      min_count=min_count)

    def transform(self, values, how, axis=0, **kwargs):
        return self._cython_operation('transform', values, how, axis, **kwargs)

    def _aggregate(self, result, counts, values, comp_ids, agg_func,
                   is_numeric, is_datetimelike, min_count=-1):
        if values.ndim > 3:
            # punting for now
            raise NotImplementedError("number of dimensions is currently "
                                      "limited to 3")
        elif values.ndim > 2:
            for i, chunk in enumerate(values.transpose(2, 0, 1)):

                chunk = chunk.squeeze()
                agg_func(result[:, :, i], counts, chunk, comp_ids,
                         min_count)
        else:
            agg_func(result, counts, values, comp_ids, min_count)

        return result

    def _transform(self, result, values, comp_ids, transform_func,
                   is_numeric, is_datetimelike, **kwargs):

        comp_ids, _, ngroups = self.group_info
        if values.ndim > 3:
            # punting for now
            raise NotImplementedError("number of dimensions is currently "
                                      "limited to 3")
        elif values.ndim > 2:
            for i, chunk in enumerate(values.transpose(2, 0, 1)):

                transform_func(result[:, :, i], values,
                               comp_ids, is_datetimelike, **kwargs)
        else:
            transform_func(result, values, comp_ids, is_datetimelike, **kwargs)

        return result

    def agg_series(self, obj, func):
        try:
            return self._aggregate_series_fast(obj, func)
        except Exception:
            return self._aggregate_series_pure_python(obj, func)

    def _aggregate_series_fast(self, obj, func):
        func = self._is_builtin_func(func)

        if obj.index._has_complex_internals:
            raise TypeError('Incompatible index for Cython grouper')

        group_index, _, ngroups = self.group_info

        # avoids object / Series creation overhead
        dummy = obj._get_values(slice(None, 0)).to_dense()
        indexer = get_group_index_sorter(group_index, ngroups)
        obj = obj._take(indexer).to_dense()
        group_index = algorithms.take_nd(
            group_index, indexer, allow_fill=False)
        grouper = reduction.SeriesGrouper(obj, func, group_index, ngroups,
                                          dummy)
        result, counts = grouper.get_result()
        return result, counts

    def _aggregate_series_pure_python(self, obj, func):

        group_index, _, ngroups = self.group_info

        counts = np.zeros(ngroups, dtype=int)
        result = None

        splitter = get_splitter(obj, group_index, ngroups, axis=self.axis)

        for label, group in splitter:
            res = func(group)
            if result is None:
                if (isinstance(res, (Series, Index, np.ndarray))):
                    raise ValueError('Function does not reduce')
                result = np.empty(ngroups, dtype='O')

            counts[label] = group.shape[0]
            result[label] = res

        result = lib.maybe_convert_objects(result, try_float=0)
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

    def __init__(self, bins, binlabels, filter_empty=False, mutated=False,
                 indexer=None):
        self.bins = ensure_int64(bins)
        self.binlabels = ensure_index(binlabels)
        self._filter_empty_groups = filter_empty
        self.mutated = mutated
        self.indexer = indexer

    @cache_readonly
    def groups(self):
        """ dict {group name -> group labels} """

        # this is mainly for compat
        # GH 3881
        result = {}
        for key, value in zip(self.binlabels, self.bins):
            if key is not NaT:
                result[key] = value
        return result

    @property
    def nkeys(self):
        return 1

    def get_iterator(self, data, axis=0):
        """
        Groupby iterator

        Returns
        -------
        Generator yielding sequence of (name, subsetted object)
        for each group
        """
        if isinstance(data, NDFrame):
            slicer = lambda start, edge: data._slice(
                slice(start, edge), axis=axis)
            length = len(data.axes[axis])
        else:
            slicer = lambda start, edge: data[slice(start, edge)]
            length = len(data)

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

        return comp_ids.astype('int64', copy=False), \
            obs_group_ids.astype('int64', copy=False), ngroups

    @cache_readonly
    def ngroups(self):
        return len(self.result_index)

    @cache_readonly
    def result_index(self):
        if len(self.binlabels) != 0 and isna(self.binlabels[0]):
            return self.binlabels[1:]

        return self.binlabels

    @property
    def levels(self):
        return [self.binlabels]

    @property
    def names(self):
        return [self.binlabels.name]

    @property
    def groupings(self):
        from pandas.core.groupby.grouper import Grouping
        return [Grouping(lvl, lvl, in_axis=False, level=None, name=name)
                for lvl, name in zip(self.levels, self.names)]

    def agg_series(self, obj, func):
        dummy = obj[:0]
        grouper = reduction.SeriesBinGrouper(obj, func, self.bins, dummy)
        return grouper.get_result()

    # ----------------------------------------------------------------------
    # cython aggregation

    _cython_functions = copy.deepcopy(BaseGrouper._cython_functions)


def _get_axes(group):
    if isinstance(group, Series):
        return [group.index]
    else:
        return group.axes


def _is_indexed_like(obj, axes):
    if isinstance(obj, Series):
        if len(axes) > 1:
            return False
        return obj.index.equals(axes[0])
    elif isinstance(obj, DataFrame):
        return obj.index.equals(axes[0])

    return False


# ----------------------------------------------------------------------
# Splitting / application


class DataSplitter(object):

    def __init__(self, data, labels, ngroups, axis=0):
        self.data = data
        self.labels = ensure_int64(labels)
        self.ngroups = ngroups

        self.axis = axis

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
            # Since I'm now compressing the group ids, it's now not "possible"
            # to produce empty slices because such groups would not be observed
            # in the data
            # if start >= end:
            #     raise AssertionError('Start %s must be less than end %s'
            #                          % (str(start), str(end)))
            yield i, self._chop(sdata, slice(start, end))

    def _get_sorted_data(self):
        return self.data._take(self.sort_idx, axis=self.axis)

    def _chop(self, sdata, slice_obj):
        return sdata.iloc[slice_obj]

    def apply(self, f):
        raise com.AbstractMethodError(self)


class SeriesSplitter(DataSplitter):

    def _chop(self, sdata, slice_obj):
        return sdata._get_values(slice_obj).to_dense()


class FrameSplitter(DataSplitter):

    def __init__(self, data, labels, ngroups, axis=0):
        super(FrameSplitter, self).__init__(data, labels, ngroups, axis=axis)

    def fast_apply(self, f, names):
        # must return keys::list, values::list, mutated::bool
        try:
            starts, ends = lib.generate_slices(self.slabels, self.ngroups)
        except Exception:
            # fails when all -1
            return [], True

        sdata = self._get_sorted_data()
        results, mutated = reduction.apply_frame_axis0(sdata, f, names,
                                                       starts, ends)

        return results, mutated

    def _chop(self, sdata, slice_obj):
        if self.axis == 0:
            return sdata.iloc[slice_obj]
        else:
            return sdata._slice(slice_obj, axis=1)  # .loc[:, slice_obj]


class NDFrameSplitter(DataSplitter):

    def __init__(self, data, labels, ngroups, axis=0):
        super(NDFrameSplitter, self).__init__(data, labels, ngroups, axis=axis)

        self.factory = data._constructor

    def _get_sorted_data(self):
        # this is the BlockManager
        data = self.data._data

        # this is sort of wasteful but...
        sorted_axis = data.axes[self.axis].take(self.sort_idx)
        sorted_data = data.reindex_axis(sorted_axis, axis=self.axis)

        return sorted_data

    def _chop(self, sdata, slice_obj):
        return self.factory(sdata.get_slice(slice_obj, axis=self.axis))


def get_splitter(data, *args, **kwargs):
    if isinstance(data, Series):
        klass = SeriesSplitter
    elif isinstance(data, DataFrame):
        klass = FrameSplitter
    else:
        klass = NDFrameSplitter

    return klass(data, *args, **kwargs)
