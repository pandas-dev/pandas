from itertools import izip
import sys
import types

import numpy as np

from pandas.core.frame import DataFrame
from pandas.core.generic import NDFrame, PandasObject
from pandas.core.index import Factor, Index, MultiIndex
from pandas.core.internals import BlockManager
from pandas.core.series import Series
from pandas.core.panel import WidePanel
from pandas.util.decorators import cache_readonly
import pandas._tseries as _tseries


def groupby(obj, grouper, **kwds):
    """
    Intercepts creation and dispatches to the appropriate class based
    on type.
    """
    if isinstance(obj, Series):
        klass = SeriesGroupBy
    elif isinstance(obj, DataFrame):
        klass = DataFrameGroupBy
    else: # pragma: no cover
        raise TypeError('invalid type: %s' % type(obj))

    return klass(obj, grouper, **kwds)

class GroupBy(object):
    """
    Class for grouping and aggregating relational data.

    Supported classes: Series, DataFrame
    """
    def __init__(self, obj, grouper=None, axis=0, level=None,
                 groupings=None, exclusions=None, name=None):
        self._name = name
        self.obj = obj
        self.axis = axis
        self.level = level

        if groupings is None:
            groupings, exclusions = _get_groupings(obj, grouper, axis=axis,
                                                   level=level)

        self.groupings = groupings
        self.exclusions = set(exclusions)

    @property
    def groups(self):
        if len(self.groupings) == 1:
            return self.primary.groups
        else:
            raise NotImplementedError

    @property
    def name(self):
        if self._name is None:
            return 'result'
        else:
            return self._name

    @property
    def _obj_with_exclusions(self):
        return self.obj

    @property
    def _group_shape(self):
        return tuple(len(ping.ids) for ping in self.groupings)

    @property
    def _agg_stride_shape(self):
        raise NotImplementedError

    def __getattribute__(self, attr):
        try:
            return object.__getattribute__(self, attr)
        except AttributeError:
            if hasattr(self.obj, attr):
                return self._make_wrapper(attr)
            raise

    def _make_wrapper(self, name):
        f = getattr(self.obj, name)
        if not isinstance(f, types.MethodType):
            return self.apply(lambda self: getattr(self, name))

        f = getattr(type(self.obj), name)

        def wrapper(*args, **kwargs):
            def curried(x):
                return f(x, *args, **kwargs)
            return self.apply(curried)

        return wrapper

    @property
    def primary(self):
        return self.groupings[0]

    def get_group(self, name, obj=None):
        if obj is None:
            obj = self.obj

        labels = self.primary.get_group_labels(name)
        axis_name = obj._get_axis_name(self.axis)
        return obj.reindex(**{axis_name : labels})

    def __iter__(self):
        """
        Groupby iterator

        Returns
        -------
        Generator yielding sequence of (name, subsetted object)
        for each group
        """
        if len(self.groupings) == 1:
            groups = self.primary.indices.keys()
            try:
                groups = sorted(groups)
            except Exception: # pragma: no cover
                pass

            for name in groups:
                yield name, self.get_group(name)
        else:
            # provide "flattened" iterator for multi-group setting
            for it in self._multi_iter():
                yield it

    def _multi_iter(self):
        tipo = type(self.obj)
        data = self.obj
        if (isinstance(self.obj, NDFrame) and
            not isinstance(self.obj, DataFrame)):
            data = self.obj._data
        elif isinstance(self.obj, Series):
            tipo = Series

        def flatten(gen, level=0):
            ids = self.groupings[level].ids
            for cat, subgen in gen:
                if isinstance(subgen, tipo):
                    yield (ids[cat],), subgen
                else:
                    for subcat, data in flatten(subgen, level=level+1):
                        if len(data) == 0:
                            continue
                        yield (ids[cat],) + subcat, data

        gen = self._generator_factory(data)

        for cats, data in flatten(gen):
            yield cats + (data,)

    def apply(self, func):
        """
        Apply function, combine results together
        """
        return self._python_apply_general(func)

    def aggregate(self, func):
        raise NotImplementedError

    def agg(self, func):
        return self.aggregate(func)

    def _get_names(self):
        axes = [ping.levels for ping in self.groupings]
        grouping_names = [ping.name for ping in self.groupings]
        shape = self._group_shape
        return zip(grouping_names, _ravel_names(axes, shape))

    def _iterate_slices(self):
        yield self.name, self.obj

    def transform(self, func):
        raise NotImplementedError

    def mean(self):
        """
        Compute mean of groups, excluding missing values
        """
        try:
            return self._cython_agg_general('mean')
        except Exception:
            return self.aggregate(np.mean)

    def sum(self):
        """
        Compute sum of values, excluding missing values
        """
        try:
            return self._cython_agg_general('add')
        except Exception:
            return self.aggregate(np.sum)

    def _cython_agg_general(self, how):
        label_list = [ping.labels for ping in self.groupings]
        shape = self._group_shape

        # TODO: address inefficiencies, like duplicating effort (should
        # aggregate all the columns at once?)

        output = {}
        cannot_agg = []
        for name, obj in self._iterate_slices():
            try:
                obj = np.asarray(obj, dtype=float)
            except ValueError:
                cannot_agg.append(name)
                continue

            result, counts =  _tseries.group_aggregate(obj, label_list,
                                                       shape, how=how)
            result = result.ravel()
            mask = counts.ravel() > 0
            output[name] = result[mask]

        return self._wrap_aggregated_output(output, mask)

    def _get_multi_index(self, mask):
        name_list = self._get_names()
        masked = [raveled[mask] for _, raveled in name_list]
        return MultiIndex.from_arrays(masked)

    def _python_agg_general(self, arg):
        group_shape = self._group_shape
        counts = np.zeros(group_shape, dtype=int)

        # want to cythonize?
        def _doit(reschunk, ctchunk, gen, shape_axis=0):
            for i, (_, subgen) in enumerate(gen):
                if isinstance(subgen, PandasObject):
                    size = subgen.shape[shape_axis]
                    ctchunk[i] = size
                    if size == 0:
                        continue
                    reschunk[i] = arg(subgen)
                else:
                    _doit(reschunk[i], ctchunk[i], subgen,
                          shape_axis=shape_axis)

        gen_factory = self._generator_factory

        try:
            stride_shape = self._agg_stride_shape
            output = np.empty(group_shape + stride_shape,
                              dtype=float)
            output.fill(np.nan)
            obj = self._obj_with_exclusions
            _doit(output, counts, gen_factory(obj),
                  shape_axis=self.axis)

            mask = counts.ravel() > 0
            output = output.reshape((np.prod(group_shape),) + stride_shape)
            output = output[mask]
        except TypeError:
            result = np.empty(group_shape, dtype=float)
            result.fill(np.nan)
            # iterate through "columns" ex exclusions to populate output dict
            output = {}
            for name, obj in self._iterate_slices():
                _doit(result, counts, gen_factory(obj))
                # TODO: same mask for every column...
                result.fill(np.nan)
                output[name] = result.ravel()

            mask = counts.ravel() > 0
            for name, result in output.iteritems():
                output[name] = result[mask]

        return self._wrap_aggregated_output(output, mask)

    def _python_apply_general(self, arg):
        result_keys = []
        result_values = []

        key_as_tuple = len(self.groupings) > 1

        not_indexed_same = False

        for data in self:
            if key_as_tuple:
                key = data[:-1]
            else:
                key = data[0]

            group = data[-1]
            group.name = key

            res = arg(group)

            if not _is_indexed_like(res, group):
               not_indexed_same = True

            result_keys.append(key)
            result_values.append(res)

        return self._wrap_applied_output(result_keys, result_values,
                                         not_indexed_same=not_indexed_same)

    def _wrap_applied_output(self, *args, **kwargs):
        raise NotImplementedError

    @property
    def _generator_factory(self):
        labels = [ping.labels for ping in self.groupings]
        shape = self._group_shape

        if isinstance(self.obj, NDFrame):
            factory = self.obj._constructor
        else:
            factory = None

        factory = None
        axis = self.axis
        return lambda obj: generate_groups(obj, labels, shape, axis=axis,
                                           factory=factory)

def _is_indexed_like(obj, other):
    if isinstance(obj, Series):
        return obj.index.equals(other.index)
    elif isinstance(obj, DataFrame):
        return obj._indexed_same(other)

class Grouping(object):

    def __init__(self, index, grouper=None, name=None, level=None):
        self.name = name
        self.level = level
        self.grouper = _convert_grouper(index, grouper)

        if level is not None:
            inds = index.labels[level]
            labels = index.levels[level].values.take(inds)

            if grouper is not None:
                self.grouper = _tseries.arrmap(labels, self.grouper)
            else:
                self.grouper = labels

        self.index = index.values

        # no level passed
        if not isinstance(self.grouper, np.ndarray):
            self.grouper = _tseries.arrmap(self.index, self.grouper)

        self.indices = _tseries.groupby_indices(self.grouper)

    def __repr__(self):
        return 'Grouping(%s)' % self.name

    def __iter__(self):
        return iter(self.indices)

    def get_group_labels(self, group):
        inds = self.indices[group]
        return self.index.take(inds)

    _labels = None
    _ids = None
    _counts = None

    @property
    def labels(self):
        if self._labels is None:
            self._make_labels()
        return self._labels

    @property
    def ids(self):
        if self._ids is None:
            self._make_labels()
        return self._ids

    @property
    def levels(self):
        return [self.ids[k] for k in sorted(self.ids)]

    @property
    def counts(self):
        if self._counts is None:
            self._make_labels()
        return self._counts

    def _make_labels(self):
        ids, labels, counts  = _tseries.group_labels(self.grouper)
        sids, slabels, scounts = sort_group_labels(ids, labels, counts)
        self._labels = slabels
        self._ids = sids
        self._counts = scounts

    _groups = None
    @property
    def groups(self):
        if self._groups is None:
            self._groups = _tseries.groupby(self.index, self.grouper)
        return self._groups

def _get_groupings(obj, grouper=None, axis=0, level=None):
    group_axis = obj._get_axis(axis)

    if level is not None and not isinstance(group_axis, MultiIndex):
        raise ValueError('can only specify level with multi-level index')

    groupings = []
    exclusions = []
    if isinstance(grouper, (tuple, list)):
        for i, arg in enumerate(grouper):
            name = 'key_%d' % i
            if isinstance(arg, basestring):
                exclusions.append(arg)
                name = arg
                arg = obj[arg]

            ping = Grouping(group_axis, arg, name=name, level=level)
            groupings.append(ping)
    else:
        name = 'key'
        if isinstance(grouper, basestring):
            exclusions.append(grouper)
            name = grouper
            grouper = obj[grouper]
        ping = Grouping(group_axis, grouper, name=name, level=level)
        groupings.append(ping)

    return groupings, exclusions

def _convert_grouper(axis, grouper):
    if isinstance(grouper, dict):
        return grouper.get
    elif isinstance(grouper, Series):
        if grouper.index.equals(axis):
            return np.asarray(grouper, dtype=object)
        else:
            return grouper.__getitem__
    elif isinstance(grouper, (list, np.ndarray)):
        assert(len(grouper) == len(axis))
        return np.asarray(grouper, dtype=object)
    else:
        return grouper

class SeriesGroupBy(GroupBy):

    _cythonized_methods = set(['add', 'mean'])

    @property
    def _agg_stride_shape(self):
        return ()

    def get_group(self, name, obj=None):
        # faster get_group for Series
        if obj is None:
            obj = self.obj

        inds = self.primary.indices[name]
        return obj.take(inds)

    def aggregate(self, arg):
        """
        See doc for DataFrame.groupby, group series using mapper (dict or key
        function, apply given function to group, return result as series).

        Main difference here is that applyfunc must return a value, so that the
        result is a sensible series.

        Parameters
        ----------
        mapper : function
            Called on each element of the Series index to determine the groups
        arg : function
            Function to use to aggregate each group

        Returns
        -------
        Series or DataFrame
        """
        if isinstance(arg, basestring):
            return getattr(self, arg)()

        if len(self.groupings) > 1:
            return self._python_agg_general(arg)

        if hasattr(arg,'__iter__'):
            ret = self._aggregate_multiple_funcs(arg)
        else:
            try:
                result = self._aggregate_simple(arg)
            except Exception:
                result = self._aggregate_named(arg)

            if len(result) > 0:
                if isinstance(result.values()[0], Series):
                    ret = DataFrame(result).T
                else:
                    ret = Series(result)
            else:
                ret = Series({})

        return ret

    def _wrap_aggregated_output(self, output, mask):
        if isinstance(output, dict):
            # sort of a kludge
            output = output[self.name]

        if len(self.groupings) > 1:
            index = self._get_multi_index(mask)
            return Series(output, index=index)
        else:
            name_list = self._get_names()
            return Series(output, index=name_list[0][1])

    def _wrap_applied_output(self, keys, values, not_indexed_same=False):
        if len(keys) == 0:
            return Series([])

        if isinstance(values[0], Series):
            if not_indexed_same:
                data_dict = dict(zip(keys, values))
                result = DataFrame(data_dict).T
                if len(self.groupings) > 1:
                    result.index = MultiIndex.from_tuples(keys)
                return result
            else:
                cat_values = np.concatenate([x.values for x in values])
                cat_index = np.concatenate([np.asarray(x.index)
                                            for x in values])
                return Series(cat_values, index=cat_index)
        else:
            if len(self.groupings) > 1:
                index = MultiIndex.from_tuples(keys)
                return Series(values, index)
            else:
                return Series(values, keys)

    def _aggregate_multiple_funcs(self, arg):
        if not isinstance(arg, dict):
            arg = dict((func.__name__, func) for func in arg)

        results = {}

        for name, func in arg.iteritems():
            try:
                result = func(self)
            except Exception:
                result = self.aggregate(func)
            results[name] = result

        return DataFrame(results)

    def _aggregate_simple(self, arg):
        values = self.obj.values
        result = {}
        for k, v in self.primary.indices.iteritems():
            result[k] = arg(values.take(v))

        return result

    def _aggregate_named(self, arg):
        result = {}

        for name in self.primary:
            grp = self.get_group(name)
            grp.name = name
            output = arg(grp)
            result[name] = output

        return result

    def transform(self, func):
        """
        For given Series, group index by given mapper function or dict, take
        the sub-Series (reindex) for this group and call apply(applyfunc)
        on this sub-Series. Return a Series of the results for each
        key.

        Parameters
        ----------
        mapper : function
            on being called on each element of the Series
            index, determines the groups.

        applyfunc : function to apply to each group

        Note
        ----
        This function does not aggregate like groupby/tgroupby,
        the results of the given function on the subSeries should be another
        Series.

        Example
        -------
        series.transform(lambda x: mapping[x],
                         lambda x: (x - x.mean()) / x.std())

        Returns
        -------
        Series standardized by each unique value of mapping
        """
        result = self.obj.copy()

        for name, group in self:
            # XXX
            group.name = name
            res = func(group)
            indexer, _ = self.obj.index.get_indexer(group.index)
            np.put(result, indexer, res)

        return result

def _ravel_names(axes, shape):
    """
    Compute labeling vector for raveled values vector
    """
    unrolled = []
    for i, ax in enumerate(axes):
        tile_shape = shape[:i] + shape[i+1:] + (1,)
        tiled = np.tile(ax, tile_shape)
        tiled = tiled.swapaxes(i, -1)
        unrolled.append(tiled.ravel())

    return unrolled

class DataFrameGroupBy(GroupBy):

    def get_group(self, name, obj=None):
        # faster get_group for Series
        if obj is None:
            obj = self.obj

        inds = self.primary.indices[name]
        return obj.take(inds, axis=self.axis)

    @property
    def _agg_stride_shape(self):
        if self.axis == 0:
            n = len(self.obj.columns)
        else:
            n = len(self.obj.index)

        n -= len(self.exclusions)

        return n,

    def __getitem__(self, key):
        if key not in self.obj:
            raise KeyError('column %s not found' % key)
        return SeriesGroupBy(self.obj[key], groupings=self.groupings,
                             exclusions=self.exclusions, name=key)

    def _iterate_slices(self):
        if self.axis == 0:
            slice_axis = self.obj.columns
            slicer = lambda x: self.obj[x]
        else:
            slice_axis = self.obj.index
            slicer = self.obj.xs

        for val in slice_axis:
            if val in self.exclusions:
                continue

            yield val, slicer(val)

    @cache_readonly
    def _obj_with_exclusions(self):
        if len(self.exclusions) > 0:
            return self.obj.drop(self.exclusions, axis=1)
        else:
            return self.obj

    def aggregate(self, arg):
        """
        Aggregate using input function or dict of {column -> function}

        Parameters
        ----------
        arg : function or dict

            Function to use for aggregating groups. If a function, must either
            work when passed a DataFrame or when passed to DataFrame.apply. If
            pass a dict, the keys must be DataFrame column names

        Returns
        -------
        aggregated : DataFrame
        """
        if isinstance(arg, basestring):
            return getattr(self, arg)()

        result = {}
        if isinstance(arg, dict):
            for col, func in arg.iteritems():
                result[col] = self[col].agg(func)

            result = DataFrame(result)
        else:
            if len(self.groupings) > 1:
                return self._python_agg_general(arg)
            result = self._aggregate_generic(arg, axis=self.axis)

        return result

    def _aggregate_generic(self, agger, axis=0):
        result = {}

        obj = self._obj_with_exclusions

        try:
            for name in self.primary:
                data = self.get_group(name, obj=obj)
                try:
                    result[name] = agger(data)
                except Exception:
                    result[name] = data.apply(agger, axis=axis)
        except Exception, e1:
            if axis == 0:
                try:
                    return self._aggregate_item_by_item(agger)
                except Exception:
                    raise e1
            else:
                raise e1

        result = DataFrame(result)
        if axis == 0:
            result = result.T

        return result

    def _aggregate_item_by_item(self, agger):
        # only for axis==0

        obj = self._obj_with_exclusions

        result = {}
        cannot_agg = []
        for item in obj:
            try:
                result[item] = self[item].agg(agger)
            except (ValueError, TypeError):
                cannot_agg.append(item)
                continue

        return DataFrame(result)

    def _wrap_aggregated_output(self, output, mask):
        if len(self.groupings) > 1:
            index = self._get_multi_index(mask)
            result = DataFrame(output, index=index)
        else:
            name_list = self._get_names()
            result = DataFrame(output, index=name_list[0][1])

        if self.axis == 1:
            result = result.T

        return result

    def _wrap_applied_output(self, keys, values, not_indexed_same=False):
        if len(keys) == 0:
            # XXX
            return DataFrame({})

        if isinstance(values[0], DataFrame):
            return _concat_frames(values)
        else:
            if len(self.groupings) > 1:
                keys = MultiIndex.from_tuples(keys)

            # obj = self._obj_with_exclusions

            if self.axis == 0:
                stacked_values = np.vstack(values)
                columns = values[0].index
                index = keys
            else:
                stacked_values = np.vstack(values)
                index = values[0].index
                columns = keys

            return DataFrame(stacked_values, index=index, columns=columns)

    def transform(self, func):
        """
        For given DataFrame, group index by given mapper function or dict, take
        the sub-DataFrame (reindex) for this group and call apply(func)
        on this sub-DataFrame. Return a DataFrame of the results for each
        key.

        Note: this function does not aggregate like groupby/tgroupby,
        the results of the given function on the subDataFrame should be another
        DataFrame.

        Parameters
        ----------
        mapper : function, dict-like, or string
            Mapping or mapping function. If string given, must be a column
            name in the frame
        func : function
            Function to apply to each subframe

        Note
        ----
        Each subframe is endowed the attribute 'name' in case you need to know
        which group you are working on.

        Example
        --------
        >>> grouped = df.groupby(lambda x: mapping[x])
        >>> grouped.transform(lambda x: (x - x.mean()) / x.std())
        """
        applied = []

        obj = self._obj_with_exclusions
        for name, group in self:
            group.name = name

            try:
                res = group.apply(func, axis=self.axis)
            except Exception: # pragma: no cover
                res = func(group)

            # broadcasting
            if isinstance(res, Series):
                if res.index is obj.index:
                    group.T.values[:] = res
                else:
                    group.values[:] = res

                applied.append(group)
            else:
                applied.append(res)

        return _concat_frames(applied, obj.index, obj.columns,
                              axis=self.axis)

def _concat_frames(frames, index=None, columns=None, axis=0):
    if axis == 0:
        all_index = [np.asarray(x.index) for x in frames]
        new_index = Index(np.concatenate(all_index))

        if columns is None:
            new_columns = frames[0].columns
        else:
            new_columns = columns
    else:
        all_columns = [np.asarray(x.columns) for x in frames]
        new_columns = Index(np.concatenate(all_columns))
        if index is None:
            new_index = frames[0].index
        else:
            new_index = index

    new_values = np.concatenate([x.values for x in frames], axis=axis)
    result = DataFrame(new_values, index=new_index, columns=new_columns)
    return result.reindex(index=index, columns=columns)

class WidePanelGroupBy(GroupBy):

    def aggregate(self, func):
        """
        For given DataFrame, group index by given mapper function or dict, take
        the sub-DataFrame (reindex) for this group and call apply(func)
        on this sub-DataFrame. Return a DataFrame of the results for each
        key.

        Parameters
        ----------
        mapper : function, dict-like, or string
            Mapping or mapping function. If string given, must be a column
            name in the frame
        func : function
            Function to use for aggregating groups

        N.B.: func must produce one value from a Series, otherwise
        an error will occur.

        Optional: provide set mapping as dictionary
        """
        return self._aggregate_generic(func, axis=self.axis)

    def _aggregate_generic(self, agger, axis=0):
        result = {}

        obj = self._obj_with_exclusions

        for name in self.primary:
            data = self.get_group(name, obj=obj)
            try:
                result[name] = agger(data)
            except Exception:
                result[name] = data.apply(agger, axis=axis)

        result = WidePanel.fromDict(result, intersect=False)

        if axis > 0:
            result = result.swapaxes(0, axis)

        return result


class NDArrayGroupBy(GroupBy):
    pass

#-------------------------------------------------------------------------------
# Grouping generator for BlockManager

def generate_groups(data, label_list, shape, axis=0, factory=lambda x: x):
    """
    Parameters
    ----------
    data : BlockManager

    Returns
    -------
    generator
    """
    sorted_data, sorted_labels = _group_reorder(data, label_list, axis=axis)

    gen = _generate_groups(sorted_data, sorted_labels, shape,
                           0, len(label_list[0]), axis=axis, which=0,
                           factory=factory)
    for key, group in gen:
        yield key, group

def _group_reorder(data, label_list, axis=0):
    indexer = np.lexsort(label_list[::-1])
    sorted_labels = [labels.take(indexer) for labels in label_list]

    if isinstance(data, BlockManager):
        # this is sort of wasteful but...
        sorted_axis = data.axes[axis].take(indexer)
        sorted_data = data.reindex_axis(sorted_axis, axis=axis)
    if isinstance(data, Series):
        sorted_axis = data.index.take(indexer)
        sorted_data = data.reindex(sorted_axis)
    elif isinstance(data, DataFrame):
        sorted_data = data.take(indexer, axis=axis)

    return sorted_data, sorted_labels

def _generate_groups(data, labels, shape, start, end, axis=0, which=0,
                     factory=lambda x: x):
    axis_labels = labels[which][start:end]
    edges = axis_labels.searchsorted(np.arange(1, shape[which] + 1),
                                     side='left')

    if isinstance(data, DataFrame):
        def slicer(data, slob):
            if axis == 0:
                return data[slob]
            else:
                return data.ix[:, slob]
    elif isinstance(data, BlockManager):
        def slicer(data, slob):
            return factory(data.get_slice(slob, axis=axis))
    else:
        def slicer(data, slob):
            return data[slob]

    do_slice = which == len(labels) - 1

    # omit -1 values at beginning-- NA values
    left = axis_labels.searchsorted(0)

    # time to actually aggregate
    for i, right in enumerate(edges):
        if do_slice:
            slob = slice(start + left, start + right)
            yield i, slicer(data, slob)
        else:
            # yield subgenerators, yikes
            yield i, _generate_groups(data, labels, shape, start + left,
                                      start + right, axis=axis,
                                      which=which + 1, factory=factory)

        left = right

#----------------------------------------------------------------------
# sorting levels...cleverly?

def sort_group_labels(ids, labels, counts):
    n = len(ids)
    rng = np.arange(n)
    values = Series(ids, index=rng, dtype=object).values
    indexer = values.argsort()

    reverse_indexer = np.empty(n, dtype=np.int32)
    reverse_indexer.put(indexer, np.arange(n))

    new_labels = reverse_indexer.take(labels)
    np.putmask(new_labels, labels == -1, -1)

    new_ids = dict(izip(rng, values.take(indexer)))
    new_counts = counts.take(indexer)

    return new_ids, new_labels, new_counts
