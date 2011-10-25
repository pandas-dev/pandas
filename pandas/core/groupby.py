from itertools import izip
import types

import numpy as np

from pandas.core.frame import DataFrame
from pandas.core.generic import NDFrame, PandasObject
from pandas.core.index import Index, Int64Index, MultiIndex
from pandas.core.internals import BlockManager
from pandas.core.series import Series
from pandas.core.panel import Panel
from pandas.util.decorators import cache_readonly
import pandas._tseries as lib


class GroupBy(object):
    """
    Class for grouping and aggregating relational data. See aggregate,
    transform, and apply functions on this object.

    It's easiest to use obj.groupby(...) to use GroupBy, but you can also do:

    ::

        grouped = groupby(obj, ...)

    Parameters
    ----------
    obj : pandas object
    axis : int, default 0
    level : int, default None
        Level of MultiIndex
    groupings : list of Grouping objects
        Most users should ignore this
    exclusions : array-like, optional
        List of columns to exclude
    name : string
        Most users should ignore this

    Notes
    -----
    After grouping, see aggregate, apply, and transform functions. Here are some
    other brief notes about usage. When grouping by multiple groups, the result
    index will be a MultiIndex (hierarhical) by default.

    Iteration produces (key, group) tuples, i.e. chunking the data by group. So
    you can write code like:

    ::

        grouped = obj.groupby(grouper, axis=axis)
        for key, group in grouped:
            # do something with the data

    Function calls on GroupBy, if not specially implemented, "dispatch" to the
    grouped data. So if you group a DataFrame and wish to invoke the std()
    method on each group, you can simply do:

    ::

        df.groupby(mapper).std()

    rather than

    ::

        df.groupby(mapper).aggregate(np.std)

    You can pass arguments to these "wrapped" functions, too.

    See the online documentation for full exposition on these topics and much
    more

    Returns
    -------
    **Attributes**
    groups : dict
        {group name -> group labels}
    len(grouped) : int
        Number of groups
    """

    def __init__(self, obj, grouper=None, axis=0, level=None,
                 groupings=None, exclusions=None, name=None, as_index=True):
        self._name = name
        self.obj = obj
        self.axis = axis
        self.level = level

        if not as_index:
            if not isinstance(obj, DataFrame):
                raise TypeError('as_index=False only valid with DataFrame')
            if axis != 0:
                raise ValueError('as_index=False only valid for axis=0')

        self.as_index = as_index

        if groupings is None:
            groupings, exclusions = _get_groupings(obj, grouper, axis=axis,
                                                   level=level)

        self.groupings = groupings
        self.exclusions = set(exclusions)

    def __len__(self):
        return len(self.groups)

    _groups = None
    @property
    def groups(self):
        if self._groups is not None:
            return self._groups

        if len(self.groupings) == 1:
            self._groups = self.primary.groups
        else:
            to_groupby = zip(*(ping.grouper for ping in self.groupings))
            to_groupby = Index(to_groupby)

            axis = self.obj._get_axis(self.axis)
            self._groups = axis.groupby(to_groupby)

        return self._groups

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

    @cache_readonly
    def use_take(self):
        group_axis = self.obj._get_axis(self.axis)
        return isinstance(group_axis, Int64Index) and len(self.groupings) == 1

    def get_group(self, name, obj=None):
        if obj is None:
            obj = self.obj

        if self.use_take:
            inds = self.primary.indices[name]
            return obj.take(inds, axis=self.axis)
        else:
            labels = self.groups[name]
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
            groups = self.primary.groups.keys()
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

        def flatten(gen, level=0, shape_axis=0):
            ids = self.groupings[level].ids
            for cat, subgen in gen:
                if isinstance(subgen, tipo):
                    yield (ids[cat],), subgen
                else:
                    for subcat, data in flatten(subgen, level=level+1,
                                                shape_axis=shape_axis):
                        if data.shape[shape_axis] == 0:
                            continue
                        yield (ids[cat],) + subcat, data

        return flatten(self._generator_factory(data), shape_axis=self.axis)

    def apply(self, func, *args, **kwargs):
        """
        Apply function and combine results together in an intelligent way. The
        split-apply-combine combination rules attempt to be as common sense
        based as possible. For example:

        case 1:
        group DataFrame
        apply aggregation function (f(chunk) -> Series)
        yield DataFrame, with group axis having group labels

        case 2:
        group DataFrame
        apply transform function ((f(chunk) -> DataFrame with same indexes)
        yield DataFrame with resulting chunks glued together

        case 3:
        group Series
        apply function with f(chunk) -> DataFrame
        yield DataFrame with result of chunks glued together

        Parameters
        ----------
        func : function

        Notes
        -----
        See online documentation for full exposition on how to use apply

        See also
        --------
        aggregate, transform

        Returns
        -------
        applied : type depending on grouped object and function
        """
        return self._python_apply_general(func, *args, **kwargs)

    def aggregate(self, func, *args, **kwargs):
        raise NotImplementedError

    def agg(self, func, *args, **kwargs):
        """
        See docstring for aggregate
        """
        return self.aggregate(func, *args, **kwargs)

    def _get_names(self):
        axes = [ping.group_index for ping in self.groupings]
        grouping_names = [ping.name for ping in self.groupings]
        shape = self._group_shape
        return zip(grouping_names, _ravel_names(axes, shape))

    def _iterate_slices(self):
        yield self.name, self.obj

    def transform(self, func, *args, **kwargs):
        raise NotImplementedError

    def mean(self):
        """
        Compute mean of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex
        """
        return self._cython_agg_general('mean')

    def size(self):
        """
        Compute group sizes
        """
        result = sorted((k, len(v)) for k, v in self.groups.iteritems())
        keys, values = zip(*result)

        if len(self.groupings) > 1:
            index = MultiIndex.from_tuples(keys)
        else:
            index = Index(keys)

        return Series(values, index=index)

    def sum(self):
        """
        Compute sum of values, excluding missing values

        For multiple groupings, the result index will be a MultiIndex
        """
        try:
            return self._cython_agg_general('add')
        except Exception:
            return self.aggregate(lambda x: np.sum(x, axis=self.axis))

    def _cython_agg_general(self, how):
        label_list = [ping.labels for ping in self.groupings]
        shape = self._group_shape

        # TODO: address inefficiencies, like duplicating effort (should
        # aggregate all the columns at once?)

        output = {}
        cannot_agg = []
        for name, obj in self._iterate_slices():
            if issubclass(obj.dtype.type, np.number):
                if obj.dtype != np.float64:
                    obj = obj.astype('f8')
            else:
                cannot_agg.append(name)
                continue

            result, counts =  lib.group_aggregate(obj, label_list,
                                                  shape, how=how)
            result = result.ravel()
            mask = counts.ravel() > 0
            output[name] = result[mask]

        return self._wrap_aggregated_output(output, mask)

    def _get_multi_index(self, mask):
        masked = [labels for _, labels in self._get_group_levels(mask)]
        names = [ping.name for ping in self.groupings]
        return MultiIndex.from_arrays(masked, names=names)

    def _get_group_levels(self, mask):
        name_list = self._get_names()
        return [(name, raveled[mask]) for name, raveled in name_list]

    def _python_agg_general(self, func, *args, **kwargs):
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
                    reschunk[i] = func(subgen, *args, **kwargs)
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
                output[name] = result.ravel().copy()
                result.fill(np.nan)

            mask = counts.ravel() > 0
            for name, result in output.iteritems():
                output[name] = result[mask]

        return self._wrap_aggregated_output(output, mask)

    def _python_apply_general(self, func, *args, **kwargs):
        result_keys = []
        result_values = []

        not_indexed_same = False
        for key, group in self:
            group.name = key

            # group might be modified
            group_axes = _get_axes(group)

            res = func(group, *args, **kwargs)

            if not _is_indexed_like(res, group_axes):
                not_indexed_same = True

            result_keys.append(key)
            result_values.append(res)

        return self._wrap_applied_output(result_keys, result_values,
                                         not_indexed_same=not_indexed_same)

    def _wrap_applied_output(self, *args, **kwargs):
        raise NotImplementedError

    def _wrap_frames(self, keys, values, not_indexed_same=False):
        if not_indexed_same:
            result = _concat_frames_hierarchical(values, keys,
                                                 self.groupings,
                                                 axis=self.axis)
        else:
            result = _concat_frames(values, self.obj.index)

        return result

    @property
    def _generator_factory(self):
        labels = [ping.labels for ping in self.groupings]
        shape = self._group_shape

        if isinstance(self.obj, NDFrame):
            factory = self.obj._constructor
        else:
            factory = None

        return lambda obj: generate_groups(obj, labels, shape,
                                           axis=self.axis,
                                           factory=factory)

def groupby(obj, by, **kwds):
    if isinstance(obj, Series):
        klass = SeriesGroupBy
    elif isinstance(obj, DataFrame):
        klass = DataFrameGroupBy
    else: # pragma: no cover
        raise TypeError('invalid type: %s' % type(obj))

    return klass(obj, by, **kwds)
groupby.__doc__ = GroupBy.__doc__

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

class Grouping(object):
    """
    Holds the grouping information for a single key

    Parameters
    ----------
    index : Index
    grouper :
    name :
    level :

    Returns
    -------
    **Attributes**:
      * indices : dict of {group -> index_list}
      * labels : ndarray, group labels
      * ids : mapping of label -> group
      * reverse_ids : mapping of group -> label
      * counts : array of group counts
      * group_index : unique groups
      * groups : dict of {group -> label_list}
    """
    def __init__(self, index, grouper=None, name=None, level=None):
        self.name = name
        self.level = level
        self.grouper = _convert_grouper(index, grouper)
        self.index = index

        if level is not None:
            if not isinstance(level, int):
                assert(level in index.names)
                level = index.names.index(level)

            inds = index.labels[level]
            labels = index.levels[level].take(inds)

            if grouper is not None:
                self.grouper = labels.map(self.grouper)
            else:
                self.grouper = labels

        # no level passed
        if not isinstance(self.grouper, np.ndarray):
            self.grouper = self.index.map(self.grouper)

    def __repr__(self):
        return 'Grouping(%s)' % self.name

    def __iter__(self):
        return iter(self.indices)

    _labels = None
    _ids = None
    _counts = None

    @cache_readonly
    def indices(self):
        return lib.groupby_indices(self.grouper)

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

    @cache_readonly
    def reverse_ids(self):
        return dict((v, k) for k, v in self.ids.iteritems())

    @property
    def counts(self):
        if self._counts is None:
            self._make_labels()
        return self._counts

    @cache_readonly
    def group_index(self):
        return Index([self.ids[i] for i in range(len(self.ids))])

    def _make_labels(self):
        ids, labels, counts  = lib.group_labels(self.grouper)
        sids, slabels, scounts = sort_group_labels(ids, labels, counts)
        self._labels = slabels
        self._ids = sids
        self._counts = scounts

    _groups = None
    @property
    def groups(self):
        if self._groups is None:
            self._groups = self.index.groupby(self.grouper)
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
            return grouper.reindex(axis).astype(object)
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

    def aggregate(self, func_or_funcs, *args, **kwargs):
        """
        Apply aggregation function or functions to groups, yielding most likely
        Series but in some cases DataFrame depending on the output of the
        aggregation function

        Parameters
        ----------
        func_or_funcs : function or list / dict of functions
            List/dict of functions will produce DataFrame with column names
            determined by the function names themselves (list) or the keys in
            the dict

        Notes
        -----
        agg is an alias for aggregate. Use it.

        Example
        -------
        >>> series
        bar    1.0
        baz    2.0
        qot    3.0
        qux    4.0

        >>> mapper = lambda x: x[0] # first letter
        >>> grouped = series.groupby(mapper)

        >>> grouped.aggregate(np.sum)
        b    3.0
        q    7.0

        >>> grouped.aggregate([np.sum, np.mean, np.std])
           mean  std  sum
        b  1.5   0.5  3
        q  3.5   0.5  7

        >>> grouped.agg({'result' : lambda x: x.mean() / x.std(),
        ...              'total' : np.sum})
           result  total
        b  2.121   3
        q  4.95    7

        See also
        --------
        apply, transform

        Returns
        -------
        Series or DataFrame
        """
        if isinstance(func_or_funcs, basestring):
            return getattr(self, func_or_funcs)(*args, **kwargs)

        if hasattr(func_or_funcs,'__iter__'):
            ret = self._aggregate_multiple_funcs(func_or_funcs)
        else:
            if len(self.groupings) > 1:
                return self._python_agg_general(func_or_funcs, *args, **kwargs)

            try:
                result = self._aggregate_simple(func_or_funcs, *args, **kwargs)
            except Exception:
                result = self._aggregate_named(func_or_funcs, *args, **kwargs)

            ret = Series(result)

        if not self.as_index:  # pragma: no cover
            print 'Warning, ignoring as_index=True'

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

        key_names = [ping.name for ping in self.groupings]

        if isinstance(values[0], Series):
            if not_indexed_same:
                data_dict = dict(zip(keys, values))
                result = DataFrame(data_dict).T
                if len(self.groupings) > 1:
                    result.index = MultiIndex.from_tuples(keys, names=key_names)
                return result
            else:
                cat_values = np.concatenate([x.values for x in values])
                cat_index = values[0].index
                if len(values) > 1:
                    cat_index = cat_index.append([x.index for x in values[1:]])
                return Series(cat_values, index=cat_index)
        elif isinstance(values[0], DataFrame):
            # possible that Series -> DataFrame by applied function
            return self._wrap_frames(keys, values,
                                     not_indexed_same=not_indexed_same)
        else:
            if len(self.groupings) > 1:
                index = MultiIndex.from_tuples(keys, names=key_names)
                return Series(values, index)
            else:
                return Series(values, keys)

    def _aggregate_multiple_funcs(self, arg):
        if not isinstance(arg, dict):
            arg = dict((func.__name__, func) for func in arg)

        results = {}

        for name, func in arg.iteritems():
            results[name] = self.aggregate(func)

        return DataFrame(results)

    def _aggregate_simple(self, func, *args, **kwargs):
        values = self.obj.values
        result = {}
        for k, v in self.primary.indices.iteritems():
            agged = func(values.take(v), *args, **kwargs)
            if isinstance(agged, np.ndarray):
                raise Exception('Must produce aggregated value')
            result[k] = agged

        return result

    def _aggregate_named(self, func, *args, **kwargs):
        result = {}

        for name in self.primary:
            grp = self.get_group(name)
            grp.name = name
            output = func(grp, *args, **kwargs)
            if isinstance(output, np.ndarray):
                raise Exception('Must produce aggregated value')
            result[name] = output

        return result

    def transform(self, func, *args, **kwargs):
        """
        Call function producing a like-indexed Series on each group and return
        a Series with the transformed values

        Parameters
        ----------
        func : function
            To apply to each group. Should return a Series with the same index

        Example
        -------
        >>> grouped.transform(lambda x: (x - x.mean()) / x.std())

        Returns
        -------
        transformed : Series
        """
        result = self.obj.copy()

        for name, group in self:
            group.name = name
            res = func(group, *args, **kwargs)
            indexer = self.obj.index.get_indexer(group.index)
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

    @property
    def _agg_stride_shape(self):
        if self.axis == 0:
            n = len(self.obj.columns)
        else:
            n = len(self.obj.index)

        n -= len(self.exclusions)

        return n,

    def __getitem__(self, key):
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

    def aggregate(self, arg, *args, **kwargs):
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
            if self.axis != 0:  # pragma: no cover
                raise ValueError('Can only pass dict with axis=0')

            for col, func in arg.iteritems():
                result[col] = self[col].agg(func)

            result = DataFrame(result)
        else:
            if len(self.groupings) > 1:
                try:
                    return self._python_agg_general(arg, *args, **kwargs)
                except Exception:
                    return self._aggregate_item_by_item(arg, *args, **kwargs)
            result = self._aggregate_generic(arg, *args, **kwargs)

        if not self.as_index:
            if isinstance(result.index, MultiIndex):
                zipped = zip(result.index.levels, result.index.labels,
                             result.index.names)
                for i, (lev, lab, name) in enumerate(zipped):
                    result.insert(i, name, lev.values.take(lab))
                result = result.consolidate()
            else:
                values = result.index.values
                name = self.groupings[0].name
                result.insert(0, name, values)
            result.index = np.arange(len(result))

        return result

    def _aggregate_generic(self, func, *args, **kwargs):
        result = {}
        axis = self.axis
        obj = self._obj_with_exclusions

        try:
            for name in self.primary:
                data = self.get_group(name, obj=obj)
                try:
                    result[name] = func(data, *args, **kwargs)
                except Exception:
                    wrapper = lambda x: func(x, *args, **kwargs)
                    result[name] = data.apply(wrapper, axis=axis)
        except Exception, e1:
            if axis == 0:
                return self._aggregate_item_by_item(func, *args, **kwargs)
            else:
                raise e1

        if result:
            if axis == 0:
                result = DataFrame(result, index=obj.columns).T
            else:
                result = DataFrame(result, index=obj.index)
        else:
            result = DataFrame(result)

        return result

    def _aggregate_item_by_item(self, func, *args, **kwargs):
        # only for axis==0

        obj = self._obj_with_exclusions

        result = {}
        cannot_agg = []
        for item in obj:
            try:
                result[item] = self[item].agg(func, *args, **kwargs)
            except (ValueError, TypeError):
                cannot_agg.append(item)
                continue

        return DataFrame(result)

    def _wrap_aggregated_output(self, output, mask):
        agg_axis = 0 if self.axis == 1 else 1
        agg_labels = self.obj._get_axis(agg_axis)
        if isinstance(output, dict):
            if len(output) == len(agg_labels):
                output_keys = agg_labels
            else:
                output_keys = sorted(output)
                try:
                    output_keys.sort()
                except Exception:  # pragma: no cover
                    pass

                if isinstance(agg_labels, MultiIndex):
                    output_keys = MultiIndex.from_tuples(output_keys,
                                                         names=agg_labels.names)
        else:
            output_keys = agg_labels

        if len(self.groupings) > 1:
            if not self.as_index:
                result = DataFrame(output, columns=output_keys)
                group_levels = self._get_group_levels(mask)
                for i, (name, labels) in enumerate(group_levels):
                    result.insert(i, name, labels)
                result = result.consolidate()
            else:
                index = self._get_multi_index(mask)
                result = DataFrame(output, index=index, columns=output_keys)
        else:
            name_list = self._get_names()
            name, labels = name_list[0]
            if not self.as_index:
                result = DataFrame(output, columns=output_keys)
                result.insert(0, name, labels)
            else:
                result = DataFrame(output, index=labels, columns=output_keys)

        if self.axis == 1:
            result = result.T

        return result

    def _wrap_applied_output(self, keys, values, not_indexed_same=False):
        if len(keys) == 0:
            # XXX
            return DataFrame({})

        key_names = [ping.name for ping in self.groupings]

        if isinstance(values[0], DataFrame):
            return self._wrap_frames(keys, values,
                                     not_indexed_same=not_indexed_same)
        else:
            if len(self.groupings) > 1:
                keys = MultiIndex.from_tuples(keys, names=key_names)

            if isinstance(values[0], np.ndarray):
                if self.axis == 0:
                    stacked_values = np.vstack([np.asarray(x)
                                                for x in values])
                    columns = values[0].index
                    index = keys
                else:
                    stacked_values = np.vstack([np.asarray(x)
                                                for x in values]).T
                    index = values[0].index
                    columns = keys
                return DataFrame(stacked_values, index=index,
                                 columns=columns)
            else:
                return Series(values, index=keys)

    def transform(self, func, *args, **kwargs):
        """
        Call function producing a like-indexed DataFrame on each group and
        return a DataFrame having the same indexes as the original object
        filled with the transformed values

        Parameters
        ----------
        f : function
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
                wrapper = lambda x: func(x, *args, **kwargs)
                res = group.apply(wrapper, axis=self.axis)
            except Exception: # pragma: no cover
                res = func(group, *args, **kwargs)

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

def _concat_frames(frames, index, columns=None, axis=0):
    if len(frames) == 1:
        return frames[0]

    if axis == 0:
        new_index = _concat_indexes([x.index for x in frames])
        if columns is None:
            new_columns = frames[0].columns
        else:
            new_columns = columns
    else:
        new_columns = _concat_indexes([x.columns for x in frames])
        new_index = index

    if frames[0]._is_mixed_type:
        new_data = {}
        for col in new_columns:
            new_data[col] = np.concatenate([x[col].values for x in frames])
        return DataFrame(new_data, index=new_index, columns=new_columns)
    else:
        new_values = np.concatenate([x.values for x in frames], axis=axis)
        result = DataFrame(new_values, index=new_index, columns=new_columns)
        return result.reindex(index=index, columns=columns)

def _concat_indexes(indexes):
    return indexes[0].append(indexes[1:])

def _concat_frames_hierarchical(frames, keys, groupings, axis=0):
    if axis == 0:
        indexes = [x.index for x in frames]
        new_index = _make_concat_multiindex(indexes, keys, groupings)
        new_columns = frames[0].columns
    else:
        all_columns = [x.columns for x in frames]
        new_columns = _make_concat_multiindex(all_columns, keys, groupings)
        new_index = frames[0].index

    if frames[0]._is_mixed_type:
        new_data = {}
        for col in new_columns:
            new_data[col] = np.concatenate([x[col].values for x in frames])
        return DataFrame(new_data, index=new_index, columns=new_columns)
    else:
        new_values = np.concatenate([x.values for x in frames], axis=axis)
        return DataFrame(new_values, index=new_index, columns=new_columns)

def _make_concat_multiindex(indexes, keys, groupings):
    if not _all_indexes_same(indexes):
        label_list = []

        # things are potentially different sizes, so compute the exact labels
        # for each level and pass those to MultiIndex.from_arrays
        if len(groupings) == 1:
            zipped = [keys]
        else:
            zipped = zip(*keys)

        for hlevel in zipped:
            to_concat = []
            for k, index in zip(hlevel, indexes):
                to_concat.append(np.repeat(k, len(index)))
            label_list.append(np.concatenate(to_concat))

        concat_index = _concat_indexes(indexes)

        # these go at the end
        if isinstance(concat_index, MultiIndex):
            for level in range(concat_index.nlevels):
                label_list.append(concat_index.get_level_values(level))
        else:
            label_list.append(concat_index.values)

        return MultiIndex.from_arrays(label_list)

    new_index = indexes[0]
    n = len(new_index)

    # do something a bit more speedy
    levels = [ping.group_index for ping in  groupings]
    levels.append(new_index)

    # construct labels
    labels = []

    for hlevel, ping in zip(zip(*keys), groupings):
        get_id = ping.reverse_ids.__getitem__
        mapped = [get_id(x) for x in hlevel]
        labels.append(np.repeat(mapped, n))

    # last labels for the new level
    labels.append(np.tile(np.arange(n), len(indexes)))
    return MultiIndex(levels=levels, labels=labels)

def _all_indexes_same(indexes):
    first = indexes[0]
    for index in indexes[1:]:
        if not first.equals(index):
            return False
    return True

class PanelGroupBy(GroupBy):

    def aggregate(self, func, *args, **kwargs):
        """
        Aggregate using input function or dict of {column -> function}

        Parameters
        ----------
        arg : function or dict
            Function to use for aggregating groups. If a function, must either
            work when passed a Panel or when passed to Panel.apply. If
            pass a dict, the keys must be DataFrame column names

        Returns
        -------
        aggregated : Panel
        """
        return self._aggregate_generic(func, *args, **kwargs)

    def _aggregate_generic(self, func, *args, **kwargs):
        result = {}

        axis = self.axis

        obj = self._obj_with_exclusions

        for name in self.primary:
            data = self.get_group(name, obj=obj)
            try:
                result[name] = func(data, *args, **kwargs)
            except Exception:
                wrapper = lambda x: func(x, *args, **kwargs)
                result[name] = data.apply(wrapper, axis=axis)

        result = Panel.fromDict(result, intersect=False)

        if axis > 0:
            result = result.swapaxes(0, axis)

        return result


class NDArrayGroupBy(GroupBy):
    pass

#----------------------------------------------------------------------
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
