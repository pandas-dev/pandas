from itertools import izip
import types

import numpy as np

from pandas.core.frame import DataFrame
from pandas.core.generic import NDFrame
from pandas.core.index import Index, MultiIndex
from pandas.core.internals import BlockManager, make_block
from pandas.core.series import Series
from pandas.core.panel import Panel
from pandas.util.decorators import cache_readonly, Appender
import pandas.core.algorithms as algos
import pandas.core.common as com
import pandas._tseries as lib


class GroupByError(Exception):
    pass


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
                 groupings=None, exclusions=None, column=None, as_index=True,
                 sort=True):
        self._column = column

        if isinstance(obj, NDFrame):
            obj._consolidate_inplace()

        self.obj = obj


        self.axis = axis
        self.level = level

        if not as_index:
            if not isinstance(obj, DataFrame):
                raise TypeError('as_index=False only valid with DataFrame')
            if axis != 0:
                raise ValueError('as_index=False only valid for axis=0')

        self.as_index = as_index
        self.grouper = grouper
        self.sort = sort

        if groupings is None:
            groupings, exclusions = _get_groupings(obj, grouper, axis=axis,
                                                   level=level, sort=sort)

        self.groupings = groupings
        self.exclusions = set(exclusions) if exclusions else set()

    def __len__(self):
        return len(self.indices)

    @cache_readonly
    def groups(self):
        if len(self.groupings) == 1:
            return self.primary.groups
        else:
            to_groupby = zip(*(ping.grouper for ping in self.groupings))
            to_groupby = Index(to_groupby)

            axis = self.obj._get_axis(self.axis)
            return axis.groupby(to_groupby)

    @cache_readonly
    def indices(self):
        if len(self.groupings) == 1:
            return self.primary.indices
        else:
            # TODO: this is massively inefficient
            to_groupby = zip(*(ping.grouper for ping in self.groupings))
            to_groupby = Index(to_groupby)
            return lib.groupby_indices(to_groupby)

    @property
    def name(self):
        if self._column is None:
            return 'result'
        else:
            return self._column

    @property
    def _obj_with_exclusions(self):
        return self.obj

    def __getattr__(self, attr):
        if hasattr(self.obj, attr) and attr != '_cache':
            return self._make_wrapper(attr)
        raise AttributeError("'%s' object has no attribute '%s'" %
                             (type(self).__name__, attr))

    def __getitem__(self, key):
        raise NotImplementedError

    def _make_wrapper(self, name):
        f = getattr(self.obj, name)
        if not isinstance(f, types.MethodType):
            return self.apply(lambda self: getattr(self, name))

        f = getattr(type(self.obj), name)

        def wrapper(*args, **kwargs):
            # a little trickery for aggregation functions that need an axis
            # argument
            kwargs_with_axis = kwargs.copy()
            if 'axis' not in kwargs_with_axis:
                kwargs_with_axis['axis'] = self.axis

            def curried_with_axis(x):
                return f(x, *args, **kwargs_with_axis)
            def curried(x):
                return f(x, *args, **kwargs)

            try:
                return self.apply(curried_with_axis)
            except Exception:
                return self.apply(curried)

        return wrapper

    @property
    def primary(self):
        return self.groupings[0]

    def get_group(self, name, obj=None):
        if obj is None:
            obj = self.obj

        inds = self.indices[name]
        return obj.take(inds, axis=self.axis)

    def __iter__(self):
        """
        Groupby iterator

        Returns
        -------
        Generator yielding sequence of (name, subsetted object)
        for each group
        """
        if len(self.groupings) == 1:
            groups = self.indices.keys()
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
        data = self.obj

        comp_ids, _, ngroups = self._group_info
        label_list = [ping.labels for ping in self.groupings]
        level_list = [ping.group_index for ping in self.groupings]
        mapper = _KeyMapper(comp_ids, ngroups, label_list, level_list)

        for label, group in self._generate_groups(data, comp_ids, ngroups):
            key = mapper.get_key(label)
            yield key, group

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

    def _iterate_slices(self):
        yield self.name, self.obj

    def transform(self, func, *args, **kwargs):
        raise NotImplementedError

    def mean(self):
        """
        Compute mean of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex
        """
        try:
            return self._cython_agg_general('mean')
        except GroupByError:
            raise
        except Exception:  # pragma: no cover
            f = lambda x: x.mean(axis=self.axis)
            return self._python_agg_general(f)

    def std(self):
        """
        Compute mean of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex
        """
        return self._cython_agg_general('std')

    def size(self):
        """
        Compute group sizes
        """
        result = sorted((k, len(v)) for k, v in self.groups.iteritems())
        keys, values = zip(*result)

        if len(self.groupings) > 1:
            names = [ping.name for ping in self.groupings]
            index = MultiIndex.from_tuples(keys, names=names)
        else:
            index = Index(keys, name=self.groupings[0].name)

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
        # TODO: address inefficiencies, like duplicating effort (should
        # aggregate all the columns at once?)

        comp_ids, obs_group_ids, max_group = self._group_info

        output = {}
        for name, obj in self._iterate_slices():
            if not issubclass(obj.dtype.type, (np.number, np.bool_)):
                continue

            obj = com._ensure_float64(obj)
            result, counts = cython_aggregate(obj, comp_ids,
                                              max_group, how=how)
            mask = counts > 0
            output[name] = result[mask]

        if len(output) == 0:
            raise GroupByError('No numeric types to aggregate')

        return self._wrap_aggregated_output(output, mask, obs_group_ids)

    def _python_agg_general(self, func, *args, **kwargs):
        agg_func = lambda x: func(x, *args, **kwargs)

        comp_ids, obs_group_ids, max_group = self._group_info

        # iterate through "columns" ex exclusions to populate output dict
        output = {}
        for name, obj in self._iterate_slices():
            try:
                result, counts = self._aggregate_series(obj, agg_func,
                                                        comp_ids, max_group)
                output[name] = result
            except TypeError:
                continue

        if len(output) == 0:
            return self._python_apply_general(func, *args, **kwargs)

        mask = counts.ravel() > 0
        for name, result in output.iteritems():
            output[name] = result[mask]

        return self._wrap_aggregated_output(output, mask, obs_group_ids)

    @property
    def _group_info(self):
        if len(self.groupings) > 1:
            all_labels = [ping.labels for ping in self.groupings]
            group_index = get_group_index(all_labels, self._group_shape)
            comp_ids, obs_group_ids = _compress_group_index(group_index)
        else:
            ping = self.groupings[0]
            group_index = ping.labels

        comp_ids, obs_group_ids = _compress_group_index(group_index)
        ngroups = len(obs_group_ids)
        comp_ids = com._ensure_int32(comp_ids)
        return comp_ids, obs_group_ids, ngroups

    @property
    def _group_shape(self):
        return tuple(ping.ngroups for ping in self.groupings)

    def _get_multi_index(self, mask, obs_ids):
        masked = [labels for _, labels in
                  self._get_group_levels(mask, obs_ids)]
        names = [ping.name for ping in self.groupings]
        return MultiIndex.from_arrays(masked, names=names)

    def _get_group_levels(self, mask, obs_ids):
        recons_labels = decons_group_index(obs_ids, self._group_shape)

        name_list = []
        for ping, labels in zip(self.groupings, recons_labels):
            labels = com._ensure_platform_int(labels)
            name_list.append((ping.name, ping.group_index.take(labels)))

        return name_list

    def _aggregate_series(self, obj, func, group_index, ngroups):
        try:
            return self._aggregate_series_fast(obj, func, group_index, ngroups)
        except Exception:
            return self._aggregate_series_pure_python(obj, func, group_index,
                                                      ngroups)

    def _aggregate_series_fast(self, obj, func, group_index, ngroups):
        if obj.index._has_complex_internals:
            raise TypeError('Incompatible index for Cython grouper')

        # avoids object / Series creation overhead
        dummy = obj[:0].copy()
        indexer = lib.groupsort_indexer(group_index, ngroups)[0]
        obj = obj.take(indexer)
        group_index = group_index.take(indexer)
        grouper = lib.SeriesGrouper(obj, func, group_index, ngroups,
                                    dummy)
        result, counts = grouper.get_result()
        return result, counts

    def _aggregate_series_pure_python(self, obj, func, group_index, ngroups):
        counts = np.zeros(ngroups, dtype=int)
        result = None

        for label, group in self._generate_groups(obj, group_index, ngroups):
            res = func(group)
            if result is None:
                try:
                    assert(not isinstance(res, np.ndarray))
                    assert(not (isinstance(res, list) and
                                len(res) == len(self.dummy)))

                    result = np.empty(ngroups, dtype='O')
                except Exception:
                    raise ValueError('function does not reduce')

            counts[label] = group.shape[0]
            result[label] = res

        result = lib.maybe_convert_objects(result, try_float=0)
        return result, counts

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
        from pandas.tools.merge import concat

        if not_indexed_same:
            group_keys = keys
            group_levels = [ping.group_index for ping in self.groupings]
            group_names = [ping.name for ping in self.groupings]
            result = concat(values, axis=self.axis, keys=group_keys,
                            levels=group_levels, names=group_names)
        else:
            result = concat(values, axis=self.axis)
            ax = self.obj._get_axis(self.axis)
            result = result.reindex_axis(ax, axis=self.axis)

        return result

    def _generate_groups(self, obj, group_index, ngroups):
        if isinstance(obj, NDFrame) and not isinstance(obj, DataFrame):
            factory = obj._constructor
            obj = obj._data
        else:
            factory = None

        return generate_groups(obj, group_index, ngroups,
                               axis=self.axis, factory=factory)

@Appender(GroupBy.__doc__)
def groupby(obj, by, **kwds):
    if isinstance(obj, Series):
        klass = SeriesGroupBy
    elif isinstance(obj, DataFrame):
        klass = DataFrameGroupBy
    else: # pragma: no cover
        raise TypeError('invalid type: %s' % type(obj))

    return klass(obj, by, **kwds)

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
      * counts : array of group counts
      * group_index : unique groups
      * groups : dict of {group -> label_list}
    """
    def __init__(self, index, grouper=None, name=None, level=None,
                 sort=True):
        self.name = name
        self.level = level
        self.grouper = _convert_grouper(index, grouper)
        self.index = index
        self.sort = sort

        # right place for this?
        if isinstance(grouper, Series) and name is None:
            self.name = grouper.name

        # pre-computed
        self._was_factor = False

        if level is not None:
            if not isinstance(level, int):
                assert(level in index.names)
                level = index.names.index(level)

            inds = index.labels[level]
            level_index = index.levels[level]

            if self.name is None:
                self.name = index.names[level]

            # XXX complete hack

            level_values = index.levels[level].take(inds)
            if grouper is not None:
                self.grouper = level_values.map(self.grouper)
            else:
                self._was_factor = True
                self._labels = inds
                self._group_index = level_index
                self.grouper = level_values
        else:
            # no level passed
            if not isinstance(self.grouper, np.ndarray):
                self.grouper = self.index.map(self.grouper)

    def __repr__(self):
        return 'Grouping(%s)' % self.name

    def __iter__(self):
        return iter(self.indices)

    _labels = None
    _counts = None
    _group_index = None

    @property
    def ngroups(self):
        return len(self.group_index)

    @cache_readonly
    def indices(self):
        return _groupby_indices(self.grouper)

    @property
    def labels(self):
        if self._labels is None:
            self._make_labels()
        return self._labels

    @property
    def counts(self):
        if self._counts is None:
            if self._was_factor:
                self._counts = lib.group_count(self.labels, self.ngroups)
            else:
                self._make_labels()
        return self._counts

    @property
    def group_index(self):
        if self._group_index is None:
            self._make_labels()
        return self._group_index

    def _make_labels(self):
        if self._was_factor:  # pragma: no cover
            raise Exception('Should not call this method grouping by level')
        else:
            labs, uniques, counts = algos.factorize(self.grouper,
                                                    sort=self.sort)
            uniques = Index(uniques, name=self.name)
            self._labels = labs
            self._group_index = uniques
            self._counts = counts

    _groups = None
    @property
    def groups(self):
        if self._groups is None:
            self._groups = self.index.groupby(self.grouper)
        return self._groups


def _get_groupings(obj, grouper=None, axis=0, level=None, sort=True):
    group_axis = obj._get_axis(axis)

    if level is not None and not isinstance(group_axis, MultiIndex):
        raise ValueError('can only specify level with multi-level index')

    if not isinstance(grouper, (tuple, list)):
        groupers = [grouper]
    else:
        groupers = grouper

    # what are we after, exactly?
    match_axis_length = len(groupers) == len(group_axis)
    any_callable = any(callable(g) or isinstance(g, dict) for g in groupers)
    any_arraylike = any(isinstance(g, (list, tuple, np.ndarray))
                        for g in groupers)

    try:
        if isinstance(obj, DataFrame):
            all_in_columns = all(g in obj.columns for g in groupers)
        else:
            all_in_columns = False
    except Exception:
        all_in_columns = False

    if (not any_callable and not all_in_columns
        and not any_arraylike and match_axis_length
        and not level):
        groupers = [com._asarray_tuplesafe(groupers)]

    if isinstance(level, (tuple, list)):
        if grouper is None:
            groupers = [None] * len(level)
        levels = level
    else:
        levels = [level] * len(groupers)

    groupings = []
    exclusions = []
    for i, (gpr, level) in enumerate(zip(groupers, levels)):
        name = None
        if _is_label_like(gpr):
            exclusions.append(gpr)
            name = gpr
            gpr = obj[gpr]
        ping = Grouping(group_axis, gpr, name=name, level=level, sort=sort)
        if ping.name is None:
            ping.name = 'key_%d' % i
        groupings.append(ping)

    return groupings, exclusions

def _is_label_like(val):
    return isinstance(val, basestring) or np.isscalar(val)

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
        return grouper
        # return np.asarray(grouper, dtype=object)
    else:
        return grouper

class SeriesGroupBy(GroupBy):

    _cythonized_methods = set(['add', 'mean'])

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
                return self._python_agg_general(func_or_funcs, *args, **kwargs)
            except Exception:
                result = self._aggregate_named(func_or_funcs, *args, **kwargs)

            index = Index(sorted(result), name=self.groupings[0].name)
            ret = Series(result, index=index)

        if not self.as_index:  # pragma: no cover
            print 'Warning, ignoring as_index=True'

        return ret

    def _aggregate_multiple_funcs(self, arg):
        if not isinstance(arg, dict):
            arg = dict((func.__name__, func) for func in arg)

        results = {}

        for name, func in arg.iteritems():
            results[name] = self.aggregate(func)

        return DataFrame(results)

    def _wrap_aggregated_output(self, output, mask, comp_ids):
        # sort of a kludge
        output = output[self.name]
        index = self._get_multi_index(mask, comp_ids)
        return Series(output, index=index)

    def _wrap_applied_output(self, keys, values, not_indexed_same=False):
        if len(keys) == 0:
            return Series([])

        key_names = [ping.name for ping in self.groupings]

        def _get_index():
            if len(self.groupings) > 1:
                index = MultiIndex.from_tuples(keys, names=key_names)
            else:
                ping = self.groupings[0]
                if len(keys) == ping.ngroups:
                    index = ping.group_index
                    index.name = key_names[0]
                else:
                    index = Index(keys, name=key_names[0])
            return index

        if isinstance(values[0], Series):
            if not_indexed_same:
                data_dict = dict(zip(keys, values))
                result = DataFrame(data_dict).T
                result.index = _get_index()
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
            return Series(values, index=_get_index())

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

class DataFrameGroupBy(GroupBy):

    def __getitem__(self, key):
        if self._column is not None:
            raise Exception('Column %s already selected' % self._column)

        if key not in self.obj:  # pragma: no cover
            raise KeyError(str(key))

        # kind of a kludge
        if self.as_index:
            return SeriesGroupBy(self.obj[key], column=key,
                                 groupings=self.groupings,
                                 exclusions=self.exclusions)
        else:
            return DataFrameGroupBy(self.obj, self.grouper, column=key,
                                    groupings=self.groupings,
                                    exclusions=self.exclusions,
                                    as_index=self.as_index)

    def _iterate_slices(self):
        if self.axis == 0:
            # kludge
            if self._column is None:
                slice_axis = self.obj.columns
            else:
                slice_axis = [self._column]
            slicer = lambda x: self.obj[x]
        else:
            slice_axis = self.obj.index
            slicer = self.obj.xs

        for val in slice_axis:
            if val in self.exclusions:
                continue

            yield val, slicer(val)

    def _cython_agg_general(self, how):

        comp_ids, obs_group_ids, max_group = self._group_info

        obj = self._obj_with_exclusions
        if self.axis == 1:
            obj = obj.T

        new_blocks = []

        for block in obj._data.blocks:
            values = block.values.T
            if not issubclass(values.dtype.type, (np.number, np.bool_)):
                continue

            values = com._ensure_float64(values)
            result, counts = cython_aggregate(values, comp_ids,
                                              max_group, how=how)

            mask = counts > 0
            if len(mask) > 0:
                result = result[mask]
            newb = make_block(result.T, block.items, block.ref_items)
            new_blocks.append(newb)

        if len(new_blocks) == 0:
            raise GroupByError('No numeric types to aggregate')

        agg_axis = 0 if self.axis == 1 else 1
        agg_labels = self._obj_with_exclusions._get_axis(agg_axis)

        if sum(len(x.items) for x in new_blocks) == len(agg_labels):
            output_keys = agg_labels
        else:
            output_keys = []
            for b in new_blocks:
                output_keys.extend(b.items)
            try:
                output_keys.sort()
            except TypeError:  # pragma: no cover
                pass

            if isinstance(agg_labels, MultiIndex):
                output_keys = MultiIndex.from_tuples(output_keys,
                                                     names=agg_labels.names)

        if not self.as_index:
            index = np.arange(new_blocks[0].values.shape[1])
            mgr = BlockManager(new_blocks, [output_keys, index])
            result = DataFrame(mgr)
            group_levels = self._get_group_levels(mask, obs_group_ids)
            for i, (name, labels) in enumerate(group_levels):
                result.insert(i, name, labels)
            result = result.consolidate()
        else:
            index = self._get_multi_index(mask, obs_group_ids)
            mgr = BlockManager(new_blocks, [output_keys, index])
            result = DataFrame(mgr)

        if self.axis == 1:
            result = result.T

        return result

    @cache_readonly
    def _obj_with_exclusions(self):
        if self._column is not None:
            return self.obj.reindex(columns=[self._column])

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

            obj = self._obj_with_exclusions
            for col, func in arg.iteritems():
                colg = SeriesGroupBy(obj[col], column=col,
                                     groupings=self.groupings)
                result[col] = colg.agg(func)

            result = DataFrame(result)
        elif isinstance(arg, list):
            return self._aggregate_multiple_funcs(arg)
        else:
            if len(self.groupings) > 1:
                return self._python_agg_general(arg, *args, **kwargs)
            else:
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

    def _aggregate_multiple_funcs(self, arg):
        from pandas.tools.merge import concat

        if self.axis != 0:
            raise NotImplementedError

        obj = self._obj_with_exclusions

        results = []
        keys = []
        for col in obj:
            try:
                colg = SeriesGroupBy(obj[col], column=col,
                                     groupings=self.groupings)
                results.append(colg.agg(arg))
                keys.append(col)
            except TypeError:
                pass

        result = concat(results, keys=keys, axis=1)

        return result

    def _aggregate_generic(self, func, *args, **kwargs):
        assert(len(self.groupings) == 1)

        axis = self.axis
        obj = self._obj_with_exclusions

        result = {}
        if axis == 0:
            try:
                for name in self.indices:
                    data = self.get_group(name, obj=obj)
                    result[name] = func(data, *args, **kwargs)
            except Exception:
                return self._aggregate_item_by_item(func, *args, **kwargs)
        else:
            for name in self.indices:
                try:
                    data = self.get_group(name, obj=obj)
                    result[name] = func(data, *args, **kwargs)
                except Exception:
                    wrapper = lambda x: func(x, *args, **kwargs)
                    result[name] = data.apply(wrapper, axis=axis)

        result_index = self.groupings[0].group_index

        if result:
            if axis == 0:
                result = DataFrame(result, index=obj.columns,
                                   columns=result_index).T
            else:
                result = DataFrame(result, index=obj.index,
                                   columns=result_index)
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
                colg = SeriesGroupBy(obj[item], column=item,
                                     groupings=self.groupings)
                result[item] = colg.agg(func, *args, **kwargs)
            except (ValueError, TypeError):
                cannot_agg.append(item)
                continue

        result_columns = obj.columns
        if cannot_agg:
            result_columns = result_columns.drop(cannot_agg)

        return DataFrame(result, columns=result_columns)

    def _wrap_aggregated_output(self, output, mask, comp_ids):
        agg_axis = 0 if self.axis == 1 else 1
        agg_labels = self._obj_with_exclusions._get_axis(agg_axis)

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

        if not self.as_index:
            result = DataFrame(output, columns=output_keys)
            group_levels = self._get_group_levels(mask, comp_ids)
            for i, (name, labels) in enumerate(group_levels):
                result.insert(i, name, labels)
            result = result.consolidate()
        else:
            index = self._get_multi_index(mask, comp_ids)
            result = DataFrame(output, index=index, columns=output_keys)

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
                key_index = MultiIndex.from_tuples(keys, names=key_names)
            else:
                ping = self.groupings[0]
                if len(keys) == ping.ngroups:
                    key_index = ping.group_index
                    key_index.name = key_names[0]

                    key_lookup = Index(keys)
                    indexer = key_lookup.get_indexer(key_index)

                    # reorder the values
                    values = [values[i] for i in indexer]
                else:
                    key_index = Index(keys, name=key_names[0])

            if isinstance(values[0], np.ndarray):
                if self.axis == 0:
                    stacked_values = np.vstack([np.asarray(x)
                                                for x in values])
                    columns = values[0].index
                    index = key_index
                else:
                    stacked_values = np.vstack([np.asarray(x)
                                                for x in values]).T
                    index = values[0].index
                    columns = key_index
                return DataFrame(stacked_values, index=index,
                                 columns=columns)
            else:
                return Series(values, index=key_index)

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
        from pandas.tools.merge import concat

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

        concat_index = obj.columns if self.axis == 0 else obj.index
        concatenated = concat(applied, join_axes=[concat_index],
                              axis=self.axis, verify_integrity=False)
        return concatenated.reindex_like(obj)

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

def generate_groups(data, group_index, ngroups, axis=0, factory=lambda x: x):
    """
    Parameters
    ----------
    data : BlockManager

    Returns
    -------
    generator
    """
    group_index = com._ensure_int32(group_index)

    indexer = lib.groupsort_indexer(group_index, ngroups)[0]
    group_index = group_index.take(indexer)

    if isinstance(data, BlockManager):
        # this is sort of wasteful but...
        sorted_axis = data.axes[axis].take(indexer)
        sorted_data = data.reindex_axis(sorted_axis, axis=axis)
    if isinstance(data, Series):
        sorted_axis = data.index.take(indexer)
        sorted_data = data.reindex(sorted_axis)
    elif isinstance(data, DataFrame):
        sorted_data = data.take(indexer, axis=axis)

    if isinstance(sorted_data, DataFrame):
        def _get_slice(slob):
            if axis == 0:
                return sorted_data[slob]
            else:
                return sorted_data.ix[:, slob]
    elif isinstance(sorted_data, BlockManager):
        def _get_slice(slob):
            return factory(sorted_data.get_slice(slob, axis=axis))
    elif isinstance(sorted_data, Series):
        def _get_slice(slob):
            return sorted_data._get_values(slob)
    else:  # pragma: no cover
        def _get_slice(slob):
            return sorted_data[slob]

    starts, ends = lib.generate_slices(group_index, ngroups)

    for i, (start, end) in enumerate(zip(starts, ends)):
        # Since I'm now compressing the group ids, it's now not "possible" to
        # produce empty slices because such groups would not be observed in the
        # data
        assert(start < end)
        yield i, _get_slice(slice(start, end))

def get_group_index(label_list, shape):
    """
    For the particular label_list, gets the offsets into the hypothetical list
    representing the totally ordered cartesian product of all possible label
    combinations.
    """
    if len(label_list) == 1:
        return label_list[0]

    n = len(label_list[0])
    group_index = np.zeros(n, dtype=int)
    mask = np.zeros(n, dtype=bool)
    for i in xrange(len(shape)):
        stride = np.prod([x for x in shape[i+1:]], dtype=int)
        group_index += label_list[i] * stride
        mask |= label_list[i] < 0

    np.putmask(group_index, mask, -1)
    return group_index

def decons_group_index(comp_labels, shape):
    # reconstruct labels
    label_list = []
    factor = 1
    y = 0
    x = comp_labels
    for i in reversed(xrange(len(shape))):
        labels = (x - y) % (factor * shape[i]) // factor
        np.putmask(labels, comp_labels < 0, -1)
        label_list.append(labels)
        y = labels * factor
        factor *= shape[i]
    return label_list[::-1]


class _KeyMapper(object):
    """
    Ease my suffering. Map compressed group id -> key tuple
    """
    def __init__(self, comp_ids, ngroups, labels, levels):
        self.levels = levels
        self.labels = labels
        self.comp_ids = comp_ids.astype('i8')

        self.k = len(labels)
        self.tables = [lib.Int64HashTable(ngroups) for _ in range(self.k)]

        self._populate_tables()

    def _populate_tables(self):
        for labs, table in zip(self.labels, self.tables):
            table.map(self.comp_ids, labs.astype('i8'))

    def get_key(self, comp_id):
        return tuple(level[table.get_item(comp_id)]
                     for table, level in zip(self.tables, self.levels))

#----------------------------------------------------------------------
# Group aggregations in Cython


def cython_aggregate(values, group_index, ngroups, how='add'):
    agg_func = _cython_functions[how]
    if values.ndim == 1:
        squeeze = True
        values = values[:, None]
        out_shape = (ngroups, 1)
    else:
        squeeze = False
        out_shape = (ngroups, values.shape[1])

    trans_func = _cython_transforms.get(how, lambda x: x)

    # will be filled in Cython function
    result = np.empty(out_shape, dtype=np.float64)
    counts = np.zeros(ngroups, dtype=np.int32)

    agg_func(result, counts, values, group_index)
    result = trans_func(result)

    if squeeze:
        result = result.squeeze()

    return result, counts

_cython_functions = {
    'add' : lib.group_add,
    'mean' : lib.group_mean,
    'var' : lib.group_var,
    'std' : lib.group_var
}

_cython_transforms = {
    'std' : np.sqrt
}

#----------------------------------------------------------------------
# sorting levels...cleverly?

def _compress_group_index(group_index, sort=True):
    """
    Group_index is offsets into cartesian product of all possible labels. This
    space can be huge, so this function compresses it, by computing offsets
    (comp_ids) into the list of unique labels (obs_group_ids).
    """

    uniques = []
    table = lib.Int64HashTable(len(group_index))

    group_index = com._ensure_int64(group_index)

    # note, group labels come out ascending (ie, 1,2,3 etc)
    comp_ids = table.get_labels_groupby(group_index, uniques)

    # these are the unique ones we observed, in the order we observed them
    obs_group_ids = np.array(uniques, dtype='i8')

    if sort and len(obs_group_ids) > 0:
        # sorter is index where elements ought to go
        sorter = obs_group_ids.argsort()

        # reverse_indexer is where elements came from
        reverse_indexer = np.empty(len(sorter), dtype='i4')
        reverse_indexer.put(sorter, np.arange(len(sorter)))

        mask = comp_ids < 0

        # move comp_ids to right locations (ie, unsort ascending labels)
        comp_ids = reverse_indexer.take(comp_ids)
        np.putmask(comp_ids, mask, -1)

        # sort observed ids
        obs_group_ids = obs_group_ids.take(sorter)

    return comp_ids, obs_group_ids

def _groupby_indices(values):
    if values.dtype != np.object_:
        values = values.astype('O')
    return lib.groupby_indices(values)

def numpy_groupby(data, labels, axis=0):
    s = np.argsort(labels)
    keys, inv = np.unique(labels, return_inverse = True)
    i = inv.take(s)
    groups_at = np.where(i != np.concatenate(([-1], i[:-1])))[0]
    ordered_data = data.take(s, axis=axis)
    group_sums = np.add.reduceat(ordered_data, groups_at, axis=axis)

    return group_sums
