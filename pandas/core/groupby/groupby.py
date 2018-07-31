"""
Provide the groupby split-apply-combine paradigm. Define the GroupBy
class providing the base-class of operations.

The SeriesGroupBy and DataFrameGroupBy sub-class
(defined in pandas.core.groupby.generic)
expose these user-facing objects to provide specific functionailty.
"""

import types
from functools import wraps, partial
import datetime
import collections
import warnings
from contextlib import contextmanager

import numpy as np

from pandas._libs import groupby as libgroupby, Timestamp
from pandas.util._validators import validate_kwargs
from pandas.util._decorators import (
    cache_readonly, Substitution, Appender)

from pandas import compat
from pandas.compat import zip, range, callable, set_function_name
from pandas.compat.numpy import function as nv

from pandas.core.dtypes.common import (
    is_numeric_dtype,
    is_scalar,
    ensure_float)
from pandas.core.dtypes.cast import maybe_downcast_to_dtype
from pandas.core.dtypes.missing import isna, notna

from pandas.core.groupby import base
from pandas.core.base import (PandasObject, SelectionMixin, GroupByError,
                              DataError, SpecificationError)
from pandas.core.index import Index, MultiIndex
from pandas.core.generic import NDFrame
from pandas.core.frame import DataFrame
from pandas.core.series import Series
from pandas.core.sorting import get_group_index_sorter
import pandas.core.common as com
import pandas.core.algorithms as algorithms
from pandas.core.config import option_context

_doc_template = """

        See also
        --------
        pandas.Series.%(name)s
        pandas.DataFrame.%(name)s
        pandas.Panel.%(name)s
"""

_apply_docs = dict(
    template="""
    Apply function `func`  group-wise and combine the results together.

    The function passed to `apply` must take a {input} as its first
    argument and return a DataFrame, Series or scalar. `apply` will
    then take care of combining the results back together into a single
    dataframe or series. `apply` is therefore a highly flexible
    grouping method.

    While `apply` is a very flexible method, its downside is that
    using it can be quite a bit slower than using more specific methods
    like `agg` or `transform`. Pandas offers a wide range of method that will
    be much faster than using `apply` for their specific purposes, so try to
    use them before reaching for `apply`.

    Parameters
    ----------
    func : callable
        A callable that takes a {input} as its first argument, and
        returns a dataframe, a series or a scalar. In addition the
        callable may take positional and keyword arguments.
    args, kwargs : tuple and dict
        Optional positional and keyword arguments to pass to `func`.

    Returns
    -------
    applied : Series or DataFrame

    Notes
    -----
    In the current implementation `apply` calls `func` twice on the
    first group to decide whether it can take a fast or slow code
    path. This can lead to unexpected behavior if `func` has
    side-effects, as they will take effect twice for the first
    group.

    Examples
    --------
    {examples}

    See also
    --------
    pipe : Apply function to the full GroupBy object instead of to each
        group.
    aggregate : Apply aggregate function to the GroupBy object.
    transform : Apply function column-by-column to the GroupBy object.
    Series.apply : Apply a function to a Series.
    DataFrame.apply : Apply a function to each row or column of a DataFrame.
    """,
    dataframe_examples="""
    >>> df = pd.DataFrame({'A': 'a a b'.split(),
                           'B': [1,2,3],
                           'C': [4,6, 5]})
    >>> g = df.groupby('A')

    Notice that ``g`` has two groups, ``a`` and ``b``.
    Calling `apply` in various ways, we can get different grouping results:

    Example 1: below the function passed to `apply` takes a DataFrame as
    its argument and returns a DataFrame. `apply` combines the result for
    each group together into a new DataFrame:

    >>> g[['B', 'C']].apply(lambda x: x / x.sum())
              B    C
    0  0.333333  0.4
    1  0.666667  0.6
    2  1.000000  1.0

    Example 2: The function passed to `apply` takes a DataFrame as
    its argument and returns a Series.  `apply` combines the result for
    each group together into a new DataFrame:

    >>> g[['B', 'C']].apply(lambda x: x.max() - x.min())
       B  C
    A
    a  1  2
    b  0  0

    Example 3: The function passed to `apply` takes a DataFrame as
    its argument and returns a scalar. `apply` combines the result for
    each group together into a Series, including setting the index as
    appropriate:

    >>> g.apply(lambda x: x.C.max() - x.B.min())
    A
    a    5
    b    2
    dtype: int64
    """,
    series_examples="""
    >>> s = pd.Series([0, 1, 2], index='a a b'.split())
    >>> g = s.groupby(s.index)

    From ``s`` above we can see that ``g`` has two groups, ``a`` and ``b``.
    Calling `apply` in various ways, we can get different grouping results:

    Example 1: The function passed to `apply` takes a Series as
    its argument and returns a Series.  `apply` combines the result for
    each group together into a new Series:

    >>> g.apply(lambda x:  x*2 if x.name == 'b' else x/2)
    0    0.0
    1    0.5
    2    4.0
    dtype: float64

    Example 2: The function passed to `apply` takes a Series as
    its argument and returns a scalar. `apply` combines the result for
    each group together into a Series, including setting the index as
    appropriate:

    >>> g.apply(lambda x: x.max() - x.min())
    a    1
    b    0
    dtype: int64
    """)

_pipe_template = """\
Apply a function `func` with arguments to this %(klass)s object and return
the function's result.

%(versionadded)s

Use `.pipe` when you want to improve readability by chaining together
functions that expect Series, DataFrames, GroupBy or Resampler objects.
Instead of writing

>>> h(g(f(df.groupby('group')), arg1=a), arg2=b, arg3=c)

You can write

>>> (df.groupby('group')
...    .pipe(f)
...    .pipe(g, arg1=a)
...    .pipe(h, arg2=b, arg3=c))

which is much more readable.

Parameters
----------
func : callable or tuple of (callable, string)
    Function to apply to this %(klass)s object or, alternatively,
    a `(callable, data_keyword)` tuple where `data_keyword` is a
    string indicating the keyword of `callable` that expects the
    %(klass)s object.
args : iterable, optional
       positional arguments passed into `func`.
kwargs : dict, optional
         a dictionary of keyword arguments passed into `func`.

Returns
-------
object : the return type of `func`.

Notes
-----
See more `here
<http://pandas.pydata.org/pandas-docs/stable/groupby.html#piping-function-calls>`_

Examples
--------
%(examples)s

See Also
--------
pandas.Series.pipe : Apply a function with arguments to a series
pandas.DataFrame.pipe: Apply a function with arguments to a dataframe
apply : Apply function to each group instead of to the
    full %(klass)s object.
"""

_transform_template = """
Call function producing a like-indexed %(klass)s on each group and
return a %(klass)s having the same indexes as the original object
filled with the transformed values

Parameters
----------
f : function
    Function to apply to each group

Notes
-----
Each group is endowed the attribute 'name' in case you need to know
which group you are working on.

The current implementation imposes three requirements on f:

* f must return a value that either has the same shape as the input
  subframe or can be broadcast to the shape of the input subframe.
  For example, f returns a scalar it will be broadcast to have the
  same shape as the input subframe.
* if this is a DataFrame, f must support application column-by-column
  in the subframe. If f also supports application to the entire subframe,
  then a fast path is used starting from the second chunk.
* f must not mutate groups. Mutation is not supported and may
  produce unexpected results.

Returns
-------
%(klass)s

See also
--------
aggregate, transform

Examples
--------

# Same shape
>>> df = pd.DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
...                           'foo', 'bar'],
...                    'B' : ['one', 'one', 'two', 'three',
...                          'two', 'two'],
...                    'C' : [1, 5, 5, 2, 5, 5],
...                    'D' : [2.0, 5., 8., 1., 2., 9.]})
>>> grouped = df.groupby('A')
>>> grouped.transform(lambda x: (x - x.mean()) / x.std())
          C         D
0 -1.154701 -0.577350
1  0.577350  0.000000
2  0.577350  1.154701
3 -1.154701 -1.000000
4  0.577350 -0.577350
5  0.577350  1.000000

# Broadcastable
>>> grouped.transform(lambda x: x.max() - x.min())
   C    D
0  4  6.0
1  3  8.0
2  4  6.0
3  3  8.0
4  4  6.0
5  3  8.0

"""


class GroupByPlot(PandasObject):
    """
    Class implementing the .plot attribute for groupby objects
    """

    def __init__(self, groupby):
        self._groupby = groupby

    def __call__(self, *args, **kwargs):
        def f(self):
            return self.plot(*args, **kwargs)
        f.__name__ = 'plot'
        return self._groupby.apply(f)

    def __getattr__(self, name):
        def attr(*args, **kwargs):
            def f(self):
                return getattr(self.plot, name)(*args, **kwargs)
            # GH-21609
            # naming the function in order to later determine
            # if it requires skipping during fast_apply
            # some functions (e.g. scatter) would
            # create spurious extra plots if fast path is taken
            f.__name__ = name
            return self._groupby.apply(f)
        return attr


@contextmanager
def _group_selection_context(groupby):
    """
    set / reset the _group_selection_context
    """
    groupby._set_group_selection()
    yield groupby
    groupby._reset_group_selection()


class _GroupBy(PandasObject, SelectionMixin):
    _group_selection = None
    _apply_whitelist = frozenset([])

    def __init__(self, obj, keys=None, axis=0, level=None,
                 grouper=None, exclusions=None, selection=None, as_index=True,
                 sort=True, group_keys=True, squeeze=False,
                 observed=False, **kwargs):

        self._selection = selection

        if isinstance(obj, NDFrame):
            obj._consolidate_inplace()

        self.level = level

        if not as_index:
            if not isinstance(obj, DataFrame):
                raise TypeError('as_index=False only valid with DataFrame')
            if axis != 0:
                raise ValueError('as_index=False only valid for axis=0')

        self.as_index = as_index
        self.keys = keys
        self.sort = sort
        self.group_keys = group_keys
        self.squeeze = squeeze
        self.observed = observed
        self.mutated = kwargs.pop('mutated', False)

        if grouper is None:
            from pandas.core.groupby.grouper import _get_grouper
            grouper, exclusions, obj = _get_grouper(obj, keys,
                                                    axis=axis,
                                                    level=level,
                                                    sort=sort,
                                                    observed=observed,
                                                    mutated=self.mutated)

        self.obj = obj
        self.axis = obj._get_axis_number(axis)
        self.grouper = grouper
        self.exclusions = set(exclusions) if exclusions else set()

        # we accept no other args
        validate_kwargs('group', kwargs, {})

    def __len__(self):
        return len(self.groups)

    def __unicode__(self):
        # TODO: Better unicode/repr for GroupBy object
        return object.__repr__(self)

    def _assure_grouper(self):
        """
        we create the grouper on instantiation
        sub-classes may have a different policy
        """
        pass

    @property
    def groups(self):
        """ dict {group name -> group labels} """
        self._assure_grouper()
        return self.grouper.groups

    @property
    def ngroups(self):
        self._assure_grouper()
        return self.grouper.ngroups

    @property
    def indices(self):
        """ dict {group name -> group indices} """
        self._assure_grouper()
        return self.grouper.indices

    def _get_indices(self, names):
        """
        safe get multiple indices, translate keys for
        datelike to underlying repr
        """

        def get_converter(s):
            # possibly convert to the actual key types
            # in the indices, could be a Timestamp or a np.datetime64
            if isinstance(s, (Timestamp, datetime.datetime)):
                return lambda key: Timestamp(key)
            elif isinstance(s, np.datetime64):
                return lambda key: Timestamp(key).asm8
            else:
                return lambda key: key

        if len(names) == 0:
            return []

        if len(self.indices) > 0:
            index_sample = next(iter(self.indices))
        else:
            index_sample = None     # Dummy sample

        name_sample = names[0]
        if isinstance(index_sample, tuple):
            if not isinstance(name_sample, tuple):
                msg = ("must supply a tuple to get_group with multiple"
                       " grouping keys")
                raise ValueError(msg)
            if not len(name_sample) == len(index_sample):
                try:
                    # If the original grouper was a tuple
                    return [self.indices[name] for name in names]
                except KeyError:
                    # turns out it wasn't a tuple
                    msg = ("must supply a a same-length tuple to get_group"
                           " with multiple grouping keys")
                    raise ValueError(msg)

            converters = [get_converter(s) for s in index_sample]
            names = [tuple(f(n) for f, n in zip(converters, name))
                     for name in names]

        else:
            converter = get_converter(index_sample)
            names = [converter(name) for name in names]

        return [self.indices.get(name, []) for name in names]

    def _get_index(self, name):
        """ safe get index, translate keys for datelike to underlying repr """
        return self._get_indices([name])[0]

    @cache_readonly
    def _selected_obj(self):

        if self._selection is None or isinstance(self.obj, Series):
            if self._group_selection is not None:
                return self.obj[self._group_selection]
            return self.obj
        else:
            return self.obj[self._selection]

    def _reset_group_selection(self):
        """
        Clear group based selection. Used for methods needing to return info on
        each group regardless of whether a group selection was previously set.
        """
        if self._group_selection is not None:
            # GH12839 clear cached selection too when changing group selection
            self._group_selection = None
            self._reset_cache('_selected_obj')

    def _set_group_selection(self):
        """
        Create group based selection. Used when selection is not passed
        directly but instead via a grouper.

        NOTE: this should be paired with a call to _reset_group_selection
        """
        grp = self.grouper
        if not (self.as_index and
                getattr(grp, 'groupings', None) is not None and
                self.obj.ndim > 1 and
                self._group_selection is None):
            return

        ax = self.obj._info_axis
        groupers = [g.name for g in grp.groupings
                    if g.level is None and g.in_axis]

        if len(groupers):
            # GH12839 clear selected obj cache when group selection changes
            self._group_selection = ax.difference(Index(groupers)).tolist()
            self._reset_cache('_selected_obj')

    def _set_result_index_ordered(self, result):
        # set the result index on the passed values object and
        # return the new object, xref 8046

        # the values/counts are repeated according to the group index
        # shortcut if we have an already ordered grouper
        if not self.grouper.is_monotonic:
            index = Index(np.concatenate(
                self._get_indices(self.grouper.result_index)))
            result.set_axis(index, axis=self.axis, inplace=True)
            result = result.sort_index(axis=self.axis)

        result.set_axis(self.obj._get_axis(self.axis), axis=self.axis,
                        inplace=True)
        return result

    def _dir_additions(self):
        return self.obj._dir_additions() | self._apply_whitelist

    def __getattr__(self, attr):
        if attr in self._internal_names_set:
            return object.__getattribute__(self, attr)
        if attr in self.obj:
            return self[attr]
        if hasattr(self.obj, attr):
            return self._make_wrapper(attr)

        raise AttributeError("%r object has no attribute %r" %
                             (type(self).__name__, attr))

    @Substitution(klass='GroupBy',
                  versionadded='.. versionadded:: 0.21.0',
                  examples="""\
>>> df = pd.DataFrame({'A': 'a b a b'.split(), 'B': [1, 2, 3, 4]})
>>> df
   A  B
0  a  1
1  b  2
2  a  3
3  b  4

To get the difference between each groups maximum and minimum value in one
pass, you can do

>>> df.groupby('A').pipe(lambda x: x.max() - x.min())
   B
A
a  2
b  2""")
    @Appender(_pipe_template)
    def pipe(self, func, *args, **kwargs):
        return com._pipe(self, func, *args, **kwargs)

    plot = property(GroupByPlot)

    def _make_wrapper(self, name):
        if name not in self._apply_whitelist:
            is_callable = callable(getattr(self._selected_obj, name, None))
            kind = ' callable ' if is_callable else ' '
            msg = ("Cannot access{0}attribute {1!r} of {2!r} objects, try "
                   "using the 'apply' method".format(kind, name,
                                                     type(self).__name__))
            raise AttributeError(msg)

        self._set_group_selection()

        # need to setup the selection
        # as are not passed directly but in the grouper
        f = getattr(self._selected_obj, name)
        if not isinstance(f, types.MethodType):
            return self.apply(lambda self: getattr(self, name))

        f = getattr(type(self._selected_obj), name)

        def wrapper(*args, **kwargs):
            # a little trickery for aggregation functions that need an axis
            # argument
            kwargs_with_axis = kwargs.copy()
            if 'axis' not in kwargs_with_axis or \
               kwargs_with_axis['axis'] is None:
                kwargs_with_axis['axis'] = self.axis

            def curried_with_axis(x):
                return f(x, *args, **kwargs_with_axis)

            def curried(x):
                return f(x, *args, **kwargs)

            # preserve the name so we can detect it when calling plot methods,
            # to avoid duplicates
            curried.__name__ = curried_with_axis.__name__ = name

            # special case otherwise extra plots are created when catching the
            # exception below
            if name in base.plotting_methods:
                return self.apply(curried)
            try:
                return self.apply(curried_with_axis)
            except Exception:
                try:
                    return self.apply(curried)
                except Exception:

                    # related to : GH3688
                    # try item-by-item
                    # this can be called recursively, so need to raise
                    # ValueError
                    # if we don't have this method to indicated to aggregate to
                    # mark this column as an error
                    try:
                        return self._aggregate_item_by_item(name,
                                                            *args, **kwargs)
                    except (AttributeError):
                        raise ValueError

        return wrapper

    def get_group(self, name, obj=None):
        """
        Constructs NDFrame from group with provided name

        Parameters
        ----------
        name : object
            the name of the group to get as a DataFrame
        obj : NDFrame, default None
            the NDFrame to take the DataFrame out of.  If
            it is None, the object groupby was called on will
            be used

        Returns
        -------
        group : same type as obj
        """
        if obj is None:
            obj = self._selected_obj

        inds = self._get_index(name)
        if not len(inds):
            raise KeyError(name)

        return obj._take(inds, axis=self.axis)

    def __iter__(self):
        """
        Groupby iterator

        Returns
        -------
        Generator yielding sequence of (name, subsetted object)
        for each group
        """
        return self.grouper.get_iterator(self.obj, axis=self.axis)

    @Appender(_apply_docs['template']
              .format(input="dataframe",
                      examples=_apply_docs['dataframe_examples']))
    def apply(self, func, *args, **kwargs):

        func = self._is_builtin_func(func)

        # this is needed so we don't try and wrap strings. If we could
        # resolve functions to their callable functions prior, this
        # wouldn't be needed
        if args or kwargs:
            if callable(func):

                @wraps(func)
                def f(g):
                    with np.errstate(all='ignore'):
                        return func(g, *args, **kwargs)
            else:
                raise ValueError('func must be a callable if args or '
                                 'kwargs are supplied')
        else:
            f = func

        # ignore SettingWithCopy here in case the user mutates
        with option_context('mode.chained_assignment', None):
            try:
                result = self._python_apply_general(f)
            except Exception:

                # gh-20949
                # try again, with .apply acting as a filtering
                # operation, by excluding the grouping column
                # This would normally not be triggered
                # except if the udf is trying an operation that
                # fails on *some* columns, e.g. a numeric operation
                # on a string grouper column

                with _group_selection_context(self):
                    return self._python_apply_general(f)

        return result

    def _python_apply_general(self, f):
        keys, values, mutated = self.grouper.apply(f, self._selected_obj,
                                                   self.axis)

        return self._wrap_applied_output(
            keys,
            values,
            not_indexed_same=mutated or self.mutated)

    def _iterate_slices(self):
        yield self._selection_name, self._selected_obj

    def transform(self, func, *args, **kwargs):
        raise com.AbstractMethodError(self)

    def _cumcount_array(self, ascending=True):
        """
        Parameters
        ----------
        ascending : bool, default True
            If False, number in reverse, from length of group - 1 to 0.

        Notes
        -----
        this is currently implementing sort=False
        (though the default is sort=True) for groupby in general
        """
        ids, _, ngroups = self.grouper.group_info
        sorter = get_group_index_sorter(ids, ngroups)
        ids, count = ids[sorter], len(ids)

        if count == 0:
            return np.empty(0, dtype=np.int64)

        run = np.r_[True, ids[:-1] != ids[1:]]
        rep = np.diff(np.r_[np.nonzero(run)[0], count])
        out = (~run).cumsum()

        if ascending:
            out -= np.repeat(out[run], rep)
        else:
            out = np.repeat(out[np.r_[run[1:], True]], rep) - out

        rev = np.empty(count, dtype=np.intp)
        rev[sorter] = np.arange(count, dtype=np.intp)
        return out[rev].astype(np.int64, copy=False)

    def _try_cast(self, result, obj, numeric_only=False):
        """
        try to cast the result to our obj original type,
        we may have roundtripped thru object in the mean-time

        if numeric_only is True, then only try to cast numerics
        and not datetimelikes

        """
        if obj.ndim > 1:
            dtype = obj.values.dtype
        else:
            dtype = obj.dtype

        if not is_scalar(result):
            if numeric_only and is_numeric_dtype(dtype) or not numeric_only:
                result = maybe_downcast_to_dtype(result, dtype)

        return result

    def _transform_should_cast(self, func_nm):
        """
        Parameters:
        -----------
        func_nm: str
            The name of the aggregation function being performed

        Returns:
        --------
        bool
            Whether transform should attempt to cast the result of aggregation
        """
        return (self.size().fillna(0) > 0).any() and (
            func_nm not in base.cython_cast_blacklist)

    def _cython_transform(self, how, numeric_only=True, **kwargs):
        output = collections.OrderedDict()
        for name, obj in self._iterate_slices():
            is_numeric = is_numeric_dtype(obj.dtype)
            if numeric_only and not is_numeric:
                continue

            try:
                result, names = self.grouper.transform(obj.values, how,
                                                       **kwargs)
            except NotImplementedError:
                continue
            except AssertionError as e:
                raise GroupByError(str(e))
            if self._transform_should_cast(how):
                output[name] = self._try_cast(result, obj)
            else:
                output[name] = result

        if len(output) == 0:
            raise DataError('No numeric types to aggregate')

        return self._wrap_transformed_output(output, names)

    def _cython_agg_general(self, how, alt=None, numeric_only=True,
                            min_count=-1):
        output = {}
        for name, obj in self._iterate_slices():
            is_numeric = is_numeric_dtype(obj.dtype)
            if numeric_only and not is_numeric:
                continue

            try:
                result, names = self.grouper.aggregate(obj.values, how,
                                                       min_count=min_count)
            except AssertionError as e:
                raise GroupByError(str(e))
            output[name] = self._try_cast(result, obj)

        if len(output) == 0:
            raise DataError('No numeric types to aggregate')

        return self._wrap_aggregated_output(output, names)

    def _python_agg_general(self, func, *args, **kwargs):
        func = self._is_builtin_func(func)
        f = lambda x: func(x, *args, **kwargs)

        # iterate through "columns" ex exclusions to populate output dict
        output = {}
        for name, obj in self._iterate_slices():
            try:
                result, counts = self.grouper.agg_series(obj, f)
                output[name] = self._try_cast(result, obj, numeric_only=True)
            except TypeError:
                continue

        if len(output) == 0:
            return self._python_apply_general(f)

        if self.grouper._filter_empty_groups:

            mask = counts.ravel() > 0
            for name, result in compat.iteritems(output):

                # since we are masking, make sure that we have a float object
                values = result
                if is_numeric_dtype(values.dtype):
                    values = ensure_float(values)

                output[name] = self._try_cast(values[mask], result)

        return self._wrap_aggregated_output(output)

    def _wrap_applied_output(self, *args, **kwargs):
        raise com.AbstractMethodError(self)

    def _concat_objects(self, keys, values, not_indexed_same=False):
        from pandas.core.reshape.concat import concat

        def reset_identity(values):
            # reset the identities of the components
            # of the values to prevent aliasing
            for v in com._not_none(*values):
                ax = v._get_axis(self.axis)
                ax._reset_identity()
            return values

        if not not_indexed_same:
            result = concat(values, axis=self.axis)
            ax = self._selected_obj._get_axis(self.axis)

            if isinstance(result, Series):
                result = result.reindex(ax)
            else:

                # this is a very unfortunate situation
                # we have a multi-index that is NOT lexsorted
                # and we have a result which is duplicated
                # we can't reindex, so we resort to this
                # GH 14776
                if isinstance(ax, MultiIndex) and not ax.is_unique:
                    indexer = algorithms.unique1d(
                        result.index.get_indexer_for(ax.values))
                    result = result.take(indexer, axis=self.axis)
                else:
                    result = result.reindex(ax, axis=self.axis)

        elif self.group_keys:

            values = reset_identity(values)
            if self.as_index:

                # possible MI return case
                group_keys = keys
                group_levels = self.grouper.levels
                group_names = self.grouper.names

                result = concat(values, axis=self.axis, keys=group_keys,
                                levels=group_levels, names=group_names,
                                sort=False)
            else:

                # GH5610, returns a MI, with the first level being a
                # range index
                keys = list(range(len(values)))
                result = concat(values, axis=self.axis, keys=keys)
        else:
            values = reset_identity(values)
            result = concat(values, axis=self.axis)

        if (isinstance(result, Series) and
                getattr(self, '_selection_name', None) is not None):

            result.name = self._selection_name

        return result

    def _apply_filter(self, indices, dropna):
        if len(indices) == 0:
            indices = np.array([], dtype='int64')
        else:
            indices = np.sort(np.concatenate(indices))
        if dropna:
            filtered = self._selected_obj.take(indices, axis=self.axis)
        else:
            mask = np.empty(len(self._selected_obj.index), dtype=bool)
            mask.fill(False)
            mask[indices.astype(int)] = True
            # mask fails to broadcast when passed to where; broadcast manually.
            mask = np.tile(mask, list(self._selected_obj.shape[1:]) + [1]).T
            filtered = self._selected_obj.where(mask)  # Fill with NaNs.
        return filtered


class GroupBy(_GroupBy):

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
    After grouping, see aggregate, apply, and transform functions. Here are
    some other brief notes about usage. When grouping by multiple groups, the
    result index will be a MultiIndex (hierarchical) by default.

    Iteration produces (key, group) tuples, i.e. chunking the data by group. So
    you can write code like:

    ::

        grouped = obj.groupby(keys, axis=axis)
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
    def _bool_agg(self, val_test, skipna):
        """Shared func to call any / all Cython GroupBy implementations"""

        def objs_to_bool(vals):
            try:
                vals = vals.astype(np.bool)
            except ValueError:  # for objects
                vals = np.array([bool(x) for x in vals])

            return vals.view(np.uint8)

        def result_to_bool(result):
            return result.astype(np.bool, copy=False)

        return self._get_cythonized_result('group_any_all', self.grouper,
                                           aggregate=True,
                                           cython_dtype=np.uint8,
                                           needs_values=True,
                                           needs_mask=True,
                                           pre_processing=objs_to_bool,
                                           post_processing=result_to_bool,
                                           val_test=val_test, skipna=skipna)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def any(self, skipna=True):
        """
        Returns True if any value in the group is truthful, else False

        Parameters
        ----------
        skipna : bool, default True
            Flag to ignore nan values during truth testing
        """
        return self._bool_agg('any', skipna)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def all(self, skipna=True):
        """Returns True if all values in the group are truthful, else False

        Parameters
        ----------
        skipna : bool, default True
            Flag to ignore nan values during truth testing
        """
        return self._bool_agg('all', skipna)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def count(self):
        """Compute count of group, excluding missing values"""

        # defined here for API doc
        raise NotImplementedError

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def mean(self, *args, **kwargs):
        """
        Compute mean of groups, excluding missing values.

        Returns
        -------
        pandas.Series or pandas.DataFrame

        Examples
        --------
        >>> df = pd.DataFrame({'A': [1, 1, 2, 1, 2],
        ...                    'B': [np.nan, 2, 3, 4, 5],
        ...                    'C': [1, 2, 1, 1, 2]}, columns=['A', 'B', 'C'])

        Groupby one column and return the mean of the remaining columns in
        each group.

        >>> df.groupby('A').mean()
        >>>
             B         C
        A
        1  3.0  1.333333
        2  4.0  1.500000

        Groupby two columns and return the mean of the remaining column.

        >>> df.groupby(['A', 'B']).mean()
        >>>
               C
        A B
        1 2.0  2
          4.0  1
        2 3.0  1
          5.0  2

        Groupby one column and return the mean of only particular column in
        the group.

        >>> df.groupby('A')['B'].mean()
        >>>
        A
        1    3.0
        2    4.0
        Name: B, dtype: float64
        """
        nv.validate_groupby_func('mean', args, kwargs, ['numeric_only'])
        try:
            return self._cython_agg_general('mean', **kwargs)
        except GroupByError:
            raise
        except Exception:  # pragma: no cover
            with _group_selection_context(self):
                f = lambda x: x.mean(axis=self.axis, **kwargs)
                return self._python_agg_general(f)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def median(self, **kwargs):
        """
        Compute median of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex
        """
        try:
            return self._cython_agg_general('median', **kwargs)
        except GroupByError:
            raise
        except Exception:  # pragma: no cover

            def f(x):
                if isinstance(x, np.ndarray):
                    x = Series(x)
                return x.median(axis=self.axis, **kwargs)
            with _group_selection_context(self):
                return self._python_agg_general(f)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def std(self, ddof=1, *args, **kwargs):
        """
        Compute standard deviation of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex

        Parameters
        ----------
        ddof : integer, default 1
            degrees of freedom
        """

        # TODO: implement at Cython level?
        nv.validate_groupby_func('std', args, kwargs)
        return np.sqrt(self.var(ddof=ddof, **kwargs))

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def var(self, ddof=1, *args, **kwargs):
        """
        Compute variance of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex

        Parameters
        ----------
        ddof : integer, default 1
            degrees of freedom
        """
        nv.validate_groupby_func('var', args, kwargs)
        if ddof == 1:
            return self._cython_agg_general('var', **kwargs)
        else:
            f = lambda x: x.var(ddof=ddof, **kwargs)
            with _group_selection_context(self):
                return self._python_agg_general(f)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def sem(self, ddof=1):
        """
        Compute standard error of the mean of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex

        Parameters
        ----------
        ddof : integer, default 1
            degrees of freedom
        """

        return self.std(ddof=ddof) / np.sqrt(self.count())

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def size(self):
        """Compute group sizes"""
        result = self.grouper.size()

        if isinstance(self.obj, Series):
            result.name = getattr(self.obj, 'name', None)
        return result

    @classmethod
    def _add_numeric_operations(cls):
        """ add numeric operations to the GroupBy generically """

        def groupby_function(name, alias, npfunc,
                             numeric_only=True, _convert=False,
                             min_count=-1):

            _local_template = "Compute %(f)s of group values"

            @Substitution(name='groupby', f=name)
            @Appender(_doc_template)
            @Appender(_local_template)
            def f(self, **kwargs):
                if 'numeric_only' not in kwargs:
                    kwargs['numeric_only'] = numeric_only
                if 'min_count' not in kwargs:
                    kwargs['min_count'] = min_count

                self._set_group_selection()
                try:
                    return self._cython_agg_general(
                        alias, alt=npfunc, **kwargs)
                except AssertionError as e:
                    raise SpecificationError(str(e))
                except Exception:
                    result = self.aggregate(
                        lambda x: npfunc(x, axis=self.axis))
                    if _convert:
                        result = result._convert(datetime=True)
                    return result

            set_function_name(f, name, cls)

            return f

        def first_compat(x, axis=0):

            def first(x):

                x = np.asarray(x)
                x = x[notna(x)]
                if len(x) == 0:
                    return np.nan
                return x[0]

            if isinstance(x, DataFrame):
                return x.apply(first, axis=axis)
            else:
                return first(x)

        def last_compat(x, axis=0):

            def last(x):

                x = np.asarray(x)
                x = x[notna(x)]
                if len(x) == 0:
                    return np.nan
                return x[-1]

            if isinstance(x, DataFrame):
                return x.apply(last, axis=axis)
            else:
                return last(x)

        cls.sum = groupby_function('sum', 'add', np.sum, min_count=0)
        cls.prod = groupby_function('prod', 'prod', np.prod, min_count=0)
        cls.min = groupby_function('min', 'min', np.min, numeric_only=False)
        cls.max = groupby_function('max', 'max', np.max, numeric_only=False)
        cls.first = groupby_function('first', 'first', first_compat,
                                     numeric_only=False)
        cls.last = groupby_function('last', 'last', last_compat,
                                    numeric_only=False)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def ohlc(self):
        """
        Compute sum of values, excluding missing values
        For multiple groupings, the result index will be a MultiIndex
        """

        return self._apply_to_column_groupbys(
            lambda x: x._cython_agg_general('ohlc'))

    @Appender(DataFrame.describe.__doc__)
    def describe(self, **kwargs):
        with _group_selection_context(self):
            result = self.apply(lambda x: x.describe(**kwargs))
            if self.axis == 1:
                return result.T
            return result.unstack()

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def resample(self, rule, *args, **kwargs):
        """
        Provide resampling when using a TimeGrouper
        Return a new grouper with our resampler appended
        """
        from pandas.core.resample import get_resampler_for_grouping
        return get_resampler_for_grouping(self, rule, *args, **kwargs)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def rolling(self, *args, **kwargs):
        """
        Return a rolling grouper, providing rolling
        functionality per group

        """
        from pandas.core.window import RollingGroupby
        return RollingGroupby(self, *args, **kwargs)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def expanding(self, *args, **kwargs):
        """
        Return an expanding grouper, providing expanding
        functionality per group

        """
        from pandas.core.window import ExpandingGroupby
        return ExpandingGroupby(self, *args, **kwargs)

    def _fill(self, direction, limit=None):
        """Shared function for `pad` and `backfill` to call Cython method

        Parameters
        ----------
        direction : {'ffill', 'bfill'}
            Direction passed to underlying Cython function. `bfill` will cause
            values to be filled backwards. `ffill` and any other values will
            default to a forward fill
        limit : int, default None
            Maximum number of consecutive values to fill. If `None`, this
            method will convert to -1 prior to passing to Cython

        Returns
        -------
        `Series` or `DataFrame` with filled values

        See Also
        --------
        pad
        backfill
        """
        # Need int value for Cython
        if limit is None:
            limit = -1

        return self._get_cythonized_result('group_fillna_indexer',
                                           self.grouper, needs_mask=True,
                                           cython_dtype=np.int64,
                                           result_is_index=True,
                                           direction=direction, limit=limit)

    @Substitution(name='groupby')
    def pad(self, limit=None):
        """
        Forward fill the values

        Parameters
        ----------
        limit : integer, optional
            limit of how many values to fill

        See Also
        --------
        Series.pad
        DataFrame.pad
        Series.fillna
        DataFrame.fillna
        """
        return self._fill('ffill', limit=limit)
    ffill = pad

    @Substitution(name='groupby')
    def backfill(self, limit=None):
        """
        Backward fill the values

        Parameters
        ----------
        limit : integer, optional
            limit of how many values to fill

        See Also
        --------
        Series.backfill
        DataFrame.backfill
        Series.fillna
        DataFrame.fillna
        """
        return self._fill('bfill', limit=limit)
    bfill = backfill

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def nth(self, n, dropna=None):
        """
        Take the nth row from each group if n is an int, or a subset of rows
        if n is a list of ints.

        If dropna, will take the nth non-null row, dropna is either
        Truthy (if a Series) or 'all', 'any' (if a DataFrame);
        this is equivalent to calling dropna(how=dropna) before the
        groupby.

        Parameters
        ----------
        n : int or list of ints
            a single nth value for the row or a list of nth values
        dropna : None or str, optional
            apply the specified dropna operation before counting which row is
            the nth row. Needs to be None, 'any' or 'all'

        Examples
        --------

        >>> df = pd.DataFrame({'A': [1, 1, 2, 1, 2],
        ...                    'B': [np.nan, 2, 3, 4, 5]}, columns=['A', 'B'])
        >>> g = df.groupby('A')
        >>> g.nth(0)
             B
        A
        1  NaN
        2  3.0
        >>> g.nth(1)
             B
        A
        1  2.0
        2  5.0
        >>> g.nth(-1)
             B
        A
        1  4.0
        2  5.0
        >>> g.nth([0, 1])
             B
        A
        1  NaN
        1  2.0
        2  3.0
        2  5.0

        Specifying `dropna` allows count ignoring ``NaN``

        >>> g.nth(0, dropna='any')
             B
        A
        1  2.0
        2  3.0

        NaNs denote group exhausted when using dropna

        >>> g.nth(3, dropna='any')
            B
        A
        1 NaN
        2 NaN

        Specifying `as_index=False` in `groupby` keeps the original index.

        >>> df.groupby('A', as_index=False).nth(1)
           A    B
        1  1  2.0
        4  2  5.0
        """

        if isinstance(n, int):
            nth_values = [n]
        elif isinstance(n, (set, list, tuple)):
            nth_values = list(set(n))
            if dropna is not None:
                raise ValueError(
                    "dropna option with a list of nth values is not supported")
        else:
            raise TypeError("n needs to be an int or a list/set/tuple of ints")

        nth_values = np.array(nth_values, dtype=np.intp)
        self._set_group_selection()

        if not dropna:
            mask = np.in1d(self._cumcount_array(), nth_values) | \
                np.in1d(self._cumcount_array(ascending=False) + 1, -nth_values)

            out = self._selected_obj[mask]
            if not self.as_index:
                return out

            ids, _, _ = self.grouper.group_info
            out.index = self.grouper.result_index[ids[mask]]

            return out.sort_index() if self.sort else out

        if dropna not in ['any', 'all']:
            if isinstance(self._selected_obj, Series) and dropna is True:
                warnings.warn("the dropna={dropna} keyword is deprecated,"
                              "use dropna='all' instead. "
                              "For a Series groupby, dropna must be "
                              "either None, 'any' or 'all'.".format(
                                  dropna=dropna),
                              FutureWarning,
                              stacklevel=2)
                dropna = 'all'
            else:
                # Note: when agg-ing picker doesn't raise this,
                # just returns NaN
                raise ValueError("For a DataFrame groupby, dropna must be "
                                 "either None, 'any' or 'all', "
                                 "(was passed %s)." % (dropna),)

        # old behaviour, but with all and any support for DataFrames.
        # modified in GH 7559 to have better perf
        max_len = n if n >= 0 else - 1 - n
        dropped = self.obj.dropna(how=dropna, axis=self.axis)

        # get a new grouper for our dropped obj
        if self.keys is None and self.level is None:

            # we don't have the grouper info available
            # (e.g. we have selected out
            # a column that is not in the current object)
            axis = self.grouper.axis
            grouper = axis[axis.isin(dropped.index)]

        else:

            # create a grouper with the original parameters, but on the dropped
            # object
            from pandas.core.groupby.grouper import _get_grouper
            grouper, _, _ = _get_grouper(dropped, key=self.keys,
                                         axis=self.axis, level=self.level,
                                         sort=self.sort,
                                         mutated=self.mutated)

        grb = dropped.groupby(grouper, as_index=self.as_index, sort=self.sort)
        sizes, result = grb.size(), grb.nth(n)
        mask = (sizes < max_len).values

        # set the results which don't meet the criteria
        if len(result) and mask.any():
            result.loc[mask] = np.nan

        # reset/reindex to the original groups
        if len(self.obj) == len(dropped) or \
           len(result) == len(self.grouper.result_index):
            result.index = self.grouper.result_index
        else:
            result = result.reindex(self.grouper.result_index)

        return result

    @Substitution(name='groupby')
    def ngroup(self, ascending=True):
        """
        Number each group from 0 to the number of groups - 1.

        This is the enumerative complement of cumcount.  Note that the
        numbers given to the groups match the order in which the groups
        would be seen when iterating over the groupby object, not the
        order they are first observed.

        .. versionadded:: 0.20.2

        Parameters
        ----------
        ascending : bool, default True
            If False, number in reverse, from number of group - 1 to 0.

        Examples
        --------

        >>> df = pd.DataFrame({"A": list("aaabba")})
        >>> df
           A
        0  a
        1  a
        2  a
        3  b
        4  b
        5  a
        >>> df.groupby('A').ngroup()
        0    0
        1    0
        2    0
        3    1
        4    1
        5    0
        dtype: int64
        >>> df.groupby('A').ngroup(ascending=False)
        0    1
        1    1
        2    1
        3    0
        4    0
        5    1
        dtype: int64
        >>> df.groupby(["A", [1,1,2,3,2,1]]).ngroup()
        0    0
        1    0
        2    1
        3    3
        4    2
        5    0
        dtype: int64

        See also
        --------
        .cumcount : Number the rows in each group.
        """

        with _group_selection_context(self):
            index = self._selected_obj.index
            result = Series(self.grouper.group_info[0], index)
            if not ascending:
                result = self.ngroups - 1 - result
            return result

    @Substitution(name='groupby')
    def cumcount(self, ascending=True):
        """
        Number each item in each group from 0 to the length of that group - 1.

        Essentially this is equivalent to

        >>> self.apply(lambda x: pd.Series(np.arange(len(x)), x.index))

        Parameters
        ----------
        ascending : bool, default True
            If False, number in reverse, from length of group - 1 to 0.

        Examples
        --------

        >>> df = pd.DataFrame([['a'], ['a'], ['a'], ['b'], ['b'], ['a']],
        ...                   columns=['A'])
        >>> df
           A
        0  a
        1  a
        2  a
        3  b
        4  b
        5  a
        >>> df.groupby('A').cumcount()
        0    0
        1    1
        2    2
        3    0
        4    1
        5    3
        dtype: int64
        >>> df.groupby('A').cumcount(ascending=False)
        0    3
        1    2
        2    1
        3    1
        4    0
        5    0
        dtype: int64

        See also
        --------
        .ngroup : Number the groups themselves.
        """

        with _group_selection_context(self):
            index = self._selected_obj.index
            cumcounts = self._cumcount_array(ascending=ascending)
            return Series(cumcounts, index)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def rank(self, method='average', ascending=True, na_option='keep',
             pct=False, axis=0):
        """
        Provides the rank of values within each group.

        Parameters
        ----------
        method : {'average', 'min', 'max', 'first', 'dense'}, default 'average'
            * average: average rank of group
            * min: lowest rank in group
            * max: highest rank in group
            * first: ranks assigned in order they appear in the array
            * dense: like 'min', but rank always increases by 1 between groups
        ascending : boolean, default True
            False for ranks by high (1) to low (N)
        na_option :  {'keep', 'top', 'bottom'}, default 'keep'
            * keep: leave NA values where they are
            * top: smallest rank if ascending
            * bottom: smallest rank if descending
        pct : boolean, default False
            Compute percentage rank of data within each group
        axis : int, default 0
            The axis of the object over which to compute the rank.

        Returns
        -----
        DataFrame with ranking of values within each group
        """
        return self._cython_transform('rank', numeric_only=False,
                                      ties_method=method, ascending=ascending,
                                      na_option=na_option, pct=pct, axis=axis)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def cumprod(self, axis=0, *args, **kwargs):
        """Cumulative product for each group"""
        nv.validate_groupby_func('cumprod', args, kwargs,
                                 ['numeric_only', 'skipna'])
        if axis != 0:
            return self.apply(lambda x: x.cumprod(axis=axis, **kwargs))

        return self._cython_transform('cumprod', **kwargs)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def cumsum(self, axis=0, *args, **kwargs):
        """Cumulative sum for each group"""
        nv.validate_groupby_func('cumsum', args, kwargs,
                                 ['numeric_only', 'skipna'])
        if axis != 0:
            return self.apply(lambda x: x.cumsum(axis=axis, **kwargs))

        return self._cython_transform('cumsum', **kwargs)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def cummin(self, axis=0, **kwargs):
        """Cumulative min for each group"""
        if axis != 0:
            return self.apply(lambda x: np.minimum.accumulate(x, axis))

        return self._cython_transform('cummin', numeric_only=False)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def cummax(self, axis=0, **kwargs):
        """Cumulative max for each group"""
        if axis != 0:
            return self.apply(lambda x: np.maximum.accumulate(x, axis))

        return self._cython_transform('cummax', numeric_only=False)

    def _get_cythonized_result(self, how, grouper, aggregate=False,
                               cython_dtype=None, needs_values=False,
                               needs_mask=False, needs_ngroups=False,
                               result_is_index=False,
                               pre_processing=None, post_processing=None,
                               **kwargs):
        """Get result for Cythonized functions

        Parameters
        ----------
        how : str, Cythonized function name to be called
        grouper : Grouper object containing pertinent group info
        aggregate : bool, default False
            Whether the result should be aggregated to match the number of
            groups
        cython_dtype : default None
            Type of the array that will be modified by the Cython call. If
            `None`, the type will be inferred from the values of each slice
        needs_values : bool, default False
            Whether the values should be a part of the Cython call
            signature
        needs_mask : bool, default False
            Whether boolean mask needs to be part of the Cython call
            signature
        needs_ngroups : bool, default False
            Whether number of groups is part of the Cython call signature
        result_is_index : bool, default False
            Whether the result of the Cython operation is an index of
            values to be retrieved, instead of the actual values themselves
        pre_processing : function, default None
            Function to be applied to `values` prior to passing to Cython
            Raises if `needs_values` is False
        post_processing : function, default None
            Function to be applied to result of Cython function
        **kwargs : dict
            Extra arguments to be passed back to Cython funcs

        Returns
        -------
        `Series` or `DataFrame`  with filled values
        """
        if result_is_index and aggregate:
            raise ValueError("'result_is_index' and 'aggregate' cannot both "
                             "be True!")
        if post_processing:
            if not callable(pre_processing):
                raise ValueError("'post_processing' must be a callable!")
        if pre_processing:
            if not callable(pre_processing):
                raise ValueError("'pre_processing' must be a callable!")
            if not needs_values:
                raise ValueError("Cannot use 'pre_processing' without "
                                 "specifying 'needs_values'!")

        labels, _, ngroups = grouper.group_info
        output = collections.OrderedDict()
        base_func = getattr(libgroupby, how)

        for name, obj in self._iterate_slices():
            if aggregate:
                result_sz = ngroups
            else:
                result_sz = len(obj.values)

            if not cython_dtype:
                cython_dtype = obj.values.dtype

            result = np.zeros(result_sz, dtype=cython_dtype)
            func = partial(base_func, result, labels)
            if needs_values:
                vals = obj.values
                if pre_processing:
                    vals = pre_processing(vals)
                func = partial(func, vals)

            if needs_mask:
                mask = isna(obj.values).view(np.uint8)
                func = partial(func, mask)

            if needs_ngroups:
                func = partial(func, ngroups)

            func(**kwargs)  # Call func to modify indexer values in place

            if result_is_index:
                result = algorithms.take_nd(obj.values, result)

            if post_processing:
                result = post_processing(result)

            output[name] = result

        if aggregate:
            return self._wrap_aggregated_output(output)
        else:
            return self._wrap_transformed_output(output)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def shift(self, periods=1, freq=None, axis=0):
        """
        Shift each group by periods observations

        Parameters
        ----------
        periods : integer, default 1
            number of periods to shift
        freq : frequency string
        axis : axis to shift, default 0
        """

        if freq is not None or axis != 0:
            return self.apply(lambda x: x.shift(periods, freq, axis))

        return self._get_cythonized_result('group_shift_indexer',
                                           self.grouper, cython_dtype=np.int64,
                                           needs_ngroups=True,
                                           result_is_index=True,
                                           periods=periods)

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def pct_change(self, periods=1, fill_method='pad', limit=None, freq=None,
                   axis=0):
        """Calculate pct_change of each value to previous entry in group"""
        if freq is not None or axis != 0:
            return self.apply(lambda x: x.pct_change(periods=periods,
                                                     fill_method=fill_method,
                                                     limit=limit, freq=freq,
                                                     axis=axis))

        filled = getattr(self, fill_method)(limit=limit).drop(
            self.grouper.names, axis=1)
        shifted = filled.shift(periods=periods, freq=freq)

        return (filled / shifted) - 1

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def head(self, n=5):
        """
        Returns first n rows of each group.

        Essentially equivalent to ``.apply(lambda x: x.head(n))``,
        except ignores as_index flag.

        Examples
        --------

        >>> df = pd.DataFrame([[1, 2], [1, 4], [5, 6]],
                              columns=['A', 'B'])
        >>> df.groupby('A', as_index=False).head(1)
           A  B
        0  1  2
        2  5  6
        >>> df.groupby('A').head(1)
           A  B
        0  1  2
        2  5  6
        """
        self._reset_group_selection()
        mask = self._cumcount_array() < n
        return self._selected_obj[mask]

    @Substitution(name='groupby')
    @Appender(_doc_template)
    def tail(self, n=5):
        """
        Returns last n rows of each group

        Essentially equivalent to ``.apply(lambda x: x.tail(n))``,
        except ignores as_index flag.

        Examples
        --------

        >>> df = pd.DataFrame([['a', 1], ['a', 2], ['b', 1], ['b', 2]],
                              columns=['A', 'B'])
        >>> df.groupby('A').tail(1)
           A  B
        1  a  2
        3  b  2
        >>> df.groupby('A').head(1)
           A  B
        0  a  1
        2  b  1
        """
        self._reset_group_selection()
        mask = self._cumcount_array(ascending=False) < n
        return self._selected_obj[mask]


GroupBy._add_numeric_operations()


@Appender(GroupBy.__doc__)
def groupby(obj, by, **kwds):
    if isinstance(obj, Series):
        from pandas.core.groupby.generic import SeriesGroupBy
        klass = SeriesGroupBy
    elif isinstance(obj, DataFrame):
        from pandas.core.groupby.generic import DataFrameGroupBy
        klass = DataFrameGroupBy
    else:  # pragma: no cover
        raise TypeError('invalid type: %s' % type(obj))

    return klass(obj, by, **kwds)
