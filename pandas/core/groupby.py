import types
from functools import wraps
import numpy as np
import datetime
import collections
import warnings
import copy

from pandas.compat import(
    zip, builtins, range, long, lzip,
    OrderedDict, callable, filter, map
)
from pandas import compat

from pandas.core.base import PandasObject
from pandas.core.categorical import Categorical
from pandas.core.frame import DataFrame
from pandas.core.generic import NDFrame
from pandas.core.index import Index, MultiIndex, CategoricalIndex, _ensure_index
from pandas.core.internals import BlockManager, make_block
from pandas.core.series import Series
from pandas.core.panel import Panel
from pandas.util.decorators import (cache_readonly, Appender, make_signature,
                                    deprecate_kwarg)
import pandas.core.algorithms as algos
import pandas.core.common as com
from pandas.core.common import(_possibly_downcast_to_dtype, isnull,
                               notnull, _DATELIKE_DTYPES, is_numeric_dtype,
                               is_timedelta64_dtype, is_datetime64_dtype,
                               is_categorical_dtype, _values_from_object,
                               is_datetime_or_timedelta_dtype, is_bool,
                               is_bool_dtype, AbstractMethodError)
from pandas.core.config import option_context
import pandas.lib as lib
from pandas.lib import Timestamp
import pandas.tslib as tslib
import pandas.algos as _algos
import pandas.hashtable as _hash

_agg_doc = """Aggregate using input function or dict of {column -> function}

Parameters
----------
arg : function or dict
    Function to use for aggregating groups. If a function, must either
    work when passed a DataFrame or when passed to DataFrame.apply. If
    passed a dict, the keys must be DataFrame column names.

Notes
-----
Numpy functions mean/median/prod/sum/std/var are special cased so the
default behavior is applying the function along axis=0
(e.g., np.mean(arr_2d, axis=0)) as opposed to
mimicking the default Numpy behavior (e.g., np.mean(arr_2d)).

Returns
-------
aggregated : DataFrame
"""


# special case to prevent duplicate plots when catching exceptions when
# forwarding methods from NDFrames
_plotting_methods = frozenset(['plot', 'boxplot', 'hist'])

_common_apply_whitelist = frozenset([
    'last', 'first',
    'head', 'tail', 'median',
    'mean', 'sum', 'min', 'max',
    'cumsum', 'cumprod', 'cummin', 'cummax', 'cumcount',
    'resample',
    'describe',
    'rank', 'quantile',
    'fillna',
    'mad',
    'any', 'all',
    'take',
    'idxmax', 'idxmin',
    'shift', 'tshift',
    'ffill', 'bfill',
    'pct_change', 'skew',
    'corr', 'cov', 'diff',
]) | _plotting_methods

_series_apply_whitelist = \
    (_common_apply_whitelist - set(['boxplot'])) | \
    frozenset(['dtype', 'unique'])

_dataframe_apply_whitelist = \
    _common_apply_whitelist | frozenset(['dtypes', 'corrwith'])

_cython_transforms = frozenset(['cumprod', 'cumsum', 'shift'])

class GroupByError(Exception):
    pass


class DataError(GroupByError):
    pass


class SpecificationError(GroupByError):
    pass


def _groupby_function(name, alias, npfunc, numeric_only=True,
                      _convert=False):
    def f(self):
        self._set_selection_from_grouper()
        try:
            return self._cython_agg_general(alias, numeric_only=numeric_only)
        except AssertionError as e:
            raise SpecificationError(str(e))
        except Exception:
            result = self.aggregate(lambda x: npfunc(x, axis=self.axis))
            if _convert:
                result = result._convert(datetime=True)
            return result

    f.__doc__ = "Compute %s of group values" % name
    f.__name__ = name

    return f


def _first_compat(x, axis=0):
    def _first(x):
        x = np.asarray(x)
        x = x[notnull(x)]
        if len(x) == 0:
            return np.nan
        return x[0]

    if isinstance(x, DataFrame):
        return x.apply(_first, axis=axis)
    else:
        return _first(x)


def _last_compat(x, axis=0):
    def _last(x):
        x = np.asarray(x)
        x = x[notnull(x)]
        if len(x) == 0:
            return np.nan
        return x[-1]

    if isinstance(x, DataFrame):
        return x.apply(_last, axis=axis)
    else:
        return _last(x)


class Grouper(object):
    """
    A Grouper allows the user to specify a groupby instruction for a target object

    This specification will select a column via the key parameter, or if the level and/or
    axis parameters are given, a level of the index of the target object.

    These are local specifications and will override 'global' settings, that is the parameters
    axis and level which are passed to the groupby itself.

    Parameters
    ----------
    key : string, defaults to None
        groupby key, which selects the grouping column of the target
    level : name/number, defaults to None
        the level for the target index
    freq : string / frequency object, defaults to None
        This will groupby the specified frequency if the target selection (via key or level) is
        a datetime-like object. For full specification of available frequencies, please see
        `here <http://pandas.pydata.org/pandas-docs/stable/timeseries.html>`_.
    axis : number/name of the axis, defaults to 0
    sort : boolean, default to False
        whether to sort the resulting labels

    additional kwargs to control time-like groupers (when freq is passed)

    closed : closed end of interval; left or right
    label : interval boundary to use for labeling; left or right
    convention : {'start', 'end', 'e', 's'}
        If grouper is PeriodIndex

    Returns
    -------
    A specification for a groupby instruction

    Examples
    --------

    Syntactic sugar for ``df.groupby('A')``

    >>> df.groupby(Grouper(key='A'))

    Specify a resample operation on the column 'date'

    >>> df.groupby(Grouper(key='date', freq='60s'))

    Specify a resample operation on the level 'date' on the columns axis
    with a frequency of 60s

    >>> df.groupby(Grouper(level='date', freq='60s', axis=1))
    """

    def __new__(cls, *args, **kwargs):
        if kwargs.get('freq') is not None:
            from pandas.tseries.resample import TimeGrouper
            cls = TimeGrouper
        return super(Grouper, cls).__new__(cls)

    def __init__(self, key=None, level=None, freq=None, axis=0, sort=False):
        self.key=key
        self.level=level
        self.freq=freq
        self.axis=axis
        self.sort=sort

        self.grouper=None
        self.obj=None
        self.indexer=None
        self.binner=None

    @property
    def ax(self):
        return self.grouper

    def _get_grouper(self, obj):

        """
        Parameters
        ----------
        obj : the subject object

        Returns
        -------
        a tuple of binner, grouper, obj (possibly sorted)
        """

        self._set_grouper(obj)
        self.grouper, exclusions, self.obj = _get_grouper(self.obj, [self.key], axis=self.axis,
                                                          level=self.level, sort=self.sort)
        return self.binner, self.grouper, self.obj

    def _set_grouper(self, obj, sort=False):
        """
        given an object and the specifications, setup the internal grouper for this particular specification

        Parameters
        ----------
        obj : the subject object

        """

        if self.key is not None and self.level is not None:
            raise ValueError("The Grouper cannot specify both a key and a level!")

        # the key must be a valid info item
        if self.key is not None:
            key = self.key
            if key not in obj._info_axis:
                raise KeyError("The grouper name {0} is not found".format(key))
            ax = Index(obj[key], name=key)

        else:
            ax = obj._get_axis(self.axis)
            if self.level is not None:
                level = self.level

                # if a level is given it must be a mi level or
                # equivalent to the axis name
                if isinstance(ax, MultiIndex):
                    level = ax._get_level_number(level)
                    ax = Index(ax.get_level_values(level), name=ax.names[level])

                else:
                    if level not in (0, ax.name):
                        raise ValueError("The level {0} is not valid".format(level))

        # possibly sort
        if (self.sort or sort) and not ax.is_monotonic:
            indexer = self.indexer = ax.argsort(kind='quicksort')
            ax = ax.take(indexer)
            obj = obj.take(indexer, axis=self.axis, convert=False, is_copy=False)

        self.obj = obj
        self.grouper = ax
        return self.grouper

    def _get_binner_for_grouping(self, obj):
        """ default to the standard binner here """
        group_axis = obj._get_axis(self.axis)
        return Grouping(group_axis, None, obj=obj, name=self.key,
                        level=self.level, sort=self.sort, in_axis=False)

    @property
    def groups(self):
        return self.grouper.groups


class GroupByPlot(PandasObject):
    """
    Class implementing the .plot attribute for groupby objects
    """
    def __init__(self, groupby):
        self._groupby = groupby

    def __call__(self, *args, **kwargs):
        def f(self, *args, **kwargs):
            return self.plot(*args, **kwargs)
        f.__name__ = 'plot'
        return self._groupby.apply(f)

    def __getattr__(self, name):
        def attr(*args, **kwargs):
            def f(self):
                return getattr(self.plot, name)(*args, **kwargs)
            return self._groupby.apply(f)
        return attr


class GroupBy(PandasObject):

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
    _apply_whitelist = _common_apply_whitelist
    _internal_names = ['_cache']
    _internal_names_set = set(_internal_names)
    _group_selection = None

    def __init__(self, obj, keys=None, axis=0, level=None,
                 grouper=None, exclusions=None, selection=None, as_index=True,
                 sort=True, group_keys=True, squeeze=False):
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

        if grouper is None:
            grouper, exclusions, obj = _get_grouper(obj, keys, axis=axis,
                                                    level=level, sort=sort)

        self.obj = obj
        self.axis = obj._get_axis_number(axis)
        self.grouper = grouper
        self.exclusions = set(exclusions) if exclusions else set()

    def __len__(self):
        return len(self.groups)

    def __unicode__(self):
        # TODO: Better unicode/repr for GroupBy object
        return object.__repr__(self)

    @property
    def groups(self):
        """ dict {group name -> group labels} """
        return self.grouper.groups

    @property
    def ngroups(self):
        return self.grouper.ngroups

    @property
    def indices(self):
        """ dict {group name -> group indices} """
        return self.grouper.indices

    def _get_indices(self, names):
        """ safe get multiple indices, translate keys for datelike to underlying repr """

        def get_converter(s):
            # possibly convert to the actual key types
            # in the indices, could be a Timestamp or a np.datetime64
            if isinstance(s, (Timestamp,datetime.datetime)):
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
            names = [tuple([f(n) for f, n in zip(converters, name)]) for name in names]

        else:
            converter = get_converter(index_sample)
            names = [converter(name) for name in names]

        return [self.indices.get(name, []) for name in names]

    def _get_index(self, name):
        """ safe get index, translate keys for datelike to underlying repr """
        return self._get_indices([name])[0]

    @property
    def name(self):
        if self._selection is None:
            return None  # 'result'
        else:
            return self._selection

    @property
    def _selection_list(self):
        if not isinstance(self._selection, (list, tuple, Series, Index, np.ndarray)):
            return [self._selection]
        return self._selection

    @cache_readonly
    def _selected_obj(self):

        if self._selection is None or isinstance(self.obj, Series):
            if self._group_selection is not None:
                return self.obj[self._group_selection]
            return self.obj
        else:
            return self.obj[self._selection]

    def _set_selection_from_grouper(self):
        """ we may need create a selection if we have non-level groupers """
        grp = self.grouper
        if self.as_index and getattr(grp,'groupings',None) is not None and self.obj.ndim > 1:
            ax = self.obj._info_axis
            groupers = [g.name for g in grp.groupings
                           if g.level is None and g.in_axis]

            if len(groupers):
                self._group_selection = ax.difference(Index(groupers)).tolist()

    def _set_result_index_ordered(self, result):
        # set the result index on the passed values object
        # return the new object
        # related 8046

        # the values/counts are repeated according to the group index
        indices = self.indices

        # shortcut of we have an already ordered grouper
        if not self.grouper.is_monotonic:
            index = Index(np.concatenate(self._get_indices(self.grouper.result_index)))
            result.index = index
            result = result.sort_index()

        result.index = self.obj.index
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

    def __getitem__(self, key):
        raise NotImplementedError('Not implemented: %s' % key)

    plot = property(GroupByPlot)

    def _make_wrapper(self, name):
        if name not in self._apply_whitelist:
            is_callable = callable(getattr(self._selected_obj, name, None))
            kind = ' callable ' if is_callable else ' '
            msg = ("Cannot access{0}attribute {1!r} of {2!r} objects, try "
                   "using the 'apply' method".format(kind, name,
                                                     type(self).__name__))
            raise AttributeError(msg)

        # need to setup the selection
        # as are not passed directly but in the grouper
        self._set_selection_from_grouper()

        f = getattr(self._selected_obj, name)
        if not isinstance(f, types.MethodType):
            return self.apply(lambda self: getattr(self, name))

        f = getattr(type(self._selected_obj), name)

        def wrapper(*args, **kwargs):
            # a little trickery for aggregation functions that need an axis
            # argument
            kwargs_with_axis = kwargs.copy()
            if 'axis' not in kwargs_with_axis or kwargs_with_axis['axis']==None:
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
            if name in _plotting_methods:
                return self.apply(curried)

            try:
                return self.apply(curried_with_axis)
            except Exception:
                try:
                    return self.apply(curried)
                except Exception:

                    # related to : GH3688
                    # try item-by-item
                    # this can be called recursively, so need to raise ValueError if
                    # we don't have this method to indicated to aggregate to
                    # mark this column as an error
                    try:
                        return self._aggregate_item_by_item(name, *args, **kwargs)
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
        group : type of obj
        """
        if obj is None:
            obj = self._selected_obj

        inds = self._get_index(name)
        if not len(inds):
            raise KeyError(name)

        return obj.take(inds, axis=self.axis, convert=False)

    def __iter__(self):
        """
        Groupby iterator

        Returns
        -------
        Generator yielding sequence of (name, subsetted object)
        for each group
        """
        return self.grouper.get_iterator(self.obj, axis=self.axis)

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
        See online documentation for full exposition on how to use apply.

        In the current implementation apply calls func twice on the
        first group to decide whether it can take a fast or slow code
        path. This can lead to unexpected behavior if func has
        side-effects, as they will take effect twice for the first
        group.


        See also
        --------
        aggregate, transform

        Returns
        -------
        applied : type depending on grouped object and function
        """
        func = _intercept_function(func)

        @wraps(func)
        def f(g):
            return func(g, *args, **kwargs)

        # ignore SettingWithCopy here in case the user mutates
        with option_context('mode.chained_assignment',None):
            return self._python_apply_general(f)

    def _python_apply_general(self, f):
        keys, values, mutated = self.grouper.apply(f, self._selected_obj,
                                                   self.axis)

        return self._wrap_applied_output(keys, values,
                                         not_indexed_same=mutated)

    def aggregate(self, func, *args, **kwargs):
        raise AbstractMethodError(self)

    @Appender(_agg_doc)
    def agg(self, func, *args, **kwargs):
        return self.aggregate(func, *args, **kwargs)

    def _iterate_slices(self):
        yield self.name, self._selected_obj

    def transform(self, func, *args, **kwargs):
        raise AbstractMethodError(self)

    def irow(self, i):
        """
        DEPRECATED. Use ``.nth(i)`` instead
        """

        # 10177
        warnings.warn("irow(i) is deprecated. Please use .nth(i)",
                      FutureWarning, stacklevel=2)
        return self.nth(i)

    def count(self):
        """ Compute count of group, excluding missing values """

        # defined here for API doc
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
            self._set_selection_from_grouper()
            f = lambda x: x.mean(axis=self.axis)
            return self._python_agg_general(f)

    def median(self):
        """
        Compute median of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex
        """
        try:
            return self._cython_agg_general('median')
        except GroupByError:
            raise
        except Exception:  # pragma: no cover

            self._set_selection_from_grouper()
            def f(x):
                if isinstance(x, np.ndarray):
                    x = Series(x)
                return x.median(axis=self.axis)
            return self._python_agg_general(f)

    def std(self, ddof=1):
        """
        Compute standard deviation of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex
        """
        # todo, implement at cython level?
        return np.sqrt(self.var(ddof=ddof))

    def var(self, ddof=1):
        """
        Compute variance of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex
        """
        if ddof == 1:
            return self._cython_agg_general('var')
        else:
            self._set_selection_from_grouper()
            f = lambda x: x.var(ddof=ddof)
            return self._python_agg_general(f)

    def sem(self, ddof=1):
        """
        Compute standard error of the mean of groups, excluding missing values

        For multiple groupings, the result index will be a MultiIndex
        """
        return self.std(ddof=ddof)/np.sqrt(self.count())

    def size(self):
        """
        Compute group sizes

        """
        return self.grouper.size()

    sum = _groupby_function('sum', 'add', np.sum)
    prod = _groupby_function('prod', 'prod', np.prod)
    min = _groupby_function('min', 'min', np.min, numeric_only=False)
    max = _groupby_function('max', 'max', np.max, numeric_only=False)
    first = _groupby_function('first', 'first', _first_compat,
                              numeric_only=False, _convert=True)
    last = _groupby_function('last', 'last', _last_compat, numeric_only=False,
                             _convert=True)

    def ohlc(self):
        """
        Compute sum of values, excluding missing values
        For multiple groupings, the result index will be a MultiIndex
        """
        return self._apply_to_column_groupbys(
            lambda x: x._cython_agg_general('ohlc'))

    def nth(self, n, dropna=None):
        """
        Take the nth row from each group if n is an int, or a subset of rows
        if n is a list of ints.

        If dropna, will take the nth non-null row, dropna is either
        Truthy (if a Series) or 'all', 'any' (if a DataFrame); this is equivalent
        to calling dropna(how=dropna) before the groupby.

        Parameters
        ----------
        n : int or list of ints
            a single nth value for the row or a list of nth values
        dropna : None or str, optional
            apply the specified dropna operation before counting which row is
            the nth row. Needs to be None, 'any' or 'all'

        Examples
        --------
        >>> df = DataFrame([[1, np.nan], [1, 4], [5, 6]], columns=['A', 'B'])
        >>> g = df.groupby('A')
        >>> g.nth(0)
           A   B
        0  1 NaN
        2  5   6
        >>> g.nth(1)
           A  B
        1  1  4
        >>> g.nth(-1)
           A  B
        1  1  4
        2  5  6
        >>> g.nth(0, dropna='any')
           B
        A
        1  4
        5  6
        >>> g.nth(1, dropna='any')  # NaNs denote group exhausted when using dropna
            B
        A
        1 NaN
        5 NaN

        """
        if isinstance(n, int):
            nth_values = [n]
        elif isinstance(n, (set, list, tuple)):
            nth_values = list(set(n))
            if dropna is not None:
                raise ValueError("dropna option with a list of nth values is not supported")
        else:
            raise TypeError("n needs to be an int or a list/set/tuple of ints")

        m = self.grouper._max_groupsize
        # filter out values that are outside [-m, m)
        pos_nth_values = [i for i in nth_values if i >= 0 and i < m]
        neg_nth_values = [i for i in nth_values if i < 0 and i >= -m]

        self._set_selection_from_grouper()
        if not dropna:  # good choice
            if not pos_nth_values and not neg_nth_values:
                # no valid nth values
                return self._selected_obj.loc[[]]

            rng = np.zeros(m, dtype=bool)
            for i in pos_nth_values:
                rng[i] = True
            is_nth = self._cumcount_array(rng)

            if neg_nth_values:
                rng = np.zeros(m, dtype=bool)
                for i in neg_nth_values:
                    rng[- i - 1] = True
                is_nth |= self._cumcount_array(rng, ascending=False)

            result = self._selected_obj[is_nth]

            # the result index
            if self.as_index:
                ax = self.obj._info_axis
                names = self.grouper.names
                if self.obj.ndim == 1:
                    # this is a pass-thru
                    pass
                elif all([ n in ax for n in names ]):
                    result.index = MultiIndex.from_arrays([self.obj[name][is_nth] for name in names]).set_names(names)
                elif self._group_selection is not None:
                    result.index = self.obj._get_axis(self.axis)[is_nth]

                result = result.sort_index()

            return result

        if (isinstance(self._selected_obj, DataFrame)
           and dropna not in ['any', 'all']):
            # Note: when agg-ing picker doesn't raise this, just returns NaN
            raise ValueError("For a DataFrame groupby, dropna must be "
                             "either None, 'any' or 'all', "
                             "(was passed %s)." % (dropna),)

        # old behaviour, but with all and any support for DataFrames.
        # modified in GH 7559 to have better perf
        max_len = n if n >= 0 else - 1 - n
        dropped = self.obj.dropna(how=dropna, axis=self.axis)

        # get a new grouper for our dropped obj
        if self.keys is None and self.level is None:

            # we don't have the grouper info available (e.g. we have selected out
            # a column that is not in the current object)
            axis = self.grouper.axis
            grouper = axis[axis.isin(dropped.index)]
            keys = self.grouper.names
        else:

            # create a grouper with the original parameters, but on the dropped object
            grouper, _, _ = _get_grouper(dropped, key=self.keys, axis=self.axis,
                                         level=self.level, sort=self.sort)

        sizes = dropped.groupby(grouper).size()
        result = dropped.groupby(grouper).nth(n)
        mask = (sizes<max_len).values

        # set the results which don't meet the criteria
        if len(result) and mask.any():
            result.loc[mask] = np.nan

        # reset/reindex to the original groups
        if len(self.obj) == len(dropped) or len(result) == len(self.grouper.result_index):
            result.index = self.grouper.result_index
        else:
            result = result.reindex(self.grouper.result_index)

        return result

    def cumcount(self, ascending=True):
        """
        Number each item in each group from 0 to the length of that group - 1.

        Essentially this is equivalent to

        >>> self.apply(lambda x: Series(np.arange(len(x)), x.index))

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

        """
        self._set_selection_from_grouper()

        index = self._selected_obj.index
        cumcounts = self._cumcount_array(ascending=ascending)
        return Series(cumcounts, index)

    def cumprod(self, axis=0):
        """
        Cumulative product for each group

        """
        if axis != 0:
            return self.apply(lambda x: x.cumprod(axis=axis))

        return self._cython_transform('cumprod')

    def cumsum(self, axis=0):
        """
        Cumulative sum for each group

        """
        if axis != 0:
            return self.apply(lambda x: x.cumprod(axis=axis))

        return self._cython_transform('cumsum')

    def shift(self, periods=1, freq=None, axis=0):
        """
        Shift each group by periods observations
        """

        if freq is not None or axis != 0:
            return self.apply(lambda x: x.shift(periods, freq, axis))

        labels, _,  ngroups = self.grouper.group_info
        # filled in by Cython
        indexer = np.zeros_like(labels)
        _algos.group_shift_indexer(indexer, labels, ngroups, periods)

        output = {}
        for name, obj in self._iterate_slices():
            output[name] = com.take_nd(obj.values, indexer)

        return self._wrap_transformed_output(output)

    def head(self, n=5):
        """
        Returns first n rows of each group.

        Essentially equivalent to ``.apply(lambda x: x.head(n))``,
        except ignores as_index flag.

        Examples
        --------

        >>> df = DataFrame([[1, 2], [1, 4], [5, 6]],
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
        obj = self._selected_obj
        in_head = self._cumcount_array() < n
        head = obj[in_head]
        return head

    def tail(self, n=5):
        """
        Returns last n rows of each group

        Essentially equivalent to ``.apply(lambda x: x.tail(n))``,
        except ignores as_index flag.

        Examples
        --------

        >>> df = DataFrame([['a', 1], ['a', 2], ['b', 1], ['b', 2]],
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
        obj = self._selected_obj
        rng = np.arange(0, -self.grouper._max_groupsize, -1, dtype='int64')
        in_tail = self._cumcount_array(rng, ascending=False) > -n
        tail = obj[in_tail]
        return tail

    def _cumcount_array(self, arr=None, ascending=True):
        """
        arr is where cumcount gets its values from

        note: this is currently implementing sort=False (though the default is sort=True)
              for groupby in general
        """
        if arr is None:
            arr = np.arange(self.grouper._max_groupsize, dtype='int64')

        len_index = len(self._selected_obj.index)
        cumcounts = np.zeros(len_index, dtype=arr.dtype)
        if not len_index:
            return cumcounts

        indices, values = [], []
        for v in self.indices.values():
            indices.append(v)

            if ascending:
                values.append(arr[:len(v)])
            else:
                values.append(arr[len(v)-1::-1])

        indices = np.concatenate(indices)
        values = np.concatenate(values)
        cumcounts[indices] = values

        return cumcounts

    def _index_with_as_index(self, b):
        """
        Take boolean mask of index to be returned from apply, if as_index=True

        """
        # TODO perf, it feels like this should already be somewhere...
        from itertools import chain
        original = self._selected_obj.index
        gp = self.grouper
        levels = chain((gp.levels[i][gp.labels[i][b]]
                        for i in range(len(gp.groupings))),
                       (original.get_level_values(i)[b]
                        for i in range(original.nlevels)))
        new = MultiIndex.from_arrays(list(levels))
        new.names = gp.names + original.names
        return new

    def _try_cast(self, result, obj):
        """
        try to cast the result to our obj original type,
        we may have roundtripped thru object in the mean-time

        """
        if obj.ndim > 1:
            dtype = obj.values.dtype
        else:
            dtype = obj.dtype

        if not np.isscalar(result):
            result = _possibly_downcast_to_dtype(result, dtype)

        return result

    def _cython_transform(self, how, numeric_only=True):
        output = {}
        for name, obj in self._iterate_slices():
            is_numeric = is_numeric_dtype(obj.dtype)
            if numeric_only and not is_numeric:
                continue

            try:
                result, names = self.grouper.transform(obj.values, how)
            except AssertionError as e:
                raise GroupByError(str(e))
            output[name] = self._try_cast(result, obj)

        if len(output) == 0:
            raise DataError('No numeric types to aggregate')

        return self._wrap_transformed_output(output, names)

    def _cython_agg_general(self, how, numeric_only=True):
        output = {}
        for name, obj in self._iterate_slices():
            is_numeric = is_numeric_dtype(obj.dtype)
            if numeric_only and not is_numeric:
                continue

            try:
                result, names = self.grouper.aggregate(obj.values, how)
            except AssertionError as e:
                raise GroupByError(str(e))
            output[name] = self._try_cast(result, obj)

        if len(output) == 0:
            raise DataError('No numeric types to aggregate')

        return self._wrap_aggregated_output(output, names)

    def _python_agg_general(self, func, *args, **kwargs):
        func = _intercept_function(func)
        f = lambda x: func(x, *args, **kwargs)

        # iterate through "columns" ex exclusions to populate output dict
        output = {}
        for name, obj in self._iterate_slices():
            try:
                result, counts = self.grouper.agg_series(obj, f)
                output[name] = self._try_cast(result, obj)
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
                    values = com.ensure_float(values)

                output[name] = self._try_cast(values[mask], result)

        return self._wrap_aggregated_output(output)

    def _wrap_applied_output(self, *args, **kwargs):
        raise AbstractMethodError(self)

    def _concat_objects(self, keys, values, not_indexed_same=False):
        from pandas.tools.merge import concat

        if not not_indexed_same:
            result = concat(values, axis=self.axis)
            ax = self._selected_obj._get_axis(self.axis)

            if isinstance(result, Series):
                result = result.reindex(ax)
            else:
                result = result.reindex_axis(ax, axis=self.axis)

        elif self.group_keys:

            if self.as_index:

                # possible MI return case
                group_keys = keys
                group_levels = self.grouper.levels
                group_names = self.grouper.names
                result = concat(values, axis=self.axis, keys=group_keys,
                                levels=group_levels, names=group_names)
            else:

                # GH5610, returns a MI, with the first level being a
                # range index
                keys = list(range(len(values)))
                result = concat(values, axis=self.axis, keys=keys)
        else:
            result = concat(values, axis=self.axis)

        return result

    def _apply_filter(self, indices, dropna):
        if len(indices) == 0:
            indices = []
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


@Appender(GroupBy.__doc__)
def groupby(obj, by, **kwds):
    if isinstance(obj, Series):
        klass = SeriesGroupBy
    elif isinstance(obj, DataFrame):
        klass = DataFrameGroupBy
    else:  # pragma: no cover
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


class BaseGrouper(object):
    """
    This is an internal Grouper class, which actually holds the generated groups
    """

    def __init__(self, axis, groupings, sort=True, group_keys=True):
        self._filter_empty_groups = self.compressed = len(groupings) != 1
        self.axis, self.groupings, self.sort, self.group_keys = \
                axis, groupings, sort, group_keys

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
            mapper = _KeyMapper(comp_ids, ngroups, self.labels, self.levels)
            return [mapper.get_key(i) for i in range(ngroups)]

    def apply(self, f, data, axis=0):
        mutated = False
        splitter = self._get_splitter(data, axis=axis)
        group_keys = self._get_group_keys()

        # oh boy
        f_name = com._get_callable_name(f)
        if (f_name not in _plotting_methods and
                hasattr(splitter, 'fast_apply') and axis == 0):
            try:
                values, mutated = splitter.fast_apply(f, group_keys)
                return group_keys, values, mutated
            except (lib.InvalidApply):
                # we detect a mutation of some kind
                # so take slow path
                pass
            except (Exception) as e:
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
            keys = [_values_from_object(ping.group_index) for ping in self.groupings]
            return _get_indices_dict(label_list, keys)

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
        ids = com._ensure_platform_int(ids)
        out = np.bincount(ids[ids != -1], minlength=ngroup)
        return Series(out, index=self.result_index, dtype='int64')

    @cache_readonly
    def _max_groupsize(self):
        '''
        Compute size of largest group

        '''
        # For many items in each group this is much faster than
        # self.size().max(), in worst case marginally slower
        if self.indices:
            return max(len(v) for v in self.indices.values())
        else:
            return 0

    @cache_readonly
    def groups(self):
        """ dict {group name -> group labels} """
        if len(self.groupings) == 1:
            return self.groupings[0].groups
        else:
            to_groupby = lzip(*(ping.grouper for ping in self.groupings))
            to_groupby = Index(to_groupby)
            return self.axis.groupby(to_groupby.values)

    @cache_readonly
    def is_monotonic(self):
        # return if my group orderings are monotonic
        return Index(self.group_info[0]).is_monotonic

    @cache_readonly
    def group_info(self):
        comp_ids, obs_group_ids = self._get_compressed_labels()

        ngroups = len(obs_group_ids)
        comp_ids = com._ensure_int64(comp_ids)
        return comp_ids, obs_group_ids, ngroups


    def _get_compressed_labels(self):
        all_labels = [ping.labels for ping in self.groupings]
        if len(all_labels) > 1:
            group_index = get_group_index(all_labels, self.shape,
                                          sort=True, xnull=True)
            return _compress_group_index(group_index, sort=self.sort)

        ping = self.groupings[0]
        return ping.labels, np.arange(len(ping.group_index))

    @cache_readonly
    def ngroups(self):
        return len(self.result_index)

    @property
    def recons_labels(self):
        comp_ids, obs_ids, _ = self.group_info
        labels = (ping.labels for ping in self.groupings)
        return decons_obs_group_ids(comp_ids,
                obs_ids, self.shape, labels, xnull=True)

    @cache_readonly
    def result_index(self):
        if not self.compressed and len(self.groupings) == 1:
            return self.groupings[0].group_index.rename(self.names[0])

        return MultiIndex(levels=[ping.group_index for ping in self.groupings],
                          labels=self.recons_labels,
                          verify_integrity=False,
                          names=self.names)

    def get_group_levels(self):
        if not self.compressed and len(self.groupings) == 1:
            return [self.groupings[0].group_index]

        name_list = []
        for ping, labels in zip(self.groupings, self.recons_labels):
            labels = com._ensure_platform_int(labels)
            levels = ping.group_index.take(labels)

            name_list.append(levels)

        return name_list

    #------------------------------------------------------------
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
                'f': lambda func, a, b, c, d: func(a, b, c, d, 1)
            },
            'last': 'group_last',
            },

        'transform': {
            'cumprod' : 'group_cumprod',
            'cumsum' : 'group_cumsum',
            }
    }

    _cython_arity = {
        'ohlc': 4,  # OHLC
    }

    _name_functions = {}

    def _get_cython_function(self, kind, how, values, is_numeric):

        dtype_str = values.dtype.name

        def get_func(fname):
            # see if there is a fused-type version of function
            # only valid for numeric
            f = getattr(_algos, fname, None)
            if f is not None and is_numeric:
                return f

            # otherwise find dtype-specific version, falling back to object
            for dt in [dtype_str, 'object']:
                f = getattr(_algos, "%s_%s" % (fname, dtype_str), None)
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
        return func, dtype_str

    def _cython_operation(self, kind, values, how, axis):
        assert kind in ['transform', 'aggregate']

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

        is_numeric = is_numeric_dtype(values.dtype)

        if is_datetime_or_timedelta_dtype(values.dtype):
            values = values.view('int64')
            is_numeric = True
        elif is_bool_dtype(values.dtype):
            values = _algos.ensure_float64(values)
        elif com.is_integer_dtype(values):
            values = values.astype('int64', copy=False)
        elif is_numeric:
            values = _algos.ensure_float64(values)
        else:
            values = values.astype(object)

        try:
            func, dtype_str = self._get_cython_function(kind, how, values, is_numeric)
        except NotImplementedError:
            if is_numeric:
                values = _algos.ensure_float64(values)
                func, dtype_str = self._get_cython_function(kind, how, values, is_numeric)
            else:
                raise

        if is_numeric:
            out_dtype = '%s%d' % (values.dtype.kind, values.dtype.itemsize)
        else:
            out_dtype = 'object'

        labels, _, _ = self.group_info

        if kind == 'aggregate':
            result = np.empty(out_shape, dtype=out_dtype)
            result.fill(np.nan)
            counts = np.zeros(self.ngroups, dtype=np.int64)
            result = self._aggregate(result, counts, values, labels, func, is_numeric)
        elif kind == 'transform':
            result = np.empty_like(values, dtype=out_dtype)
            result.fill(np.nan)
            # temporary storange for running-total type tranforms
            accum = np.empty(out_shape, dtype=out_dtype)
            result = self._transform(result, accum, values, labels, func, is_numeric)

        if com.is_integer_dtype(result):
            if len(result[result == tslib.iNaT]) > 0:
                result = result.astype('float64')
                result[result == tslib.iNaT] = np.nan

        if kind == 'aggregate' and self._filter_empty_groups and not counts.all():
            if result.ndim == 2:
                try:
                    result = lib.row_bool_subset(
                        result, (counts > 0).view(np.uint8))
                except ValueError:
                    result = lib.row_bool_subset_object(
                                    com._ensure_object(result),
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

    def aggregate(self, values, how, axis=0):
        return self._cython_operation('aggregate', values, how, axis)

    def transform(self, values, how, axis=0):
        return self._cython_operation('transform', values, how, axis)

    def _aggregate(self, result, counts, values, comp_ids, agg_func, is_numeric):
        if values.ndim > 3:
            # punting for now
            raise NotImplementedError("number of dimensions is currently "
                                      "limited to 3")
        elif values.ndim > 2:
            for i, chunk in enumerate(values.transpose(2, 0, 1)):

                chunk = chunk.squeeze()
                agg_func(result[:, :, i], counts, chunk, comp_ids)
        else:
            agg_func(result, counts, values, comp_ids)

        return result

    def _transform(self, result, accum, values, comp_ids, transform_func, is_numeric):
        comp_ids, _, ngroups = self.group_info
        if values.ndim > 3:
            # punting for now
            raise NotImplementedError("number of dimensions is currently "
                                      "limited to 3")
        elif values.ndim > 2:
            for i, chunk in enumerate(values.transpose(2, 0, 1)):

                chunk = chunk.squeeze()
                agg_func(result[:, :, i], values, comp_ids, accum)
        else:
            transform_func(result, values, comp_ids, accum)

        return result

    def agg_series(self, obj, func):
        try:
            return self._aggregate_series_fast(obj, func)
        except Exception:
            return self._aggregate_series_pure_python(obj, func)

    def _aggregate_series_fast(self, obj, func):
        func = _intercept_function(func)

        if obj.index._has_complex_internals:
            raise TypeError('Incompatible index for Cython grouper')

        group_index, _, ngroups = self.group_info

        # avoids object / Series creation overhead
        dummy = obj._get_values(slice(None, 0)).to_dense()
        indexer = _get_group_index_sorter(group_index, ngroups)
        obj = obj.take(indexer, convert=False)
        group_index = com.take_nd(group_index, indexer, allow_fill=False)
        grouper = lib.SeriesGrouper(obj, func, group_index, ngroups,
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
                if (isinstance(res, (Series, Index, np.ndarray)) or
                        isinstance(res, list)):
                    raise ValueError('Function does not reduce')
                result = np.empty(ngroups, dtype='O')

            counts[label] = group.shape[0]
            result[label] = res

        result = lib.maybe_convert_objects(result, try_float=0)
        return result, counts


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

class BinGrouper(BaseGrouper):

    def __init__(self, bins, binlabels, filter_empty=False):
        self.bins = com._ensure_int64(bins)
        self.binlabels = _ensure_index(binlabels)
        self._filter_empty_groups = filter_empty

    @cache_readonly
    def groups(self):
        """ dict {group name -> group labels} """

        # this is mainly for compat
        # GH 3881
        result = {}
        for key, value in zip(self.binlabels, self.bins):
            if key is not tslib.NaT:
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
            slicer = lambda start,edge: data._slice(slice(start,edge),axis=axis)
            length = len(data.axes[axis])
        else:
            slicer = lambda start,edge: data[slice(start,edge)]
            length = len(data)

        start = 0
        for edge, label in zip(self.bins, self.binlabels):
            if label is not tslib.NaT:
                yield label, slicer(start,edge)
            start = edge

        if start < length:
            yield self.binlabels[-1], slicer(start,None)

    def apply(self, f, data, axis=0):
        result_keys = []
        result_values = []
        mutated = False
        for key, group in self.get_iterator(data, axis=axis):
            object.__setattr__(group, 'name', key)

            # group might be modified
            group_axes = _get_axes(group)
            res = f(group)

            if not _is_indexed_like(res, group_axes):
                mutated = True

            result_keys.append(key)
            result_values.append(res)

        return result_keys, result_values, mutated

    @cache_readonly
    def indices(self):
        indices = collections.defaultdict(list)

        i = 0
        for label, bin in zip(self.binlabels, self.bins):
            if i < bin:
                if label is not tslib.NaT:
                    indices[label] = list(range(i, bin))
                i = bin
        return indices

    @cache_readonly
    def group_info(self):
        ngroups = self.ngroups
        obs_group_ids = np.arange(ngroups)
        rep = np.diff(np.r_[0, self.bins])

        rep = com._ensure_platform_int(rep)
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
        if len(self.binlabels) != 0 and isnull(self.binlabels[0]):
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
        # for compat
        return None

    #----------------------------------------------------------------------
    # cython aggregation

    _cython_functions = copy.deepcopy(BaseGrouper._cython_functions)
    _cython_functions['aggregate']['ohlc'] = 'group_ohlc'
    _cython_functions['aggregate'].pop('median')

    _name_functions = {
        'ohlc': lambda *args: ['open', 'high', 'low', 'close']
    }

    def agg_series(self, obj, func):
        dummy = obj[:0]
        grouper = lib.SeriesBinGrouper(obj, func, self.bins, dummy)
        return grouper.get_result()


class Grouping(object):

    """
    Holds the grouping information for a single key

    Parameters
    ----------
    index : Index
    grouper :
    obj :
    name :
    level :
    in_axis : if the Grouping is a column in self.obj and hence among
        Groupby.exclusions list

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

    def __init__(self, index, grouper=None, obj=None, name=None, level=None,
                 sort=True, in_axis=False):

        self.name = name
        self.level = level
        self.grouper = _convert_grouper(index, grouper)
        self.index = index
        self.sort = sort
        self.obj = obj
        self.in_axis = in_axis

        # right place for this?
        if isinstance(grouper, (Series, Index)) and name is None:
            self.name = grouper.name

        if isinstance(grouper, MultiIndex):
            self.grouper = grouper.values

        # pre-computed
        self._should_compress = True

        # we have a single grouper which may be a myriad of things, some of which are
        # dependent on the passing in level
        #

        if level is not None:
            if not isinstance(level, int):
                if level not in index.names:
                    raise AssertionError('Level %s not in index' % str(level))
                level = index.names.index(level)

            inds = index.labels[level]
            level_index = index.levels[level]

            if self.name is None:
                self.name = index.names[level]

            # XXX complete hack

            if grouper is not None:
                level_values = index.levels[level].take(inds)
                self.grouper = level_values.map(self.grouper)
            else:
                # all levels may not be observed
                labels, uniques = algos.factorize(inds, sort=True)

                if len(uniques) > 0 and uniques[0] == -1:
                    # handle NAs
                    mask = inds != -1
                    ok_labels, uniques = algos.factorize(inds[mask], sort=True)

                    labels = np.empty(len(inds), dtype=inds.dtype)
                    labels[mask] = ok_labels
                    labels[~mask] = -1

                if len(uniques) < len(level_index):
                    level_index = level_index.take(uniques)

                self._labels = labels
                self._group_index = level_index
                self.grouper = level_index.take(labels)
        else:
            if isinstance(self.grouper, (list, tuple)):
                self.grouper = com._asarray_tuplesafe(self.grouper)

            # a passed Categorical
            elif is_categorical_dtype(self.grouper):

                # must have an ordered categorical
                if self.sort:
                    if not self.grouper.ordered:

                        # technically we cannot group on an unordered Categorical
                        # but this a user convenience to do so; the ordering
                        # is preserved and if it's a reduction it doesn't make any difference
                        pass

                # fix bug #GH8868 sort=False being ignored in categorical groupby
                else:
                    cat = self.grouper.unique()
                    self.grouper = self.grouper.reorder_categories(cat.categories)

                # we make a CategoricalIndex out of the cat grouper
                # preserving the categories / ordered attributes
                self._labels = self.grouper.codes

                c = self.grouper.categories
                self._group_index = CategoricalIndex(Categorical.from_codes(np.arange(len(c)),
                                                     categories=c,
                                                     ordered=self.grouper.ordered))

            # a passed Grouper like
            elif isinstance(self.grouper, Grouper):

                # get the new grouper
                grouper = self.grouper._get_binner_for_grouping(self.obj)
                self.obj = self.grouper.obj
                self.grouper = grouper
                if self.name is None:
                    self.name = grouper.name

            # we are done
            if isinstance(self.grouper, Grouping):
                self.grouper = self.grouper.grouper

            # no level passed
            elif not isinstance(self.grouper, (Series, Index, Categorical, np.ndarray)):
                if getattr(self.grouper, 'ndim', 1) != 1:
                    t = self.name or str(type(self.grouper))
                    raise ValueError("Grouper for '%s' not 1-dimensional" % t)
                self.grouper = self.index.map(self.grouper)
                if not (hasattr(self.grouper, "__len__") and
                        len(self.grouper) == len(self.index)):
                    errmsg = ('Grouper result violates len(labels) == '
                              'len(data)\nresult: %s' %
                              com.pprint_thing(self.grouper))
                    self.grouper = None  # Try for sanity
                    raise AssertionError(errmsg)

        # if we have a date/time-like grouper, make sure that we have Timestamps like
        if getattr(self.grouper,'dtype',None) is not None:
            if is_datetime64_dtype(self.grouper):
                from pandas import to_datetime
                self.grouper = to_datetime(self.grouper)
            elif is_timedelta64_dtype(self.grouper):
                from pandas import to_timedelta
                self.grouper = to_timedelta(self.grouper)

    def __repr__(self):
        return 'Grouping(%s)' % self.name

    def __iter__(self):
        return iter(self.indices)

    _labels = None
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
    def group_index(self):
        if self._group_index is None:
            self._make_labels()
        return self._group_index

    def _make_labels(self):
        if self._labels is None or self._group_index is None:
            labels, uniques = algos.factorize(self.grouper, sort=self.sort)
            uniques = Index(uniques, name=self.name)
            self._labels = labels
            self._group_index = uniques

    @cache_readonly
    def groups(self):
        return self.index.groupby(self.grouper)

def _get_grouper(obj, key=None, axis=0, level=None, sort=True):
    """
    create and return a BaseGrouper, which is an internal
    mapping of how to create the grouper indexers.
    This may be composed of multiple Grouping objects, indicating
    multiple groupers

    Groupers are ultimately index mappings. They can originate as:
    index mappings, keys to columns, functions, or Groupers

    Groupers enable local references to axis,level,sort, while
    the passed in axis, level, and sort are 'global'.

    This routine tries to figure out what the passing in references
    are and then creates a Grouping for each one, combined into
    a BaseGrouper.

    """

    group_axis = obj._get_axis(axis)

    # validate that the passed level is compatible with the passed
    # axis of the object
    if level is not None:
        if not isinstance(group_axis, MultiIndex):
            if isinstance(level, compat.string_types):
                if obj.index.name != level:
                    raise ValueError('level name %s is not the name of the '
                                     'index' % level)
            elif level > 0:
                raise ValueError('level > 0 only valid with MultiIndex')

            level = None
            key = group_axis

    # a passed-in Grouper, directly convert
    if isinstance(key, Grouper):
        binner, grouper, obj = key._get_grouper(obj)
        if key.key is None:
            return grouper, [], obj
        else:
            return grouper, set([key.key]), obj

    # already have a BaseGrouper, just return it
    elif isinstance(key, BaseGrouper):
        return key, [], obj

    if not isinstance(key, (tuple, list)):
        keys = [key]
    else:
        keys = key

    # what are we after, exactly?
    match_axis_length = len(keys) == len(group_axis)
    any_callable = any(callable(g) or isinstance(g, dict) for g in keys)
    any_groupers = any(isinstance(g, Grouper) for g in keys)
    any_arraylike = any(isinstance(g, (list, tuple, Series, Index, np.ndarray))
                        for g in keys)

    try:
        if isinstance(obj, DataFrame):
            all_in_columns = all(g in obj.columns for g in keys)
        else:
            all_in_columns = False
    except Exception:
        all_in_columns = False

    if (not any_callable and not all_in_columns
        and not any_arraylike and not any_groupers
            and match_axis_length
            and level is None):
        keys = [com._asarray_tuplesafe(keys)]

    if isinstance(level, (tuple, list)):
        if key is None:
            keys = [None] * len(level)
        levels = level
    else:
        levels = [level] * len(keys)

    groupings = []
    exclusions = []

    # if the actual grouper should be obj[key]
    def is_in_axis(key):
        if not _is_label_like(key):
            try:
                obj._data.items.get_loc(key)
            except Exception:
                return False

        return True

    # if the the grouper is obj[name]
    def is_in_obj(gpr):
        try:
            return id(gpr) == id(obj[gpr.name])
        except Exception:
            return False

    for i, (gpr, level) in enumerate(zip(keys, levels)):

        if is_in_obj(gpr):  # df.groupby(df['name'])
            in_axis, name = True, gpr.name
            exclusions.append(name)

        elif is_in_axis(gpr):  # df.groupby('name')
            in_axis, name, gpr = True, gpr, obj[gpr]
            exclusions.append(name)

        else:
            in_axis, name = False, None

        if is_categorical_dtype(gpr) and len(gpr) != len(obj):
            raise ValueError("Categorical dtype grouper must have len(grouper) == len(data)")

        ping = Grouping(group_axis, gpr, obj=obj, name=name,
                        level=level, sort=sort, in_axis=in_axis)

        groupings.append(ping)

    if len(groupings) == 0:
        raise ValueError('No group keys passed!')

    # create the internals grouper
    grouper = BaseGrouper(group_axis, groupings, sort=sort)

    return grouper, exclusions, obj


def _is_label_like(val):
    return isinstance(val, compat.string_types) or np.isscalar(val)


def _convert_grouper(axis, grouper):
    if isinstance(grouper, dict):
        return grouper.get
    elif isinstance(grouper, Series):
        if grouper.index.equals(axis):
            return grouper._values
        else:
            return grouper.reindex(axis)._values
    elif isinstance(grouper, (list, Series, Index, np.ndarray)):
        if len(grouper) != len(axis):
            raise AssertionError('Grouper and axis must be same length')
        return grouper
    else:
        return grouper

def _whitelist_method_generator(klass, whitelist) :
    """
    Yields all GroupBy member defs for DataFrame/Series names in _whitelist.

    Parameters
    ----------
    klass - class where members are defined.  Should be Series or DataFrame

    whitelist - list of names of klass methods to be constructed

    Returns
    -------
    The generator yields a sequence of strings, each suitable for exec'ing,
    that define implementations of the named methods for DataFrameGroupBy
    or SeriesGroupBy.

    Since we don't want to override methods explicitly defined in the
    base class, any such name is skipped.
    """

    method_wrapper_template = \
    """def %(name)s(%(sig)s) :
    \"""
    %(doc)s
    \"""
    f = %(self)s.__getattr__('%(name)s')
    return f(%(args)s)"""
    property_wrapper_template = \
    """@property
def %(name)s(self) :
    \"""
    %(doc)s
    \"""
    return self.__getattr__('%(name)s')"""
    for name in whitelist :
        # don't override anything that was explicitly defined
        # in the base class
        if hasattr(GroupBy,name) :
            continue
        # ugly, but we need the name string itself in the method.
        f = getattr(klass,name)
        doc = f.__doc__
        doc = doc if type(doc)==str else ''
        if type(f) == types.MethodType :
            wrapper_template = method_wrapper_template
            decl, args = make_signature(f)
            # pass args by name to f because otherwise
            # GroupBy._make_wrapper won't know whether
            # we passed in an axis parameter.
            args_by_name = ['{0}={0}'.format(arg) for arg in args[1:]]
            params = {'name':name,
                      'doc':doc,
                      'sig':','.join(decl),
                      'self':args[0],
                      'args':','.join(args_by_name)}
        else :
            wrapper_template = property_wrapper_template
            params = {'name':name, 'doc':doc}
        yield wrapper_template % params

class SeriesGroupBy(GroupBy):
    #
    # Make class defs of attributes on SeriesGroupBy whitelist
    _apply_whitelist = _series_apply_whitelist
    for _def_str in _whitelist_method_generator(Series,_series_apply_whitelist) :
        exec(_def_str)

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

        Examples
        --------
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
        if isinstance(func_or_funcs, compat.string_types):
            return getattr(self, func_or_funcs)(*args, **kwargs)

        if hasattr(func_or_funcs, '__iter__'):
            ret = self._aggregate_multiple_funcs(func_or_funcs)
        else:
            cyfunc = _intercept_cython(func_or_funcs)
            if cyfunc and not args and not kwargs:
                return getattr(self, cyfunc)()

            if self.grouper.nkeys > 1:
                return self._python_agg_general(func_or_funcs, *args, **kwargs)

            try:
                return self._python_agg_general(func_or_funcs, *args, **kwargs)
            except Exception:
                result = self._aggregate_named(func_or_funcs, *args, **kwargs)

            index = Index(sorted(result), name=self.grouper.names[0])
            ret = Series(result, index=index)

        if not self.as_index:  # pragma: no cover
            print('Warning, ignoring as_index=True')

        return ret

    def _aggregate_multiple_funcs(self, arg):
        if isinstance(arg, dict):
            columns = list(arg.keys())
            arg = list(arg.items())
        elif any(isinstance(x, (tuple, list)) for x in arg):
            arg = [(x, x) if not isinstance(x, (tuple, list)) else x
                   for x in arg]

            # indicated column order
            columns = lzip(*arg)[0]
        else:
            # list of functions / function names
            columns = []
            for f in arg:
                if isinstance(f, compat.string_types):
                    columns.append(f)
                else:
                    # protect against callables without names
                    columns.append(com._get_callable_name(f))
            arg = lzip(columns, arg)

        results = {}
        for name, func in arg:
            if name in results:
                raise SpecificationError('Function names must be unique, '
                                         'found multiple named %s' % name)

            results[name] = self.aggregate(func)

        return DataFrame(results, columns=columns)

    def _wrap_output(self, output, index, names=None):
        """ common agg/transform wrapping logic """
        output = output[self.name]

        if names is not None:
            return DataFrame(output, index=index, columns=names)
        else:
            name = self.name
            if name is None:
                name = self._selected_obj.name
            return Series(output, index=index, name=name)

    def _wrap_aggregated_output(self, output, names=None):
        return self._wrap_output(output=output,
                                 index=self.grouper.result_index,
                                 names=names)

    def _wrap_transformed_output(self, output, names=None):
        return self._wrap_output(output=output,
                                 index=self.obj.index,
                                 names=names)

    def _wrap_applied_output(self, keys, values, not_indexed_same=False):
        if len(keys) == 0:
            # GH #6265
            return Series([], name=self.name)

        def _get_index():
            if self.grouper.nkeys > 1:
                index = MultiIndex.from_tuples(keys, names=self.grouper.names)
            else:
                index = Index(keys, name=self.grouper.names[0])
            return index

        if isinstance(values[0], dict):
            # GH #823
            index = _get_index()
            return DataFrame(values, index=index).stack()

        if isinstance(values[0], (Series, dict)):
            return self._concat_objects(keys, values,
                                        not_indexed_same=not_indexed_same)
        elif isinstance(values[0], DataFrame):
            # possible that Series -> DataFrame by applied function
            return self._concat_objects(keys, values,
                                        not_indexed_same=not_indexed_same)
        else:
            # GH #6265
            return Series(values, index=_get_index(), name=self.name)

    def _aggregate_named(self, func, *args, **kwargs):
        result = {}

        for name, group in self:
            group.name = name
            output = func(group, *args, **kwargs)
            if isinstance(output, (Series, Index, np.ndarray)):
                raise Exception('Must produce aggregated value')
            result[name] = self._try_cast(output, group)

        return result

    def transform(self, func, *args, **kwargs):
        """
        Call function producing a like-indexed Series on each group and return
        a Series with the transformed values

        Parameters
        ----------
        func : function
            To apply to each group. Should return a Series with the same index

        Examples
        --------
        >>> grouped.transform(lambda x: (x - x.mean()) / x.std())

        Returns
        -------
        transformed : Series
        """

        func = _intercept_cython(func) or func

        # if string function
        if isinstance(func, compat.string_types):
            if func in _cython_transforms:
                # cythonized transform
                return getattr(self, func)(*args, **kwargs)
            else:
                # cythonized aggregation and merge
                return self._transform_fast(lambda : getattr(self, func)(*args, **kwargs))

        # reg transform
        dtype = self._selected_obj.dtype
        result = self._selected_obj.values.copy()

        wrapper = lambda x: func(x, *args, **kwargs)
        for i, (name, group) in enumerate(self):
            object.__setattr__(group, 'name', name)
            res = wrapper(group)

            if hasattr(res, 'values'):
                res = res.values

            # may need to astype
            try:
                common_type = np.common_type(np.array(res), result)
                if common_type != result.dtype:
                    result = result.astype(common_type)
            except:
                pass

            indexer = self._get_index(name)
            result[indexer] = res

        result = _possibly_downcast_to_dtype(result, dtype)
        return self._selected_obj.__class__(result,
                                            index=self._selected_obj.index,
                                            name=self._selected_obj.name)

    def _transform_fast(self, func):
        """
        fast version of transform, only applicable to builtin/cythonizable functions
        """
        if isinstance(func, compat.string_types):
            func = getattr(self,func)

        ids, _, ngroup = self.grouper.group_info
        mask = ids != -1

        out = func().values[ids]
        if not mask.all():
            out = np.where(mask, out, np.nan)

        obs = np.zeros(ngroup, dtype='bool')
        obs[ids[mask]] = True
        if not obs.all():
            out = self._try_cast(out, self._selected_obj)

        return Series(out, index=self.obj.index)

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
        >>> grouped.filter(lambda x: x.mean() > 0)

        Returns
        -------
        filtered : Series
        """
        if isinstance(func, compat.string_types):
            wrapper = lambda x: getattr(x, func)(*args, **kwargs)
        else:
            wrapper = lambda x: func(x, *args, **kwargs)

        # Interpret np.nan as False.
        def true_and_notnull(x, *args, **kwargs):
            b = wrapper(x, *args, **kwargs)
            return b and notnull(b)

        try:
            indices = [self._get_index(name) for name, group in self
                       if true_and_notnull(group)]
        except ValueError:
            raise TypeError("the filter must return a boolean result")
        except TypeError:
            raise TypeError("the filter must return a boolean result")

        filtered = self._apply_filter(indices, dropna)
        return filtered

    def nunique(self, dropna=True):
        ids, _, _ = self.grouper.group_info
        val = self.obj.get_values()

        try:
            sorter = np.lexsort((val, ids))
        except TypeError:  # catches object dtypes
            assert val.dtype == object, \
                'val.dtype must be object, got %s' % val.dtype
            val, _ = algos.factorize(val, sort=False)
            sorter = np.lexsort((val, ids))
            isnull = lambda a: a == -1
        else:
            isnull = com.isnull

        ids, val = ids[sorter], val[sorter]

        # group boundaries are where group ids change
        # unique observations are where sorted values change
        idx = np.r_[0, 1 + np.nonzero(ids[1:] != ids[:-1])[0]]
        inc = np.r_[1, val[1:] != val[:-1]]

        # 1st item of each group is a new unique observation
        mask = isnull(val)
        if dropna:
            inc[idx] = 1
            inc[mask] = 0
        else:
            inc[mask & np.r_[False, mask[:-1]]] = 0
            inc[idx] = 1

        out = np.add.reduceat(inc, idx).astype('int64', copy=False)
        return Series(out if ids[0] != -1 else out[1:],
                      index=self.grouper.result_index,
                      name=self.name)

    @deprecate_kwarg('take_last', 'keep', mapping={True: 'last', False: 'first'})
    @Appender(Series.nlargest.__doc__)
    def nlargest(self, n=5, keep='first'):
        # ToDo: When we remove deprecate_kwargs, we can remote these methods
        # and include nlargest and nsmallest to _series_apply_whitelist
        return self.apply(lambda x: x.nlargest(n=n, keep=keep))


    @deprecate_kwarg('take_last', 'keep', mapping={True: 'last', False: 'first'})
    @Appender(Series.nsmallest.__doc__)
    def nsmallest(self, n=5, keep='first'):
        return self.apply(lambda x: x.nsmallest(n=n, keep=keep))

    def value_counts(self, normalize=False, sort=True, ascending=False,
                     bins=None, dropna=True):

        from functools import partial
        from pandas.tools.tile import cut
        from pandas.tools.merge import _get_join_indexers

        if bins is not None and not np.iterable(bins):
            # scalar bins cannot be done at top level
            # in a backward compatible way
            return self.apply(Series.value_counts,
                              normalize=normalize,
                              sort=sort,
                              ascending=ascending,
                              bins=bins)

        ids, _, _ = self.grouper.group_info
        val = self.obj.get_values()

        # groupby removes null keys from groupings
        mask = ids != -1
        ids, val = ids[mask], val[mask]

        if bins is None:
            lab, lev = algos.factorize(val, sort=True)
        else:
            cat, bins = cut(val, bins, retbins=True)
            # bins[:-1] for backward compat;
            # o.w. cat.categories could be better
            lab, lev, dropna = cat.codes, bins[:-1], False

        sorter = np.lexsort((lab, ids))
        ids, lab = ids[sorter], lab[sorter]

        # group boundaries are where group ids change
        idx = np.r_[0, 1 + np.nonzero(ids[1:] != ids[:-1])[0]]

        # new values are where sorted labels change
        inc = np.r_[True, lab[1:] != lab[:-1]]
        inc[idx] = True  # group boundaries are also new values
        out = np.diff(np.nonzero(np.r_[inc, True])[0]) # value counts

        # num. of times each group should be repeated
        rep = partial(np.repeat, repeats=np.add.reduceat(inc, idx))

        # multi-index components
        labels = list(map(rep, self.grouper.recons_labels)) + [lab[inc]]
        levels = [ping.group_index for ping in self.grouper.groupings] + [lev]
        names = self.grouper.names + [self.name]

        if dropna:
            mask = labels[-1] != -1
            if mask.all():
                dropna = False
            else:
                out, labels = out[mask], [label[mask] for label in labels]

        if normalize:
            out = out.astype('float')
            acc = rep(np.diff(np.r_[idx, len(ids)]))
            out /= acc[mask] if dropna else acc

        if sort and bins is None:
            cat = ids[inc][mask] if dropna else ids[inc]
            sorter = np.lexsort((out if ascending else -out, cat))
            out, labels[-1] = out[sorter], labels[-1][sorter]

        if bins is None:
            mi = MultiIndex(levels=levels, labels=labels, names=names,
                            verify_integrity=False)

            if com.is_integer_dtype(out):
                out = com._ensure_int64(out)
            return Series(out, index=mi)

        # for compat. with algos.value_counts need to ensure every
        # bin is present at every index level, null filled with zeros
        diff = np.zeros(len(out), dtype='bool')
        for lab in labels[:-1]:
            diff |= np.r_[True, lab[1:] != lab[:-1]]

        ncat, nbin = diff.sum(), len(levels[-1])

        left = [np.repeat(np.arange(ncat), nbin),
                np.tile(np.arange(nbin), ncat)]

        right = [diff.cumsum() - 1, labels[-1]]

        _, idx = _get_join_indexers(left, right, sort=False, how='left')
        out = np.where(idx != -1, out[idx], 0)

        if sort:
            sorter = np.lexsort((out if ascending else -out, left[0]))
            out, left[-1] = out[sorter], left[-1][sorter]

        # build the multi-index w/ full levels
        labels = list(map(lambda lab: np.repeat(lab[diff], nbin), labels[:-1]))
        labels.append(left[-1])

        mi = MultiIndex(levels=levels, labels=labels, names=names,
                        verify_integrity=False)

        if com.is_integer_dtype(out):
            out = com._ensure_int64(out)
        return Series(out, index=mi)

    def count(self):
        """ Compute count of group, excluding missing values """
        ids, _, ngroups = self.grouper.group_info
        val = self.obj.get_values()

        mask = (ids != -1) & ~isnull(val)
        ids = com._ensure_platform_int(ids)
        out = np.bincount(ids[mask], minlength=ngroups) if ngroups != 0 else []

        return Series(out, index=self.grouper.result_index, name=self.name, dtype='int64')

    def _apply_to_column_groupbys(self, func):
        """ return a pass thru """
        return func(self)

class NDFrameGroupBy(GroupBy):

    def _iterate_slices(self):
        if self.axis == 0:
            # kludge
            if self._selection is None:
                slice_axis = self.obj.columns
            else:
                slice_axis = self._selection_list
            slicer = lambda x: self.obj[x]
        else:
            slice_axis = self.obj.index
            slicer = self.obj.xs

        for val in slice_axis:
            if val in self.exclusions:
                continue
            yield val, slicer(val)

    def _cython_agg_general(self, how, numeric_only=True):
        new_items, new_blocks = self._cython_agg_blocks(how, numeric_only=numeric_only)
        return self._wrap_agged_blocks(new_items, new_blocks)

    def _wrap_agged_blocks(self, items, blocks):
        obj = self._obj_with_exclusions

        new_axes = list(obj._data.axes)

        # more kludge
        if self.axis == 0:
            new_axes[0], new_axes[1] = new_axes[1], self.grouper.result_index
        else:
            new_axes[self.axis] = self.grouper.result_index

        # Make sure block manager integrity check passes.
        assert new_axes[0].equals(items)
        new_axes[0] = items

        mgr = BlockManager(blocks, new_axes)

        new_obj = type(obj)(mgr)

        return self._post_process_cython_aggregate(new_obj)

    _block_agg_axis = 0

    def _cython_agg_blocks(self, how, numeric_only=True):
        data, agg_axis = self._get_data_to_aggregate()

        new_blocks = []

        if numeric_only:
            data = data.get_numeric_data(copy=False)

        for block in data.blocks:

            result, _ = self.grouper.aggregate(block.values, how, axis=agg_axis)

            # see if we can cast the block back to the original dtype
            result = block._try_coerce_and_cast_result(result)

            newb = make_block(result, placement=block.mgr_locs)
            new_blocks.append(newb)

        if len(new_blocks) == 0:
            raise DataError('No numeric types to aggregate')

        return data.items, new_blocks

    def _get_data_to_aggregate(self):
        obj = self._obj_with_exclusions
        if self.axis == 0:
            return obj.swapaxes(0, 1)._data, 1
        else:
            return obj._data, self.axis

    def _post_process_cython_aggregate(self, obj):
        # undoing kludge from below
        if self.axis == 0:
            obj = obj.swapaxes(0, 1)
        return obj

    @cache_readonly
    def _obj_with_exclusions(self):
        if self._selection is not None:
            return self.obj.reindex(columns=self._selection_list)

        if len(self.exclusions) > 0:
            return self.obj.drop(self.exclusions, axis=1)
        else:
            return self.obj

    @Appender(_agg_doc)
    def aggregate(self, arg, *args, **kwargs):
        if isinstance(arg, compat.string_types):
            return getattr(self, arg)(*args, **kwargs)

        result = OrderedDict()
        if isinstance(arg, dict):
            if self.axis != 0:  # pragma: no cover
                raise ValueError('Can only pass dict with axis=0')

            obj = self._selected_obj

            if any(isinstance(x, (list, tuple, dict)) for x in arg.values()):
                new_arg = OrderedDict()
                for k, v in compat.iteritems(arg):
                    if not isinstance(v, (tuple, list, dict)):
                        new_arg[k] = [v]
                    else:
                        new_arg[k] = v
                arg = new_arg

            keys = []
            if self._selection is not None:
                subset = obj
                if isinstance(subset, DataFrame):
                    raise NotImplementedError("Aggregating on a DataFrame is "
                                              "not supported")

                for fname, agg_how in compat.iteritems(arg):
                    colg = SeriesGroupBy(subset, selection=self._selection,
                                         grouper=self.grouper)
                    result[fname] = colg.aggregate(agg_how)
                    keys.append(fname)
            else:
                for col, agg_how in compat.iteritems(arg):
                    colg = SeriesGroupBy(obj[col], selection=col,
                                         grouper=self.grouper)
                    result[col] = colg.aggregate(agg_how)
                    keys.append(col)

            if isinstance(list(result.values())[0], DataFrame):
                from pandas.tools.merge import concat
                result = concat([result[k] for k in keys], keys=keys, axis=1)
            else:
                result = DataFrame(result)
        elif isinstance(arg, list):
            return self._aggregate_multiple_funcs(arg)
        else:
            cyfunc = _intercept_cython(arg)
            if cyfunc and not args and not kwargs:
                return getattr(self, cyfunc)()

            if self.grouper.nkeys > 1:
                return self._python_agg_general(arg, *args, **kwargs)
            else:

                # try to treat as if we are passing a list
                try:
                    assert not args and not kwargs
                    result = self._aggregate_multiple_funcs([arg])
                    result.columns = Index(result.columns.levels[0],
                                           name=self._selected_obj.columns.name)
                except:
                    result = self._aggregate_generic(arg, *args, **kwargs)

        if not self.as_index:
            self._insert_inaxis_grouper_inplace(result)
            result.index = np.arange(len(result))

        return result._convert(datetime=True)

    def _aggregate_multiple_funcs(self, arg):
        from pandas.tools.merge import concat

        if self.axis != 0:
            raise NotImplementedError("axis other than 0 is not supported")

        obj = self._obj_with_exclusions

        results = []
        keys = []
        for col in obj:
            try:
                colg = SeriesGroupBy(obj[col], selection=col,
                                     grouper=self.grouper)
                results.append(colg.aggregate(arg))
                keys.append(col)
            except (TypeError, DataError):
                pass
            except SpecificationError:
                raise
        result = concat(results, keys=keys, axis=1)

        return result

    def _aggregate_generic(self, func, *args, **kwargs):
        if self.grouper.nkeys != 1:
            raise AssertionError('Number of keys must be 1')

        axis = self.axis
        obj = self._obj_with_exclusions

        result = {}
        if axis != obj._info_axis_number:
            try:
                for name, data in self:
                    result[name] = self._try_cast(func(data, *args, **kwargs),
                                                  data)
            except Exception:
                return self._aggregate_item_by_item(func, *args, **kwargs)
        else:
            for name in self.indices:
                try:
                    data = self.get_group(name, obj=obj)
                    result[name] = self._try_cast(func(data, *args, **kwargs),
                                                  data)
                except Exception:
                    wrapper = lambda x: func(x, *args, **kwargs)
                    result[name] = data.apply(wrapper, axis=axis)

        return self._wrap_generic_output(result, obj)

    def _wrap_aggregated_output(self, output, names=None):
        raise AbstractMethodError(self)

    def _aggregate_item_by_item(self, func, *args, **kwargs):
        # only for axis==0

        obj = self._obj_with_exclusions
        result = {}
        cannot_agg = []
        errors=None
        for item in obj:
            try:
                data = obj[item]
                colg = SeriesGroupBy(data, selection=item,
                                     grouper=self.grouper)
                result[item] = self._try_cast(
                    colg.aggregate(func, *args, **kwargs), data)
            except ValueError:
                cannot_agg.append(item)
                continue
            except TypeError as e:
                cannot_agg.append(item)
                errors=e
                continue

        result_columns = obj.columns
        if cannot_agg:
            result_columns = result_columns.drop(cannot_agg)

            # GH6337
            if not len(result_columns) and errors is not None:
                raise errors

        return DataFrame(result, columns=result_columns)

    def _decide_output_index(self, output, labels):
        if len(output) == len(labels):
            output_keys = labels
        else:
            output_keys = sorted(output)
            try:
                output_keys.sort()
            except Exception:  # pragma: no cover
                pass

            if isinstance(labels, MultiIndex):
                output_keys = MultiIndex.from_tuples(output_keys,
                                                     names=labels.names)

        return output_keys

    def _wrap_applied_output(self, keys, values, not_indexed_same=False):
        from pandas.core.index import _all_indexes_same

        if len(keys) == 0:
            # XXX
            return DataFrame({})

        key_names = self.grouper.names

        if isinstance(values[0], DataFrame):
            return self._concat_objects(keys, values,
                                        not_indexed_same=not_indexed_same)
        elif self.grouper.groupings is not None:
            if len(self.grouper.groupings) > 1:
                key_index = MultiIndex.from_tuples(keys, names=key_names)

            else:
                ping = self.grouper.groupings[0]
                if len(keys) == ping.ngroups:
                    key_index = ping.group_index
                    key_index.name = key_names[0]

                    key_lookup = Index(keys)
                    indexer = key_lookup.get_indexer(key_index)

                    # reorder the values
                    values = [values[i] for i in indexer]
                else:

                    key_index = Index(keys, name=key_names[0])

                # don't use the key indexer
                if not self.as_index:
                    key_index = None

            # make Nones an empty object
            if com._count_not_none(*values) != len(values):
                try:
                    v = next(v for v in values if v is not None)
                except StopIteration:
                    # If all values are None, then this will throw an error.
                    # We'd prefer it return an empty dataframe.
                    return DataFrame()
                if v is None:
                    return DataFrame()
                elif isinstance(v, NDFrame):
                    values = [
                        x if x is not None else
                        v._constructor(**v._construct_axes_dict())
                        for x in values
                        ]

            v = values[0]

            if isinstance(v, (np.ndarray, Index, Series)):
                if isinstance(v, Series):
                    applied_index = self._selected_obj._get_axis(self.axis)
                    all_indexed_same = _all_indexes_same([
                        x.index for x in values
                    ])
                    singular_series = (len(values) == 1 and
                                       applied_index.nlevels == 1)

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
                            from pandas.tools.merge import concat
                            return concat(values)

                    if not all_indexed_same:
                        # GH 8467
                        return self._concat_objects(
                            keys, values, not_indexed_same=True,
                        )

                try:
                    if self.axis == 0:
                        # GH6124 if the list of Series have a consistent name,
                        # then propagate that name to the result.
                        index = v.index.copy()
                        if index.name is None:
                            # Only propagate the series name to the result
                            # if all series have a consistent name.  If the
                            # series do not have a consistent name, do
                            # nothing.
                            names = set(v.name for v in values)
                            if len(names) == 1:
                                index.name = list(names)[0]

                        # normally use vstack as its faster than concat
                        # and if we have mi-columns
                        if isinstance(v.index, MultiIndex) or key_index is None:
                            stacked_values = np.vstack(map(np.asarray, values))
                            result = DataFrame(stacked_values, index=key_index,
                                               columns=index)
                        else:
                            # GH5788 instead of stacking; concat gets the dtypes correct
                            from pandas.tools.merge import concat
                            result = concat(values, keys=key_index,
                                            names=key_index.names,
                                            axis=self.axis).unstack()
                            result.columns = index
                    else:
                        stacked_values = np.vstack(map(np.asarray, values))
                        result = DataFrame(stacked_values.T, index=v.index,
                                           columns=key_index)

                except (ValueError, AttributeError):
                    # GH1738: values is list of arrays of unequal lengths fall
                    # through to the outer else caluse
                    return Series(values, index=key_index)

                # if we have date/time like in the original, then coerce dates
                # as we are stacking can easily have object dtypes here
                if (self._selected_obj.ndim == 2 and
                        self._selected_obj.dtypes.isin(_DATELIKE_DTYPES).any()):
                    result = result._convert(numeric=True)
                    date_cols = self._selected_obj.select_dtypes(
                        include=list(_DATELIKE_DTYPES)).columns
                    date_cols = date_cols.intersection(result.columns)
                    result[date_cols] = (result[date_cols]
                                         ._convert(datetime=True,
                                                          coerce=True))
                else:
                    result = result._convert(datetime=True)

                return self._reindex_output(result)

            else:
                # only coerce dates if we find at least 1 datetime
                coerce = True if any([ isinstance(v,Timestamp) for v in values ]) else False
                return (Series(values, index=key_index)
                        ._convert(datetime=True,
                                         coerce=coerce))

        else:
            # Handle cases like BinGrouper
            return self._concat_objects(keys, values,
                                        not_indexed_same=not_indexed_same)

    def _transform_general(self, func, *args, **kwargs):
        from pandas.tools.merge import concat

        applied = []

        obj = self._obj_with_exclusions
        gen = self.grouper.get_iterator(obj, axis=self.axis)
        fast_path, slow_path = self._define_paths(func, *args, **kwargs)

        path = None
        for name, group in gen:
            object.__setattr__(group, 'name', name)

            if path is None:
                # Try slow path and fast path.
                try:
                    path, res = self._choose_path(fast_path, slow_path, group)
                except TypeError:
                    return self._transform_item_by_item(obj, fast_path)
                except Exception:  # pragma: no cover
                    res = fast_path(group)
                    path = fast_path
            else:
                res = path(group)

            # broadcasting
            if isinstance(res, Series):
                if res.index.is_(obj.index):
                    group.T.values[:] = res
                else:
                    group.values[:] = res

                applied.append(group)
            else:
                applied.append(res)

        concat_index = obj.columns if self.axis == 0 else obj.index
        concatenated = concat(applied, join_axes=[concat_index],
                              axis=self.axis, verify_integrity=False)
        return self._set_result_index_ordered(concatenated)

    def transform(self, func, *args, **kwargs):
        """
        Call function producing a like-indexed DataFrame on each group and
        return a DataFrame having the same indexes as the original object
        filled with the transformed values

        Parameters
        ----------
        f : function
            Function to apply to each subframe

        Notes
        -----
        Each subframe is endowed the attribute 'name' in case you need to know
        which group you are working on.

        Examples
        --------
        >>> grouped = df.groupby(lambda x: mapping[x])
        >>> grouped.transform(lambda x: (x - x.mean()) / x.std())
        """

        # optimized transforms
        func = _intercept_cython(func) or func
        if isinstance(func, compat.string_types):
            if func in _cython_transforms:
                # cythonized transform
                return getattr(self, func)(*args, **kwargs)
            else:
                # cythonized aggregation and merge
                result = getattr(self, func)(*args, **kwargs)
        else:
            return self._transform_general(func, *args, **kwargs)

        # a reduction transform
        if not isinstance(result, DataFrame):
            return self._transform_general(func, *args, **kwargs)

        obj = self._obj_with_exclusions
        # nuiscance columns
        if not result.columns.equals(obj.columns):
            return self._transform_general(func, *args, **kwargs)

        results = np.empty_like(obj.values, result.values.dtype)
        indices = self.indices
        for (name, group), (i, row) in zip(self, result.iterrows()):
            indexer = self._get_index(name)
            if len(indexer) > 0:
                results[indexer] = np.tile(row.values,len(indexer)).reshape(len(indexer),-1)

        counts = self.size().fillna(0).values
        if any(counts == 0):
            results = self._try_cast(results, obj[result.columns])

        return (DataFrame(results,columns=result.columns,index=obj.index)
                ._convert(datetime=True))

    def _define_paths(self, func, *args, **kwargs):
        if isinstance(func, compat.string_types):
            fast_path = lambda group: getattr(group, func)(*args, **kwargs)
            slow_path = lambda group: group.apply(
                lambda x: getattr(x, func)(*args, **kwargs), axis=self.axis)
        else:
            fast_path = lambda group: func(group, *args, **kwargs)
            slow_path = lambda group: group.apply(
                lambda x: func(x, *args, **kwargs), axis=self.axis)
        return fast_path, slow_path

    def _choose_path(self, fast_path, slow_path, group):
        path = slow_path
        res = slow_path(group)

        # if we make it here, test if we can use the fast path
        try:
            res_fast = fast_path(group)

            # compare that we get the same results
            if res.shape == res_fast.shape:
                res_r = res.values.ravel()
                res_fast_r = res_fast.values.ravel()
                mask = notnull(res_r)
            if (res_r[mask] == res_fast_r[mask]).all():
                path = fast_path

        except:
            pass
        return path, res

    def _transform_item_by_item(self, obj, wrapper):
        # iterate through columns
        output = {}
        inds = []
        for i, col in enumerate(obj):
            try:
                output[col] = self[col].transform(wrapper)
                inds.append(i)
            except Exception:
                pass

        if len(output) == 0:  # pragma: no cover
            raise TypeError('Transform function invalid for data types')

        columns = obj.columns
        if len(output) < len(obj.columns):
            columns = columns.take(inds)

        return DataFrame(output, index=obj.index, columns=columns)

    def filter(self, func, dropna=True, *args, **kwargs):
        """
        Return a copy of a DataFrame excluding elements from groups that
        do not satisfy the boolean criterion specified by func.

        Parameters
        ----------
        f : function
            Function to apply to each subframe. Should return True or False.
        dropna : Drop groups that do not pass the filter. True by default;
            if False, groups that evaluate False are filled with NaNs.

        Notes
        -----
        Each subframe is endowed the attribute 'name' in case you need to know
        which group you are working on.

        Examples
        --------
        >>> grouped = df.groupby(lambda x: mapping[x])
        >>> grouped.filter(lambda x: x['A'].sum() + x['B'].sum() > 0)
        """

        indices = []

        obj = self._selected_obj
        gen = self.grouper.get_iterator(obj, axis=self.axis)

        for name, group in gen:
            object.__setattr__(group, 'name', name)

            res = func(group, *args, **kwargs)

            try:
                res = res.squeeze()
            except AttributeError:  # allow e.g., scalars and frames to pass
                pass

            # interpret the result of the filter
            if is_bool(res) or (lib.isscalar(res) and isnull(res)):
                if res and notnull(res):
                    indices.append(self._get_index(name))
            else:
                # non scalars aren't allowed
                raise TypeError("filter function returned a %s, "
                                "but expected a scalar bool" %
                                type(res).__name__)

        return self._apply_filter(indices, dropna)


class DataFrameGroupBy(NDFrameGroupBy):
    _apply_whitelist = _dataframe_apply_whitelist
    #
    # Make class defs of attributes on DataFrameGroupBy whitelist.
    for _def_str in _whitelist_method_generator(DataFrame,_apply_whitelist) :
        exec(_def_str)

    _block_agg_axis = 1

    def __getitem__(self, key):
        if self._selection is not None:
            raise Exception('Column(s) %s already selected' % self._selection)

        if isinstance(key, (list, tuple, Series, Index, np.ndarray)):
            if len(self.obj.columns.intersection(key)) != len(key):
                bad_keys = list(set(key).difference(self.obj.columns))
                raise KeyError("Columns not found: %s"
                               % str(bad_keys)[1:-1])
            return DataFrameGroupBy(self.obj, self.grouper, selection=key,
                                    grouper=self.grouper,
                                    exclusions=self.exclusions,
                                    as_index=self.as_index)

        elif not self.as_index:
            if key not in self.obj.columns:
                raise KeyError("Column not found: %s" % key)
            return DataFrameGroupBy(self.obj, self.grouper, selection=key,
                                    grouper=self.grouper,
                                    exclusions=self.exclusions,
                                    as_index=self.as_index)

        else:
            if key not in self.obj:
                raise KeyError("Column not found: %s" % key)
            # kind of a kludge
            return SeriesGroupBy(self.obj[key], selection=key,
                                 grouper=self.grouper,
                                 exclusions=self.exclusions)

    def _wrap_generic_output(self, result, obj):
        result_index = self.grouper.levels[0]

        if result:
            if self.axis == 0:
                result = DataFrame(result, index=obj.columns,
                                   columns=result_index).T
            else:
                result = DataFrame(result, index=obj.index,
                                   columns=result_index)
        else:
            result = DataFrame(result)

        return result

    def _get_data_to_aggregate(self):
        obj = self._obj_with_exclusions
        if self.axis == 1:
            return obj.T._data, 1
        else:
            return obj._data, 1

    def _insert_inaxis_grouper_inplace(self, result):
        # zip in reverse so we can always insert at loc 0
        izip = zip(* map(reversed, (
            self.grouper.names,
            self.grouper.get_group_levels(),
            [grp.in_axis for grp in self.grouper.groupings])))

        for name, lev, in_axis in izip:
            if in_axis:
                result.insert(0, name, lev)

    def _wrap_aggregated_output(self, output, names=None):
        agg_axis = 0 if self.axis == 1 else 1
        agg_labels = self._obj_with_exclusions._get_axis(agg_axis)

        output_keys = self._decide_output_index(output, agg_labels)

        if not self.as_index:
            result = DataFrame(output, columns=output_keys)
            self._insert_inaxis_grouper_inplace(result)
            result = result.consolidate()
        else:
            index = self.grouper.result_index
            result = DataFrame(output, index=index, columns=output_keys)

        if self.axis == 1:
            result = result.T

        return self._reindex_output(result)._convert(datetime=True)

    def _wrap_transformed_output(self, output, names=None):
        return DataFrame(output, index=self.obj.index)

    def _wrap_agged_blocks(self, items, blocks):
        if not self.as_index:
            index = np.arange(blocks[0].values.shape[1])
            mgr = BlockManager(blocks, [items, index])
            result = DataFrame(mgr)

            self._insert_inaxis_grouper_inplace(result)
            result = result.consolidate()
        else:
            index = self.grouper.result_index
            mgr = BlockManager(blocks, [items, index])
            result = DataFrame(mgr)

        if self.axis == 1:
            result = result.T

        return self._reindex_output(result)._convert(datetime=True)

    def _reindex_output(self, result):
        """
        if we have categorical groupers, then we want to make sure that
        we have a fully reindex-output to the levels. These may have not participated in
        the groupings (e.g. may have all been nan groups)

        This can re-expand the output space
        """
        groupings = self.grouper.groupings
        if groupings is None:
            return result
        elif len(groupings) == 1:
            return result
        elif not any([isinstance(ping.grouper, (Categorical, CategoricalIndex))
                      for ping in groupings]):
            return result

        levels_list = [ ping.group_index for ping in groupings ]
        index = MultiIndex.from_product(levels_list, names=self.grouper.names)
        d = { self.obj._get_axis_name(self.axis) : index, 'copy' : False }
        return result.reindex(**d).sortlevel(axis=self.axis)

    def _iterate_column_groupbys(self):
        for i, colname in enumerate(self._selected_obj.columns):
            yield colname, SeriesGroupBy(self._selected_obj.iloc[:, i],
                                         selection=colname,
                                         grouper=self.grouper,
                                         exclusions=self.exclusions)

    def _apply_to_column_groupbys(self, func):
        from pandas.tools.merge import concat
        return concat(
            (func(col_groupby) for _, col_groupby
             in self._iterate_column_groupbys()),
            keys=self._selected_obj.columns, axis=1)

    def count(self):
        """ Compute count of group, excluding missing values """
        from functools import partial
        from pandas.lib import count_level_2d
        from pandas.core.common import _isnull_ndarraylike as isnull

        data, _ = self._get_data_to_aggregate()
        ids, _, ngroups = self.grouper.group_info
        mask = ids != -1

        val = ((mask & ~isnull(blk.get_values())) for blk in data.blocks)
        loc = (blk.mgr_locs for blk in data.blocks)

        counter = partial(count_level_2d, labels=ids, max_bin=ngroups, axis=1)
        blk = map(make_block, map(counter, val), loc)

        return self._wrap_agged_blocks(data.items, list(blk))


from pandas.tools.plotting import boxplot_frame_groupby
DataFrameGroupBy.boxplot = boxplot_frame_groupby


class PanelGroupBy(NDFrameGroupBy):

    def _iterate_slices(self):
        if self.axis == 0:
            # kludge
            if self._selection is None:
                slice_axis = self._selected_obj.items
            else:
                slice_axis = self._selection_list
            slicer = lambda x: self._selected_obj[x]
        else:
            raise NotImplementedError("axis other than 0 is not supported")

        for val in slice_axis:
            if val in self.exclusions:
                continue

            yield val, slicer(val)

    def aggregate(self, arg, *args, **kwargs):
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
        if isinstance(arg, compat.string_types):
            return getattr(self, arg)(*args, **kwargs)

        return self._aggregate_generic(arg, *args, **kwargs)

    def _wrap_generic_output(self, result, obj):
        if self.axis == 0:
            new_axes = list(obj.axes)
            new_axes[0] = self.grouper.result_index
        elif self.axis == 1:
            x, y, z = obj.axes
            new_axes = [self.grouper.result_index, z, x]
        else:
            x, y, z = obj.axes
            new_axes = [self.grouper.result_index, y, x]

        result = Panel._from_axes(result, new_axes)

        if self.axis == 1:
            result = result.swapaxes(0, 1).swapaxes(0, 2)
        elif self.axis == 2:
            result = result.swapaxes(0, 2)

        return result

    def _aggregate_item_by_item(self, func, *args, **kwargs):
        obj = self._obj_with_exclusions
        result = {}

        if self.axis > 0:
            for item in obj:
                try:
                    itemg = DataFrameGroupBy(obj[item],
                                             axis=self.axis - 1,
                                             grouper=self.grouper)
                    result[item] = itemg.aggregate(func, *args, **kwargs)
                except (ValueError, TypeError):
                    raise
            new_axes = list(obj.axes)
            new_axes[self.axis] = self.grouper.result_index
            return Panel._from_axes(result, new_axes)
        else:
            raise ValueError("axis value must be greater than 0")

    def _wrap_aggregated_output(self, output, names=None):
        raise AbstractMethodError(self)


class NDArrayGroupBy(GroupBy):
    pass


#----------------------------------------------------------------------
# Splitting / application


class DataSplitter(object):

    def __init__(self, data, labels, ngroups, axis=0):
        self.data = data
        self.labels = com._ensure_int64(labels)
        self.ngroups = ngroups

        self.axis = axis

    @cache_readonly
    def slabels(self):
        # Sorted labels
        return com.take_nd(self.labels, self.sort_idx, allow_fill=False)

    @cache_readonly
    def sort_idx(self):
        # Counting sort indexer
        return _get_group_index_sorter(self.labels, self.ngroups)

    def __iter__(self):
        sdata = self._get_sorted_data()

        if self.ngroups == 0:
            raise StopIteration

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
        return self.data.take(self.sort_idx, axis=self.axis, convert=False)

    def _chop(self, sdata, slice_obj):
        return sdata.iloc[slice_obj]

    def apply(self, f):
        raise AbstractMethodError(self)


class ArraySplitter(DataSplitter):
    pass


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
        except:
            # fails when all -1
            return [], True

        sdata = self._get_sorted_data()
        results, mutated = lib.apply_frame_axis0(sdata, f, names, starts, ends)

        return results, mutated

    def _chop(self, sdata, slice_obj):
        if self.axis == 0:
            return sdata.iloc[slice_obj]
        else:
            return sdata._slice(slice_obj, axis=1)  # ix[:, slice_obj]


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


#----------------------------------------------------------------------
# Misc utilities


def get_group_index(labels, shape, sort, xnull):
    """
    For the particular label_list, gets the offsets into the hypothetical list
    representing the totally ordered cartesian product of all possible label
    combinations, *as long as* this space fits within int64 bounds;
    otherwise, though group indices identify unique combinations of
    labels, they cannot be deconstructed.
    - If `sort`, rank of returned ids preserve lexical ranks of labels.
      i.e. returned id's can be used to do lexical sort on labels;
    - If `xnull` nulls (-1 labels) are passed through.

    Parameters
    ----------
    labels: sequence of arrays
        Integers identifying levels at each location
    shape: sequence of ints same length as labels
        Number of unique levels at each location
    sort: boolean
        If the ranks of returned ids should match lexical ranks of labels
    xnull: boolean
        If true nulls are excluded. i.e. -1 values in the labels are
        passed through
    Returns
    -------
    An array of type int64 where two elements are equal if their corresponding
    labels are equal at all location.
    """
    def _int64_cut_off(shape):
        acc = long(1)
        for i, mul in enumerate(shape):
            acc *= long(mul)
            if not acc < _INT64_MAX:
                return i
        return len(shape)

    def loop(labels, shape):
        # how many levels can be done without overflow:
        nlev = _int64_cut_off(shape)

        # compute flat ids for the first `nlev` levels
        stride = np.prod(shape[1:nlev], dtype='i8')
        out = stride * labels[0].astype('i8', subok=False, copy=False)

        for i in range(1, nlev):
            stride //= shape[i]
            out += labels[i] * stride

        if xnull: # exclude nulls
            mask = labels[0] == -1
            for lab in labels[1:nlev]:
                mask |= lab == -1
            out[mask] = -1

        if nlev == len(shape):  # all levels done!
            return out

        # compress what has been done so far in order to avoid overflow
        # to retain lexical ranks, obs_ids should be sorted
        comp_ids, obs_ids = _compress_group_index(out, sort=sort)

        labels = [comp_ids] + labels[nlev:]
        shape = [len(obs_ids)] + shape[nlev:]

        return loop(labels, shape)

    def maybe_lift(lab, size):  # pormote nan values
        return (lab + 1, size + 1) if (lab == -1).any() else (lab, size)

    labels = map(com._ensure_int64, labels)
    if not xnull:
        labels, shape = map(list, zip(*map(maybe_lift, labels, shape)))

    return loop(list(labels), list(shape))


_INT64_MAX = np.iinfo(np.int64).max


def _int64_overflow_possible(shape):
    the_prod = long(1)
    for x in shape:
        the_prod *= long(x)

    return the_prod >= _INT64_MAX


def decons_group_index(comp_labels, shape):
    # reconstruct labels
    if _int64_overflow_possible(shape):
        # at some point group indices are factorized,
        # and may not be deconstructed here! wrong path!
        raise ValueError('cannot deconstruct factorized group indices!')

    label_list = []
    factor = 1
    y = 0
    x = comp_labels
    for i in reversed(range(len(shape))):
        labels = (x - y) % (factor * shape[i]) // factor
        np.putmask(labels, comp_labels < 0, -1)
        label_list.append(labels)
        y = labels * factor
        factor *= shape[i]
    return label_list[::-1]


def decons_obs_group_ids(comp_ids, obs_ids, shape, labels, xnull):
    """
    reconstruct labels from observed group ids

    Parameters
    ----------
    xnull: boolean,
        if nulls are excluded; i.e. -1 labels are passed through
    """
    from pandas.hashtable import unique_label_indices

    if not xnull:
        lift = np.fromiter(((a == -1).any() for a in labels), dtype='i8')
        shape = np.asarray(shape, dtype='i8') + lift

    if not _int64_overflow_possible(shape):
        # obs ids are deconstructable! take the fast route!
        out = decons_group_index(obs_ids, shape)
        return out if xnull or not lift.any() \
                else [x - y for x, y in zip(out, lift)]

    i = unique_label_indices(comp_ids)
    i8copy = lambda a: a.astype('i8', subok=False, copy=True)
    return [i8copy(lab[i]) for lab in labels]


def _indexer_from_factorized(labels, shape, compress=True):
    ids = get_group_index(labels, shape, sort=True, xnull=False)

    if not compress:
        ngroups = (ids.size and ids.max()) + 1
    else:
        ids, obs = _compress_group_index(ids, sort=True)
        ngroups = len(obs)

    return _get_group_index_sorter(ids, ngroups)


def _lexsort_indexer(keys, orders=None, na_position='last'):
    labels = []
    shape = []
    if isinstance(orders, bool):
        orders = [orders] * len(keys)
    elif orders is None:
        orders = [True] * len(keys)

    for key, order in zip(keys, orders):

        # we are already a Categorical
        if is_categorical_dtype(key):
            c = key

        # create the Categorical
        else:
            c = Categorical(key,ordered=True)

        if na_position not in ['last','first']:
            raise ValueError('invalid na_position: {!r}'.format(na_position))

        n = len(c.categories)
        codes = c.codes.copy()

        mask = (c.codes == -1)
        if order: # ascending
            if na_position == 'last':
                codes = np.where(mask, n, codes)
            elif na_position == 'first':
                codes += 1
        else: # not order means descending
            if na_position == 'last':
                codes = np.where(mask, n, n-codes-1)
            elif na_position == 'first':
                codes = np.where(mask, 0, n-codes)
        if mask.any():
            n += 1

        shape.append(n)
        labels.append(codes)

    return _indexer_from_factorized(labels, shape)

def _nargsort(items, kind='quicksort', ascending=True, na_position='last'):
    """
    This is intended to be a drop-in replacement for np.argsort which handles NaNs
    It adds ascending and na_position parameters.
    GH #6399, #5231
    """

    # specially handle Categorical
    if is_categorical_dtype(items):
        return items.argsort(ascending=ascending)

    items = np.asanyarray(items)
    idx = np.arange(len(items))
    mask = isnull(items)
    non_nans = items[~mask]
    non_nan_idx = idx[~mask]
    nan_idx = np.nonzero(mask)[0]
    if not ascending:
        non_nans = non_nans[::-1]
        non_nan_idx = non_nan_idx[::-1]
    indexer = non_nan_idx[non_nans.argsort(kind=kind)]
    if not ascending:
        indexer = indexer[::-1]
    # Finally, place the NaNs at the end or the beginning according to na_position
    if na_position == 'last':
        indexer = np.concatenate([indexer, nan_idx])
    elif na_position == 'first':
        indexer = np.concatenate([nan_idx, indexer])
    else:
        raise ValueError('invalid na_position: {!r}'.format(na_position))
    return indexer


class _KeyMapper(object):

    """
    Ease my suffering. Map compressed group id -> key tuple
    """

    def __init__(self, comp_ids, ngroups, labels, levels):
        self.levels = levels
        self.labels = labels
        self.comp_ids = comp_ids.astype(np.int64)

        self.k = len(labels)
        self.tables = [_hash.Int64HashTable(ngroups) for _ in range(self.k)]

        self._populate_tables()

    def _populate_tables(self):
        for labs, table in zip(self.labels, self.tables):
            table.map(self.comp_ids, labs.astype(np.int64))

    def get_key(self, comp_id):
        return tuple(level[table.get_item(comp_id)]
                     for table, level in zip(self.tables, self.levels))


def _get_indices_dict(label_list, keys):
    shape = list(map(len, keys))

    group_index = get_group_index(label_list, shape, sort=True, xnull=True)
    ngroups = ((group_index.size and group_index.max()) + 1) \
              if _int64_overflow_possible(shape) \
              else np.prod(shape, dtype='i8')

    sorter = _get_group_index_sorter(group_index, ngroups)

    sorted_labels = [lab.take(sorter) for lab in label_list]
    group_index = group_index.take(sorter)

    return lib.indices_fast(sorter, group_index, keys, sorted_labels)


#----------------------------------------------------------------------
# sorting levels...cleverly?

def _get_group_index_sorter(group_index, ngroups):
    """
    _algos.groupsort_indexer implements `counting sort` and it is at least
    O(ngroups), where
        ngroups = prod(shape)
        shape = map(len, keys)
    that is, linear in the number of combinations (cartesian product) of unique
    values of groupby keys. This can be huge when doing multi-key groupby.
    np.argsort(kind='mergesort') is O(count x log(count)) where count is the
    length of the data-frame;
    Both algorithms are `stable` sort and that is necessary for correctness of
    groupby operations. e.g. consider:
        df.groupby(key)[col].transform('first')
    """
    count = len(group_index)
    alpha = 0.0  # taking complexities literally; there may be
    beta  = 1.0  # some room for fine-tuning these parameters
    if alpha + beta * ngroups < count * np.log(count):
        sorter, _ = _algos.groupsort_indexer(com._ensure_int64(group_index),
                                             ngroups)
        return com._ensure_platform_int(sorter)
    else:
        return group_index.argsort(kind='mergesort')


def _compress_group_index(group_index, sort=True):
    """
    Group_index is offsets into cartesian product of all possible labels. This
    space can be huge, so this function compresses it, by computing offsets
    (comp_ids) into the list of unique labels (obs_group_ids).
    """

    size_hint = min(len(group_index), _hash._SIZE_HINT_LIMIT)
    table = _hash.Int64HashTable(size_hint)

    group_index = com._ensure_int64(group_index)

    # note, group labels come out ascending (ie, 1,2,3 etc)
    comp_ids, obs_group_ids = table.get_labels_groupby(group_index)

    if sort and len(obs_group_ids) > 0:
        obs_group_ids, comp_ids = _reorder_by_uniques(obs_group_ids, comp_ids)

    return comp_ids, obs_group_ids


def _reorder_by_uniques(uniques, labels):
    # sorter is index where elements ought to go
    sorter = uniques.argsort()

    # reverse_indexer is where elements came from
    reverse_indexer = np.empty(len(sorter), dtype=np.int64)
    reverse_indexer.put(sorter, np.arange(len(sorter)))

    mask = labels < 0

    # move labels to right locations (ie, unsort ascending labels)
    labels = com.take_nd(reverse_indexer, labels, allow_fill=False)
    np.putmask(labels, mask, -1)

    # sort observed ids
    uniques = com.take_nd(uniques, sorter, allow_fill=False)

    return uniques, labels


_func_table = {
    builtins.sum: np.sum,
    builtins.max: np.max,
    builtins.min: np.min
}


_cython_table = {
    builtins.sum: 'sum',
    builtins.max: 'max',
    builtins.min: 'min',
    np.sum: 'sum',
    np.mean: 'mean',
    np.prod: 'prod',
    np.std: 'std',
    np.var: 'var',
    np.median: 'median',
    np.max: 'max',
    np.min: 'min',
    np.cumprod: 'cumprod',
    np.cumsum: 'cumsum'
}


def _intercept_function(func):
    return _func_table.get(func, func)


def _intercept_cython(func):
    return _cython_table.get(func)


def _groupby_indices(values):
    return _algos.groupby_indices(_values_from_object(com._ensure_object(values)))


def numpy_groupby(data, labels, axis=0):
    s = np.argsort(labels)
    keys, inv = np.unique(labels, return_inverse=True)
    i = inv.take(s)
    groups_at = np.where(i != np.concatenate(([-1], i[:-1])))[0]
    ordered_data = data.take(s, axis=axis)
    group_sums = np.add.reduceat(ordered_data, groups_at, axis=axis)

    return group_sums
