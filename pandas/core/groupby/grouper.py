"""
Provide user facing operators for doing the split part of the
split-apply-combine paradigm.
"""

import warnings
import numpy as np

from pandas.util._decorators import cache_readonly

from pandas import compat
from pandas.compat import zip, callable

from pandas.core.dtypes.generic import ABCSeries
from pandas.core.arrays import ExtensionArray, Categorical
from pandas.core.index import (
    Index, MultiIndex, CategoricalIndex)
from pandas.core.dtypes.common import (
    ensure_categorical,
    is_hashable,
    is_list_like,
    is_timedelta64_dtype,
    is_datetime64_dtype,
    is_categorical_dtype,
    is_scalar)
from pandas.core.series import Series
from pandas.core.frame import DataFrame
import pandas.core.common as com
from pandas.core.groupby.ops import BaseGrouper
import pandas.core.algorithms as algorithms
from pandas.io.formats.printing import pprint_thing


class Grouper(object):
    """
    A Grouper allows the user to specify a groupby instruction for a target
    object

    This specification will select a column via the key parameter, or if the
    level and/or axis parameters are given, a level of the index of the target
    object.

    These are local specifications and will override 'global' settings,
    that is the parameters axis and level which are passed to the groupby
    itself.

    Parameters
    ----------
    key : string, defaults to None
        groupby key, which selects the grouping column of the target
    level : name/number, defaults to None
        the level for the target index
    freq : string / frequency object, defaults to None
        This will groupby the specified frequency if the target selection
        (via key or level) is a datetime-like object. For full specification
        of available frequencies, please see `here
        <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`_.
    axis : number/name of the axis, defaults to 0
    sort : boolean, default to False
        whether to sort the resulting labels

    additional kwargs to control time-like groupers (when `freq` is passed)

    closed : closed end of interval; 'left' or 'right'
    label : interval boundary to use for labeling; 'left' or 'right'
    convention : {'start', 'end', 'e', 's'}
        If grouper is PeriodIndex
    base, loffset

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
    _attributes = ('key', 'level', 'freq', 'axis', 'sort')

    def __new__(cls, *args, **kwargs):
        if kwargs.get('freq') is not None:
            from pandas.core.resample import TimeGrouper
            cls = TimeGrouper
        return super(Grouper, cls).__new__(cls)

    def __init__(self, key=None, level=None, freq=None, axis=0, sort=False):
        self.key = key
        self.level = level
        self.freq = freq
        self.axis = axis
        self.sort = sort

        self.grouper = None
        self.obj = None
        self.indexer = None
        self.binner = None
        self._grouper = None

    @property
    def ax(self):
        return self.grouper

    def _get_grouper(self, obj, validate=True):
        """
        Parameters
        ----------
        obj : the subject object
        validate : boolean, default True
            if True, validate the grouper

        Returns
        -------
        a tuple of binner, grouper, obj (possibly sorted)
        """

        self._set_grouper(obj)
        self.grouper, exclusions, self.obj = _get_grouper(self.obj, [self.key],
                                                          axis=self.axis,
                                                          level=self.level,
                                                          sort=self.sort,
                                                          validate=validate)
        return self.binner, self.grouper, self.obj

    def _set_grouper(self, obj, sort=False):
        """
        given an object and the specifications, setup the internal grouper
        for this particular specification

        Parameters
        ----------
        obj : the subject object
        sort : bool, default False
            whether the resulting grouper should be sorted
        """

        if self.key is not None and self.level is not None:
            raise ValueError(
                "The Grouper cannot specify both a key and a level!")

        # Keep self.grouper value before overriding
        if self._grouper is None:
            self._grouper = self.grouper

        # the key must be a valid info item
        if self.key is not None:
            key = self.key
            # The 'on' is already defined
            if getattr(self.grouper, 'name', None) == key and \
                    isinstance(obj, ABCSeries):
                ax = self._grouper.take(obj.index)
            else:
                if key not in obj._info_axis:
                    raise KeyError(
                        "The grouper name {0} is not found".format(key))
                ax = Index(obj[key], name=key)

        else:
            ax = obj._get_axis(self.axis)
            if self.level is not None:
                level = self.level

                # if a level is given it must be a mi level or
                # equivalent to the axis name
                if isinstance(ax, MultiIndex):
                    level = ax._get_level_number(level)
                    ax = Index(ax._get_level_values(level),
                               name=ax.names[level])

                else:
                    if level not in (0, ax.name):
                        raise ValueError(
                            "The level {0} is not valid".format(level))

        # possibly sort
        if (self.sort or sort) and not ax.is_monotonic:
            # use stable sort to support first, last, nth
            indexer = self.indexer = ax.argsort(kind='mergesort')
            ax = ax.take(indexer)
            obj = obj._take(indexer, axis=self.axis, is_copy=False)

        self.obj = obj
        self.grouper = ax
        return self.grouper

    @property
    def groups(self):
        return self.grouper.groups

    def __repr__(self):
        attrs_list = ["{}={!r}".format(attr_name, getattr(self, attr_name))
                      for attr_name in self._attributes
                      if getattr(self, attr_name) is not None]
        attrs = ", ".join(attrs_list)
        cls_name = self.__class__.__name__
        return "{}({})".format(cls_name, attrs)


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
    observed : boolean, default False
        If we are a Categorical, use the observed values
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
                 sort=True, observed=False, in_axis=False):

        self.name = name
        self.level = level
        self.grouper = _convert_grouper(index, grouper)
        self.all_grouper = None
        self.index = index
        self.sort = sort
        self.obj = obj
        self.observed = observed
        self.in_axis = in_axis

        # right place for this?
        if isinstance(grouper, (Series, Index)) and name is None:
            self.name = grouper.name

        if isinstance(grouper, MultiIndex):
            self.grouper = grouper.values

        # we have a single grouper which may be a myriad of things,
        # some of which are dependent on the passing in level

        if level is not None:
            if not isinstance(level, int):
                if level not in index.names:
                    raise AssertionError('Level %s not in index' % str(level))
                level = index.names.index(level)

            if self.name is None:
                self.name = index.names[level]

            self.grouper, self._labels, self._group_index = \
                index._get_grouper_for_level(self.grouper, level)

        # a passed Grouper like, directly get the grouper in the same way
        # as single grouper groupby, use the group_info to get labels
        elif isinstance(self.grouper, Grouper):
            # get the new grouper; we already have disambiguated
            # what key/level refer to exactly, don't need to
            # check again as we have by this point converted these
            # to an actual value (rather than a pd.Grouper)
            _, grouper, _ = self.grouper._get_grouper(self.obj, validate=False)
            if self.name is None:
                self.name = grouper.result_index.name
            self.obj = self.grouper.obj
            self.grouper = grouper

        else:
            if self.grouper is None and self.name is not None:
                self.grouper = self.obj[self.name]

            elif isinstance(self.grouper, (list, tuple)):
                self.grouper = com.asarray_tuplesafe(self.grouper)

            # a passed Categorical
            elif is_categorical_dtype(self.grouper):

                from pandas.core.groupby.categorical import recode_for_groupby
                self.grouper, self.all_grouper = recode_for_groupby(
                    self.grouper, self.sort, observed)
                categories = self.grouper.categories

                # we make a CategoricalIndex out of the cat grouper
                # preserving the categories / ordered attributes
                self._labels = self.grouper.codes
                if observed:
                    codes = algorithms.unique1d(self.grouper.codes)
                else:
                    codes = np.arange(len(categories))

                self._group_index = CategoricalIndex(
                    Categorical.from_codes(
                        codes=codes,
                        categories=categories,
                        ordered=self.grouper.ordered))

            # we are done
            if isinstance(self.grouper, Grouping):
                self.grouper = self.grouper.grouper

            # no level passed
            elif not isinstance(self.grouper,
                                (Series, Index, ExtensionArray, np.ndarray)):
                if getattr(self.grouper, 'ndim', 1) != 1:
                    t = self.name or str(type(self.grouper))
                    raise ValueError("Grouper for '%s' not 1-dimensional" % t)
                self.grouper = self.index.map(self.grouper)
                if not (hasattr(self.grouper, "__len__") and
                        len(self.grouper) == len(self.index)):
                    errmsg = ('Grouper result violates len(labels) == '
                              'len(data)\nresult: %s' %
                              pprint_thing(self.grouper))
                    self.grouper = None  # Try for sanity
                    raise AssertionError(errmsg)

        # if we have a date/time-like grouper, make sure that we have
        # Timestamps like
        if getattr(self.grouper, 'dtype', None) is not None:
            if is_datetime64_dtype(self.grouper):
                from pandas import to_datetime
                self.grouper = to_datetime(self.grouper)
            elif is_timedelta64_dtype(self.grouper):
                from pandas import to_timedelta
                self.grouper = to_timedelta(self.grouper)

    def __repr__(self):
        return 'Grouping({0})'.format(self.name)

    def __iter__(self):
        return iter(self.indices)

    _labels = None
    _group_index = None

    @property
    def ngroups(self):
        return len(self.group_index)

    @cache_readonly
    def indices(self):
        # we have a list of groupers
        if isinstance(self.grouper, BaseGrouper):
            return self.grouper.indices

        values = ensure_categorical(self.grouper)
        return values._reverse_indexer()

    @property
    def labels(self):
        if self._labels is None:
            self._make_labels()
        return self._labels

    @cache_readonly
    def result_index(self):
        if self.all_grouper is not None:
            from pandas.core.groupby.categorical import recode_from_groupby
            return recode_from_groupby(self.all_grouper,
                                       self.sort, self.group_index)
        return self.group_index

    @property
    def group_index(self):
        if self._group_index is None:
            self._make_labels()
        return self._group_index

    def _make_labels(self):
        if self._labels is None or self._group_index is None:
            # we have a list of groupers
            if isinstance(self.grouper, BaseGrouper):
                labels = self.grouper.label_info
                uniques = self.grouper.result_index
            else:
                labels, uniques = algorithms.factorize(
                    self.grouper, sort=self.sort)
                uniques = Index(uniques, name=self.name)
            self._labels = labels
            self._group_index = uniques

    @cache_readonly
    def groups(self):
        return self.index.groupby(Categorical.from_codes(self.labels,
                                                         self.group_index))


def _get_grouper(obj, key=None, axis=0, level=None, sort=True,
                 observed=False, mutated=False, validate=True):
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

    If observed & we have a categorical grouper, only show the observed
    values

    If validate, then check for key/level overlaps

    """
    group_axis = obj._get_axis(axis)

    # validate that the passed single level is compatible with the passed
    # axis of the object
    if level is not None:
        # TODO: These if-block and else-block are almost same.
        # MultiIndex instance check is removable, but it seems that there are
        # some processes only for non-MultiIndex in else-block,
        # eg. `obj.index.name != level`. We have to consider carefully whether
        # these are applicable for MultiIndex. Even if these are applicable,
        # we need to check if it makes no side effect to subsequent processes
        # on the outside of this condition.
        # (GH 17621)
        if isinstance(group_axis, MultiIndex):
            if is_list_like(level) and len(level) == 1:
                level = level[0]

            if key is None and is_scalar(level):
                # Get the level values from group_axis
                key = group_axis.get_level_values(level)
                level = None

        else:
            # allow level to be a length-one list-like object
            # (e.g., level=[0])
            # GH 13901
            if is_list_like(level):
                nlevels = len(level)
                if nlevels == 1:
                    level = level[0]
                elif nlevels == 0:
                    raise ValueError('No group keys passed!')
                else:
                    raise ValueError('multiple levels only valid with '
                                     'MultiIndex')

            if isinstance(level, compat.string_types):
                if obj.index.name != level:
                    raise ValueError('level name %s is not the name of the '
                                     'index' % level)
            elif level > 0 or level < -1:
                raise ValueError('level > 0 or level < -1 only valid with '
                                 ' MultiIndex')

            # NOTE: `group_axis` and `group_axis.get_level_values(level)`
            # are same in this section.
            level = None
            key = group_axis

    # a passed-in Grouper, directly convert
    if isinstance(key, Grouper):
        binner, grouper, obj = key._get_grouper(obj, validate=False)
        if key.key is None:
            return grouper, [], obj
        else:
            return grouper, set([key.key]), obj

    # already have a BaseGrouper, just return it
    elif isinstance(key, BaseGrouper):
        return key, [], obj

    # In the future, a tuple key will always mean an actual key,
    # not an iterable of keys. In the meantime, we attempt to provide
    # a warning. We can assume that the user wanted a list of keys when
    # the key is not in the index. We just have to be careful with
    # unhashble elements of `key`. Any unhashable elements implies that
    # they wanted a list of keys.
    # https://github.com/pandas-dev/pandas/issues/18314
    is_tuple = isinstance(key, tuple)
    all_hashable = is_tuple and is_hashable(key)

    if is_tuple:
        if ((all_hashable and key not in obj and set(key).issubset(obj))
                or not all_hashable):
            # column names ('a', 'b') -> ['a', 'b']
            # arrays like (a, b) -> [a, b]
            msg = ("Interpreting tuple 'by' as a list of keys, rather than "
                   "a single key. Use 'by=[...]' instead of 'by=(...)'. In "
                   "the future, a tuple will always mean a single key.")
            warnings.warn(msg, FutureWarning, stacklevel=5)
            key = list(key)

    if not isinstance(key, list):
        keys = [key]
        match_axis_length = False
    else:
        keys = key
        match_axis_length = len(keys) == len(group_axis)

    # what are we after, exactly?
    any_callable = any(callable(g) or isinstance(g, dict) for g in keys)
    any_groupers = any(isinstance(g, Grouper) for g in keys)
    any_arraylike = any(isinstance(g, (list, tuple, Series, Index, np.ndarray))
                        for g in keys)

    try:
        if isinstance(obj, DataFrame):
            all_in_columns_index = all(g in obj.columns or g in obj.index.names
                                       for g in keys)
        else:
            all_in_columns_index = False
    except Exception:
        all_in_columns_index = False

    if not any_callable and not all_in_columns_index and \
       not any_arraylike and not any_groupers and \
       match_axis_length and level is None:
        keys = [com.asarray_tuplesafe(keys)]

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

    # if the grouper is obj[name]
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
            if gpr in obj:
                if validate:
                    stacklevel = 5  # Number of stack levels from df.groupby
                    obj._check_label_or_level_ambiguity(
                        gpr, stacklevel=stacklevel)
                in_axis, name, gpr = True, gpr, obj[gpr]
                exclusions.append(name)
            elif obj._is_level_reference(gpr):
                in_axis, name, level, gpr = False, None, gpr, None
            else:
                raise KeyError(gpr)
        elif isinstance(gpr, Grouper) and gpr.key is not None:
            # Add key to exclusions
            exclusions.append(gpr.key)
            in_axis, name = False, None
        else:
            in_axis, name = False, None

        if is_categorical_dtype(gpr) and len(gpr) != obj.shape[axis]:
            raise ValueError(
                ("Length of grouper ({len_gpr}) and axis ({len_axis})"
                 " must be same length"
                 .format(len_gpr=len(gpr), len_axis=obj.shape[axis])))

        # create the Grouping
        # allow us to passing the actual Grouping as the gpr
        ping = Grouping(group_axis,
                        gpr,
                        obj=obj,
                        name=name,
                        level=level,
                        sort=sort,
                        observed=observed,
                        in_axis=in_axis) \
            if not isinstance(gpr, Grouping) else gpr

        groupings.append(ping)

    if len(groupings) == 0:
        raise ValueError('No group keys passed!')

    # create the internals grouper
    grouper = BaseGrouper(group_axis, groupings, sort=sort, mutated=mutated)
    return grouper, exclusions, obj


def _is_label_like(val):
    return (isinstance(val, (compat.string_types, tuple)) or
            (val is not None and is_scalar(val)))


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
            raise ValueError('Grouper and axis must be same length')
        return grouper
    else:
        return grouper
