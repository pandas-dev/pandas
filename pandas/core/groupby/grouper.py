"""
Provide user facing operators for doing the split part of the
split-apply-combine paradigm.
"""

from typing import Dict, Hashable, List, Optional, Tuple

import numpy as np

from pandas._typing import FrameOrSeries
from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.common import (
    ensure_categorical,
    is_categorical_dtype,
    is_datetime64_dtype,
    is_list_like,
    is_scalar,
    is_timedelta64_dtype,
)
from pandas.core.dtypes.generic import ABCSeries

import pandas.core.algorithms as algorithms
from pandas.core.arrays import Categorical, ExtensionArray
import pandas.core.common as com
from pandas.core.frame import DataFrame
from pandas.core.groupby import ops
from pandas.core.groupby.categorical import recode_for_groupby, recode_from_groupby
from pandas.core.indexes.api import CategoricalIndex, Index, MultiIndex
from pandas.core.series import Series

from pandas.io.formats.printing import pprint_thing


class Grouper:
    """
    A Grouper allows the user to specify a groupby instruction for an object.

    This specification will select a column via the key parameter, or if the
    level and/or axis parameters are given, a level of the index of the target
    object.

    If `axis` and/or `level` are passed as keywords to both `Grouper` and
    `groupby`, the values passed to `Grouper` take precedence.

    Parameters
    ----------
    key : str, defaults to None
        Groupby key, which selects the grouping column of the target.
    level : name/number, defaults to None
        The level for the target index.
    freq : str / frequency object, defaults to None
        This will groupby the specified frequency if the target selection
        (via key or level) is a datetime-like object. For full specification
        of available frequencies, please see `here
        <http://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases>`_.
    axis : str, int, defaults to 0
        Number/name of the axis.
    sort : bool, default to False
        Whether to sort the resulting labels.
    closed : {'left' or 'right'}
        Closed end of interval. Only when `freq` parameter is passed.
    label : {'left' or 'right'}
        Interval boundary to use for labeling.
        Only when `freq` parameter is passed.
    convention : {'start', 'end', 'e', 's'}
        If grouper is PeriodIndex and `freq` parameter is passed.
    base : int, default 0
        Only when `freq` parameter is passed.
    loffset : str, DateOffset, timedelta object
        Only when `freq` parameter is passed.

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

    _attributes: Tuple[str, ...] = ("key", "level", "freq", "axis", "sort")

    def __new__(cls, *args, **kwargs):
        if kwargs.get("freq") is not None:
            from pandas.core.resample import TimeGrouper

            cls = TimeGrouper
        return super().__new__(cls)

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

    def _get_grouper(self, obj, validate: bool = True):
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
        self.grouper, _, self.obj = get_grouper(
            self.obj,
            [self.key],
            axis=self.axis,
            level=self.level,
            sort=self.sort,
            validate=validate,
        )
        return self.binner, self.grouper, self.obj

    def _set_grouper(self, obj: FrameOrSeries, sort: bool = False):
        """
        given an object and the specifications, setup the internal grouper
        for this particular specification

        Parameters
        ----------
        obj : Series or DataFrame
        sort : bool, default False
            whether the resulting grouper should be sorted
        """
        assert obj is not None

        if self.key is not None and self.level is not None:
            raise ValueError("The Grouper cannot specify both a key and a level!")

        # Keep self.grouper value before overriding
        if self._grouper is None:
            self._grouper = self.grouper

        # the key must be a valid info item
        if self.key is not None:
            key = self.key
            # The 'on' is already defined
            if getattr(self.grouper, "name", None) == key and isinstance(
                obj, ABCSeries
            ):
                ax = self._grouper.take(obj.index)
            else:
                if key not in obj._info_axis:
                    raise KeyError(f"The grouper name {key} is not found")
                ax = Index(obj[key], name=key)

        else:
            ax = obj._get_axis(self.axis)
            if self.level is not None:
                level = self.level

                # if a level is given it must be a mi level or
                # equivalent to the axis name
                if isinstance(ax, MultiIndex):
                    level = ax._get_level_number(level)
                    ax = Index(ax._get_level_values(level), name=ax.names[level])

                else:
                    if level not in (0, ax.name):
                        raise ValueError(f"The level {level} is not valid")

        # possibly sort
        if (self.sort or sort) and not ax.is_monotonic:
            # use stable sort to support first, last, nth
            indexer = self.indexer = ax.argsort(kind="mergesort")
            ax = ax.take(indexer)
            obj = obj.take(indexer, axis=self.axis, is_copy=False)

        self.obj = obj
        self.grouper = ax
        return self.grouper

    @property
    def groups(self):
        return self.grouper.groups

    def __repr__(self) -> str:
        attrs_list = (
            f"{attr_name}={repr(getattr(self, attr_name))}"
            for attr_name in self._attributes
            if getattr(self, attr_name) is not None
        )
        attrs = ", ".join(attrs_list)
        cls_name = type(self).__name__
        return f"{cls_name}({attrs})"


class Grouping:
    """
    Holds the grouping information for a single key

    Parameters
    ----------
    index : Index
    grouper :
    obj Union[DataFrame, Series]:
    name :
    level :
    observed : bool, default False
        If we are a Categorical, use the observed values
    in_axis : if the Grouping is a column in self.obj and hence among
        Groupby.exclusions list

    Returns
    -------
    **Attributes**:
      * indices : dict of {group -> index_list}
      * codes : ndarray, group codes
      * group_index : unique groups
      * groups : dict of {group -> label_list}
    """

    def __init__(
        self,
        index: Index,
        grouper=None,
        obj: Optional[FrameOrSeries] = None,
        name=None,
        level=None,
        sort: bool = True,
        observed: bool = False,
        in_axis: bool = False,
    ):
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
                    raise AssertionError(f"Level {level} not in index")
                level = index.names.index(level)

            if self.name is None:
                self.name = index.names[level]

            (
                self.grouper,
                self._codes,
                self._group_index,
            ) = index._get_grouper_for_level(self.grouper, level)

        # a passed Grouper like, directly get the grouper in the same way
        # as single grouper groupby, use the group_info to get codes
        elif isinstance(self.grouper, Grouper):
            # get the new grouper; we already have disambiguated
            # what key/level refer to exactly, don't need to
            # check again as we have by this point converted these
            # to an actual value (rather than a pd.Grouper)
            _, grouper, _ = self.grouper._get_grouper(self.obj, validate=False)
            if self.name is None:
                self.name = grouper.result_index.name
            self.obj = self.grouper.obj
            self.grouper = grouper._get_grouper()

        else:
            if self.grouper is None and self.name is not None and self.obj is not None:
                self.grouper = self.obj[self.name]

            elif isinstance(self.grouper, (list, tuple)):
                self.grouper = com.asarray_tuplesafe(self.grouper)

            # a passed Categorical
            elif is_categorical_dtype(self.grouper):

                self.grouper, self.all_grouper = recode_for_groupby(
                    self.grouper, self.sort, observed
                )
                categories = self.grouper.categories

                # we make a CategoricalIndex out of the cat grouper
                # preserving the categories / ordered attributes
                self._codes = self.grouper.codes
                if observed:
                    codes = algorithms.unique1d(self.grouper.codes)
                    codes = codes[codes != -1]
                    if sort or self.grouper.ordered:
                        codes = np.sort(codes)
                else:
                    codes = np.arange(len(categories))

                self._group_index = CategoricalIndex(
                    Categorical.from_codes(
                        codes=codes, categories=categories, ordered=self.grouper.ordered
                    ),
                    name=self.name,
                )

            # we are done
            if isinstance(self.grouper, Grouping):
                self.grouper = self.grouper.grouper

            # no level passed
            elif not isinstance(
                self.grouper, (Series, Index, ExtensionArray, np.ndarray)
            ):
                if getattr(self.grouper, "ndim", 1) != 1:
                    t = self.name or str(type(self.grouper))
                    raise ValueError(f"Grouper for '{t}' not 1-dimensional")
                self.grouper = self.index.map(self.grouper)
                if not (
                    hasattr(self.grouper, "__len__")
                    and len(self.grouper) == len(self.index)
                ):
                    grper = pprint_thing(self.grouper)
                    errmsg = (
                        "Grouper result violates len(labels) == "
                        f"len(data)\nresult: {grper}"
                    )
                    self.grouper = None  # Try for sanity
                    raise AssertionError(errmsg)

        # if we have a date/time-like grouper, make sure that we have
        # Timestamps like
        if getattr(self.grouper, "dtype", None) is not None:
            if is_datetime64_dtype(self.grouper):
                self.grouper = self.grouper.astype("datetime64[ns]")
            elif is_timedelta64_dtype(self.grouper):

                self.grouper = self.grouper.astype("timedelta64[ns]")

    def __repr__(self) -> str:
        return f"Grouping({self.name})"

    def __iter__(self):
        return iter(self.indices)

    _codes: Optional[np.ndarray] = None
    _group_index: Optional[Index] = None

    @property
    def ngroups(self) -> int:
        return len(self.group_index)

    @cache_readonly
    def indices(self):
        # we have a list of groupers
        if isinstance(self.grouper, ops.BaseGrouper):
            return self.grouper.indices

        values = ensure_categorical(self.grouper)
        return values._reverse_indexer()

    @property
    def codes(self) -> np.ndarray:
        if self._codes is None:
            self._make_codes()
        return self._codes

    @cache_readonly
    def result_index(self) -> Index:
        if self.all_grouper is not None:
            return recode_from_groupby(self.all_grouper, self.sort, self.group_index)
        return self.group_index

    @property
    def group_index(self) -> Index:
        if self._group_index is None:
            self._make_codes()
        assert self._group_index is not None
        return self._group_index

    def _make_codes(self) -> None:
        if self._codes is None or self._group_index is None:
            # we have a list of groupers
            if isinstance(self.grouper, ops.BaseGrouper):
                codes = self.grouper.codes_info
                uniques = self.grouper.result_index
            else:
                codes, uniques = algorithms.factorize(self.grouper, sort=self.sort)
                uniques = Index(uniques, name=self.name)
            self._codes = codes
            self._group_index = uniques

    @cache_readonly
    def groups(self) -> Dict[Hashable, np.ndarray]:
        return self.index.groupby(Categorical.from_codes(self.codes, self.group_index))


def get_grouper(
    obj: FrameOrSeries,
    key=None,
    axis: int = 0,
    level=None,
    sort: bool = True,
    observed: bool = False,
    mutated: bool = False,
    validate: bool = True,
) -> "Tuple[ops.BaseGrouper, List[Hashable], FrameOrSeries]":
    """
    Create and return a BaseGrouper, which is an internal
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
    values.

    If validate, then check for key/level overlaps.

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
                    raise ValueError("No group keys passed!")
                else:
                    raise ValueError("multiple levels only valid with MultiIndex")

            if isinstance(level, str):
                if obj._get_axis(axis).name != level:
                    raise ValueError(
                        f"level name {level} is not the name "
                        f"of the {obj._get_axis_name(axis)}"
                    )
            elif level > 0 or level < -1:
                raise ValueError("level > 0 or level < -1 only valid with MultiIndex")

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
            return grouper, [key.key], obj

    # already have a BaseGrouper, just return it
    elif isinstance(key, ops.BaseGrouper):
        return key, [], obj

    if not isinstance(key, list):
        keys = [key]
        match_axis_length = False
    else:
        keys = key
        match_axis_length = len(keys) == len(group_axis)

    # what are we after, exactly?
    any_callable = any(callable(g) or isinstance(g, dict) for g in keys)
    any_groupers = any(isinstance(g, Grouper) for g in keys)
    any_arraylike = any(
        isinstance(g, (list, tuple, Series, Index, np.ndarray)) for g in keys
    )

    # is this an index replacement?
    if (
        not any_callable
        and not any_arraylike
        and not any_groupers
        and match_axis_length
        and level is None
    ):
        if isinstance(obj, DataFrame):
            all_in_columns_index = all(
                g in obj.columns or g in obj.index.names for g in keys
            )
        else:
            assert isinstance(obj, Series)
            all_in_columns_index = all(g in obj.index.names for g in keys)

        if not all_in_columns_index:
            keys = [com.asarray_tuplesafe(keys)]

    if isinstance(level, (tuple, list)):
        if key is None:
            keys = [None] * len(level)
        levels = level
    else:
        levels = [level] * len(keys)

    groupings: List[Grouping] = []
    exclusions: List[Hashable] = []

    # if the actual grouper should be obj[key]
    def is_in_axis(key) -> bool:
        if not _is_label_like(key):
            items = obj._data.items
            try:
                items.get_loc(key)
            except (KeyError, TypeError):
                # TypeError shows up here if we pass e.g. Int64Index
                return False

        return True

    # if the grouper is obj[name]
    def is_in_obj(gpr) -> bool:
        if not hasattr(gpr, "name"):
            return False
        try:
            return gpr is obj[gpr.name]
        except (KeyError, IndexError):
            return False

    for i, (gpr, level) in enumerate(zip(keys, levels)):

        if is_in_obj(gpr):  # df.groupby(df['name'])
            in_axis, name = True, gpr.name
            exclusions.append(name)

        elif is_in_axis(gpr):  # df.groupby('name')
            if gpr in obj:
                if validate:
                    obj._check_label_or_level_ambiguity(gpr, axis=axis)
                in_axis, name, gpr = True, gpr, obj[gpr]
                exclusions.append(name)
            elif obj._is_level_reference(gpr, axis=axis):
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
                f"Length of grouper ({len(gpr)}) and axis ({obj.shape[axis]})"
                " must be same length"
            )

        # create the Grouping
        # allow us to passing the actual Grouping as the gpr
        ping = (
            Grouping(
                group_axis,
                gpr,
                obj=obj,
                name=name,
                level=level,
                sort=sort,
                observed=observed,
                in_axis=in_axis,
            )
            if not isinstance(gpr, Grouping)
            else gpr
        )

        groupings.append(ping)

    if len(groupings) == 0 and len(obj):
        raise ValueError("No group keys passed!")
    elif len(groupings) == 0:
        groupings.append(Grouping(Index([], dtype="int"), np.array([], dtype=np.intp)))

    # create the internals grouper
    grouper = ops.BaseGrouper(group_axis, groupings, sort=sort, mutated=mutated)
    return grouper, exclusions, obj


def _is_label_like(val) -> bool:
    return isinstance(val, (str, tuple)) or (val is not None and is_scalar(val))


def _convert_grouper(axis: Index, grouper):
    if isinstance(grouper, dict):
        return grouper.get
    elif isinstance(grouper, Series):
        if grouper.index.equals(axis):
            return grouper._values
        else:
            return grouper.reindex(axis)._values
    elif isinstance(grouper, (list, Series, Index, np.ndarray)):
        if len(grouper) != len(axis):
            raise ValueError("Grouper and axis must be same length")
        return grouper
    else:
        return grouper
