import operator

import numpy as np
import pandas as pd

from pandas.core.base import PandasObject, IndexOpsMixin
from pandas.core.common import (_values_from_object, _ensure_platform_int,
                                notnull, is_datetime_or_timedelta_dtype,
                                is_integer_dtype, is_float_dtype)
from pandas.core.index import (Index, _ensure_index, default_pprint,
                               InvalidIndexError, MultiIndex)
from pandas.lib import (Interval, IntervalMixin, IntervalTree,
                        interval_bounds_to_intervals,
                        intervals_to_interval_bounds)
from pandas.util.decorators import cache_readonly
import pandas.core.common as com


_VALID_CLOSED = set(['left', 'right', 'both', 'neither'])


def _get_next_label(label):
    dtype = getattr(label, 'dtype', type(label))
    if isinstance(label, (pd.Timestamp, pd.Timedelta)):
        dtype = 'datetime64'
    if is_datetime_or_timedelta_dtype(dtype):
        return label + np.timedelta64(1, 'ns')
    elif is_integer_dtype(dtype):
        return label + 1
    elif is_float_dtype(dtype):
        return np.nextafter(label, np.infty)
    else:
        raise TypeError('cannot determine next label for type %r'
                        % type(label))


def _get_prev_label(label):
    dtype = getattr(label, 'dtype', type(label))
    if isinstance(label, (pd.Timestamp, pd.Timedelta)):
        dtype = 'datetime64'
    if is_datetime_or_timedelta_dtype(dtype):
        return label - np.timedelta64(1, 'ns')
    elif is_integer_dtype(dtype):
        return label - 1
    elif is_float_dtype(dtype):
        return np.nextafter(label, -np.infty)
    else:
        raise TypeError('cannot determine next label for type %r'
                        % type(label))


def _get_interval_closed_bounds(interval):
    """
    Given an Interval or IntervalIndex, return the corresponding interval with
    closed bounds.
    """
    left, right = interval.left, interval.right
    if interval.open_left:
        left = _get_next_label(left)
    if interval.open_right:
        right = _get_prev_label(right)
    return left, right


class IntervalIndex(IntervalMixin, Index):
    """
    Immutable Index implementing an ordered, sliceable set. IntervalIndex
    represents an Index of intervals that are all closed on the same side.

    .. versionadded:: 0.18

    Properties
    ----------
    left, right : array-like (1-dimensional)
        Left and right bounds for each interval.
    closed : {'left', 'right', 'both', 'neither'}, optional
        Whether the intervals are closed on the left-side, right-side, both or
        neither. Defaults to 'right'.
    name : object, optional
        Name to be stored in the index.
    """
    _typ = 'intervalindex'
    _comparables = ['name']
    _attributes = ['name', 'closed']
    _allow_index_ops = True
    _engine = None # disable it

    def __new__(cls, left, right, closed='right', name=None, fastpath=False):
        # TODO: validation
        result = IntervalMixin.__new__(cls)
        result._left = _ensure_index(left)
        result._right = _ensure_index(right)
        result._closed = closed
        result.name = name
        if not fastpath:
            result._validate()
        result._reset_identity()
        return result

    def _validate(self):
        """Verify that the IntervalIndex is valid.
        """
        # TODO: exclude periods?
        if self.closed not in _VALID_CLOSED:
            raise ValueError("invalid options for 'closed': %s" % self.closed)
        if len(self.left) != len(self.right):
            raise ValueError('left and right must have the same length')
        left_valid = notnull(self.left)
        right_valid = notnull(self.right)
        if not (left_valid == right_valid).all():
            raise ValueError('missing values must be missing in the same '
                             'location both left and right sides')
        if not (self.left[left_valid] <= self.right[left_valid]).all():
            raise ValueError('left side of interval must be <= right side')

    def _simple_new(cls, values, name=None, **kwargs):
        # ensure we don't end up here (this is a superclass method)
        raise NotImplementedError

    def _cleanup(self):
        pass

    @property
    def _engine(self):
        raise NotImplementedError

    @cache_readonly
    def _tree(self):
        return IntervalTree(self.left, self.right, closed=self.closed)

    @property
    def _constructor(self):
        return type(self).from_intervals

    @classmethod
    def from_breaks(cls, breaks, closed='right', name=None):
        """
        Construct an IntervalIndex from an array of splits

        Parameters
        ----------
        breaks : array-like (1-dimensional)
            Left and right bounds for each interval.
        closed : {'left', 'right', 'both', 'neither'}, optional
            Whether the intervals are closed on the left-side, right-side, both
            or neither. Defaults to 'right'.
        name : object, optional
            Name to be stored in the index.

        Examples
        --------

        >>> IntervalIndex.from_breaks([0, 1, 2, 3])
        IntervalIndex(left=[0, 1, 2],
                      right=[1, 2, 3],
                      closed='right')
        """
        return cls(breaks[:-1], breaks[1:], closed, name)

    @classmethod
    def from_intervals(cls, data, name=None):
        """
        Construct an IntervalIndex from a 1d array of Interval objects

        Parameters
        ----------
        data : array-like (1-dimensional)
            Array of Interval objects. All intervals must be closed on the same
            sides.
        name : object, optional
            Name to be stored in the index.

        Examples
        --------

        >>> IntervalIndex.from_intervals([Interval(0, 1), Interval(1, 2)])
        IntervalIndex(left=[0, 1],
                      right=[1, 2],
                      closed='right')

        The generic Index constructor work identically when it infers an array
        of all intervals:

        >>> Index([Interval(0, 1), Interval(1, 2)])
        IntervalIndex(left=[0, 1],
                      right=[1, 2],
                      closed='right')
        """
        data = np.asarray(data)
        left, right, closed = intervals_to_interval_bounds(data)
        return cls(left, right, closed, name)

    @classmethod
    def from_tuples(cls, data, closed='right', name=None):
        left = []
        right = []
        for l, r in data:
            left.append(l)
            right.append(r)
        return cls(np.array(left), np.array(right), closed, name)

    def to_tuples(self):
        return Index(com._asarray_tuplesafe(zip(self.left, self.right)))

    @cache_readonly
    def _multiindex(self):
        return MultiIndex.from_arrays([self.left, self.right],
                                      names=['left', 'right'])

    @property
    def left(self):
        return self._left

    @property
    def right(self):
        return self._right

    @property
    def closed(self):
        return self._closed

    def __len__(self):
        return len(self.left)

    @cache_readonly
    def values(self):
        """Returns the IntervalIndex's data as a numpy array of Interval
        objects (with dtype='object')
        """
        left = np.asarray(self.left)
        right = np.asarray(self.right)
        return interval_bounds_to_intervals(left, right, self.closed)

    def __array__(self, result=None):
        """ the array interface, return my values """
        return self.values

    def __array_wrap__(self, result, context=None):
        # we don't want the superclass implementation
        return result

    def _array_values(self):
        return self.values

    def __reduce__(self):
        return self.__class__, (self.left, self.right, self.closed, self.name)

    def _shallow_copy(self, values=None, name=None):
        name = name if name is not None else self.name
        if values is not None:
            return type(self).from_intervals(values, name=name)
        else:
            return self.copy(name=name)

    def copy(self, deep=False, name=None):
        left = self.left.copy(deep=True) if deep else self.left
        right = self.right.copy(deep=True) if deep else self.right
        name = name if name is not None else self.name
        return type(self)(left, right, closed=self.closed, name=name,
                          fastpath=True)

    @cache_readonly
    def dtype(self):
        return np.dtype('O')

    @cache_readonly
    def mid(self):
        """Returns the mid-point of each interval in the index as an array
        """
        try:
            return Index(0.5 * (self.left.values + self.right.values))
        except TypeError:
            # datetime safe version
            delta = self.right.values - self.left.values
            return Index(self.left.values + 0.5 * delta)

    @cache_readonly
    def is_monotonic_increasing(self):
        return self._multiindex.is_monotonic_increasing

    @cache_readonly
    def is_monotonic_decreasing(self):
        return self._multiindex.is_monotonic_decreasing

    @cache_readonly
    def is_unique(self):
        return self._multiindex.is_unique

    @cache_readonly
    def is_non_overlapping_monotonic(self):
        # must be increasing  (e.g., [0, 1), [1, 2), [2, 3), ... )
        # or decreasing (e.g., [-1, 0), [-2, -1), [-3, -2), ...)
        # we already require left <= right
        return ((self.right[:-1] <= self.left[1:]).all() or
                (self.left[:-1] >= self.right[1:]).all())

    def _convert_scalar_indexer(self, key, kind=None):
        return key

    def _maybe_cast_slice_bound(self, label, side, kind):
        return getattr(self, side)._maybe_cast_slice_bound(label, side, kind)

    def _convert_list_indexer(self, keyarr, kind=None):
        """
        we are passed a list indexer.
        Return our indexer or raise if all of the values are not included in the categories
        """
        locs = self.get_indexer(keyarr)
        # TODO: handle keyarr if it includes intervals
        if (locs == -1).any():
            raise KeyError("a list-indexer must only include existing intervals")

        return locs

    def _check_method(self, method):
        if method is not None:
            raise NotImplementedError(
                'method %r not yet implemented for IntervalIndex' % method)

    def _searchsorted_monotonic(self, label, side, exclude_label=False):
        if not self.is_non_overlapping_monotonic:
            raise KeyError('can only get slices from an IntervalIndex if '
                           'bounds are non-overlapping and all monotonic '
                           'increasing or decreasing')

        if isinstance(label, IntervalMixin):
            raise NotImplementedError

        if ((side == 'left' and self.left.is_monotonic_increasing) or
                (side == 'right' and self.left.is_monotonic_decreasing)):
            sub_idx = self.right
            if self.open_right or exclude_label:
                label = _get_next_label(label)
        else:
            sub_idx = self.left
            if self.open_left or exclude_label:
                label = _get_prev_label(label)

        return sub_idx._searchsorted_monotonic(label, side)

    def _get_loc_only_exact_matches(self, key):
        return self._multiindex._tuple_index.get_loc(key)

    def _find_non_overlapping_monotonic_bounds(self, key):
        if isinstance(key, IntervalMixin):
            start = self._searchsorted_monotonic(
                key.left, 'left', exclude_label=key.open_left)
            stop = self._searchsorted_monotonic(
                key.right, 'right', exclude_label=key.open_right)
        else:
            # scalar
            start = self._searchsorted_monotonic(key, 'left')
            stop = self._searchsorted_monotonic(key, 'right')
        return start, stop

    def get_loc(self, key, method=None):
        self._check_method(method)

        original_key = key

        if self.is_non_overlapping_monotonic:
            if isinstance(key, Interval):
                left = self._maybe_cast_slice_bound(key.left, 'left', None)
                right = self._maybe_cast_slice_bound(key.right, 'right', None)
                key = Interval(left, right, key.closed)
            else:
                key = self._maybe_cast_slice_bound(key, 'left', None)

            start, stop = self._find_non_overlapping_monotonic_bounds(key)

            if start + 1 == stop:
                return start
            elif start < stop:
                return slice(start, stop)
            else:
                raise KeyError(original_key)

        else:
            # use the interval tree
            if isinstance(key, Interval):
                left, right = _get_interval_closed_bounds(key)
                return self._tree.get_loc_interval(left, right)
            else:
                return self._tree.get_loc(key)

    def get_value(self, series, key):
        # this method seems necessary for Series.__getitem__ but I have no idea
        # what it should actually do here...
        loc = self.get_loc(key)  # nb. this can't handle slice objects
        return series.iloc[loc]

    def get_indexer(self, target, method=None, limit=None, tolerance=None):
        self._check_method(method)
        target = _ensure_index(target)

        if self.is_non_overlapping_monotonic:
            start, stop = self._find_non_overlapping_monotonic_bounds(target)

            start_plus_one = start + 1
            if (start_plus_one < stop).any():
                raise ValueError('indexer corresponds to non-unique elements')
            return np.where(start_plus_one == stop, start, -1)

        else:
            if isinstance(target, IntervalIndex):
                raise NotImplementedError(
                    'have not yet implemented get_indexer '
                    'for IntervalIndex indexers')
            else:
                return self._tree.get_indexer(target)

    def delete(self, loc):
        new_left = self.left.delete(loc)
        new_right = self.right.delete(loc)
        return type(self)(new_left, new_right, self.closed, self.name,
                          fastpath=True)

    def insert(self, loc, item):
        if not isinstance(item, Interval):
            raise ValueError('can only insert Interval objects into an '
                             'IntervalIndex')
        if not item.closed == self.closed:
            raise ValueError('inserted item must be closed on the same side '
                             'as the index')
        new_left = self.left.insert(loc, item.left)
        new_right = self.right.insert(loc, item.right)
        return type(self)(new_left, new_right, self.closed, self.name,
                          fastpath=True)

    def _as_like_interval_index(self, other, error_msg):
        self._assert_can_do_setop(other)
        other = _ensure_index(other)
        if (not isinstance(other, IntervalIndex) or
                self.closed != other.closed):
            raise ValueError(error_msg)
        return other

    def append(self, other):
        msg = ('can only append two IntervalIndex objects that are closed on '
               'the same side')
        other = self._as_like_interval_index(other, msg)
        new_left = self.left.append(other.left)
        new_right = self.right.append(other.right)
        if other.name is not None and other.name != self.name:
            name = None
        else:
            name = self.name
        return type(self)(new_left, new_right, self.closed, name,
                          fastpath=True)

    def take(self, indexer, axis=0):
        indexer = com._ensure_platform_int(indexer)
        new_left = self.left.take(indexer)
        new_right = self.right.take(indexer)
        return type(self)(new_left, new_right, self.closed, self.name,
                          fastpath=True)

    def __contains__(self, key):
        try:
            self.get_loc(key)
            return True
        except KeyError:
            return False

    def __getitem__(self, value):
        left = self.left[value]
        right = self.right[value]
        if not isinstance(left, Index):
            return Interval(left, right, self.closed)
        else:
            return type(self)(left, right, self.closed, self.name)

    # __repr__ associated methods are based on MultiIndex

    def _format_attrs(self):
        attrs = [('left', default_pprint(self.left)),
                 ('right', default_pprint(self.right)),
                 ('closed', repr(self.closed))]
        if self.name is not None:
            attrs.append(('name', default_pprint(self.name)))
        return attrs

    def _format_space(self):
        return "\n%s" % (' ' * (len(self.__class__.__name__) + 1))

    def _format_data(self):
        return None

    def argsort(self, *args, **kwargs):
        return np.lexsort((self.right, self.left))

    def equals(self, other):
        if self.is_(other):
            return True
        try:
            return (self.left.equals(other.left)
                    and self.right.equals(other.right)
                    and self.closed == other.closed)
        except AttributeError:
            return False

    def _setop(op_name):
        def func(self, other):
            msg = ('can only do set operations between two IntervalIndex '
                   'objects that are closed on the same side')
            other = self._as_like_interval_index(other, msg)
            result = getattr(self._multiindex, op_name)(other._multiindex)
            result_name = self.name if self.name == other.name else None
            return type(self).from_tuples(result.values, closed=self.closed,
                                          name=result_name)
        return func

    union = _setop('union')
    intersection = _setop('intersection')
    difference = _setop('difference')
    sym_diff = _setop('sym_diff')

    # TODO: arithmetic operations


IntervalIndex._add_logical_methods_disabled()
