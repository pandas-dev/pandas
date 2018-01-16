import numpy as np

from pandas._libs.interval import (Interval, IntervalMixin,
                                   intervals_to_interval_bounds)
from pandas.compat.numpy import function as nv
from pandas.core.common import _all_not_none
from pandas.core.config import get_option
from pandas.core.dtypes.cast import maybe_convert_platform
from pandas.core.dtypes.common import (_ensure_platform_int,
                                       is_categorical_dtype, is_float_dtype,
                                       is_integer_dtype, is_interval_dtype,
                                       is_scalar, is_string_dtype)
from pandas.core.dtypes.dtypes import IntervalDtype
from pandas.core.dtypes.generic import (ABCDatetimeIndex, ABCPeriodIndex,
                                        ABCSeries)
from pandas.core.dtypes.missing import isna, notna
from pandas.core.extensions import ExtensionArray
from pandas.core.indexes.base import Index, _ensure_index

_VALID_CLOSED = set(['left', 'right', 'both', 'neither'])


class ScalarDataError(TypeError):
    # XXX: this is a "hack" to get the right class name in the error
    # message.
    pass


class IntervalArray(IntervalMixin, ExtensionArray):
    dtype = IntervalDtype()
    ndim = 1
    can_hold_na = True
    _na_value = fill_value = np.nan

    def __new__(cls, data, closed=None, copy=False, dtype=None,
                fastpath=False, verify_integrity=True):

        from pandas.core.indexes.interval import IntervalIndex

        if fastpath:
            return cls._simple_new(data.left, data.right, closed,
                                   copy=copy, verify_integrity=False)

        if isinstance(data, ABCSeries) and is_interval_dtype(data):
            data = data.values
        if isinstance(data, (cls, IntervalIndex)):
            left = data.left
            right = data.right
            closed = data.closed
        else:

            # don't allow scalars
            if is_scalar(data):
                cls._scalar_data_error(data)

            data = maybe_convert_platform_interval(data)
            left, right, infer_closed = intervals_to_interval_bounds(data)

            if _all_not_none(closed, infer_closed) and closed != infer_closed:
                # GH 18421
                msg = ("conflicting values for closed: constructor got "
                       "'{closed}', inferred from data '{infer_closed}'"
                       .format(closed=closed, infer_closed=infer_closed))
                raise ValueError(msg)

            closed = closed or infer_closed

        return cls._simple_new(left, right, closed,
                               copy=copy, verify_integrity=verify_integrity)

    @classmethod
    def _simple_new(cls, left, right, closed=None,
                    copy=False, verify_integrity=True):
        result = IntervalMixin.__new__(cls)

        if closed is None:
            closed = 'right'
        left = _ensure_index(left, copy=copy)
        right = _ensure_index(right, copy=copy)

        # coerce dtypes to match if needed
        if is_float_dtype(left) and is_integer_dtype(right):
            right = right.astype(left.dtype)
        elif is_float_dtype(right) and is_integer_dtype(left):
            left = left.astype(right.dtype)

        if type(left) != type(right):
            msg = ('must not have differing left [{ltype}] and right '
                   '[{rtype}] types')
            raise ValueError(msg.format(ltype=type(left).__name__,
                                        rtype=type(right).__name__))
        elif is_categorical_dtype(left.dtype) or is_string_dtype(left.dtype):
            # GH 19016
            msg = ('category, object, and string subtypes are not supported '
                   'for IntervalIndex')
            raise TypeError(msg)
        elif isinstance(left, ABCPeriodIndex):
            msg = 'Period dtypes are not supported, use a PeriodIndex instead'
            raise ValueError(msg)
        elif (isinstance(left, ABCDatetimeIndex) and
                str(left.tz) != str(right.tz)):
            msg = ("left and right must have the same time zone, got "
                   "'{left_tz}' and '{right_tz}'")
            raise ValueError(msg.format(left_tz=left.tz, right_tz=right.tz))

        result._left = left
        result._right = right
        result._closed = closed
        if verify_integrity:
            result._validate()
        return result

    @classmethod
    def from_breaks(cls, breaks, closed='right', copy=False):
        """
        Construct an IntervalIndex from an array of splits

        Parameters
        ----------
        breaks : array-like (1-dimensional)
            Left and right bounds for each interval.
        closed : {'left', 'right', 'both', 'neither'}, default 'right'
            Whether the intervals are closed on the left-side, right-side, both
            or neither.
        copy : boolean, default False
            copy the data

        Examples
        --------
        >>> pd.IntervalIndex.from_breaks([0, 1, 2, 3])
        IntervalIndex([(0, 1], (1, 2], (2, 3]]
                      closed='right',
                      dtype='interval[int64]')

        See Also
        --------
        interval_range : Function to create a fixed frequency IntervalIndex
        IntervalIndex.from_arrays : Construct an IntervalIndex from a left and
                                    right array
        IntervalIndex.from_intervals : Construct an IntervalIndex from an array
                                       of Interval objects
        IntervalIndex.from_tuples : Construct an IntervalIndex from a
                                    list/array of tuples
        """
        breaks = maybe_convert_platform_interval(breaks)

        return cls.from_arrays(breaks[:-1], breaks[1:], closed, copy=copy)

    @classmethod
    def from_arrays(cls, left, right, closed='right', copy=False):
        """
        Construct an IntervalIndex from a a left and right array

        Parameters
        ----------
        left : array-like (1-dimensional)
            Left bounds for each interval.
        right : array-like (1-dimensional)
            Right bounds for each interval.
        closed : {'left', 'right', 'both', 'neither'}, default 'right'
            Whether the intervals are closed on the left-side, right-side, both
            or neither.
        copy : boolean, default False
            copy the data

        Examples
        --------
        >>> pd.IntervalIndex.from_arrays([0, 1, 2], [1, 2, 3])
        IntervalIndex([(0, 1], (1, 2], (2, 3]]
                      closed='right',
                      dtype='interval[int64]')

        See Also
        --------
        interval_range : Function to create a fixed frequency IntervalIndex
        IntervalIndex.from_breaks : Construct an IntervalIndex from an array of
                                    splits
        IntervalIndex.from_intervals : Construct an IntervalIndex from an array
                                       of Interval objects
        IntervalIndex.from_tuples : Construct an IntervalIndex from a
                                    list/array of tuples
        """
        left = maybe_convert_platform_interval(left)
        right = maybe_convert_platform_interval(right)

        return cls._simple_new(left, right, closed, copy=copy,
                               verify_integrity=True)

    @classmethod
    def from_intervals(cls, data, copy=False):
        """
        Construct an IntervalIndex from a 1d array of Interval objects

        Parameters
        ----------
        data : array-like (1-dimensional)
            Array of Interval objects. All intervals must be closed on the same
            sides.
        copy : boolean, default False
            by-default copy the data, this is compat only and ignored

        Examples
        --------
        >>> pd.IntervalIndex.from_intervals([pd.Interval(0, 1),
        ...                                  pd.Interval(1, 2)])
        IntervalIndex([(0, 1], (1, 2]]
                      closed='right', dtype='interval[int64]')

        The generic Index constructor work identically when it infers an array
        of all intervals:

        >>> pd.Index([pd.Interval(0, 1), pd.Interval(1, 2)])
        IntervalIndex([(0, 1], (1, 2]]
                      closed='right', dtype='interval[int64]')

        See Also
        --------
        interval_range : Function to create a fixed frequency IntervalIndex
        IntervalIndex.from_arrays : Construct an IntervalIndex from a left and
                                    right array
        IntervalIndex.from_breaks : Construct an IntervalIndex from an array of
                                    splits
        IntervalIndex.from_tuples : Construct an IntervalIndex from a
                                    list/array of tuples
        """
        from pandas.core.indexes.interval import IntervalIndex

        if isinstance(data, (cls, IntervalIndex)):
            left, right, closed = data.left, data.right, data.closed
        else:
            data = maybe_convert_platform_interval(data)
            left, right, closed = intervals_to_interval_bounds(data)
        return cls.from_arrays(left, right, closed, copy=False)

    @classmethod
    def from_tuples(cls, data, closed='right', copy=False):
        """
        Construct an IntervalIndex from a list/array of tuples

        Parameters
        ----------
        data : array-like (1-dimensional)
            Array of tuples
        closed : {'left', 'right', 'both', 'neither'}, default 'right'
            Whether the intervals are closed on the left-side, right-side, both
            or neither.
        copy : boolean, default False
            by-default copy the data, this is compat only and ignored

        Examples
        --------
        >>>  pd.IntervalIndex.from_tuples([(0, 1), (1,2)])
        IntervalIndex([(0, 1], (1, 2]],
                      closed='right', dtype='interval[int64]')

        See Also
        --------
        interval_range : Function to create a fixed frequency IntervalIndex
        IntervalIndex.from_arrays : Construct an IntervalIndex from a left and
                                    right array
        IntervalIndex.from_breaks : Construct an IntervalIndex from an array of
                                    splits
        IntervalIndex.from_intervals : Construct an IntervalIndex from an array
                                       of Interval objects
        """
        if len(data):
            left, right = [], []
        else:
            left = right = data

        for d in data:
            if isna(d):
                lhs = rhs = np.nan
            else:
                lhs, rhs = d
            left.append(lhs)
            right.append(rhs)

        # TODO
        # if we have nulls and we previous had *only*
        # integer data, then we have changed the dtype

        return cls.from_arrays(left, right, closed, copy=False)

    def _validate(self):
        """
        Verify that the IntervalIndex is valid.
        """
        if self.closed not in _VALID_CLOSED:
            raise ValueError("invalid options for 'closed': {closed}"
                             .format(closed=self.closed))
        if len(self.left) != len(self.right):
            raise ValueError('left and right must have the same length')
        left_mask = notna(self.left)
        right_mask = notna(self.right)
        if not (left_mask == right_mask).all():
            raise ValueError('missing values must be missing in the same '
                             'location both left and right sides')
        if not (self.left[left_mask] <= self.right[left_mask]).all():
            raise ValueError('left side of interval must be <= right side')
        self._mask = ~left_mask

    # ---------
    # Interface
    # ---------
    def __iter__(self):
        return iter(self.values)

    def __len__(self):
        return len(self.left)

    def __getitem__(self, value):
        mask = self.isna()
        if is_scalar(mask) and mask:
            return self.fill_value

        left = self.left[value]
        right = self.right[value]

        # scalar
        if not isinstance(left, Index):
            return Interval(left, right, self.closed)

        return self._shallow_copy(left, right)

    def _shallow_copy(self, left=None, right=None):
        from pandas.core.indexes.interval import IntervalIndex

        if left is None:

            # no values passed
            # XXX: is ^ right? Or does that mean just left wasn't passed?
            left, right = self.left, self.right

        elif right is None:

            # only single value passed, could be an IntervalIndex
            # or array of Intervals
            if not isinstance(left, (type(self), IntervalIndex)):
                left = type(self).from_intervals(left)

            left, right = left.left, left.right
        else:

            # both left and right are values
            pass

        return self._simple_new(left, right, closed=self.closed,
                                verify_integrity=False)

    @classmethod
    def concat_same_type(cls, to_concat):
        closed = set(interval.closed for interval in to_concat)
        if len(closed) != 1:
            raise ValueError("Intervals must all be closed on the same side.")
        closed = closed.pop()

        # TODO: avoid intermediate list
        left = np.concatenate([interval.left for interval in to_concat])
        right = np.concatenate([interval.right for interval in to_concat])
        return cls._simple_new(left, right, closed=closed, copy=False)

    # TODO: doc
    def copy(self, deep=False):
        left = self.left.copy(deep=True) if deep else self.left
        right = self.right.copy(deep=True) if deep else self.right
        closed = self.closed
        return type(self).from_arrays(left, right, closed=closed)

    def formatting_values(self):
        return self.values

    def get_values(self):
        return self.values

    def isna(self):
        return isna(self.left)

    def nbytes(self):
        # XXX: https://github.com/pandas-dev/pandas/issues/19209
        return self.values.nbytes

    def take(self, indices, axis=0, allow_fill=True, fill_value=None,
             **kwargs):
        nv.validate_take(tuple(), kwargs)
        indices = _ensure_platform_int(indices)
        left, right = self.left, self.right

        if fill_value is None:
            fill_value = self._na_value
        mask = indices == -1

        if not mask.any():
            # we won't change dtype here in this case
            # if we don't need
            allow_fill = False

        taker = lambda x: x.take(indices, allow_fill=allow_fill,
                                 fill_value=fill_value)

        try:
            new_left = taker(left)
            new_right = taker(right)
        except ValueError:

            # we need to coerce; migth have NA's in an
            # integer dtype
            new_left = taker(left.astype(float))
            new_right = taker(right.astype(float))

        return self._shallow_copy(new_left, new_right)

    take_nd = take

    def _format_data(self):

        # TODO: integrate with categorical and make generic
        # name argument is unused here; just for compat with base / categorical
        n = len(self)
        max_seq_items = min((get_option(
            'display.max_seq_items') or n) // 10, 10)

        formatter = str

        if n == 0:
            summary = '[]'
        elif n == 1:
            first = formatter(self[0])
            summary = '[{first}]'.format(first=first)
        elif n == 2:
            first = formatter(self[0])
            last = formatter(self[-1])
            summary = '[{first}, {last}]'.format(first=first, last=last)
        else:

            if n > max_seq_items:
                n = min(max_seq_items // 2, 10)
                head = [formatter(x) for x in self[:n]]
                tail = [formatter(x) for x in self[-n:]]
                summary = '[{head} ... {tail}]'.format(
                    head=', '.join(head), tail=', '.join(tail))
            else:
                head = []
                tail = [formatter(x) for x in self]
                summary = '[{tail}]'.format(tail=', '.join(tail))

        return summary

    def _format_space(self):
        space = ' ' * (len(self.__class__.__name__) + 1)
        return "\n{space}".format(space=space)

    @property
    def left(self):
        """
        Return the left endpoints of each Interval in the IntervalIndex as
        an Index
        """
        return self._left

    @property
    def right(self):
        """
        Return the right endpoints of each Interval in the IntervalIndex as
        an Index
        """
        return self._right

    @property
    def closed(self):
        """
        Whether the intervals are closed on the left-side, right-side, both or
        neither
        """
        return self._closed

    @property
    def length(self):
        """
        Return an Index with entries denoting the length of each Interval in
        the IntervalIndex
        """
        try:
            return self.right - self.left
        except TypeError:
            # length not defined for some types, e.g. string
            msg = ('IntervalIndex contains Intervals without defined length, '
                   'e.g. Intervals with string endpoints')
            raise TypeError(msg)

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__, self._format_data())

    @property
    def values(self):
        """
        Return the IntervalIndex's data as a numpy array of Interval
        objects (with dtype='object')
        """
        left = self.left
        right = self.right
        mask = self.isna()
        closed = self._closed

        result = np.empty(len(left), dtype=object)
        for i in range(len(left)):
            if mask[i]:
                result[i] = np.nan
            else:
                result[i] = Interval(left[i], right[i], closed)
        return result

    @classmethod
    def _scalar_data_error(cls, data):
        # TODO: array-mixin
        raise ScalarDataError(
            '{0}(...) must be called with a collection of some '
            'kind, {1} was passed'.format(cls.__name__, repr(data))
        )

    def slice(self, slicer):
        left = self.left[slicer]
        right = self.right[slicer]
        return self._simple_new(left, right, closed=self.closed,
                                verify_integrity=False)


def maybe_convert_platform_interval(values):
    """
    Try to do platform conversion, with special casing for IntervalIndex.
    Wrapper around maybe_convert_platform that alters the default return
    dtype in certain cases to be compatible with IntervalIndex.  For example,
    empty lists return with integer dtype instead of object dtype, which is
    prohibited for IntervalIndex.

    Parameters
    ----------
    values : array-like

    Returns
    -------
    array
    """
    if isinstance(values, (list, tuple)) and len(values) == 0:
        # GH 19016
        # empty lists/tuples get object dtype by default, but this is not
        # prohibited for IntervalIndex, so coerce to integer instead
        return np.array([], dtype=np.int64)
    return maybe_convert_platform(values)
