# -*- coding: utf-8 -*-
import numbers

from cpython.object cimport (Py_EQ, Py_NE, Py_GT, Py_LT, Py_GE, Py_LE,
                             PyObject_RichCompare)

cimport cython
from cython cimport Py_ssize_t

import numpy as np
from numpy cimport ndarray


cimport util
util.import_array()

from tslibs import Timestamp
from tslibs.timezones cimport tz_compare


_VALID_CLOSED = frozenset(['left', 'right', 'both', 'neither'])


cdef class IntervalMixin(object):

    @property
    def closed_left(self):
        """
        Check if the interval is closed on the left side.

        For the meaning of `closed` and `open` see :class:`~pandas.Interval`.

        Returns
        -------
        bool
            ``True`` if the Interval is closed on the left-side, else
            ``False``.
        """
        return self.closed in ('left', 'both')

    @property
    def closed_right(self):
        """
        Check if the interval is closed on the right side.

        For the meaning of `closed` and `open` see :class:`~pandas.Interval`.

        Returns
        -------
        bool
            ``True`` if the Interval is closed on the left-side, else
            ``False``.
        """
        return self.closed in ('right', 'both')

    @property
    def open_left(self):
        """
        Check if the interval is open on the left side.

        For the meaning of `closed` and `open` see :class:`~pandas.Interval`.

        Returns
        -------
        bool
            ``True`` if the Interval is closed on the left-side, else
            ``False``.
        """
        return not self.closed_left

    @property
    def open_right(self):
        """
        Check if the interval is open on the right side.

        For the meaning of `closed` and `open` see :class:`~pandas.Interval`.

        Returns
        -------
        bool
            ``True`` if the Interval is closed on the left-side, else
            ``False``.
        """
        return not self.closed_right

    @property
    def mid(self):
        """
        Return the midpoint of the Interval
        """
        try:
            return 0.5 * (self.left + self.right)
        except TypeError:
            # datetime safe version
            return self.left + 0.5 * self.length

    @property
    def length(self):
        """Return the length of the Interval"""
        try:
            return self.right - self.left
        except TypeError:
            # length not defined for some types, e.g. string
            msg = 'cannot compute length between {left!r} and {right!r}'
            raise TypeError(msg.format(left=self.left, right=self.right))

    def _check_closed_matches(self, other, name='other'):
        """Check if the closed attribute of `other` matches.

        Note that 'left' and 'right' are considered different from 'both'.

        Parameters
        ----------
        other : Interval, IntervalIndex, IntervalArray
        name : str
            Name to use for 'other' in the error message.

        Raises
        ------
        ValueError
            When `other` is not closed exactly the same as self.
        """
        if self.closed != other.closed:
            msg = "'{}.closed' is '{}', expected '{}'."
            raise ValueError(msg.format(name, other.closed, self.closed))


cdef _interval_like(other):
    return (hasattr(other, 'left')
            and hasattr(other, 'right')
            and hasattr(other, 'closed'))


cdef class Interval(IntervalMixin):
    """
    Immutable object implementing an Interval, a bounded slice-like interval.

    .. versionadded:: 0.20.0

    Parameters
    ----------
    left : orderable scalar
        Left bound for the interval.
    right : orderable scalar
        Right bound for the interval.
    closed : {'left', 'right', 'both', 'neither'}, default 'right'
        Whether the interval is closed on the left-side, right-side, both or
        neither.
    closed : {'right', 'left', 'both', 'neither'}, default 'right'
        Whether the interval is closed on the left-side, right-side, both or
        neither. See the Notes for more detailed explanation.

    Notes
    -----
    The parameters `left` and `right` must be from the same type, you must be
    able to compare them and they must satisfy ``left <= right``.

    A closed interval (in mathematics denoted by square brackets) contains
    its endpoints, i.e. the closed interval ``[0, 5]`` is characterized by the
    conditions ``0 <= x <= 5``. This is what ``closed='both'`` stands for.
    An open interval (in mathematics denoted by parentheses) does not contain
    its endpoints, i.e. the open interval ``(0, 5)`` is characterized by the
    conditions ``0 < x < 5``. This is what ``closed='neither'`` stands for.
    Intervals can also be half-open or half-closed, i.e. ``[0, 5)`` is
    described by ``0 <= x < 5`` (``closed='left'``) and ``(0, 5]`` is
    described by ``0 < x <= 5`` (``closed='right'``).

    Examples
    --------
    It is possible to build Intervals of different types, like numeric ones:

    >>> iv = pd.Interval(left=0, right=5)
    >>> iv
    Interval(0, 5, closed='right')

    You can check if an element belongs to it

    >>> 2.5 in iv
    True

    You can test the bounds (``closed='right'``, so ``0 < x <= 5``):

    >>> 0 in iv
    False
    >>> 5 in iv
    True
    >>> 0.0001 in iv
    True

    Calculate its length

    >>> iv.length
    5

    You can operate with `+` and `*` over an Interval and the operation
    is applied to each of its bounds, so the result depends on the type
    of the bound elements

    >>> shifted_iv = iv + 3
    >>> shifted_iv
    Interval(3, 8, closed='right')
    >>> extended_iv = iv * 10.0
    >>> extended_iv
    Interval(0.0, 50.0, closed='right')

    To create a time interval you can use Timestamps as the bounds

    >>> year_2017 = pd.Interval(pd.Timestamp('2017-01-01 00:00:00'),
    ...                         pd.Timestamp('2018-01-01 00:00:00'),
    ...                         closed='left')
    >>> pd.Timestamp('2017-01-01 00:00') in year_2017
    True
    >>> year_2017.length
    Timedelta('365 days 00:00:00')

    And also you can create string intervals

    >>> volume_1 = pd.Interval('Ant', 'Dog', closed='both')
    >>> 'Bee' in volume_1
    True

    See Also
    --------
    IntervalIndex : An Index of Interval objects that are all closed on the
        same side.
    cut : Convert continuous data into discrete bins (Categorical
        of Interval objects).
    qcut : Convert continuous data into bins (Categorical of Interval objects)
        based on quantiles.
    Period : Represents a period of time.
    """
    _typ = "interval"

    cdef readonly object left
    """Left bound for the interval"""

    cdef readonly object right
    """Right bound for the interval"""

    cdef readonly str closed
    """
    Whether the interval is closed on the left-side, right-side, both or
    neither
    """

    def __init__(self, left, right, str closed='right'):
        # note: it is faster to just do these checks than to use a special
        # constructor (__cinit__/__new__) to avoid them
        if closed not in _VALID_CLOSED:
            msg = "invalid option for 'closed': {closed}".format(closed=closed)
            raise ValueError(msg)
        if not left <= right:
            raise ValueError('left side of interval must be <= right side')
        if (isinstance(left, Timestamp) and
                not tz_compare(left.tzinfo, right.tzinfo)):
            # GH 18538
            msg = ("left and right must have the same time zone, got "
                   "'{left_tz}' and '{right_tz}'")
            raise ValueError(msg.format(left_tz=left.tzinfo,
                                        right_tz=right.tzinfo))
        self.left = left
        self.right = right
        self.closed = closed

    def __hash__(self):
        return hash((self.left, self.right, self.closed))

    def __contains__(self, key):
        if _interval_like(key):
            raise TypeError('__contains__ not defined for two intervals')
        return ((self.left < key if self.open_left else self.left <= key) and
                (key < self.right if self.open_right else key <= self.right))

    def __richcmp__(self, other, int op):
        if hasattr(other, 'ndim'):
            # let numpy (or IntervalIndex) handle vectorization
            return NotImplemented

        if _interval_like(other):
            self_tuple = (self.left, self.right, self.closed)
            other_tuple = (other.left, other.right, other.closed)
            return PyObject_RichCompare(self_tuple, other_tuple, op)

        # nb. could just return NotImplemented now, but handling this
        # explicitly allows us to opt into the Python 3 behavior, even on
        # Python 2.
        if op == Py_EQ or op == Py_NE:
            return NotImplemented
        else:
            name = type(self).__name__
            other = type(other).__name__
            op_str = {Py_LT: '<', Py_LE: '<=', Py_GT: '>', Py_GE: '>='}[op]
            raise TypeError('unorderable types: {name}() {op} {other}()'
                            .format(name=name, op=op_str, other=other))

    def __reduce__(self):
        args = (self.left, self.right, self.closed)
        return (type(self), args)

    def _repr_base(self):
        left = self.left
        right = self.right

        # TODO: need more general formatting methodology here
        if isinstance(left, Timestamp) and isinstance(right, Timestamp):
            left = left._short_repr
            right = right._short_repr

        return left, right

    def __repr__(self):

        left, right = self._repr_base()
        name = type(self).__name__
        repr_str = '{name}({left!r}, {right!r}, closed={closed!r})'.format(
            name=name, left=left, right=right, closed=self.closed)
        return repr_str

    def __str__(self):

        left, right = self._repr_base()
        start_symbol = '[' if self.closed_left else '('
        end_symbol = ']' if self.closed_right else ')'
        return '{start}{left}, {right}{end}'.format(
            start=start_symbol, left=left, right=right, end=end_symbol)

    def __add__(self, y):
        if isinstance(y, numbers.Number):
            return Interval(self.left + y, self.right + y)
        elif isinstance(y, Interval) and isinstance(self, numbers.Number):
            return Interval(y.left + self, y.right + self)
        return NotImplemented

    def __sub__(self, y):
        if isinstance(y, numbers.Number):
            return Interval(self.left - y, self.right - y)
        return NotImplemented

    def __mul__(self, y):
        if isinstance(y, numbers.Number):
            return Interval(self.left * y, self.right * y)
        elif isinstance(y, Interval) and isinstance(self, numbers.Number):
            return Interval(y.left * self, y.right * self)
        return NotImplemented

    def __div__(self, y):
        if isinstance(y, numbers.Number):
            return Interval(self.left / y, self.right / y)
        return NotImplemented

    def __truediv__(self, y):
        if isinstance(y, numbers.Number):
            return Interval(self.left / y, self.right / y)
        return NotImplemented

    def __floordiv__(self, y):
        if isinstance(y, numbers.Number):
            return Interval(self.left // y, self.right // y)
        return NotImplemented


@cython.wraparound(False)
@cython.boundscheck(False)
def intervals_to_interval_bounds(ndarray intervals,
                                 bint validate_closed=True):
    """
    Parameters
    ----------
    intervals : ndarray
        object array of Intervals / nulls

    validate_closed: boolean, default True
        boolean indicating if all intervals must be closed on the same side.
        Mismatching closed will raise if True, else return None for closed.

    Returns
    -------
    tuples (left: ndarray object array,
            right: ndarray object array,
            closed: str)

    """

    cdef:
        object closed = None, interval
        int64_t n = len(intervals)
        ndarray left, right
        bint seen_closed = False

    left = np.empty(n, dtype=intervals.dtype)
    right = np.empty(n, dtype=intervals.dtype)

    for i in range(len(intervals)):
        interval = intervals[i]
        if interval is None or util.is_nan(interval):
            left[i] = np.nan
            right[i] = np.nan
            continue

        if not isinstance(interval, Interval):
            raise TypeError("type {typ} with value {iv} is not an interval"
                            .format(typ=type(interval), iv=interval))

        left[i] = interval.left
        right[i] = interval.right
        if not seen_closed:
            seen_closed = True
            closed = interval.closed
        elif closed != interval.closed:
            closed = None
            if validate_closed:
                msg = 'intervals must all be closed on the same side'
                raise ValueError(msg)

    return left, right, closed

include "intervaltree.pxi"
