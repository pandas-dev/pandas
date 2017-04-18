cimport numpy as np
import numpy as np
import pandas as pd

cimport util
cimport cython
import cython
from numpy cimport *
from tslib import Timestamp

from cpython.object cimport (Py_EQ, Py_NE, Py_GT, Py_LT, Py_GE, Py_LE,
                             PyObject_RichCompare)

import numbers
_VALID_CLOSED = frozenset(['left', 'right', 'both', 'neither'])

cdef class IntervalMixin:
    property closed_left:
        def __get__(self):
            return self.closed == 'left' or self.closed == 'both'

    property closed_right:
        def __get__(self):
            return self.closed == 'right' or self.closed == 'both'

    property open_left:
        def __get__(self):
            return not self.closed_left

    property open_right:
        def __get__(self):
            return not self.closed_right

    property mid:
        def __get__(self):
            try:
                return 0.5 * (self.left + self.right)
            except TypeError:
                # datetime safe version
                return self.left + 0.5 * (self.right - self.left)


cdef _interval_like(other):
    return (hasattr(other, 'left')
            and hasattr(other, 'right')
            and hasattr(other, 'closed'))


cdef class Interval(IntervalMixin):
    """
    Immutable object implementing an Interval, a bounded slice-like interval.

    .. versionadded:: 0.20.0

    Attributes
    ----------
    left, right : values
        Left and right bounds for each interval.
    closed : {'left', 'right', 'both', 'neither'}
        Whether the interval is closed on the left-side, right-side, both or
        neither. Defaults to 'right'.
    """

    cdef readonly object left, right
    cdef readonly str closed

    def __init__(self, left, right, str closed='right'):
        # note: it is faster to just do these checks than to use a special
        # constructor (__cinit__/__new__) to avoid them
        if closed not in _VALID_CLOSED:
            raise ValueError("invalid option for 'closed': %s" % closed)
        if not left <= right:
            raise ValueError('left side of interval must be <= right side')
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
            op_str = {Py_LT: '<', Py_LE: '<=', Py_GT: '>', Py_GE: '>='}[op]
            raise TypeError(
                'unorderable types: %s() %s %s()' %
                (type(self).__name__, op_str, type(other).__name__))

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
        return ('%s(%r, %r, closed=%r)' %
                (type(self).__name__, left, right, self.closed))

    def __str__(self):

        left, right = self._repr_base()
        start_symbol = '[' if self.closed_left else '('
        end_symbol = ']' if self.closed_right else ')'
        return '%s%s, %s%s' % (start_symbol, left, right, end_symbol)

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
cpdef intervals_to_interval_bounds(ndarray intervals):
    """
    Parameters
    ----------
    intervals: ndarray object array of Intervals / nulls

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

    left = np.empty(n, dtype=object)
    right = np.empty(n, dtype=object)

    for i in range(len(intervals)):
        interval = intervals[i]
        if util._checknull(interval):
            left[i] = np.nan
            right[i] = np.nan
            continue

        if not isinstance(interval, Interval):
            raise TypeError("type {} with value {} is not an interval".format(
                type(interval), interval))

        left[i] = interval.left
        right[i] = interval.right
        if closed is None:
            closed = interval.closed
        elif closed != interval.closed:
            raise ValueError('intervals must all be closed on the same side')

    return left, right, closed

include "intervaltree.pxi"
