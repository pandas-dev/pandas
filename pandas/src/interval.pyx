cimport numpy as np
import numpy as np
import pandas as pd

cimport cython
import cython

from cpython.object cimport (Py_EQ, Py_NE, Py_GT, Py_LT, Py_GE, Py_LE,
                             PyObject_RichCompare)


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
            raise TypeError('unorderable types: %s() %s %s()' %
                            (type(self).__name__, op_str, type(other).__name__))

    def __reduce__(self):
        args = (self.left, self.right, self.closed)
        return (type(self), args)

    def __repr__(self):
        return ('%s(%r, %r, closed=%r)' %
                (type(self).__name__, self.left, self.right, self.closed))

    def __str__(self):
        start_symbol = '[' if self.closed_left else '('
        end_symbol = ']' if self.closed_right else ')'
        return '%s%s, %s%s' % (start_symbol, self.left, self.right, end_symbol)


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef interval_bounds_to_intervals(np.ndarray left, np.ndarray right,
                                   str closed):
    result = np.empty(len(left), dtype=object)
    nulls = pd.isnull(left) | pd.isnull(right)
    result[nulls] = np.nan
    for i in np.flatnonzero(~nulls):
        result[i] = Interval(left[i], right[i], closed)
    return result


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef intervals_to_interval_bounds(np.ndarray intervals):
    left = np.empty(len(intervals), dtype=object)
    right = np.empty(len(intervals), dtype=object)
    cdef str closed = None
    for i in range(len(intervals)):
        interval = intervals[i]
        left[i] = interval.left
        right[i] = interval.right
        if closed is None:
            closed = interval.closed
        elif closed != interval.closed:
            raise ValueError('intervals must all be closed on the same side')
    return left, right, closed
