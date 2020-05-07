"""
We define base classes that will be inherited by Timestamp, Timedelta, etc
in order to allow for fast isinstance checks without circular dependency issues.

This is analogous to core.dtypes.generic.
"""

from cpython.datetime cimport datetime, timedelta


cdef class ABCTimedelta(timedelta):
    pass


cdef class ABCTimestamp(datetime):
    pass


cdef class ABCPeriod:
    pass


cdef class ABCTick:
    pass


cdef bint is_tick_object(object obj):
    return isinstance(obj, ABCTick)


cdef bint is_period_object(object obj):
    return isinstance(obj, ABCPeriod)
