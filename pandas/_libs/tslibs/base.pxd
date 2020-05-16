from cpython.datetime cimport datetime, timedelta

cdef class ABCTimedelta(timedelta):
    pass


cdef class ABCTimestamp(datetime):
    pass


cdef class ABCTick:
    pass


cdef bint is_tick_object(object obj)
