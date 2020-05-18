from cpython.datetime cimport datetime, timedelta

cdef class ABCTimedelta(timedelta):
    pass


cdef class ABCTimestamp(datetime):
    pass


cdef class ABCPeriod:
    pass


cdef bint is_period_object(object obj)
