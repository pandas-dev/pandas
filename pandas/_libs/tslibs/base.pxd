from cpython.datetime cimport datetime, timedelta

cdef class ABCTimedelta(timedelta):
    pass


cdef class ABCTimestamp(datetime):
    pass
