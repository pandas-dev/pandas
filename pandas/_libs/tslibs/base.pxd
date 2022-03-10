from cpython.datetime cimport datetime


cdef class ABCTimestamp(datetime):
    pass


cdef bytes AM_LOCAL
cdef bytes PM_LOCAL
