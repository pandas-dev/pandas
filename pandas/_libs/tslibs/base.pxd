from cpython.datetime cimport datetime


cdef class ABCTimestamp(datetime):
    pass


cdef str AM_LOCAL
cdef str PM_LOCAL
