# -*- coding: utf-8 -*-

from numpy cimport ndarray, uint8_t

cpdef bint checknull(object val)
cpdef bint checknull_old(object val)
cpdef ndarray[uint8_t] isnaobj(ndarray arr)

cdef bint is_null_datetime64(v)
cdef bint is_null_timedelta64(v)
cdef bint is_null_period(v)

cdef class C_NAType:
    pass

cdef C_NAType C_NA
