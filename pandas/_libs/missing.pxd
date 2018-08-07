# -*- coding: utf-8 -*-
# cython: profile=False

from tslibs.nattype cimport is_null_datetimelike

cpdef bint checknull(object val)
cpdef bint checknull_old(object val)

cdef bint is_null_datetime64(v)
cdef bint is_null_timedelta64(v)
cdef bint is_null_period(v)
