# -*- coding: utf-8 -*-
# cython: profile=False

from numpy cimport int64_t

cdef int64_t gmtime(object date)
cpdef object to_datetime(int64_t timestamp)
cpdef object to_timestamp(object dt)

cdef _to_i8(object val)
