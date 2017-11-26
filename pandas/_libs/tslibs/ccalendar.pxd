# -*- coding: utf-8 -*-
# cython: profile=False

from cython cimport Py_ssize_t

from numpy cimport int64_t, int32_t


cpdef monthrange(int64_t year, Py_ssize_t month)

cdef int dayofweek(int y, int m, int m) nogil
cdef bint is_leapyear(int64_t year) nogil
cdef int32_t get_days_in_month(int year, Py_ssize_t month) nogil
