# -*- coding: utf-8 -*-

from cython cimport Py_ssize_t

from numpy cimport int64_t, int32_t


cdef int dayofweek(int y, int m, int d) nogil
cdef bint is_leapyear(int64_t year) nogil
cpdef int32_t get_days_in_month(int year, Py_ssize_t month) nogil
cpdef int32_t get_week_of_year(int year, int month, int day) nogil
cpdef int32_t get_day_of_year(int year, int month, int day) nogil
