# -*- coding: utf-8 -*-
# cython: profile=False
# cython: boundscheck=False
"""
Cython implementations of functions resembling the stdlib calendar module
"""

cimport cython
from cython cimport Py_ssize_t

import numpy as np
cimport numpy as np
from numpy cimport int64_t, int32_t
np.import_array()


# ----------------------------------------------------------------------
# Constants

# Slightly more performant cython lookups than a 2D table
cdef int32_t* days_per_month_array = [
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 
    31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

# ----------------------------------------------------------------------


@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline int32_t get_days_in_month(int year, int month) nogil:
    return days_per_month_array[12 * is_leapyear(year) + month - 1]


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef monthrange(int64_t year, Py_ssize_t month):
    cdef:
        int32_t days

    if month < 1 or month > 12:
        raise ValueError("bad month number 0; must be 1-12")

    days = get_days_in_month(year, month)
    return (dayofweek(year, month, 1), days)


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision
cdef int dayofweek(int y, int m, int d) nogil:
    """Sakamoto's method, from wikipedia"""
    cdef:
        int day
        int* sakamoto_arr = [0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4]

    y -= m < 3
    day = (y + y / 4 - y / 100 + y / 400 + sakamoto_arr[m - 1] + d) % 7
    # convert to python day
    return (day + 6) % 7


cdef int is_leapyear(int64_t year) nogil:
    """Returns 1 if the given year is a leap year, 0 otherwise."""
    return ((year & 0x3) == 0 and  # year % 4 == 0 
            ((year % 100) != 0 or (year % 400) == 0))
