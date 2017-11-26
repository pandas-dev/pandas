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
# The first 12 entries correspond to month lengths for non-leap years.
# The remaining 12 entries give month lengths for leap years
cdef int32_t* days_per_month_array = [
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 
    31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

# ----------------------------------------------------------------------


@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline int32_t get_days_in_month(int year, Py_ssize_t month) nogil:
    """Return the number of days in the given month of the given year.

    Parameters
    ----------
    year : int
    month : int

    Returns
    -------
    days_in_month : int

    Notes
    -----
    Assumes that the arguments are valid.  Passing a month not between 1 and 12
    risks a segfault.
    """
    return days_per_month_array[12 * is_leapyear(year) + month - 1]


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef monthrange(int64_t year, Py_ssize_t month):
    """
    Return tuple containing the weekday of the first day of the month and
    the number of days in the month.

    Parameters
    ----------
    year : int
    month : int

    Returns
    -------
    weekday : int
    days_in_month : int

    Raises
    ------
    ValueError if month is invalid
    """
    cdef:
        int32_t days

    if month < 1 or month > 12:
        raise ValueError("bad month number 0; must be 1-12")

    days = get_days_in_month(year, month)
    return (dayofweek(year, month, 1), days)


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision
cdef int dayofweek(int year, int month, int day) nogil:
    """Find the day of week for the date described by the Y/M/D triple y, m, d
    using Sakamoto's method, from wikipedia.

    https://en.wikipedia.org/wiki/\
    Determination_of_the_day_of_the_week#Sakamoto.27s_methods

    Parameters
    ----------
    year : int
    month : int
    day : int

    Returns
    -------
    weekday : int

    Notes
    -----
    Assumes that y, m, d, represents a valid date.
    """
    cdef:
        int day
        int* sakamoto_arr = [0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4]

    y -= m < 3
    day = (y + y / 4 - y / 100 + y / 400 + sakamoto_arr[m - 1] + d) % 7
    # convert to python day
    return (day + 6) % 7


cdef bint is_leapyear(int64_t year) nogil:
    """Returns 1 if the given year is a leap year, 0 otherwise.

    Parameters
    ----------
    year : int

    Returns
    -------
    is_leap : bool
    """
    return ((year & 0x3) == 0 and  # year % 4 == 0
            ((year % 100) != 0 or (year % 400) == 0))
