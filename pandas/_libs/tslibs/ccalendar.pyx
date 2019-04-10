# -*- coding: utf-8 -*-
# cython: boundscheck=False
"""
Cython implementations of functions resembling the stdlib calendar module
"""

import cython

from numpy cimport int64_t, int32_t

from locale import LC_TIME

from pandas._config.localization import set_locale
from pandas._libs.tslibs.strptime import LocaleTime

# ----------------------------------------------------------------------
# Constants

# Slightly more performant cython lookups than a 2D table
# The first 12 entries correspond to month lengths for non-leap years.
# The remaining 12 entries give month lengths for leap years
cdef int32_t* days_per_month_array = [
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
    31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

cdef int* sakamoto_arr = [0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4]

# The first 13 entries give the month days elapsed as of the first of month N
# (or the total number of days in the year for N=13) in non-leap years.
# The remaining 13 entries give the days elapsed in leap years.
cdef int32_t* _month_offset = [
    0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365,
    0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]

# Canonical location for other modules to find name constants
MONTHS = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL',
          'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
# The first blank line is consistent with calendar.month_name in the calendar
# standard library
MONTHS_FULL = ['', 'January', 'February', 'March', 'April', 'May', 'June',
               'July', 'August', 'September', 'October', 'November',
               'December']
MONTH_NUMBERS = {name: num for num, name in enumerate(MONTHS)}
MONTH_ALIASES = {(num + 1): name for num, name in enumerate(MONTHS)}
MONTH_TO_CAL_NUM = {name: num + 1 for num, name in enumerate(MONTHS)}

DAYS = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']
DAYS_FULL = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday',
             'Saturday', 'Sunday']
int_to_weekday = {num: name for num, name in enumerate(DAYS)}
weekday_to_int = {int_to_weekday[key]: key for key in int_to_weekday}

DAY_SECONDS = 86400
HOUR_SECONDS = 3600

# ----------------------------------------------------------------------


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef int32_t get_days_in_month(int year, Py_ssize_t month) nogil:
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
@cython.cdivision
cdef int dayofweek(int y, int m, int d) nogil:
    """Find the day of week for the date described by the Y/M/D triple y, m, d
    using Sakamoto's method, from wikipedia.

    0 represents Monday.  See [1]_.

    Parameters
    ----------
    y : int
    m : int
    d : int

    Returns
    -------
    weekday : int

    Notes
    -----
    Assumes that y, m, d, represents a valid date.

    See Also
    --------
    [1] https://docs.python.org/3/library/calendar.html#calendar.weekday

    [2] https://en.wikipedia.org/wiki/\
    Determination_of_the_day_of_the_week#Sakamoto.27s_methods
    """
    cdef:
        int day

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


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef int32_t get_week_of_year(int year, int month, int day) nogil:
    """Return the ordinal week-of-year for the given day.

    Parameters
    ----------
    year : int
    month : int
    day : int

    Returns
    -------
    week_of_year : int32_t

    Notes
    -----
    Assumes the inputs describe a valid date.
    """
    cdef:
        int32_t doy, dow
        int woy

    doy = get_day_of_year(year, month, day)
    dow = dayofweek(year, month, day)

    # estimate
    woy = (doy - 1) - dow + 3
    if woy >= 0:
        woy = woy // 7 + 1

    # verify
    if woy < 0:
        if (woy > -2) or (woy == -2 and is_leapyear(year - 1)):
            woy = 53
        else:
            woy = 52
    elif woy == 53:
        if 31 - day + dow < 3:
            woy = 1

    return woy


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef int32_t get_day_of_year(int year, int month, int day) nogil:
    """Return the ordinal day-of-year for the given day.

    Parameters
    ----------
    year : int
    month : int
    day : int

    Returns
    -------
    day_of_year : int32_t

    Notes
    -----
    Assumes the inputs describe a valid date.
    """
    cdef:
        bint isleap
        int32_t mo_off
        int day_of_year

    isleap = is_leapyear(year)

    mo_off = _month_offset[isleap * 13 + month - 1]

    day_of_year = mo_off + day
    return day_of_year


def get_locale_names(name_type: object, locale: object=None):
    """Returns an array of localized day or month names

    Parameters
    ----------
    name_type : string, attribute of LocaleTime() in which to return localized
        names
    locale : string

    Returns
    -------
    list of locale names
    """
    with set_locale(locale, LC_TIME):
        return getattr(LocaleTime(), name_type)
