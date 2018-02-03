# -*- coding: utf-8 -*-
cimport cython

import numpy as np
from numpy cimport int64_t

from util cimport INT32_MIN
from libc.limits cimport INT_MAX

from period_conversion cimport get_freq_group

cdef double SECONDS_PER_DAY = 86400
cdef int BASE_YEAR = 1970
cdef int BASE_WEEK_TO_DAY_OFFSET = 1  # diff between day 0 and end of week
cdef int DAYS_PER_WEEK = 7
cdef int BUSINESS_DAYS_PER_WEEK = 5

# ----------------------------------------------------------------------
# ccalendar-like functions

@cython.cdivision
cdef inline int monthToQuarter(int month) nogil:
    return ((month - 1) / 3) + 1


@cython.cdivision
cdef int dInfoCalc_YearOffset(int64_t year, int calendar) nogil except? -1:
    """
    Return the year offset, that is the absolute date of the day
    31.12.(year-1) in the given calendar.

    Note:
       For the Julian calendar we shift the absdate (which is measured
       using the Gregorian Epoch) value by two days because the Epoch
       (0001-01-01) in the Julian calendar lies 2 days before the Epoch in
       the Gregorian calendar.
    """
    year -= 1
    if calendar == GREGORIAN_CALENDAR:
        if (year >= 0 or -1 / 4 == -1):  # TODO: DOES THIS CONDITION MAKE SENSE
            return year * 365 + year / 4 - year / 100 + year / 400
        else:
            return (year * 365 + (year - 3) / 4 -
                    (year - 99) / 100 + (year - 399) / 400)
    elif calendar == JULIAN_CALENDAR:
        if (year >= 0 or -1 / 4 == -1):  # TODO: DOES THIS CONDITION MAKE SENSE
            return year * 365 + year / 4 - 2
        else:
            return year * 365 + (year - 3) / 4 - 2
    else:
        return -1
        # raise ValueError("unknown calendar")


@cython.cdivision
cdef int dInfoCalc_DayOfWeek(int64_t absdate) nogil:
    """Return the day of the week for the given absolute date"""
    cdef:
        int day_of_week

    if absdate >= 1:
        day_of_week = (absdate - 1) % 7
    else:
        day_of_week = 6 - ((-absdate) % 7)
    return day_of_week


@cython.cdivision
cdef bint dInfoCalc_Leapyear(int64_t year, int calendar) nogil:
    """ Return 1/0 iff year points to a leap year in calendar."""
    if calendar == GREGORIAN_CALENDAR:
        return (year % 4 == 0) and ((year % 100 != 0) or (year % 400 == 0))
    else:
        return (year % 4 == 0)


# ----------------------------------------------------------------------

@cython.cdivision
cdef int _ISOWeek(date_info *dinfo):
    cdef:
        int week

    # Estimate
    week = (dinfo.day_of_year - 1) - dinfo.day_of_week + 3
    if week >= 0:
        week = week / 7 + 1

    # Verify
    if week < 0:
        # The day lies in last week of the previous year
        if (week > -2) or (week == -2 and
                           dInfoCalc_Leapyear(dinfo.year - 1, dinfo.calendar)):
            week = 53
        else:
            week = 52
    elif week == 53:
        # Check if the week belongs to year or year+1
        if (31 - dinfo.day + dinfo.day_of_week < 3):
            week = 1

    return week


cdef int dInfoCalc_SetFromAbsDateTime(date_info *dinfo,
                                      int64_t absdate, double abstime,
                                      int calendar) nogil except -1:
    """
    Set the instance's value using the given date and time. calendar
    may be set to the flags: GREGORIAN_CALENDAR, JULIAN_CALENDAR to
    indicate the calendar to be used.
    """
    # Bounds check
    if not (abstime >= 0.0 and abstime <= SECONDS_PER_DAY):
        return -1
    # Py_AssertWithArg(abstime >= 0.0 and abstime <= SECONDS_PER_DAY,
    #                 PyExc_ValueError,
    #                 "abstime out of range (0.0 - 86400.0): %f", abstime);

    # Calculate the date
    if dInfoCalc_SetFromAbsDate(dinfo, absdate, calendar):
        return -1

    # Calculate the time
    if dInfoCalc_SetFromAbsTime(dinfo, abstime):
        return -1

    return 0


@cython.cdivision
cdef int dInfoCalc_SetFromAbsTime(date_info *dinfo, double abstime) nogil:
    """Sets the time part of the DateTime object."""
    cdef:
        int inttime
        int hour, minute
        double second

    inttime = <int>abstime
    hour = inttime / 3600
    minute = (inttime % 3600) / 60
    second = abstime - <double>(hour * 3600 + minute * 60)

    dinfo.hour = hour
    dinfo.minute = minute
    dinfo.second = second

    dinfo.abstime = abstime

    return 0


@cython.boundscheck(False)
@cython.cdivision
cdef int dInfoCalc_SetFromAbsDate(date_info *dinfo,
                                  int64_t absdate, int calendar) nogil:
    """
    Sets the date part of the date_info struct using the indicated
    calendar.

    XXX This could also be done using some integer arithmetics rather
       than with this iterative approach...
    """
    cdef:
        int64_t year
        int64_t yearoffset
        int leap, dayoffset, month
        int64_t[:] monthoffset

    # Approximate year
    if calendar == GREGORIAN_CALENDAR:
        year = <int64_t>((<double>absdate) / 365.2425)
    elif calendar == JULIAN_CALENDAR:
        year = <int64_t>((<double>absdate) / 365.25)
    # else:
    #    Py_Error(PyExc_ValueError, "unknown calendar")

    if absdate > 0:
        year += 1

    # Apply corrections to reach the correct year
    while True:
        # Calculate the year offset
        yearoffset = dInfoCalc_YearOffset(year, calendar)
        if yearoffset == INT32_MIN:
            return INT32_MIN

        # Backward correction: absdate must be greater than the yearoffset
        if yearoffset >= absdate:
            year -= 1
            continue

        dayoffset = absdate - yearoffset
        leap = dInfoCalc_Leapyear(year, calendar)

        # Forward correction: non leap years only have 365 days
        if dayoffset > 365 and not leap:
            year += 1
            continue

        break

    dinfo.year = year
    dinfo.calendar = calendar

    # Now iterate to find the month
    monthoffset = month_offset[leap]

    for month in range(1, 13):
        if monthoffset[month] >= dayoffset:
            break

    dinfo.month = month
    dinfo.quarter = monthToQuarter(month)
    dinfo.day = dayoffset - month_offset[leap][month - 1]

    dinfo.day_of_week = dInfoCalc_DayOfWeek(absdate)
    dinfo.day_of_year = dayoffset
    dinfo.absdate = absdate

    return 0


@cython.boundscheck(False)
@cython.cdivision
cdef int dInfoCalc_SetFromDateAndTime(date_info *dinfo, int year,
                                      int month, int day, int hour,
                                      int minute, double second,
                                      int calendar) nogil:
    """
    Set the instance's value using the given date and time. calendar may be set
    to the flags: GREGORIAN_CALENDAR, JULIAN_CALENDAR to indicate the calendar
    to be used. */
    """
    # Calculate the absolute date
    cdef:
        bint leap
        int64_t absdate
        int yearoffset

    # Range check
    if not year > -(INT_MAX / 366) and year < (INT_MAX / 366):
        return 1
        # raise ValueError("year out of range: %i" % year)

    # Is it a leap year?
    leap = dInfoCalc_Leapyear(year, calendar)

    # Negative month values indicate months relative to the years end
    if month < 0:
        month += 13

    if not (month >= 1 and month <= 12):
        return 1
        # raise ValueError("month out of range (1-12): %i" % month)

    # Negative values indicate days relative to the months end
    if day < 0:
        day += <int>days_in_month[leap][month - 1] + 1

    if not (day >= 1 and day <= days_in_month[leap][month - 1]):
        return 1
        # raise ValueError("day out of range: %i" % day)

    yearoffset = dInfoCalc_YearOffset(year, calendar)
    if yearoffset == INT32_MIN:
        return INT32_MIN

    absdate = day + month_offset[leap][month - 1] + yearoffset

    dinfo.absdate = absdate

    dinfo.year = year;
    dinfo.month = month
    dinfo.quarter = ((month - 1) / 3) + 1
    dinfo.day = day

    dinfo.day_of_week = dInfoCalc_DayOfWeek(absdate)
    dinfo.day_of_year = <short>(absdate - yearoffset)

    dinfo.calendar = calendar

    # Calculate the absolute time
    if not (hour >= 0 and hour <= 23):
        return 1
        # raise ValueError("hour out of range (0-23): %i" % hour)
    if not (minute >= 0 and minute <= 59):
        return 1
        # raise ValueError("minute out of range (0-59): %i" % minute)
    if not (second >= <double>0.0 and
            (second < <double>60.0 or
             (hour == 23 and minute == 59 and second < <double>61.0))):
        return 1
        # raise ValueError("second out of range (0.0 - <60.0; <61.0 for "
        #                  "23:59): %f" % second)

    dinfo.abstime = <double>(hour * 3600 + minute * 60) + second

    dinfo.hour = hour
    dinfo.minute = minute
    dinfo.second = second
    return 0


cdef int64_t absdate_from_ymd(int y, int m, int d) nogil:
    cdef:
        date_info tempDate

    if dInfoCalc_SetFromDateAndTime(&tempDate, y, m, d, 0, 0, 0,
                                    GREGORIAN_CALENDAR):
        return INT32_MIN

    return tempDate.absdate


@cython.cdivision
cdef int64_t get_period_ordinal(int year, int month, int day,
                                int hour, int minute, int second,
                                int microseconds, int picoseconds,
                                int freq) nogil except INT32_MIN:
    """generate an ordinal in period space"""
    cdef:
        int64_t absdays, delta, seconds
        int64_t weeks, days
        int64_t ordinal, day_adj
        int freq_group, fmonth, mdiff

    freq_group = get_freq_group(freq)

    if freq == FR_SEC or freq == FR_MS or freq == FR_US or freq == FR_NS:
        absdays = absdate_from_ymd(year, month, day)
        delta = absdays - ORD_OFFSET
        seconds = <int64_t>(delta * 86400 + hour * 3600 + minute * 60 + second)

        if freq == FR_MS:
            return seconds * 1000 + microseconds / 1000

        elif freq == FR_US:
            return seconds * 1000000 + microseconds

        elif freq == FR_NS:
            return (seconds * 1000000000 +
                    microseconds * 1000 + picoseconds / 1000)

        return seconds

    if freq == FR_MIN:
        absdays = absdate_from_ymd(year, month, day)
        delta = absdays - ORD_OFFSET
        return <int64_t>(delta * 1440 + hour * 60 + minute)

    if freq == FR_HR:
        absdays = absdate_from_ymd(year, month, day)
        if absdays == INT32_MIN:
            return INT32_MIN

        delta = (absdays - ORD_OFFSET)
        return <int64_t>(delta * 24 + hour)

    if freq == FR_DAY:
        return <int64_t>(absdate_from_ymd(year, month, day) - ORD_OFFSET)

    if freq == FR_UND:
        return <int64_t>(absdate_from_ymd(year, month, day) - ORD_OFFSET)

    if freq == FR_BUS:
        days = absdate_from_ymd(year, month, day)
        if days == INT32_MIN:
            return INT32_MIN

        # calculate the current week assuming sunday as last day of a week
        weeks = (days - BASE_WEEK_TO_DAY_OFFSET) / DAYS_PER_WEEK
        # calculate the current weekday (in range 1 .. 7)
        delta = (days - BASE_WEEK_TO_DAY_OFFSET) % DAYS_PER_WEEK + 1
        # return the number of business days in full weeks plus the business
        # days in the last - possible partial - week
        if delta <= BUSINESS_DAYS_PER_WEEK:
            return (<int64_t>(weeks * BUSINESS_DAYS_PER_WEEK) +
                    delta - BDAY_OFFSET)
        else:
            return (<int64_t>(weeks * BUSINESS_DAYS_PER_WEEK) +
                    BUSINESS_DAYS_PER_WEEK + 1 - BDAY_OFFSET)

    if freq_group == FR_WK:
        ordinal = <int64_t>absdate_from_ymd(year, month, day)
        if ordinal == INT32_MIN:
            return INT32_MIN

        day_adj = freq - FR_WK
        return (ordinal - (1 + day_adj)) / 7 + 1 - WEEK_OFFSET

    if freq == FR_MTH:
        return (year - BASE_YEAR) * 12 + month - 1

    if freq_group == FR_QTR:
        fmonth = freq - FR_QTR
        if fmonth == 0:
            fmonth = 12

        mdiff = month - fmonth
        if mdiff < 0:
            mdiff += 12
        if month >= fmonth:
            mdiff += 12

        return (year - BASE_YEAR) * 4 + (mdiff - 1) / 3

    if freq_group == FR_ANN:
        fmonth = freq - FR_ANN
        if fmonth == 0:
            fmonth = 12
        if month <= fmonth:
            return year - BASE_YEAR
        else:
            return year - BASE_YEAR + 1

    # Py_Error(PyExc_RuntimeError, "Unable to generate frequency ordinal")
