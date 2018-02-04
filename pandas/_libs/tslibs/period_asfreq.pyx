# -*- coding: utf-8 -*-
cimport cython

import numpy as np
from numpy cimport int64_t

from util cimport INT32_MIN

from period_info cimport (dInfoCalc_SetFromAbsDateTime,
                          dInfoCalc_SetFromAbsDate,
                          dInfoCalc_Leapyear,
                          absdate_from_ymd, monthToQuarter, _ISOWeek)
from period_conversion cimport (get_daytime_conversion_factor, max_value,
                                get_abs_time,
                                get_freq_group, get_freq_group_index)

# ----------------------------------------------------------------------
# Constants

cdef int BASE_YEAR = 1970

cdef enum CALENDARS:
    GREGORIAN_CALENDAR = 1
    JULIAN_CALENDAR = 2

cdef enum OFFSETS:
    ORD_OFFSET = 719163LL   # days until 1970-01-01
    BDAY_OFFSET = 513689LL  # days until 1970-01-01
    WEEK_OFFSET = 102737LL

cdef enum FREQS:
    FR_ANN = 1000      # Annual
    FR_QTR = 2000      # Quarterly - December year end (default quarterly)
    FR_MTH = 3000      # Monthly
    FR_WK = 4000       # Weekly
    FR_BUS = 5000      # Business days
    FR_DAY = 6000      # Daily
    FR_HR = 7000       # Hourly
    FR_MIN = 8000      # Minutely
    FR_SEC = 9000      # Secondly
    FR_MS = 10000      # Millisecondly
    FR_US = 11000      # Microsecondly
    FR_NS = 12000      # Nanosecondly
    FR_UND = -10000    # Undefined

# Table of number of days in a month (0-based, without and with leap)
cdef int64_t[:, :] days_in_month = np.array(
    # Windows builds seem to require super-explicit casting
    [[<int64_t>val for val in row] for row in
     [[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
      [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]]],
    dtype=np.int64)


# ---------------------------------------------------------------
# Code derived from scikits.timeseries

@cython.cdivision
cdef int mod_compat(int x, int m) nogil:
    cdef:
        int result = x % m
    if result < 0:
        return result + m
    return result


@cython.cdivision
cdef int floordiv(int x, int divisor) nogil:
    if x < 0:
        if mod_compat(x, divisor):
            return x / divisor - 1
        else:
            return x / divisor
    else:
        return x / divisor


# --------------------------------------------------------------------
# Date Info Construction

cdef int get_date_info(int64_t ordinal, int freq,
                       date_info *dinfo) nogil except -1:
    cdef:
        int64_t absdate = get_python_ordinal(ordinal, freq)
        double abstime = get_abs_time(freq, absdate - ORD_OFFSET, ordinal)

    while abstime < 0:
        abstime += 86400
        absdate -= 1

    while abstime >= 86400:
        abstime -= 86400
        absdate += 1

    if dInfoCalc_SetFromAbsDateTime(dinfo, absdate, abstime,
                                    GREGORIAN_CALENDAR):
        return -1

    return 0

# ----------------------------------------------------------------------

cdef int64_t get_python_ordinal(int64_t period_ordinal, int freq) nogil:
    """
    Returns the proleptic Gregorian ordinal of the date, as an integer.
    This corresponds to the number of days since Jan., 1st, 1AD.
    When the instance has a frequency less than daily, the proleptic date
    is calculated for the last day of the period.
    """
    cdef:
        asfreq_info af_info
        freq_conv_func toDaily = NULL

    if freq == FR_DAY:
        return period_ordinal + ORD_OFFSET

    toDaily = get_asfreq_func(freq, FR_DAY)
    get_asfreq_info(freq, FR_DAY, &af_info)

    return toDaily(period_ordinal, 'E', &af_info) + ORD_OFFSET


cdef void get_asfreq_info(int fromFreq, int toFreq,
                          asfreq_info *af_info) nogil:
    cdef:
        int fromGroup = get_freq_group(fromFreq)
        int toGroup = get_freq_group(toFreq)

    af_info.intraday_conversion_factor = get_daytime_conversion_factor(
        get_freq_group_index(max_value(fromGroup, FR_DAY)),
        get_freq_group_index(max_value(toGroup, FR_DAY)))

    if fromGroup == FR_WK:
        af_info.from_week_end = calc_week_end(fromFreq, fromGroup)
    elif fromGroup == FR_ANN:
        af_info.from_a_year_end = calc_a_year_end(fromFreq, fromGroup)
    elif fromGroup == FR_QTR:
        af_info.from_q_year_end = calc_a_year_end(fromFreq, fromGroup)

    if toGroup == FR_WK:
        af_info.to_week_end = calc_week_end(toFreq, toGroup)
    elif toGroup == FR_ANN:
        af_info.to_a_year_end = calc_a_year_end(toFreq, toGroup)
    elif toGroup == FR_QTR:
        af_info.to_q_year_end = calc_a_year_end(toFreq, toGroup)


cdef int calc_week_end(int freq, int group) nogil:
    return freq - group


@cython.cdivision
cdef int calc_a_year_end(int freq, int group) nogil:
    cdef:
        int result = (freq - group) % 12
    if result == 0:
        return 12
    else:
        return result

# ----------------------------------------------------------------------

cdef int64_t asfreq(int64_t ordinal, int freq1, int freq2, char relation):
    cdef:
        int64_t val
        freq_conv_func func
        asfreq_info finfo

    func = get_asfreq_func(freq1, freq2)

    get_asfreq_info(freq1, freq2, &finfo)
    val = func(ordinal, relation, &finfo)

    if val == INT32_MIN:
        # // Py_Error(PyExc_ValueError, "Unable to convert to desired
        # // frequency.");
        return INT32_MIN

    return val


cdef freq_conv_func get_asfreq_func(int fromFreq, int toFreq) nogil:
    cdef:
        int fromGroup = get_freq_group(fromFreq)
        int toGroup = get_freq_group(toFreq)

    if fromGroup == FR_UND:
        fromGroup = FR_DAY

    if fromGroup == FR_ANN:
        if toGroup == FR_ANN:
            return asfreq_AtoA
        elif toGroup == FR_QTR:
            return asfreq_AtoQ
        elif toGroup == FR_MTH:
            return asfreq_AtoM
        elif toGroup == FR_WK:
            return asfreq_AtoW
        elif toGroup == FR_BUS:
            return asfreq_AtoB
        elif toGroup in [FR_DAY, FR_HR, FR_MIN, FR_SEC, FR_MS, FR_US, FR_NS]:
            return asfreq_AtoDT
        else:
            return <freq_conv_func>nofunc

    elif fromGroup == FR_QTR:
        if toGroup == FR_ANN:
            return asfreq_QtoA
        elif toGroup == FR_QTR:
            return asfreq_QtoQ
        elif toGroup == FR_MTH:
            return asfreq_QtoM
        elif toGroup == FR_WK:
            return asfreq_QtoW
        elif toGroup == FR_BUS:
            return asfreq_QtoB
        elif toGroup in [FR_DAY, FR_HR, FR_MIN, FR_SEC, FR_MS, FR_US, FR_NS]:
            return asfreq_QtoDT
        else:
            return <freq_conv_func>nofunc

    elif fromGroup == FR_MTH:
        if toGroup == FR_ANN:
            return asfreq_MtoA
        elif toGroup == FR_QTR:
            return asfreq_MtoQ
        elif toGroup == FR_MTH:
            return <freq_conv_func>no_op
        elif toGroup == FR_WK:
            return asfreq_MtoW
        elif toGroup == FR_BUS:
            return asfreq_MtoB
        elif toGroup in [FR_DAY, FR_HR, FR_MIN, FR_SEC, FR_MS, FR_US, FR_NS]:
            return asfreq_MtoDT
        else:
            return <freq_conv_func>nofunc

    elif fromGroup == FR_WK:
        if toGroup == FR_ANN:
            return asfreq_WtoA
        elif toGroup == FR_QTR:
            return asfreq_WtoQ
        elif toGroup == FR_MTH:
            return asfreq_WtoM
        elif toGroup == FR_WK:
            return asfreq_WtoW
        elif toGroup == FR_BUS:
            return asfreq_WtoB
        elif toGroup in [FR_DAY, FR_HR, FR_MIN, FR_SEC, FR_MS, FR_US, FR_NS]:
            return asfreq_WtoDT
        else:
            return <freq_conv_func>nofunc

    elif fromGroup == FR_BUS:
        if toGroup == FR_ANN:
            return asfreq_BtoA
        elif toGroup == FR_QTR:
            return asfreq_BtoQ
        elif toGroup == FR_MTH:
            return asfreq_BtoM
        elif toGroup == FR_WK:
            return asfreq_BtoW
        elif toGroup == FR_BUS:
            return <freq_conv_func>no_op
        elif toGroup in [FR_DAY, FR_HR, FR_MIN, FR_SEC, FR_MS, FR_US, FR_NS]:
            return asfreq_BtoDT
        else:
            return <freq_conv_func>nofunc

    elif fromGroup in [FR_DAY, FR_HR, FR_MIN, FR_SEC, FR_MS, FR_US, FR_NS]:
        if toGroup == FR_ANN:
            return asfreq_DTtoA
        elif toGroup == FR_QTR:
            return asfreq_DTtoQ
        elif toGroup == FR_MTH:
            return asfreq_DTtoM
        elif toGroup == FR_WK:
            return asfreq_DTtoW
        elif toGroup == FR_BUS:
            return asfreq_DTtoB
        elif toGroup in [FR_DAY, FR_HR, FR_MIN, FR_SEC, FR_MS, FR_US, FR_NS]:
            if fromGroup > toGroup:
                return asfreq_DownsampleWithinDay
            else:
                return asfreq_UpsampleWithinDay
        else:
            return <freq_conv_func>nofunc

    else:
        return <freq_conv_func>nofunc


cdef int64_t nofunc(int64_t ordinal, char relation, asfreq_info *af_info):
    return INT32_MIN


cdef int64_t no_op(int64_t ordinal, char relation, asfreq_info *af_info):
    return ordinal


# ---------------------------------------------------------------

cdef int64_t DtoQ_yq(int64_t ordinal, asfreq_info *af_info, int *year,
                     int *quarter) nogil:
    cdef:
        date_info dinfo

    if dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET,
                                GREGORIAN_CALENDAR):
        return INT32_MIN

    if af_info.to_q_year_end != 12:
        dinfo.month -= af_info.to_q_year_end
        if dinfo.month <= 0:
            dinfo.month += 12
        else:
            dinfo.year += 1

        dinfo.quarter = monthToQuarter(dinfo.month)

    year[0] = dinfo.year
    quarter[0] = dinfo.quarter
    return 0


cdef inline int64_t transform_via_day(int64_t ordinal, char relation,
                                      asfreq_info *af_info,
                                      freq_conv_func first_func,
                                      freq_conv_func second_func) nogil:
    cdef:
        int64_t result

    result = (first_func)(ordinal, relation, af_info)
    result = (second_func)(result, relation, af_info)

    return result


cdef inline int64_t upsample_daytime(int64_t ordinal,
                                     asfreq_info *af_info, int atEnd) nogil:
    if atEnd:
        return (ordinal + 1) * af_info.intraday_conversion_factor - 1
    else:
        return ordinal * af_info.intraday_conversion_factor


@cython.cdivision
cdef inline int64_t downsample_daytime(int64_t ordinal,
                                       asfreq_info *af_info, int atEnd) nogil:
    return ordinal / af_info.intraday_conversion_factor


# ----------------------------------------------------------------------
# From Annual

@cython.cdivision
cdef int64_t asfreq_AtoDT(int64_t year, char relation,
                          asfreq_info *af_info) nogil:
    cdef:
        int64_t absdate
        int month = (af_info.from_a_year_end) % 12

    # start from 1970
    year += BASE_YEAR

    month += 1

    if af_info.from_a_year_end != 12:
        year -= 1

    if relation == 'E':
        year += 1

    absdate = absdate_from_ymd(year, month, 1)

    if absdate == INT32_MIN:
        return INT32_MIN

    if relation == 'E':
        absdate -= 1

    return upsample_daytime(absdate - ORD_OFFSET, af_info, relation != 'S')


cdef int64_t asfreq_AtoA(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_AtoDT,
                             asfreq_DTtoA)


cdef int64_t asfreq_AtoQ(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_AtoDT,
                             asfreq_DTtoQ)


cdef int64_t asfreq_AtoM(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_AtoDT,
                             asfreq_DTtoM)


cdef int64_t asfreq_AtoW(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_AtoDT,
                             asfreq_DTtoW)


cdef int64_t asfreq_AtoB(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    cdef:
        date_info dinfo

    if dInfoCalc_SetFromAbsDate(&dinfo,
                                asfreq_AtoDT(ordinal, relation,
                                             af_info) + ORD_OFFSET,
                                GREGORIAN_CALENDAR):
        return INT32_MIN

    if relation == 'S':
        return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week)
    else:
        return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week)


@cython.cdivision
cdef int64_t DtoB_weekday(int64_t absdate) nogil:
    return (((absdate) / 7) * 5) + (absdate) % 7 - BDAY_OFFSET


cdef int64_t DtoB_WeekendToMonday(int64_t absdate, int day_of_week) nogil:
    if day_of_week > 4:
        # change to Monday after weekend
        absdate += (7 - day_of_week)
    return DtoB_weekday(absdate)


cdef int64_t DtoB_WeekendToFriday(int64_t absdate, int day_of_week) nogil:
    if day_of_week > 4:
        # change to friday before weekend
        absdate -= (day_of_week - 4)

    return DtoB_weekday(absdate)


# ----------------------------------------------------------------------
# From Quarterly

cdef void QtoD_ym(int64_t ordinal, int *y, int *m, asfreq_info *af_info) nogil:
    y[0] = floordiv(ordinal, 4) + BASE_YEAR
    m[0] = mod_compat(ordinal, 4) * 3 + 1

    if af_info.from_q_year_end != 12:
        m[0] += af_info.from_q_year_end
        if m[0] > 12:
            m[0] -= 12
        else:
            y[0] -= 1


cdef int64_t asfreq_QtoDT(int64_t ordinal, char relation,
                          asfreq_info *af_info) nogil:
    cdef:
        int64_t absdate
        int y, m

    if relation == 'E':
        ordinal += 1

    QtoD_ym(ordinal, &y, &m, af_info)

    absdate = absdate_from_ymd(y, m, 1)
    if absdate == INT32_MIN:
        return INT32_MIN

    if relation == 'E':
        absdate -= 1

    return upsample_daytime(absdate - ORD_OFFSET, af_info, relation != 'S')


cdef int64_t asfreq_QtoQ(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_QtoDT,
                             asfreq_DTtoQ)


cdef int64_t asfreq_QtoA(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_QtoDT,
                             asfreq_DTtoA)


cdef int64_t asfreq_QtoM(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_QtoDT,
                             asfreq_DTtoM)


cdef int64_t asfreq_QtoW(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_QtoDT,
                             asfreq_DTtoW)


cdef int64_t asfreq_QtoB(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    cdef:
        date_info dinfo

    if dInfoCalc_SetFromAbsDate(&dinfo,
                                asfreq_QtoDT(ordinal, relation,
                                             af_info) + ORD_OFFSET,
                                GREGORIAN_CALENDAR):
        return INT32_MIN

    if relation == 'S':
        return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week)
    else:
        return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week)


# ----------------------------------------------------------------------
# From Monthly

cdef void MtoD_ym(int64_t ordinal, int *y, int *m) nogil:
    y[0] = floordiv(ordinal, 12) + BASE_YEAR
    m[0] = mod_compat(ordinal, 12) + 1


cdef int64_t asfreq_MtoDT(int64_t ordinal, char relation,
                          asfreq_info *af_info) nogil:
    cdef:
        int64_t absdate
        int y, m

    if relation == 'E':
        ordinal += 1

    MtoD_ym(ordinal, &y, &m)
    absdate = absdate_from_ymd(y, m, 1)
    if absdate == INT32_MIN:
        return INT32_MIN

    ordinal = absdate - ORD_OFFSET

    if relation == 'E':
        ordinal -= 1

    return upsample_daytime(ordinal, af_info, relation != 'S')


cdef int64_t asfreq_MtoA(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_MtoDT,
                             asfreq_DTtoA)


cdef int64_t asfreq_MtoQ(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_MtoDT,
                             asfreq_DTtoQ)


cdef int64_t asfreq_MtoW(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_MtoDT,
                             asfreq_DTtoW)


cdef int64_t asfreq_MtoB(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    cdef:
        date_info dinfo

    if dInfoCalc_SetFromAbsDate(&dinfo,
                                asfreq_MtoDT(ordinal, relation,
                                             af_info) + ORD_OFFSET,
                                GREGORIAN_CALENDAR):
        return INT32_MIN

    if relation == 'S':
        return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week)
    else:
        return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week)


# ----------------------------------------------------------------------
# From Weekly

cdef int64_t asfreq_WtoDT(int64_t ordinal, char relation,
                          asfreq_info *af_info) nogil:
    ordinal += WEEK_OFFSET
    if relation != 'S':
        ordinal += 1

    ordinal = ordinal * 7 - 6 + af_info.from_week_end - ORD_OFFSET

    if relation != 'S':
        ordinal -= 1

    return upsample_daytime(ordinal, af_info, relation != 'S')


cdef int64_t asfreq_WtoA(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_WtoDT,
                             asfreq_DTtoA)


cdef int64_t asfreq_WtoQ(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_WtoDT,
                             asfreq_DTtoQ)


cdef int64_t asfreq_WtoM(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_WtoDT,
                             asfreq_DTtoM)


cdef int64_t asfreq_WtoW(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_WtoDT,
                             asfreq_DTtoW)


cdef int64_t asfreq_WtoB(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    cdef:
        date_info dinfo

    if dInfoCalc_SetFromAbsDate(&dinfo,
                                asfreq_WtoDT(ordinal, relation,
                                             af_info) + ORD_OFFSET,
                                GREGORIAN_CALENDAR):
        return INT32_MIN

    if relation == 'S':
        return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week)
    else:
        return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week)


# ----------------------------------------------------------------------
# From Business-Freq

@cython.cdivision
cdef int64_t asfreq_BtoDT(int64_t ordinal, char relation,
                          asfreq_info *af_info) nogil:
    ordinal += BDAY_OFFSET
    ordinal = (((ordinal - 1) / 5) * 7 + mod_compat(ordinal - 1, 5) +
               1 - ORD_OFFSET)

    return upsample_daytime(ordinal, af_info, relation != 'S')


cdef int64_t asfreq_BtoA(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_BtoDT,
                             asfreq_DTtoA)


cdef int64_t asfreq_BtoQ(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_BtoDT,
                             asfreq_DTtoQ)


cdef int64_t asfreq_BtoM(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_BtoDT,
                             asfreq_DTtoM)


cdef int64_t asfreq_BtoW(int64_t ordinal, char relation,
                         asfreq_info *af_info) nogil:
    return transform_via_day(ordinal, relation, af_info, asfreq_BtoDT,
                             asfreq_DTtoW)


# ----------------------------------------------------------------------
# From Daily

cdef int64_t asfreq_DTtoA(int64_t ordinal, char relation,
                          asfreq_info *af_info) nogil:
    cdef:
        date_info dinfo

    ordinal = downsample_daytime(ordinal, af_info, 0)
    if dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET,
                                GREGORIAN_CALENDAR):
        return INT32_MIN

    if dinfo.month > af_info.to_a_year_end:
        return <int64_t>(dinfo.year + 1 - BASE_YEAR)
    else:
        return <int64_t>(dinfo.year - BASE_YEAR)


cdef int64_t asfreq_DTtoQ(int64_t ordinal, char relation,
                          asfreq_info *af_info) nogil:
    cdef:
        int year, quarter

    ordinal = downsample_daytime(ordinal, af_info, 0)

    if DtoQ_yq(ordinal, af_info, &year, &quarter) == INT32_MIN:
        return INT32_MIN

    return <int64_t>((year - BASE_YEAR) * 4 + quarter - 1)


cdef int64_t asfreq_DTtoM(int64_t ordinal, char relation,
                          asfreq_info *af_info) nogil:
    cdef:
        date_info dinfo

    ordinal = downsample_daytime(ordinal, af_info, 0)

    if dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET,
                                GREGORIAN_CALENDAR):
        return INT32_MIN
    return <int64_t>((dinfo.year - BASE_YEAR) * 12 + dinfo.month - 1)


@cython.cdivision
cdef int64_t asfreq_DTtoW(int64_t ordinal, char relation,
                          asfreq_info *af_info) nogil:
    ordinal = downsample_daytime(ordinal, af_info, 0)
    return ((ordinal + ORD_OFFSET - (1 + af_info.to_week_end)) / 7 +
            1 - WEEK_OFFSET)


cdef int64_t asfreq_DTtoB(int64_t ordinal, char relation,
                          asfreq_info *af_info) nogil:
    cdef:
        date_info dinfo

    ordinal = downsample_daytime(ordinal, af_info, 0)

    if dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET,
                                GREGORIAN_CALENDAR):
        return INT32_MIN

    if relation == 'S':
        return DtoB_WeekendToFriday(dinfo.absdate, dinfo.day_of_week)
    else:
        return DtoB_WeekendToMonday(dinfo.absdate, dinfo.day_of_week)


# all intra day calculations are now done within one function
cdef int64_t asfreq_DownsampleWithinDay(int64_t ordinal, char relation,
                                        asfreq_info *af_info) nogil:
    return downsample_daytime(ordinal, af_info, relation == 'E')


cdef int64_t asfreq_UpsampleWithinDay(int64_t ordinal, char relation,
                                      asfreq_info *af_info) nogil:
    return upsample_daytime(ordinal, af_info, relation == 'E')


# ----------------------------------------------------------------------
# Period Accessors

cdef int pqyear(int64_t ordinal, int freq):
    cdef:
        int year, quarter
    if _quarter_year(ordinal, freq, &year, &quarter) == INT32_MIN:
        return INT32_MIN
    return year


cdef int pquarter(int64_t ordinal, int freq):
    cdef:
        int year, quarter
    if _quarter_year(ordinal, freq, &year, &quarter) == INT32_MIN:
        return INT32_MIN
    return quarter


cdef int pday_of_year(int64_t ordinal, int freq):
    cdef:
        date_info dinfo

    if get_date_info(ordinal, freq, &dinfo) == INT32_MIN:
        return INT32_MIN
    return dinfo.day_of_year


cdef int pweek(int64_t ordinal, int freq):
    cdef:
        date_info dinfo

    if get_date_info(ordinal, freq, &dinfo) == INT32_MIN:
        return INT32_MIN
    return _ISOWeek(&dinfo)


cdef int pweekday(int64_t ordinal, int freq):
    cdef:
        date_info dinfo

    if get_date_info(ordinal, freq, &dinfo) == INT32_MIN:
        return INT32_MIN
    return dinfo.day_of_week


cdef int pyear(int64_t ordinal, int freq):
    cdef:
        date_info dinfo

    get_date_info(ordinal, freq, &dinfo)
    return dinfo.year


cdef int pmonth(int64_t ordinal, int freq):
    cdef:
        date_info dinfo

    if get_date_info(ordinal, freq, &dinfo) == INT32_MIN:
        return INT32_MIN
    return dinfo.month


cdef int pday(int64_t ordinal, int freq):
    cdef:
        date_info dinfo

    if get_date_info(ordinal, freq, &dinfo) == INT32_MIN:
        return INT32_MIN
    return dinfo.day


cdef int phour(int64_t ordinal, int freq):
    cdef:
        date_info dinfo

    if get_date_info(ordinal, freq, &dinfo) == INT32_MIN:
        return INT32_MIN
    return dinfo.hour


cdef int pminute(int64_t ordinal, int freq):
    cdef:
        date_info dinfo

    if get_date_info(ordinal, freq, &dinfo) == INT32_MIN:
        return INT32_MIN
    return dinfo.minute


cdef int psecond(int64_t ordinal, int freq):
    cdef:
        date_info dinfo

    if get_date_info(ordinal, freq, &dinfo) == INT32_MIN:
        return INT32_MIN
    return <int>dinfo.second


@cython.boundscheck(False)
cdef int pdays_in_month(int64_t ordinal, int freq):
    cdef:
        date_info dinfo
        int days
        Py_ssize_t leap

    if get_date_info(ordinal, freq, &dinfo) == INT32_MIN:
        return INT32_MIN

    leap = <Py_ssize_t>dInfoCalc_Leapyear(dinfo.year, dinfo.calendar)
    days = days_in_month[leap][dinfo.month - 1]
    return days


@cython.cdivision
cdef int _quarter_year(int64_t ordinal, int freq, int *year, int *quarter):
    cdef:
        asfreq_info af_info
        int qtr_freq

    ordinal = get_python_ordinal(ordinal, freq) - ORD_OFFSET

    if get_freq_group(freq) == FR_QTR:
        qtr_freq = freq
    else:
        qtr_freq = FR_QTR

    get_asfreq_info(FR_DAY, qtr_freq, &af_info)

    if DtoQ_yq(ordinal, &af_info, year, quarter) == INT32_MIN:
        return INT32_MIN

    if (qtr_freq % 1000) > 12:
        year[0] -= 1

    return 0
