# -*- coding: utf-8 -*-
from datetime import datetime, date

from cpython cimport (
    PyUnicode_Check,
    PyObject_RichCompareBool,
    Py_EQ, Py_NE)

from numpy cimport int64_t, import_array, ndarray
import numpy as np
import_array()

from libc.stdlib cimport free, malloc
from libc.time cimport strftime, tm
from libc.string cimport strlen, memset

cimport cython

from cpython.datetime cimport (PyDateTime_Check, PyDelta_Check,
                               PyDateTime_IMPORT)
# import datetime C API
PyDateTime_IMPORT

from np_datetime cimport (npy_datetimestruct, dtstruct_to_dt64,
                          dt64_to_dtstruct,
                          pandas_datetime_to_datetimestruct,
                          NPY_DATETIMEUNIT, NPY_FR_D)

cdef extern from "../src/datetime/np_datetime.h":
    int64_t npy_datetimestruct_to_datetime(NPY_DATETIMEUNIT fr,
                                           npy_datetimestruct *d) nogil

cimport util
from util cimport is_period_object, is_string_object, INT32_MIN

from timestamps import Timestamp
from timezones cimport is_utc, is_tzlocal, get_dst_info
from timedeltas import Timedelta
from timedeltas cimport delta_to_nanoseconds

cimport ccalendar
from ccalendar cimport dayofweek, get_day_of_year, is_leapyear
from ccalendar import MONTH_NUMBERS
from conversion cimport tz_convert_utc_to_tzlocal
from frequencies cimport (get_freq_code, get_base_alias,
                          get_to_timestamp_base, get_freq_str,
                          get_rule_month)
from parsing import parse_time_string
from resolution import Resolution
from nattype import nat_strings, NaT, iNaT
from nattype cimport _nat_scalar_rules, NPY_NAT, is_null_datetimelike
from offsets cimport to_offset
from offsets import _Tick

cdef bint PY2 = str == bytes


ctypedef struct asfreq_info:
    int64_t intraday_conversion_factor
    int is_end
    int to_end
    int from_end

ctypedef int64_t (*freq_conv_func)(int64_t, asfreq_info*) nogil


cdef extern from *:
    """
    /*** FREQUENCY CONSTANTS ***/

    #define FR_ANN 1000      /* Annual */
    #define FR_ANNDEC FR_ANN /* Annual - December year end*/
    #define FR_ANNJAN 1001   /* Annual - January year end*/
    #define FR_ANNFEB 1002   /* Annual - February year end*/
    #define FR_ANNMAR 1003   /* Annual - March year end*/
    #define FR_ANNAPR 1004   /* Annual - April year end*/
    #define FR_ANNMAY 1005   /* Annual - May year end*/
    #define FR_ANNJUN 1006   /* Annual - June year end*/
    #define FR_ANNJUL 1007   /* Annual - July year end*/
    #define FR_ANNAUG 1008   /* Annual - August year end*/
    #define FR_ANNSEP 1009   /* Annual - September year end*/
    #define FR_ANNOCT 1010   /* Annual - October year end*/
    #define FR_ANNNOV 1011   /* Annual - November year end*/

    /* The standard quarterly frequencies with various fiscal year ends
       eg, Q42005 for Q@OCT runs Aug 1, 2005 to Oct 31, 2005 */
    #define FR_QTR 2000      /* Quarterly - December year end (default Q) */
    #define FR_QTRDEC FR_QTR /* Quarterly - December year end */
    #define FR_QTRJAN 2001   /* Quarterly - January year end */
    #define FR_QTRFEB 2002   /* Quarterly - February year end */
    #define FR_QTRMAR 2003   /* Quarterly - March year end */
    #define FR_QTRAPR 2004   /* Quarterly - April year end */
    #define FR_QTRMAY 2005   /* Quarterly - May year end */
    #define FR_QTRJUN 2006   /* Quarterly - June year end */
    #define FR_QTRJUL 2007   /* Quarterly - July year end */
    #define FR_QTRAUG 2008   /* Quarterly - August year end */
    #define FR_QTRSEP 2009   /* Quarterly - September year end */
    #define FR_QTROCT 2010   /* Quarterly - October year end */
    #define FR_QTRNOV 2011   /* Quarterly - November year end */

    #define FR_MTH 3000 /* Monthly */

    #define FR_WK 4000     /* Weekly */
    #define FR_WKSUN FR_WK /* Weekly - Sunday end of week */
    #define FR_WKMON 4001  /* Weekly - Monday end of week */
    #define FR_WKTUE 4002  /* Weekly - Tuesday end of week */
    #define FR_WKWED 4003  /* Weekly - Wednesday end of week */
    #define FR_WKTHU 4004  /* Weekly - Thursday end of week */
    #define FR_WKFRI 4005  /* Weekly - Friday end of week */
    #define FR_WKSAT 4006  /* Weekly - Saturday end of week */

    #define FR_BUS 5000 /* Business days */
    #define FR_DAY 6000 /* Daily */
    #define FR_HR 7000  /* Hourly */
    #define FR_MIN 8000 /* Minutely */
    #define FR_SEC 9000 /* Secondly */
    #define FR_MS 10000 /* Millisecondly */
    #define FR_US 11000 /* Microsecondly */
    #define FR_NS 12000 /* Nanosecondly */

    #define FR_UND -10000 /* Undefined */

    static int64_t daytime_conversion_factor_matrix[7][7] = {
        {1, 24, 1440, 86400, 86400000, 86400000000, 86400000000000},
        {0,  1,   60,  3600,  3600000,  3600000000,  3600000000000},
        {0,  0,   1,     60,    60000,    60000000,    60000000000},
        {0,  0,   0,      1,     1000,     1000000,     1000000000},
        {0,  0,   0,      0,        1,        1000,        1000000},
        {0,  0,   0,      0,        0,           1,           1000},
        {0,  0,   0,      0,        0,           0,              1}};

    int max_value(int a, int b) { return a > b ? a : b; }

    static int min_value(int a, int b) { return a < b ? a : b; }

    npy_int64 get_daytime_conversion_factor(int from_index, int to_index) {
        int row = min_value(from_index, to_index);
        int col = max_value(from_index, to_index);
        // row or col < 6 means frequency strictly lower than Daily, which
        // do not use daytime_conversion_factors
        if (row < 6) {
            return 0;
        } else if (col < 6) {
            return 0;
        }
        return daytime_conversion_factor_matrix[row - 6][col - 6];
    }
    """
    int64_t get_daytime_conversion_factor(int from_index, int to_index) nogil
    int max_value(int left, int right) nogil
    int FR_ANN
    int FR_QTR
    int FR_MTH
    int FR_WK
    int FR_DAY
    int FR_HR
    int FR_MIN
    int FR_SEC
    int FR_MS
    int FR_US
    int FR_NS
    int FR_BUS
    int FR_UND


cdef int64_t nofunc(int64_t ordinal, asfreq_info *af_info):
    return np.iinfo(np.int32).min


cdef int64_t no_op(int64_t ordinal, asfreq_info *af_info):
    return ordinal


cdef freq_conv_func get_asfreq_func(int from_freq, int to_freq) nogil:
    cdef:
        int from_group = get_freq_group(from_freq)
        int to_group = get_freq_group(to_freq)

    if from_group == FR_UND:
        from_group = FR_DAY

    if from_group == FR_BUS:
        if to_group == FR_ANN:
            return <freq_conv_func>asfreq_BtoA
        elif to_group == FR_QTR:
            return <freq_conv_func>asfreq_BtoQ
        elif to_group == FR_MTH:
            return <freq_conv_func>asfreq_BtoM
        elif to_group == FR_WK:
            return <freq_conv_func>asfreq_BtoW
        elif to_group == FR_BUS:
            return <freq_conv_func>no_op
        elif to_group  in [FR_DAY, FR_HR, FR_MIN, FR_SEC, FR_MS, FR_US, FR_NS]:
            return <freq_conv_func>asfreq_BtoDT
        else:
            return <freq_conv_func>nofunc

    elif to_group == FR_BUS:
        if from_group == FR_ANN:
            return <freq_conv_func>asfreq_AtoB
        elif from_group == FR_QTR:
            return <freq_conv_func>asfreq_QtoB
        elif from_group == FR_MTH:
            return <freq_conv_func>asfreq_MtoB
        elif from_group == FR_WK:
            return <freq_conv_func>asfreq_WtoB
        elif from_group in [FR_DAY, FR_HR, FR_MIN, FR_SEC,
                            FR_MS, FR_US, FR_NS]:
            return <freq_conv_func>asfreq_DTtoB
        else:
            return <freq_conv_func>nofunc

    elif from_group == FR_ANN:
        if to_group == FR_ANN:
            return <freq_conv_func>asfreq_AtoA
        elif to_group == FR_QTR:
            return <freq_conv_func>asfreq_AtoQ
        elif to_group == FR_MTH:
            return <freq_conv_func>asfreq_AtoM
        elif to_group == FR_WK:
            return <freq_conv_func>asfreq_AtoW
        elif to_group in [FR_DAY, FR_HR, FR_MIN, FR_SEC, FR_MS, FR_US, FR_NS]:
            return <freq_conv_func>asfreq_AtoDT
        else:
            return <freq_conv_func>nofunc

    elif from_group == FR_QTR:
        if to_group == FR_ANN:
            return <freq_conv_func>asfreq_QtoA
        elif to_group == FR_QTR:
            return <freq_conv_func>asfreq_QtoQ
        elif to_group == FR_MTH:
            return <freq_conv_func>asfreq_QtoM
        elif to_group == FR_WK:
            return <freq_conv_func>asfreq_QtoW
        elif to_group in [FR_DAY, FR_HR, FR_MIN, FR_SEC, FR_MS, FR_US, FR_NS]:
            return <freq_conv_func>asfreq_QtoDT
        else:
            return <freq_conv_func>nofunc

    elif from_group == FR_MTH:
        if to_group == FR_ANN:
            return <freq_conv_func>asfreq_MtoA
        elif to_group == FR_QTR:
            return <freq_conv_func>asfreq_MtoQ
        elif to_group == FR_MTH:
            return <freq_conv_func>no_op
        elif to_group == FR_WK:
            return <freq_conv_func>asfreq_MtoW
        elif to_group in [FR_DAY, FR_HR, FR_MIN, FR_SEC, FR_MS, FR_US, FR_NS]:
            return <freq_conv_func>asfreq_MtoDT
        else:
            return <freq_conv_func>nofunc

    elif from_group == FR_WK:
        if to_group == FR_ANN:
            return <freq_conv_func>asfreq_WtoA
        elif to_group == FR_QTR:
            return <freq_conv_func>asfreq_WtoQ
        elif to_group == FR_MTH:
            return <freq_conv_func>asfreq_WtoM
        elif to_group == FR_WK:
            return <freq_conv_func>asfreq_WtoW
        elif to_group in [FR_DAY, FR_HR, FR_MIN, FR_SEC, FR_MS, FR_US, FR_NS]:
            return <freq_conv_func>asfreq_WtoDT
        else:
            return <freq_conv_func>nofunc

    elif from_group in [FR_DAY, FR_HR, FR_MIN, FR_SEC, FR_MS, FR_US, FR_NS]:
        if to_group == FR_ANN:
            return <freq_conv_func>asfreq_DTtoA
        elif to_group == FR_QTR:
            return <freq_conv_func>asfreq_DTtoQ
        elif to_group == FR_MTH:
            return <freq_conv_func>asfreq_DTtoM
        elif to_group == FR_WK:
            return <freq_conv_func>asfreq_DTtoW
        elif to_group in [FR_DAY, FR_HR, FR_MIN, FR_SEC, FR_MS, FR_US, FR_NS]:
            if from_group > to_group:
                return <freq_conv_func>downsample_daytime
            else:
                return <freq_conv_func>upsample_daytime

        else:
            return <freq_conv_func>nofunc

    else:
        return <freq_conv_func>nofunc


# --------------------------------------------------------------------
# Frequency Conversion Helpers

cdef int64_t DtoB_weekday(int64_t unix_date) nogil:
    return ((unix_date + 4) // 7) * 5 + ((unix_date + 4) % 7) - 4


cdef int64_t DtoB(npy_datetimestruct *dts, int roll_back, int64_t unix_date):
    cdef:
        int day_of_week = dayofweek(dts.year, dts.month, dts.day)

    if roll_back == 1:
        if day_of_week > 4:
            # change to friday before weekend
            unix_date -= (day_of_week - 4)
    else:
        if day_of_week > 4:
            # change to Monday after weekend
            unix_date += (7 - day_of_week)

    return DtoB_weekday(unix_date)


cdef inline int64_t upsample_daytime(int64_t ordinal, asfreq_info *af_info):
    if (af_info.is_end):
        return (ordinal + 1) * af_info.intraday_conversion_factor - 1
    else:
        return ordinal * af_info.intraday_conversion_factor


cdef inline int64_t downsample_daytime(int64_t ordinal, asfreq_info *af_info):
    return ordinal // (af_info.intraday_conversion_factor)


cdef inline int64_t transform_via_day(int64_t ordinal,
                                      asfreq_info *af_info,
                                      freq_conv_func first_func,
                                      freq_conv_func second_func):
    cdef:
        int64_t result

    result = first_func(ordinal, af_info)
    result = second_func(result, af_info)
    return result

# --------------------------------------------------------------------
# Conversion _to_ Daily Freq

cdef void AtoD_ym(int64_t ordinal, int64_t *year,
                  int *month, asfreq_info *af_info):
    year[0] = ordinal + 1970
    month[0] = 1

    if af_info.from_end != 12:
        month[0] += af_info.from_end
        if month[0] > 12:
            #  This case is never reached, but is kept for symmetry
            # with QtoD_ym
            month[0] -= 12
        else:
            year[0] -= 1


cdef int64_t asfreq_AtoDT(int64_t ordinal, asfreq_info *af_info):
    cdef:
        int64_t unix_date, year
        int month

    ordinal += af_info.is_end
    AtoD_ym(ordinal, &year, &month, af_info)

    unix_date = unix_date_from_ymd(year, month, 1)
    unix_date -= af_info.is_end
    return upsample_daytime(unix_date, af_info)


cdef void QtoD_ym(int64_t ordinal, int *year,
                  int *month, asfreq_info *af_info):
    year[0] = ordinal // 4 + 1970
    month[0] = (ordinal % 4) * 3 + 1

    if af_info.from_end != 12:
        month[0] += af_info.from_end
        if month[0] > 12:
            month[0] -= 12
        else:
            year[0] -= 1


cdef int64_t asfreq_QtoDT(int64_t ordinal, asfreq_info *af_info):
    cdef:
        int64_t unix_date
        int year, month

    ordinal += af_info.is_end
    QtoD_ym(ordinal, &year, &month, af_info)

    unix_date = unix_date_from_ymd(year, month, 1)
    unix_date -= af_info.is_end
    return upsample_daytime(unix_date, af_info)


cdef void MtoD_ym(int64_t ordinal, int *year, int *month):
    year[0] = ordinal // 12 + 1970
    month[0] = ordinal % 12 + 1


cdef int64_t asfreq_MtoDT(int64_t ordinal, asfreq_info *af_info):
    cdef:
        int64_t unix_date
        int year, month

    ordinal += af_info.is_end
    MtoD_ym(ordinal, &year, &month)

    unix_date = unix_date_from_ymd(year, month, 1)
    unix_date -= af_info.is_end
    return upsample_daytime(unix_date, af_info)


cdef int64_t asfreq_WtoDT(int64_t ordinal, asfreq_info *af_info):
    ordinal = (ordinal * 7 + af_info.from_end - 4 +
               (7 - 1) * (af_info.is_end - 1))
    return upsample_daytime(ordinal, af_info)


# --------------------------------------------------------------------
# Conversion _to_ BusinessDay Freq

cdef int64_t asfreq_AtoB(int64_t ordinal, asfreq_info *af_info):
    cdef:
        int roll_back
        npy_datetimestruct dts
        int64_t unix_date = asfreq_AtoDT(ordinal, af_info)

    pandas_datetime_to_datetimestruct(unix_date, NPY_FR_D, &dts)
    roll_back = af_info.is_end
    return DtoB(&dts, roll_back, unix_date)


cdef int64_t asfreq_QtoB(int64_t ordinal, asfreq_info *af_info):
    cdef:
        int roll_back
        npy_datetimestruct dts
        int64_t unix_date = asfreq_QtoDT(ordinal, af_info)

    pandas_datetime_to_datetimestruct(unix_date, NPY_FR_D, &dts)
    roll_back = af_info.is_end
    return DtoB(&dts, roll_back, unix_date)


cdef int64_t asfreq_MtoB(int64_t ordinal, asfreq_info *af_info):
    cdef:
        int roll_back
        npy_datetimestruct dts
        int64_t unix_date = asfreq_MtoDT(ordinal, af_info)

    pandas_datetime_to_datetimestruct(unix_date, NPY_FR_D, &dts)
    roll_back = af_info.is_end
    return DtoB(&dts, roll_back, unix_date)


cdef int64_t asfreq_WtoB(int64_t ordinal, asfreq_info *af_info):
    cdef:
        int roll_back
        npy_datetimestruct dts
        int64_t unix_date = asfreq_WtoDT(ordinal, af_info)

    pandas_datetime_to_datetimestruct(unix_date, NPY_FR_D, &dts)
    roll_back = af_info.is_end
    return DtoB(&dts, roll_back, unix_date)


cdef int64_t asfreq_DTtoB(int64_t ordinal, asfreq_info *af_info):
    cdef:
        int roll_back
        npy_datetimestruct dts
        int64_t unix_date = downsample_daytime(ordinal, af_info)

    pandas_datetime_to_datetimestruct(unix_date, NPY_FR_D, &dts)
    # This usage defines roll_back the opposite way from the others
    roll_back = 1 - af_info.is_end
    return DtoB(&dts, roll_back, unix_date)


# ----------------------------------------------------------------------
# Conversion _from_ Daily Freq

cdef int64_t asfreq_DTtoA(int64_t ordinal, asfreq_info *af_info):
    cdef:
        npy_datetimestruct dts

    ordinal = downsample_daytime(ordinal, af_info)
    pandas_datetime_to_datetimestruct(ordinal, NPY_FR_D, &dts)
    if dts.month > af_info.to_end:
        return <int64_t>(dts.year + 1 - 1970)
    else:
        return <int64_t>(dts.year - 1970)


cdef int DtoQ_yq(int64_t ordinal, asfreq_info *af_info, int *year):
    cdef:
        npy_datetimestruct dts
        int quarter

    pandas_datetime_to_datetimestruct(ordinal, NPY_FR_D, &dts)
    # TODO: Another version of this function used
    # date_info_from_days_and_time(&dts, unix_date, 0)
    # instead of pandas_datetime_to_datetimestruct; is one more performant?
    if af_info.to_end != 12:
        dts.month -= af_info.to_end
        if dts.month <= 0:
            dts.month += 12
        else:
            dts.year += 1

    year[0] = dts.year
    quarter = month_to_quarter(dts.month)
    return quarter


cdef int64_t asfreq_DTtoQ(int64_t ordinal, asfreq_info *af_info):
    cdef:
        int year, quarter

    ordinal = downsample_daytime(ordinal, af_info)

    quarter = DtoQ_yq(ordinal, af_info, &year)
    return <int64_t>((year - 1970) * 4 + quarter - 1)


cdef int64_t asfreq_DTtoM(int64_t ordinal, asfreq_info *af_info):
    cdef:
        npy_datetimestruct dts

    ordinal = downsample_daytime(ordinal, af_info)
    pandas_datetime_to_datetimestruct(ordinal, NPY_FR_D, &dts)
    return <int64_t>((dts.year - 1970) * 12 + dts.month - 1)


cdef int64_t asfreq_DTtoW(int64_t ordinal, asfreq_info *af_info):
    ordinal = downsample_daytime(ordinal, af_info)
    return (ordinal + 3 - af_info.to_end) // 7 + 1


# --------------------------------------------------------------------
# Conversion _from_ BusinessDay Freq

cdef int64_t asfreq_BtoDT(int64_t ordinal, asfreq_info *af_info):
    ordinal = ((ordinal + 3) // 5) * 7 + (ordinal + 3) % 5 -3
    return upsample_daytime(ordinal, af_info)


cdef int64_t asfreq_BtoA(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_BtoDT,
                             <freq_conv_func>asfreq_DTtoA)


cdef int64_t asfreq_BtoQ(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_BtoDT,
                             <freq_conv_func>asfreq_DTtoQ)


cdef int64_t asfreq_BtoM(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_BtoDT,
                             <freq_conv_func>asfreq_DTtoM)


cdef int64_t asfreq_BtoW(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_BtoDT,
                             <freq_conv_func>asfreq_DTtoW)


# ----------------------------------------------------------------------
# Conversion _from_ Annual Freq

cdef int64_t asfreq_AtoA(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_AtoDT,
                             <freq_conv_func>asfreq_DTtoA)


cdef int64_t asfreq_AtoQ(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_AtoDT,
                             <freq_conv_func>asfreq_DTtoQ);


cdef int64_t asfreq_AtoM(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_AtoDT,
                             <freq_conv_func>asfreq_DTtoM)


cdef int64_t asfreq_AtoW(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_AtoDT,
                             <freq_conv_func>asfreq_DTtoW)


# ----------------------------------------------------------------------
# Conversion _from_ Quarterly Freq

cdef int64_t asfreq_QtoQ(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_QtoDT,
                             <freq_conv_func>asfreq_DTtoQ)


cdef int64_t asfreq_QtoA(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_QtoDT,
                             <freq_conv_func>asfreq_DTtoA)


cdef int64_t asfreq_QtoM(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_QtoDT,
                             <freq_conv_func>asfreq_DTtoM)


cdef int64_t asfreq_QtoW(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_QtoDT,
                             <freq_conv_func>asfreq_DTtoW)


# ----------------------------------------------------------------------
# Conversion _from_ Monthly Freq

cdef int64_t asfreq_MtoA(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_MtoDT,
                             <freq_conv_func>asfreq_DTtoA)


cdef int64_t asfreq_MtoQ(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_MtoDT,
                             <freq_conv_func>asfreq_DTtoQ)


cdef int64_t asfreq_MtoW(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_MtoDT,
                             <freq_conv_func>asfreq_DTtoW)


# ----------------------------------------------------------------------
# Conversion _from_ Weekly Freq

cdef int64_t asfreq_WtoA(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_WtoDT,
                             <freq_conv_func>asfreq_DTtoA)


cdef int64_t asfreq_WtoQ(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_WtoDT,
                             <freq_conv_func>asfreq_DTtoQ)


cdef int64_t asfreq_WtoM(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_WtoDT,
                             <freq_conv_func>asfreq_DTtoM)


cdef int64_t asfreq_WtoW(int64_t ordinal, asfreq_info *af_info):
    return transform_via_day(ordinal, af_info,
                             <freq_conv_func>asfreq_WtoDT,
                             <freq_conv_func>asfreq_DTtoW)


# ----------------------------------------------------------------------

@cython.cdivision
cdef char* c_strftime(npy_datetimestruct *dts, char *fmt):
    """
    Generate a nice string representation of the period
    object, originally from DateObject_strftime

    Parameters
    ----------
    dts : npy_datetimestruct*
    fmt : char*

    Returns
    -------
    result : char*
    """
    cdef:
        tm c_date
        char *result
        int result_len = strlen(fmt) + 50

    c_date.tm_sec = dts.sec
    c_date.tm_min = dts.min
    c_date.tm_hour = dts.hour
    c_date.tm_mday = dts.day
    c_date.tm_mon = dts.month - 1
    c_date.tm_year = dts.year - 1900
    c_date.tm_wday = (dayofweek(dts.year, dts.month, dts.day) + 1) % 7
    c_date.tm_yday = get_day_of_year(dts.year, dts.month, dts.day) - 1
    c_date.tm_isdst = -1

    result = <char*>malloc(result_len * sizeof(char))

    strftime(result, result_len, fmt, &c_date)

    return result


# ----------------------------------------------------------------------
# Conversion between date_info and npy_datetimestruct

cdef inline int get_freq_group(int freq) nogil:
    return (freq // 1000) * 1000


cdef inline int get_freq_group_index(int freq) nogil:
    return freq // 1000


# Find the unix_date (days elapsed since datetime(1970, 1, 1)
# for the given year/month/day.
# Assumes GREGORIAN_CALENDAR */
cdef int64_t unix_date_from_ymd(int year, int month, int day) nogil:
    # Calculate the absolute date
    cdef:
        npy_datetimestruct dts
        int64_t unix_date

    memset(&dts, 0, sizeof(npy_datetimestruct))
    dts.year = year
    dts.month = month
    dts.day = day
    unix_date = npy_datetimestruct_to_datetime(NPY_FR_D, &dts)
    return unix_date


# specifically _dont_ use cdvision or else ordinals near -1 are assigned to
# incorrect dates GH#19643
@cython.cdivision(False)
cdef int64_t get_period_ordinal(npy_datetimestruct *dts, int freq) nogil:
    """
    Generate an ordinal in period space

    Parameters
    ----------
    dts: npy_datetimestruct*
    freq : int

    Returns
    -------
    period_ordinal : int64_t
    """
    cdef:
        int64_t unix_date, seconds, delta
        int64_t weeks
        int64_t day_adj
        int freq_group, fmonth, mdiff

    freq_group = get_freq_group(freq)

    if freq_group == FR_ANN:
        fmonth = freq - FR_ANN
        if fmonth == 0:
            fmonth = 12

        mdiff = dts.month - fmonth
        if mdiff <= 0:
            return dts.year - 1970
        else:
            return dts.year - 1970 + 1

    elif freq_group == FR_QTR:
        fmonth = freq - FR_QTR
        if fmonth == 0:
            fmonth = 12

        mdiff = dts.month - fmonth
        # TODO: Aren't the next two conditions equivalent to
        # unconditional incrementing?
        if mdiff < 0:
            mdiff += 12
        if dts.month >= fmonth:
            mdiff += 12

        return (dts.year - 1970) * 4 + (mdiff - 1) // 3

    elif freq == FR_MTH:
        return (dts.year - 1970) * 12 + dts.month - 1

    unix_date = npy_datetimestruct_to_datetime(NPY_FR_D, dts)

    if freq >= FR_SEC:
        seconds = unix_date * 86400 + dts.hour * 3600 + dts.min * 60 + dts.sec

        if freq == FR_MS:
            return seconds * 1000 + dts.us // 1000

        elif freq == FR_US:
            return seconds * 1000000 + dts.us

        elif freq == FR_NS:
            return (seconds * 1000000000 +
                    dts.us * 1000 + dts.ps // 1000)

        else:
            return seconds

    elif freq == FR_MIN:
        return unix_date * 1440 + dts.hour * 60 + dts.min

    elif freq == FR_HR:
        return unix_date * 24 + dts.hour

    elif freq == FR_DAY:
        return unix_date

    elif freq == FR_UND:
        return unix_date

    elif freq == FR_BUS:
        # calculate the current week (counting from 1970-01-01) treating
        # sunday as last day of a week
        weeks = (unix_date + 3) // 7
        # calculate the current weekday (in range 1 .. 7)
        delta = (unix_date + 3) % 7 + 1
        # return the number of business days in full weeks plus the business
        # days in the last - possible partial - week
        if delta <= 5:
            return (5 * weeks) + delta - 4
        else:
            return (5 * weeks) + (5 + 1) - 4

    elif freq_group == FR_WK:
        day_adj = freq - FR_WK
        return (unix_date + 3 - day_adj) // 7 + 1

    # raise ValueError


cdef void get_date_info(int64_t ordinal, int freq,
                        npy_datetimestruct *dts) nogil:
    cdef:
        int64_t unix_date
        double abstime

    unix_date = get_unix_date(ordinal, freq)
    abstime = get_abs_time(freq, unix_date, ordinal)

    while abstime < 0:
        abstime += 86400
        unix_date -= 1

    while abstime >= 86400:
        abstime -= 86400
        unix_date += 1

    date_info_from_days_and_time(dts, unix_date, abstime)


cdef int64_t get_unix_date(int64_t period_ordinal, int freq) nogil:
    """
    Returns the proleptic Gregorian ordinal of the date, as an integer.
    This corresponds to the number of days since Jan., 1st, 1970 AD.
    When the instance has a frequency less than daily, the proleptic date
    is calculated for the last day of the period.

    Parameters
    ----------
    period_ordinal : int64_t
    freq : int

    Returns
    -------
    unix_date : int64_t number of days since datetime(1970, 1, 1)
    """
    cdef:
        asfreq_info af_info
        freq_conv_func toDaily = NULL

    if freq == FR_DAY:
        return period_ordinal

    toDaily = get_asfreq_func(freq, FR_DAY)
    get_asfreq_info(freq, FR_DAY, True, &af_info)
    return toDaily(period_ordinal, &af_info)


@cython.cdivision
cdef void date_info_from_days_and_time(npy_datetimestruct *dts,
                                       int64_t unix_date,
                                       double abstime) nogil:
    """
    Set the instance's value using the given date and time.

    Parameters
    ----------
    dts : npy_datetimestruct*
    unix_date : int64_t
        days elapsed since datetime(1970, 1, 1)
    abstime : double
        seconds elapsed since beginning of day described by unix_date

    Notes
    -----
    Updates dts inplace
    """
    cdef:
        int inttime
        int hour, minute
        double second, subsecond_fraction

    # Bounds check
    # The calling function is responsible for ensuring that
    # abstime >= 0.0 and abstime <= 86400

    # Calculate the date
    pandas_datetime_to_datetimestruct(unix_date, NPY_FR_D, dts)

    # Calculate the time
    inttime = <int>abstime
    hour = inttime / 3600
    minute = (inttime % 3600) / 60
    second = abstime - <double>(hour * 3600 + minute * 60)

    dts.hour = hour
    dts.min = minute
    dts.sec = <int>second

    subsecond_fraction = second - dts.sec
    dts.us = int((subsecond_fraction) * 1e6)
    dts.ps = int(((subsecond_fraction) * 1e6 - dts.us) * 1e6)


@cython.cdivision
cdef double get_abs_time(int freq, int64_t unix_date, int64_t ordinal) nogil:
    cdef:
        int freq_index, day_index, base_index
        int64_t per_day, start_ord
        double unit, result

    if freq <= FR_DAY:
        return 0

    freq_index = freq // 1000
    day_index = FR_DAY // 1000
    base_index = FR_SEC // 1000

    per_day = get_daytime_conversion_factor(day_index, freq_index)
    unit = get_daytime_conversion_factor(freq_index, base_index)

    if base_index < freq_index:
        unit = 1 / unit

    start_ord = unix_date * per_day
    result = <double>(unit * (ordinal - start_ord))
    return result


cdef int get_yq(int64_t ordinal, int freq, int *quarter, int *year):
    """
    Find the year and quarter of a Period with the given ordinal and frequency

    Parameters
    ----------
    ordinal : int64_t
    freq : int
    quarter : *int
    year : *int

    Returns
    -------
    qtr_freq : int
        describes the implied quarterly frequency associated with `freq`

    Notes
    -----
    Sets quarter and year inplace
    """
    cdef:
        asfreq_info af_info
        int qtr_freq
        int64_t unix_date

    unix_date = get_unix_date(ordinal, freq)

    if get_freq_group(freq) == FR_QTR:
        qtr_freq = freq
    else:
        qtr_freq = FR_QTR

    assert (qtr_freq % 1000) <= 12
    get_asfreq_info(FR_DAY, qtr_freq, True, &af_info)

    quarter[0] = DtoQ_yq(unix_date, &af_info, year)
    return qtr_freq


cdef inline int month_to_quarter(int month):
    return (month - 1) // 3 + 1


# ----------------------------------------------------------------------
# Period logic


@cython.wraparound(False)
@cython.boundscheck(False)
def dt64arr_to_periodarr(int64_t[:] dtarr, int freq, tz=None):
    """
    Convert array of datetime64 values (passed in as 'i8' dtype) to a set of
    periods corresponding to desired frequency, per period convention.
    """
    cdef:
        int64_t[:] out
        Py_ssize_t i, l
        npy_datetimestruct dts

    l = len(dtarr)

    out = np.empty(l, dtype='i8')

    if tz is None:
        with nogil:
            for i in range(l):
                if dtarr[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue
                dt64_to_dtstruct(dtarr[i], &dts)
                out[i] = get_period_ordinal(&dts, freq)
    else:
        out = localize_dt64arr_to_period(dtarr, freq, tz)
    return out.base  # .base to access underlying np.ndarray


@cython.wraparound(False)
@cython.boundscheck(False)
def periodarr_to_dt64arr(int64_t[:] periodarr, int freq):
    """
    Convert array to datetime64 values from a set of ordinals corresponding to
    periods per period convention.
    """
    cdef:
        int64_t[:] out
        Py_ssize_t i, l

    l = len(periodarr)

    out = np.empty(l, dtype='i8')

    with nogil:
        for i in range(l):
            if periodarr[i] == NPY_NAT:
                out[i] = NPY_NAT
                continue
            out[i] = period_ordinal_to_dt64(periodarr[i], freq)

    return out.base  # .base to access underlying np.ndarray


cpdef int64_t period_asfreq(int64_t ordinal, int freq1, int freq2, bint end):
    """
    Convert period ordinal from one frequency to another, and if upsampling,
    choose to use start ('S') or end ('E') of period.
    """
    cdef:
        int64_t retval
        freq_conv_func func
        asfreq_info af_info

    if ordinal == iNaT:
        return iNaT

    func = get_asfreq_func(freq1, freq2)
    get_asfreq_info(freq1, freq2, end, &af_info)
    retval = func(ordinal, &af_info)

    if retval == INT32_MIN:
        raise ValueError('Frequency conversion failed')

    return retval


cdef void get_asfreq_info(int from_freq, int to_freq,
                          bint is_end, asfreq_info *af_info) nogil:
    """
    Construct the `asfreq_info` object used to convert an ordinal from
    `from_freq` to `to_freq`.

    Parameters
    ----------
    from_freq : int
    to_freq int
    is_end : bool
    af_info : *asfreq_info
    """
    cdef:
        int from_group = get_freq_group(from_freq)
        int to_group = get_freq_group(to_freq)

    af_info.is_end = is_end

    af_info.intraday_conversion_factor = get_daytime_conversion_factor(
        get_freq_group_index(max_value(from_group, FR_DAY)),
        get_freq_group_index(max_value(to_group, FR_DAY)))

    if from_group == FR_WK:
        af_info.from_end = calc_week_end(from_freq, from_group)
    elif from_group == FR_ANN:
        af_info.from_end = calc_a_year_end(from_freq, from_group)
    elif from_group == FR_QTR:
        af_info.from_end = calc_a_year_end(from_freq, from_group)

    if to_group == FR_WK:
        af_info.to_end = calc_week_end(to_freq, to_group)
    elif to_group == FR_ANN:
        af_info.to_end = calc_a_year_end(to_freq, to_group)
    elif to_group == FR_QTR:
        af_info.to_end = calc_a_year_end(to_freq, to_group)


@cython.cdivision
cdef int calc_a_year_end(int freq, int group) nogil:
    cdef:
        int result = (freq - group) % 12
    if result == 0:
        return 12
    else:
        return result


cdef inline int calc_week_end(int freq, int group) nogil:
    return freq - group


def period_asfreq_arr(ndarray[int64_t] arr, int freq1, int freq2, bint end):
    """
    Convert int64-array of period ordinals from one frequency to another, and
    if upsampling, choose to use start ('S') or end ('E') of period.
    """
    cdef:
        int64_t[:] result
        Py_ssize_t i, n
        freq_conv_func func
        asfreq_info af_info
        int64_t val

    n = len(arr)
    result = np.empty(n, dtype=np.int64)

    func = get_asfreq_func(freq1, freq2)
    get_asfreq_info(freq1, freq2, end, &af_info)

    mask = arr == iNaT
    if mask.any():      # NaT process
        for i in range(n):
            val = arr[i]
            if val != iNaT:
                val = func(val, &af_info)
                if val == INT32_MIN:
                    raise ValueError("Unable to convert to desired frequency.")
            result[i] = val
    else:
        for i in range(n):
            val = func(arr[i], &af_info)
            if val == INT32_MIN:
                raise ValueError("Unable to convert to desired frequency.")
            result[i] = val

    return result.base  # .base to access underlying np.ndarray


cpdef int64_t period_ordinal(int y, int m, int d, int h, int min,
                             int s, int us, int ps, int freq):
    """
    Find the ordinal representation of the given datetime components at the
    frequency `freq`.

    Parameters
    ----------
    y : int
    m : int
    d : int
    h : int
    min : int
    s : int
    us : int
    ps : int

    Returns
    -------
    ordinal : int64_t
    """
    cdef:
        npy_datetimestruct dts
    dts.year = y
    dts.month = m
    dts.day = d
    dts.hour = h
    dts.min = min
    dts.sec = s
    dts.us = us
    dts.ps = ps
    return get_period_ordinal(&dts, freq)


cpdef int64_t period_ordinal_to_dt64(int64_t ordinal, int freq) nogil:
    cdef:
        npy_datetimestruct dts

    if ordinal == NPY_NAT:
        return NPY_NAT

    get_date_info(ordinal, freq, &dts)
    return dtstruct_to_dt64(&dts)


def period_format(int64_t value, int freq, object fmt=None):
    cdef:
        int freq_group

    if value == iNaT:
        return repr(NaT)

    if fmt is None:
        freq_group = get_freq_group(freq)
        if freq_group == 1000:    # FR_ANN
            fmt = b'%Y'
        elif freq_group == 2000:  # FR_QTR
            fmt = b'%FQ%q'
        elif freq_group == 3000:  # FR_MTH
            fmt = b'%Y-%m'
        elif freq_group == 4000:  # WK
            left = period_asfreq(value, freq, 6000, 0)
            right = period_asfreq(value, freq, 6000, 1)
            return '%s/%s' % (period_format(left, 6000),
                              period_format(right, 6000))
        elif (freq_group == 5000      # BUS
              or freq_group == 6000):  # DAY
            fmt = b'%Y-%m-%d'
        elif freq_group == 7000:   # HR
            fmt = b'%Y-%m-%d %H:00'
        elif freq_group == 8000:   # MIN
            fmt = b'%Y-%m-%d %H:%M'
        elif freq_group == 9000:   # SEC
            fmt = b'%Y-%m-%d %H:%M:%S'
        elif freq_group == 10000:  # MILLISEC
            fmt = b'%Y-%m-%d %H:%M:%S.%l'
        elif freq_group == 11000:  # MICROSEC
            fmt = b'%Y-%m-%d %H:%M:%S.%u'
        elif freq_group == 12000:  # NANOSEC
            fmt = b'%Y-%m-%d %H:%M:%S.%n'
        else:
            raise ValueError('Unknown freq: %d' % freq)

    return _period_strftime(value, freq, fmt)


cdef list extra_fmts = [(b"%q", b"^`AB`^"),
                        (b"%f", b"^`CD`^"),
                        (b"%F", b"^`EF`^"),
                        (b"%l", b"^`GH`^"),
                        (b"%u", b"^`IJ`^"),
                        (b"%n", b"^`KL`^")]

cdef list str_extra_fmts = ["^`AB`^", "^`CD`^", "^`EF`^",
                            "^`GH`^", "^`IJ`^", "^`KL`^"]

cdef object _period_strftime(int64_t value, int freq, object fmt):
    cdef:
        Py_ssize_t i
        npy_datetimestruct dts
        char *formatted
        object pat, repl, result
        list found_pat = [False] * len(extra_fmts)
        int year, quarter

    if PyUnicode_Check(fmt):
        fmt = fmt.encode('utf-8')

    get_date_info(value, freq, &dts)
    for i in range(len(extra_fmts)):
        pat = extra_fmts[i][0]
        repl = extra_fmts[i][1]
        if pat in fmt:
            fmt = fmt.replace(pat, repl)
            found_pat[i] = True

    formatted = c_strftime(&dts, <char*> fmt)

    result = util.char_to_string(formatted)
    free(formatted)

    for i in range(len(extra_fmts)):
        if found_pat[i]:
            if get_yq(value, freq, &quarter, &year) < 0:
                raise ValueError('Unable to get quarter and year')

            if i == 0:
                repl = '%d' % quarter
            elif i == 1:  # %f, 2-digit year
                repl = '%.2d' % (year % 100)
            elif i == 2:
                repl = '%d' % year
            elif i == 3:
                repl = '%03d' % (value % 1000)
            elif i == 4:
                repl = '%06d' % (value % 1000000)
            elif i == 5:
                repl = '%09d' % (value % 1000000000)

            result = result.replace(str_extra_fmts[i], repl)

    if PY2:
        result = result.decode('utf-8', 'ignore')

    return result


# ----------------------------------------------------------------------
# period accessors

ctypedef int (*accessor)(int64_t ordinal, int freq) except INT32_MIN


cdef int pyear(int64_t ordinal, int freq):
    cdef:
        npy_datetimestruct dts
    get_date_info(ordinal, freq, &dts)
    return dts.year


@cython.cdivision
cdef int pqyear(int64_t ordinal, int freq):
    cdef:
        int year, quarter
    get_yq(ordinal, freq, &quarter, &year)
    return year


cdef int pquarter(int64_t ordinal, int freq):
    cdef:
        int year, quarter
    get_yq(ordinal, freq, &quarter, &year)
    return quarter


cdef int pmonth(int64_t ordinal, int freq):
    cdef:
        npy_datetimestruct dts
    get_date_info(ordinal, freq, &dts)
    return dts.month


cdef int pday(int64_t ordinal, int freq):
    cdef:
        npy_datetimestruct dts
    get_date_info(ordinal, freq, &dts)
    return dts.day


cdef int pweekday(int64_t ordinal, int freq):
    cdef:
        npy_datetimestruct dts
    get_date_info(ordinal, freq, &dts)
    return dayofweek(dts.year, dts.month, dts.day)


cdef int pday_of_year(int64_t ordinal, int freq):
    cdef:
        npy_datetimestruct dts
    get_date_info(ordinal, freq, &dts)
    return get_day_of_year(dts.year, dts.month, dts.day)


cdef int pweek(int64_t ordinal, int freq):
    cdef:
        npy_datetimestruct dts
    get_date_info(ordinal, freq, &dts)
    return ccalendar.get_week_of_year(dts.year, dts.month, dts.day)


cdef int phour(int64_t ordinal, int freq):
    cdef:
        npy_datetimestruct dts
    get_date_info(ordinal, freq, &dts)
    return dts.hour


cdef int pminute(int64_t ordinal, int freq):
    cdef:
        npy_datetimestruct dts
    get_date_info(ordinal, freq, &dts)
    return dts.min


cdef int psecond(int64_t ordinal, int freq):
    cdef:
        npy_datetimestruct dts
    get_date_info(ordinal, freq, &dts)
    return <int>dts.sec


cdef int pdays_in_month(int64_t ordinal, int freq):
    cdef:
        npy_datetimestruct dts
    get_date_info(ordinal, freq, &dts)
    return ccalendar.get_days_in_month(dts.year, dts.month)


def get_period_field_arr(int code, int64_t[:] arr, int freq):
    cdef:
        Py_ssize_t i, sz
        int64_t[:] out
        accessor f

    func = _get_accessor_func(code)
    if func is NULL:
        raise ValueError('Unrecognized period code: %d' % code)

    sz = len(arr)
    out = np.empty(sz, dtype=np.int64)

    for i in range(sz):
        if arr[i] == iNaT:
            out[i] = -1
            continue
        out[i] = func(arr[i], freq)

    return out.base  # .base to access underlying np.ndarray


cdef accessor _get_accessor_func(int code):
    if code == 0:
        return <accessor>pyear
    elif code == 1:
        return <accessor>pqyear
    elif code == 2:
        return <accessor>pquarter
    elif code == 3:
        return <accessor>pmonth
    elif code == 4:
        return <accessor>pday
    elif code == 5:
        return <accessor>phour
    elif code == 6:
        return <accessor>pminute
    elif code == 7:
        return <accessor>psecond
    elif code == 8:
        return <accessor>pweek
    elif code == 9:
        return <accessor>pday_of_year
    elif code == 10:
        return <accessor>pweekday
    elif code == 11:
        return <accessor>pdays_in_month
    return NULL


def extract_ordinals(object[:] values, freq):
    cdef:
        Py_ssize_t i, n = len(values)
        int64_t[:] ordinals = np.empty(n, dtype=np.int64)
        object p

    freqstr = Period._maybe_convert_freq(freq).freqstr

    for i in range(n):
        p = values[i]

        if is_null_datetimelike(p):
            ordinals[i] = iNaT
        else:
            try:
                ordinals[i] = p.ordinal

                if p.freqstr != freqstr:
                    msg = DIFFERENT_FREQ_INDEX.format(freqstr, p.freqstr)
                    raise IncompatibleFrequency(msg)

            except AttributeError:
                p = Period(p, freq=freq)
                if p is NaT:
                    # input may contain NaT-like string
                    ordinals[i] = iNaT
                else:
                    ordinals[i] = p.ordinal

    return ordinals.base  # .base to access underlying np.ndarray


def extract_freq(object[:] values):
    cdef:
        Py_ssize_t i, n = len(values)
        object p

    for i in range(n):
        p = values[i]

        try:
            # now Timestamp / NaT has freq attr
            if is_period_object(p):
                return p.freq
        except AttributeError:
            pass

    raise ValueError('freq not specified and cannot be inferred')


# -----------------------------------------------------------------------
# period helpers

@cython.wraparound(False)
@cython.boundscheck(False)
cdef int64_t[:] localize_dt64arr_to_period(int64_t[:] stamps,
                                           int freq, object tz):
    cdef:
        Py_ssize_t n = len(stamps)
        int64_t[:] result = np.empty(n, dtype=np.int64)
        ndarray[int64_t] trans
        int64_t[:] deltas
        Py_ssize_t[:] pos
        npy_datetimestruct dts
        int64_t local_val

    if is_utc(tz) or tz is None:
        with nogil:
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                dt64_to_dtstruct(stamps[i], &dts)
                result[i] = get_period_ordinal(&dts, freq)

    elif is_tzlocal(tz):
        for i in range(n):
            if stamps[i] == NPY_NAT:
                result[i] = NPY_NAT
                continue
            local_val = tz_convert_utc_to_tzlocal(stamps[i], tz)
            dt64_to_dtstruct(local_val, &dts)
            result[i] = get_period_ordinal(&dts, freq)
    else:
        # Adjust datetime64 timestamp, recompute datetimestruct
        trans, deltas, typ = get_dst_info(tz)

        if typ not in ['pytz', 'dateutil']:
            # static/fixed; in this case we know that len(delta) == 1
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                dt64_to_dtstruct(stamps[i] + deltas[0], &dts)
                result[i] = get_period_ordinal(&dts, freq)
        else:
            pos = trans.searchsorted(stamps, side='right') - 1

            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                dt64_to_dtstruct(stamps[i] + deltas[pos[i]], &dts)
                result[i] = get_period_ordinal(&dts, freq)

    return result


_DIFFERENT_FREQ = "Input has different freq={1} from Period(freq={0})"
DIFFERENT_FREQ_INDEX = ("Input has different freq={1} "
                        "from PeriodIndex(freq={0})")


class IncompatibleFrequency(ValueError):
    pass


cdef class _Period(object):

    cdef readonly:
        int64_t ordinal
        object freq

    _typ = 'period'

    def __cinit__(self, ordinal, freq):
        self.ordinal = ordinal
        self.freq = freq

    @classmethod
    def _maybe_convert_freq(cls, object freq):

        if isinstance(freq, (int, tuple)):
            code, stride = get_freq_code(freq)
            freq = get_freq_str(code, stride)

        freq = to_offset(freq)

        if freq.n <= 0:
            raise ValueError('Frequency must be positive, because it'
                             ' represents span: {0}'.format(freq.freqstr))

        return freq

    @classmethod
    def _from_ordinal(cls, ordinal, freq):
        """
        Fast creation from an ordinal and freq that are already validated!
        """
        if ordinal == iNaT:
            return NaT
        else:
            freq = cls._maybe_convert_freq(freq)
            self = _Period.__new__(cls, ordinal, freq)
            return self

    def __richcmp__(self, other, op):
        if is_period_object(other):
            if other.freq != self.freq:
                msg = _DIFFERENT_FREQ.format(self.freqstr, other.freqstr)
                raise IncompatibleFrequency(msg)
            return PyObject_RichCompareBool(self.ordinal, other.ordinal, op)
        elif other is NaT:
            return _nat_scalar_rules[op]
        # index/series like
        elif hasattr(other, '_typ'):
            return NotImplemented
        else:
            if op == Py_EQ:
                return NotImplemented
            elif op == Py_NE:
                return NotImplemented
            raise TypeError('Cannot compare type %r with type %r' %
                            (type(self).__name__, type(other).__name__))

    def __hash__(self):
        return hash((self.ordinal, self.freqstr))

    def _add_delta(self, other):
        cdef:
            int64_t nanos, offset_nanos

        if (PyDelta_Check(other) or util.is_timedelta64_object(other) or
                isinstance(other, _Tick)):
            offset = to_offset(self.freq.rule_code)
            if isinstance(offset, _Tick):
                nanos = delta_to_nanoseconds(other)
                offset_nanos = delta_to_nanoseconds(offset)
                if nanos % offset_nanos == 0:
                    ordinal = self.ordinal + (nanos // offset_nanos)
                    return Period(ordinal=ordinal, freq=self.freq)
            msg = 'Input cannot be converted to Period(freq={0})'
            raise IncompatibleFrequency(msg.format(self.freqstr))
        elif util.is_offset_object(other):
            freqstr = other.rule_code
            base = get_base_alias(freqstr)
            if base == self.freq.rule_code:
                ordinal = self.ordinal + other.n
                return Period(ordinal=ordinal, freq=self.freq)
            msg = _DIFFERENT_FREQ.format(self.freqstr, other.freqstr)
            raise IncompatibleFrequency(msg)
        else:  # pragma no cover
            return NotImplemented

    def __add__(self, other):
        if is_period_object(self):
            if (PyDelta_Check(other) or util.is_timedelta64_object(other) or
                    util.is_offset_object(other)):
                return self._add_delta(other)
            elif other is NaT:
                return NaT
            elif util.is_integer_object(other):
                ordinal = self.ordinal + other * self.freq.n
                return Period(ordinal=ordinal, freq=self.freq)
            elif (PyDateTime_Check(other) or
                  is_period_object(other) or util.is_datetime64_object(other)):
                # can't add datetime-like
                # GH#17983
                sname = type(self).__name__
                oname = type(other).__name__
                raise TypeError("unsupported operand type(s) for +: '{self}' "
                                "and '{other}'".format(self=sname,
                                                       other=oname))
            else:  # pragma: no cover
                return NotImplemented
        elif is_period_object(other):
            # this can be reached via __radd__ because of cython rules
            return other + self
        else:
            return NotImplemented

    def __sub__(self, other):
        if is_period_object(self):
            if (PyDelta_Check(other) or util.is_timedelta64_object(other) or
                    util.is_offset_object(other)):
                neg_other = -other
                return self + neg_other
            elif util.is_integer_object(other):
                ordinal = self.ordinal - other * self.freq.n
                return Period(ordinal=ordinal, freq=self.freq)
            elif is_period_object(other):
                if other.freq != self.freq:
                    msg = _DIFFERENT_FREQ.format(self.freqstr, other.freqstr)
                    raise IncompatibleFrequency(msg)
                return (self.ordinal - other.ordinal) * self.freq
            elif getattr(other, '_typ', None) == 'periodindex':
                # GH#21314 PeriodIndex - Period returns an object-index
                # of DateOffset objects, for which we cannot use __neg__
                # directly, so we have to apply it pointwise
                return other.__sub__(self).map(lambda x: -x)
            else:  # pragma: no cover
                return NotImplemented
        elif is_period_object(other):
            if self is NaT:
                return NaT
            return NotImplemented
        else:
            return NotImplemented

    def asfreq(self, freq, how='E'):
        """
        Convert Period to desired frequency, either at the start or end of the
        interval

        Parameters
        ----------
        freq : string
        how : {'E', 'S', 'end', 'start'}, default 'end'
            Start or end of the timespan

        Returns
        -------
        resampled : Period
        """
        freq = self._maybe_convert_freq(freq)
        how = _validate_end_alias(how)
        base1, mult1 = get_freq_code(self.freq)
        base2, mult2 = get_freq_code(freq)

        # mult1 can't be negative or 0
        end = how == 'E'
        if end:
            ordinal = self.ordinal + mult1 - 1
        else:
            ordinal = self.ordinal
        ordinal = period_asfreq(ordinal, base1, base2, end)

        return Period(ordinal=ordinal, freq=freq)

    @property
    def start_time(self):
        """
        Get the Timestamp for the start of the period.

        Returns
        -------
        Timestamp

        See also
        --------
        Period.end_time : Return the end Timestamp.
        Period.dayofyear : Return the day of year.
        Period.daysinmonth : Return the days in that month.
        Period.dayofweek : Return the day of the week.

        Examples
        --------
        >>> period = pd.Period('2012-1-1', freq='D')
        >>> period
        Period('2012-01-01', 'D')

        >>> period.start_time
        Timestamp('2012-01-01 00:00:00')

        >>> period.end_time
        Timestamp('2012-01-01 23:59:59.999999999')
        """
        return self.to_timestamp(how='S')

    @property
    def end_time(self):
        # freq.n can't be negative or 0
        # ordinal = (self + self.freq.n).start_time.value - 1
        ordinal = (self + 1).start_time.value - 1
        return Timestamp(ordinal)

    def to_timestamp(self, freq=None, how='start', tz=None):
        """
        Return the Timestamp representation of the Period at the target
        frequency at the specified end (how) of the Period

        Parameters
        ----------
        freq : string or DateOffset
            Target frequency. Default is 'D' if self.freq is week or
            longer and 'S' otherwise
        how: str, default 'S' (start)
            'S', 'E'. Can be aliased as case insensitive
            'Start', 'Finish', 'Begin', 'End'

        Returns
        -------
        Timestamp
        """
        if freq is not None:
            freq = self._maybe_convert_freq(freq)
        how = _validate_end_alias(how)

        end = how == 'E'
        if end:
            return (self + 1).to_timestamp(how='start') - Timedelta(1, 'ns')

        if freq is None:
            base, mult = get_freq_code(self.freq)
            freq = get_to_timestamp_base(base)

        base, mult = get_freq_code(freq)
        val = self.asfreq(freq, how)

        dt64 = period_ordinal_to_dt64(val.ordinal, base)
        return Timestamp(dt64, tz=tz)

    @property
    def year(self):
        base, mult = get_freq_code(self.freq)
        return pyear(self.ordinal, base)

    @property
    def month(self):
        base, mult = get_freq_code(self.freq)
        return pmonth(self.ordinal, base)

    @property
    def day(self):
        """
        Get day of the month that a Period falls on.

        Returns
        -------
        int

        See Also
        --------
        Period.dayofweek : Get the day of the week

        Period.dayofyear : Get the day of the year

        Examples
        --------
        >>> p = pd.Period("2018-03-11", freq='H')
        >>> p.day
        11
        """
        base, mult = get_freq_code(self.freq)
        return pday(self.ordinal, base)

    @property
    def hour(self):
        """
        Get the hour of the day component of the Period.

        Returns
        -------
        int
            The hour as an integer, between 0 and 23.

        See Also
        --------
        Period.second : Get the second component of the Period.
        Period.minute : Get the minute component of the Period.

        Examples
        --------
        >>> p = pd.Period("2018-03-11 13:03:12.050000")
        >>> p.hour
        13

        Period longer than a day

        >>> p = pd.Period("2018-03-11", freq="M")
        >>> p.hour
        0
        """
        base, mult = get_freq_code(self.freq)
        return phour(self.ordinal, base)

    @property
    def minute(self):
        """
        Get minute of the hour component of the Period.

        Returns
        -------
        int
            The minute as an integer, between 0 and 59.

        See Also
        --------
        Period.hour : Get the hour component of the Period.
        Period.second : Get the second component of the Period.

        Examples
        --------
        >>> p = pd.Period("2018-03-11 13:03:12.050000")
        >>> p.minute
        3
        """
        base, mult = get_freq_code(self.freq)
        return pminute(self.ordinal, base)

    @property
    def second(self):
        """
        Get the second component of the Period.

        Returns
        -------
        int
            The second of the Period (ranges from 0 to 59).

        See Also
        --------
        Period.hour : Get the hour component of the Period.
        Period.minute : Get the minute component of the Period.

        Examples
        --------
        >>> p = pd.Period("2018-03-11 13:03:12.050000")
        >>> p.second
        12
        """
        base, mult = get_freq_code(self.freq)
        return psecond(self.ordinal, base)

    @property
    def weekofyear(self):
        base, mult = get_freq_code(self.freq)
        return pweek(self.ordinal, base)

    @property
    def week(self):
        """
        Get the week of the year on the given Period.

        Returns
        -------
        int

        See Also
        --------
        Period.dayofweek : Get the day component of the Period.
        Period.weekday : Get the day component of the Period.

        Examples
        --------
        >>> p = pd.Period("2018-03-11", "H")
        >>> p.week
        10

        >>> p = pd.Period("2018-02-01", "D")
        >>> p.week
        5

        >>> p = pd.Period("2018-01-06", "D")
        >>> p.week
        1
        """
        return self.weekofyear

    @property
    def dayofweek(self):
        """
        Day of the week the period lies in, with Monday=0 and Sunday=6.

        If the period frequency is lower than daily (e.g. hourly), and the
        period spans over multiple days, the day at the start of the period is
        used.

        If the frequency is higher than daily (e.g. monthly), the last day
        of the period is used.

        Returns
        -------
        int
            Day of the week.

        See Also
        --------
        Period.dayofweek : Day of the week the period lies in.
        Period.weekday : Alias of Period.dayofweek.
        Period.day : Day of the month.
        Period.dayofyear : Day of the year.

        Examples
        --------
        >>> per = pd.Period('2017-12-31 22:00', 'H')
        >>> per.dayofweek
        6

        For periods that span over multiple days, the day at the beginning of
        the period is returned.

        >>> per = pd.Period('2017-12-31 22:00', '4H')
        >>> per.dayofweek
        6
        >>> per.start_time.dayofweek
        6

        For periods with a frequency higher than days, the last day of the
        period is returned.

        >>> per = pd.Period('2018-01', 'M')
        >>> per.dayofweek
        2
        >>> per.end_time.dayofweek
        2
        """
        base, mult = get_freq_code(self.freq)
        return pweekday(self.ordinal, base)

    @property
    def weekday(self):
        """
        Day of the week the period lies in, with Monday=0 and Sunday=6.

        If the period frequency is lower than daily (e.g. hourly), and the
        period spans over multiple days, the day at the start of the period is
        used.

        If the frequency is higher than daily (e.g. monthly), the last day
        of the period is used.

        Returns
        -------
        int
            Day of the week.

        See Also
        --------
        Period.dayofweek : Day of the week the period lies in.
        Period.weekday : Alias of Period.dayofweek.
        Period.day : Day of the month.
        Period.dayofyear : Day of the year.

        Examples
        --------
        >>> per = pd.Period('2017-12-31 22:00', 'H')
        >>> per.dayofweek
        6

        For periods that span over multiple days, the day at the beginning of
        the period is returned.

        >>> per = pd.Period('2017-12-31 22:00', '4H')
        >>> per.dayofweek
        6
        >>> per.start_time.dayofweek
        6

        For periods with a frequency higher than days, the last day of the
        period is returned.

        >>> per = pd.Period('2018-01', 'M')
        >>> per.dayofweek
        2
        >>> per.end_time.dayofweek
        2
        """
        # Docstring is a duplicate from dayofweek. Reusing docstrings with
        # Appender doesn't work for properties in Cython files, and setting
        # the __doc__ attribute is also not possible.
        return self.dayofweek

    @property
    def dayofyear(self):
        """
        Return the day of the year.

        This attribute returns the day of the year on which the particular
        date occurs. The return value ranges between 1 to 365 for regular
        years and 1 to 366 for leap years.

        Returns
        -------
        int
            The day of year.

        See Also
        --------
        Period.day : Return the day of the month.
        Period.dayofweek : Return the day of week.
        PeriodIndex.dayofyear : Return the day of year of all indexes.

        Examples
        --------
        >>> period = pd.Period("2015-10-23", freq='H')
        >>> period.dayofyear
        296
        >>> period = pd.Period("2012-12-31", freq='D')
        >>> period.dayofyear
        366
        >>> period = pd.Period("2013-01-01", freq='D')
        >>> period.dayofyear
        1
        """
        base, mult = get_freq_code(self.freq)
        return pday_of_year(self.ordinal, base)

    @property
    def quarter(self):
        base, mult = get_freq_code(self.freq)
        return pquarter(self.ordinal, base)

    @property
    def qyear(self):
        """
        Fiscal year the Period lies in according to its starting-quarter.

        The `year` and the `qyear` of the period will be the same if the fiscal
        and calendar years are the same. When they are not, the fiscal year
        can be different from the calendar year of the period.

        Returns
        -------
        int
            The fiscal year of the period.

        See Also
        --------
        Period.year : Return the calendar year of the period.

        Examples
        --------
        If the natural and fiscal year are the same, `qyear` and `year` will
        be the same.

        >>> per = pd.Period('2018Q1', freq='Q')
        >>> per.qyear
        2018
        >>> per.year
        2018

        If the fiscal year starts in April (`Q-MAR`), the first quarter of
        2018 will start in April 2017. `year` will then be 2018, but `qyear`
        will be the fiscal year, 2018.

        >>> per = pd.Period('2018Q1', freq='Q-MAR')
        >>> per.start_time
        Timestamp('2017-04-01 00:00:00')
        >>> per.qyear
        2018
        >>> per.year
        2017
        """
        base, mult = get_freq_code(self.freq)
        return pqyear(self.ordinal, base)

    @property
    def days_in_month(self):
        """
        Get the total number of days in the month that this period falls on.

        Returns
        -------
        int

        See Also
        --------
        Period.daysinmonth : Gets the number of days in the month.
        DatetimeIndex.daysinmonth : Gets the number of days in the month.
        calendar.monthrange : Returns a tuple containing weekday
            (0-6 ~ Mon-Sun) and number of days (28-31).

        Examples
        --------
        >>> p = pd.Period('2018-2-17')
        >>> p.days_in_month
        28

        >>> pd.Period('2018-03-01').days_in_month
        31

        Handles the leap year case as well:

        >>> p = pd.Period('2016-2-17')
        >>> p.days_in_month
        29
        """
        base, mult = get_freq_code(self.freq)
        return pdays_in_month(self.ordinal, base)

    @property
    def daysinmonth(self):
        """
        Get the total number of days of the month that the Period falls in.

        Returns
        -------
        int

        See Also
        --------
        Period.days_in_month : Return the days of the month
        Period.dayofyear : Return the day of the year

        Examples
        --------
        >>> p = pd.Period("2018-03-11", freq='H')
        >>> p.daysinmonth
        31
        """
        return self.days_in_month

    @property
    def is_leap_year(self):
        return bool(is_leapyear(self.year))

    @classmethod
    def now(cls, freq=None):
        return Period(datetime.now(), freq=freq)

    # HACK IT UP AND YOU BETTER FIX IT SOON
    def __str__(self):
        return self.__unicode__()

    @property
    def freqstr(self):
        return self.freq.freqstr

    def __repr__(self):
        base, mult = get_freq_code(self.freq)
        formatted = period_format(self.ordinal, base)
        return "Period('%s', '%s')" % (formatted, self.freqstr)

    def __unicode__(self):
        """
        Return a string representation for a particular DataFrame

        Invoked by unicode(df) in py2 only. Yields a Unicode String in both
        py2/py3.
        """
        base, mult = get_freq_code(self.freq)
        formatted = period_format(self.ordinal, base)
        value = ("%s" % formatted)
        return value

    def __setstate__(self, state):
        self.freq = state[1]
        self.ordinal = state[2]

    def __reduce__(self):
        object_state = None, self.freq, self.ordinal
        return (Period, object_state)

    def strftime(self, fmt):
        """
        Returns the string representation of the :class:`Period`, depending
        on the selected ``fmt``. ``fmt`` must be a string
        containing one or several directives.  The method recognizes the same
        directives as the :func:`time.strftime` function of the standard Python
        distribution, as well as the specific additional directives ``%f``,
        ``%F``, ``%q``. (formatting & docs originally from scikits.timeries)

        +-----------+--------------------------------+-------+
        | Directive | Meaning                        | Notes |
        +===========+================================+=======+
        | ``%a``    | Locale's abbreviated weekday   |       |
        |           | name.                          |       |
        +-----------+--------------------------------+-------+
        | ``%A``    | Locale's full weekday name.    |       |
        +-----------+--------------------------------+-------+
        | ``%b``    | Locale's abbreviated month     |       |
        |           | name.                          |       |
        +-----------+--------------------------------+-------+
        | ``%B``    | Locale's full month name.      |       |
        +-----------+--------------------------------+-------+
        | ``%c``    | Locale's appropriate date and  |       |
        |           | time representation.           |       |
        +-----------+--------------------------------+-------+
        | ``%d``    | Day of the month as a decimal  |       |
        |           | number [01,31].                |       |
        +-----------+--------------------------------+-------+
        | ``%f``    | 'Fiscal' year without a        | \(1)  |
        |           | century  as a decimal number   |       |
        |           | [00,99]                        |       |
        +-----------+--------------------------------+-------+
        | ``%F``    | 'Fiscal' year with a century   | \(2)  |
        |           | as a decimal number            |       |
        +-----------+--------------------------------+-------+
        | ``%H``    | Hour (24-hour clock) as a      |       |
        |           | decimal number [00,23].        |       |
        +-----------+--------------------------------+-------+
        | ``%I``    | Hour (12-hour clock) as a      |       |
        |           | decimal number [01,12].        |       |
        +-----------+--------------------------------+-------+
        | ``%j``    | Day of the year as a decimal   |       |
        |           | number [001,366].              |       |
        +-----------+--------------------------------+-------+
        | ``%m``    | Month as a decimal number      |       |
        |           | [01,12].                       |       |
        +-----------+--------------------------------+-------+
        | ``%M``    | Minute as a decimal number     |       |
        |           | [00,59].                       |       |
        +-----------+--------------------------------+-------+
        | ``%p``    | Locale's equivalent of either  | \(3)  |
        |           | AM or PM.                      |       |
        +-----------+--------------------------------+-------+
        | ``%q``    | Quarter as a decimal number    |       |
        |           | [01,04]                        |       |
        +-----------+--------------------------------+-------+
        | ``%S``    | Second as a decimal number     | \(4)  |
        |           | [00,61].                       |       |
        +-----------+--------------------------------+-------+
        | ``%U``    | Week number of the year        | \(5)  |
        |           | (Sunday as the first day of    |       |
        |           | the week) as a decimal number  |       |
        |           | [00,53].  All days in a new    |       |
        |           | year preceding the first       |       |
        |           | Sunday are considered to be in |       |
        |           | week 0.                        |       |
        +-----------+--------------------------------+-------+
        | ``%w``    | Weekday as a decimal number    |       |
        |           | [0(Sunday),6].                 |       |
        +-----------+--------------------------------+-------+
        | ``%W``    | Week number of the year        | \(5)  |
        |           | (Monday as the first day of    |       |
        |           | the week) as a decimal number  |       |
        |           | [00,53].  All days in a new    |       |
        |           | year preceding the first       |       |
        |           | Monday are considered to be in |       |
        |           | week 0.                        |       |
        +-----------+--------------------------------+-------+
        | ``%x``    | Locale's appropriate date      |       |
        |           | representation.                |       |
        +-----------+--------------------------------+-------+
        | ``%X``    | Locale's appropriate time      |       |
        |           | representation.                |       |
        +-----------+--------------------------------+-------+
        | ``%y``    | Year without century as a      |       |
        |           | decimal number [00,99].        |       |
        +-----------+--------------------------------+-------+
        | ``%Y``    | Year with century as a decimal |       |
        |           | number.                        |       |
        +-----------+--------------------------------+-------+
        | ``%Z``    | Time zone name (no characters  |       |
        |           | if no time zone exists).       |       |
        +-----------+--------------------------------+-------+
        | ``%%``    | A literal ``'%'`` character.   |       |
        +-----------+--------------------------------+-------+

        Notes
        -----

        (1)
            The ``%f`` directive is the same as ``%y`` if the frequency is
            not quarterly.
            Otherwise, it corresponds to the 'fiscal' year, as defined by
            the :attr:`qyear` attribute.

        (2)
            The ``%F`` directive is the same as ``%Y`` if the frequency is
            not quarterly.
            Otherwise, it corresponds to the 'fiscal' year, as defined by
            the :attr:`qyear` attribute.

        (3)
            The ``%p`` directive only affects the output hour field
            if the ``%I`` directive is used to parse the hour.

        (4)
            The range really is ``0`` to ``61``; this accounts for leap
            seconds and the (very rare) double leap seconds.

        (5)
            The ``%U`` and ``%W`` directives are only used in calculations
            when the day of the week and the year are specified.

        Examples
        --------

        >>> a = Period(freq='Q-JUL', year=2006, quarter=1)
        >>> a.strftime('%F-Q%q')
        '2006-Q1'
        >>> # Output the last month in the quarter of this date
        >>> a.strftime('%b-%Y')
        'Oct-2005'
        >>>
        >>> a = Period(freq='D', year=2001, month=1, day=1)
        >>> a.strftime('%d-%b-%Y')
        '01-Jan-2006'
        >>> a.strftime('%b. %d, %Y was a %A')
        'Jan. 01, 2001 was a Monday'
        """
        base, mult = get_freq_code(self.freq)
        return period_format(self.ordinal, base, fmt)


class Period(_Period):
    """
    Represents a period of time

    Parameters
    ----------
    value : Period or compat.string_types, default None
        The time period represented (e.g., '4Q2005')
    freq : str, default None
        One of pandas period strings or corresponding objects
    year : int, default None
    month : int, default 1
    quarter : int, default None
    day : int, default 1
    hour : int, default 0
    minute : int, default 0
    second : int, default 0
    """

    def __new__(cls, value=None, freq=None, ordinal=None,
                year=None, month=None, quarter=None, day=None,
                hour=None, minute=None, second=None):
        # freq points to a tuple (base, mult);  base is one of the defined
        # periods such as A, Q, etc. Every five minutes would be, e.g.,
        # ('T', 5) but may be passed in as a string like '5T'

        # ordinal is the period offset from the gregorian proleptic epoch

        cdef _Period self

        if freq is not None:
            freq = cls._maybe_convert_freq(freq)

        if ordinal is not None and value is not None:
            raise ValueError(("Only value or ordinal but not both should be "
                              "given but not both"))
        elif ordinal is not None:
            if not util.is_integer_object(ordinal):
                raise ValueError("Ordinal must be an integer")
            if freq is None:
                raise ValueError('Must supply freq for ordinal value')

        elif value is None:
            if (year is None and month is None and
                    quarter is None and day is None and
                    hour is None and minute is None and second is None):
                ordinal = iNaT
            else:
                if freq is None:
                    raise ValueError("If value is None, freq cannot be None")

                # set defaults
                month = 1 if month is None else month
                day = 1 if day is None else day
                hour = 0 if hour is None else hour
                minute = 0 if minute is None else minute
                second = 0 if second is None else second

                ordinal = _ordinal_from_fields(year, month, quarter, day,
                                               hour, minute, second, freq)

        elif is_period_object(value):
            other = value
            if freq is None or get_freq_code(
                    freq) == get_freq_code(other.freq):
                ordinal = other.ordinal
                freq = other.freq
            else:
                converted = other.asfreq(freq)
                ordinal = converted.ordinal

        elif is_null_datetimelike(value) or value in nat_strings:
            ordinal = iNaT

        elif is_string_object(value) or util.is_integer_object(value):
            if util.is_integer_object(value):
                value = str(value)
            value = value.upper()
            dt, _, reso = parse_time_string(value, freq)
            if dt is NaT:
                ordinal = iNaT

            if freq is None:
                try:
                    freq = Resolution.get_freq(reso)
                except KeyError:
                    raise ValueError(
                        "Invalid frequency or could not infer: %s" % reso)

        elif isinstance(value, datetime):
            dt = value
            if freq is None:
                raise ValueError('Must supply freq for datetime value')
        elif util.is_datetime64_object(value):
            dt = Timestamp(value)
            if freq is None:
                raise ValueError('Must supply freq for datetime value')
        elif isinstance(value, date):
            dt = datetime(year=value.year, month=value.month, day=value.day)
            if freq is None:
                raise ValueError('Must supply freq for datetime value')
        else:
            msg = "Value must be Period, string, integer, or datetime"
            raise ValueError(msg)

        if ordinal is None:
            base, mult = get_freq_code(freq)
            ordinal = period_ordinal(dt.year, dt.month, dt.day,
                                     dt.hour, dt.minute, dt.second,
                                     dt.microsecond, 0, base)

        return cls._from_ordinal(ordinal, freq)


cdef int64_t _ordinal_from_fields(int year, int month, quarter, int day,
                                  int hour, int minute, int second, freq):
    base, mult = get_freq_code(freq)
    if quarter is not None:
        year, month = quarter_to_myear(year, quarter, freq)

    return period_ordinal(year, month, day, hour,
                          minute, second, 0, 0, base)


def quarter_to_myear(int year, int quarter, freq):
    """
    A quarterly frequency defines a "year" which may not coincide with
    the calendar-year.  Find the calendar-year and calendar-month associated
    with the given year and quarter under the `freq`-derived calendar.

    Parameters
    ----------
    year : int
    quarter : int
    freq : DateOffset

    Returns
    -------
    year : int
    month : int

    See Also
    --------
    Period.qyear
    """
    if quarter <= 0 or quarter > 4:
        raise ValueError('Quarter must be 1 <= q <= 4')

    mnum = MONTH_NUMBERS[get_rule_month(freq)] + 1
    month = (mnum + (quarter - 1) * 3) % 12 + 1
    if month > mnum:
        year -= 1

    return year, month


def _validate_end_alias(how):
    how_dict = {'S': 'S', 'E': 'E',
                'START': 'S', 'FINISH': 'E',
                'BEGIN': 'S', 'END': 'E'}
    how = how_dict.get(str(how).upper())
    if how not in {'S', 'E'}:
        raise ValueError('How must be one of S or E')
    return how
