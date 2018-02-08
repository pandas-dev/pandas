/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.

Borrowed and derived code from scikits.timeseries that we will expose via
Cython to pandas. This primarily concerns interval representation and
frequency conversion routines.

See end of file for stuff pandas uses (search for 'pandas').
*/

#include "period_helper.h"
#include "../datetime/np_datetime.h"

/* ------------------------------------------------------------------
 * Code derived from scikits.timeseries
 * ------------------------------------------------------------------*/

static int mod_compat(int x, int m) {
    int result = x % m;
    if (result < 0) return result + m;
    return result;
}

static int floordiv(int x, int divisor) {
    if (x < 0) {
        if (mod_compat(x, divisor)) {
            return x / divisor - 1;
        } else {
            return x / divisor;
        }
    } else {
        return x / divisor;
    }
}


static int monthToQuarter(int month) { return ((month - 1) / 3) + 1; }


/* Find the absdate (days elapsed since datetime(1, 1, 1)
 * for the given year/month/day.
 * Assumes GREGORIAN_CALENDAR */
static npy_int64 dInfoCalc_SetFromDateAndTime(int year, int month, int day) {
    /* Calculate the absolute date */
    pandas_datetimestruct dts;
    npy_int64 unix_date;

    memset(&dts, 0, sizeof(pandas_datetimestruct));
    dts.year = year;
    dts.month = month;
    dts.day = day;
    unix_date = pandas_datetimestruct_to_datetime(PANDAS_FR_D, &dts);
    return ORD_OFFSET + unix_date;
}

/* Sets the date part of the date_info struct
   Assumes GREGORIAN_CALENDAR */
static int dInfoCalc_SetFromAbsDate(register struct date_info *dinfo,
                                    npy_int64 absdate) {
    pandas_datetimestruct dts;

    pandas_datetime_to_datetimestruct(absdate - ORD_OFFSET, PANDAS_FR_D, &dts);
    dinfo->year = dts.year;
    dinfo->month = dts.month;
    dinfo->day = dts.day;

    dinfo->absdate = absdate;
    return 0;
}

///////////////////////////////////////////////

// frequency specific conversion routines
// each function must take an integer fromDate and
// a char relation ('S' or 'E' for 'START' or 'END')
///////////////////////////////////////////////////////////////////////

// helpers for frequency conversion routines //

static npy_int64 daytime_conversion_factor_matrix[7][7] = {
    {1, 24, 1440, 86400, 86400000, 86400000000, 86400000000000},
    {0,  1,   60,  3600,  3600000,  3600000000,  3600000000000},
    {0,  0,   1,     60,    60000,    60000000,    60000000000},
    {0,  0,   0,      1,     1000,     1000000,     1000000000},
    {0,  0,   0,      0,        1,        1000,        1000000},
    {0,  0,   0,      0,        0,           1,           1000},
    {0,  0,   0,      0,        0,           0,              1}};

PANDAS_INLINE int max_value(int a, int b) { return a > b ? a : b; }

PANDAS_INLINE int min_value(int a, int b) { return a < b ? a : b; }

PANDAS_INLINE int get_freq_group(int freq) { return (freq / 1000) * 1000; }

PANDAS_INLINE int get_freq_group_index(int freq) { return freq / 1000; }


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

PANDAS_INLINE npy_int64 upsample_daytime(npy_int64 ordinal,
                                         asfreq_info *af_info) {
    if (af_info->is_end) {
        return (ordinal + 1) * af_info->intraday_conversion_factor - 1;
    } else {
        return ordinal * af_info->intraday_conversion_factor;
    }
}

PANDAS_INLINE npy_int64 downsample_daytime(npy_int64 ordinal,
                                           asfreq_info *af_info) {
    return ordinal / (af_info->intraday_conversion_factor);
}

PANDAS_INLINE npy_int64 transform_via_day(npy_int64 ordinal,
                                          asfreq_info *af_info,
                                          freq_conv_func first_func,
                                          freq_conv_func second_func) {
    npy_int64 result;

    result = (*first_func)(ordinal, af_info);
    result = (*second_func)(result, af_info);

    return result;
}

static npy_int64 DtoB_weekday(npy_int64 absdate) {
    return (((absdate) / 7) * 5) + (absdate) % 7 - BDAY_OFFSET;
}

static npy_int64 DtoB(struct date_info *dinfo, int roll_back) {
    int day_of_week = dayofweek(dinfo->year, dinfo->month, dinfo->day);
    npy_int64 absdate = dinfo->absdate;

    if (roll_back == 1) {
        if (day_of_week > 4) {
            // change to friday before weekend
            absdate -= (day_of_week - 4);
        }
    } else {
        if (day_of_week > 4) {
            // change to Monday after weekend
            absdate += (7 - day_of_week);
        }
    }
    return DtoB_weekday(absdate);
}

static npy_int64 absdate_from_ymd(int y, int m, int d) {
    return dInfoCalc_SetFromDateAndTime(y, m, d);
}

//************ FROM DAILY ***************

static npy_int64 asfreq_DTtoA(npy_int64 ordinal, asfreq_info *af_info) {
    struct date_info dinfo;
    ordinal = downsample_daytime(ordinal, af_info);
    dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET);
    if (dinfo.month > af_info->to_a_year_end) {
        return (npy_int64)(dinfo.year + 1 - BASE_YEAR);
    } else {
        return (npy_int64)(dinfo.year - BASE_YEAR);
    }
}

static npy_int64 DtoQ_yq(npy_int64 ordinal, asfreq_info *af_info, int *year,
                         int *quarter) {
    struct date_info dinfo;
    dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET);
    if (af_info->to_q_year_end != 12) {
        dinfo.month -= af_info->to_q_year_end;
        if (dinfo.month <= 0) {
            dinfo.month += 12;
        } else {
            dinfo.year += 1;
        }
    }

    *year = dinfo.year;
    *quarter = monthToQuarter(dinfo.month);

    return 0;
}

static npy_int64 asfreq_DTtoQ(npy_int64 ordinal, asfreq_info *af_info) {
    int year, quarter;

    ordinal = downsample_daytime(ordinal, af_info);

    DtoQ_yq(ordinal, af_info, &year, &quarter);
    return (npy_int64)((year - BASE_YEAR) * 4 + quarter - 1);
}

static npy_int64 asfreq_DTtoM(npy_int64 ordinal, asfreq_info *af_info) {
    struct date_info dinfo;

    ordinal = downsample_daytime(ordinal, af_info);

    dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET);
    return (npy_int64)((dinfo.year - BASE_YEAR) * 12 + dinfo.month - 1);
}

static npy_int64 asfreq_DTtoW(npy_int64 ordinal, asfreq_info *af_info) {
    ordinal = downsample_daytime(ordinal, af_info);
    return (ordinal + ORD_OFFSET - (1 + af_info->to_week_end)) / 7 + 1 -
           WEEK_OFFSET;
}

static npy_int64 asfreq_DTtoB(npy_int64 ordinal, asfreq_info *af_info) {
    struct date_info dinfo;
    int roll_back;

    ordinal = downsample_daytime(ordinal, af_info);

    dInfoCalc_SetFromAbsDate(&dinfo, ordinal + ORD_OFFSET);

    // This usage defines roll_back the opposite way from the others
    roll_back = 1 - af_info->is_end;
    return DtoB(&dinfo, roll_back);
}

// all intra day calculations are now done within one function
static npy_int64 asfreq_DownsampleWithinDay(npy_int64 ordinal,
                                            asfreq_info *af_info) {
    return downsample_daytime(ordinal, af_info);
}

static npy_int64 asfreq_UpsampleWithinDay(npy_int64 ordinal,
                                          asfreq_info *af_info) {
    return upsample_daytime(ordinal, af_info);
}
//************ FROM BUSINESS ***************

static npy_int64 asfreq_BtoDT(npy_int64 ordinal, asfreq_info *af_info) {
    ordinal += BDAY_OFFSET;
    ordinal =
        (((ordinal - 1) / 5) * 7 + mod_compat(ordinal - 1, 5) + 1 - ORD_OFFSET);

    return upsample_daytime(ordinal, af_info);
}

static npy_int64 asfreq_BtoA(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_BtoDT, asfreq_DTtoA);
}

static npy_int64 asfreq_BtoQ(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_BtoDT, asfreq_DTtoQ);
}

static npy_int64 asfreq_BtoM(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_BtoDT, asfreq_DTtoM);
}

static npy_int64 asfreq_BtoW(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_BtoDT, asfreq_DTtoW);
}

//************ FROM WEEKLY ***************

static npy_int64 asfreq_WtoDT(npy_int64 ordinal, asfreq_info *af_info) {
    ordinal = (ordinal + WEEK_OFFSET) * 7 +
               af_info->from_week_end - ORD_OFFSET +
               (7 - 1) * (af_info->is_end - 1);
    return upsample_daytime(ordinal, af_info);
}

static npy_int64 asfreq_WtoA(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_WtoDT, asfreq_DTtoA);
}

static npy_int64 asfreq_WtoQ(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_WtoDT, asfreq_DTtoQ);
}

static npy_int64 asfreq_WtoM(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_WtoDT, asfreq_DTtoM);
}

static npy_int64 asfreq_WtoW(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_WtoDT, asfreq_DTtoW);
}

static npy_int64 asfreq_WtoB(npy_int64 ordinal, asfreq_info *af_info) {
    struct date_info dinfo;
    int roll_back = af_info->is_end;
    dInfoCalc_SetFromAbsDate(
            &dinfo, asfreq_WtoDT(ordinal, af_info) + ORD_OFFSET);

    return DtoB(&dinfo, roll_back);
}

//************ FROM MONTHLY ***************
static void MtoD_ym(npy_int64 ordinal, int *y, int *m) {
    *y = floordiv(ordinal, 12) + BASE_YEAR;
    *m = mod_compat(ordinal, 12) + 1;
}

static npy_int64 asfreq_MtoDT(npy_int64 ordinal, asfreq_info *af_info) {
    npy_int64 absdate;
    int y, m;

    ordinal += af_info->is_end;
    MtoD_ym(ordinal, &y, &m);
    absdate = absdate_from_ymd(y, m, 1);
    ordinal = absdate - ORD_OFFSET;

    ordinal -= af_info->is_end;
    return upsample_daytime(ordinal, af_info);
}

static npy_int64 asfreq_MtoA(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_MtoDT, asfreq_DTtoA);
}

static npy_int64 asfreq_MtoQ(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_MtoDT, asfreq_DTtoQ);
}

static npy_int64 asfreq_MtoW(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_MtoDT, asfreq_DTtoW);
}

static npy_int64 asfreq_MtoB(npy_int64 ordinal, asfreq_info *af_info) {
    struct date_info dinfo;
    int roll_back = af_info->is_end;

    dInfoCalc_SetFromAbsDate(
            &dinfo, asfreq_MtoDT(ordinal, af_info) + ORD_OFFSET);

    return DtoB(&dinfo, roll_back);
}

//************ FROM QUARTERLY ***************

static void QtoD_ym(npy_int64 ordinal, int *y, int *m, asfreq_info *af_info) {
    *y = floordiv(ordinal, 4) + BASE_YEAR;
    *m = mod_compat(ordinal, 4) * 3 + 1;

    if (af_info->from_q_year_end != 12) {
        *m += af_info->from_q_year_end;
        if (*m > 12) {
            *m -= 12;
        } else {
            *y -= 1;
        }
    }
}

static npy_int64 asfreq_QtoDT(npy_int64 ordinal, asfreq_info *af_info) {
    npy_int64 absdate;
    int y, m;

    ordinal += af_info->is_end;
    QtoD_ym(ordinal, &y, &m, af_info);

    absdate = absdate_from_ymd(y, m, 1);

    absdate -= af_info->is_end;
    return upsample_daytime(absdate - ORD_OFFSET, af_info);
}

static npy_int64 asfreq_QtoQ(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_QtoDT, asfreq_DTtoQ);
}

static npy_int64 asfreq_QtoA(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_QtoDT, asfreq_DTtoA);
}

static npy_int64 asfreq_QtoM(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_QtoDT, asfreq_DTtoM);
}

static npy_int64 asfreq_QtoW(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_QtoDT, asfreq_DTtoW);
}

static npy_int64 asfreq_QtoB(npy_int64 ordinal, asfreq_info *af_info) {
    struct date_info dinfo;
    int roll_back = af_info->is_end;

    dInfoCalc_SetFromAbsDate(
            &dinfo, asfreq_QtoDT(ordinal, af_info) + ORD_OFFSET);

    return DtoB(&dinfo, roll_back);
}

//************ FROM ANNUAL ***************

static npy_int64 asfreq_AtoDT(npy_int64 ordinal, asfreq_info *af_info) {
    npy_int64 absdate;

    // start from 1970
    npy_int64 year = ordinal + BASE_YEAR;

    int month = (af_info->from_a_year_end % 12) + 1;
    if (af_info->from_a_year_end != 12) {
        year -= 1;
    }

    year += af_info->is_end;
    absdate = absdate_from_ymd(year, month, 1);

    absdate -= af_info->is_end;
    return upsample_daytime(absdate - ORD_OFFSET, af_info);
}

static npy_int64 asfreq_AtoA(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_AtoDT, asfreq_DTtoA);
}

static npy_int64 asfreq_AtoQ(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_AtoDT, asfreq_DTtoQ);
}

static npy_int64 asfreq_AtoM(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_AtoDT, asfreq_DTtoM);
}

static npy_int64 asfreq_AtoW(npy_int64 ordinal, asfreq_info *af_info) {
    return transform_via_day(ordinal, af_info, asfreq_AtoDT, asfreq_DTtoW);
}

static npy_int64 asfreq_AtoB(npy_int64 ordinal, asfreq_info *af_info) {
    struct date_info dinfo;
    int roll_back = af_info->is_end;
    dInfoCalc_SetFromAbsDate(
            &dinfo, asfreq_AtoDT(ordinal, af_info) + ORD_OFFSET);

    return DtoB(&dinfo, roll_back);
}

static npy_int64 nofunc(npy_int64 ordinal, asfreq_info *af_info) {
    return INT_ERR_CODE;
}
static npy_int64 no_op(npy_int64 ordinal, asfreq_info *af_info) {
    return ordinal;
}

// end of frequency specific conversion routines

static int calc_a_year_end(int freq, int group) {
    int result = (freq - group) % 12;
    if (result == 0) {
        return 12;
    } else {
        return result;
    }
}

static int calc_week_end(int freq, int group) { return freq - group; }

void get_asfreq_info(int fromFreq, int toFreq, char relation,
                     asfreq_info *af_info) {
    int fromGroup = get_freq_group(fromFreq);
    int toGroup = get_freq_group(toFreq);

    if (relation == 'E') {
        af_info->is_end = 1;
    } else {
        af_info->is_end = 0;
    }

    af_info->intraday_conversion_factor = get_daytime_conversion_factor(
        get_freq_group_index(max_value(fromGroup, FR_DAY)),
        get_freq_group_index(max_value(toGroup, FR_DAY)));

    switch (fromGroup) {
        case FR_WK:
            af_info->from_week_end = calc_week_end(fromFreq, fromGroup);
            break;
        case FR_ANN:
            af_info->from_a_year_end = calc_a_year_end(fromFreq, fromGroup);
            break;
        case FR_QTR:
            af_info->from_q_year_end = calc_a_year_end(fromFreq, fromGroup);
            break;
    }

    switch (toGroup) {
        case FR_WK:
            af_info->to_week_end = calc_week_end(toFreq, toGroup);
            break;
        case FR_ANN:
            af_info->to_a_year_end = calc_a_year_end(toFreq, toGroup);
            break;
        case FR_QTR:
            af_info->to_q_year_end = calc_a_year_end(toFreq, toGroup);
            break;
    }
}

freq_conv_func get_asfreq_func(int fromFreq, int toFreq) {
    int fromGroup = get_freq_group(fromFreq);
    int toGroup = get_freq_group(toFreq);

    if (fromGroup == FR_UND) {
        fromGroup = FR_DAY;
    }

    switch (fromGroup) {
        case FR_ANN:
            switch (toGroup) {
                case FR_ANN:
                    return &asfreq_AtoA;
                case FR_QTR:
                    return &asfreq_AtoQ;
                case FR_MTH:
                    return &asfreq_AtoM;
                case FR_WK:
                    return &asfreq_AtoW;
                case FR_BUS:
                    return &asfreq_AtoB;
                case FR_DAY:
                case FR_HR:
                case FR_MIN:
                case FR_SEC:
                case FR_MS:
                case FR_US:
                case FR_NS:
                    return &asfreq_AtoDT;

                default:
                    return &nofunc;
            }

        case FR_QTR:
            switch (toGroup) {
                case FR_ANN:
                    return &asfreq_QtoA;
                case FR_QTR:
                    return &asfreq_QtoQ;
                case FR_MTH:
                    return &asfreq_QtoM;
                case FR_WK:
                    return &asfreq_QtoW;
                case FR_BUS:
                    return &asfreq_QtoB;
                case FR_DAY:
                case FR_HR:
                case FR_MIN:
                case FR_SEC:
                case FR_MS:
                case FR_US:
                case FR_NS:
                    return &asfreq_QtoDT;
                default:
                    return &nofunc;
            }

        case FR_MTH:
            switch (toGroup) {
                case FR_ANN:
                    return &asfreq_MtoA;
                case FR_QTR:
                    return &asfreq_MtoQ;
                case FR_MTH:
                    return &no_op;
                case FR_WK:
                    return &asfreq_MtoW;
                case FR_BUS:
                    return &asfreq_MtoB;
                case FR_DAY:
                case FR_HR:
                case FR_MIN:
                case FR_SEC:
                case FR_MS:
                case FR_US:
                case FR_NS:
                    return &asfreq_MtoDT;
                default:
                    return &nofunc;
            }

        case FR_WK:
            switch (toGroup) {
                case FR_ANN:
                    return &asfreq_WtoA;
                case FR_QTR:
                    return &asfreq_WtoQ;
                case FR_MTH:
                    return &asfreq_WtoM;
                case FR_WK:
                    return &asfreq_WtoW;
                case FR_BUS:
                    return &asfreq_WtoB;
                case FR_DAY:
                case FR_HR:
                case FR_MIN:
                case FR_SEC:
                case FR_MS:
                case FR_US:
                case FR_NS:
                    return &asfreq_WtoDT;
                default:
                    return &nofunc;
            }

        case FR_BUS:
            switch (toGroup) {
                case FR_ANN:
                    return &asfreq_BtoA;
                case FR_QTR:
                    return &asfreq_BtoQ;
                case FR_MTH:
                    return &asfreq_BtoM;
                case FR_WK:
                    return &asfreq_BtoW;
                case FR_BUS:
                    return &no_op;
                case FR_DAY:
                case FR_HR:
                case FR_MIN:
                case FR_SEC:
                case FR_MS:
                case FR_US:
                case FR_NS:
                    return &asfreq_BtoDT;
                default:
                    return &nofunc;
            }

        case FR_DAY:
        case FR_HR:
        case FR_MIN:
        case FR_SEC:
        case FR_MS:
        case FR_US:
        case FR_NS:
            switch (toGroup) {
                case FR_ANN:
                    return &asfreq_DTtoA;
                case FR_QTR:
                    return &asfreq_DTtoQ;
                case FR_MTH:
                    return &asfreq_DTtoM;
                case FR_WK:
                    return &asfreq_DTtoW;
                case FR_BUS:
                    return &asfreq_DTtoB;
                case FR_DAY:
                case FR_HR:
                case FR_MIN:
                case FR_SEC:
                case FR_MS:
                case FR_US:
                case FR_NS:
                    if (fromGroup > toGroup) {
                        return &asfreq_DownsampleWithinDay;
                    } else {
                        return &asfreq_UpsampleWithinDay;
                    }
                default:
                    return &nofunc;
            }

        default:
            return &nofunc;
    }
}

double get_abs_time(int freq, npy_int64 date_ordinal, npy_int64 ordinal) {
    int freq_index, day_index, base_index;
    npy_int64 per_day, start_ord;
    double unit, result;

    if (freq <= FR_DAY) {
        return 0;
    }

    freq_index = get_freq_group_index(freq);
    day_index = get_freq_group_index(FR_DAY);
    base_index = get_freq_group_index(FR_SEC);

    per_day = get_daytime_conversion_factor(day_index, freq_index);
    unit = get_daytime_conversion_factor(freq_index, base_index);

    if (base_index < freq_index) {
        unit = 1 / unit;
    }

    start_ord = date_ordinal * per_day;
    result = (double)(unit * (ordinal - start_ord));
    return result;
}

/* Sets the time part of the DateTime object. */
static int dInfoCalc_SetFromAbsTime(struct date_info *dinfo, double abstime) {
    int inttime;
    int hour, minute;
    double second;

    inttime = (int)abstime;
    hour = inttime / 3600;
    minute = (inttime % 3600) / 60;
    second = abstime - (double)(hour * 3600 + minute * 60);

    dinfo->hour = hour;
    dinfo->minute = minute;
    dinfo->second = second;
    return 0;
}

/* Set the instance's value using the given date and time.
   Assumes GREGORIAN_CALENDAR. */
static int dInfoCalc_SetFromAbsDateTime(struct date_info *dinfo,
                                        npy_int64 absdate, double abstime) {
    /* Bounds check */
    // The calling function is responsible for ensuring that
    // abstime >= 0.0 && abstime <= 86400

    /* Calculate the date */
    dInfoCalc_SetFromAbsDate(dinfo, absdate);

    /* Calculate the time */
    dInfoCalc_SetFromAbsTime(dinfo, abstime);

    return 0;
}

/* ------------------------------------------------------------------
 * New pandas API-helper code, to expose to cython
 * ------------------------------------------------------------------*/

npy_int64 asfreq(npy_int64 period_ordinal, int freq1, int freq2,
                 char relation) {
    npy_int64 val;
    freq_conv_func func;
    asfreq_info finfo;

    func = get_asfreq_func(freq1, freq2);

    get_asfreq_info(freq1, freq2, relation, &finfo);
    val = (*func)(period_ordinal, &finfo);
    return val;
}

/* generate an ordinal in period space */
npy_int64 get_period_ordinal(int year, int month, int day, int hour, int minute,
                             int second, int microseconds, int picoseconds,
                             int freq) {
    npy_int64 absdays, delta, seconds;
    npy_int64 weeks, days;
    npy_int64 ordinal, day_adj;
    int freq_group, fmonth, mdiff;
    freq_group = get_freq_group(freq);

    if (freq == FR_SEC || freq == FR_MS || freq == FR_US || freq == FR_NS) {
        absdays = absdate_from_ymd(year, month, day);
        delta = (absdays - ORD_OFFSET);
        seconds =
            (npy_int64)(delta * 86400 + hour * 3600 + minute * 60 + second);

        switch (freq) {
            case FR_MS:
                return seconds * 1000 + microseconds / 1000;

            case FR_US:
                return seconds * 1000000 + microseconds;

            case FR_NS:
                return seconds * 1000000000 + microseconds * 1000 +
                       picoseconds / 1000;
        }

        return seconds;
    }

    if (freq == FR_MIN) {
        absdays = absdate_from_ymd(year, month, day);
        delta = (absdays - ORD_OFFSET);
        return (npy_int64)(delta * 1440 + hour * 60 + minute);
    }

    if (freq == FR_HR) {
        absdays = absdate_from_ymd(year, month, day);
        delta = (absdays - ORD_OFFSET);
        return (npy_int64)(delta * 24 + hour);
    }

    if (freq == FR_DAY) {
        return (npy_int64)(absdate_from_ymd(year, month, day) - ORD_OFFSET);
    }

    if (freq == FR_UND) {
        return (npy_int64)(absdate_from_ymd(year, month, day) - ORD_OFFSET);
    }

    if (freq == FR_BUS) {
        days = absdate_from_ymd(year, month, day);
        // calculate the current week assuming sunday as last day of a week
        weeks = (days - BASE_WEEK_TO_DAY_OFFSET) / DAYS_PER_WEEK;
        // calculate the current weekday (in range 1 .. 7)
        delta = (days - BASE_WEEK_TO_DAY_OFFSET) % DAYS_PER_WEEK + 1;
        // return the number of business days in full weeks plus the business
        // days in the last - possible partial - week
        return (npy_int64)(weeks * BUSINESS_DAYS_PER_WEEK) +
               (delta <= BUSINESS_DAYS_PER_WEEK ? delta
                                                : BUSINESS_DAYS_PER_WEEK + 1) -
               BDAY_OFFSET;
    }

    if (freq_group == FR_WK) {
        ordinal = (npy_int64)absdate_from_ymd(year, month, day);
        day_adj = freq - FR_WK;
        return (ordinal - (1 + day_adj)) / 7 + 1 - WEEK_OFFSET;
    }

    if (freq == FR_MTH) {
        return (year - BASE_YEAR) * 12 + month - 1;
    }

    if (freq_group == FR_QTR) {
        fmonth = freq - FR_QTR;
        if (fmonth == 0) fmonth = 12;

        mdiff = month - fmonth;
        if (mdiff < 0) mdiff += 12;
        if (month >= fmonth) mdiff += 12;

        return (year - BASE_YEAR) * 4 + (mdiff - 1) / 3;
    }

    if (freq_group == FR_ANN) {
        fmonth = freq - FR_ANN;
        if (fmonth == 0) fmonth = 12;
        if (month <= fmonth) {
            return year - BASE_YEAR;
        } else {
            return year - BASE_YEAR + 1;
        }
    }

    Py_Error(PyExc_RuntimeError, "Unable to generate frequency ordinal");

onError:
    return INT_ERR_CODE;
}

/*
   Returns the proleptic Gregorian ordinal of the date, as an integer.
   This corresponds to the number of days since Jan., 1st, 1AD.
   When the instance has a frequency less than daily, the proleptic date
   is calculated for the last day of the period.
 */

npy_int64 get_python_ordinal(npy_int64 period_ordinal, int freq) {
    asfreq_info af_info;
    freq_conv_func toDaily = NULL;

    if (freq == FR_DAY) return period_ordinal + ORD_OFFSET;

    toDaily = get_asfreq_func(freq, FR_DAY);
    get_asfreq_info(freq, FR_DAY, 'E', &af_info);

    return toDaily(period_ordinal, &af_info) + ORD_OFFSET;
}


int get_yq(npy_int64 ordinal, int freq, int *quarter, int *year) {
    asfreq_info af_info;
    int qtr_freq;
    npy_int64 daily_ord;

    daily_ord = get_python_ordinal(ordinal, freq) - ORD_OFFSET;

    if (get_freq_group(freq) == FR_QTR) {
        qtr_freq = freq;
    } else {
        qtr_freq = FR_QTR;
    }
    get_asfreq_info(FR_DAY, qtr_freq, 'E', &af_info);

    DtoQ_yq(daily_ord, &af_info, year, quarter);
    return qtr_freq;
}

int _quarter_year(npy_int64 ordinal, int freq, int *year, int *quarter) {
    asfreq_info af_info;
    int qtr_freq;

    ordinal = get_python_ordinal(ordinal, freq) - ORD_OFFSET;

    if (get_freq_group(freq) == FR_QTR)
        qtr_freq = freq;
    else
        qtr_freq = FR_QTR;

    get_asfreq_info(FR_DAY, qtr_freq, 'E', &af_info);

    DtoQ_yq(ordinal, &af_info, year, quarter);

    if ((qtr_freq % 1000) > 12) *year -= 1;

    return 0;
}


int get_date_info(npy_int64 ordinal, int freq, struct date_info *dinfo) {
    npy_int64 absdate = get_python_ordinal(ordinal, freq);
    double abstime = get_abs_time(freq, absdate - ORD_OFFSET, ordinal);

    while (abstime < 0) {
        abstime += 86400;
        absdate -= 1;
    }
    while (abstime >= 86400) {
        abstime -= 86400;
        absdate += 1;
    }

    dInfoCalc_SetFromAbsDateTime(dinfo, absdate, abstime);
    return 0;
}
