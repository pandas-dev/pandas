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
npy_int64 absdate_from_ymd(int year, int month, int day) {
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
    return floordiv(absdate, 7) * 5 + mod_compat(absdate, 7) - BDAY_OFFSET;
}

static npy_int64 DtoB(struct date_info *dinfo,
                      int roll_back, npy_int64 absdate) {
    int day_of_week = dayofweek(dinfo->year, dinfo->month, dinfo->day);

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
    npy_int64 absdate;
    int roll_back;

    ordinal = downsample_daytime(ordinal, af_info);
    absdate = ordinal + ORD_OFFSET;
    dInfoCalc_SetFromAbsDate(&dinfo, absdate);

    // This usage defines roll_back the opposite way from the others
    roll_back = 1 - af_info->is_end;
    return DtoB(&dinfo, roll_back, absdate);
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
        (floordiv(ordinal - 1, 5) * 7 + mod_compat(ordinal - 1, 5) + 1 - ORD_OFFSET);

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
    npy_int64 absdate = asfreq_WtoDT(ordinal, af_info) + ORD_OFFSET;
    int roll_back = af_info->is_end;
    dInfoCalc_SetFromAbsDate(&dinfo, absdate);

    return DtoB(&dinfo, roll_back, absdate);
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
    npy_int64 absdate = asfreq_MtoDT(ordinal, af_info) + ORD_OFFSET;
    int roll_back = af_info->is_end;

    dInfoCalc_SetFromAbsDate(&dinfo, absdate);

    return DtoB(&dinfo, roll_back, absdate);
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
    npy_int64 absdate = asfreq_QtoDT(ordinal, af_info) + ORD_OFFSET;
    int roll_back = af_info->is_end;

    dInfoCalc_SetFromAbsDate(&dinfo, absdate);

    return DtoB(&dinfo, roll_back, absdate);
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
    npy_int64 absdate = asfreq_AtoDT(ordinal, af_info) + ORD_OFFSET;
    int roll_back = af_info->is_end;
    dInfoCalc_SetFromAbsDate(&dinfo, absdate);

    return DtoB(&dinfo, roll_back, absdate);
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
