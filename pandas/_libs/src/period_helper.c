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


/* Find the unix_date (days elapsed since datetime(1970, 1, 1)
 * for the given year/month/day.
 * Assumes GREGORIAN_CALENDAR */
npy_int64 unix_date_from_ymd(int year, int month, int day) {
    /* Calculate the absolute date */
    npy_datetimestruct dts;
    npy_int64 unix_date;

    memset(&dts, 0, sizeof(npy_datetimestruct));
    dts.year = year;
    dts.month = month;
    dts.day = day;
    unix_date = npy_datetimestruct_to_datetime(NPY_FR_D, &dts);
    return unix_date;
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

int max_value(int a, int b) { return a > b ? a : b; }

PANDAS_INLINE int min_value(int a, int b) { return a < b ? a : b; }

PANDAS_INLINE int get_freq_group(int freq) { return (freq / 1000) * 1000; }


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

static npy_int64 DtoB_weekday(npy_int64 unix_date) {
    return floordiv(unix_date + 4, 7) * 5 + mod_compat(unix_date + 4, 7) - 4;
}

static npy_int64 DtoB(npy_datetimestruct *dts,
                      int roll_back, npy_int64 unix_date) {
    int day_of_week = dayofweek(dts->year, dts->month, dts->day);

    if (roll_back == 1) {
        if (day_of_week > 4) {
            // change to friday before weekend
            unix_date -= (day_of_week - 4);
        }
    } else {
        if (day_of_week > 4) {
            // change to Monday after weekend
            unix_date += (7 - day_of_week);
        }
    }
    return DtoB_weekday(unix_date);
}


//************ FROM DAILY ***************

static npy_int64 asfreq_DTtoA(npy_int64 ordinal, asfreq_info *af_info) {
    npy_datetimestruct dts;
    ordinal = downsample_daytime(ordinal, af_info);
    pandas_datetime_to_datetimestruct(ordinal, NPY_FR_D, &dts);
    if (dts.month > af_info->to_end) {
        return (npy_int64)(dts.year + 1 - 1970);
    } else {
        return (npy_int64)(dts.year - 1970);
    }
}

static int DtoQ_yq(npy_int64 ordinal, asfreq_info *af_info, int *year) {
    npy_datetimestruct dts;
    int quarter;

    pandas_datetime_to_datetimestruct(ordinal, NPY_FR_D, &dts);
    if (af_info->to_end != 12) {
        dts.month -= af_info->to_end;
        if (dts.month <= 0) {
            dts.month += 12;
        } else {
            dts.year += 1;
        }
    }

    *year = dts.year;
    quarter = monthToQuarter(dts.month);
    return quarter;
}

static npy_int64 asfreq_DTtoQ(npy_int64 ordinal, asfreq_info *af_info) {
    int year, quarter;

    ordinal = downsample_daytime(ordinal, af_info);

    quarter = DtoQ_yq(ordinal, af_info, &year);
    return (npy_int64)((year - 1970) * 4 + quarter - 1);
}

static npy_int64 asfreq_DTtoM(npy_int64 ordinal, asfreq_info *af_info) {
    npy_datetimestruct dts;

    ordinal = downsample_daytime(ordinal, af_info);

    pandas_datetime_to_datetimestruct(ordinal, NPY_FR_D, &dts);
    return (npy_int64)((dts.year - 1970) * 12 + dts.month - 1);
}

static npy_int64 asfreq_DTtoW(npy_int64 ordinal, asfreq_info *af_info) {
    ordinal = downsample_daytime(ordinal, af_info);
    return floordiv(ordinal + 3 - af_info->to_end, 7) + 1;
}

static npy_int64 asfreq_DTtoB(npy_int64 ordinal, asfreq_info *af_info) {
    int roll_back;
    npy_datetimestruct dts;
    npy_int64 unix_date = downsample_daytime(ordinal, af_info);
    pandas_datetime_to_datetimestruct(unix_date, NPY_FR_D, &dts);

    // This usage defines roll_back the opposite way from the others
    roll_back = 1 - af_info->is_end;
    return DtoB(&dts, roll_back, unix_date);
}

//************ FROM BUSINESS ***************

static npy_int64 asfreq_BtoDT(npy_int64 ordinal, asfreq_info *af_info) {
    ordinal = floordiv(ordinal + 3, 5) * 7 + mod_compat(ordinal + 3, 5) - 3;

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
    ordinal = ordinal * 7 + af_info->from_end - 4 +
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
    int roll_back;
    npy_datetimestruct dts;
    npy_int64 unix_date = asfreq_WtoDT(ordinal, af_info);

    pandas_datetime_to_datetimestruct(unix_date, NPY_FR_D, &dts);
    roll_back = af_info->is_end;
    return DtoB(&dts, roll_back, unix_date);
}

//************ FROM MONTHLY ***************
static void MtoD_ym(npy_int64 ordinal, int *year, int *month) {
    *year = floordiv(ordinal, 12) + 1970;
    *month = mod_compat(ordinal, 12) + 1;
}

static npy_int64 asfreq_MtoDT(npy_int64 ordinal, asfreq_info *af_info) {
    npy_int64 unix_date;
    int year, month;

    ordinal += af_info->is_end;
    MtoD_ym(ordinal, &year, &month);

    unix_date = unix_date_from_ymd(year, month, 1);
    unix_date -= af_info->is_end;
    return upsample_daytime(unix_date, af_info);
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
    int roll_back;
    npy_datetimestruct dts;
    npy_int64 unix_date = asfreq_MtoDT(ordinal, af_info);

    pandas_datetime_to_datetimestruct(unix_date, NPY_FR_D, &dts);
    roll_back = af_info->is_end;
    return DtoB(&dts, roll_back, unix_date);
}

//************ FROM QUARTERLY ***************

static void QtoD_ym(npy_int64 ordinal, int *year, int *month,
                    asfreq_info *af_info) {
    *year = floordiv(ordinal, 4) + 1970;
    *month = mod_compat(ordinal, 4) * 3 + 1;

    if (af_info->from_end != 12) {
        *month += af_info->from_end;
        if (*month > 12) {
            *month -= 12;
        } else {
            *year -= 1;
        }
    }
}

static npy_int64 asfreq_QtoDT(npy_int64 ordinal, asfreq_info *af_info) {
    npy_int64 unix_date;
    int year, month;

    ordinal += af_info->is_end;
    QtoD_ym(ordinal, &year, &month, af_info);

    unix_date = unix_date_from_ymd(year, month, 1);
    unix_date -= af_info->is_end;
    return upsample_daytime(unix_date, af_info);
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
    int roll_back;
    npy_datetimestruct dts;
    npy_int64 unix_date = asfreq_QtoDT(ordinal, af_info);

    pandas_datetime_to_datetimestruct(unix_date, NPY_FR_D, &dts);
    roll_back = af_info->is_end;
    return DtoB(&dts, roll_back, unix_date);
}

//************ FROM ANNUAL ***************

static void AtoD_ym(npy_int64 ordinal, npy_int64 *year, int *month,
                    asfreq_info *af_info) {
    *year = ordinal + 1970;
    *month = 1;

    if (af_info->from_end != 12) {
        *month += af_info->from_end;
        if (*month > 12) {
            // This case is never reached, but is kept for symmetry
            // with QtoD_ym
            *month -= 12;
        } else {
            *year -= 1;
        }
    }
}

static npy_int64 asfreq_AtoDT(npy_int64 ordinal, asfreq_info *af_info) {
    npy_int64 unix_date, year;
    int month;

    ordinal += af_info->is_end;
    AtoD_ym(ordinal, &year, &month, af_info);

    unix_date = unix_date_from_ymd(year, month, 1);
    unix_date -= af_info->is_end;
    return upsample_daytime(unix_date, af_info);
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
    int roll_back;
    npy_datetimestruct dts;
    npy_int64 unix_date = asfreq_AtoDT(ordinal, af_info);

    pandas_datetime_to_datetimestruct(unix_date, NPY_FR_D, &dts);
    roll_back = af_info->is_end;
    return DtoB(&dts, roll_back, unix_date);
}

static npy_int64 nofunc(npy_int64 ordinal, asfreq_info *af_info) {
    return INT_ERR_CODE;
}
static npy_int64 no_op(npy_int64 ordinal, asfreq_info *af_info) {
    return ordinal;
}

// end of frequency specific conversion routines

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
                        return &downsample_daytime;
                    } else {
                        return &upsample_daytime;
                    }
                default:
                    return &nofunc;
            }

        default:
            return &nofunc;
    }
}
