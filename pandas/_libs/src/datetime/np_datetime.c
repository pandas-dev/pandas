/*

Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.

Copyright (c) 2005-2011, NumPy Developers
All rights reserved.

This file is derived from NumPy 1.7. See NUMPY_LICENSE.txt

*/

#define NO_IMPORT

#include <Python.h>
#include <datetime.h>

#include <numpy/arrayobject.h>
#include <numpy/arrayscalars.h>
#include "np_datetime.h"

#if PY_MAJOR_VERSION >= 3
#define PyInt_AsLong PyLong_AsLong
#endif

const pandas_datetimestruct _NS_MIN_DTS = {
    1677, 9, 21, 0, 12, 43, 145225, 0, 0};
const pandas_datetimestruct _NS_MAX_DTS = {
    2262, 4, 11, 23, 47, 16, 854775, 807000, 0};


const int days_per_month_table[2][12] = {
    {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}};

/*
 * Returns 1 if the given year is a leap year, 0 otherwise.
 */
int is_leapyear(npy_int64 year) {
    return (year & 0x3) == 0 && /* year % 4 == 0 */
           ((year % 100) != 0 || (year % 400) == 0);
}

/*
 * Sakamoto's method, from wikipedia
 */
int dayofweek(int y, int m, int d) {
    int day;
    static const int t[] = {0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4};
    y -= m < 3;
    day = (y + y / 4 - y / 100 + y / 400 + t[m - 1] + d) % 7;
    // convert to python day
    return (day + 6) % 7;
}

/*
 * Adjusts a datetimestruct based on a minutes offset. Assumes
 * the current values are valid.g
 */
void add_minutes_to_datetimestruct(pandas_datetimestruct *dts, int minutes) {
    int isleap;

    /* MINUTES */
    dts->min += minutes;
    while (dts->min < 0) {
        dts->min += 60;
        dts->hour--;
    }
    while (dts->min >= 60) {
        dts->min -= 60;
        dts->hour++;
    }

    /* HOURS */
    while (dts->hour < 0) {
        dts->hour += 24;
        dts->day--;
    }
    while (dts->hour >= 24) {
        dts->hour -= 24;
        dts->day++;
    }

    /* DAYS */
    if (dts->day < 1) {
        dts->month--;
        if (dts->month < 1) {
            dts->year--;
            dts->month = 12;
        }
        isleap = is_leapyear(dts->year);
        dts->day += days_per_month_table[isleap][dts->month - 1];
    } else if (dts->day > 28) {
        isleap = is_leapyear(dts->year);
        if (dts->day > days_per_month_table[isleap][dts->month - 1]) {
            dts->day -= days_per_month_table[isleap][dts->month - 1];
            dts->month++;
            if (dts->month > 12) {
                dts->year++;
                dts->month = 1;
            }
        }
    }
}

/*
 * Calculates the days offset from the 1970 epoch.
 */
npy_int64 get_datetimestruct_days(const pandas_datetimestruct *dts) {
    int i, month;
    npy_int64 year, days = 0;
    const int *month_lengths;

    year = dts->year - 1970;
    days = year * 365;

    /* Adjust for leap years */
    if (days >= 0) {
        /*
         * 1968 is the closest leap year before 1970.
         * Exclude the current year, so add 1.
         */
        year += 1;
        /* Add one day for each 4 years */
        days += year / 4;
        /* 1900 is the closest previous year divisible by 100 */
        year += 68;
        /* Subtract one day for each 100 years */
        days -= year / 100;
        /* 1600 is the closest previous year divisible by 400 */
        year += 300;
        /* Add one day for each 400 years */
        days += year / 400;
    } else {
        /*
         * 1972 is the closest later year after 1970.
         * Include the current year, so subtract 2.
         */
        year -= 2;
        /* Subtract one day for each 4 years */
        days += year / 4;
        /* 2000 is the closest later year divisible by 100 */
        year -= 28;
        /* Add one day for each 100 years */
        days -= year / 100;
        /* 2000 is also the closest later year divisible by 400 */
        /* Subtract one day for each 400 years */
        days += year / 400;
    }

    month_lengths = days_per_month_table[is_leapyear(dts->year)];
    month = dts->month - 1;

    /* Add the months */
    for (i = 0; i < month; ++i) {
        days += month_lengths[i];
    }

    /* Add the days */
    days += dts->day - 1;

    return days;
}

/*
 * Modifies '*days_' to be the day offset within the year,
 * and returns the year.
 */
static npy_int64 days_to_yearsdays(npy_int64 *days_) {
    const npy_int64 days_per_400years = (400 * 365 + 100 - 4 + 1);
    /* Adjust so it's relative to the year 2000 (divisible by 400) */
    npy_int64 days = (*days_) - (365 * 30 + 7);
    npy_int64 year;

    /* Break down the 400 year cycle to get the year and day within the year */
    if (days >= 0) {
        year = 400 * (days / days_per_400years);
        days = days % days_per_400years;
    } else {
        year = 400 * ((days - (days_per_400years - 1)) / days_per_400years);
        days = days % days_per_400years;
        if (days < 0) {
            days += days_per_400years;
        }
    }

    /* Work out the year/day within the 400 year cycle */
    if (days >= 366) {
        year += 100 * ((days - 1) / (100 * 365 + 25 - 1));
        days = (days - 1) % (100 * 365 + 25 - 1);
        if (days >= 365) {
            year += 4 * ((days + 1) / (4 * 365 + 1));
            days = (days + 1) % (4 * 365 + 1);
            if (days >= 366) {
                year += (days - 1) / 365;
                days = (days - 1) % 365;
            }
        }
    }

    *days_ = days;
    return year + 2000;
}

/*
 * Adjusts a datetimestruct based on a seconds offset. Assumes
 * the current values are valid.
 */
NPY_NO_EXPORT void add_seconds_to_datetimestruct(pandas_datetimestruct *dts,
                                                 int seconds) {
    int minutes;

    dts->sec += seconds;
    if (dts->sec < 0) {
        minutes = dts->sec / 60;
        dts->sec = dts->sec % 60;
        if (dts->sec < 0) {
            --minutes;
            dts->sec += 60;
        }
        add_minutes_to_datetimestruct(dts, minutes);
    } else if (dts->sec >= 60) {
        minutes = dts->sec / 60;
        dts->sec = dts->sec % 60;
        add_minutes_to_datetimestruct(dts, minutes);
    }
}

/*
 * Fills in the year, month, day in 'dts' based on the days
 * offset from 1970.
 */
static void set_datetimestruct_days(npy_int64 days,
                                    pandas_datetimestruct *dts) {
    const int *month_lengths;
    int i;

    dts->year = days_to_yearsdays(&days);
    month_lengths = days_per_month_table[is_leapyear(dts->year)];

    for (i = 0; i < 12; ++i) {
        if (days < month_lengths[i]) {
            dts->month = i + 1;
            dts->day = days + 1;
            return;
        } else {
            days -= month_lengths[i];
        }
    }
}

/*
 * Compares two pandas_datetimestruct objects chronologically
 */
int cmp_pandas_datetimestruct(const pandas_datetimestruct *a,
                              const pandas_datetimestruct *b) {
    if (a->year > b->year) {
        return 1;
    } else if (a->year < b->year) {
        return -1;
    }

    if (a->month > b->month) {
        return 1;
    } else if (a->month < b->month) {
        return -1;
    }

    if (a->day > b->day) {
        return 1;
    } else if (a->day < b->day) {
        return -1;
    }

    if (a->hour > b->hour) {
        return 1;
    } else if (a->hour < b->hour) {
        return -1;
    }

    if (a->min > b->min) {
        return 1;
    } else if (a->min < b->min) {
        return -1;
    }

    if (a->sec > b->sec) {
        return 1;
    } else if (a->sec < b->sec) {
        return -1;
    }

    if (a->us > b->us) {
        return 1;
    } else if (a->us < b->us) {
        return -1;
    }

    if (a->ps > b->ps) {
        return 1;
    } else if (a->ps < b->ps) {
        return -1;
    }

    if (a->as > b->as) {
        return 1;
    } else if (a->as < b->as) {
        return -1;
    }

    return 0;
}

/*
 *
 * Tests for and converts a Python datetime.datetime or datetime.date
 * object into a NumPy pandas_datetimestruct.
 *
 * While the C API has PyDate_* and PyDateTime_* functions, the following
 * implementation just asks for attributes, and thus supports
 * datetime duck typing. The tzinfo time zone conversion would require
 * this style of access anyway.
 *
 * 'out_bestunit' gives a suggested unit based on whether the object
 *      was a datetime.date or datetime.datetime object.
 *
 * If 'apply_tzinfo' is 1, this function uses the tzinfo to convert
 * to UTC time, otherwise it returns the struct with the local time.
 *
 * Returns -1 on error, 0 on success, and 1 (with no error set)
 * if obj doesn't have the neeeded date or datetime attributes.
 */
int convert_pydatetime_to_datetimestruct(PyObject *obj,
                                         pandas_datetimestruct *out,
                                         PANDAS_DATETIMEUNIT *out_bestunit,
                                         int apply_tzinfo) {
    PyObject *tmp;
    int isleap;

    /* Initialize the output to all zeros */
    memset(out, 0, sizeof(pandas_datetimestruct));
    out->month = 1;
    out->day = 1;

    /* Need at least year/month/day attributes */
    if (!PyObject_HasAttrString(obj, "year") ||
        !PyObject_HasAttrString(obj, "month") ||
        !PyObject_HasAttrString(obj, "day")) {
        return 1;
    }

    /* Get the year */
    tmp = PyObject_GetAttrString(obj, "year");
    if (tmp == NULL) {
        return -1;
    }
    out->year = PyInt_AsLong(tmp);
    if (out->year == -1 && PyErr_Occurred()) {
        Py_DECREF(tmp);
        return -1;
    }
    Py_DECREF(tmp);

    /* Get the month */
    tmp = PyObject_GetAttrString(obj, "month");
    if (tmp == NULL) {
        return -1;
    }
    out->month = PyInt_AsLong(tmp);
    if (out->month == -1 && PyErr_Occurred()) {
        Py_DECREF(tmp);
        return -1;
    }
    Py_DECREF(tmp);

    /* Get the day */
    tmp = PyObject_GetAttrString(obj, "day");
    if (tmp == NULL) {
        return -1;
    }
    out->day = PyInt_AsLong(tmp);
    if (out->day == -1 && PyErr_Occurred()) {
        Py_DECREF(tmp);
        return -1;
    }
    Py_DECREF(tmp);

    /* Validate that the month and day are valid for the year */
    if (out->month < 1 || out->month > 12) {
        goto invalid_date;
    }
    isleap = is_leapyear(out->year);
    if (out->day < 1 ||
        out->day > days_per_month_table[isleap][out->month - 1]) {
        goto invalid_date;
    }

    /* Check for time attributes (if not there, return success as a date) */
    if (!PyObject_HasAttrString(obj, "hour") ||
        !PyObject_HasAttrString(obj, "minute") ||
        !PyObject_HasAttrString(obj, "second") ||
        !PyObject_HasAttrString(obj, "microsecond")) {
        /* The best unit for date is 'D' */
        if (out_bestunit != NULL) {
            *out_bestunit = PANDAS_FR_D;
        }
        return 0;
    }

    /* Get the hour */
    tmp = PyObject_GetAttrString(obj, "hour");
    if (tmp == NULL) {
        return -1;
    }
    out->hour = PyInt_AsLong(tmp);
    if (out->hour == -1 && PyErr_Occurred()) {
        Py_DECREF(tmp);
        return -1;
    }
    Py_DECREF(tmp);

    /* Get the minute */
    tmp = PyObject_GetAttrString(obj, "minute");
    if (tmp == NULL) {
        return -1;
    }
    out->min = PyInt_AsLong(tmp);
    if (out->min == -1 && PyErr_Occurred()) {
        Py_DECREF(tmp);
        return -1;
    }
    Py_DECREF(tmp);

    /* Get the second */
    tmp = PyObject_GetAttrString(obj, "second");
    if (tmp == NULL) {
        return -1;
    }
    out->sec = PyInt_AsLong(tmp);
    if (out->sec == -1 && PyErr_Occurred()) {
        Py_DECREF(tmp);
        return -1;
    }
    Py_DECREF(tmp);

    /* Get the microsecond */
    tmp = PyObject_GetAttrString(obj, "microsecond");
    if (tmp == NULL) {
        return -1;
    }
    out->us = PyInt_AsLong(tmp);
    if (out->us == -1 && PyErr_Occurred()) {
        Py_DECREF(tmp);
        return -1;
    }
    Py_DECREF(tmp);

    if (out->hour < 0 || out->hour >= 24 || out->min < 0 || out->min >= 60 ||
        out->sec < 0 || out->sec >= 60 || out->us < 0 || out->us >= 1000000) {
        goto invalid_time;
    }

    /* Apply the time zone offset if it exists */
    if (apply_tzinfo && PyObject_HasAttrString(obj, "tzinfo")) {
        tmp = PyObject_GetAttrString(obj, "tzinfo");
        if (tmp == NULL) {
            return -1;
        }
        if (tmp == Py_None) {
            Py_DECREF(tmp);
        } else {
            PyObject *offset;
            int seconds_offset, minutes_offset;

            /* The utcoffset function should return a timedelta */
            offset = PyObject_CallMethod(tmp, "utcoffset", "O", obj);
            if (offset == NULL) {
                Py_DECREF(tmp);
                return -1;
            }
            Py_DECREF(tmp);

            /*
             * The timedelta should have a function "total_seconds"
             * which contains the value we want.
             */
            tmp = PyObject_CallMethod(offset, "total_seconds", "");
            if (tmp == NULL) {
                return -1;
            }
            seconds_offset = PyInt_AsLong(tmp);
            if (seconds_offset == -1 && PyErr_Occurred()) {
                Py_DECREF(tmp);
                return -1;
            }
            Py_DECREF(tmp);

            /* Convert to a minutes offset and apply it */
            minutes_offset = seconds_offset / 60;

            add_minutes_to_datetimestruct(out, -minutes_offset);
        }
    }

    /* The resolution of Python's datetime is 'us' */
    if (out_bestunit != NULL) {
        *out_bestunit = PANDAS_FR_us;
    }

    return 0;

invalid_date:
    PyErr_Format(PyExc_ValueError,
                 "Invalid date (%d,%d,%d) when converting to NumPy datetime",
                 (int)out->year, (int)out->month, (int)out->day);
    return -1;

invalid_time:
    PyErr_Format(PyExc_ValueError,
                 "Invalid time (%d,%d,%d,%d) when converting "
                 "to NumPy datetime",
                 (int)out->hour, (int)out->min, (int)out->sec, (int)out->us);
    return -1;
}

npy_datetime pandas_datetimestruct_to_datetime(PANDAS_DATETIMEUNIT fr,
                                               pandas_datetimestruct *d) {
    pandas_datetime_metadata meta;
    npy_datetime result = PANDAS_DATETIME_NAT;

    meta.base = fr;
    meta.num = 1;

    convert_datetimestruct_to_datetime(&meta, d, &result);
    return result;
}

void pandas_datetime_to_datetimestruct(npy_datetime val, PANDAS_DATETIMEUNIT fr,
                                       pandas_datetimestruct *result) {
    pandas_datetime_metadata meta;

    meta.base = fr;
    meta.num = 1;

    convert_datetime_to_datetimestruct(&meta, val, result);
}

void pandas_timedelta_to_timedeltastruct(npy_timedelta val,
                                         PANDAS_DATETIMEUNIT fr,
                                         pandas_timedeltastruct *result) {
  pandas_datetime_metadata meta;

  meta.base = fr;
  meta.num = 1;

  convert_timedelta_to_timedeltastruct(&meta, val, result);
}


/*
 * Converts a datetime from a datetimestruct to a datetime based
 * on some metadata. The date is assumed to be valid.
 *
 * TODO: If meta->num is really big, there could be overflow
 *
 * Returns 0 on success, -1 on failure.
 */
int convert_datetimestruct_to_datetime(pandas_datetime_metadata *meta,
                                       const pandas_datetimestruct *dts,
                                       npy_datetime *out) {
    npy_datetime ret;
    PANDAS_DATETIMEUNIT base = meta->base;

    if (base == PANDAS_FR_Y) {
        /* Truncate to the year */
        ret = dts->year - 1970;
    } else if (base == PANDAS_FR_M) {
        /* Truncate to the month */
        ret = 12 * (dts->year - 1970) + (dts->month - 1);
    } else {
        /* Otherwise calculate the number of days to start */
        npy_int64 days = get_datetimestruct_days(dts);

        switch (base) {
            case PANDAS_FR_W:
                /* Truncate to weeks */
                if (days >= 0) {
                    ret = days / 7;
                } else {
                    ret = (days - 6) / 7;
                }
                break;
            case PANDAS_FR_D:
                ret = days;
                break;
            case PANDAS_FR_h:
                ret = days * 24 + dts->hour;
                break;
            case PANDAS_FR_m:
                ret = (days * 24 + dts->hour) * 60 + dts->min;
                break;
            case PANDAS_FR_s:
                ret = ((days * 24 + dts->hour) * 60 + dts->min) * 60 + dts->sec;
                break;
            case PANDAS_FR_ms:
                ret = (((days * 24 + dts->hour) * 60 + dts->min) * 60 +
                       dts->sec) *
                          1000 +
                      dts->us / 1000;
                break;
            case PANDAS_FR_us:
                ret = (((days * 24 + dts->hour) * 60 + dts->min) * 60 +
                       dts->sec) *
                          1000000 +
                      dts->us;
                break;
            case PANDAS_FR_ns:
                ret = ((((days * 24 + dts->hour) * 60 + dts->min) * 60 +
                        dts->sec) *
                           1000000 +
                       dts->us) *
                          1000 +
                      dts->ps / 1000;
                break;
            case PANDAS_FR_ps:
                ret = ((((days * 24 + dts->hour) * 60 + dts->min) * 60 +
                        dts->sec) *
                           1000000 +
                       dts->us) *
                          1000000 +
                      dts->ps;
                break;
            case PANDAS_FR_fs:
                /* only 2.6 hours */
                ret = (((((days * 24 + dts->hour) * 60 + dts->min) * 60 +
                         dts->sec) *
                            1000000 +
                        dts->us) *
                           1000000 +
                       dts->ps) *
                          1000 +
                      dts->as / 1000;
                break;
            case PANDAS_FR_as:
                /* only 9.2 secs */
                ret = (((((days * 24 + dts->hour) * 60 + dts->min) * 60 +
                         dts->sec) *
                            1000000 +
                        dts->us) *
                           1000000 +
                       dts->ps) *
                          1000000 +
                      dts->as;
                break;
            default:
                /* Something got corrupted */
                PyErr_SetString(
                    PyExc_ValueError,
                    "NumPy datetime metadata with corrupt unit value");
                return -1;
        }
    }

    /* Divide by the multiplier */
    if (meta->num > 1) {
        if (ret >= 0) {
            ret /= meta->num;
        } else {
            ret = (ret - meta->num + 1) / meta->num;
        }
    }

    *out = ret;

    return 0;
}

/*
 * Converts a datetime based on the given metadata into a datetimestruct
 */
int convert_datetime_to_datetimestruct(pandas_datetime_metadata *meta,
                                       npy_datetime dt,
                                       pandas_datetimestruct *out) {
    npy_int64 perday;

    /* Initialize the output to all zeros */
    memset(out, 0, sizeof(pandas_datetimestruct));
    out->year = 1970;
    out->month = 1;
    out->day = 1;

    /* TODO: Change to a mechanism that avoids the potential overflow */
    dt *= meta->num;

    /*
     * Note that care must be taken with the / and % operators
     * for negative values.
     */
    switch (meta->base) {
        case PANDAS_FR_Y:
            out->year = 1970 + dt;
            break;

        case PANDAS_FR_M:
            if (dt >= 0) {
                out->year = 1970 + dt / 12;
                out->month = dt % 12 + 1;
            } else {
                out->year = 1969 + (dt + 1) / 12;
                out->month = 12 + (dt + 1) % 12;
            }
            break;

        case PANDAS_FR_W:
            /* A week is 7 days */
            set_datetimestruct_days(dt * 7, out);
            break;

        case PANDAS_FR_D:
            set_datetimestruct_days(dt, out);
            break;

        case PANDAS_FR_h:
            perday = 24LL;

            if (dt >= 0) {
                set_datetimestruct_days(dt / perday, out);
                dt = dt % perday;
            } else {
                set_datetimestruct_days(
                    dt / perday - (dt % perday == 0 ? 0 : 1), out);
                dt = (perday - 1) + (dt + 1) % perday;
            }
            out->hour = dt;
            break;

        case PANDAS_FR_m:
            perday = 24LL * 60;

            if (dt >= 0) {
                set_datetimestruct_days(dt / perday, out);
                dt = dt % perday;
            } else {
                set_datetimestruct_days(
                    dt / perday - (dt % perday == 0 ? 0 : 1), out);
                dt = (perday - 1) + (dt + 1) % perday;
            }
            out->hour = dt / 60;
            out->min = dt % 60;
            break;

        case PANDAS_FR_s:
            perday = 24LL * 60 * 60;

            if (dt >= 0) {
                set_datetimestruct_days(dt / perday, out);
                dt = dt % perday;
            } else {
                set_datetimestruct_days(
                    dt / perday - (dt % perday == 0 ? 0 : 1), out);
                dt = (perday - 1) + (dt + 1) % perday;
            }
            out->hour = dt / (60 * 60);
            out->min = (dt / 60) % 60;
            out->sec = dt % 60;
            break;

        case PANDAS_FR_ms:
            perday = 24LL * 60 * 60 * 1000;

            if (dt >= 0) {
                set_datetimestruct_days(dt / perday, out);
                dt = dt % perday;
            } else {
                set_datetimestruct_days(
                    dt / perday - (dt % perday == 0 ? 0 : 1), out);
                dt = (perday - 1) + (dt + 1) % perday;
            }
            out->hour = dt / (60 * 60 * 1000LL);
            out->min = (dt / (60 * 1000LL)) % 60;
            out->sec = (dt / 1000LL) % 60;
            out->us = (dt % 1000LL) * 1000;
            break;

        case PANDAS_FR_us:
            perday = 24LL * 60LL * 60LL * 1000LL * 1000LL;

            if (dt >= 0) {
                set_datetimestruct_days(dt / perday, out);
                dt = dt % perday;
            } else {
                set_datetimestruct_days(
                    dt / perday - (dt % perday == 0 ? 0 : 1), out);
                dt = (perday - 1) + (dt + 1) % perday;
            }
            out->hour = dt / (60 * 60 * 1000000LL);
            out->min = (dt / (60 * 1000000LL)) % 60;
            out->sec = (dt / 1000000LL) % 60;
            out->us = dt % 1000000LL;
            break;

        case PANDAS_FR_ns:
            perday = 24LL * 60LL * 60LL * 1000LL * 1000LL * 1000LL;

            if (dt >= 0) {
                set_datetimestruct_days(dt / perday, out);
                dt = dt % perday;
            } else {
                set_datetimestruct_days(
                    dt / perday - (dt % perday == 0 ? 0 : 1), out);
                dt = (perday - 1) + (dt + 1) % perday;
            }
            out->hour = dt / (60 * 60 * 1000000000LL);
            out->min = (dt / (60 * 1000000000LL)) % 60;
            out->sec = (dt / 1000000000LL) % 60;
            out->us = (dt / 1000LL) % 1000000LL;
            out->ps = (dt % 1000LL) * 1000;
            break;

        case PANDAS_FR_ps:
            perday = 24LL * 60 * 60 * 1000 * 1000 * 1000 * 1000;

            if (dt >= 0) {
                set_datetimestruct_days(dt / perday, out);
                dt = dt % perday;
            } else {
                set_datetimestruct_days(
                    dt / perday - (dt % perday == 0 ? 0 : 1), out);
                dt = (perday - 1) + (dt + 1) % perday;
            }
            out->hour = dt / (60 * 60 * 1000000000000LL);
            out->min = (dt / (60 * 1000000000000LL)) % 60;
            out->sec = (dt / 1000000000000LL) % 60;
            out->us = (dt / 1000000LL) % 1000000LL;
            out->ps = dt % 1000000LL;
            break;

        case PANDAS_FR_fs:
            /* entire range is only +- 2.6 hours */
            if (dt >= 0) {
                out->hour = dt / (60 * 60 * 1000000000000000LL);
                out->min = (dt / (60 * 1000000000000000LL)) % 60;
                out->sec = (dt / 1000000000000000LL) % 60;
                out->us = (dt / 1000000000LL) % 1000000LL;
                out->ps = (dt / 1000LL) % 1000000LL;
                out->as = (dt % 1000LL) * 1000;
            } else {
                npy_datetime minutes;

                minutes = dt / (60 * 1000000000000000LL);
                dt = dt % (60 * 1000000000000000LL);
                if (dt < 0) {
                    dt += (60 * 1000000000000000LL);
                    --minutes;
                }
                /* Offset the negative minutes */
                add_minutes_to_datetimestruct(out, minutes);
                out->sec = (dt / 1000000000000000LL) % 60;
                out->us = (dt / 1000000000LL) % 1000000LL;
                out->ps = (dt / 1000LL) % 1000000LL;
                out->as = (dt % 1000LL) * 1000;
            }
            break;

        case PANDAS_FR_as:
            /* entire range is only +- 9.2 seconds */
            if (dt >= 0) {
                out->sec = (dt / 1000000000000000000LL) % 60;
                out->us = (dt / 1000000000000LL) % 1000000LL;
                out->ps = (dt / 1000000LL) % 1000000LL;
                out->as = dt % 1000000LL;
            } else {
                npy_datetime seconds;

                seconds = dt / 1000000000000000000LL;
                dt = dt % 1000000000000000000LL;
                if (dt < 0) {
                    dt += 1000000000000000000LL;
                    --seconds;
                }
                /* Offset the negative seconds */
                add_seconds_to_datetimestruct(out, seconds);
                out->us = (dt / 1000000000000LL) % 1000000LL;
                out->ps = (dt / 1000000LL) % 1000000LL;
                out->as = dt % 1000000LL;
            }
            break;

        default:
            PyErr_SetString(PyExc_RuntimeError,
                            "NumPy datetime metadata is corrupted with invalid "
                            "base unit");
            return -1;
    }

    return 0;
}

/*
 * Converts a timedelta from a timedeltastruct to a timedelta based
 * on some metadata. The timedelta is assumed to be valid.
 *
 * Returns 0 on success, -1 on failure.
 */
int convert_timedelta_to_timedeltastruct(pandas_timedelta_metadata *meta,
                                         npy_timedelta td,
                                         pandas_timedeltastruct *out) {
    npy_int64 frac;
    npy_int64 sfrac;
    npy_int64 ifrac;
    int sign;
    npy_int64 DAY_NS = 86400000000000LL;

    /* Initialize the output to all zeros */
    memset(out, 0, sizeof(pandas_timedeltastruct));

    switch (meta->base) {
        case PANDAS_FR_ns:

        // put frac in seconds
        if (td < 0 && td % (1000LL * 1000LL * 1000LL) != 0)
            frac = td / (1000LL * 1000LL * 1000LL) - 1;
        else
            frac = td / (1000LL * 1000LL * 1000LL);

        if (frac < 0) {
            sign = -1;

            // even fraction
            if ((-frac % 86400LL) != 0) {
              out->days = -frac / 86400LL + 1;
              frac += 86400LL * out->days;
            } else {
              frac = -frac;
            }
        } else {
            sign = 1;
            out->days = 0;
        }

        if (frac >= 86400) {
            out->days += frac / 86400LL;
            frac -= out->days * 86400LL;
        }

        if (frac >= 3600) {
            out->hrs = frac / 3600LL;
            frac -= out->hrs * 3600LL;
        } else {
            out->hrs = 0;
        }

        if (frac >= 60) {
            out->min = frac / 60LL;
            frac -= out->min * 60LL;
        } else {
            out->min = 0;
        }

        if (frac >= 0) {
            out->sec = frac;
            frac -= out->sec;
        } else {
            out->sec = 0;
        }

        sfrac = (out->hrs * 3600LL + out->min * 60LL
                 + out->sec) * (1000LL * 1000LL * 1000LL);

        if (sign < 0)
            out->days = -out->days;

        ifrac = td - (out->days * DAY_NS + sfrac);

        if (ifrac != 0) {
            out->ms = ifrac / (1000LL * 1000LL);
            ifrac -= out->ms * 1000LL * 1000LL;
            out->us = ifrac / 1000LL;
            ifrac -= out->us * 1000LL;
            out->ns = ifrac;
        } else {
            out->ms = 0;
            out->us = 0;
            out->ns = 0;
        }

        out->seconds = out->hrs * 3600 + out->min * 60 + out->sec;
        out->microseconds = out->ms * 1000 + out->us;
        out->nanoseconds = out->ns;
        break;

        default:
            PyErr_SetString(PyExc_RuntimeError,
                            "NumPy timedelta metadata is corrupted with "
                            "invalid base unit");
            return -1;
    }

    return 0;
}
