/*
 * NB: This is derived from numpy 1.7 datetime.c, just enough code to
 * do what we need in numpy 1.6.
 */

#include <Python.h>
#include <datetime.h>

#include <time.h>

#include <numpy/arrayobject.h>

#include "numpy/npy_3kcompat.h"

#include "numpy/arrayscalars.h"
#include "datetime_helper.h"

/* Days per month, regular year and leap year */
NPY_NO_EXPORT int _days_per_month_table[2][12] = {
    { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 },
    { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }
};

/*
 * Returns 1 if the given year is a leap year, 0 otherwise.
 */
NPY_NO_EXPORT int
is_leapyear(npy_int64 year)
{
    return (year & 0x3) == 0 && /* year % 4 == 0 */
           ((year % 100) != 0 ||
            (year % 400) == 0);
}

/*
 * Adjusts a datetimestruct based on a minutes offset. Assumes
 * the current values are valid.
 */
NPY_NO_EXPORT void
add_minutes_to_datetimestruct(npy_datetimestruct *dts, int minutes)
{
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
        dts->day += _days_per_month_table[isleap][dts->month-1];
    }
    else if (dts->day > 28) {
        isleap = is_leapyear(dts->year);
        if (dts->day > _days_per_month_table[isleap][dts->month-1]) {
            dts->day -= _days_per_month_table[isleap][dts->month-1];
            dts->month++;
            if (dts->month > 12) {
                dts->year++;
                dts->month = 1;
            }
        }
    }
}

/*
 *
 * Tests for and converts a Python datetime.datetime or datetime.date
 * object into a NumPy npy_datetimestruct.
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
int
convert_pydatetime_to_datetimestruct(PyObject *obj, npy_datetimestruct *out,
                                     NPY_DATETIMEUNIT *out_bestunit,
                                     int apply_tzinfo)
{
    PyObject *tmp;
    int isleap;

    /* Initialize the output to all zeros */
    memset(out, 0, sizeof(npy_datetimestruct));
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
                out->day > _days_per_month_table[isleap][out->month-1]) {
        goto invalid_date;
    }

    /* Check for time attributes (if not there, return success as a date) */
    if (!PyObject_HasAttrString(obj, "hour") ||
            !PyObject_HasAttrString(obj, "minute") ||
            !PyObject_HasAttrString(obj, "second") ||
            !PyObject_HasAttrString(obj, "microsecond")) {
        /* The best unit for date is 'D' */
        if (out_bestunit != NULL) {
            *out_bestunit = NPY_FR_D;
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

    if (out->hour < 0 || out->hour >= 24 ||
            out->min < 0 || out->min >= 60 ||
            out->sec < 0 || out->sec >= 60 ||
            out->us < 0 || out->us >= 1000000) {
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
        }
        else {
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
        *out_bestunit = NPY_FR_us;
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
