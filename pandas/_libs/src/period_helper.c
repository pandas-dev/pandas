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


/* ------------------------------------------------------------------
 * New pandas API-helper code, to expose to cython
 * ------------------------------------------------------------------*/


// function to generate a nice string representation of the period
// object, originally from DateObject_strftime

char *c_strftime(struct date_info *tmp, char *fmt) {
    struct tm c_date;
    char *result;
    struct date_info dinfo = *tmp;
    int result_len = strlen(fmt) + 50;

    c_date.tm_sec = (int)dinfo.second;
    c_date.tm_min = dinfo.minute;
    c_date.tm_hour = dinfo.hour;
    c_date.tm_mday = dinfo.day;
    c_date.tm_mon = dinfo.month - 1;
    c_date.tm_year = dinfo.year - 1900;
    c_date.tm_wday = (dinfo.day_of_week + 1) % 7;
    c_date.tm_yday = dinfo.day_of_year - 1;
    c_date.tm_isdst = -1;

    result = malloc(result_len * sizeof(char));

    strftime(result, result_len, fmt, &c_date);

    return result;
}
