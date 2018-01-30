/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.

Borrowed and derived code from scikits.timeseries that we will expose via
Cython to pandas. This primarily concerns interval representation and
frequency conversion routines.
*/

#ifndef PANDAS__LIBS_SRC_PERIOD_HELPER_H_
#define PANDAS__LIBS_SRC_PERIOD_HELPER_H_

#include <Python.h>
#include "headers/stdint.h"
#include "helper.h"
#include "limits.h"
#include "numpy/ndarraytypes.h"

/*
 * declarations from period here
 */

typedef struct date_info {
    npy_int64 absdate;
    double abstime;

    double second;
    int minute;
    int hour;
    int day;
    int month;
    int quarter;
    int year;
    int day_of_week;
    int day_of_year;
    int calendar;
} date_info;

/*
 * new pandas API helper functions here
 */

char *c_strftime(struct date_info *dinfo, char *fmt);

#endif  // PANDAS__LIBS_SRC_PERIOD_HELPER_H_
