/*

Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.

Copyright (c) 2005-2011, NumPy Developers
All rights reserved.

This file is derived from NumPy 1.7. See NUMPY_LICENSE.txt

*/

#ifndef PANDAS__LIBS_SRC_DATETIME_NP_DATETIME_H_
#define PANDAS__LIBS_SRC_DATETIME_NP_DATETIME_H_

#include <numpy/ndarraytypes.h>

typedef enum {
        PANDAS_FR_Y = 0,  // Years
        PANDAS_FR_M = 1,  // Months
        PANDAS_FR_W = 2,  // Weeks
        // Gap where NPY_FR_B was
        PANDAS_FR_D = 4,  // Days
        PANDAS_FR_h = 5,  // hours
        PANDAS_FR_m = 6,  // minutes
        PANDAS_FR_s = 7,  // seconds
        PANDAS_FR_ms = 8,  // milliseconds
        PANDAS_FR_us = 9,  // microseconds
        PANDAS_FR_ns = 10,  // nanoseconds
        PANDAS_FR_ps = 11,  // picoseconds
        PANDAS_FR_fs = 12,  // femtoseconds
        PANDAS_FR_as = 13,  // attoseconds
        PANDAS_FR_GENERIC = 14  // Generic, unbound units, can
                                // convert to anything
} PANDAS_DATETIMEUNIT;

#define PANDAS_DATETIME_NUMUNITS 13

#define PANDAS_DATETIME_MAX_ISO8601_STRLEN (21+3*5+1+3*6+6+1)

#define PANDAS_DATETIME_NAT NPY_MIN_INT64

typedef struct {
        npy_int64 year;
        npy_int32 month, day, hour, min, sec, us, ps, as;
} pandas_datetimestruct;

typedef struct {
        npy_int64 days;
        npy_int32 hrs, min, sec, ms, us, ns, seconds, microseconds, nanoseconds;
} pandas_timedeltastruct;

typedef struct {
    PANDAS_DATETIMEUNIT base;
    int num;
} pandas_datetime_metadata;

typedef pandas_datetime_metadata pandas_timedelta_metadata;

extern const pandas_datetimestruct _NS_MIN_DTS;
extern const pandas_datetimestruct _NS_MAX_DTS;

// stuff pandas needs
// ----------------------------------------------------------------------------

int convert_pydatetime_to_datetimestruct(PyObject *obj,
                                         pandas_datetimestruct *out,
                                         PANDAS_DATETIMEUNIT *out_bestunit,
                                         int apply_tzinfo);

npy_datetime pandas_datetimestruct_to_datetime(PANDAS_DATETIMEUNIT fr,
                                               pandas_datetimestruct *d);

void pandas_datetime_to_datetimestruct(npy_datetime val, PANDAS_DATETIMEUNIT fr,
                                       pandas_datetimestruct *result);

void pandas_timedelta_to_timedeltastruct(npy_timedelta val,
                                         PANDAS_DATETIMEUNIT fr,
                                         pandas_timedeltastruct *result);

int dayofweek(int y, int m, int d);

extern const int days_per_month_table[2][12];

// stuff numpy-derived code needs in header
// ----------------------------------------------------------------------------

int is_leapyear(npy_int64 year);

/*
 * Converts a datetime from a datetimestruct to a datetime based
 * on some metadata. The date is assumed to be valid.
 *
 * TODO: If meta->num is really big, there could be overflow
 *
 * Returns 0 on success, -1 on failure.
 */
int
convert_datetimestruct_to_datetime(pandas_datetime_metadata *meta,
                                   const pandas_datetimestruct *dts,
                                   npy_datetime *out);

/*
 * Calculates the days offset from the 1970 epoch.
 */
npy_int64
get_datetimestruct_days(const pandas_datetimestruct *dts);


/*
 * Compares two pandas_datetimestruct objects chronologically
 */
int cmp_pandas_datetimestruct(const pandas_datetimestruct *a,
                              const pandas_datetimestruct *b);


/*
 * Adjusts a datetimestruct based on a minutes offset. Assumes
 * the current values are valid.
 */
void
add_minutes_to_datetimestruct(pandas_datetimestruct *dts, int minutes);


int
convert_datetime_to_datetimestruct(pandas_datetime_metadata *meta,
                                   npy_datetime dt,
                                   pandas_datetimestruct *out);

int
convert_timedelta_to_timedeltastruct(pandas_timedelta_metadata *meta,
                                     npy_timedelta td,
                                     pandas_timedeltastruct *out);


#endif  // PANDAS__LIBS_SRC_DATETIME_NP_DATETIME_H_
