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


#define PANDAS_DATETIME_NUMUNITS 13

#define PANDAS_DATETIME_MAX_ISO8601_STRLEN (21+3*5+1+3*6+6+1)

#define PANDAS_DATETIME_NAT NPY_MIN_INT64


typedef struct {
        npy_int64 days;
        npy_int32 hrs, min, sec, ms, us, ns, seconds, microseconds, nanoseconds;
} pandas_timedeltastruct;

typedef struct {
    NPY_DATETIMEUNIT base;
    int num;
} pandas_datetime_metadata;

typedef pandas_datetime_metadata pandas_timedelta_metadata;

extern const npy_datetimestruct _NS_MIN_DTS;
extern const npy_datetimestruct _NS_MAX_DTS;

// stuff pandas needs
// ----------------------------------------------------------------------------

int convert_pydatetime_to_datetimestruct(PyObject *obj,
                                         npy_datetimestruct *out,
                                         NPY_DATETIMEUNIT *out_bestunit,
                                         int apply_tzinfo);

npy_datetime pandas_datetimestruct_to_datetime(NPY_DATETIMEUNIT fr,
                                               npy_datetimestruct *d);

void pandas_datetime_to_datetimestruct(npy_datetime val, NPY_DATETIMEUNIT fr,
                                       npy_datetimestruct *result);

void pandas_timedelta_to_timedeltastruct(npy_timedelta val,
                                         NPY_DATETIMEUNIT fr,
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
                                   const npy_datetimestruct *dts,
                                   npy_datetime *out);

/*
 * Calculates the days offset from the 1970 epoch.
 */
npy_int64
get_datetimestruct_days(const npy_datetimestruct *dts);


/*
 * Compares two npy_datetimestruct objects chronologically
 */
int cmp_pandas_datetimestruct(const npy_datetimestruct *a,
                              const npy_datetimestruct *b);


/*
 * Adjusts a datetimestruct based on a minutes offset. Assumes
 * the current values are valid.
 */
void
add_minutes_to_datetimestruct(npy_datetimestruct *dts, int minutes);

/*
 * This provides the casting rules for the TIMEDELTA data type units.
 *
 * Notably, there is a barrier between the nonlinear years and
 * months units, and all the other units.
 */
npy_bool
can_cast_datetime64_units(NPY_DATETIMEUNIT src_unit,
                          NPY_DATETIMEUNIT dst_unit,
                          NPY_CASTING casting);


int
convert_datetime_to_datetimestruct(pandas_datetime_metadata *meta,
                                   npy_datetime dt,
                                   npy_datetimestruct *out);

int
convert_timedelta_to_timedeltastruct(pandas_timedelta_metadata *meta,
                                     npy_timedelta td,
                                     pandas_timedeltastruct *out);


NPY_DATETIMEUNIT get_datetime64_unit(PyObject *obj);


#endif  // PANDAS__LIBS_SRC_DATETIME_NP_DATETIME_H_
