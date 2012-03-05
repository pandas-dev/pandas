/*
 * This is derived from numpy 1.7
 * See NP_LICENSE.TXT
 */

#ifndef _PANDAS_DATETIME_H_
#define _PANDAS_DATETIME_H_

#define NPY_DATETIME_MAX_ISO8601_STRLEN (21+3*5+1+3*6+6+1)

// stuff pandas needs
// ----------------------------------------------------------------------------

int convert_pydatetime_to_datetimestruct(PyObject *obj, npy_datetimestruct *out,
                                         NPY_DATETIMEUNIT *out_bestunit,
                                         int apply_tzinfo);

int dayofweek(int y, int m, int d);

static int _days_per_month_table[2][12] = {
    { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 },
    { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }
};

// stuff numpy needs in header
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
convert_datetimestruct_to_datetime(PyArray_DatetimeMetaData *meta,
                                   const npy_datetimestruct *dts,
                                   npy_datetime *out);

/*
 * Calculates the days offset from the 1970 epoch.
 */
npy_int64
get_datetimestruct_days(const npy_datetimestruct *dts);

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
//npy_bool
//can_cast_timedelta64_units(NPY_DATETIMEUNIT src_unit,
//                          NPY_DATETIMEUNIT dst_unit,
//                          NPY_CASTING casting);

npy_bool
can_cast_datetime64_units(NPY_DATETIMEUNIT src_unit,
                          NPY_DATETIMEUNIT dst_unit,
                          NPY_CASTING casting);


int
convert_datetime_to_datetimestruct(PyArray_DatetimeMetaData *meta,
                                    npy_datetime dt,
                                    npy_datetimestruct *out);

#endif
