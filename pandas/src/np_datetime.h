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

/* Exported as DATETIMEUNITS in multiarraymodule.c */
static char *_datetime_strings[NPY_DATETIME_NUMUNITS] = {
    NPY_STR_Y,
    NPY_STR_M,
    NPY_STR_W,
    NPY_STR_D,
    NPY_STR_h,
    NPY_STR_m,
    NPY_STR_s,
    NPY_STR_ms,
    NPY_STR_us,
    NPY_STR_ns,
    NPY_STR_ps,
    NPY_STR_fs,
    NPY_STR_as,
    "generic"
};

static int _days_per_month_table[2][12] = {
    { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 },
    { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }
};

static int _month_offset[2][13] = {
    { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 },
    { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 }
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
