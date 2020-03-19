#ifndef PANDAS__LIBS_SRC_UJSON_DATE_CONVERSIONS
#define PANDAS__LIBS_SRC_UJSON_DATE_CONVERSIONS

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/ndarraytypes.h>
#include "datetime.h"

// Scales value inplace from nanosecond resolution to unit resolution
int scaleNanosecToUnit(npy_int64 *value, NPY_DATETIMEUNIT unit);

// Converts an int64 object representing a date to ISO format
// up to precision `base` e.g. base="s" yields 2020-01-03T00:00:00Z
// while base="ns" yields "2020-01-01T00:00:00.000000000Z"
// len is mutated to save the length of the returned string
char *int64ToIso(int64_t value, NPY_DATETIMEUNIT base, size_t *len);

// TODO: this function doesn't do a lot; should augment or replace with
// scaleNanosecToUnit
npy_datetime NpyDateTimeToEpoch(npy_datetime dt, NPY_DATETIMEUNIT base);

// Converts a Python object representing a Date / Datetime to ISO format
// up to precision `base` e.g. base="s" yields 2020-01-03T00:00:00Z
// while base="ns" yields "2020-01-01T00:00:00.000000000Z"
// len is mutated to save the length of the returned string
char *PyDateTimeToIso(PyDateTime_Date *obj, NPY_DATETIMEUNIT base, size_t *len);

// Convert a Python Date/Datetime to Unix epoch with resolution base
npy_datetime PyDateTimeToEpoch(PyDateTime_Date *dt, NPY_DATETIMEUNIT base);

char *int64ToIsoDuration(int64_t value, size_t *len);

#endif
