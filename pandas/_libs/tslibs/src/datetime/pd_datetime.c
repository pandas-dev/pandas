/*

Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.

Copyright (c) 2005-2011, NumPy Developers
All rights reserved.

This file is derived from NumPy 1.7. See NUMPY_LICENSE.txt

*/

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "datetime.h"
#include "pd_datetime.h"

static const int days_per_month_table[2][12] = {
    {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}};

/*
 * Returns 1 if the given year is a leap year, 0 otherwise.
 */
static int is_leapyear(npy_int64 year) {
    return (year & 0x3) == 0 && /* year % 4 == 0 */
           ((year % 100) != 0 || (year % 400) == 0);
}

/*
 * Calculates the days offset from the 1970 epoch.
 */
static npy_int64 get_datetimestruct_days(const npy_datetimestruct *dts) {
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
 * Converts a datetime from a datetimestruct to a datetime based
 * on a metadata unit. The date is assumed to be valid.
 */
static npy_datetime npy_datetimestruct_to_datetime(NPY_DATETIMEUNIT base,
                                            const npy_datetimestruct *dts) {
    npy_datetime ret;

    if (base == NPY_FR_Y) {
        /* Truncate to the year */
        ret = dts->year - 1970;
    } else if (base == NPY_FR_M) {
        /* Truncate to the month */
        ret = 12 * (dts->year - 1970) + (dts->month - 1);
    } else {
        /* Otherwise calculate the number of days to start */
        npy_int64 days = get_datetimestruct_days(dts);

        switch (base) {
            case NPY_FR_W:
                /* Truncate to weeks */
                if (days >= 0) {
                    ret = days / 7;
                } else {
                    ret = (days - 6) / 7;
                }
                break;
            case NPY_FR_D:
                ret = days;
                break;
            case NPY_FR_h:
                ret = days * 24 + dts->hour;
                break;
            case NPY_FR_m:
                ret = (days * 24 + dts->hour) * 60 + dts->min;
                break;
            case NPY_FR_s:
                ret = ((days * 24 + dts->hour) * 60 + dts->min) * 60 + dts->sec;
                break;
            case NPY_FR_ms:
                ret = (((days * 24 + dts->hour) * 60 + dts->min) * 60 +
                       dts->sec) *
                          1000 +
                      dts->us / 1000;
                break;
            case NPY_FR_us:
                ret = (((days * 24 + dts->hour) * 60 + dts->min) * 60 +
                       dts->sec) *
                          1000000 +
                      dts->us;
                break;
            case NPY_FR_ns:
                ret = ((((days * 24 + dts->hour) * 60 + dts->min) * 60 +
                        dts->sec) *
                           1000000 +
                       dts->us) *
                          1000 +
                      dts->ps / 1000;
                break;
            case NPY_FR_ps:
                ret = ((((days * 24 + dts->hour) * 60 + dts->min) * 60 +
                        dts->sec) *
                           1000000 +
                       dts->us) *
                          1000000 +
                      dts->ps;
                break;
            case NPY_FR_fs:
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
            case NPY_FR_as:
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
    return ret;
}

static void
pandas_datetime_destructor(PyObject *op)
{
  void *ptr = PyCapsule_GetPointer(op, PandasDateTime_CAPSULE_NAME);
  PyMem_Free(ptr);
}

static PyMethodDef module_methods[] = {
  {NULL, NULL}
};

static struct PyModuleDef pandas_datetimemodule = {
  PyModuleDef_HEAD_INIT,
  .m_name = "_pandas_datetime",
  
  .m_doc = "Internal module with datetime support for other extensions",
  .m_size = 01,
  .m_methods = module_methods
};

PyMODINIT_FUNC
PyInit_pandas_datetime(void)
{
  PyObject *module = PyModule_Create(&pandas_datetimemodule);
  if (module == NULL)
    return NULL;

  PyDateTime_IMPORT;

  PandasDateTime_CAPI *capi = PyMem_Malloc(sizeof(PandasDateTime_CAPI));
  if (capi == NULL) {
    PyErr_NoMemory();
    return NULL;
  }
  capi->npy_datetimestruct_to_datetime = npy_datetimestruct_to_datetime;
  
  PyObject *capsule = PyCapsule_New(capi, PandasDateTime_CAPSULE_NAME, pandas_datetime_destructor);
  if (capsule == NULL) {
    PyMem_Free(capi);
    return NULL;
  }

  if (PyModule_AddObject(module, "datetime_CAPI", capsule) < 0) {
    Py_DECREF(capsule);
    // cpython doesn't PyMem_Free(capi) here - mistake or intentional?
    return NULL;
  }
    
  return module;
}
