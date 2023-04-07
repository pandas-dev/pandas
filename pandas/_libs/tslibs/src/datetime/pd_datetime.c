/*

Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.

Copyright (c) 2005-2011, NumPy Developers
All rights reserved.

This file is derived from NumPy 1.7. See NUMPY_LICENSE.txt

*/

#define _PANDAS_DATETIME_IMPL

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "datetime.h"
#include "pd_datetime.h"


static void pandas_datetime_destructor(PyObject *op) {
  void *ptr = PyCapsule_GetPointer(op, PandasDateTime_CAPSULE_NAME);
  PyMem_Free(ptr);
}

static int pandas_datetime_exec(PyObject *module) {
  PyDateTime_IMPORT;
  PandasDateTime_CAPI *capi = PyMem_Malloc(sizeof(PandasDateTime_CAPI));
  if (capi == NULL) {
    PyErr_NoMemory();
    return -1;
  }
  capi->npy_datetimestruct_to_datetime = npy_datetimestruct_to_datetime;
  capi->scaleNanosecToUnit = scaleNanosecToUnit;
  capi->int64ToIso = int64ToIso;
  capi->NpyDateTimeToEpoch = NpyDateTimeToEpoch;
  capi->PyDateTimeToIso = PyDateTimeToIso;
  capi->PyDateTimeToEpoch = PyDateTimeToEpoch;
  capi->int64ToIsoDuration = int64ToIsoDuration;
  capi->pandas_datetime_to_datetimestruct = pandas_datetime_to_datetimestruct;
  capi->pandas_timedelta_to_timedeltastruct =
      pandas_timedelta_to_timedeltastruct;
  capi->convert_pydatetime_to_datetimestruct =
      convert_pydatetime_to_datetimestruct;
  capi->cmp_npy_datetimestruct = cmp_npy_datetimestruct;
  capi->get_datetime_metadata_from_dtype = get_datetime_metadata_from_dtype;
  capi->parse_iso_8601_datetime = parse_iso_8601_datetime;
  capi->get_datetime_iso_8601_strlen = get_datetime_iso_8601_strlen;
  capi->make_iso_8601_datetime = make_iso_8601_datetime;
  capi->make_iso_8601_timedelta = make_iso_8601_timedelta;

  PyObject *capsule = PyCapsule_New(capi, PandasDateTime_CAPSULE_NAME,
                                    pandas_datetime_destructor);
  if (capsule == NULL) {
    PyMem_Free(capi);
    return -1;
  }

  // Monkeypatch the top level pandas module to have an attribute for the
  // C-API. This is required because Python capsules do not support setting
  // this attribute on anything but the top level package. Ideally not
  // done when cpython gh-6898 gets implemented
  PyObject *pandas = PyImport_ImportModule("pandas");
  if (!pandas) {
    PyErr_SetString(PyExc_ImportError,
                    "pd_datetime.c could not import module pandas");
    Py_DECREF(capsule);
    return -1;
  }

  if (PyModule_AddObject(pandas, "_pandas_datetime_CAPI", capsule) < 0) {
    Py_DECREF(capsule);
    return -1;
  }

  return 0;
}

static PyModuleDef_Slot pandas_datetime_slots[] = {
    {Py_mod_exec, pandas_datetime_exec}, {0, NULL}};

static struct PyModuleDef pandas_datetimemodule = {
    PyModuleDef_HEAD_INIT,
    .m_name = "pandas._libs.pandas_datetime",

    .m_doc = "Internal module with datetime support for other extensions",
    .m_size = 0,
    .m_methods = NULL,
    .m_slots = pandas_datetime_slots};

PyMODINIT_FUNC PyInit_pandas_datetime(void) {
  return PyModuleDef_Init(&pandas_datetimemodule);
}
