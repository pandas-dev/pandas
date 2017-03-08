/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#ifndef PANDAS__LIBS_SRC_DATETIME_HELPER_H_
#define PANDAS__LIBS_SRC_DATETIME_HELPER_H_

#include <stdio.h>
#include "datetime.h"
#include "numpy/arrayobject.h"
#include "numpy/arrayscalars.h"

npy_int64 get_long_attr(PyObject *o, const char *attr) {
  npy_int64 long_val;
  PyObject *value = PyObject_GetAttrString(o, attr);
  long_val = (PyLong_Check(value) ?
              PyLong_AsLongLong(value) : PyInt_AS_LONG(value));
  Py_DECREF(value);
  return long_val;
}

npy_float64 total_seconds(PyObject *td) {
  // Python 2.6 compat
  npy_int64 microseconds = get_long_attr(td, "microseconds");
  npy_int64 seconds = get_long_attr(td, "seconds");
  npy_int64 days = get_long_attr(td, "days");
  npy_int64 days_in_seconds = days * 24LL * 3600LL;
  return (microseconds + (seconds + days_in_seconds) * 1000000.0) / 1000000.0;
}

#endif  // PANDAS__LIBS_SRC_DATETIME_HELPER_H_
