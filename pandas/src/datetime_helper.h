#include "datetime.h"
#include "numpy/arrayobject.h"
#include "numpy/arrayscalars.h"
#include <stdio.h>

#if PY_MAJOR_VERSION >= 3
#define PyInt_AS_LONG PyLong_AsLong
#endif

void mangle_nat(PyObject *val) {
  PyDateTime_GET_MONTH(val) = -1;
  PyDateTime_GET_DAY(val) = -1;
}

npy_int64 get_long_attr(PyObject *o, const char *attr) {
  PyObject *value = PyObject_GetAttrString(o, attr);
  return PyLong_Check(value) ? PyLong_AsLongLong(value) : PyInt_AS_LONG(value);
}

npy_float64 total_seconds(PyObject *td) {
  // Python 2.6 compat
  npy_int64 microseconds = get_long_attr(td, "microseconds");
  npy_int64 seconds = get_long_attr(td, "seconds");
  npy_int64 days = get_long_attr(td, "days");
  npy_int64 days_in_seconds = days * 24LL * 3600LL;
  return (microseconds + (seconds + days_in_seconds) * 1000000.0) / 1000000.0;
}
