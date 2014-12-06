#include "datetime.h"

#if PY_MAJOR_VERSION >= 3
#define PyInt_AS_LONG PyLong_AsLong
#endif

void mangle_nat(PyObject *val) {
  PyDateTime_GET_MONTH(val) = -1;
  PyDateTime_GET_DAY(val) = -1;
}

long get_long_attr(PyObject *o, const char *attr) {
  return PyInt_AS_LONG(PyObject_GetAttrString(o, attr));
}

double total_seconds(PyObject *td) {
  // Python 2.6 compat
  long microseconds = get_long_attr(td, "microseconds");
  long seconds = get_long_attr(td, "seconds");
  long days = get_long_attr(td, "days");
  long days_in_seconds = days * 24 * 3600;
  return (microseconds + (seconds + days_in_seconds) * 1000000.0) / 1000000.0;
}
