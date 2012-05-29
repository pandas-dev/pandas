#include "datetime.h"

void mangle_nat(PyObject *val) {
  PyDateTime_GET_MONTH(val) = -1;
  PyDateTime_GET_DAY(val) = -1;
}
