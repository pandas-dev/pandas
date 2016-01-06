// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#include <Python.h>

#include "pandas/pytypes.h"

namespace pandas {

namespace py {

PyObject* NAType;
PyObject* NA;

void init_natype(PyObject* na_type, PyObject* na_singleton) {
  Py_INCREF(na_type);
  Py_INCREF(na_singleton);
  NAType = na_type;
  NA = na_singleton;
}

} // namespace py

} // namespace pandas
