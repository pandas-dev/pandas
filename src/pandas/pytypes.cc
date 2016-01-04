// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#include <Python.h>

#include "pandas/pytypes.h"

namespace pandas {

namespace py {

void init_natype(PyObject* na_type) {
  Py_INCREF(na_type);
  NAType = na_type;
}

} // namespace py

} // namespace pandas
