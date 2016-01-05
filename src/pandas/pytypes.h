// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_PYTYPES_H
#define PANDAS_PYTYPES_H

#include <Python.h>

namespace pandas {

namespace py {

void init_natype(PyObject* na_type, PyObject* na_singleton);

extern PyObject* NAType;
extern PyObject* NA;

static inline bool is_na(PyObject* obj) {
  return (obj == NA) || PyObject_IsInstance(obj, NAType);
}

} // namespace py

} // namespace pandas

#endif // PANDAS_PYTYPES_H
