#include "numpy/ndarrayobject.h"

inline int
is_integer_object(PyObject* obj) {
  return PyArray_IsIntegerScalar(obj);
}
