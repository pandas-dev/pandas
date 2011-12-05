#include "numpy/ndarrayobject.h"

inline int
is_integer_object(PyObject* obj) {
  return PyArray_IsIntegerScalar(obj);
}

inline int
is_float_object(PyObject* obj) {
  return (PyFloat_Check(obj)
          || PyObject_TypeCheck(obj, &PyFloatingArrType_Type));
}
