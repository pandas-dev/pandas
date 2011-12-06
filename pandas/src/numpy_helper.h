#include "Python.h"
#include "numpy/ndarrayobject.h"

inline int
is_integer_object(PyObject* obj) {
  return PyArray_IsIntegerScalar(obj);
}

inline int
is_float_object(PyObject* obj) {
  return (PyFloat_Check(obj) || PyArray_IsScalar(obj, Floating));
}

inline int
is_bool_object(PyObject* obj) {
  return (PyBool_Check(obj) || PyArray_IsScalar(obj, Bool));
}

inline int
is_string_object(PyObject* obj) {
  return (PyString_Check(obj) || PyUnicode_Check(obj));
}

inline int
assign_value_1d(PyArrayObject* ap, Py_ssize_t _i, PyObject* v) {
  char *item;
  npy_intp i = (npy_intp) _i;
  item = PyArray_DATA(ap) + i * PyArray_STRIDE(ap, 0);
  return PyArray_DESCR(ap)->f->setitem(v, item, ap);
}
