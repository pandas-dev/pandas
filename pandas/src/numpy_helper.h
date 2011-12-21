#include "Python.h"
#include "numpy/ndarrayobject.h"

#ifndef PANDAS_INLINE
  #if defined(__GNUC__)
    #define PANDAS_INLINE __inline__
  #elif defined(_MSC_VER)
    #define PANDAS_INLINE __inline
  #elif defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
    #define PANDAS_INLINE inline
  #else
    #define PANDAS_INLINE
  #endif
#endif

PANDAS_INLINE int
is_integer_object(PyObject* obj) {
  return PyArray_IsIntegerScalar(obj);
}

PANDAS_INLINE int
is_float_object(PyObject* obj) {
  return (PyFloat_Check(obj) || PyArray_IsScalar(obj, Floating));
}

PANDAS_INLINE int
is_bool_object(PyObject* obj) {
  return (PyBool_Check(obj) || PyArray_IsScalar(obj, Bool));
}

PANDAS_INLINE int
is_string_object(PyObject* obj) {
  return (PyString_Check(obj) || PyUnicode_Check(obj));
}

PANDAS_INLINE int
assign_value_1d(PyArrayObject* ap, Py_ssize_t _i, PyObject* v) {
  npy_intp i = (npy_intp) _i;
  char *item = (char *) PyArray_DATA(ap) + i * PyArray_STRIDE(ap, 0);
  return PyArray_DESCR(ap)->f->setitem(v, item, ap);
}

PANDAS_INLINE PyObject*
get_value_1d(PyArrayObject* ap, Py_ssize_t i) {
  char *item = (char *) PyArray_DATA(ap) + i * PyArray_STRIDE(ap, 0);
  return PyArray_Scalar(item, PyArray_DESCR(ap), (PyObject*) ap);
}
