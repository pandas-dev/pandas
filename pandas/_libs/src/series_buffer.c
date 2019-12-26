#include "Python.h"

static int series_getbuffer(PyObject *obj, Py_buffer *view, int flags) {
  PyObject *arr = PyObject_CallMethod(obj, "to_numpy", NULL);
  if (arr == NULL) {
    return NULL;
  }
  
  return PyObject_GetBuffer(arr, view, flags);
}

static int series_releasebuffer(PyObject *obj, Py_buffer, *view) {
  PyObject *arr = PyObject_CallMethod(obj, "to_numpy", NULL);
  if (arr == NULL) {
    return NULL;
  }
  
  return PyObject_ReleaseBuffer(arr, view);
}
