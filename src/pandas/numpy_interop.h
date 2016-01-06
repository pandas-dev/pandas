// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_NUMPY_INTEROP_H
#define PANDAS_NUMPY_INTEROP_H

#include <Python.h>

#include <numpy/numpyconfig.h>

// Don't use the deprecated Numpy functions
#ifdef NPY_1_7_API_VERSION
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#else
#define NPY_ARRAY_NOTSWAPPED NPY_NOTSWAPPED
#define NPY_ARRAY_ALIGNED NPY_ALIGNED
#define NPY_ARRAY_WRITEABLE NPY_WRITEABLE
#define NPY_ARRAY_UPDATEIFCOPY NPY_UPDATEIFCOPY
#endif

// This is required to be able to access the NumPy C API properly in C++ files
// other than this main one
#define PY_ARRAY_UNIQUE_SYMBOL pandas_ARRAY_API
#ifndef NUMPY_IMPORT_ARRAY
#define NO_IMPORT_ARRAY
#endif

#include <numpy/arrayobject.h>
#include <numpy/ufuncobject.h>

#include "pandas/array.h"
#include "pandas/types.h"

namespace pandas {

inline int import_numpy() {
#ifdef NUMPY_IMPORT_ARRAY
  import_array1(-1);
  import_umath1(-1);
#endif

  return 0;
}

class Status;

Status numpy_type_num_to_pandas(int type_num, TypeEnum* pandas_type);

Status array_from_numpy(PyObject* arr, Array** out);
Status array_from_masked_numpy(PyObject* arr, PyObject* mask, Array** out);


// Container for strided (but contiguous) data contained in a NumPy array
class NumPyBuffer {
 public:
  NumPyBuffer() : arr_(nullptr) {}
  virtual ~NumPyBuffer();
  Status Init(PyObject* arr);

  size_t size();
  int stride();

  PyArrayObject* array() {
    return reinterpret_cast<PyArrayObject*>(arr_);
  }

  PyArray_Descr* dtype() {
    return PyArray_DESCR(array());
  }

  char* item(size_t i) {
    char* data = reinterpret_cast<char*>(PyArray_DATA(array()));
    return data + i * (PyArray_STRIDE(array(), 0));
  }

  PyObject* GetItem(size_t i);
  void PySetItem(size_t i, PyObject* val);

 protected:
  PyObject* arr_;
};

} // namespace pandas

#endif // PANDAS_NUMPY_INTEROP_H
