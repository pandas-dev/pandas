// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#include "pandas/numpy_interop.h"

#include <numpy/arrayobject.h>

#include <memory>

#include "pandas/types/boolean.h"
#include "pandas/types/integer.h"
#include "pandas/types/floating.h"
#include "pandas/status.h"

namespace pandas {

#define TYPE_MAP_CASE(NP_NAME, PD_NAME)         \
  case NPY_##NP_NAME:                           \
    *pandas_type = TypeEnum::PD_NAME;           \
    break;

template <int NPY_TYPE>
struct numpy_traits {
};


#define NUMPY_TRAITS_DECL(NPY_TYPE, PandasArrayType)    \
  template <>                                           \
  struct numpy_traits<NPY_TYPE> {                       \
    typedef PandasArrayType ArrayType;                  \
  }


NUMPY_TRAITS_DECL(NPY_INT8, Int8Array);
NUMPY_TRAITS_DECL(NPY_INT16, Int16Array);
NUMPY_TRAITS_DECL(NPY_INT32, Int32Array);
NUMPY_TRAITS_DECL(NPY_INT64, Int64Array);
NUMPY_TRAITS_DECL(NPY_UINT8, UInt8Array);
NUMPY_TRAITS_DECL(NPY_UINT16, UInt16Array);
NUMPY_TRAITS_DECL(NPY_UINT32, UInt32Array);
NUMPY_TRAITS_DECL(NPY_UINT64, UInt64Array);
NUMPY_TRAITS_DECL(NPY_BOOL, BooleanArray);
NUMPY_TRAITS_DECL(NPY_FLOAT32, FloatArray);
NUMPY_TRAITS_DECL(NPY_FLOAT64, DoubleArray);
// NUMPY_TRAITS_DECL(NPY_OBJECT, PyObjectArray);


Status numpy_type_num_to_pandas(int type_num, TypeEnum* pandas_type) {
  switch (type_num) {
    TYPE_MAP_CASE(INT8, INT8);
    TYPE_MAP_CASE(INT16, INT16);
    TYPE_MAP_CASE(INT32, INT32);
    TYPE_MAP_CASE(INT64, INT64);
    TYPE_MAP_CASE(UINT8, UINT8);
    TYPE_MAP_CASE(UINT16, UINT16);
    TYPE_MAP_CASE(UINT32, UINT32);
    TYPE_MAP_CASE(UINT64, UINT64);
    TYPE_MAP_CASE(FLOAT32, FLOAT);
    TYPE_MAP_CASE(FLOAT64, DOUBLE);
    TYPE_MAP_CASE(BOOL, BOOL);
    TYPE_MAP_CASE(OBJECT, PYOBJECT);
    default:
      return Status::NotImplemented("unsupported numpy type");
  }
  return Status::OK();
}


template <int NPY_TYPE>
static Status convert_numpy_array(PyObject* arr, Array** out) {
  typedef typename numpy_traits<NPY_TYPE>::ArrayType ArrayType;

  ArrayType* inst = new ArrayType();
  RETURN_NOT_OK(inst->InitFromNumpy(arr));
  *out = static_cast<Array*>(inst);
  return Status::OK();
}


#define NUMPY_CONVERTER_CASE(NP_NAME, PD_NAME)                      \
  case NPY_##NP_NAME:                                               \
    RETURN_NOT_OK(convert_numpy_array<NPY_##NP_NAME>(arr, out));    \
    break;


Status array_from_numpy(PyObject* arr, Array** out) {
  int type_num = PyArray_TYPE((PyArrayObject*) arr);
  switch (type_num) {
    NUMPY_CONVERTER_CASE(INT8, Int8Array);
    NUMPY_CONVERTER_CASE(INT16, Int16Array);
    NUMPY_CONVERTER_CASE(INT32, INT32);
    NUMPY_CONVERTER_CASE(INT64, INT64);
    NUMPY_CONVERTER_CASE(UINT8, UINT8);
    NUMPY_CONVERTER_CASE(UINT16, UINT16);
    NUMPY_CONVERTER_CASE(UINT32, UINT32);
    NUMPY_CONVERTER_CASE(UINT64, UINT64);
    NUMPY_CONVERTER_CASE(FLOAT32, FLOAT);
    NUMPY_CONVERTER_CASE(FLOAT64, DOUBLE);
    NUMPY_CONVERTER_CASE(BOOL, BOOL);
    // NUMPY_CONVERTER_CASE(OBJECT, PYOBJECT);
    default:
      return Status::NotImplemented("unsupported numpy type");
  }
  return Status::OK();
}


// Convert a NumPy array to a pandas::Array with appropriate missing values set
// according to the passed uint8 dtype mask array
Status array_from_masked_numpy(PyObject* arr, PyObject* mask, Array** out) {
  return Status::NotImplemented();
}


// ----------------------------------------------------------------------
// NumPy array container

NumPyBuffer::~NumPyBuffer() {
  if (arr_ != nullptr) {
    Py_DECREF(arr_);
  }
}

Status NumPyBuffer::Init(PyObject* arr) {
  if (PyArray_NDIM((PyArrayObject*) arr) != 1) {
    return Status::Invalid("Only support 1-dimensional NumPy arrays for now");
  }
  Py_INCREF(arr);
  arr_ = arr;
  return Status::OK();
}

size_t NumPyBuffer::size() {
  return PyArray_SIZE(array());
}

int NumPyBuffer::stride() {
  return static_cast<int>(PyArray_STRIDES(array())[0]);
}

  PyObject* GetItem(size_t i);
  void PySetItem(size_t i, PyObject* val);

} // namespace pandas
