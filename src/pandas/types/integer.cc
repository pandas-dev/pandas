// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#include "pandas/types/integer.h"

#include <numpy/arrayobject.h>

#include "pandas/pytypes.h"

namespace pandas {

template <typename T>
PyObject* IntegerArrayImpl<T>::GetValue(size_t i) {
  if (have_nulls_ && nulls_.IsSet(i)) {
    Py_INCREF(py::NA);
    return py::NA;
  }
  return PyArray_Scalar(numpy_array_.item(i), numpy_array_.dtype(),
      reinterpret_cast<PyObject*>(numpy_array_.array()));
}

template <typename T>
void IntegerArrayImpl<T>::SetValue(size_t i, PyObject* val) {
  if (py::is_na(val)) {
    if (!have_nulls_) {
      // TODO: raise Python exception on error status
      nulls_.Init(this->length());
      have_nulls_ = true;
    }
    nulls_.Set(i);
  } else {
    if (have_nulls_) {
      nulls_.Unset(i);
    }
    // Overflow issues
    PyArrayObject* ap = numpy_array_.array();
    PyArray_DESCR((ap))->f->setitem(val, numpy_array_.item(i), ap);
  }
}


// Instantiate templates
template class IntegerArrayImpl<UInt8Type>;
template class IntegerArrayImpl<Int8Type>;
template class IntegerArrayImpl<UInt16Type>;
template class IntegerArrayImpl<Int16Type>;
template class IntegerArrayImpl<UInt32Type>;
template class IntegerArrayImpl<Int32Type>;
template class IntegerArrayImpl<UInt64Type>;
template class IntegerArrayImpl<Int64Type>;

} // namespace pandas
