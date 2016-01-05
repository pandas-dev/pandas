// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_TYPES_FLOATING_H
#define PANDAS_TYPES_FLOATING_H

#include "pandas/array.h"
#include "pandas/numpy_interop.h"
#include "pandas/status.h"
#include "pandas/types.h"

namespace pandas {

template <typename TypeClass>
class FloatingArrayImpl : public Array {
 public:
  typedef typename TypeClass::c_type T;

  FloatingArrayImpl() : Array() {}

  Status InitFromNumpy(PyObject* arr) {
    TypePtr type(new TypeClass());
    RETURN_NOT_OK(numpy_array_.Init(arr));
    RETURN_NOT_OK(Array::Init(type, numpy_array_.size()));
    return Status::OK();
  }

  virtual PyObject* GetValue(size_t i);
  virtual void SetValue(size_t i, PyObject* val);

 protected:
  NumPyBuffer numpy_array_;
};

typedef FloatingArrayImpl<FloatType> FloatArray;
typedef FloatingArrayImpl<DoubleType> DoubleArray;

} // namespace pandas

#endif // PANDAS_TYPES_FLOATING_H
