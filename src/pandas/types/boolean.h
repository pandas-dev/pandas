// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_TYPES_BOOLEAN_H
#define PANDAS_TYPES_BOOLEAN_H

#include <Python.h>

#include "pandas/array.h"
#include "pandas/numpy_interop.h"
#include "pandas/status.h"
#include "pandas/types.h"
#include "pandas/util/bitarray.h"

namespace pandas {

class BooleanArray : public Array {
 public:

  Status InitFromNumpy(PyObject* arr) {
    TypePtr type(new BooleanType());
    RETURN_NOT_OK(numpy_array_.Init(arr));
    RETURN_NOT_OK(Array::Init(type, numpy_array_.size()));
    return Status::OK();
  }

 protected:
  NumPyBuffer numpy_array_;
  BitArray nulls_;
};

} // namespace pandas

#endif // PANDAS_TYPES_BOOLEAN_H
