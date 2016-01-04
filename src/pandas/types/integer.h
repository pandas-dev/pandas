// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_TYPES_INTEGER_H
#define PANDAS_TYPES_INTEGER_H

#include "pandas/array.h"
#include "pandas/numpy_interop.h"
#include "pandas/status.h"
#include "pandas/types.h"
#include "pandas/util/bitarray.h"

namespace pandas {

template <typename TypeClass>
class IntegerArrayImpl : public Array {
 public:
  typedef typename TypeClass::c_type T;

  IntegerArrayImpl() : Array() {}

  Status InitFromNumpy(PyObject* arr) {
    TypePtr type(new TypeClass());
    RETURN_NOT_OK(numpy_array_.Init(arr));
    RETURN_NOT_OK(Array::Init(type, numpy_array_.size()));
    return Status::OK();
  }

  size_t null_count() {
    return nulls_.set_count();
  }

 protected:
  NumPyBuffer numpy_array_;
  BitArray nulls_;
};

typedef IntegerArrayImpl<UInt8Type> UInt8Array;
typedef IntegerArrayImpl<Int8Type> Int8Array;

typedef IntegerArrayImpl<UInt16Type> UInt16Array;
typedef IntegerArrayImpl<Int16Type> Int16Array;

typedef IntegerArrayImpl<UInt32Type> UInt32Array;
typedef IntegerArrayImpl<Int32Type> Int32Array;

typedef IntegerArrayImpl<UInt64Type> UInt64Array;
typedef IntegerArrayImpl<Int64Type> Int64Array;

} // namespace pandas

#endif // PANDAS_TYPES_INTEGER_H
