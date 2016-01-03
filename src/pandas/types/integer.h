// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_TYPES_INTEGER_H
#define PANDAS_TYPES_INTEGER_H

#include "pandas/array.h"
#include "pandas/status.h"
#include "pandas/types.h"

namespace pandas {

template <typename TypeClass>
class IntegerArrayImpl : public NumPyArray {
 public:
  typedef typename TypeClass::c_type T;

  IntegerArrayImpl() : NumPyArray() {}
  Status Init(PyObject* arr) {
    TypePtr type(new TypeClass());
    RETURN_NOT_OK(NumPyArray::Init(type, length));
    return Status::OK();
  }
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
