// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_TYPES_FLOATING_H
#define PANDAS_TYPES_FLOATING_H

#include "pandas/array.h"
#include "pandas/status.h"
#include "pandas/types.h"

namespace pandas {

template <typename TypeClass>
class FloatingArrayImpl : public NumPyArray {
 public:
  typedef typename TypeClass::c_type T;

  FloatingArrayImpl() : NumPyArray() {}
  Status Init(size_t length) {
    TypePtr type(new TypeClass());
    RETURN_NOT_OK(NumPyArray::Init(type, length));
    return Status::OK();
  }
};

typedef FloatingArrayImpl<FloatType> FloatArray;
typedef FloatingArrayImpl<DoubleType> DoubleArray;

} // namespace pandas

#endif // PANDAS_TYPES_FLOATING_H
