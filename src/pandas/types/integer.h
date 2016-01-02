// Copyright 2015 Cloudera Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PANDAS_TYPES_INTEGER_H
#define PANDAS_TYPES_INTEGER_H

#include "pandas/array.h"
#include "pandas/types.h"


template <typename TypeClass>
class IntegerArrayImpl : public NumPyArray {
 public:
  typedef typename TypeClass::c_type T;

  IntegerArrayImpl() : PrimitiveArray() {}
  IntegerArrayImpl(size_t length) {
    Init(length);
  }

  void Init(PyObject* arr) {
    TypePtr type(new TypeClass());
    NumPyArray::Init(type, length);
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


#endif // PANDAS_TYPES_INTEGER_H
