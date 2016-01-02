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

#ifndef PANDAS_TYPES_FLOATING_H
#define PANDAS_TYPES_FLOATING_H

#include "pandas/array.h"
#include "pandas/types.h"

template <typename TypeClass>
class FloatingArrayImpl : public NumPyArray {
 public:
  typedef typename TypeClass::c_type T;

  FloatingArrayImpl() : PrimitiveArray() {}
  FloatingArrayImpl(size_t length) {
    Init(length);
  }

  void Init(size_t length) {
    TypePtr type(new TypeClass());
    NumPyArray::Init(type, length);
  }
};

typedef FloatingArrayImpl<FloatType> FloatArray;
typedef FloatingArrayImpl<DoubleType> DoubleArray;

#endif // PANDAS_TYPES_FLOATING_H
