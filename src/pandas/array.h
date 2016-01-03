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

#ifndef PANDAS_ARRAY_H
#define PANDAS_ARRAY_H

#include <Python.h>

#include <memory>
#include <string>
#include <vector>

#include "pandas/status.h"
#include "pandas/types.h"
#include "pandas/util.h"

namespace pandas {

class Array {
 public:

  Array() {}
  virtual ~Array() {}

  Status Init(const TypePtr& type, size_t length) {
    type_ = type;
    length_ = length;
    return Status::OK();
  }

  size_t length() const { return length_;}
  const TypePtr& type() const { return type_;}
  TypeEnum type_enum() const { return type_->type;}

 protected:
  TypePtr type_;
  size_t length_;

 private:
  DISALLOW_COPY_AND_ASSIGN(Array);
};


typedef std::shared_ptr<Array> ArrayPtr;


// Base class for arrays whose data is stored in a NumPy array
class NumPyArray : public Array {
 public:
  explicit NumPyArray() : Array(), numpy_array_(nullptr) {}
  virtual ~NumPyArray() {
    if (numpy_array_ != nullptr) {
      Py_DECREF(numpy_array_);
    }
  }

  void Init(const TypePtr& type, size_t length, PyObject* arr) {
    Array::Init(type, length);
    Py_INCREF(arr);
    numpy_array_ = arr;
  }

 protected:
  PyObject* numpy_array_;
};

} // namespace pandas

#endif // PANDAS_ARRAY_H
