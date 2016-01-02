// copyright 2015 Cloudera Inc.
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

#ifndef PANDAS_TYPES_CATEGORY_H
#define PANDAS_TYPES_CATEGORY_H

#include <memory>

#include "pandas/array.h"
#include "pandas/types.h"

class CategoryArray : public Array {
 public:
  void Init(PyObject* arr) {
    TypePtr type(new TypeClass());
    NumPyArray::Init(type, length);
  }

  std::shared_ptr<Array> codes() {
    return codes_;
  }

  std::shared_ptr<Array> categories() {
    return categories_;
  }

 private:
  std::shared_ptr<Array> codes_;
  std::shared_ptr<Array> categories_;
};

#endif // PANDAS_TYPES_CATEGORY_H
