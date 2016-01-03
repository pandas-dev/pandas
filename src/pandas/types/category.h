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

#ifndef PANDAS_TYPES_CATEGORY_H
#define PANDAS_TYPES_CATEGORY_H

#include <string>

#include "pandas/array.h"
#include "pandas/types.h"

namespace pandas {

struct CategoryType : public DataType {

  explicit CategoryType(const ArrayPtr& categories)
      : DataType(TypeEnum::CATEGORY) {
    categories_ = categories;
  }

  virtual std::string ToString() {
    std::stringstream s;
    s << "category<" << category_type()->ToString() << ">";
    return s.str();
  }

  TypePtr category_type() {
    return categories_->type();
  }

  ArrayPtr categories() {
    return categories_;
  }

 protected:
  ArrayPtr categories_;
};


class CategoryArray : public Array {
 public:
  Status Init(const TypePtr& type, PyObject* codes) {
    // return Array::Init(type, length);
    return Status::NotImplemented();
  }

  ArrayPtr codes() {
    return codes_;
  }

  ArrayPtr categories() {
    return static_cast<CategoryType*>(type_.get())->categories();
  }

 private:
  ArrayPtr codes_;
};

} // namespace pandas

#endif // PANDAS_TYPES_CATEGORY_H
