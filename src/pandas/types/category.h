// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_TYPES_CATEGORY_H
#define PANDAS_TYPES_CATEGORY_H

#include <string>

#include "pandas/array.h"
#include "pandas/status.h"
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

  virtual PyObject* GetValue(size_t i);
  virtual void SetValue(size_t i, PyObject* val);

 private:
  ArrayPtr codes_;
};

} // namespace pandas

#endif // PANDAS_TYPES_CATEGORY_H
