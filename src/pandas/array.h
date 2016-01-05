// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_ARRAY_H
#define PANDAS_ARRAY_H

#include <Python.h>

#include <memory>
#include <string>
#include <vector>

#include "pandas/types.h"
#include "pandas/util.h"

namespace pandas {

// Forward declarations
class Status;

class Array {
 public:

  Array() {}
  virtual ~Array() {}
  Status Init(const TypePtr& type, size_t length);
  size_t length() const { return length_;}
  const TypePtr& type() const { return type_;}
  TypeEnum type_enum() const { return type_->type;}

  // Dynamic dispatch array API
  virtual PyObject* GetValue(size_t i) = 0;
  virtual void SetValue(size_t i, PyObject* val) = 0;

 protected:
  TypePtr type_;
  size_t length_;

 private:
  DISALLOW_COPY_AND_ASSIGN(Array);
};


typedef std::shared_ptr<Array> ArrayPtr;

} // namespace pandas

#endif // PANDAS_ARRAY_H
