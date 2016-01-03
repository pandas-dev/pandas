// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_TYPES_BOOLEAN_H
#define PANDAS_TYPES_BOOLEAN_H

namespace pandas {

class BooleanArray : public NumPyArray {
 public:
  Status Init(const TypePtr& type, PyObject* codes) {
    // return Array::Init(type, length);
    return Status::NotImplemented();
  }
};

} // namespace pandas

#endif // PANDAS_TYPES_BOOLEAN_H
