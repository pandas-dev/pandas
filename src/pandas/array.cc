// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#include "pandas/array.h"

#include "pandas/status.h"

namespace pandas {

// ----------------------------------------------------------------------
// Base array class

Status Array::Init(const TypePtr& type, size_t length) {
  type_ = type;
  length_ = length;
  return Status::OK();
}

} // namespace pandas
