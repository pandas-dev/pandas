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

#include "pandas/dispatch.h"

namespace pandas {

#define MAKE_TYPE_CASE(NAME, CapName)           \
  case NAME:                                    \
    *out = new CapName##Type();                 \
    break;

Status primitive_type_from_enum(TypeEnum tp_enum, DataType** out) {
  switch (tp_enum) {
    MAKE_TYPE_CASE(INT8, Int8);
    MAKE_TYPE_CASE(INT16, Int16);
    MAKE_TYPE_CASE(INT32, Int32);
    MAKE_TYPE_CASE(INT64, Int64);
    MAKE_TYPE_CASE(UINT8, UInt8);
    MAKE_TYPE_CASE(UINT16, UInt16);
    MAKE_TYPE_CASE(UINT32, UInt32);
    MAKE_TYPE_CASE(UINT64, UInt64);
    MAKE_TYPE_CASE(FLOAT, Float);
    MAKE_TYPE_CASE(DOUBLE, Double);
    MAKE_TYPE_CASE(BOOL, Boolean);
    MAKE_TYPE_CASE(PYOBJECT, PyObject);
    default:
      return Status::NotImplemented("Not a primitive type");
  }
  return Status::OK();
}

} // namespace pandas
