// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#include "pandas/dispatch.h"

#include "pandas/util/status.h"

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
