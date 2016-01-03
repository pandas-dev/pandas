// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#include <Python.h>
#include <numpy/arrayobject.h>

#include "pandas/numpy_interop.h"

namespace pandas {

#define TYPE_MAP_CASE(NP_NAME, PD_NAME)         \
  case NPY_##NP_NAME:                           \
    *pandas_type = TypeEnum::PD_NAME;           \
    break;

Status numpy_type_num_to_pandas(int type_num, TypeEnum* pandas_type) {
  switch (type_num) {
    TYPE_MAP_CASE(INT8, INT8);
    TYPE_MAP_CASE(INT16, INT16);
    TYPE_MAP_CASE(INT32, INT32);
    TYPE_MAP_CASE(INT64, INT64);
    TYPE_MAP_CASE(UINT8, UINT8);
    TYPE_MAP_CASE(UINT16, UINT16);
    TYPE_MAP_CASE(UINT32, UINT32);
    TYPE_MAP_CASE(UINT64, UINT64);
    TYPE_MAP_CASE(FLOAT32, FLOAT);
    TYPE_MAP_CASE(FLOAT64, DOUBLE);
    TYPE_MAP_CASE(BOOL, BOOL);
    TYPE_MAP_CASE(OBJECT, PYOBJECT);
    default:
      return Status::NotImplemented("unsupported numpy type");
  }
  return Status::OK();
}

} // namespace pandas
