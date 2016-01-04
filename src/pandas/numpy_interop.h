// This file is a part of pandas. See LICENSE for details about reuse and
// copyright holders

#ifndef PANDAS_NUMPY_INTEROP_H
#define PANDAS_NUMPY_INTEROP_H

#include <Python.h>

#include "pandas/array.h"
#include "pandas/types.h"

namespace pandas {

class Status;

Status numpy_type_num_to_pandas(int type_num, TypeEnum* pandas_type);

Status array_from_numpy(PyObject* arr, Array** out);
Status array_from_masked_numpy(PyObject* arr, PyObject* mask, Array** out);


// Container for strided (but contiguous) data contained in a NumPy array
class NumPyBuffer {
 public:
  explicit NumPyBuffer() : arr_(nullptr) {}
  virtual ~NumPyBuffer();
  Status Init(PyObject* arr);

  size_t size();
  int stride();

 protected:
  PyObject* arr_;
};

void init_numpy();

} // namespace pandas

#endif // PANDAS_NUMPY_INTEROP_H
