/*
Copyright (c) 2025, Pandas Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#include "pandas/frame_utils.h"

/* Return whether or not the object is a local in the caller's frame.*/
int is_local_in_caller_frame_impl(PyObject *object) {
  PyFrameObject *frame = PyEval_GetFrame();
  if (frame == NULL) {
    return 0;
  }

  // Get the caller's frame (skip the current frame)
  PyFrameObject *caller_frame = PyFrame_GetBack(frame);
  if (caller_frame == NULL) {
    return 0;
  }

  // Get local variables of caller's frame and check if the object is in it
  PyObject *locals_dict = PyFrame_GetLocals(caller_frame);
  if (locals_dict == NULL) {
    Py_DECREF(caller_frame);
    return 0;
  }

  int result = 0;
  Py_ssize_t pos = 0;
  PyObject *values = PyMapping_Values(locals_dict);
  if (values == NULL) {
    Py_DECREF(locals_dict);
    Py_DECREF(caller_frame);
    return 0;
  }
  Py_ssize_t num_values = PyList_Size(values);
  for (Py_ssize_t i = 0; i < num_values; i++) {
    if (PyList_GetItem(values, i) == object) {
      result = 1;
      break;
    }
  }

  Py_DECREF(values);
  Py_DECREF(locals_dict);
  Py_DECREF(caller_frame);
  return result;
}
