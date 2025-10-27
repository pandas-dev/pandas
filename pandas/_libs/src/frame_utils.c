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
  PyObject *key, *value;
  Py_ssize_t pos = 0;
  while (PyDict_Next(locals_dict, &pos, &key, &value)) {
    if (object == value) {
      result = 1;
      break;
    }
  }

  Py_DECREF(locals_dict);
  Py_DECREF(caller_frame);
  return result;
}
