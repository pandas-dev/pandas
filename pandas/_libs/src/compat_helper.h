/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#ifndef PANDAS__LIBS_SRC_COMPAT_HELPER_H_
#define PANDAS__LIBS_SRC_COMPAT_HELPER_H_

#include "Python.h"
#include "numpy_helper.h"

/*
PySlice_GetIndicesEx changes signature in PY3
but 3.6.1 in particular changes the behavior of this function slightly
https://bugs.python.org/issue27867
*/

PANDAS_INLINE int slice_get_indices(PyObject *s, Py_ssize_t length,
                                    Py_ssize_t *start, Py_ssize_t *stop, Py_ssize_t *step,
                                    Py_ssize_t *slicelength) {
#if PY_VERSION_HEX >= 0x03060000
  return PySlice_GetIndicesEx(s, length, start, stop, step, slicelength);
#else
  return PySlice_GetIndicesEx(<PySliceObject *>s, length, start, stop, step, slicelength);
#endif
}

#endif  // PANDAS__LIBS_SRC_COMPAT_HELPER_H_
