/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#ifndef PANDAS__LIBS_SRC_COMPAT_HELPER_H_
#define PANDAS__LIBS_SRC_COMPAT_HELPER_H_

#include "Python.h"
#include "inline_helper.h"

/*
PySlice_GetIndicesEx changes signature in PY3
but 3.6.1 in particular changes the behavior of this function slightly
https://bugs.python.org/issue27867


In 3.6.1 PySlice_GetIndicesEx was changed to a macro
inadvertently breaking ABI compat.  For now, undefing
the macro, which restores compat.
https://github.com/pandas-dev/pandas/issues/15961
https://bugs.python.org/issue29943
*/

#ifndef PYPY_VERSION
# if PY_VERSION_HEX < 0x03070000 && defined(PySlice_GetIndicesEx)
#   undef PySlice_GetIndicesEx
# endif
#endif

PANDAS_INLINE int slice_get_indices(PyObject *s,
                                    Py_ssize_t length,
                                    Py_ssize_t *start,
                                    Py_ssize_t *stop,
                                    Py_ssize_t *step,
                                    Py_ssize_t *slicelength) {
#if PY_VERSION_HEX >= 0x03000000
  return PySlice_GetIndicesEx(s, length, start, stop,
                              step, slicelength);
#else
  return PySlice_GetIndicesEx((PySliceObject *)s, length, start,
                              stop, step, slicelength);
#endif
}

#endif  // PANDAS__LIBS_SRC_COMPAT_HELPER_H_
