/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#ifndef PANDAS__LIBS_SRC_NUMPY_HELPER_H_
#define PANDAS__LIBS_SRC_NUMPY_HELPER_H_

#include "Python.h"
#include "inline_helper.h"
#include "numpy/arrayobject.h"
#include "numpy/arrayscalars.h"


PANDAS_INLINE int assign_value_1d(PyArrayObject* ap, Py_ssize_t _i,
                                  PyObject* v) {
    npy_intp i = (npy_intp)_i;
    char* item = (char*)PyArray_DATA(ap) + i * PyArray_STRIDE(ap, 0);
    return PyArray_DESCR(ap)->f->setitem(v, item, ap);
}

PANDAS_INLINE PyObject* get_value_1d(PyArrayObject* ap, Py_ssize_t i) {
    char* item = (char*)PyArray_DATA(ap) + i * PyArray_STRIDE(ap, 0);
    return PyArray_Scalar(item, PyArray_DESCR(ap), (PyObject*)ap);
}

#endif  // PANDAS__LIBS_SRC_NUMPY_HELPER_H_
