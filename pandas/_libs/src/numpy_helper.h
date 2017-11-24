/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#ifndef PANDAS__LIBS_SRC_NUMPY_HELPER_H_
#define PANDAS__LIBS_SRC_NUMPY_HELPER_H_

#include "Python.h"
#include "helper.h"
#include "numpy/arrayobject.h"
#include "numpy/arrayscalars.h"


PANDAS_INLINE npy_int64 get_nat(void) { return NPY_MIN_INT64; }

PANDAS_INLINE int is_integer_object(PyObject* obj) {
    return (!PyBool_Check(obj)) && PyArray_IsIntegerScalar(obj);
}

PANDAS_INLINE int is_float_object(PyObject* obj) {
    return (PyFloat_Check(obj) || PyArray_IsScalar(obj, Floating));
}
PANDAS_INLINE int is_complex_object(PyObject* obj) {
    return (PyComplex_Check(obj) || PyArray_IsScalar(obj, ComplexFloating));
}

PANDAS_INLINE int is_bool_object(PyObject* obj) {
    return (PyBool_Check(obj) || PyArray_IsScalar(obj, Bool));
}

PANDAS_INLINE int is_string_object(PyObject* obj) {
    return (PyString_Check(obj) || PyUnicode_Check(obj));
}

PANDAS_INLINE int is_datetime64_object(PyObject* obj) {
    return PyArray_IsScalar(obj, Datetime);
}

PANDAS_INLINE int is_timedelta64_object(PyObject* obj) {
    return PyArray_IsScalar(obj, Timedelta);
}

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

// returns ASCII or UTF8 (py3) view on python str
// python object owns memory, should not be freed
PANDAS_INLINE char* get_c_string(PyObject* obj) {
#if PY_VERSION_HEX >= 0x03000000
    return PyUnicode_AsUTF8(obj);
#else
    return PyString_AsString(obj);
#endif
}

PANDAS_INLINE PyObject* char_to_string(char* data) {
#if PY_VERSION_HEX >= 0x03000000
    return PyUnicode_FromString(data);
#else
    return PyString_FromString(data);
#endif
}

PyObject* sarr_from_data(PyArray_Descr* descr, int length, void* data) {
    PyArrayObject* result;
    npy_intp dims[1] = {length};
    Py_INCREF(descr);  // newfromdescr steals a reference to descr
    result = (PyArrayObject*)PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims,
                                                  NULL, data, 0, NULL);

    // Returned array doesn't own data by default
    result->flags |= NPY_OWNDATA;

    return (PyObject*)result;
}

void transfer_object_column(char* dst, char* src, size_t stride,
                            size_t length) {
    size_t i;
    size_t sz = sizeof(PyObject*);

    for (i = 0; i < length; ++i) {
        // uninitialized data

        // Py_XDECREF(*((PyObject**) dst));

        memcpy(dst, src, sz);
        Py_INCREF(*((PyObject**)dst));
        src += sz;
        dst += stride;
    }
}


void set_array_not_contiguous(PyArrayObject* ao) {
    ao->flags &= ~(NPY_C_CONTIGUOUS | NPY_F_CONTIGUOUS);
}

// If arr is zerodim array, return a proper array scalar (e.g. np.int64).
// Otherwise, return arr as is.
PANDAS_INLINE PyObject* unbox_if_zerodim(PyObject* arr) {
    if (PyArray_IsZeroDim(arr)) {
        PyObject* ret;
        ret = PyArray_ToScalar(PyArray_DATA(arr), arr);
        return ret;
    } else {
        Py_INCREF(arr);
        return arr;
    }
}

#endif  // PANDAS__LIBS_SRC_NUMPY_HELPER_H_
