# -*- coding: utf-8 -*-
cimport cpython
from cpython cimport PyTypeObject

cdef extern from "Python.h":
    # Note: importing extern-style allows us to declare these as nogil
    # functions, whereas `from cpython cimport` does not.
    bint PyUnicode_Check(object obj) nogil
    bint PyString_Check(object obj) nogil
    bint PyBool_Check(object obj) nogil
    bint PyFloat_Check(object obj) nogil
    bint PyComplex_Check(object obj) nogil
    bint PyObject_TypeCheck(object obj, PyTypeObject* type) nogil

import numpy as np
cimport numpy as cnp
from numpy cimport ndarray, NPY_C_CONTIGUOUS, NPY_F_CONTIGUOUS
cnp.import_array()

cdef extern from "numpy/arrayobject.h":
    PyTypeObject PyFloatingArrType_Type

cdef extern from "numpy/ndarrayobject.h":
    PyTypeObject PyTimedeltaArrType_Type
    PyTypeObject PyDatetimeArrType_Type
    PyTypeObject PyComplexFloatingArrType_Type
    PyTypeObject PyBoolArrType_Type

    bint PyArray_IsIntegerScalar(obj) nogil
    bint PyArray_Check(obj) nogil


cdef int64_t get_nat():
    return np.datetime64('NaT').astype(np.int64)

# --------------------------------------------------------------------
# Type Checking

cdef bint is_string_object(object obj) nogil:
    return PyString_Check(obj) or PyUnicode_Check(obj)


cdef bint is_integer_object(object obj) nogil:
    return not PyBool_Check(obj) and PyArray_IsIntegerScalar(obj)


cdef bint is_float_object(object obj) nogil:
    return (PyFloat_Check(obj) or
            (PyObject_TypeCheck(obj, &PyFloatingArrType_Type)))


cdef bint is_complex_object(object obj) nogil:
    return (PyComplex_Check(obj) or
            PyObject_TypeCheck(obj, &PyComplexFloatingArrType_Type))


cdef bint is_bool_object(object obj) nogil:
    return (PyBool_Check(obj) or
            PyObject_TypeCheck(obj, &PyBoolArrType_Type))


cdef bint is_timedelta64_object(object obj) nogil:
    return PyObject_TypeCheck(obj, &PyTimedeltaArrType_Type)


cdef bint is_datetime64_object(object obj) nogil:
    return PyObject_TypeCheck(obj, &PyDatetimeArrType_Type)


cdef bint is_array(object o):
    return cnp.PyArray_Check(o)


cdef bint is_period_object(object val):
    return getattr(val, '_typ', '_typ') == 'period'


# --------------------------------------------------------------------


cdef inline bint _checknull(object val):
    try:
        return val is None or (cpython.PyFloat_Check(val) and val != val)
    except ValueError:
        return False


cdef inline bint _checknan(object val):
    return not cnp.PyArray_Check(val) and val != val
