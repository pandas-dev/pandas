# -*- coding: utf-8 -*-

from cython cimport Py_ssize_t

cimport numpy as cnp
from numpy cimport ndarray

cdef extern from "src/headers/stdint.h":
    enum: UINT8_MAX
    enum: UINT16_MAX
    enum: UINT32_MAX
    enum: UINT64_MAX
    enum: INT8_MIN
    enum: INT8_MAX
    enum: INT16_MIN
    enum: INT16_MAX
    enum: INT32_MAX
    enum: INT32_MIN
    enum: INT64_MAX
    enum: INT64_MIN

ctypedef fused numeric:
    cnp.int8_t
    cnp.int16_t
    cnp.int32_t
    cnp.int64_t

    cnp.uint8_t
    cnp.uint16_t
    cnp.uint32_t
    cnp.uint64_t

    cnp.float32_t
    cnp.float64_t


cdef extern from "src/numpy_helper.h":
    void set_array_not_contiguous(ndarray ao)

    int assign_value_1d(ndarray, Py_ssize_t, object) except -1
    cnp.int64_t get_nat()
    object get_value_1d(ndarray, Py_ssize_t)
    char *get_c_string(object) except NULL
    object char_to_string(char*)


cdef bint is_string_object(object obj) nogil
cdef bint is_integer_object(object obj) nogil
cdef bint is_float_object(object obj) nogil
cdef bint is_complex_object(object obj) nogil
cdef bint is_bool_object(object obj) nogil
cdef bint is_timedelta64_object(object obj) nogil
cdef bint is_datetime64_object(object obj) nogil

cdef bint is_array(object o)
cdef bint is_period_object(object val)

cdef bint _checknull(object val)
cdef bint _checknan(object val)

cdef object unbox_if_zerodim(object arr)

cdef set_value_at(ndarray arr, object loc, object value)
cdef set_value_at_unsafe(ndarray arr, object loc, object value)
cdef object get_value_at(ndarray arr, object loc)
