# -*- coding: utf-8 -*-

from cython cimport Py_ssize_t

cimport numpy as cnp
from numpy cimport ndarray, int64_t


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

cdef int64_t get_nat()
