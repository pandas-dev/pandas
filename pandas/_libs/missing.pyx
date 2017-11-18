# -*- coding: utf-8 -*-
# cython: profile=False

from cpython cimport PyFloat_Check, PyComplex_Check

cimport cython
from cython cimport Py_ssize_t

import numpy as np
cimport numpy as np
from numpy cimport ndarray, uint8_t
np.import_array()

cimport util

from tslibs.np_datetime cimport get_timedelta64_value, get_datetime64_value
from tslibs.nattype import NaT, iNaT

cdef double INF = <double> np.inf
cdef double NEGINF = -INF


cdef inline bint is_null_datetimelike(v):
    # determine if we have a null for a timedelta/datetime (or integer
    # versions)
    if util._checknull(v):
        return True
    elif v is NaT:
        return True
    elif util.is_timedelta64_object(v):
        return v.view('int64') == iNaT
    elif util.is_datetime64_object(v):
        return v.view('int64') == iNaT
    elif util.is_integer_object(v):
        return v == iNaT
    return False


cdef inline bint _check_all_nulls(object val):
    """ utility to check if a value is any type of null """
    cdef bint res
    if PyFloat_Check(val) or PyComplex_Check(val):
        res = val != val
    elif val is NaT:
        res = 1
    elif val is None:
        res = 1
    elif util.is_datetime64_object(val):
        res = get_datetime64_value(val) == NPY_NAT
    elif util.is_timedelta64_object(val):
        res = get_timedelta64_value(val) == NPY_NAT
    else:
        res = 0
    return res


cpdef bint checknull(object val):
    if util.is_float_object(val) or util.is_complex_object(val):
        return val != val  # and val != INF and val != NEGINF
    elif util.is_datetime64_object(val):
        return get_datetime64_value(val) == NPY_NAT
    elif val is NaT:
        return True
    elif util.is_timedelta64_object(val):
        return get_timedelta64_value(val) == NPY_NAT
    elif util.is_array(val):
        return False
    else:
        return util._checknull(val)


cpdef bint checknull_old(object val):
    if util.is_float_object(val) or util.is_complex_object(val):
        return val != val or val == INF or val == NEGINF
    elif util.is_datetime64_object(val):
        return get_datetime64_value(val) == NPY_NAT
    elif val is NaT:
        return True
    elif util.is_timedelta64_object(val):
        return get_timedelta64_value(val) == NPY_NAT
    elif util.is_array(val):
        return False
    else:
        return util._checknull(val)


@cython.wraparound(False)
@cython.boundscheck(False)
def isnaobj(ndarray arr):
    cdef:
        Py_ssize_t i, n
        object val
        ndarray[uint8_t] result

    assert arr.ndim == 1, "'arr' must be 1-D."

    n = len(arr)
    result = np.empty(n, dtype=np.uint8)
    for i from 0 <= i < n:
        val = arr[i]
        result[i] = _check_all_nulls(val)
    return result.view(np.bool_)


@cython.wraparound(False)
@cython.boundscheck(False)
def isnaobj_old(ndarray arr):
    cdef:
        Py_ssize_t i, n
        object val
        ndarray[uint8_t] result

    assert arr.ndim == 1, "'arr' must be 1-D."

    n = len(arr)
    result = np.zeros(n, dtype=np.uint8)
    for i from 0 <= i < n:
        val = arr[i]
        result[i] = val is NaT or util._checknull_old(val)
    return result.view(np.bool_)


@cython.wraparound(False)
@cython.boundscheck(False)
def isnaobj2d(ndarray arr):
    cdef:
        Py_ssize_t i, j, n, m
        object val
        ndarray[uint8_t, ndim=2] result

    assert arr.ndim == 2, "'arr' must be 2-D."

    n, m = (<object> arr).shape
    result = np.zeros((n, m), dtype=np.uint8)
    for i from 0 <= i < n:
        for j from 0 <= j < m:
            val = arr[i, j]
            if checknull(val):
                result[i, j] = 1
    return result.view(np.bool_)


@cython.wraparound(False)
@cython.boundscheck(False)
def isnaobj2d_old(ndarray arr):
    cdef:
        Py_ssize_t i, j, n, m
        object val
        ndarray[uint8_t, ndim=2] result

    assert arr.ndim == 2, "'arr' must be 2-D."

    n, m = (<object> arr).shape
    result = np.zeros((n, m), dtype=np.uint8)
    for i from 0 <= i < n:
        for j from 0 <= j < m:
            val = arr[i, j]
            if checknull_old(val):
                result[i, j] = 1
    return result.view(np.bool_)


cpdef bint isposinf_scalar(object val):
    if util.is_float_object(val) and val == INF:
        return True
    else:
        return False


cpdef bint isneginf_scalar(object val):
    if util.is_float_object(val) and val == NEGINF:
        return True
    else:
        return False
