# -*- coding: utf-8 -*-
# cython: profile=False

from cpython cimport PyFloat_Check, PyComplex_Check

cimport cython
from cython cimport Py_ssize_t

import numpy as np
cimport numpy as cnp
from numpy cimport ndarray, int64_t, uint8_t
cnp.import_array()

cimport util

from tslibs.np_datetime cimport get_timedelta64_value, get_datetime64_value
from tslibs.nattype import NaT

cdef double INF = <double> np.inf
cdef double NEGINF = -INF

cdef int64_t NPY_NAT = util.get_nat()


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
    """
    Return boolean describing of the input is NA-like, defined here as any
    of:
     - None
     - nan
     - NaT
     - np.datetime64 representation of NaT
     - np.timedelta64 representation of NaT

    Parameters
    ----------
    val : object

    Returns
    -------
    result : bool

    Notes
    -----
    The difference between `checknull` and `checknull_old` is that `checknull`
    does *not* consider INF or NEGINF to be NA.
    """
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
        return val is None or util.is_nan(val)


cpdef bint checknull_old(object val):
    """
    Return boolean describing of the input is NA-like, defined here as any
    of:
     - None
     - nan
     - INF
     - NEGINF
     - NaT
     - np.datetime64 representation of NaT
     - np.timedelta64 representation of NaT

    Parameters
    ----------
    val : object

    Returns
    -------
    result : bool

    Notes
    -----
    The difference between `checknull` and `checknull_old` is that `checknull`
    does *not* consider INF or NEGINF to be NA.
    """
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
        return val is None or util.is_nan(val)


cdef inline bint _check_none_nan_inf_neginf(object val):
    try:
        return val is None or (PyFloat_Check(val) and
                               (val != val or val == INF or val == NEGINF))
    except ValueError:
        return False


@cython.wraparound(False)
@cython.boundscheck(False)
def isnaobj(ndarray arr):
    """
    Return boolean mask denoting which elements of a 1-D array are na-like,
    according to the criteria defined in `_check_all_nulls`:
     - None
     - nan
     - NaT
     - np.datetime64 representation of NaT
     - np.timedelta64 representation of NaT

    Parameters
    ----------
    arr : ndarray

    Returns
    -------
    result : ndarray (dtype=np.bool_)
    """
    cdef:
        Py_ssize_t i, n
        object val
        ndarray[uint8_t] result

    assert arr.ndim == 1, "'arr' must be 1-D."

    n = len(arr)
    result = np.empty(n, dtype=np.uint8)
    for i in range(n):
        val = arr[i]
        result[i] = _check_all_nulls(val)
    return result.view(np.bool_)


@cython.wraparound(False)
@cython.boundscheck(False)
def isnaobj_old(ndarray arr):
    """
    Return boolean mask denoting which elements of a 1-D array are na-like,
    defined as being any of:
     - None
     - nan
     - INF
     - NEGINF
     - NaT

    Parameters
    ----------
    arr : ndarray

    Returns
    -------
    result : ndarray (dtype=np.bool_)
    """
    cdef:
        Py_ssize_t i, n
        object val
        ndarray[uint8_t] result

    assert arr.ndim == 1, "'arr' must be 1-D."

    n = len(arr)
    result = np.zeros(n, dtype=np.uint8)
    for i in range(n):
        val = arr[i]
        result[i] = val is NaT or _check_none_nan_inf_neginf(val)
    return result.view(np.bool_)


@cython.wraparound(False)
@cython.boundscheck(False)
def isnaobj2d(ndarray arr):
    """
    Return boolean mask denoting which elements of a 2-D array are na-like,
    according to the criteria defined in `checknull`:
     - None
     - nan
     - NaT
     - np.datetime64 representation of NaT
     - np.timedelta64 representation of NaT

    Parameters
    ----------
    arr : ndarray

    Returns
    -------
    result : ndarray (dtype=np.bool_)

    Notes
    -----
    The difference between `isnaobj2d` and `isnaobj2d_old` is that `isnaobj2d`
    does *not* consider INF or NEGINF to be NA.
    """
    cdef:
        Py_ssize_t i, j, n, m
        object val
        ndarray[uint8_t, ndim=2] result

    assert arr.ndim == 2, "'arr' must be 2-D."

    n, m = (<object> arr).shape
    result = np.zeros((n, m), dtype=np.uint8)
    for i in range(n):
        for j in range(m):
            val = arr[i, j]
            if checknull(val):
                result[i, j] = 1
    return result.view(np.bool_)


@cython.wraparound(False)
@cython.boundscheck(False)
def isnaobj2d_old(ndarray arr):
    """
    Return boolean mask denoting which elements of a 2-D array are na-like,
    according to the criteria defined in `checknull_old`:
     - None
     - nan
     - INF
     - NEGINF
     - NaT
     - np.datetime64 representation of NaT
     - np.timedelta64 representation of NaT

    Parameters
    ----------
    arr : ndarray

    Returns
    -------
    result : ndarray (dtype=np.bool_)

    Notes
    -----
    The difference between `isnaobj2d` and `isnaobj2d_old` is that `isnaobj2d`
    does *not* consider INF or NEGINF to be NA.
    """
    cdef:
        Py_ssize_t i, j, n, m
        object val
        ndarray[uint8_t, ndim=2] result

    assert arr.ndim == 2, "'arr' must be 2-D."

    n, m = (<object> arr).shape
    result = np.zeros((n, m), dtype=np.uint8)
    for i in range(n):
        for j in range(m):
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


cdef inline bint is_null_datetime64(v):
    # determine if we have a null for a datetime (or integer versions),
    # excluding np.timedelta64('nat')
    if v is None or util.is_nan(v):
        return True
    elif v is NaT:
        return True
    elif util.is_datetime64_object(v):
        return v.view('int64') == NPY_NAT
    return False


cdef inline bint is_null_timedelta64(v):
    # determine if we have a null for a timedelta (or integer versions),
    # excluding np.datetime64('nat')
    if v is None or util.is_nan(v):
        return True
    elif v is NaT:
        return True
    elif util.is_timedelta64_object(v):
        return v.view('int64') == NPY_NAT
    return False


cdef inline bint is_null_period(v):
    # determine if we have a null for a Period (or integer versions),
    # excluding np.datetime64('nat') and np.timedelta64('nat')
    if v is None or util.is_nan(v):
        return True
    elif v is NaT:
        return True
    return False
