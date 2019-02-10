# -*- coding: utf-8 -*-

import cython
from cython import Py_ssize_t

import numpy as np
cimport numpy as cnp
from numpy cimport ndarray, int64_t, uint8_t, float64_t
cnp.import_array()

cimport pandas._libs.util as util

from pandas._libs.tslibs.np_datetime cimport (
    get_timedelta64_value, get_datetime64_value)
from pandas._libs.tslibs.nattype cimport (
    checknull_with_nat, c_NaT as NaT, is_null_datetimelike)


cdef float64_t INF = <float64_t>np.inf
cdef float64_t NEGINF = -INF

cdef int64_t NPY_NAT = util.get_nat()


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
    return is_null_datetimelike(val, inat_is_null=False)


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
    if checknull(val):
        return True
    elif util.is_float_object(val) or util.is_complex_object(val):
        return val == INF or val == NEGINF
    return False


cdef inline bint _check_none_nan_inf_neginf(object val):
    try:
        return val is None or (isinstance(val, float) and
                               (val != val or val == INF or val == NEGINF))
    except ValueError:
        return False


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef ndarray[uint8_t] isnaobj(ndarray arr):
    """
    Return boolean mask denoting which elements of a 1-D array are na-like,
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
        result[i] = checknull(val)
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

    n, m = (<object>arr).shape
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

    n, m = (<object>arr).shape
    result = np.zeros((n, m), dtype=np.uint8)
    for i in range(n):
        for j in range(m):
            val = arr[i, j]
            if checknull_old(val):
                result[i, j] = 1
    return result.view(np.bool_)


def isposinf_scalar(val: object) -> bool:
    if util.is_float_object(val) and val == INF:
        return True
    else:
        return False


def isneginf_scalar(val: object) -> bool:
    if util.is_float_object(val) and val == NEGINF:
        return True
    else:
        return False


cdef inline bint is_null_datetime64(v):
    # determine if we have a null for a datetime (or integer versions),
    # excluding np.timedelta64('nat')
    if checknull_with_nat(v):
        return True
    elif util.is_datetime64_object(v):
        return get_datetime64_value(v) == NPY_NAT
    return False


cdef inline bint is_null_timedelta64(v):
    # determine if we have a null for a timedelta (or integer versions),
    # excluding np.datetime64('nat')
    if checknull_with_nat(v):
        return True
    elif util.is_timedelta64_object(v):
        return get_timedelta64_value(v) == NPY_NAT
    return False


cdef inline bint is_null_period(v):
    # determine if we have a null for a Period (or integer versions),
    # excluding np.datetime64('nat') and np.timedelta64('nat')
    return checknull_with_nat(v)
