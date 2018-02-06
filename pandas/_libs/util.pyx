# -*- coding: utf-8 -*-
cimport cpython
from cpython cimport PyTypeObject

cimport numpy as cnp
from numpy cimport ndarray, NPY_C_CONTIGUOUS, NPY_F_CONTIGUOUS
cnp.import_array()


cdef object get_value_at(ndarray arr, object loc):
    cdef:
        Py_ssize_t i, sz
        int casted

    if is_float_object(loc):
        casted = int(loc)
        if casted == loc:
            loc = casted
    i = <Py_ssize_t> loc
    sz = cnp.PyArray_SIZE(arr)

    if i < 0 and sz > 0:
        i += sz
    elif i >= sz or sz == 0:
        raise IndexError('index out of bounds')

    return get_value_1d(arr, i)


cdef set_value_at_unsafe(ndarray arr, object loc, object value):
    """Sets a value into the array without checking the writeable flag.

    This should be used when setting values in a loop, check the writeable
    flag above the loop and then eschew the check on each iteration.
    """
    cdef:
        Py_ssize_t i, sz
    if is_float_object(loc):
        casted = int(loc)
        if casted == loc:
            loc = casted
    i = <Py_ssize_t> loc
    sz = cnp.PyArray_SIZE(arr)

    if i < 0:
        i += sz
    elif i >= sz:
        raise IndexError('index out of bounds')

    assign_value_1d(arr, i, value)


cdef set_value_at(ndarray arr, object loc, object value):
    """Sets a value into the array after checking that the array is mutable.
    """
    if not cnp.PyArray_ISWRITEABLE(arr):
        raise ValueError('assignment destination is read-only')

    set_value_at_unsafe(arr, loc, value)


cdef inline object unbox_if_zerodim(object arr):
    """
    If arr is zerodim array, return a proper array scalar (e.g. np.int64).
    Otherwise, return arr as is.

    Parameters
    ----------
    arr : object

    Returns
    -------
    result : object
    """
    if cnp.PyArray_IsZeroDim(arr):
        return cnp.PyArray_ToScalar(cnp.PyArray_DATA(arr), arr)
    else:
        return arr
