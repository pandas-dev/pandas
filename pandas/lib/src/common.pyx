include "datetime.pxi"
include "Python.pxi"

cimport numpy as np
cimport cython
from numpy cimport *

from python_dict cimport (PyDict_New, PyDict_GetItem, PyDict_SetItem,
                          PyDict_Contains, PyDict_Keys)
from python_float cimport PyFloat_Check

import numpy as np
isnan = np.isnan
cdef double NaN = <double> np.NaN

from datetime import datetime as pydatetime

cdef inline object trycall(object func, object arg):
    try:
        return func(arg)
    except:
        raise Exception('Error calling func on index %s' % arg)

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a >= b else b

ctypedef unsigned char UChar

cdef int is_contiguous(ndarray arr):
    return np.PyArray_CHKFLAGS(arr, np.NPY_C_CONTIGUOUS)

cdef int _contiguous_check(ndarray arr):
    if not is_contiguous(arr):
        raise ValueError('Tried to use data field on non-contiguous array!')

cdef int16_t *get_int16_ptr(ndarray arr):
    _contiguous_check(arr)

    return <int16_t *> arr.data

cdef int32_t *get_int32_ptr(ndarray arr):
    _contiguous_check(arr)

    return <int32_t *> arr.data

cdef double_t *get_double_ptr(ndarray arr):
    _contiguous_check(arr)

    return <double_t *> arr.data

cpdef map_indices(ndarray index):
    '''
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    '''
    cdef int i, length
    cdef flatiter iter
    cdef dict result
    cdef object idx

    result = {}

    iter = <flatiter> PyArray_IterNew(index)
    length = PyArray_SIZE(index)

    for i from 0 <= i < length:
        idx = PyArray_GETITEM(index, PyArray_ITER_DATA(iter))
        result[idx] = i
        PyArray_ITER_NEXT(iter)

    return result

cdef extern from "math.h":
    double sqrt(double x)


# initialize numpy
import_array()

