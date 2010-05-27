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

cdef int64_t *get_int64_ptr(ndarray arr):
    _contiguous_check(arr)

    return <int64_t *> arr.data

cdef double_t *get_double_ptr(ndarray arr):
    _contiguous_check(arr)

    return <double_t *> arr.data

cdef extern from "math.h":
    double sqrt(double x)

cdef extern from "cobject.h":
    pass # for datetime API

cdef extern from "datetime.h":

    ctypedef class datetime.datetime [object PyDateTime_DateTime]:
        # cdef int *data
        # cdef long hashcode
        # cdef char hastzinfo
        pass

    int PyDateTime_GET_YEAR(datetime o)
    int PyDateTime_GET_MONTH(datetime o)
    int PyDateTime_GET_DAY(datetime o)
    int PyDateTime_DATE_GET_HOUR(datetime o)
    int PyDateTime_DATE_GET_MINUTE(datetime o)
    int PyDateTime_DATE_GET_SECOND(datetime o)
    int PyDateTime_DATE_GET_MICROSECOND(datetime o)
    int PyDateTime_TIME_GET_HOUR(datetime o)
    int PyDateTime_TIME_GET_MINUTE(datetime o)
    int PyDateTime_TIME_GET_SECOND(datetime o)
    int PyDateTime_TIME_GET_MICROSECOND(datetime o)
    bint PyDateTime_Check(object o)
    void PyDateTime_IMPORT()

# import datetime C API
PyDateTime_IMPORT

# initialize numpy
import_array()

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

def isAllDates(ndarray index):
    cdef int i, length
    cdef flatiter iter
    cdef object date

    iter = <flatiter> PyArray_IterNew(index)
    length = PyArray_SIZE(index)

    if length == 0:
        return False

    for i from 0 <= i < length:
        date = PyArray_GETITEM(index, PyArray_ITER_DATA(iter))

        if not PyDateTime_Check(date):
            return False
        PyArray_ITER_NEXT(iter)

    return True

def isAllDates2(ndarray[object, ndim=1] arr):
    '''
    cannot use
    '''

    cdef int i, size = len(arr)
    cdef object date

    if size == 0:
        return False

    for i from 0 <= i < size:
        date = arr[i]

        if not PyDateTime_Check(date):
            return False

    return True
