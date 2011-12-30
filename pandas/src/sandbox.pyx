# cython: wraparound=False
# cython: boundscheck=False

from numpy cimport *
cimport numpy as cnp
import numpy as np

from cpython cimport *
cimport cpython

cnp.import_array()

cdef class SeriesIterator:

    def __init__(self, arr):
        pass

    def next(self):
        pass

def foo(object o):
    cdef int64_t bar = o
    return bar

def foo2():
    print sizeof(PyObject*)

def bench_dict():
    cdef:
        # Py_ssize_t i
        dict d = {}

    for i in range(1000000):
        d[i] = i

from cpython cimport PyObject

cdef extern from "numpy/arrayobject.h":
    bint PyArray_Check(PyObject*)

cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def bench_typecheck1(ndarray[object] arr):
    cdef Py_ssize_t i, n
    n = cnp.PyArray_SIZE(arr)
    for i in range(n):
        cpython.PyFloat_Check(arr[i])

def bench_typecheck2(ndarray[object] arr):
    cdef Py_ssize_t i, n
    cdef PyObject** buf = <PyObject**> arr.data
    n = cnp.PyArray_SIZE(arr)
    for i in range(n):
        PyArray_Check(buf[i])


def foo(object _chunk, object _arr):
    cdef:
        char* dummy_buf
        ndarray arr, result, chunk

    arr = _arr
    chunk = _chunk

    dummy_buf = chunk.data
    chunk.data = arr.data

    shape = chunk.shape
    group_size = 0
    n = len(arr)

    inc = arr.dtype.itemsize

    # chunk.shape[0] = 100
    return chunk


from skiplist cimport *

def sl_test():
    cdef int ret

    np.random.seed(12345)
    n = 100

    cdef skiplist_t* skp = skiplist_init(n)

    arr = np.random.randn(n)

    for i in range(n):
        print i
        skiplist_insert(skp, arr[i])
        # val = skiplist_get(skp, 0, &ret)
        # if ret == 0:
        #     raise ValueError('%d out of bounds' % i)

        if i >= 20:
            skiplist_remove(skp, arr[i-20])

        # skiplist_remove(skp, arr[i])
        # print 'Skiplist begin: %s' % skiplist_get(skp, 0)
        # print 'Actual begin: %s' % sorted(arr[:i+1])[0]
        data = arr[max(i-19, 0):i+1]
        print 'Skiplist middle: %s' % skiplist_get(skp, len(data) // 2, &ret)
        print 'Actual middle: %s' % sorted(data)[len(data) // 2]

    skiplist_destroy(skp)

cdef double NaN = np.NaN

def _check_minp(minp, N):
    if minp > N:
        minp = N + 1
    elif minp == 0:
        minp = 1
    elif minp < 0:
        raise ValueError('min_periods must be >= 0')
    return minp

def roll_median(ndarray[float64_t] arg, int win, int minp):
    cdef double val, res, prev
    cdef:
        int ret
        skiplist_t *sl
        Py_ssize_t midpoint, nobs = 0, i


    cdef Py_ssize_t N = len(arg)
    cdef ndarray[double_t] output = np.empty(N, dtype=float)

    sl = skiplist_init(win)

    minp = _check_minp(minp, N)

    for i from 0 <= i < minp - 1:
        val = arg[i]

        # Not NaN
        if val == val:
            nobs += 1
            skiplist_insert(sl, val)

        output[i] = NaN

    for i from minp - 1 <= i < N:
        val = arg[i]

        if i > win - 1:
            prev = arg[i - win]

            if prev == prev:
                skiplist_remove(sl, prev)
                nobs -= 1

        if val == val:
            nobs += 1
            skiplist_insert(sl, val)

        if nobs >= minp:
            midpoint = nobs / 2
            if nobs % 2:
                res = skiplist_get(sl, midpoint, &ret)
            else:
                res = (skiplist_get(sl, midpoint, &ret) +
                       skiplist_get(sl, (midpoint - 1), &ret)) / 2
        else:
            res = NaN

        output[i] = res

    skiplist_destroy(sl)

    return output

include "hashtable.pyx"

include "join.pyx"
