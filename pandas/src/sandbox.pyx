from numpy cimport ndarray, int64_t
cimport numpy as cnp
import numpy as np

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

