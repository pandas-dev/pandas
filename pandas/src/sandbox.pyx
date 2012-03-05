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

# cdef extern from "numpy/arrayobject.h":
#     bint PyArray_Check(PyObject*)

cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def bench_typecheck1(ndarray[object] arr):
    cdef Py_ssize_t i, n
    n = cnp.PyArray_SIZE(arr)
    for i in range(n):
        cpython.PyFloat_Check(arr[i])

# def bench_typecheck2(ndarray[object] arr):
#     cdef Py_ssize_t i, n
#     cdef PyObject** buf = <PyObject**> arr.data
#     n = cnp.PyArray_SIZE(arr)
#     for i in range(n):
#         PyArray_Check(buf[i])


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

cdef extern from "Python.h":
    bint PyDict_Contains(object, PyObject*)
    PyObject* PyDict_GetItem(object, PyObject*)
    long PyInt_AS_LONG(PyObject*)

def get_indexer(ndarray values, dict mapping):
    cdef:
        Py_ssize_t i, length
        ndarray fill_vec
        PyObject **buf
        int32_t *resbuf
        PyObject* val

    length = len(values)
    buf = <PyObject**> values.data
    fill_vec = np.empty(length, dtype='i4')
    resbuf = <int32_t*> fill_vec.data

    for i in range(length):
        val = buf[i]
        if PyDict_Contains(mapping, val):
            resbuf[i] = PyInt_AS_LONG(PyDict_GetItem(mapping, val))
        else:
            resbuf[i] = -1
    return fill_vec


@cython.wraparound(False)
@cython.boundscheck(False)
def fancy_inc(ndarray[int64_t, ndim=2] values,
              ndarray[int64_t] iarr, ndarray[int64_t] jarr, int64_t inc):
    cdef:
        Py_ssize_t i, n = len(iarr)

    for i in range(n):
        values[iarr[i], jarr[i]] += inc



# def foo2(o):
#     return util.is_integer_object(o)

# def foo3(o):
#     return util.get_base_ndarray(o)


cimport util

from khash cimport *

cdef class Int64HashTable:

    cdef:
        kh_int64_t *table

    def __init__(self, size_hint=1):
        if size_hint is not None:
            kh_resize_int64(self.table, size_hint)

    def __cinit__(self):
        self.table = kh_init_int64()

    def __dealloc__(self):
        kh_destroy_int64(self.table)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def get_labels(self, ndarray[int64_t] values):
        cdef:
            Py_ssize_t i, n = len(values)
            ndarray[int32_t] labels
            Py_ssize_t idx, count = 0
            int ret = 0
            int64_t val
            khiter_t k

        labels = np.empty(n, dtype=np.int32)

        for i in range(n):
            val = values[i]
            k = kh_get_int64(self.table, val)
            if k != self.table.n_buckets:
                idx = self.table.vals[k]
                labels[i] = idx
            else:
                k = kh_put_int64(self.table, val, &ret)
                self.table.vals[k] = count
                labels[i] = count
                count += 1

        return labels
