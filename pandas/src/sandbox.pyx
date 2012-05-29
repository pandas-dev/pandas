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

#----------------------------------------------------------------------
# isnull / notnull related

cdef double INF = <double> np.inf
cdef double NEGINF = -INF

cdef inline bint _checknull(object val):
    return not np.PyArray_Check(val) and (val is None or val != val)

cdef inline bint _checknan(object val):
    return not np.PyArray_Check(val) and val != val

cpdef checknull(object val):
    if util.is_float_object(val):
        return val != val or val == INF or val == NEGINF
    elif is_array(val):
        return False
    else:
        return _checknull(val)

@cython.wraparound(False)
@cython.boundscheck(False)
def isnullobj(ndarray[object] arr):
    cdef Py_ssize_t i, n
    cdef object val
    cdef ndarray[uint8_t] result

    n = len(arr)
    result = np.zeros(n, dtype=np.uint8)
    for i from 0 <= i < n:
        result[i] = _checknull(arr[i])
    return result.view(np.bool_)

@cython.wraparound(False)
@cython.boundscheck(False)
def isnullobj2d(ndarray[object, ndim=2] arr):
    cdef Py_ssize_t i, j, n, m
    cdef object val
    cdef ndarray[uint8_t, ndim=2] result

    n, m = (<object> arr).shape
    result = np.zeros((n, m), dtype=np.uint8)
    for i from 0 <= i < n:
        for j from 0 <= j < m:
            val = arr[i, j]
            if checknull(val):
                result[i, j] = 1
    return result.view(np.bool_)

from util cimport is_array

from numpy import nan

cdef extern from "math.h":
    double sqrt(double x)
    double fabs(double)

cdef float64_t FP_ERR = 1e-13

cimport util

cdef:
    int TIEBREAK_AVERAGE = 0
    int TIEBREAK_MIN = 1
    int TIEBREAK_MAX = 2
    int TIEBREAK_FIRST = 3

tiebreakers = {
    'average' : TIEBREAK_AVERAGE,
    'min' : TIEBREAK_MIN,
    'max' : TIEBREAK_MAX,
    'first' : TIEBREAK_FIRST
}

from khash cimport *

def test(ndarray arr, Py_ssize_t size_hint):
    cdef:
        kh_pymap_t *table
        int ret = 0
        khiter_t k
        PyObject **data
        Py_ssize_t i, n
        ndarray[Py_ssize_t] indexer

    table = kh_init_pymap()
    kh_resize_pymap(table, size_hint)

    data = <PyObject**> arr.data
    n = len(arr)

    indexer = np.empty(n, dtype=np.int_)

    for i in range(n):
        k = kh_put_pymap(table, data[i], &ret)

        # if not ret:
        #     kh_del_pymap(table, k)

        table.vals[k] = i

    for i in range(n):
        k = kh_get_pymap(table, data[i])
        indexer[i] = table.vals[k]

    kh_destroy_pymap(table)

    return indexer


def test_str(ndarray arr, Py_ssize_t size_hint):
    cdef:
        kh_str_t *table
        kh_cstr_t val
        int ret = 0
        khiter_t k
        PyObject **data
        Py_ssize_t i, n
        ndarray[Py_ssize_t] indexer

    table = kh_init_str()
    kh_resize_str(table, size_hint)

    data = <PyObject**> arr.data
    n = len(arr)

    indexer = np.empty(n, dtype=np.int_)

    for i in range(n):
        k = kh_put_str(table, util.get_c_string(<object> data[i]), &ret)

        # if not ret:
        #     kh_del_str(table, k)

        table.vals[k] = i

    # for i in range(n):
    #     k = kh_get_str(table, PyString_AsString(<object> data[i]))
    #     indexer[i] = table.vals[k]

    kh_destroy_str(table)

    return indexer

# def test2(ndarray[object] arr):
#     cdef:
#         dict table
#         object obj
#         Py_ssize_t i, loc, n
#         ndarray[Py_ssize_t] indexer

#     n = len(arr)
#     indexer = np.empty(n, dtype=np.int_)

#     table = {}
#     for i in range(n):
#         table[arr[i]] = i

#     for i in range(n):
#         indexer[i] =  table[arr[i]]

#     return indexer

def obj_unique(ndarray[object] arr):
    cdef:
        kh_pyset_t *table
        # PyObject *obj
        object obj
        PyObject **data
        int ret = 0
        khiter_t k
        Py_ssize_t i, n
        list uniques

    n = len(arr)
    uniques = []

    table = kh_init_pyset()

    data = <PyObject**> arr.data

    # size hint
    kh_resize_pyset(table, n // 10)

    for i in range(n):
        obj = arr[i]

        k = kh_get_pyset(table, <PyObject*> obj)
        if not kh_exist_pyset(table, k):
            k = kh_put_pyset(table, <PyObject*> obj, &ret)
            # uniques.append(obj)
            # Py_INCREF(<object> obj)

    kh_destroy_pyset(table)

    return None

def int64_unique(ndarray[int64_t] arr):
    cdef:
        kh_int64_t *table
        # PyObject *obj
        int64_t obj
        PyObject **data
        int ret = 0
        khiter_t k
        Py_ssize_t i, j, n
        ndarray[int64_t] uniques

    n = len(arr)
    uniques = np.empty(n, dtype='i8')

    table = kh_init_int64()
    kh_resize_int64(table, n)

    j = 0

    for i in range(n):
        obj = arr[i]

        k = kh_get_int64(table, obj)
        if not kh_exist_int64(table, k):
            k = kh_put_int64(table, obj, &ret)
            uniques[j] = obj
            j += 1
            # Py_INCREF(<object> obj)

    kh_destroy_int64(table)

    return np.sort(uniques[:j])


# cdef extern from "kvec.h":

#     ctypedef struct kv_int64_t:
#         size_t n, m
#         int64_t *a


def test_foo(ndarray[int64_t] values):
    cdef int64_t val

    val = values[0]
    print val

# cdef extern from "foo.h":
#     double add_things(double *a, double *b, double *c, int n)


# def cython_test(ndarray a, ndarray b, ndarray c):
#     return add_things(<double*> a.data,
#                       <double*> b.data,
#                       <double*> c.data, len(a))


# def cython_test2(ndarray[float64_t] a, ndarray[float64_t] b,
#                  ndarray[float64_t] c):
#     cdef:
#         Py_ssize_t i, n = len(a)
#         float64_t result = 0

#     for i in range(n):
#         result += a[i] + b[i] + c[i]

#     return result

@cython.boundscheck(False)
@cython.wraparound(False)
def inner(ndarray[float64_t] x, ndarray[float64_t] y):
    cdef Py_ssize_t i, n = len(x)
    cdef float64_t result = 0
    for i in range(n):
        result += x[i] * y[i]
    return result

def indices_fast(ndarray[int64_t] labels, list keys,
                 list sorted_labels):
    cdef:
        Py_ssize_t i, j, k, lab, cur, start, n = len(labels)
        dict result = {}
        object tup

    index = np.arange(n)

    k = len(keys)

    if n == 0:
        return result

    start = 0
    cur = labels[0]
    for i in range(1, n):
        lab = labels[i]

        if lab != cur:
            if lab != -1:
                tup = PyTuple_New(k)
                for j in range(k):
                    val = util.get_value_at(keys[j],
                                            sorted_labels[j][cur])
                    PyTuple_SET_ITEM(tup, j, val)
                    Py_INCREF(val)

                result[tup] = index[start:i]
            start = i
        cur = lab

    return result
