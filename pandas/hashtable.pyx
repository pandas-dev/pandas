# cython: profile=False

from cpython cimport PyObject, Py_INCREF, PyList_Check, PyTuple_Check

from khash cimport *
from numpy cimport *
from cpython cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

from util cimport _checknan
cimport util

import numpy as np
nan = np.nan

cdef extern from "numpy/npy_math.h":
    double NAN "NPY_NAN"

cimport cython
cimport numpy as cnp

cnp.import_array()
cnp.import_ufunc()

cdef int64_t iNaT = util.get_nat()
_SIZE_HINT_LIMIT = (1 << 20) + 7

cdef extern from "datetime.h":
    bint PyDateTime_Check(object o)
    void PyDateTime_IMPORT()

PyDateTime_IMPORT

cdef extern from "Python.h":
    int PySlice_Check(object)

cdef size_t _INIT_VEC_CAP = 32


include "hashtable_class_helper.pxi"
include "hashtable_func_helper.pxi"

cdef class Factorizer:
    cdef public PyObjectHashTable table
    cdef public ObjectVector uniques
    cdef public Py_ssize_t count

    def __init__(self, size_hint):
        self.table = PyObjectHashTable(size_hint)
        self.uniques = ObjectVector()
        self.count = 0

    def get_count(self):
        return self.count

    def factorize(self, ndarray[object] values, sort=False, na_sentinel=-1,
                  check_null=True):
        """
        Factorize values with nans replaced by na_sentinel
        >>> factorize(np.array([1,2,np.nan], dtype='O'), na_sentinel=20)
        array([ 0,  1, 20])
        """
        labels = self.table.get_labels(values, self.uniques,
                                       self.count, na_sentinel, check_null)
        mask = (labels == na_sentinel)
        # sort on
        if sort:
            if labels.dtype != np.intp:
                labels = labels.astype(np.intp)
            sorter = self.uniques.to_array().argsort()
            reverse_indexer = np.empty(len(sorter), dtype=np.intp)
            reverse_indexer.put(sorter, np.arange(len(sorter)))
            labels = reverse_indexer.take(labels, mode='clip')
            labels[mask] = na_sentinel
        self.count = len(self.uniques)
        return labels

    def unique(self, ndarray[object] values):
        # just for fun
        return self.table.unique(values)


cdef class Int64Factorizer:
    cdef public Int64HashTable table
    cdef public Int64Vector uniques
    cdef public Py_ssize_t count

    def __init__(self, size_hint):
        self.table = Int64HashTable(size_hint)
        self.uniques = Int64Vector()
        self.count = 0

    def get_count(self):
        return self.count

    def factorize(self, int64_t[:] values, sort=False,
                  na_sentinel=-1, check_null=True):
        labels = self.table.get_labels(values, self.uniques,
                                       self.count, na_sentinel,
                                       check_null)

        # sort on
        if sort:
            if labels.dtype != np.intp:
                labels = labels.astype(np.intp)

            sorter = self.uniques.to_array().argsort()
            reverse_indexer = np.empty(len(sorter), dtype=np.intp)
            reverse_indexer.put(sorter, np.arange(len(sorter)))

            labels = reverse_indexer.take(labels)

        self.count = len(self.uniques)
        return labels


@cython.wraparound(False)
@cython.boundscheck(False)
cdef build_count_table_object(ndarray[object] values,
                              ndarray[uint8_t, cast=True] mask,
                              kh_pymap_t *table):
    cdef:
        khiter_t k
        Py_ssize_t i, n = len(values)
        int ret = 0

    kh_resize_pymap(table, n // 10)

    for i in range(n):
        if mask[i]:
            continue

        val = values[i]
        k = kh_get_pymap(table, <PyObject*> val)
        if k != table.n_buckets:
            table.vals[k] += 1
        else:
            k = kh_put_pymap(table, <PyObject*> val, &ret)
            table.vals[k] = 1


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef value_count_object(ndarray[object] values,
                         ndarray[uint8_t, cast=True] mask):
    cdef:
        Py_ssize_t i
        kh_pymap_t *table
        int k

    table = kh_init_pymap()
    build_count_table_object(values, mask, table)

    i = 0
    result_keys = np.empty(table.n_occupied, dtype=object)
    result_counts = np.zeros(table.n_occupied, dtype=np.int64)
    for k in range(table.n_buckets):
        if kh_exist_pymap(table, k):
            result_keys[i] = <object> table.keys[k]
            result_counts[i] = table.vals[k]
            i += 1
    kh_destroy_pymap(table)

    return result_keys, result_counts


@cython.wraparound(False)
@cython.boundscheck(False)
def mode_object(ndarray[object] values, ndarray[uint8_t, cast=True] mask):
    cdef:
        int count, max_count = 2
        int j = -1 # so you can do +=
        int k
        ndarray[object] modes
        kh_pymap_t *table

    table = kh_init_pymap()
    build_count_table_object(values, mask, table)

    modes = np.empty(table.n_buckets, dtype=np.object_)
    for k in range(table.n_buckets):
        if kh_exist_pymap(table, k):
            count = table.vals[k]

            if count == max_count:
                j += 1
            elif count > max_count:
                max_count = count
                j = 0
            else:
                continue
            modes[j] = <object> table.keys[k]

    kh_destroy_pymap(table)

    return modes[:j + 1]


@cython.wraparound(False)
@cython.boundscheck(False)
def mode_int64(int64_t[:] values):
    cdef:
        int count, max_count = 2
        int j = -1 # so you can do +=
        int k
        kh_int64_t *table
        ndarray[int64_t] modes

    table = kh_init_int64()

    build_count_table_int64(values, table, 0)

    modes = np.empty(table.n_buckets, dtype=np.int64)

    with nogil:
        for k in range(table.n_buckets):
            if kh_exist_int64(table, k):
                count = table.vals[k]

                if count == max_count:
                    j += 1
                elif count > max_count:
                    max_count = count
                    j = 0
                else:
                    continue
                modes[j] = table.keys[k]

    kh_destroy_int64(table)

    return modes[:j + 1]


@cython.wraparound(False)
@cython.boundscheck(False)
def duplicated_object(ndarray[object] values, object keep='first'):
    cdef:
        Py_ssize_t i, n
        dict seen = dict()
        object row

    n = len(values)
    cdef ndarray[uint8_t] result = np.zeros(n, dtype=np.uint8)

    if keep == 'last':
        for i from n > i >= 0:
            row = values[i]
            if row in seen:
                result[i] = 1
            else:
                seen[row] = i
                result[i] = 0
    elif keep == 'first':
        for i from 0 <= i < n:
            row = values[i]
            if row in seen:
                result[i] = 1
            else:
                seen[row] = i
                result[i] = 0
    elif keep is False:
        for i from 0 <= i < n:
            row = values[i]
            if row in seen:
                result[i] = 1
                result[seen[row]] = 1
            else:
                seen[row] = i
                result[i] = 0
    else:
        raise ValueError('keep must be either "first", "last" or False')

    return result.view(np.bool_)


@cython.wraparound(False)
@cython.boundscheck(False)
def unique_label_indices(ndarray[int64_t, ndim=1] labels):
    """
    indices of the first occurrences of the unique labels
    *excluding* -1. equivelent to:
        np.unique(labels, return_index=True)[1]
    """
    cdef:
        int ret = 0
        Py_ssize_t i, n = len(labels)
        kh_int64_t * table = kh_init_int64()
        Int64Vector idx = Int64Vector()
        ndarray[int64_t, ndim=1] arr
        Int64VectorData *ud = idx.data

    kh_resize_int64(table, min(n, _SIZE_HINT_LIMIT))

    with nogil:
        for i in range(n):
            kh_put_int64(table, labels[i], &ret)
            if ret != 0:
                if needs_resize(ud):
                    with gil:
                        idx.resize()
                append_data_int64(ud, i)

    kh_destroy_int64(table)

    arr = idx.to_array()
    arr = arr[labels[arr].argsort()]

    return arr[1:] if arr.size != 0 and labels[arr[0]] == -1 else arr
