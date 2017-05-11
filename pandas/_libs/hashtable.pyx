# cython: profile=False

from cpython cimport PyObject, Py_INCREF, PyList_Check, PyTuple_Check

from khash cimport *
from numpy cimport *

from libc.stdlib cimport malloc, free
from cpython cimport (PyMem_Malloc, PyMem_Realloc, PyMem_Free,
                      PyString_Check, PyBytes_Check,
                      PyUnicode_Check)

from util cimport _checknan
cimport util

import numpy as np
nan = np.nan

cdef extern from "numpy/npy_math.h":
    double NAN "NPY_NAN"

cimport cython
cimport numpy as cnp

from pandas._libs.lib import checknull

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

cdef size_t _INIT_VEC_CAP = 128

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
        if self.uniques.external_view_exists:
            uniques = ObjectVector()
            uniques.extend(self.uniques.to_array())
            self.uniques = uniques
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
        """
        Factorize values with nans replaced by na_sentinel
        >>> factorize(np.array([1,2,np.nan], dtype='O'), na_sentinel=20)
        array([ 0,  1, 20])
        """
        if self.uniques.external_view_exists:
            uniques = Int64Vector()
            uniques.extend(self.uniques.to_array())
            self.uniques = uniques
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
