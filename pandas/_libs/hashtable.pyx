cimport cython

from cpython.ref cimport PyObject, Py_INCREF
from cpython.mem cimport PyMem_Malloc, PyMem_Free

from libc.stdlib cimport malloc, free

import numpy as np
cimport numpy as cnp
from numpy cimport ndarray, uint8_t, uint32_t, float64_t
cnp.import_array()

cdef extern from "numpy/npy_math.h":
    float64_t NAN "NPY_NAN"

from pandas._libs.khash cimport (
    khiter_t,
    kh_str_t,
    kh_init_str,
    kh_put_str,
    kh_exist_str,
    kh_get_str,
    kh_destroy_str,
    kh_resize_str,
    kh_put_strbox,
    kh_get_strbox,
    kh_init_strbox,
    kh_int64_t,
    kh_init_int64,
    kh_resize_int64,
    kh_destroy_int64,
    kh_get_int64,
    kh_exist_int64,
    kh_put_int64,
    kh_float64_t,
    kh_exist_float64,
    kh_put_float64,
    kh_init_float64,
    kh_get_float64,
    kh_destroy_float64,
    kh_resize_float64,
    kh_resize_uint64,
    kh_exist_uint64,
    kh_destroy_uint64,
    kh_put_uint64,
    kh_get_uint64,
    kh_init_uint64,
    kh_destroy_pymap,
    kh_exist_pymap,
    kh_init_pymap,
    kh_get_pymap,
    kh_put_pymap,
    kh_resize_pymap,
)


cimport pandas._libs.util as util

from pandas._libs.missing cimport checknull


cdef int64_t NPY_NAT = util.get_nat()
_SIZE_HINT_LIMIT = (1 << 20) + 7


cdef Py_ssize_t _INIT_VEC_CAP = 128

include "hashtable_class_helper.pxi"
include "hashtable_func_helper.pxi"

cdef class Factorizer:
    cdef public:
        PyObjectHashTable table
        ObjectVector uniques
        Py_ssize_t count

    def __init__(self, size_hint):
        self.table = PyObjectHashTable(size_hint)
        self.uniques = ObjectVector()
        self.count = 0

    def get_count(self):
        return self.count

    def factorize(
        self, ndarray[object] values, sort=False, na_sentinel=-1, na_value=None
    ):
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
                                       self.count, na_sentinel, na_value)
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
    cdef public:
        Int64HashTable table
        Int64Vector uniques
        Py_ssize_t count

    def __init__(self, size_hint):
        self.table = Int64HashTable(size_hint)
        self.uniques = Int64Vector()
        self.count = 0

    def get_count(self):
        return self.count

    def factorize(self, const int64_t[:] values, sort=False,
                  na_sentinel=-1, na_value=None):
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
                                       na_value=na_value)

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
def unique_label_indices(const int64_t[:] labels):
    """
    Indices of the first occurrences of the unique labels
    *excluding* -1. equivalent to:
        np.unique(labels, return_index=True)[1]
    """
    cdef:
        int ret = 0
        Py_ssize_t i, n = len(labels)
        kh_int64_t *table = kh_init_int64()
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
    arr = arr[np.asarray(labels)[arr].argsort()]

    return arr[1:] if arr.size != 0 and labels[arr[0]] == -1 else arr
