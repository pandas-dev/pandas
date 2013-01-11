from cpython cimport PyObject, Py_INCREF, PyList_Check, PyTuple_Check

from khash cimport *
from numpy cimport *

from util cimport _checknan
cimport util

import numpy as np

ONAN = np.nan

cimport cython
cimport numpy as cnp

cnp.import_array()
cnp.import_ufunc()

cdef int64_t iNaT = util.get_nat()

cdef extern from "datetime.h":
    bint PyDateTime_Check(object o)
    void PyDateTime_IMPORT()

PyDateTime_IMPORT

cdef extern from "Python.h":
    int PySlice_Check(object)


def list_to_object_array(list obj):
    '''
    Convert list to object ndarray. Seriously can't believe I had to write this
    function
    '''
    cdef:
        Py_ssize_t i, n
        ndarray[object] arr

    n = len(obj)
    arr = np.empty(n, dtype=object)

    for i from 0 <= i < n:
        arr[i] = obj[i]

    return arr


cdef size_t _INIT_VEC_CAP = 32

cdef class ObjectVector:

    cdef:
        size_t n, m
        ndarray ao
        PyObject **data

    def __cinit__(self):
        self.n = 0
        self.m = _INIT_VEC_CAP
        self.ao = np.empty(_INIT_VEC_CAP, dtype=object)
        self.data = <PyObject**> self.ao.data

    def __len__(self):
        return self.n

    def to_array(self):
        self.ao.resize(self.n)
        return self.ao

    cdef inline append(self, object o):
        if self.n == self.m:
            self.m = self.m * 2
            self.ao.resize(self.m)
            self.data = <PyObject**> self.ao.data

        Py_INCREF(o)
        self.data[self.n] = <PyObject*> o
        self.n += 1


cdef class Int64Vector:

    cdef:
        size_t n, m
        ndarray ao
        int64_t *data

    def __cinit__(self):
        self.n = 0
        self.m = _INIT_VEC_CAP
        self.ao = np.empty(_INIT_VEC_CAP, dtype=np.int64)
        self.data = <int64_t*> self.ao.data

    def __len__(self):
        return self.n

    def to_array(self):
        self.ao.resize(self.n)
        return self.ao

    cdef inline append(self, int64_t x):
        if self.n == self.m:
            self.m = self.m * 2
            self.ao.resize(self.m)
            self.data = <int64_t*> self.ao.data

        self.data[self.n] = x
        self.n += 1

cdef class Float64Vector:

    cdef:
        size_t n, m
        ndarray ao
        float64_t *data

    def __cinit__(self):
        self.n = 0
        self.m = _INIT_VEC_CAP
        self.ao = np.empty(_INIT_VEC_CAP, dtype=np.float64)
        self.data = <float64_t*> self.ao.data

    def __len__(self):
        return self.n

    def to_array(self):
        self.ao.resize(self.n)
        return self.ao

    cdef inline append(self, float64_t x):
        if self.n == self.m:
            self.m = self.m * 2
            self.ao.resize(self.m)
            self.data = <float64_t*> self.ao.data

        self.data[self.n] = x
        self.n += 1


cdef class HashTable:
    pass


cdef class StringHashTable(HashTable):
    cdef kh_str_t *table

    # def __init__(self, size_hint=1):
    #     if size_hint is not None:
    #         kh_resize_str(self.table, size_hint)

    def __cinit__(self, int size_hint=1):
        self.table = kh_init_str()
        if size_hint is not None:
            kh_resize_str(self.table, size_hint)

    def __dealloc__(self):
        kh_destroy_str(self.table)

    cdef inline int check_type(self, object val):
        return util.is_string_object(val)

    cpdef get_item(self, object val):
        cdef khiter_t k
        k = kh_get_str(self.table, util.get_c_string(val))
        if k != self.table.n_buckets:
            return self.table.vals[k]
        else:
            raise KeyError(val)

    def get_iter_test(self, object key, Py_ssize_t iterations):
        cdef Py_ssize_t i, val
        for i in range(iterations):
            k = kh_get_str(self.table, util.get_c_string(key))
            if k != self.table.n_buckets:
                val = self.table.vals[k]

    cpdef set_item(self, object key, Py_ssize_t val):
        cdef:
            khiter_t k
            int ret = 0
            char* buf

        buf = util.get_c_string(key)

        k = kh_put_str(self.table, buf, &ret)
        self.table.keys[k] = key
        if kh_exist_str(self.table, k):
            self.table.vals[k] = val
        else:
            raise KeyError(key)

    def get_indexer(self, ndarray[object] values):
        cdef:
            Py_ssize_t i, n = len(values)
            ndarray[int64_t] labels = np.empty(n, dtype=np.int64)
            char *buf
            int64_t *resbuf = <int64_t*> labels.data
            khiter_t k
            kh_str_t *table = self.table

        for i in range(n):
            buf = util.get_c_string(values[i])
            k = kh_get_str(table, buf)
            if k != table.n_buckets:
                resbuf[i] = table.vals[k]
            else:
                resbuf[i] = -1
        return labels

    def unique(self, ndarray[object] values):
        cdef:
            Py_ssize_t i, n = len(values)
            Py_ssize_t idx, count = 0
            int ret = 0
            object val
            char *buf
            khiter_t k
            ObjectVector uniques = ObjectVector()

        for i in range(n):
            val = values[i]
            buf = util.get_c_string(val)
            k = kh_get_str(self.table, buf)
            if k == self.table.n_buckets:
                k = kh_put_str(self.table, buf, &ret)
                # print 'putting %s, %s' % (val, count)
                count += 1
                uniques.append(val)

        # return None
        return uniques.to_array()

    def factorize(self, ndarray[object] values):
        cdef:
            Py_ssize_t i, n = len(values)
            ndarray[int64_t] labels = np.empty(n, dtype=np.int64)
            dict reverse = {}
            Py_ssize_t idx, count = 0
            int ret = 0
            object val
            char *buf
            khiter_t k

        for i in range(n):
            val = values[i]
            buf = util.get_c_string(val)
            k = kh_get_str(self.table, buf)
            if k != self.table.n_buckets:
                idx = self.table.vals[k]
                labels[i] = idx
            else:
                k = kh_put_str(self.table, buf, &ret)
                # print 'putting %s, %s' % (val, count)

                self.table.vals[k] = count
                reverse[count] = val
                labels[i] = count
                count += 1

        # return None
        return reverse, labels

cdef class Int32HashTable(HashTable):
    cdef kh_int32_t *table

    def __init__(self, size_hint=1):
        if size_hint is not None:
            kh_resize_int32(self.table, size_hint)

    def __cinit__(self):
        self.table = kh_init_int32()

    def __dealloc__(self):
        kh_destroy_int32(self.table)

    cdef inline int check_type(self, object val):
        return util.is_string_object(val)

    cpdef get_item(self, int32_t val):
        cdef khiter_t k
        k = kh_get_int32(self.table, val)
        if k != self.table.n_buckets:
            return self.table.vals[k]
        else:
            raise KeyError(val)

    def get_iter_test(self, int32_t key, Py_ssize_t iterations):
        cdef Py_ssize_t i, val=0
        for i in range(iterations):
            k = kh_get_int32(self.table, val)
            if k != self.table.n_buckets:
                val = self.table.vals[k]

    cpdef set_item(self, int32_t key, Py_ssize_t val):
        cdef:
            khiter_t k
            int ret = 0

        k = kh_put_int32(self.table, key, &ret)
        self.table.keys[k] = key
        if kh_exist_int32(self.table, k):
            self.table.vals[k] = val
        else:
            raise KeyError(key)

    def map_locations(self, ndarray[int32_t] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            int32_t val
            khiter_t k

        for i in range(n):
            val = values[i]
            k = kh_put_int32(self.table, val, &ret)
            self.table.vals[k] = i

    def lookup(self, ndarray[int32_t] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            int32_t val
            khiter_t k
            ndarray[int32_t] locs = np.empty(n, dtype=np.int64)

        for i in range(n):
            val = values[i]
            k = kh_get_int32(self.table, val)
            if k != self.table.n_buckets:
                locs[i] = self.table.vals[k]
            else:
                locs[i] = -1

        return locs

    def factorize(self, ndarray[int32_t] values):
        cdef:
            Py_ssize_t i, n = len(values)
            ndarray[int64_t] labels = np.empty(n, dtype=np.int64)
            dict reverse = {}
            Py_ssize_t idx, count = 0
            int ret = 0
            int32_t val
            khiter_t k

        for i in range(n):
            val = values[i]
            k = kh_get_int32(self.table, val)
            if k != self.table.n_buckets:
                idx = self.table.vals[k]
                labels[i] = idx
            else:
                k = kh_put_int32(self.table, val, &ret)
                self.table.vals[k] = count
                reverse[count] = val
                labels[i] = count
                count += 1

        # return None
        return reverse, labels

cdef class Int64HashTable: #(HashTable):
    # cdef kh_int64_t *table

    def __cinit__(self, size_hint=1):
        self.table = kh_init_int64()
        if size_hint is not None:
            kh_resize_int64(self.table, size_hint)

    def __dealloc__(self):
        kh_destroy_int64(self.table)

    def __contains__(self, object key):
        cdef khiter_t k
        k = kh_get_int64(self.table, key)
        return k != self.table.n_buckets

    def __len__(self):
        return self.table.size

    cpdef get_item(self, int64_t val):
        cdef khiter_t k
        k = kh_get_int64(self.table, val)
        if k != self.table.n_buckets:
            return self.table.vals[k]
        else:
            raise KeyError(val)

    def get_iter_test(self, int64_t key, Py_ssize_t iterations):
        cdef Py_ssize_t i, val=0
        for i in range(iterations):
            k = kh_get_int64(self.table, val)
            if k != self.table.n_buckets:
                val = self.table.vals[k]

    cpdef set_item(self, int64_t key, Py_ssize_t val):
        cdef:
            khiter_t k
            int ret = 0

        k = kh_put_int64(self.table, key, &ret)
        self.table.keys[k] = key
        if kh_exist_int64(self.table, k):
            self.table.vals[k] = val
        else:
            raise KeyError(key)

    def map(self, ndarray[int64_t] keys, ndarray[int64_t] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            int64_t key
            khiter_t k

        for i in range(n):
            key = keys[i]
            k = kh_put_int64(self.table, key, &ret)
            self.table.vals[k] = <Py_ssize_t> values[i]

    def map_locations(self, ndarray[int64_t] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            int64_t val
            khiter_t k

        for i in range(n):
            val = values[i]
            k = kh_put_int64(self.table, val, &ret)
            self.table.vals[k] = i

    def lookup(self, ndarray[int64_t] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            int64_t val
            khiter_t k
            ndarray[int64_t] locs = np.empty(n, dtype=np.int64)

        for i in range(n):
            val = values[i]
            k = kh_get_int64(self.table, val)
            if k != self.table.n_buckets:
                locs[i] = self.table.vals[k]
            else:
                locs[i] = -1

        return locs

    def factorize(self, ndarray[object] values):
        reverse = {}
        labels = self.get_labels(values, reverse, 0)
        return reverse, labels

    def get_labels(self, ndarray[int64_t] values, Int64Vector uniques,
                   Py_ssize_t count_prior, Py_ssize_t na_sentinel):
        cdef:
            Py_ssize_t i, n = len(values)
            ndarray[int64_t] labels
            Py_ssize_t idx, count = count_prior
            int ret = 0
            int64_t val
            khiter_t k

        labels = np.empty(n, dtype=np.int64)

        for i in range(n):
            val = values[i]
            k = kh_get_int64(self.table, val)
            if k != self.table.n_buckets:
                idx = self.table.vals[k]
                labels[i] = idx
            else:
                k = kh_put_int64(self.table, val, &ret)
                self.table.vals[k] = count
                uniques.append(val)
                labels[i] = count
                count += 1

        return labels

    def get_labels_groupby(self, ndarray[int64_t] values):
        cdef:
            Py_ssize_t i, n = len(values)
            ndarray[int64_t] labels
            Py_ssize_t idx, count = 0
            int ret = 0
            int64_t val
            khiter_t k
            Int64Vector uniques = Int64Vector()

        labels = np.empty(n, dtype=np.int64)

        for i in range(n):
            val = values[i]

            # specific for groupby
            if val < 0:
                labels[i] = -1
                continue

            k = kh_get_int64(self.table, val)
            if k != self.table.n_buckets:
                idx = self.table.vals[k]
                labels[i] = idx
            else:
                k = kh_put_int64(self.table, val, &ret)
                self.table.vals[k] = count
                uniques.append(val)
                labels[i] = count
                count += 1

        arr_uniques = uniques.to_array()

        return labels, arr_uniques

    def unique(self, ndarray[int64_t] values):
        cdef:
            Py_ssize_t i, n = len(values)
            Py_ssize_t idx, count = 0
            int ret = 0
            ndarray result
            int64_t val
            khiter_t k
            Int64Vector uniques = Int64Vector()

        for i in range(n):
            val = values[i]
            k = kh_get_int64(self.table, val)
            if k == self.table.n_buckets:
                k = kh_put_int64(self.table, val, &ret)
                uniques.append(val)
                count += 1

        result = uniques.to_array()

        return result


cdef class Float64HashTable(HashTable):
    # cdef kh_float64_t *table

    def __cinit__(self, size_hint=1):
        self.table = kh_init_float64()
        if size_hint is not None:
            kh_resize_float64(self.table, size_hint)

    def __len__(self):
        return self.table.size

    def __dealloc__(self):
        kh_destroy_float64(self.table)

    def factorize(self, ndarray[float64_t] values):
        uniques = Float64Vector()
        labels = self.get_labels(values, uniques, 0, -1)
        return uniques.to_array(), labels

    def get_labels(self, ndarray[float64_t] values,
                     Float64Vector uniques,
                     Py_ssize_t count_prior, int64_t na_sentinel):
        cdef:
            Py_ssize_t i, n = len(values)
            ndarray[int64_t] labels
            Py_ssize_t idx, count = count_prior
            int ret = 0
            float64_t val
            khiter_t k

        labels = np.empty(n, dtype=np.int64)

        for i in range(n):
            val = values[i]

            if val != val:
                labels[i] = na_sentinel
                continue

            k = kh_get_float64(self.table, val)
            if k != self.table.n_buckets:
                idx = self.table.vals[k]
                labels[i] = idx
            else:
                k = kh_put_float64(self.table, val, &ret)
                self.table.vals[k] = count
                uniques.append(val)
                labels[i] = count
                count += 1

        return labels

    def map_locations(self, ndarray[float64_t] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            khiter_t k

        for i in range(n):
            k = kh_put_float64(self.table, values[i], &ret)
            self.table.vals[k] = i

    def lookup(self, ndarray[float64_t] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            float64_t val
            khiter_t k
            ndarray[int64_t] locs = np.empty(n, dtype=np.int64)

        for i in range(n):
            val = values[i]
            k = kh_get_float64(self.table, val)
            if k != self.table.n_buckets:
                locs[i] = self.table.vals[k]
            else:
                locs[i] = -1

        return locs

    def unique(self, ndarray[float64_t] values):
        cdef:
            Py_ssize_t i, n = len(values)
            Py_ssize_t idx, count = 0
            int ret = 0
            float64_t val
            khiter_t k
            Float64Vector uniques = Float64Vector()
            bint seen_na = 0

        for i in range(n):
            val = values[i]

            if val == val:
                k = kh_get_float64(self.table, val)
                if k == self.table.n_buckets:
                    k = kh_put_float64(self.table, val, &ret)
                    uniques.append(val)
                    count += 1
            elif not seen_na:
                seen_na = 1
                uniques.append(ONAN)

        return uniques.to_array()

cdef class PyObjectHashTable(HashTable):
    # cdef kh_pymap_t *table

    def __init__(self, size_hint=1):
        self.table = kh_init_pymap()
        kh_resize_pymap(self.table, size_hint)

    def __dealloc__(self):
        if self.table is not NULL:
            self.destroy()

    def __len__(self):
        return self.table.size

    def __contains__(self, object key):
        cdef khiter_t k
        hash(key)
        k = kh_get_pymap(self.table, <PyObject*>key)
        return k != self.table.n_buckets

    def destroy(self):
        kh_destroy_pymap(self.table)
        self.table = NULL

    cpdef get_item(self, object val):
        cdef khiter_t k
        k = kh_get_pymap(self.table, <PyObject*>val)
        if k != self.table.n_buckets:
            return self.table.vals[k]
        else:
            raise KeyError(val)

    def get_iter_test(self, object key, Py_ssize_t iterations):
        cdef Py_ssize_t i, val
        for i in range(iterations):
            k = kh_get_pymap(self.table, <PyObject*>key)
            if k != self.table.n_buckets:
                val = self.table.vals[k]

    cpdef set_item(self, object key, Py_ssize_t val):
        cdef:
            khiter_t k
            int ret = 0
            char* buf

        hash(key)
        k = kh_put_pymap(self.table, <PyObject*>key, &ret)
        # self.table.keys[k] = key
        if kh_exist_pymap(self.table, k):
            self.table.vals[k] = val
        else:
            raise KeyError(key)

    def map_locations(self, ndarray[object] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            object val
            khiter_t k

        for i in range(n):
            val = values[i]
            hash(val)
            k = kh_put_pymap(self.table, <PyObject*>val, &ret)
            self.table.vals[k] = i

    def lookup(self, ndarray[object] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            object val
            khiter_t k
            ndarray[int64_t] locs = np.empty(n, dtype=np.int64)

        for i in range(n):
            val = values[i]
            hash(val)
            k = kh_get_pymap(self.table, <PyObject*>val)
            if k != self.table.n_buckets:
                locs[i] = self.table.vals[k]
            else:
                locs[i] = -1

        return locs

    def lookup2(self, ndarray[object] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            object val
            khiter_t k
            long hval
            ndarray[int64_t] locs = np.empty(n, dtype=np.int64)

        # for i in range(n):
        #     val = values[i]
            # hval = PyObject_Hash(val)
            # k = kh_get_pymap(self.table, <PyObject*>val)

        return locs

    def unique(self, ndarray[object] values):
        cdef:
            Py_ssize_t i, n = len(values)
            Py_ssize_t idx, count = 0
            int ret = 0
            object val
            ndarray result
            khiter_t k
            ObjectVector uniques = ObjectVector()
            bint seen_na = 0

        for i in range(n):
            val = values[i]
            hash(val)
            if not _checknan(val):
                k = kh_get_pymap(self.table, <PyObject*>val)
                if k == self.table.n_buckets:
                    k = kh_put_pymap(self.table, <PyObject*>val, &ret)
                    uniques.append(val)
            elif not seen_na:
                seen_na = 1
                uniques.append(ONAN)

        result = uniques.to_array()

        return result

    def get_labels(self, ndarray[object] values, ObjectVector uniques,
                     Py_ssize_t count_prior, int64_t na_sentinel):
        cdef:
            Py_ssize_t i, n = len(values)
            ndarray[int64_t] labels
            Py_ssize_t idx, count = count_prior
            int ret = 0
            object val
            khiter_t k

        labels = np.empty(n, dtype=np.int64)

        for i in range(n):
            val = values[i]
            hash(val)

            if val != val or val is None:
                labels[i] = na_sentinel
                continue

            k = kh_get_pymap(self.table, <PyObject*>val)
            if k != self.table.n_buckets:
                idx = self.table.vals[k]
                labels[i] = idx
            else:
                k = kh_put_pymap(self.table, <PyObject*>val, &ret)
                self.table.vals[k] = count
                uniques.append(val)
                labels[i] = count
                count += 1

        return labels


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

    def factorize(self, ndarray[object] values, sort=False, na_sentinel=-1):
        labels = self.table.get_labels(values, self.uniques,
                                       self.count, na_sentinel)

        # sort on
        if sort:
            if labels.dtype != np.int_:
                labels = labels.astype(np.int_)

            sorter = self.uniques.to_array().argsort()
            reverse_indexer = np.empty(len(sorter), dtype=np.int_)
            reverse_indexer.put(sorter, np.arange(len(sorter)))

            labels = reverse_indexer.take(labels)

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

    def factorize(self, ndarray[int64_t] values, sort=False,
                  na_sentinel=-1):
        labels = self.table.get_labels(values, self.uniques,
                                       self.count, na_sentinel)

        # sort on
        if sort:
            if labels.dtype != np.int_:
                labels = labels.astype(np.int_)

            sorter = self.uniques.to_array().argsort()
            reverse_indexer = np.empty(len(sorter), dtype=np.int_)
            reverse_indexer.put(sorter, np.arange(len(sorter)))

            labels = reverse_indexer.take(labels)

        self.count = len(self.uniques)
        return labels



def value_count_int64(ndarray[int64_t] values):
    cdef:
        Py_ssize_t i, n = len(values)
        kh_int64_t *table
        int ret = 0
        list uniques = []

    table = kh_init_int64()
    kh_resize_int64(table, n)

    for i in range(n):
        val = values[i]
        k = kh_get_int64(table, val)
        if k != table.n_buckets:
            table.vals[k] += 1
        else:
            k = kh_put_int64(table, val, &ret)
            table.vals[k] = 1

    # for (k = kh_begin(h); k != kh_end(h); ++k)
    # 	if (kh_exist(h, k)) kh_value(h, k) = 1;
    i = 0
    result_keys = np.empty(table.n_occupied, dtype=np.int64)
    result_counts = np.zeros(table.n_occupied, dtype=np.int64)
    for k in range(table.n_buckets):
        if kh_exist_int64(table, k):
            result_keys[i] = table.keys[k]
            result_counts[i] = table.vals[k]
            i += 1
    kh_destroy_int64(table)

    return result_keys, result_counts

def value_count_object(ndarray[object] values,
                       ndarray[uint8_t, cast=True] mask):
    cdef:
        Py_ssize_t i, n = len(values)
        kh_pymap_t *table
        int ret = 0
        list uniques = []

    table = kh_init_pymap()
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

