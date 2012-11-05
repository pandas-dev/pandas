from cpython cimport PyObject

from khash cimport *
from numpy cimport *

from util cimport _checknan
cimport util

import numpy as np

ONAN = np.nan


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


cdef extern from "kvec.h":

    ctypedef struct kv_int64_t:
        size_t n, m
        int64_t* a

    ctypedef struct kv_object_t:
        size_t n, m
        PyObject** a

    inline void kv_object_push(kv_object_t *v, PyObject* x)
    inline void kv_object_destroy(kv_object_t *v)
    inline void kv_int64_push(kv_int64_t *v, int64_t x)


cdef class ObjectVector:

    cdef:
        bint owndata
        kv_object_t vec

    def __cinit__(self):
        self.owndata = 1

    def to_array(self, xfer_data=True):
        """ Here we use the __array__ method, that is called when numpy
            tries to get an array from the object."""
        cdef:
            npy_intp shape[1]
            ndarray result

        shape[0] = <npy_intp> self.vec.n

        # Create a 1D array, of length 'size'
        result = PyArray_SimpleNewFromData(1, shape,
                                           np.NPY_OBJECT, self.vec.a)

        # urgh, mingw32 barfs because of this

        # if xfer_data:
        #     self.owndata = 0
        #     util.set_array_owndata(result)
        # return result

        return result.copy()

    cdef inline append(self, object o):
        kv_object_push(&self.vec, <PyObject*> o)

    def __dealloc__(self):
        if self.owndata:
            kv_object_destroy(&self.vec)


cdef class Int64Vector:

    cdef:
        bint owndata
        kv_int64_t vec

    def __cinit__(self):
        self.owndata = 1

    def to_array(self, xfer_data=True):
        """ Here we use the __array__ method, that is called when numpy
            tries to get an array from the object."""
        cdef:
            npy_intp shape[1]
            ndarray result

        shape[0] = <npy_intp> self.vec.n

        # Create a 1D array, of length 'size'
        result = PyArray_SimpleNewFromData(1, shape, np.NPY_INT64,
                                           self.vec.a)

        if xfer_data:
            self.owndata = 0
            util.set_array_owndata(result)

        return result

    cdef inline append(self, int64_t x):
        kv_int64_push(&self.vec, x)

    def __dealloc__(self):
        if self.owndata:
            free(self.vec.a)


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
            list uniques = []

        for i in range(n):
            val = values[i]
            buf = util.get_c_string(val)
            k = kh_get_str(self.table, buf)
            if k == self.table.n_buckets:
                k = kh_put_str(self.table, buf, &ret)
                # print 'putting %s, %s' % (val, count)
                if not ret:
                    kh_del_str(self.table, k)
                count += 1
                uniques.append(val)

        # return None
        return uniques

    def factorize(self, ndarray[object] values):
        cdef:
            Py_ssize_t i, n = len(values)
            ndarray[int64_t] labels = np.empty(n, dtype=np.int64)
            ndarray[int64_t] counts = np.empty(n, dtype=np.int64)
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
                counts[idx] = counts[idx] + 1
            else:
                k = kh_put_str(self.table, buf, &ret)
                # print 'putting %s, %s' % (val, count)
                if not ret:
                    kh_del_str(self.table, k)

                self.table.vals[k] = count
                reverse[count] = val
                labels[i] = count
                counts[count] = 1
                count += 1

        # return None
        return reverse, labels, counts[:count].copy()

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
            ndarray[int64_t] counts = np.empty(n, dtype=np.int64)
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
                counts[idx] = counts[idx] + 1
            else:
                k = kh_put_int32(self.table, val, &ret)
                if not ret:
                    kh_del_int32(self.table, k)
                self.table.vals[k] = count
                reverse[count] = val
                labels[i] = count
                counts[count] = 1
                count += 1

        # return None
        return reverse, labels, counts[:count].copy()

cdef class Int64HashTable(HashTable):
    cdef kh_int64_t *table

    def __init__(self, size_hint=1):
        if size_hint is not None:
            kh_resize_int64(self.table, size_hint)

    def __cinit__(self):
        self.table = kh_init_int64()

    def __dealloc__(self):
        kh_destroy_int64(self.table)

    def __contains__(self, object key):
        cdef khiter_t k
        k = kh_get_int64(self.table, key)
        return k != self.table.n_buckets

    def __len__(self):
        return self.table.size

    cdef inline bint has_key(self, int64_t val):
        cdef khiter_t k
        k = kh_get_int64(self.table, val)
        return k != self.table.n_buckets

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

    def lookup_i4(self, ndarray[int64_t] values):
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
        labels, counts = self.get_labels(values, reverse, 0)
        return reverse, labels, counts

    def get_labels(self, ndarray[int64_t] values, list uniques,
                   Py_ssize_t count_prior, Py_ssize_t na_sentinel):
        cdef:
            Py_ssize_t i, n = len(values)
            ndarray[int64_t] labels
            ndarray[int64_t] counts
            Py_ssize_t idx, count = count_prior
            int ret = 0
            int64_t val
            khiter_t k

        labels = np.empty(n, dtype=np.int64)
        counts = np.empty(count_prior + n, dtype=np.int64)

        for i in range(n):
            val = values[i]
            k = kh_get_int64(self.table, val)
            if k != self.table.n_buckets:
                idx = self.table.vals[k]
                labels[i] = idx
                counts[idx] = counts[idx] + 1
            else:
                k = kh_put_int64(self.table, val, &ret)
                self.table.vals[k] = count
                uniques.append(val)
                labels[i] = count
                counts[count] = 1
                count += 1

        return labels, counts[:count].copy()

    def get_labels_groupby(self, ndarray[int64_t] values, list uniques):
        cdef:
            Py_ssize_t i, n = len(values)
            ndarray[int64_t] labels
            Py_ssize_t idx, count = 0
            int ret = 0
            int64_t val
            khiter_t k

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

        return labels

    def unique(self, ndarray[int64_t] values):
        cdef:
            Py_ssize_t i, n = len(values)
            Py_ssize_t idx, count = 0
            int ret = 0
            ndarray result
            int64_t val
            khiter_t k
            Int64Vector uniques = Int64Vector()

        # TODO: kvec

        for i in range(n):
            val = values[i]
            k = kh_get_int64(self.table, val)
            if k == self.table.n_buckets:
                k = kh_put_int64(self.table, val, &ret)
                uniques.append(val)
                count += 1

        result = uniques.to_array(xfer_data=True)

        # result = np.array(uniques, copy=False)
        # result.base = <PyObject*> uniques
        # Py_INCREF(uniques)

        return result

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

cdef class Float64HashTable(HashTable):
    cdef kh_float64_t *table

    def __init__(self, size_hint=1):
        if size_hint is not None:
            kh_resize_float64(self.table, size_hint)

    def __cinit__(self):
        self.table = kh_init_float64()

    def __len__(self):
        return self.table.size

    def __dealloc__(self):
        kh_destroy_float64(self.table)

    def factorize(self, ndarray[float64_t] values):
        uniques = []
        labels, counts = self.get_labels(values, uniques, 0, -1)
        return uniques, labels, counts

    cpdef get_labels(self, ndarray[float64_t] values, list uniques,
                     Py_ssize_t count_prior, int64_t na_sentinel):
        cdef:
            Py_ssize_t i, n = len(values)
            ndarray[int64_t] labels
            ndarray[int64_t] counts
            Py_ssize_t idx, count = count_prior
            int ret = 0
            float64_t val
            khiter_t k

        labels = np.empty(n, dtype=np.int64)
        counts = np.empty(count_prior + n, dtype=np.int64)

        for i in range(n):
            val = values[i]

            if val != val:
                labels[i] = na_sentinel
                continue

            k = kh_get_float64(self.table, val)
            if k != self.table.n_buckets:
                idx = self.table.vals[k]
                labels[i] = idx
                counts[idx] = counts[idx] + 1
            else:
                k = kh_put_float64(self.table, val, &ret)
                self.table.vals[k] = count
                uniques.append(val)
                labels[i] = count
                counts[count] = 1
                count += 1

        return labels, counts[:count].copy()

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
            list uniques = []
            bint seen_na = 0

        # TODO: kvec

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

        return uniques

cdef class PyObjectHashTable(HashTable):
    cdef kh_pymap_t *table

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

    cpdef destroy(self):
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

        result = uniques.to_array(xfer_data=True)

        # result = np.array(uniques, copy=False)
        # result.base = <PyObject*> uniques
        # Py_INCREF(uniques)

        return result

    cpdef get_labels(self, ndarray[object] values, list uniques,
                     Py_ssize_t count_prior, int64_t na_sentinel):
        cdef:
            Py_ssize_t i, n = len(values)
            ndarray[int64_t] labels
            ndarray[int64_t] counts
            Py_ssize_t idx, count = count_prior
            int ret = 0
            object val
            khiter_t k

        labels = np.empty(n, dtype=np.int64)
        counts = np.empty(count_prior + n, dtype=np.int64)

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
                counts[idx] = counts[idx] + 1
            else:
                k = kh_put_pymap(self.table, <PyObject*>val, &ret)
                self.table.vals[k] = count
                uniques.append(val)
                labels[i] = count
                counts[count] = 1
                count += 1

        return labels, counts[:count].copy()


cdef class Factorizer:
    cdef public PyObjectHashTable table
    cdef public uniques
    cdef public Py_ssize_t count

    def __init__(self, size_hint):
        self.table = PyObjectHashTable(size_hint)
        self.uniques = []
        self.count = 0

    def get_count(self):
        return self.count

    def factorize(self, ndarray[object] values, sort=False, na_sentinel=-1):
        labels, counts = self.table.get_labels(values, self.uniques,
                                               self.count, na_sentinel)

        # sort on
        if sort:
            if labels.dtype != np.int_:
                labels = labels.astype(np.int_)

            sorter = list_to_object_array(self.uniques).argsort()
            reverse_indexer = np.empty(len(sorter), dtype=np.int_)
            reverse_indexer.put(sorter, np.arange(len(sorter)))

            labels = reverse_indexer.take(labels)
            counts = counts.take(sorter)

        self.count = len(counts)
        return labels, counts

    def unique(self, ndarray[object] values):
        # just for fun
        return self.table.unique(values)


cdef class Int64Factorizer:
    cdef public Int64HashTable table
    cdef public list uniques
    cdef public Py_ssize_t count

    def __init__(self, size_hint):
        self.table = Int64HashTable(size_hint)
        self.uniques = []
        self.count = 0

    def get_count(self):
        return self.count

    def factorize(self, ndarray[int64_t] values, sort=False,
                  na_sentinel=-1):
        labels, counts = self.table.get_labels(values, self.uniques,
                                               self.count, na_sentinel)

        # sort on
        if sort:
            if labels.dtype != np.int_:
                labels = labels.astype(np.int_)

            sorter = list_to_object_array(self.uniques).argsort()
            reverse_indexer = np.empty(len(sorter), dtype=np.int_)
            reverse_indexer.put(sorter, np.arange(len(sorter)))

            labels = reverse_indexer.take(labels)
            counts = counts.take(sorter)

        self.count = len(counts)
        return labels, counts


cdef class DictFactorizer:

    cdef public:
        dict table
        list uniques
        Py_ssize_t count

    def __init__(self, table=None, uniques=None):
        if table is None:
            self.table = {}
        else:
            self.table = table

        if uniques is None:
            self.uniques = []
            self.count = 0
        else:
            self.uniques = uniques
            self.count = len(uniques)

    def get_count(self):
        return self.count

    def get_labels(self, ndarray[object] values):
        cdef:
            Py_ssize_t i, n = len(values)
            ndarray[int64_t] labels
            ndarray[int64_t] counts
            Py_ssize_t idx, count = self.count
            int ret = 0
            object val

        labels = np.empty(n, dtype=np.int64)
        counts = np.empty(count + n, dtype=np.int64)

        for i in range(n):
            val = values[i]

            if val in self.table:
                idx = self.table[val]
                labels[i] = idx
                counts[idx] = counts[idx] + 1
            else:
                self.table[val] = count
                self.uniques.append(val)
                labels[i] = count
                counts[count] = 1
                count += 1

        return labels, counts[:count].copy()

    def factorize(self, ndarray[object] values, sort=False):
        labels, counts = self.get_labels(values)

        # sort on
        if sort:
            if labels.dtype != np.int_:
                labels = labels.astype(np.int_)

            sorter = list_to_object_array(self.uniques).argsort()
            reverse_indexer = np.empty(len(sorter), dtype=np.int_)
            reverse_indexer.put(sorter, np.arange(len(sorter)))

            labels = reverse_indexer.take(labels)
            counts = counts.take(sorter)

        self.count = len(counts)
        return labels, counts

    def unique(self, ndarray[object] values):
        cdef:
            Py_ssize_t i, n = len(values)
            Py_ssize_t idx, count = self.count
            object val

        for i in range(n):
            val = values[i]
            if val not in self.table:
                self.table[val] = count
                self.uniques.append(val)
                count += 1
        return self.uniques


    def unique_int64(self, ndarray[int64_t] values):
        cdef:
            Py_ssize_t i, n = len(values)
            Py_ssize_t idx, count = self.count
            int64_t val

        for i in range(n):
            val = values[i]
            if val not in self.table:
                self.table[val] = count
                self.uniques.append(val)
                count += 1
        return self.uniques

def lookup2(ndarray[object] values):
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

