"""
Template for each `dtype` helper function for hashtable

WARNING: DO NOT edit .pxi FILE directly, .pxi is generated from .pxi.in
"""

#----------------------------------------------------------------------
# VectorData
#----------------------------------------------------------------------


ctypedef struct Float64VectorData:
    float64_t *data
    size_t n, m


@cython.wraparound(False)
@cython.boundscheck(False)
cdef void append_data_float64(Float64VectorData *data,
                                float64_t x) nogil:

    data.data[data.n] = x
    data.n += 1


ctypedef struct Int64VectorData:
    int64_t *data
    size_t n, m


@cython.wraparound(False)
@cython.boundscheck(False)
cdef void append_data_int64(Int64VectorData *data,
                                int64_t x) nogil:

    data.data[data.n] = x
    data.n += 1

ctypedef fused vector_data:
    Int64VectorData
    Float64VectorData

cdef bint needs_resize(vector_data *data) nogil:
    return data.n == data.m

#----------------------------------------------------------------------
# Vector
#----------------------------------------------------------------------

cdef class Float64Vector:

    cdef:
        Float64VectorData *data
        ndarray ao

    def __cinit__(self):
        self.data = <Float64VectorData *>PyMem_Malloc(
            sizeof(Float64VectorData))
        if not self.data:
            raise MemoryError()
        self.data.n = 0
        self.data.m = _INIT_VEC_CAP
        self.ao = np.empty(self.data.m, dtype=np.float64)
        self.data.data = <float64_t*> self.ao.data

    cdef resize(self):
        self.data.m = max(self.data.m * 4, _INIT_VEC_CAP)
        self.ao.resize(self.data.m)
        self.data.data = <float64_t*> self.ao.data

    def __dealloc__(self):
        PyMem_Free(self.data)

    def __len__(self):
        return self.data.n

    def to_array(self):
        self.ao.resize(self.data.n)
        self.data.m = self.data.n
        return self.ao

    cdef inline void append(self, float64_t x):

        if needs_resize(self.data):
            self.resize()

        append_data_float64(self.data, x)

cdef class Int64Vector:

    cdef:
        Int64VectorData *data
        ndarray ao

    def __cinit__(self):
        self.data = <Int64VectorData *>PyMem_Malloc(
            sizeof(Int64VectorData))
        if not self.data:
            raise MemoryError()
        self.data.n = 0
        self.data.m = _INIT_VEC_CAP
        self.ao = np.empty(self.data.m, dtype=np.int64)
        self.data.data = <int64_t*> self.ao.data

    cdef resize(self):
        self.data.m = max(self.data.m * 4, _INIT_VEC_CAP)
        self.ao.resize(self.data.m)
        self.data.data = <int64_t*> self.ao.data

    def __dealloc__(self):
        PyMem_Free(self.data)

    def __len__(self):
        return self.data.n

    def to_array(self):
        self.ao.resize(self.data.n)
        self.data.m = self.data.n
        return self.ao

    cdef inline void append(self, int64_t x):

        if needs_resize(self.data):
            self.resize()

        append_data_int64(self.data, x)


cdef class ObjectVector:

    cdef:
        PyObject **data
        size_t n, m
        ndarray ao

    def __cinit__(self):
        self.n = 0
        self.m = _INIT_VEC_CAP
        self.ao = np.empty(_INIT_VEC_CAP, dtype=object)
        self.data = <PyObject**> self.ao.data

    def __len__(self):
        return self.n

    cdef inline append(self, object o):
        if self.n == self.m:
            self.m = max(self.m * 2, _INIT_VEC_CAP)
            self.ao.resize(self.m)
            self.data = <PyObject**> self.ao.data

        Py_INCREF(o)
        self.data[self.n] = <PyObject*> o
        self.n += 1

    def to_array(self):
        self.ao.resize(self.n)
        self.m = self.n
        return self.ao


#----------------------------------------------------------------------
# HashTable
#----------------------------------------------------------------------


cdef class HashTable:
    pass

cdef class Float64HashTable(HashTable):

    def __cinit__(self, size_hint=1):
        self.table = kh_init_float64()
        if size_hint is not None:
            kh_resize_float64(self.table, size_hint)

    def __len__(self):
        return self.table.size

    def __dealloc__(self):
        kh_destroy_float64(self.table)

    def __contains__(self, object key):
        cdef khiter_t k
        k = kh_get_float64(self.table, key)
        return k != self.table.n_buckets

    cpdef get_item(self, float64_t val):
        cdef khiter_t k
        k = kh_get_float64(self.table, val)
        if k != self.table.n_buckets:
            return self.table.vals[k]
        else:
            raise KeyError(val)

    def get_iter_test(self, float64_t key, Py_ssize_t iterations):
        cdef Py_ssize_t i, val=0
        for i in range(iterations):
            k = kh_get_float64(self.table, val)
            if k != self.table.n_buckets:
                val = self.table.vals[k]

    cpdef set_item(self, float64_t key, Py_ssize_t val):
        cdef:
            khiter_t k
            int ret = 0

        k = kh_put_float64(self.table, key, &ret)
        self.table.keys[k] = key
        if kh_exist_float64(self.table, k):
            self.table.vals[k] = val
        else:
            raise KeyError(key)

    @cython.boundscheck(False)
    def map(self, float64_t[:] keys, int64_t[:] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            float64_t key
            khiter_t k

        with nogil:
            for i in range(n):
                key = keys[i]
                k = kh_put_float64(self.table, key, &ret)
                self.table.vals[k] = <Py_ssize_t> values[i]

    @cython.boundscheck(False)
    def map_locations(self, ndarray[float64_t, ndim=1] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            float64_t val
            khiter_t k

        with nogil:
            for i in range(n):
                val = values[i]
                k = kh_put_float64(self.table, val, &ret)
                self.table.vals[k] = i

    @cython.boundscheck(False)
    def lookup(self, float64_t[:] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            float64_t val
            khiter_t k
            int64_t[:] locs = np.empty(n, dtype=np.int64)

        with nogil:
            for i in range(n):
                val = values[i]
                k = kh_get_float64(self.table, val)
                if k != self.table.n_buckets:
                    locs[i] = self.table.vals[k]
                else:
                    locs[i] = -1

        return np.asarray(locs)

    def factorize(self, float64_t values):
        uniques = Float64Vector()
        labels = self.get_labels(values, uniques, 0, 0)
        return uniques.to_array(), labels

    @cython.boundscheck(False)
    def get_labels(self, float64_t[:] values, Float64Vector uniques,
                   Py_ssize_t count_prior, Py_ssize_t na_sentinel,
                   bint check_null=True):
        cdef:
            Py_ssize_t i, n = len(values)
            int64_t[:] labels
            Py_ssize_t idx, count = count_prior
            int ret = 0
            float64_t val
            khiter_t k
            Float64VectorData *ud

        labels = np.empty(n, dtype=np.int64)
        ud = uniques.data

        with nogil:
            for i in range(n):
                val = values[i]

                if check_null and val != val:
                    labels[i] = na_sentinel
                    continue

                k = kh_get_float64(self.table, val)

                if k != self.table.n_buckets:
                    idx = self.table.vals[k]
                    labels[i] = idx
                else:
                    k = kh_put_float64(self.table, val, &ret)
                    self.table.vals[k] = count

                    if needs_resize(ud):
                        with gil:
                            uniques.resize()
                    append_data_float64(ud, val)
                    labels[i] = count
                    count += 1

        return np.asarray(labels)

    @cython.boundscheck(False)
    def get_labels_groupby(self, float64_t[:] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int64_t[:] labels
            Py_ssize_t idx, count = 0
            int ret = 0
            float64_t val
            khiter_t k
            Float64Vector uniques = Float64Vector()
            Float64VectorData *ud

        labels = np.empty(n, dtype=np.int64)
        ud = uniques.data

        with nogil:
            for i in range(n):
                val = values[i]

                # specific for groupby
                if val < 0:
                    labels[i] = -1
                    continue

                k = kh_get_float64(self.table, val)
                if k != self.table.n_buckets:
                    idx = self.table.vals[k]
                    labels[i] = idx
                else:
                    k = kh_put_float64(self.table, val, &ret)
                    self.table.vals[k] = count

                    if needs_resize(ud):
                        with gil:
                            uniques.resize()
                    append_data_float64(ud, val)
                    labels[i] = count
                    count += 1

        arr_uniques = uniques.to_array()

        return np.asarray(labels), arr_uniques

    @cython.boundscheck(False)
    def unique(self, float64_t[:] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            float64_t val
            khiter_t k
            bint seen_na = 0
            Float64Vector uniques = Float64Vector()
            Float64VectorData *ud

        ud = uniques.data

        with nogil:
            for i in range(n):
                val = values[i]

                if val == val:
                    k = kh_get_float64(self.table, val)
                    if k == self.table.n_buckets:
                        kh_put_float64(self.table, val, &ret)
                        if needs_resize(ud):
                            with gil:
                                uniques.resize()
                        append_data_float64(ud, val)
                elif not seen_na:
                    seen_na = 1
                    if needs_resize(ud):
                        with gil:
                            uniques.resize()
                    append_data_float64(ud, NAN)

        return uniques.to_array()

cdef class Int64HashTable(HashTable):

    def __cinit__(self, size_hint=1):
        self.table = kh_init_int64()
        if size_hint is not None:
            kh_resize_int64(self.table, size_hint)

    def __len__(self):
        return self.table.size

    def __dealloc__(self):
        kh_destroy_int64(self.table)

    def __contains__(self, object key):
        cdef khiter_t k
        k = kh_get_int64(self.table, key)
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

    @cython.boundscheck(False)
    def map(self, int64_t[:] keys, int64_t[:] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            int64_t key
            khiter_t k

        with nogil:
            for i in range(n):
                key = keys[i]
                k = kh_put_int64(self.table, key, &ret)
                self.table.vals[k] = <Py_ssize_t> values[i]

    @cython.boundscheck(False)
    def map_locations(self, ndarray[int64_t, ndim=1] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            int64_t val
            khiter_t k

        with nogil:
            for i in range(n):
                val = values[i]
                k = kh_put_int64(self.table, val, &ret)
                self.table.vals[k] = i

    @cython.boundscheck(False)
    def lookup(self, int64_t[:] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            int64_t val
            khiter_t k
            int64_t[:] locs = np.empty(n, dtype=np.int64)

        with nogil:
            for i in range(n):
                val = values[i]
                k = kh_get_int64(self.table, val)
                if k != self.table.n_buckets:
                    locs[i] = self.table.vals[k]
                else:
                    locs[i] = -1

        return np.asarray(locs)

    def factorize(self, int64_t values):
        uniques = Int64Vector()
        labels = self.get_labels(values, uniques, 0, 0)
        return uniques.to_array(), labels

    @cython.boundscheck(False)
    def get_labels(self, int64_t[:] values, Int64Vector uniques,
                   Py_ssize_t count_prior, Py_ssize_t na_sentinel,
                   bint check_null=True):
        cdef:
            Py_ssize_t i, n = len(values)
            int64_t[:] labels
            Py_ssize_t idx, count = count_prior
            int ret = 0
            int64_t val
            khiter_t k
            Int64VectorData *ud

        labels = np.empty(n, dtype=np.int64)
        ud = uniques.data

        with nogil:
            for i in range(n):
                val = values[i]

                if check_null and val == iNaT:
                    labels[i] = na_sentinel
                    continue

                k = kh_get_int64(self.table, val)

                if k != self.table.n_buckets:
                    idx = self.table.vals[k]
                    labels[i] = idx
                else:
                    k = kh_put_int64(self.table, val, &ret)
                    self.table.vals[k] = count

                    if needs_resize(ud):
                        with gil:
                            uniques.resize()
                    append_data_int64(ud, val)
                    labels[i] = count
                    count += 1

        return np.asarray(labels)

    @cython.boundscheck(False)
    def get_labels_groupby(self, int64_t[:] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int64_t[:] labels
            Py_ssize_t idx, count = 0
            int ret = 0
            int64_t val
            khiter_t k
            Int64Vector uniques = Int64Vector()
            Int64VectorData *ud

        labels = np.empty(n, dtype=np.int64)
        ud = uniques.data

        with nogil:
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

                    if needs_resize(ud):
                        with gil:
                            uniques.resize()
                    append_data_int64(ud, val)
                    labels[i] = count
                    count += 1

        arr_uniques = uniques.to_array()

        return np.asarray(labels), arr_uniques

    @cython.boundscheck(False)
    def unique(self, int64_t[:] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            int64_t val
            khiter_t k
            bint seen_na = 0
            Int64Vector uniques = Int64Vector()
            Int64VectorData *ud

        ud = uniques.data

        with nogil:
            for i in range(n):
                val = values[i]

                k = kh_get_int64(self.table, val)
                if k == self.table.n_buckets:
                    kh_put_int64(self.table, val, &ret)
                    if needs_resize(ud):
                        with gil:
                            uniques.resize()
                    append_data_int64(ud, val)

        return uniques.to_array()


cdef class StringHashTable(HashTable):
    cdef kh_str_t *table

    def __cinit__(self, int size_hint=1):
        self.table = kh_init_str()
        if size_hint is not None:
            kh_resize_str(self.table, size_hint)

    def __dealloc__(self):
        kh_destroy_str(self.table)

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
                kh_put_str(self.table, buf, &ret)
                uniques.append(val)

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

        return reverse, labels


na_sentinel = object

cdef class PyObjectHashTable(HashTable):

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
        if key != key or key is None:
            key = na_sentinel
        k = kh_get_pymap(self.table, <PyObject*>key)
        return k != self.table.n_buckets

    def destroy(self):
        kh_destroy_pymap(self.table)
        self.table = NULL

    cpdef get_item(self, object val):
        cdef khiter_t k
        if val != val or val is None:
            val = na_sentinel
        k = kh_get_pymap(self.table, <PyObject*>val)
        if k != self.table.n_buckets:
            return self.table.vals[k]
        else:
            raise KeyError(val)

    def get_iter_test(self, object key, Py_ssize_t iterations):
        cdef Py_ssize_t i, val
        if key != key or key is None:
            key = na_sentinel
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
        if key != key or key is None:
            key = na_sentinel
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
            if val != val or val is None:
                val = na_sentinel

            k = kh_put_pymap(self.table, <PyObject*>val, &ret)
            self.table.vals[k] = i

    def lookup(self, ndarray[object] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            object val
            khiter_t k
            int64_t[:] locs = np.empty(n, dtype=np.int64)

        for i in range(n):
            val = values[i]
            hash(val)
            if val != val or val is None:
                val = na_sentinel

            k = kh_get_pymap(self.table, <PyObject*>val)
            if k != self.table.n_buckets:
                locs[i] = self.table.vals[k]
            else:
                locs[i] = -1

        return np.asarray(locs)

    def unique(self, ndarray[object] values):
        cdef:
            Py_ssize_t i, n = len(values)
            int ret = 0
            object val
            khiter_t k
            ObjectVector uniques = ObjectVector()
            bint seen_na = 0

        for i in range(n):
            val = values[i]
            hash(val)
            if not _checknan(val):
                k = kh_get_pymap(self.table, <PyObject*>val)
                if k == self.table.n_buckets:
                    kh_put_pymap(self.table, <PyObject*>val, &ret)
                    uniques.append(val)
            elif not seen_na:
                seen_na = 1
                uniques.append(nan)

        return uniques.to_array()

    def get_labels(self, ndarray[object] values, ObjectVector uniques,
                   Py_ssize_t count_prior, int64_t na_sentinel,
                   bint check_null=True):
        cdef:
            Py_ssize_t i, n = len(values)
            int64_t[:] labels
            Py_ssize_t idx, count = count_prior
            int ret = 0
            object val
            khiter_t k

        labels = np.empty(n, dtype=np.int64)

        for i in range(n):
            val = values[i]
            hash(val)

            if check_null and val != val or val is None:
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

        return np.asarray(labels)