from pandas._libs.khash cimport (
    kh_int64_t, kh_uint64_t, kh_float64_t, kh_pymap_t, kh_str_t, uint64_t,
    int64_t, float64_t)
from numpy cimport ndarray

# prototypes for sharing

cdef class HashTable:
    pass

cdef class UInt64HashTable(HashTable):
    cdef kh_uint64_t *table

    cpdef get_item(self, uint64_t val)
    cpdef set_item(self, uint64_t key, Py_ssize_t val)

cdef class Int64HashTable(HashTable):
    cdef kh_int64_t *table

    cpdef get_item(self, int64_t val)
    cpdef set_item(self, int64_t key, Py_ssize_t val)

cdef class Float64HashTable(HashTable):
    cdef kh_float64_t *table

    cpdef get_item(self, float64_t val)
    cpdef set_item(self, float64_t key, Py_ssize_t val)

cdef class PyObjectHashTable(HashTable):
    cdef kh_pymap_t *table

    cpdef get_item(self, object val)
    cpdef set_item(self, object key, Py_ssize_t val)


cdef class StringHashTable(HashTable):
    cdef kh_str_t *table

    cpdef get_item(self, object val)
    cpdef set_item(self, object key, Py_ssize_t val)

cdef struct Int64VectorData:
    int64_t *data
    size_t n, m

cdef class Int64Vector:
    cdef Int64VectorData *data
    cdef ndarray ao
    cdef bint external_view_exists

    cdef resize(self)
    cpdef to_array(self)
    cdef inline void append(self, int64_t x)
    cdef extend(self, int64_t[:] x)
