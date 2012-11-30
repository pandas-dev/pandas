from khash cimport *

# prototypes for sharing

cdef class HashTable:
    pass

cdef class Int64HashTable(HashTable):
    cdef kh_int64_t *table

    cpdef get_item(self, int64_t val)
    cpdef set_item(self, int64_t key, Py_ssize_t val)


cdef class Float64HashTable(HashTable):
    cdef kh_float64_t *table

    # cpdef get_item(self, float64_t val)
    # cpdef set_item(self, float64_t key, Py_ssize_t val)

cdef class PyObjectHashTable(HashTable):
    cdef kh_pymap_t *table

    cpdef get_item(self, object val)
    cpdef set_item(self, object key, Py_ssize_t val)
