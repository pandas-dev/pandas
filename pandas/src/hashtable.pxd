from khash cimport *

# prototypes for sharing

# cdef class StringHashTable:
#     cdef kh_str_t *table

#     cdef inline int check_type(self, object)
#     cpdef get_item(self, object)
#     cpdef set_item(self, object, Py_ssize_t)

# cdef class Int32HashTable:
#     cdef kh_int32_t *table

#     cdef inline int check_type(self, object)
#     cpdef get_item(self, int32_t)
#     cpdef set_item(self, int32_t, Py_ssize_t)

# cdef class Int64HashTable:
#     cdef kh_int64_t *table

#     cdef inline bint has_key(self, int64_t)
#     cpdef get_item(self, int64_t)
#     cpdef set_item(self, int64_t, Py_ssize_t)


# cdef class Float64HashTable:
#     cdef kh_float64_t *table

#     cpdef get_labels(self, ndarray, list, Py_ssize_t, int32_t)


# cdef class PyObjectHashTable:
#     cdef kh_pymap_t *table

#     cdef destroy(self)
#     cpdef get_item(self, object)
#     cpdef set_item(self, object, Py_ssize_t)
#     cpdef get_labels(self, ndarray, list, Py_ssize_t, int32_t)


# cdef class Factorizer:
#     cdef public PyObjectHashTable table
#     cdef public uniques
#     cdef public Py_ssize_t count


# cdef class Int64Factorizer:
#     cdef public Int64HashTable table
#     cdef public list uniques
#     cdef public Py_ssize_t count
