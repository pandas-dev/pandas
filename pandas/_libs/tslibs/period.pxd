from numpy cimport int64_t

cdef class _Period:
    cdef readonly:
        int64_t ordinal
        object freq
