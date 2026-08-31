from numpy cimport int64_t


cdef extern from "pandas/portable.h":
    int checked_add(int64_t a, int64_t b, int64_t *res) noexcept nogil
    int checked_sub(int64_t a, int64_t b, int64_t *res) noexcept nogil
    int checked_mul(int64_t a, int64_t b, int64_t *res) noexcept nogil
