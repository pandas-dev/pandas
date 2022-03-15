from cpython.datetime cimport tzinfo
from numpy cimport int64_t


cdef int64_t tz_convert_utc_to_tzlocal(
    int64_t utc_val, tzinfo tz, bint* fold=*
) except? -1
cpdef int64_t tz_convert_from_utc_single(int64_t val, tzinfo tz)
cdef int64_t tz_localize_to_utc_single(
    int64_t val, tzinfo tz, object ambiguous=*, object nonexistent=*
) except? -1

cdef Py_ssize_t bisect_right_i8(int64_t *data, int64_t val, Py_ssize_t n)
