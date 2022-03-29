from cpython.datetime cimport tzinfo
from numpy cimport (
    int64_t,
    intp_t,
)


cdef int64_t localize_tzinfo_api(
    int64_t utc_val, tzinfo tz, bint* fold=*
) except? -1
cdef int64_t tz_convert_from_utc_single(
    int64_t utc_val, tzinfo tz, bint* fold=?, Py_ssize_t* outpos=?
) except? -1
cdef int64_t tz_localize_to_utc_single(
    int64_t val, tzinfo tz, object ambiguous=*, object nonexistent=*
) except? -1

cdef Py_ssize_t bisect_right_i8(int64_t *data, int64_t val, Py_ssize_t n)

cdef bint infer_dateutil_fold(
    int64_t value,
    const int64_t[::1] trans,
    const int64_t[::1] deltas,
    intp_t pos,
)
