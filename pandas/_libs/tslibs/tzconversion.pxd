from cpython.datetime cimport tzinfo
from numpy cimport (
    int64_t,
    ndarray,
)


cdef int64_t localize_tzinfo_api(
    int64_t utc_val, tzinfo tz, bint* fold=*
) except? -1
cpdef int64_t tz_convert_from_utc_single(int64_t val, tzinfo tz)
cdef int64_t tz_localize_to_utc_single(
    int64_t val, tzinfo tz, object ambiguous=*, object nonexistent=*
) except? -1

cdef Py_ssize_t bisect_right_i8(int64_t *data, int64_t val, Py_ssize_t n)


cdef class Localizer:
    cdef readonly:
        tzinfo tz
        bint use_utc, use_fixed, use_tzlocal, use_dst, use_pytz, use_dateutil
        ndarray trans
        Py_ssize_t ntrans
        const int64_t[::1] deltas
        int64_t delta

    cdef inline int64_t utc_val_to_local_val(
        self, int64_t utc_val, bint* fold=*
    )
    cdef tzinfo adjust_pytz_tzinfo(self, int64_t utc_val)
