from cpython.datetime cimport tzinfo
from numpy cimport (
    int64_t,
    intp_t,
    ndarray,
)


cdef int64_t tz_convert_utc_to_tzlocal(
    int64_t utc_val, tzinfo tz, bint* fold=*
) except? -1
cpdef int64_t tz_convert_from_utc_single(int64_t val, tzinfo tz)
cdef int64_t tz_localize_to_utc_single(
    int64_t val, tzinfo tz, object ambiguous=*, object nonexistent=*
) except? -1


cdef class Localizer:
    cdef:
        tzinfo tz
        bint use_utc
        bint use_fixed
        bint use_tzlocal
        bint use_dst
        bint use_pytz
        ndarray trans
        int64_t[:] deltas
        int64_t delta
        str typ

    cdef intp_t* prepare1(self, int64_t utc_val)
    cdef intp_t* prepare(self, const int64_t[:] stamps)
    cdef int64_t utc_val_to_local_val(self, int64_t utc_val, intp_t* pos, Py_ssize_t i)
