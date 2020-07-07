from cpython.datetime cimport tzinfo
from numpy cimport int64_t, intp_t, ndarray


cdef int64_t tz_convert_utc_to_tzlocal(int64_t utc_val, tzinfo tz, bint* fold=*)
cpdef int64_t tz_convert_single(int64_t val, tzinfo tz1, tzinfo tz2)
cdef int64_t tz_localize_to_utc_single(
    int64_t val, tzinfo tz, object ambiguous=*, object nonexistent=*
) except? -1


cdef class Localizer:
    cdef:
        bint use_utc, use_tzlocal, use_fixed, use_pytz
        int noffsets
        int64_t* utcoffsets
        intp_t* positions
        ndarray positions_arr  # needed to avoid segfault
        int64_t delta
        tzinfo tz

    cdef inline int64_t get_local_timestamp(self, int64_t utc_value, Py_ssize_t i)
