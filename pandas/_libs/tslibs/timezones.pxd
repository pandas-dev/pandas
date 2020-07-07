from cpython.datetime cimport datetime, timedelta, tzinfo

from numpy cimport int64_t, intp_t, ndarray

cdef tzinfo utc_pytz

cpdef bint is_utc(tzinfo tz)
cdef bint is_tzlocal(tzinfo tz)

cdef bint treat_tz_as_pytz(tzinfo tz)

cpdef bint tz_compare(tzinfo start, tzinfo end)
cpdef object get_timezone(tzinfo tz)
cpdef object maybe_get_tz(object tz)

cdef timedelta get_utcoffset(tzinfo tz, datetime obj)
cdef bint is_fixed_offset(tzinfo tz)

cdef object get_dst_info(tzinfo tz)


cdef class TZ:
    cdef:
        bint use_utc, use_tzlocal, use_fixed, use_pytz
        int noffsets
        int64_t* utcoffsets
        intp_t* positions
        ndarray positions_arr  # needed to avoid segfault
        int64_t delta
        tzinfo tz

    cdef inline int64_t get_local_timestamp(self, int64_t utc_value, Py_ssize_t i)
