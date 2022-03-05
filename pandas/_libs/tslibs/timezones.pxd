from cpython.datetime cimport (
    datetime,
    timedelta,
    tzinfo,
)
from numpy cimport (
    int64_t,
    intp_t,
    ndarray,
)


cdef tzinfo utc_pytz

cpdef bint is_utc(tzinfo tz)
cdef bint is_tzlocal(tzinfo tz)

cdef bint treat_tz_as_pytz(tzinfo tz)

cpdef bint tz_compare(tzinfo start, tzinfo end)
cpdef object get_timezone(tzinfo tz)
cpdef tzinfo maybe_get_tz(object tz)

cdef timedelta get_utcoffset(tzinfo tz, datetime obj)
cdef bint is_fixed_offset(tzinfo tz)

cdef object get_dst_info(tzinfo tz)


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

    cdef int64_t prepare1(self, int64_t utc_val)
    cdef intp_t* prepare(self, const int64_t[:] stamps)
