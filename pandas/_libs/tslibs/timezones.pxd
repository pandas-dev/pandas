from cpython.datetime cimport tzinfo

from numpy cimport int64_t, intp_t, ndarray

cdef tzinfo utc_pytz

cpdef bint is_utc(tzinfo tz)
cdef bint is_tzlocal(tzinfo tz)

cdef bint treat_tz_as_pytz(tzinfo tz)

cpdef bint tz_compare(object start, object end)
cpdef object get_timezone(object tz)
cpdef object maybe_get_tz(object tz)

cdef get_utcoffset(tzinfo tz, obj)
cdef bint is_fixed_offset(tzinfo tz)

cdef object get_dst_info(tzinfo tz)


ctypedef struct TZConvertInfo:
    bint use_utc
    bint use_tzlocal
    bint use_fixed
    ndarray[int64_t, ndim=1]* utcoffsets
    ndarray[intp_t, ndim=1]* positions
    int64_t delta
    int noffsets

cdef TZConvertInfo get_tzconverter(tzinfo tz, const int64_t[:] values)
