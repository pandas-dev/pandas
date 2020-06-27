from cpython.datetime cimport tzinfo

cdef tzinfo utc_pytz

cpdef bint is_utc(object tz)
cdef bint is_tzlocal(object tz)

cdef bint treat_tz_as_pytz(object tz)

cpdef bint tz_compare(object start, object end)
cpdef object get_timezone(object tz)
cpdef object maybe_get_tz(object tz)

cdef get_utcoffset(tzinfo tz, obj)
cdef bint is_fixed_offset(tzinfo tz)

cdef object get_dst_info(object tz)
