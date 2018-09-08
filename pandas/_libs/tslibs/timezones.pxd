# -*- coding: utf-8 -*-

cdef bint is_utc(object tz)
cdef bint is_tzlocal(object tz)

cdef bint treat_tz_as_pytz(object tz)
cdef bint treat_tz_as_dateutil(object tz)

cpdef bint tz_compare(object start, object end)
cpdef object get_timezone(object tz)
cpdef object maybe_get_tz(object tz)

cdef get_utcoffset(tzinfo, obj)
cdef bint is_fixed_offset(object tz)

cdef object get_dst_info(object tz, bint dst)
