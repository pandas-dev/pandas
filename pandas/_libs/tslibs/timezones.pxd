# -*- coding: utf-8 -*-
# cython: profile=False

cdef bint is_utc(object tz)
cdef bint is_tzlocal(object tz)

cdef bint treat_tz_as_pytz(object tz)
cdef bint treat_tz_as_dateutil(object tz)

cpdef object get_timezone(object tz)

cpdef get_utcoffset(tzinfo, obj)
