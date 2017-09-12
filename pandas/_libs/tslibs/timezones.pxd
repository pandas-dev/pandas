# -*- coding: utf-8 -*-
# cython: profile=False

cdef bint _is_utc(object tz)
cdef bint _is_tzlocal(object tz)

cdef bint _treat_tz_as_pytz(object tz)
cdef bint _treat_tz_as_dateutil(object tz)

cdef object _get_zone(object tz)

cpdef _get_utcoffset(tzinfo, obj)
