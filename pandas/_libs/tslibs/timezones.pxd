#!/usr/bin/env python
# -*- coding: utf-8 -*-
# cython: profile=False

from numpy cimport ndarray, float64_t

cdef object _get_zone(object tz)
cdef object _tz_cache_key(object tz)
cdef bint _is_utc(object tz)
cdef bint _is_tzlocal(object tz)
cdef bint _treat_tz_as_pytz(object tz)
cdef bint _treat_tz_as_dateutil(object tz)
cpdef object maybe_get_tz(object tz)

cpdef _get_utcoffset(tzinfo, obj)
cpdef ndarray _unbox_utcoffsets(object transinfo)
cdef bint _is_fixed_offset(object tz)
cdef object _get_utc_trans_times_from_dateutil_tz(object tz)

cpdef object _get_dst_info(object tz)

cdef float64_t total_seconds(object td)
