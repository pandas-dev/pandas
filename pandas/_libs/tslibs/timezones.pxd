# -*- coding: utf-8 -*-
# cython: profile=False

from numpy cimport ndarray, int64_t

cdef bint is_utc(object tz)
cdef bint is_tzlocal(object tz)

cdef bint treat_tz_as_pytz(object tz)
cdef bint treat_tz_as_dateutil(object tz)

cpdef object get_timezone(object tz)
cpdef object maybe_get_tz(object tz)

cpdef get_utcoffset(tzinfo, obj)
cdef bint is_fixed_offset(object tz)

cdef object get_dst_info(object tz)

cdef ndarray[int64_t] _infer_dst(ndarray[int64_t] vals,
                                 ndarray[int64_t] result_a,
                                 ndarray[int64_t] result_b)
