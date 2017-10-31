# -*- coding: utf-8 -*-
# cython: profile=False

from numpy cimport int64_t

from datetime cimport pandas_datetimestruct


cdef class _TSObject:
    cdef:
        pandas_datetimestruct dts      # pandas_datetimestruct
        int64_t value               # numpy dt64
        object tzinfo

cdef void _localize_tso(_TSObject obj, object tz)

cpdef int64_t tz_convert_single(int64_t val, object tz1, object tz2)
