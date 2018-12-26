# -*- coding: utf-8 -*-

from cpython.datetime cimport datetime, tzinfo

from numpy cimport int64_t, int32_t

from pandas._libs.tslibs.np_datetime cimport npy_datetimestruct


cdef class _TSObject:
    cdef:
        npy_datetimestruct dts      # npy_datetimestruct
        int64_t value               # numpy dt64
        object tzinfo


cdef convert_to_tsobject(object ts, object tz, object unit,
                         bint dayfirst, bint yearfirst,
                         int32_t nanos=*)

cdef _TSObject convert_datetime_to_tsobject(datetime ts, object tz,
                                            int32_t nanos=*)

cpdef int64_t tz_convert_single(int64_t val, object tz1, object tz2)

cdef int64_t get_datetime64_nanos(object val) except? -1

cpdef int64_t pydt_to_i8(object pydt) except? -1

cdef maybe_datetimelike_to_i8(object val)

cdef int64_t tz_convert_utc_to_tzlocal(int64_t utc_val, tzinfo tz)

cpdef datetime localize_pydatetime(datetime dt, object tz)
