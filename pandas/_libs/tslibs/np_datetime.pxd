# -*- coding: utf-8 -*-
# cython: profile=False

from numpy cimport int64_t, int32_t


cdef extern from "../src/datetime/np_datetime.h":
    ctypedef struct pandas_datetimestruct:
        int64_t year
        int32_t month, day, hour, min, sec, us, ps, as


cdef check_dts_bounds(pandas_datetimestruct *dts)

cdef int64_t dtstruct_to_dt64(pandas_datetimestruct* dts) nogil
cdef void dt64_to_dtstruct(int64_t dt64, pandas_datetimestruct* out) nogil
