# -*- coding: utf-8 -*-
# cython: profile=False

from cpython.datetime cimport date, datetime

from numpy cimport int64_t, int32_t


cdef extern from "../src/datetime/np_datetime.h":
    ctypedef struct pandas_datetimestruct:
        int64_t year
        int32_t month, day, hour, min, sec, us, ps, as


cdef int reverse_ops[6]

cdef bint cmp_scalar(int64_t lhs, int64_t rhs, int op) except -1

cdef check_dts_bounds(pandas_datetimestruct *dts)

cdef int64_t dtstruct_to_dt64(pandas_datetimestruct* dts) nogil
cdef void dt64_to_dtstruct(int64_t dt64, pandas_datetimestruct* out) nogil

cdef int64_t pydatetime_to_dt64(datetime val, pandas_datetimestruct *dts)
cdef int64_t pydate_to_dt64(date val, pandas_datetimestruct *dts)
