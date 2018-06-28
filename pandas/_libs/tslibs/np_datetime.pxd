# -*- coding: utf-8 -*-
# cython: profile=False

from cpython.datetime cimport date, datetime

from numpy cimport int64_t, int32_t

cdef extern from "numpy/ndarrayobject.h":
    ctypedef int64_t npy_timedelta
    ctypedef int64_t npy_datetime

cdef extern from "numpy/ndarraytypes.h":
    ctypedef struct PyArray_DatetimeMetaData:
        PANDAS_DATETIMEUNIT base
        int64_t num

cdef extern from "numpy/arrayscalars.h":
    ctypedef struct PyDatetimeScalarObject:
        # PyObject_HEAD
        npy_datetime obval
        PyArray_DatetimeMetaData obmeta

    ctypedef struct PyTimedeltaScalarObject:
        # PyObject_HEAD
        npy_timedelta obval
        PyArray_DatetimeMetaData obmeta

cdef extern from "../src/datetime/np_datetime.h":
    ctypedef struct pandas_datetimestruct:
        int64_t year
        int32_t month, day, hour, min, sec, us, ps, as

    ctypedef struct pandas_timedeltastruct:
        int64_t days
        int32_t hrs, min, sec, ms, us, ns, seconds, microseconds, nanoseconds

    ctypedef enum PANDAS_DATETIMEUNIT:
        PANDAS_FR_Y
        PANDAS_FR_M
        PANDAS_FR_W
        PANDAS_FR_D
        PANDAS_FR_B
        PANDAS_FR_h
        PANDAS_FR_m
        PANDAS_FR_s
        PANDAS_FR_ms
        PANDAS_FR_us
        PANDAS_FR_ns
        PANDAS_FR_ps
        PANDAS_FR_fs
        PANDAS_FR_as

    void pandas_datetime_to_datetimestruct(npy_datetime val,
                                           PANDAS_DATETIMEUNIT fr,
                                           pandas_datetimestruct *result) nogil


cdef int reverse_ops[6]

cdef bint cmp_scalar(int64_t lhs, int64_t rhs, int op) except -1

cdef check_dts_bounds(pandas_datetimestruct *dts)

cdef int64_t dtstruct_to_dt64(pandas_datetimestruct* dts) nogil
cdef void dt64_to_dtstruct(int64_t dt64, pandas_datetimestruct* out) nogil
cdef void td64_to_tdstruct(int64_t td64, pandas_timedeltastruct* out) nogil

cdef int64_t pydatetime_to_dt64(datetime val, pandas_datetimestruct *dts)
cdef int64_t pydate_to_dt64(date val, pandas_datetimestruct *dts)

cdef npy_datetime get_datetime64_value(object obj) nogil
cdef npy_timedelta get_timedelta64_value(object obj) nogil
cdef PANDAS_DATETIMEUNIT get_datetime64_unit(object obj) nogil

cdef int _string_to_dts(object val, pandas_datetimestruct* dts,
                        int* out_local, int* out_tzoffset) except? -1
