# -*- coding: utf-8 -*-
# cython: profile=False

from cpython.datetime cimport date, datetime

from numpy cimport int64_t, int32_t

cdef extern from "numpy/ndarrayobject.h":
    ctypedef int64_t npy_timedelta
    ctypedef int64_t npy_datetime

cdef extern from "numpy/ndarraytypes.h":
    # ctypedef struct npy_datetimestruct:
    #    int64_t year
    #    int32_t month, day, hour, min, sec, us, ps, as

    ctypedef enum NPY_DATETIMEUNIT:
        NPY_FR_Y = 0    # Years
        NPY_FR_M = 1    # Months
        NPY_FR_W = 2    # Weeks
        # Gap where 1.6 NPY_FR_B (value 3) was
        NPY_FR_D = 4    # Days
        NPY_FR_h = 5    # hours
        NPY_FR_m = 6    # minutes
        NPY_FR_s = 7    # seconds
        NPY_FR_ms = 8   # milliseconds
        NPY_FR_us = 9   # microseconds
        NPY_FR_ns = 10  # nanoseconds
        NPY_FR_ps = 11  # picoseconds
        NPY_FR_fs = 12  # femtoseconds
        NPY_FR_as = 13  # attoseconds
        NPY_FR_GENERIC = 14  # Generic, unbound units, can convert to anything

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

cdef check_dts_bounds(pandas_datetimestruct *dts)

cdef int64_t dtstruct_to_dt64(pandas_datetimestruct* dts) nogil
cdef void dt64_to_dtstruct(int64_t dt64, pandas_datetimestruct* out) nogil

cdef int64_t pydatetime_to_dt64(datetime val, pandas_datetimestruct *dts)
cdef int64_t pydate_to_dt64(date val, pandas_datetimestruct *dts)

cdef npy_datetime get_datetime64_value(object obj) nogil
cdef npy_timedelta get_timedelta64_value(object obj) nogil
cdef PANDAS_DATETIMEUNIT get_datetime64_unit(object obj) nogil
