# -*- coding: utf-8 -*-
# cython: profile=False

from cpython.datetime cimport (datetime, date,
                               PyDateTime_IMPORT,
                               PyDateTime_GET_YEAR, PyDateTime_GET_MONTH,
                               PyDateTime_GET_DAY, PyDateTime_DATE_GET_HOUR,
                               PyDateTime_DATE_GET_MINUTE,
                               PyDateTime_DATE_GET_SECOND,
                               PyDateTime_DATE_GET_MICROSECOND)
PyDateTime_IMPORT

from numpy cimport int64_t

cdef extern from "numpy/ndarrayobject.h":
    ctypedef int64_t npy_timedelta
    ctypedef int64_t npy_datetime

cdef extern from "../src/datetime/np_datetime.h":
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

    int cmp_pandas_datetimestruct(pandas_datetimestruct *a,
                                  pandas_datetimestruct *b)

    npy_datetime pandas_datetimestruct_to_datetime(PANDAS_DATETIMEUNIT fr,
                                                   pandas_datetimestruct *d
                                                   ) nogil

    void pandas_datetime_to_datetimestruct(npy_datetime val,
                                           PANDAS_DATETIMEUNIT fr,
                                           pandas_datetimestruct *result) nogil

    pandas_datetimestruct _NS_MIN_DTS, _NS_MAX_DTS

# ----------------------------------------------------------------------


class OutOfBoundsDatetime(ValueError):
    pass


cdef inline check_dts_bounds(pandas_datetimestruct *dts):
    """Raises OutOfBoundsDatetime if the given date is outside the range that
    can be represented by nanosecond-resolution 64-bit integers."""
    cdef:
        bint error = False

    if (dts.year <= 1677 and
            cmp_pandas_datetimestruct(dts, &_NS_MIN_DTS) == -1):
        error = True
    elif (dts.year >= 2262 and
          cmp_pandas_datetimestruct(dts, &_NS_MAX_DTS) == 1):
        error = True

    if error:
        fmt = '%d-%.2d-%.2d %.2d:%.2d:%.2d' % (dts.year, dts.month,
                                               dts.day, dts.hour,
                                               dts.min, dts.sec)
        raise OutOfBoundsDatetime(
            'Out of bounds nanosecond timestamp: {fmt}'.format(fmt=fmt))


# ----------------------------------------------------------------------
# Conversion

cdef inline int64_t dtstruct_to_dt64(pandas_datetimestruct* dts) nogil:
    """Convenience function to call pandas_datetimestruct_to_datetime
    with the by-far-most-common frequency PANDAS_FR_ns"""
    return pandas_datetimestruct_to_datetime(PANDAS_FR_ns, dts)


cdef inline void dt64_to_dtstruct(int64_t dt64,
                                  pandas_datetimestruct* out) nogil:
    """Convenience function to call pandas_datetime_to_datetimestruct
    with the by-far-most-common frequency PANDAS_FR_ns"""
    pandas_datetime_to_datetimestruct(dt64, PANDAS_FR_ns, out)
    return


cdef inline int64_t pydatetime_to_dt64(datetime val,
                                       pandas_datetimestruct *dts):
    dts.year = PyDateTime_GET_YEAR(val)
    dts.month = PyDateTime_GET_MONTH(val)
    dts.day = PyDateTime_GET_DAY(val)
    dts.hour = PyDateTime_DATE_GET_HOUR(val)
    dts.min = PyDateTime_DATE_GET_MINUTE(val)
    dts.sec = PyDateTime_DATE_GET_SECOND(val)
    dts.us = PyDateTime_DATE_GET_MICROSECOND(val)
    dts.ps = dts.as = 0
    return dtstruct_to_dt64(dts)


cdef inline int64_t pydate_to_dt64(date val,
                                   pandas_datetimestruct *dts):
    dts.year = PyDateTime_GET_YEAR(val)
    dts.month = PyDateTime_GET_MONTH(val)
    dts.day = PyDateTime_GET_DAY(val)
    dts.hour = dts.min = dts.sec = dts.us = 0
    dts.ps = dts.as = 0
    return dtstruct_to_dt64(dts)
