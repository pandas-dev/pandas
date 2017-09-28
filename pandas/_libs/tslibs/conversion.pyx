# -*- coding: utf-8 -*-
# cython: profile=False


from datetime import (
	date as pydate,
	datetime as pydatetime)

from cpython.datetime cimport (
	PyDateTime_Check,
    PyDateTime_GET_YEAR,
    PyDateTime_GET_MONTH,
    PyDateTime_GET_DAY,
    PyDateTime_DATE_GET_HOUR,
    PyDateTime_DATE_GET_MINUTE,
    PyDateTime_DATE_GET_SECOND,
    PyDateTime_IMPORT)
PyDateTime_IMPORT

from datetime cimport (get_datetime64_value, _pydatetime_to_dts,
                       pandas_datetimestruct)

import numpy as np
cimport numpy as cnp
from numpy cimport int64_t, ndarray
cnp.import_array()

cimport util

from timezones cimport get_utcoffset, is_utc

# ----------------------------------------------------------------------
# Constants
cdef int _EPOCH_ORD = 719163

# ----------------------------------------------------------------------
# Non-pandas-specific

cpdef object to_datetime(int64_t timestamp):
    return pydatetime.utcfromtimestamp(timestamp / 1000.0)


cdef inline int64_t gmtime(object date):
    cdef int y, m, d, h, mn, s, days

    y = PyDateTime_GET_YEAR(date)
    m = PyDateTime_GET_MONTH(date)
    d = PyDateTime_GET_DAY(date)
    h = PyDateTime_DATE_GET_HOUR(date)
    mn = PyDateTime_DATE_GET_MINUTE(date)
    s = PyDateTime_DATE_GET_SECOND(date)

    days = pydate(y, m, 1).toordinal() - _EPOCH_ORD + d - 1
    return ((<int64_t> (((days * 24 + h) * 60 + mn))) * 60 + s) * 1000


cpdef object to_timestamp(object dt):
    return gmtime(dt)


def array_to_timestamp(ndarray[object, ndim=1] arr):
    cdef int i, n
    cdef ndarray[int64_t, ndim=1] result

    n = len(arr)
    result = np.empty(n, dtype=np.int64)

    for i in range(n):
        result[i] = gmtime(arr[i])

    return result


def time64_to_datetime(ndarray[int64_t, ndim=1] arr):
    cdef int i, n
    cdef ndarray[object, ndim=1] result

    n = len(arr)
    result = np.empty(n, dtype=object)

    for i in range(n):
        result[i] = to_datetime(arr[i])

    return result


# ----------------------------------------------------------------------

cdef inline _to_i8(object val):
    cdef pandas_datetimestruct dts
    try:
        return val.value
    except AttributeError:
        if util.is_datetime64_object(val):
            return get_datetime64_value(val)
        elif PyDateTime_Check(val):
            tzinfo = getattr(val, 'tzinfo', None)
            # Save the original date value so we can get the utcoffset from it.
            ival = _pydatetime_to_dts(val, &dts)
            if tzinfo is not None and not is_utc(tzinfo):
                offset = get_utcoffset(tzinfo, val)
                ival -= int(offset.total_seconds() * 1e9)
            return ival
        return val
