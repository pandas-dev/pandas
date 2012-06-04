# cython: profile=False
cimport numpy as np
import numpy as np

from numpy cimport int32_t, int64_t, import_array, ndarray
from cpython cimport *

# this is our datetime.pxd
from datetime cimport *
from util cimport is_integer_object, is_datetime64_object

from datetime import timedelta
from dateutil.parser import parse as parse_date
cimport util

import cython

# initialize numpy
import_array()

# import datetime C API
PyDateTime_IMPORT

# Period logic
#----------------------------------------------------------------------

cdef inline int64_t apply_mult(int64_t period_ord, int64_t mult):
    """
    Get freq+multiple ordinal value from corresponding freq-only ordinal value.
    For example, 5min ordinal will be 1/5th the 1min ordinal (rounding down to
    integer).
    """
    if mult == 1:
        return period_ord

    return (period_ord - 1) // mult

cdef inline int64_t remove_mult(int64_t period_ord_w_mult, int64_t mult):
    """
    Get freq-only ordinal value from corresponding freq+multiple ordinal.
    """
    if mult == 1:
        return period_ord_w_mult

    return period_ord_w_mult * mult + 1;

def dt64arr_to_periodarr(ndarray[int64_t] dtarr, int freq):
    """
    Convert array of datetime64 values (passed in as 'i8' dtype) to a set of
    periods corresponding to desired frequency, per period convention.
    """
    cdef:
        ndarray[int64_t] out
        Py_ssize_t i, l
        pandas_datetimestruct dts

    l = len(dtarr)

    out = np.empty(l, dtype='i8')

    for i in range(l):
        pandas_datetime_to_datetimestruct(dtarr[i], PANDAS_FR_ns, &dts)
        out[i] = get_period_ordinal(dts.year, dts.month, dts.day,
                                    dts.hour, dts.min, dts.sec, freq)
    return out

def periodarr_to_dt64arr(ndarray[int64_t] periodarr, int freq):
    """
    Convert array to datetime64 values from a set of ordinals corresponding to
    periods per period convention.
    """
    cdef:
        ndarray[int64_t] out
        Py_ssize_t i, l

    l = len(periodarr)

    out = np.empty(l, dtype='i8')

    for i in range(l):
        out[i] = period_ordinal_to_dt64(periodarr[i], freq)

    return out

cdef char START = 'S'
cdef char END = 'E'

cpdef int64_t period_asfreq(int64_t period_ordinal, int freq1, int freq2,
                            bint end):
    """
    Convert period ordinal from one frequency to another, and if upsampling,
    choose to use start ('S') or end ('E') of period.
    """
    cdef:
        int64_t retval

    if end:
        retval = asfreq(period_ordinal, freq1, freq2, END)
    else:
        retval = asfreq(period_ordinal, freq1, freq2, START)

    if retval == INT32_MIN:
        raise ValueError('Frequency conversion failed')

    return retval

def period_asfreq_arr(ndarray[int64_t] arr, int freq1, int freq2, bint end):
    """
    Convert int64-array of period ordinals from one frequency to another, and
    if upsampling, choose to use start ('S') or end ('E') of period.
    """
    cdef:
        ndarray[int64_t] result
        Py_ssize_t i, n
        freq_conv_func func
        asfreq_info finfo
        int64_t val, ordinal
        char relation

    n = len(arr)
    result = np.empty(n, dtype=np.int64)

    func = get_asfreq_func(freq1, freq2)
    get_asfreq_info(freq1, freq2, &finfo)

    if end:
        relation = END
    else:
        relation = START

    for i in range(n):
        val = func(arr[i], relation, &finfo)
        if val == -1:
            raise ValueError("Unable to convert to desired frequency.")
        result[i] = val

    return result

def period_ordinal(int y, int m, int d, int h, int min, int s, int freq):
    cdef:
        int64_t ordinal

    return get_period_ordinal(y, m, d, h, min, s, freq)


cpdef int64_t period_ordinal_to_dt64(int64_t ordinal, int freq):
    cdef:
        pandas_datetimestruct dts
        date_info dinfo

    get_date_info(ordinal, freq, &dinfo)

    dts.year = dinfo.year
    dts.month = dinfo.month
    dts.day = dinfo.day
    dts.hour = dinfo.hour
    dts.min = dinfo.minute
    dts.sec = int(dinfo.second)
    dts.us = dts.ps = 0

    return pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &dts)

def period_ordinal_to_string(int64_t value, int freq):
    cdef:
        char *ptr

    ptr = period_to_string(value, freq)

    if ptr == NULL:
        raise ValueError("Could not create string from ordinal '%s'" % value)

    return <object> ptr

def period_strftime(int64_t value, int freq, object fmt):
    cdef:
        char *ptr

    ptr = period_to_string2(value, freq, <char*>fmt)

    if ptr == NULL:
        raise ValueError("Could not create string with fmt '%s'" % fmt)

    return <object> ptr

# period accessors

ctypedef int (*accessor)(int64_t ordinal, int freq) except INT32_MIN

def get_period_field(int code, int64_t value, int freq):
    cdef accessor f = _get_accessor_func(code)
    return f(value, freq)

def get_period_field_arr(int code, ndarray[int64_t] arr, int freq):
    cdef:
        Py_ssize_t i, sz
        ndarray[int64_t] out
        accessor f

    f = _get_accessor_func(code)

    sz = len(arr)
    out = np.empty(sz, dtype=np.int64)

    for i in range(sz):
        out[i] = f(arr[i], freq)

    return out



cdef accessor _get_accessor_func(int code):
    if code == 0:
        return &pyear
    elif code == 1:
        return &pqyear
    elif code == 2:
        return &pquarter
    elif code == 3:
        return &pmonth
    elif code == 4:
        return &pday
    elif code == 5:
        return &phour
    elif code == 6:
        return &pminute
    elif code == 7:
        return &psecond
    elif code == 8:
        return &pweek
    elif code == 9:
        return &pday_of_year
    elif code == 10:
        return &pweekday
    else:
        raise ValueError('Unrecognized code: %s' % code)

