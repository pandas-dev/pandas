# cython: profile=False

cimport numpy as np
import numpy as np

from numpy cimport int64_t as i8, import_array, ndarray
from cpython cimport *

from libc.stdlib cimport free

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


cdef extern from "period.h":
    ctypedef struct date_info:
        i8 absdate
        i8 abstime

        i8 nanosecond
        i8 microsecond
        i8 second
        i8 minute
        i8 hour
        i8 day
        i8 month
        i8 quarter
        i8 year
        i8 day_of_week
        i8 day_of_year
        i8 calendar

    ctypedef struct asfreq_info:
        i8 from_week_end
        i8 to_week_end

        i8 from_a_year_end
        i8 to_a_year_end

        i8 from_q_year_end
        i8 to_q_year_end

    ctypedef i8 (*freq_conv_func)(i8 ordinal, char relation,
                                  asfreq_info* af_info)

    i8 asfreq(i8 dtordinal, i8 freq1, i8 freq2, char relation) except INT64_MIN
    freq_conv_func get_asfreq_func(i8 fromFreq, i8 toFreq)
    void get_asfreq_info(i8 fromFreq, i8 toFreq, asfreq_info *af_info)

    i8 get_period_ordinal(i8 year, i8 month, i8 day, i8 hour, i8 minute,
                          i8 second, i8 microsecond, i8 freq) except INT64_MIN

    i8 get_python_ordinal(i8 period_ordinal, i8 freq) except INT64_MIN

    i8 get_date_info(i8 ordinal, i8 freq, date_info *dinfo) except INT64_MIN

    i8 pyear(i8 ordinal, i8 freq) except INT64_MIN
    i8 pqyear(i8 ordinal, i8 freq) except INT64_MIN
    i8 pquarter(i8 ordinal, i8 freq) except INT64_MIN
    i8 pmonth(i8 ordinal, i8 freq) except INT64_MIN
    i8 pday(i8 ordinal, i8 freq) except INT64_MIN
    i8 pweekday(i8 ordinal, i8 freq) except INT64_MIN
    i8 pday_of_week(i8 ordinal, i8 freq) except INT64_MIN
    i8 pday_of_year(i8 ordinal, i8 freq) except INT64_MIN
    i8 pweek(i8 ordinal, i8 freq) except INT64_MIN
    i8 phour(i8 ordinal, i8 freq) except INT64_MIN
    i8 pminute(i8 ordinal, i8 freq) except INT64_MIN
    i8 psecond(i8 ordinal, i8 freq) except INT64_MIN
    i8 pmicrosecond(i8 ordinal, i8 freq) except INT64_MIN
    char *c_strftime(date_info *dinfo, char *fmt) except NULL
    i8 get_yq(i8 ordinal, i8 freq, i8 *quarter, i8 *year) except -1


# Period logic
#----------------------------------------------------------------------
cdef inline i8 apply_mult(i8 period_ord, i8 mult):
    """
    Get freq+multiple ordinal value from corresponding freq-only ordinal value.
    For example, 5min ordinal will be 1/5th the 1min ordinal (rounding down to
    integer).
    """
    if mult == 1:
        return period_ord

    return (period_ord - 1) // mult


cdef inline i8 remove_mult(i8 period_ord_w_mult, i8 mult):
    """
    Get freq-only ordinal value from corresponding freq+multiple ordinal.
    """
    if mult == 1:
        return period_ord_w_mult

    return period_ord_w_mult * mult + 1;


cpdef ndarray[i8] dt64arr_to_periodarr(i8[:] dtarr, i8 freq):
    """
    Convert array of datetime64 values (passed in as 'i8' dtype) to a set of
    periods corresponding to desired frequency, per period convention.
    """
    cdef:
        Py_ssize_t i, l = len(dtarr)
        ndarray[i8] out = np.empty(l, dtype='i8')
        pandas_datetimestruct dts

    for i in range(l):
        pandas_datetime_to_datetimestruct(dtarr[i], PANDAS_FR_ns, &dts)
        out[i] = get_period_ordinal(dts.year, dts.month, dts.day, dts.hour,
                                    dts.min, dts.sec, dts.us, freq)
    return out


cpdef ndarray[i8] periodarr_to_dt64arr(i8[:] periodarr, i8 freq):
    """
    Convert array to datetime64 values from a set of ordinals corresponding to
    periods per period convention.
    """
    cdef:
        Py_ssize_t i, l = len(periodarr)
        ndarray[i8] out = np.empty(l, dtype='i8')

    for i in range(l):
        out[i] = period_ordinal_to_dt64(periodarr[i], freq)

    return out


cdef char START = 'S'
cdef char END = 'E'


cpdef i8 period_asfreq(i8 period_ordinal, i8 freq1, i8 freq2,
                       i8 end) except INT64_MIN:
    """
    Convert period ordinal from one frequency to another, and if upsampling,
    choose to use start ('S') or end ('E') of period.
    """
    cdef:
        char how = END if end else START
        i8 retval = asfreq(period_ordinal, freq1, freq2, how)

    if retval == INT64_MIN:
        raise ValueError('Frequency conversion failed')

    return retval

cpdef ndarray[i8] period_asfreq_arr(i8[:] arr, i8 freq1, i8 freq2, i8 end):
    """
    Convert int64-array of period ordinals from one frequency to another, and
    if upsampling, choose to use start ('S') or end ('E') of period.
    """
    cdef:
        Py_ssize_t i, n = len(arr)
        ndarray[i8] result = np.empty(n, dtype=np.int64)
        freq_conv_func func = get_asfreq_func(freq1, freq2)
        asfreq_info finfo
        i8 val, ordinal
        char relation = END if end else START

    get_asfreq_info(freq1, freq2, &finfo)


    for i in range(n):
        val = func(arr[i], relation, &finfo)

        if val == INT64_MIN:
            raise ValueError("Unable to convert to desired frequency.")

        result[i] = val

    return result


cpdef i8 period_ordinal(i8 y, i8 m, i8 d, i8 h, i8 min, i8 s, i8 us,
                        i8 freq) except INT64_MIN:
    cdef i8 ordinal = get_period_ordinal(y, m, d, h, min, s, us, freq)

    if ordinal == INT64_MIN:
        raise ValueError('Unable to retrieve ordinal')

    return ordinal


cpdef i8 period_ordinal_to_dt64(i8 ordinal, i8 freq) except INT64_MIN:
    cdef:
        pandas_datetimestruct dts
        date_info dinfo

    get_date_info(ordinal, freq, &dinfo)

    dts.year = dinfo.year
    dts.month = dinfo.month
    dts.day = dinfo.day
    dts.hour = dinfo.hour
    dts.min = dinfo.minute
    dts.sec = dinfo.second
    dts.us = dinfo.microsecond
    dts.ps = 0

    return pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &dts)


cpdef str period_format(i8 value, i8 freq, bytes fmt=None):
    cdef i8 freq_group

    if fmt is None:
        freq_group = (freq // 1000) * 1000

        if freq_group == 1000: # FR_ANN
            fmt = b'%Y'
        elif freq_group == 2000: # FR_QTR
            fmt = b'%FQ%q'
        elif freq_group == 3000: # FR_MTH
            fmt = b'%Y-%m'
        elif freq_group == 4000: # WK
            left = period_asfreq(value, freq, 6000, 0)
            right = period_asfreq(value, freq, 6000, 1)
            return '%s/%s' % (period_format(left, 6000),
                              period_format(right, 6000))
        elif (freq_group == 5000 # BUS
              or freq_group == 6000): # DAY
            fmt = b'%Y-%m-%d'
        elif freq_group == 7000: # HR
            fmt = b'%Y-%m-%d %H:00'
        elif freq_group == 8000: # MIN
            fmt = b'%Y-%m-%d %H:%M'
        elif freq_group == 9000: # SEC
            fmt = b'%Y-%m-%d %H:%M:%S'
        elif freq_group == 11000:
            fmt = b'%Y-%m-%d %H:%M:%S.%%06u'
        else:
            raise ValueError('Unknown freq: %d' % freq)

    return _period_strftime(value, freq, fmt)


cdef list extra_fmts = [(b"%q", b"^`AB`^"),
                        (b"%f", b"^`CD`^"),
                        (b"%F", b"^`EF`^")]

cdef list str_extra_fmts = ["^`AB`^", "^`CD`^", "^`EF`^"]


cdef str _period_strftime(i8 value, i8 freq, bytes fmt):
    cdef:
        Py_ssize_t i
        date_info dinfo
        char *formatted
        object pat, repl, result
        list found_pat = [False] * len(extra_fmts)
        i8 year, quarter

    if PyUnicode_Check(fmt):
        fmt = fmt.encode('utf-8')

    get_date_info(value, freq, &dinfo)
    for i in range(len(extra_fmts)):
        pat = extra_fmts[i][0]
        repl = extra_fmts[i][1]

        if pat in fmt:
            fmt = fmt.replace(pat, repl)
            found_pat[i] = True

    formatted = c_strftime(&dinfo, <char*> fmt)

    result = util.char_to_string(formatted)
    free(formatted)

    for i in range(len(extra_fmts)):

        if found_pat[i]:

            if get_yq(value, freq, &quarter, &year) < 0:
                raise ValueError('Unable to get quarter and year')

            if i == 0:
                repl = '%d' % quarter
            elif i == 1:  # %f, 2-digit year
                repl = '%.2d' % (year % 100)
            elif i == 2:
                repl = '%d' % year

            result = result.replace(str_extra_fmts[i], repl)

    # Py3?
    if not PyString_Check(result):
        result = str(result)

    return result

# period accessors
ctypedef i8 (*accessor)(i8 ordinal, i8 freq) except INT64_MIN


cpdef i8 get_period_field(i8 code, i8 value, i8 freq) except INT64_MIN:
    cdef:
        accessor f = _get_accessor_func(code)
        i8 r = f(value, freq)

    if code == 11:
        print 'ordinal:', value
        print 'us:     ', r

    if r == INT64_MIN:
        raise ValueError('Unable to retrieve property')

    return r


cpdef ndarray[i8] get_period_field_arr(i8 code, i8[:] arr, i8 freq):
    cdef:
        Py_ssize_t i, sz
        i8 v
        ndarray[i8] out
        accessor f

    f = _get_accessor_func(code)

    sz = len(arr)
    out = np.empty(sz, dtype=np.int64)

    for i in range(sz):
        v = f(arr[i], freq)

        if v == INT64_MIN:
            raise ValueError('Cannot get freq %i' % freq)

        out[i] = v

    return out


cdef accessor _get_accessor_func(i8 code):
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
    elif code == 11:
        return &pmicrosecond
    else:
        raise ValueError('Unrecognized code: %s' % code)
