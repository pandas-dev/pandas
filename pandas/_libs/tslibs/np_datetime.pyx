# -*- coding: utf-8 -*-
# cython: profile=False
cimport cython
from cython cimport Py_ssize_t

from cpython cimport (Py_EQ, Py_NE, Py_GE, Py_GT, Py_LT, Py_LE,
                      PyUnicode_Check, PyUnicode_AsASCIIString)

from cpython.datetime cimport (datetime, date,
                               PyDateTime_IMPORT,
                               PyDateTime_GET_YEAR, PyDateTime_GET_MONTH,
                               PyDateTime_GET_DAY, PyDateTime_DATE_GET_HOUR,
                               PyDateTime_DATE_GET_MINUTE,
                               PyDateTime_DATE_GET_SECOND,
                               PyDateTime_DATE_GET_MICROSECOND)
PyDateTime_IMPORT

import numpy as np
from numpy cimport int64_t

cdef extern from "../src/datetime/np_datetime.h":
    int cmp_pandas_datetimestruct(pandas_datetimestruct *a,
                                  pandas_datetimestruct *b)

    npy_datetime pandas_datetimestruct_to_datetime(PANDAS_DATETIMEUNIT fr,
                                                   pandas_datetimestruct *d
                                                   ) nogil

    void pandas_datetime_to_datetimestruct(npy_datetime val,
                                           PANDAS_DATETIMEUNIT fr,
                                           pandas_datetimestruct *result) nogil

    void pandas_timedelta_to_timedeltastruct(npy_timedelta val,
                                             PANDAS_DATETIMEUNIT fr,
                                             pandas_timedeltastruct *result
                                            ) nogil

    pandas_datetimestruct _NS_MIN_DTS, _NS_MAX_DTS
    void set_datetimestruct_days(int64_t days,
                                 pandas_datetimestruct *dts) nogil


cdef extern from "../src/datetime/np_datetime_strings.h":
    int parse_iso_8601_datetime(char *str, int len,
                                pandas_datetimestruct *out,
                                int *out_local, int *out_tzoffset)

# ----------------------------------------------------------------------
# numpy object inspection

cdef inline npy_datetime get_datetime64_value(object obj) nogil:
    """
    returns the int64 value underlying scalar numpy datetime64 object

    Note that to interpret this as a datetime, the corresponding unit is
    also needed.  That can be found using `get_datetime64_unit`.
    """
    return (<PyDatetimeScalarObject*>obj).obval


cdef inline npy_timedelta get_timedelta64_value(object obj) nogil:
    """
    returns the int64 value underlying scalar numpy timedelta64 object
    """
    return (<PyTimedeltaScalarObject*>obj).obval


cdef inline PANDAS_DATETIMEUNIT get_datetime64_unit(object obj) nogil:
    """
    returns the unit part of the dtype for a numpy datetime64 object.
    """
    return <PANDAS_DATETIMEUNIT>(<PyDatetimeScalarObject*>obj).obmeta.base

# ----------------------------------------------------------------------
# Comparison

cdef int reverse_ops[6]

reverse_ops[Py_LT] = Py_GT
reverse_ops[Py_LE] = Py_GE
reverse_ops[Py_EQ] = Py_EQ
reverse_ops[Py_NE] = Py_NE
reverse_ops[Py_GT] = Py_LT
reverse_ops[Py_GE] = Py_LE


cdef inline bint cmp_scalar(int64_t lhs, int64_t rhs, int op) except -1:
    """
    cmp_scalar is a more performant version of PyObject_RichCompare
    typed for int64_t arguments.
    """
    if op == Py_EQ:
        return lhs == rhs
    elif op == Py_NE:
        return lhs != rhs
    elif op == Py_LT:
        return lhs < rhs
    elif op == Py_LE:
        return lhs <= rhs
    elif op == Py_GT:
        return lhs > rhs
    elif op == Py_GE:
        return lhs >= rhs


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

cdef inline void td64_to_tdstruct(int64_t td64,
                                  pandas_timedeltastruct* out) nogil:
    """Convenience function to call pandas_timedelta_to_timedeltastruct
    with the by-far-most-common frequency PANDAS_FR_ns"""
    pandas_timedelta_to_timedeltastruct(td64, PANDAS_FR_ns, out)
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


cdef inline int _string_to_dts(object val, pandas_datetimestruct* dts,
                               int* out_local, int* out_tzoffset) except? -1:
    cdef:
        int result
        char *tmp

    if PyUnicode_Check(val):
        val = PyUnicode_AsASCIIString(val)

    tmp = val
    result = _cstring_to_dts(tmp, len(val), dts, out_local, out_tzoffset)

    if result == -1:
        raise ValueError('Unable to parse %s' % str(val))
    return result


cdef inline int _cstring_to_dts(char *val, int length,
                                pandas_datetimestruct* dts,
                                int* out_local, int* out_tzoffset) except? -1:
    # Note: without this "extra layer" between _string_to_dts
    # and parse_iso_8601_datetime, calling _string_to_dts raises
    # `SystemError: <class 'str'> returned a result with an error set`
    # in Python3
    cdef:
        int result

    result = parse_iso_8601_datetime(val, length,
                                     dts, out_local, out_tzoffset)
    return result


# ----------------------------------------------------------------------
# Unit Conversion
cdef datetime EPOCH = datetime(1970, 1, 1)

cdef int64_t* _coeffs = [0,                                   # PANDAS_FR_Y
                         0,                                   # PANDAS_FR_M
                         7 * 24 * 3600 * 1000 * 1000 * 1000,  # PANDAS_FR_W
                         0,                                   # NPY_FR_B dummy
                         24 * 3600 * 1000 * 1000 * 1000,      # PANDAS_FR_D
                         3600 * 1000 * 1000 * 1000L,          # PANDAS_FR_h
                         60 * 1000 * 1000 * 1000L,            # PANDAS_FR_m
                         1000 * 1000 * 1000,                  # PANDAS_FR_s
                         1000 * 1000L,                        # PANDAS_FR_ms
                         1000,                                # PANDAS_FR_us
                         1,                                   # PANDAS_FR_ns
                         # From here down we divide instead of multiply
                         1000,                                # PANDAS_FR_ps
                         1000 * 1000,                         # PANDAS_FR_fs
                         1000 * 1000 * 1000]                  # PANDAS_FR_as

# The largest absolute values these can take _without_ raising.
cdef int64_t* _bounds = [292,                             # PANDAS_FR_Y dummy
                         3507,                            # PANDAS_FR_M dummy
                         15250,                           # PANDAS_FR_W
                         0,                               # NPY_FR_B dummy
                         106751,                          # PANDAS_FR_D
                         2562047,                         # PANDAS_FR_h
                         153722867,                       # PANDAS_FR_m
                         9223372036,                      # PANDAS_FR_s
                         9223372036854,                   # PANDAS_FR_ms
                         9223372036854775,                # PANDAS_FR_us
                         9223372036854775807,             # PANDAS_FR_ns
                         9223372036854775807,             # PANDAS_FR_ps
                         9223372036854775807,             # PANDAS_FR_fs
                         9223372036854775807]             # PANDAS_FR_as

# Type names for the np.datetime64 types that are liable to overflow;
# used so we can render the correct exception message
cdef dict type_names = {PANDAS_FR_Y: 'Y', PANDAS_FR_M: 'M', PANDAS_FR_W: 'W',
                        PANDAS_FR_D: 'D', PANDAS_FR_h: 'h', PANDAS_FR_m: 'm',
                        PANDAS_FR_s: 's', PANDAS_FR_ms: 'ms',
                        PANDAS_FR_us: 'us'}


cdef int64_t convert_to_ns(int64_t val, PANDAS_DATETIMEUNIT unit) except? -1:
    """Convert the int64_t representation of a timestamp with the given unit
    to a representation using PANDAS_FR_ns.
    """
    cdef:
        datetime dt
        int64_t year, month
        int64_t coeff, bound

    bound = _bounds[<Py_ssize_t>unit]
    if abs(val) > bound:
        unit_name = type_names[unit]
        val_ns = np.datetime64(val, unit_name).astype('datetime64[ns]')
        fmt = str(val_ns).replace('T', ' ')
        raise OutOfBoundsDatetime('Out of bounds nanosecond timestamp: '
                                  '{fmt}'.format(fmt=fmt))

    if unit == PANDAS_FR_Y:
        dt = datetime(1970 + val, 1, 1)
        return int((dt - EPOCH).total_seconds() * 1e9)

    elif unit == PANDAS_FR_M:
        if val >= 0:
            year = 1970 + val // 12
            month = val % 12 + 1
        else:
            year = 1969 + (val + 1) // 12
            month = 12 + (val + 1) % 12

        dt = datetime(year, month, 1)
        return int((dt - EPOCH).total_seconds() * 1e9)

    elif unit < PANDAS_FR_ns:
        coeff = _coeffs[<Py_ssize_t>unit]
        return val * coeff

    elif unit > PANDAS_FR_ns:
        # no risk of overflows
        coeff = _coeffs[<Py_ssize_t>unit]
        return val // coeff


@cython.cdivision
cdef int convert_datetime_to_dtstruct(int64_t dt, pandas_datetimestruct *out):
    """
    convert a nanosecond (PANDAS_FR_ns) int64_t timestamp to
    a pandas_datetimestruct

    Parameters
    ----------
    dt : int64_t
    out : pandas_datetimestruct*

    Returns
    -------
    code : 0 on success
    """
    cdef:
        int64_t perday = 24LL * 60LL * 60LL * 1000LL * 1000LL * 1000LL

    # Note that care must be taken with the / and % operators
    # for negative values.

    if dt >= 0:
        set_datetimestruct_days(dt / perday, out)
        dt = dt % perday;
    else:
        if dt % perday == 0:
            set_datetimestruct_days(dt / perday - 0, out)
        else:
            set_datetimestruct_days(dt / perday - 1, out)
        dt = (perday - 1) + (dt + 1) % perday

    out.hour = dt / (60 * 60 * 1000000000LL)
    out.min = (dt / (60 * 1000000000LL)) % 60
    out.sec = (dt / 1000000000LL) % 60
    out.us = (dt / 1000LL) % 1000000LL
    out.ps = (dt % 1000LL) * 1000

    return 0
