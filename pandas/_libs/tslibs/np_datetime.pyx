# -*- coding: utf-8 -*-
# cython: profile=False

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

from numpy cimport int64_t

cdef extern from "../src/datetime/np_datetime.h":
    int cmp_npy_datetimestruct(npy_datetimestruct *a,
                               npy_datetimestruct *b)

    npy_datetime npy_datetimestruct_to_datetime(NPY_DATETIMEUNIT fr,
                                                npy_datetimestruct *d) nogil

    void pandas_datetime_to_datetimestruct(npy_datetime val,
                                           NPY_DATETIMEUNIT fr,
                                           npy_datetimestruct *result) nogil

    void pandas_timedelta_to_timedeltastruct(npy_timedelta val,
                                             NPY_DATETIMEUNIT fr,
                                             pandas_timedeltastruct *result
                                            ) nogil

    npy_datetimestruct _NS_MIN_DTS, _NS_MAX_DTS

cdef extern from "../src/datetime/np_datetime_strings.h":
    int parse_iso_8601_datetime(char *str, int len,
                                npy_datetimestruct *out,
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


cdef inline NPY_DATETIMEUNIT get_datetime64_unit(object obj) nogil:
    """
    returns the unit part of the dtype for a numpy datetime64 object.
    """
    return <NPY_DATETIMEUNIT>(<PyDatetimeScalarObject*>obj).obmeta.base

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


cdef inline check_dts_bounds(npy_datetimestruct *dts):
    """Raises OutOfBoundsDatetime if the given date is outside the range that
    can be represented by nanosecond-resolution 64-bit integers."""
    cdef:
        bint error = False

    if (dts.year <= 1677 and
            cmp_npy_datetimestruct(dts, &_NS_MIN_DTS) == -1):
        error = True
    elif (dts.year >= 2262 and
          cmp_npy_datetimestruct(dts, &_NS_MAX_DTS) == 1):
        error = True

    if error:
        fmt = '%d-%.2d-%.2d %.2d:%.2d:%.2d' % (dts.year, dts.month,
                                               dts.day, dts.hour,
                                               dts.min, dts.sec)
        raise OutOfBoundsDatetime(
            'Out of bounds nanosecond timestamp: {fmt}'.format(fmt=fmt))


# ----------------------------------------------------------------------
# Conversion

cdef inline int64_t dtstruct_to_dt64(npy_datetimestruct* dts) nogil:
    """Convenience function to call npy_datetimestruct_to_datetime
    with the by-far-most-common frequency NPY_FR_ns"""
    return npy_datetimestruct_to_datetime(NPY_FR_ns, dts)


cdef inline void dt64_to_dtstruct(int64_t dt64,
                                  npy_datetimestruct* out) nogil:
    """Convenience function to call pandas_datetime_to_datetimestruct
    with the by-far-most-common frequency NPY_FR_ns"""
    pandas_datetime_to_datetimestruct(dt64, NPY_FR_ns, out)
    return

cdef inline void td64_to_tdstruct(int64_t td64,
                                  pandas_timedeltastruct* out) nogil:
    """Convenience function to call pandas_timedelta_to_timedeltastruct
    with the by-far-most-common frequency NPY_FR_ns"""
    pandas_timedelta_to_timedeltastruct(td64, NPY_FR_ns, out)
    return


cdef inline int64_t pydatetime_to_dt64(datetime val,
                                       npy_datetimestruct *dts):
    """
    Note we are assuming that the datetime object is timezone-naive.
    """
    dts.year = PyDateTime_GET_YEAR(val)
    dts.month = PyDateTime_GET_MONTH(val)
    dts.day = PyDateTime_GET_DAY(val)
    dts.hour = PyDateTime_DATE_GET_HOUR(val)
    dts.min = PyDateTime_DATE_GET_MINUTE(val)
    dts.sec = PyDateTime_DATE_GET_SECOND(val)
    dts.us = PyDateTime_DATE_GET_MICROSECOND(val)
    dts.ps = dts.as = 0
    return dtstruct_to_dt64(dts)


cdef inline int64_t pydate_to_dt64(date val, npy_datetimestruct *dts):
    dts.year = PyDateTime_GET_YEAR(val)
    dts.month = PyDateTime_GET_MONTH(val)
    dts.day = PyDateTime_GET_DAY(val)
    dts.hour = dts.min = dts.sec = dts.us = 0
    dts.ps = dts.as = 0
    return dtstruct_to_dt64(dts)


cdef inline int _string_to_dts(object val, npy_datetimestruct* dts,
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
                                npy_datetimestruct* dts,
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
