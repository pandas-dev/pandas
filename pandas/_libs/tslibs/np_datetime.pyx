from cpython.datetime cimport (
    PyDateTime_DATE_GET_HOUR,
    PyDateTime_DATE_GET_MICROSECOND,
    PyDateTime_DATE_GET_MINUTE,
    PyDateTime_DATE_GET_SECOND,
    PyDateTime_GET_DAY,
    PyDateTime_GET_MONTH,
    PyDateTime_GET_YEAR,
    import_datetime,
)
from cpython.object cimport (
    Py_EQ,
    Py_GE,
    Py_GT,
    Py_LE,
    Py_LT,
    Py_NE,
)

import_datetime()

cimport numpy as cnp

cnp.import_array()
from numpy cimport (
    int64_t,
    ndarray,
)

from pandas._libs.tslibs.util cimport get_c_string_buf_and_size


cdef extern from "src/datetime/np_datetime.h":
    int cmp_npy_datetimestruct(npy_datetimestruct *a,
                               npy_datetimestruct *b)

    # AS, FS, PS versions exist but are not imported because they are not used.
    npy_datetimestruct _NS_MIN_DTS, _NS_MAX_DTS
    npy_datetimestruct _US_MIN_DTS, _US_MAX_DTS
    npy_datetimestruct _MS_MIN_DTS, _MS_MAX_DTS
    npy_datetimestruct _S_MIN_DTS, _S_MAX_DTS
    npy_datetimestruct _M_MIN_DTS, _M_MAX_DTS

    PyArray_DatetimeMetaData get_datetime_metadata_from_dtype(cnp.PyArray_Descr *dtype);

cdef extern from "src/datetime/np_datetime_strings.h":
    int parse_iso_8601_datetime(const char *str, int len, int want_exc,
                                npy_datetimestruct *out,
                                NPY_DATETIMEUNIT *out_bestunit,
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


cdef NPY_DATETIMEUNIT get_unit_from_dtype(cnp.dtype dtype):
    # NB: caller is responsible for ensuring this is *some* datetime64 or
    #  timedelta64 dtype, otherwise we can segfault
    cdef:
        cnp.PyArray_Descr* descr = <cnp.PyArray_Descr*>dtype
        PyArray_DatetimeMetaData meta
    meta = get_datetime_metadata_from_dtype(descr)
    return meta.base


def py_get_unit_from_dtype(dtype):
    # for testing get_unit_from_dtype; adds 896 bytes to the .so file.
    return get_unit_from_dtype(dtype)


# ----------------------------------------------------------------------
# Comparison


cdef bint cmp_dtstructs(
    npy_datetimestruct* left, npy_datetimestruct* right, int op
):
    cdef:
        int cmp_res

    cmp_res = cmp_npy_datetimestruct(left, right)
    if op == Py_EQ:
        return cmp_res == 0
    if op == Py_NE:
        return cmp_res != 0
    if op == Py_GT:
        return cmp_res == 1
    if op == Py_LT:
        return cmp_res == -1
    if op == Py_GE:
        return cmp_res == 1 or cmp_res == 0
    else:
        # i.e. op == Py_LE
        return cmp_res == -1 or cmp_res == 0


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
    """
    Raised when the datetime is outside the range that
    can be represented.
    """
    pass


class OutOfBoundsTimedelta(ValueError):
    """
    Raised when encountering a timedelta value that cannot be represented
    as a timedelta64[ns].
    """
    # Timedelta analogue to OutOfBoundsDatetime
    pass


cdef check_dts_bounds(npy_datetimestruct *dts, NPY_DATETIMEUNIT unit=NPY_FR_ns):
    """Raises OutOfBoundsDatetime if the given date is outside the range that
    can be represented by nanosecond-resolution 64-bit integers."""
    cdef:
        bint error = False
        npy_datetimestruct cmp_upper, cmp_lower

    if unit == NPY_FR_ns:
        cmp_upper = _NS_MAX_DTS
        cmp_lower = _NS_MIN_DTS
    elif unit == NPY_FR_us:
        cmp_upper = _US_MAX_DTS
        cmp_lower = _US_MIN_DTS
    elif unit == NPY_FR_ms:
        cmp_upper = _MS_MAX_DTS
        cmp_lower = _MS_MIN_DTS
    elif unit == NPY_FR_s:
        cmp_upper = _S_MAX_DTS
        cmp_lower = _S_MIN_DTS
    elif unit == NPY_FR_m:
        cmp_upper = _M_MAX_DTS
        cmp_lower = _M_MIN_DTS
    else:
        raise NotImplementedError(unit)

    if cmp_npy_datetimestruct(dts, &cmp_lower) == -1:
        error = True
    elif cmp_npy_datetimestruct(dts, &cmp_upper) == 1:
        error = True

    if error:
        fmt = (f'{dts.year}-{dts.month:02d}-{dts.day:02d} '
               f'{dts.hour:02d}:{dts.min:02d}:{dts.sec:02d}')
        # TODO: "nanosecond" in the message assumes NPY_FR_ns
        raise OutOfBoundsDatetime(f'Out of bounds nanosecond timestamp: {fmt}')


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


# just exposed for testing at the moment
def py_td64_to_tdstruct(int64_t td64, NPY_DATETIMEUNIT unit):
    cdef:
        pandas_timedeltastruct tds
    pandas_timedelta_to_timedeltastruct(td64, unit, &tds)
    return tds  # <- returned as a dict to python


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


cdef inline void pydate_to_dtstruct(date val, npy_datetimestruct *dts):
    dts.year = PyDateTime_GET_YEAR(val)
    dts.month = PyDateTime_GET_MONTH(val)
    dts.day = PyDateTime_GET_DAY(val)
    dts.hour = dts.min = dts.sec = dts.us = 0
    dts.ps = dts.as = 0
    return

cdef inline int64_t pydate_to_dt64(date val, npy_datetimestruct *dts):
    pydate_to_dtstruct(val, dts)
    return dtstruct_to_dt64(dts)


cdef inline int string_to_dts(
    str val,
    npy_datetimestruct* dts,
    NPY_DATETIMEUNIT* out_bestunit,
    int* out_local,
    int* out_tzoffset,
    bint want_exc,
) except? -1:
    cdef:
        Py_ssize_t length
        const char* buf

    buf = get_c_string_buf_and_size(val, &length)
    return parse_iso_8601_datetime(buf, length, want_exc,
                                   dts, out_bestunit, out_local, out_tzoffset)


cpdef ndarray astype_overflowsafe(
    ndarray values,
    cnp.dtype dtype,
    bint copy=True,
):
    """
    Convert an ndarray with datetime64[X] to datetime64[Y], raising on overflow.
    """
    if values.descr.type_num != cnp.NPY_DATETIME:
        # aka values.dtype.kind != "M"
        raise TypeError("astype_overflowsafe values must have datetime64 dtype")
    if dtype.type_num != cnp.NPY_DATETIME:
        raise TypeError("astype_overflowsafe dtype must be datetime64")

    cdef:
        NPY_DATETIMEUNIT from_unit = get_unit_from_dtype(values.dtype)
        NPY_DATETIMEUNIT to_unit = get_unit_from_dtype(dtype)

    if (
        from_unit == NPY_DATETIMEUNIT.NPY_FR_GENERIC
        or to_unit == NPY_DATETIMEUNIT.NPY_FR_GENERIC
    ):
        # without raising explicitly here, we end up with a SystemError
        # built-in function [...] returned a result with an error
        raise ValueError("datetime64 values and dtype must have a unit specified")

    if from_unit == to_unit:
        # Check this before allocating result for perf, might save some memory
        if copy:
            return values.copy()
        return values

    cdef:
        ndarray i8values = values.view("i8")

        # equiv: result = np.empty((<object>values).shape, dtype="i8")
        ndarray iresult = cnp.PyArray_EMPTY(
            values.ndim, values.shape, cnp.NPY_INT64, 0
        )

        cnp.broadcast mi = cnp.PyArray_MultiIterNew2(iresult, i8values)
        Py_ssize_t i, N = values.size
        int64_t value, new_value
        npy_datetimestruct dts

    for i in range(N):
        # Analogous to: item = values[i]
        value = (<int64_t*>cnp.PyArray_MultiIter_DATA(mi, 1))[0]

        if value == NPY_DATETIME_NAT:
            new_value = NPY_DATETIME_NAT
        else:
            pandas_datetime_to_datetimestruct(value, from_unit, &dts)
            check_dts_bounds(&dts, to_unit)
            new_value = npy_datetimestruct_to_datetime(to_unit, &dts)

        # Analogous to: iresult[i] = new_value
        (<int64_t*>cnp.PyArray_MultiIter_DATA(mi, 0))[0] = new_value

        cnp.PyArray_MultiIter_NEXT(mi)

    return iresult.view(dtype)
