import cython

from cpython.datetime cimport (
    date,
    datetime,
    time,
    tzinfo,
)

import numpy as np

from numpy cimport (
    int64_t,
    intp_t,
    ndarray,
)

from .conversion cimport normalize_i8_stamp

from .dtypes import Resolution

from .nattype cimport (
    NPY_NAT,
    c_NaT as NaT,
)
from .np_datetime cimport (
    dt64_to_dtstruct,
    npy_datetimestruct,
)
from .offsets cimport BaseOffset
from .period cimport get_period_ordinal
from .timestamps cimport create_timestamp_from_ts
from .tzconversion cimport Localizer

# -------------------------------------------------------------------------

cdef inline object create_datetime_from_ts(
    int64_t value,
    npy_datetimestruct dts,
    tzinfo tz,
    object freq,
    bint fold,
):
    """
    Convenience routine to construct a datetime.datetime from its parts.
    """
    return datetime(
        dts.year, dts.month, dts.day, dts.hour, dts.min, dts.sec, dts.us,
        tz, fold=fold,
    )


cdef inline object create_date_from_ts(
    int64_t value,
    npy_datetimestruct dts,
    tzinfo tz,
    object freq,
    bint fold
):
    """
    Convenience routine to construct a datetime.date from its parts.
    """
    # GH#25057 add fold argument to match other func_create signatures
    return date(dts.year, dts.month, dts.day)


cdef inline object create_time_from_ts(
    int64_t value,
    npy_datetimestruct dts,
    tzinfo tz,
    object freq,
    bint fold
):
    """
    Convenience routine to construct a datetime.time from its parts.
    """
    return time(dts.hour, dts.min, dts.sec, dts.us, tz, fold=fold)


@cython.wraparound(False)
@cython.boundscheck(False)
def ints_to_pydatetime(
    const int64_t[:] stamps,
    tzinfo tz=None,
    BaseOffset freq=None,
    bint fold=False,
    str box="datetime"
) -> np.ndarray:
    """
    Convert an i8 repr to an ndarray of datetimes, date, time or Timestamp.

    Parameters
    ----------
    stamps : array of i8
    tz : str, optional
         convert to this timezone
    freq : BaseOffset, optional
         freq to convert
    fold : bint, default is 0
        Due to daylight saving time, one wall clock time can occur twice
        when shifting from summer to winter time; fold describes whether the
        datetime-like corresponds  to the first (0) or the second time (1)
        the wall clock hits the ambiguous time

        .. versionadded:: 1.1.0
    box : {'datetime', 'timestamp', 'date', 'time'}, default 'datetime'
        * If datetime, convert to datetime.datetime
        * If date, convert to datetime.date
        * If time, convert to datetime.time
        * If Timestamp, convert to pandas.Timestamp

    Returns
    -------
    ndarray[object] of type specified by box
    """
    cdef:
        Py_ssize_t i, n = len(stamps)
        intp_t* pos
        npy_datetimestruct dts
        object new_tz
        int64_t value, local_val
        ndarray[object] result = np.empty(n, dtype=object)
        object (*func_create)(int64_t, npy_datetimestruct, tzinfo, object, bint)
        Localizer info = Localizer(tz)

    if box == "date":
        assert (tz is None), "tz should be None when converting to date"

        func_create = create_date_from_ts
    elif box == "timestamp":
        func_create = create_timestamp_from_ts
    elif box == "time":
        func_create = create_time_from_ts
    elif box == "datetime":
        func_create = create_datetime_from_ts
    else:
        raise ValueError(
            "box must be one of 'datetime', 'date', 'time' or 'timestamp'"
        )

    pos = info.prepare(stamps)

    for i in range(n):
        new_tz = tz
        value = stamps[i]

        if value == NPY_NAT:
            result[i] = <object>NaT
            continue

        local_val = info.utc_val_to_local_val(value, pos, i)

        if info.use_pytz:
            # find right representation of dst etc in pytz timezone
            new_tz = tz._tzinfos[tz._transition_info[pos[i]]]

        dt64_to_dtstruct(local_val, &dts)
        result[i] = func_create(value, dts, new_tz, freq, fold)

    return result


# -------------------------------------------------------------------------

cdef:
    int RESO_NS = 0
    int RESO_US = 1
    int RESO_MS = 2
    int RESO_SEC = 3
    int RESO_MIN = 4
    int RESO_HR = 5
    int RESO_DAY = 6
    int RESO_MTH = 7
    int RESO_QTR = 8
    int RESO_YR = 9


cdef inline int _reso_stamp(npy_datetimestruct *dts):
    if dts.us != 0:
        if dts.us % 1000 == 0:
            return RESO_MS
        return RESO_US
    elif dts.sec != 0:
        return RESO_SEC
    elif dts.min != 0:
        return RESO_MIN
    elif dts.hour != 0:
        return RESO_HR
    return RESO_DAY


@cython.wraparound(False)
@cython.boundscheck(False)
def get_resolution(const int64_t[:] stamps, tzinfo tz=None) -> Resolution:
    cdef:
        Py_ssize_t i, n = len(stamps)
        npy_datetimestruct dts
        int reso = RESO_DAY, curr_reso
        intp_t* pos
        int64_t local_val
        Localizer info = Localizer(tz)

    pos = info.prepare(stamps)

    for i in range(n):
        if stamps[i] == NPY_NAT:
            continue

        local_val = info.utc_val_to_local_val(stamps[i], pos, i)

        dt64_to_dtstruct(local_val, &dts)
        curr_reso = _reso_stamp(&dts)
        if curr_reso < reso:
            reso = curr_reso

    return Resolution(reso)


# -------------------------------------------------------------------------

@cython.wraparound(False)
@cython.boundscheck(False)
cpdef ndarray[int64_t] normalize_i8_timestamps(const int64_t[:] stamps, tzinfo tz):
    """
    Normalize each of the (nanosecond) timezone aware timestamps in the given
    array by rounding down to the beginning of the day (i.e. midnight).
    This is midnight for timezone, `tz`.

    Parameters
    ----------
    stamps : int64 ndarray
    tz : tzinfo or None

    Returns
    -------
    result : int64 ndarray of converted of normalized nanosecond timestamps
    """
    cdef:
        Py_ssize_t i, n = len(stamps)
        int64_t[:] result = np.empty(n, dtype=np.int64)
        intp_t* pos
        int64_t local_val
        Localizer info = Localizer(tz)

    pos = info.prepare(stamps)

    for i in range(n):
        if stamps[i] == NPY_NAT:
            result[i] = NPY_NAT
            continue

        local_val = info.utc_val_to_local_val(stamps[i], pos, i)

        result[i] = normalize_i8_stamp(local_val)

    return result.base  # `.base` to access underlying ndarray


@cython.wraparound(False)
@cython.boundscheck(False)
def is_date_array_normalized(const int64_t[:] stamps, tzinfo tz=None) -> bool:
    """
    Check if all of the given (nanosecond) timestamps are normalized to
    midnight, i.e. hour == minute == second == 0.  If the optional timezone
    `tz` is not None, then this is midnight for this timezone.

    Parameters
    ----------
    stamps : int64 ndarray
    tz : tzinfo or None

    Returns
    -------
    is_normalized : bool True if all stamps are normalized
    """
    cdef:
        Py_ssize_t i, n = len(stamps)
        int64_t local_val
        int64_t day_nanos = 24 * 3600 * 1_000_000_000
        intp_t* pos
        Localizer info = Localizer(tz)

    return info.is_date_array_normalized(stamps)


# -------------------------------------------------------------------------


@cython.wraparound(False)
@cython.boundscheck(False)
def dt64arr_to_periodarr(const int64_t[:] stamps, int freq, tzinfo tz):
    cdef:
        Py_ssize_t i, n = len(stamps)
        int64_t[:] result = np.empty(n, dtype=np.int64)
        intp_t* pos
        npy_datetimestruct dts
        int64_t local_val
        Localizer info = Localizer(tz)

    pos = info.prepare(stamps)

    for i in range(n):
        if stamps[i] == NPY_NAT:
            result[i] = NPY_NAT
            continue

        local_val = info.utc_val_to_local_val(stamps[i], pos, i)

        dt64_to_dtstruct(local_val, &dts)
        result[i] = get_period_ordinal(&dts, freq)

    return result.base  # .base to get underlying ndarray
