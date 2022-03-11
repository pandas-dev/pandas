import cython

cimport numpy as cnp
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

cnp.import_array()

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
from .timezones cimport (
    get_dst_info,
    is_tzlocal,
    is_utc,
)
from .tzconversion cimport tz_convert_utc_to_tzlocal


cdef const int64_t[::1] _deltas_placeholder = np.array([], dtype=np.int64)

ctypedef int64_t (*localizer_func)(Localizer, int64_t, intp_t*, intp_t)


@cython.freelist(16)
@cython.final
cdef class Localizer:
    cdef readonly:
        tzinfo tz
        bint use_utc, use_fixed, use_tzlocal, use_dst, use_pytz
        ndarray trans
        const int64_t[::1] deltas
        int64_t delta
        str typ

    # TODO: report cython bug; this declaration works if on its own line
    #  but raises at compile-time if included in the 'cdef readonly'
    #  declarations above.
    # cdef int64_t (*func)(Localizer, int64_t, intp_t*, intp_t)
    # cdef localizer_func func

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    def __cinit__(self, tzinfo tz):
        self.tz = tz
        self.use_utc = self.use_tzlocal = self.use_fixed = self.use_dst = self.use_pytz = False
        self.delta = -1  # placeholder
        self.deltas = _deltas_placeholder

        if is_utc(tz) or tz is None:
            self.use_utc = True
            # self.func = self.func_use_utc

        elif is_tzlocal(tz):
            self.use_tzlocal = True
            # self.func = self.func_use_tzlocal

        else:
            trans, deltas, typ = get_dst_info(tz)
            self.trans = trans
            self.deltas = deltas
            self.typ = typ

            if typ not in ["pytz", "dateutil"]:
                # static/fixed; in this case we know that len(delta) == 1
                self.use_fixed = True
                self.delta = deltas[0]
                # self.func = self.func_use_fixed
            else:
                self.use_dst = True
                # self.func = self.func_use_dst
                if typ == "pytz":
                    self.use_pytz = True

    cdef int64_t func_use_utc(self, int64_t utc_val, intp_t* pos, intp_t i):
        return utc_val

    cdef int64_t func_use_tzlocal(self, int64_t utc_val, intp_t* pos, intp_t i):
        return tz_convert_utc_to_tzlocal(utc_val, self.tz)

    cdef int64_t func_use_fixed(self, int64_t utc_val, intp_t* pos, intp_t i):
        return utc_val + self.delta

    cdef int64_t func_use_dst(self, int64_t utc_val, intp_t* pos, intp_t i):
        return utc_val + self.deltas[pos[i]]


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
        Localizer info = Localizer(tz)
        const int64_t[::1] deltas = info.deltas
        int64_t value, local_val, delta = info.delta
        intp_t* pos
        Py_ssize_t i, n = len(stamps)
        npy_datetimestruct dts

        tzinfo new_tz
        ndarray[object] result = np.empty(n, dtype=object)
        object (*func_create)(int64_t, npy_datetimestruct, tzinfo, object, bint)

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

    if info.use_dst:
        pos = <intp_t*>cnp.PyArray_DATA(info.trans.searchsorted(stamps, side="right") - 1)

    for i in range(n):
        new_tz = tz
        value = stamps[i]

        if value == NPY_NAT:
            result[i] = <object>NaT
            continue

        if info.use_utc:
            local_val = value
        elif info.use_tzlocal:
            local_val = tz_convert_utc_to_tzlocal(value, tz)
        elif info.use_fixed:
            local_val = value + delta
        else:
            local_val = value + info.deltas[pos[i]]

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


def get_resolution(const int64_t[:] stamps, tzinfo tz=None) -> Resolution:
    cdef:
        Localizer info = Localizer(tz)
        const int64_t[::1] deltas = info.deltas
        int64_t local_val, delta = info.delta
        intp_t* pos
        Py_ssize_t i, n = len(stamps)
        npy_datetimestruct dts

        int reso = RESO_DAY, curr_reso

    if info.use_dst:
        pos = <intp_t*>cnp.PyArray_DATA(info.trans.searchsorted(stamps, side="right") - 1)

    for i in range(n):
        if stamps[i] == NPY_NAT:
            continue

        if info.use_utc:
            local_val = stamps[i]
        elif info.use_tzlocal:
            local_val = tz_convert_utc_to_tzlocal(stamps[i], tz)
        elif info.use_fixed:
            local_val = stamps[i] + info.delta
        else:
            local_val = stamps[i] + deltas[pos[i]]

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
        Localizer info = Localizer(tz)
        const int64_t[::1] deltas = info.deltas
        int64_t local_val, delta = info.delta
        intp_t* pos
        Py_ssize_t i, n = len(stamps)

        int64_t[:] result = np.empty(n, dtype=np.int64)

    if info.use_dst:
        pos = <intp_t*>cnp.PyArray_DATA(info.trans.searchsorted(stamps, side="right") - 1)

    for i in range(n):
        if stamps[i] == NPY_NAT:
            result[i] = NPY_NAT
            continue

        if info.use_utc:
            local_val = stamps[i]
        elif info.use_tzlocal:
            local_val = tz_convert_utc_to_tzlocal(stamps[i], tz)
        elif info.use_fixed:
            local_val = stamps[i] + delta
        else:
            local_val = stamps[i] + info.deltas[pos[i]]

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
        Localizer info = Localizer(tz)
        const int64_t[::1] deltas = info.deltas
        int64_t local_val, delta = info.delta
        intp_t* pos
        Py_ssize_t i, n = len(stamps)

        int64_t day_nanos = 24 * 3600 * 1_000_000_000

    if info.use_dst:
        pos = <intp_t*>cnp.PyArray_DATA(info.trans.searchsorted(stamps, side="right") - 1)

    for i in range(n):
        if info.use_utc:
            local_val = stamps[i]
        elif info.use_tzlocal:
            local_val = tz_convert_utc_to_tzlocal(stamps[i], tz)
        elif info.use_fixed:
            local_val = stamps[i] + info.delta
        else:
            local_val = stamps[i] + deltas[pos[i]]

        if local_val % day_nanos != 0:
            return False

    return True


# -------------------------------------------------------------------------


@cython.wraparound(False)
@cython.boundscheck(False)
def dt64arr_to_periodarr(const int64_t[:] stamps, int freq, tzinfo tz):
    cdef:
        Localizer info = Localizer(tz)
        const int64_t[::1] deltas = info.deltas
        int64_t local_val, delta = info.delta
        intp_t* pos
        Py_ssize_t i, n = len(stamps)
        npy_datetimestruct dts

        int64_t[:] result = np.empty(n, dtype=np.int64)

    if info.use_dst:
        pos = <intp_t*>cnp.PyArray_DATA(info.trans.searchsorted(stamps, side="right") - 1)

    for i in range(n):
        if stamps[i] == NPY_NAT:
            result[i] = NPY_NAT
            continue

        if info.use_utc:
            local_val = stamps[i]
        elif info.use_tzlocal:
            local_val = tz_convert_utc_to_tzlocal(stamps[i], tz)
        elif info.use_fixed:
            local_val = stamps[i] + delta
        else:
            local_val = stamps[i] + info.deltas[pos[i]]

        dt64_to_dtstruct(local_val, &dts)
        result[i] = get_period_ordinal(&dts, freq)

    return result.base  # .base to get underlying ndarray
