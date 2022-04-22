import cython

from cpython.datetime cimport (
    date,
    datetime,
    time,
    tzinfo,
)

import numpy as np

cimport numpy as cnp
from numpy cimport (
    int64_t,
    intp_t,
    ndarray,
)

cnp.import_array()

from .dtypes import Resolution

from .ccalendar cimport DAY_NANOS
from .dtypes cimport c_Resolution
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
    is_zoneinfo,
)
from .tzconversion cimport (
    bisect_right_i8,
    localize_tzinfo_api,
)


cdef const int64_t[::1] _deltas_placeholder = np.array([], dtype=np.int64)


@cython.freelist(16)
@cython.internal
@cython.final
cdef class Localizer:
    cdef:
        tzinfo tz
        bint use_utc, use_fixed, use_tzlocal, use_dst, use_pytz
        ndarray trans
        Py_ssize_t ntrans
        const int64_t[::1] deltas
        int64_t delta

    @cython.initializedcheck(False)
    @cython.boundscheck(False)
    def __cinit__(self, tzinfo tz):
        self.tz = tz
        self.use_utc = self.use_tzlocal = self.use_fixed = False
        self.use_dst = self.use_pytz = False
        self.ntrans = -1  # placeholder
        self.delta = -1  # placeholder
        self.deltas = _deltas_placeholder

        if is_utc(tz) or tz is None:
            self.use_utc = True

        elif is_tzlocal(tz) or is_zoneinfo(tz):
            self.use_tzlocal = True

        else:
            trans, deltas, typ = get_dst_info(tz)
            self.trans = trans
            self.ntrans = trans.shape[0]
            self.deltas = deltas

            if typ != "pytz" and typ != "dateutil":
                # static/fixed; in this case we know that len(delta) == 1
                self.use_fixed = True
                self.delta = deltas[0]
            else:
                self.use_dst = True
                if typ == "pytz":
                    self.use_pytz = True


# -------------------------------------------------------------------------


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
        int64_t utc_val, local_val
        Py_ssize_t pos, i, n = stamps.shape[0]
        int64_t* tdata = NULL

        npy_datetimestruct dts
        tzinfo new_tz
        ndarray[object] result = np.empty(n, dtype=object)
        bint use_date = False, use_time = False, use_ts = False, use_pydt = False

    if box == "date":
        assert (tz is None), "tz should be None when converting to date"
        use_date = True
    elif box == "timestamp":
        use_ts = True
    elif box == "time":
        use_time = True
    elif box == "datetime":
        use_pydt = True
    else:
        raise ValueError(
            "box must be one of 'datetime', 'date', 'time' or 'timestamp'"
        )

    if info.use_dst:
        tdata = <int64_t*>cnp.PyArray_DATA(info.trans)

    for i in range(n):
        utc_val = stamps[i]
        new_tz = tz

        if utc_val == NPY_NAT:
            result[i] = <object>NaT
            continue

        if info.use_utc:
            local_val = utc_val
        elif info.use_tzlocal:
            local_val = utc_val + localize_tzinfo_api(utc_val, tz)
        elif info.use_fixed:
            local_val = utc_val + info.delta
        else:
            pos = bisect_right_i8(tdata, utc_val, info.ntrans) - 1
            local_val = utc_val + info.deltas[pos]

            if info.use_pytz:
                # find right representation of dst etc in pytz timezone
                new_tz = tz._tzinfos[tz._transition_info[pos]]

        dt64_to_dtstruct(local_val, &dts)

        if use_ts:
            result[i] = create_timestamp_from_ts(utc_val, dts, new_tz, freq, fold)
        elif use_pydt:
            result[i] = datetime(
                dts.year, dts.month, dts.day, dts.hour, dts.min, dts.sec, dts.us,
                new_tz, fold=fold,
            )
        elif use_date:
            result[i] = date(dts.year, dts.month, dts.day)
        else:
            result[i] = time(dts.hour, dts.min, dts.sec, dts.us, new_tz, fold=fold)

    return result


# -------------------------------------------------------------------------


cdef inline c_Resolution _reso_stamp(npy_datetimestruct *dts):
    if dts.us != 0:
        if dts.us % 1000 == 0:
            return c_Resolution.RESO_MS
        return c_Resolution.RESO_US
    elif dts.sec != 0:
        return c_Resolution.RESO_SEC
    elif dts.min != 0:
        return c_Resolution.RESO_MIN
    elif dts.hour != 0:
        return c_Resolution.RESO_HR
    return c_Resolution.RESO_DAY


@cython.wraparound(False)
@cython.boundscheck(False)
def get_resolution(const int64_t[:] stamps, tzinfo tz=None) -> Resolution:
    cdef:
        Localizer info = Localizer(tz)
        int64_t utc_val, local_val
        Py_ssize_t pos, i, n = stamps.shape[0]
        int64_t* tdata = NULL

        npy_datetimestruct dts
        c_Resolution reso = c_Resolution.RESO_DAY, curr_reso

    if info.use_dst:
        tdata = <int64_t*>cnp.PyArray_DATA(info.trans)

    for i in range(n):
        utc_val = stamps[i]
        if utc_val == NPY_NAT:
            continue

        if info.use_utc:
            local_val = utc_val
        elif info.use_tzlocal:
            local_val = utc_val + localize_tzinfo_api(utc_val, tz)
        elif info.use_fixed:
            local_val = utc_val + info.delta
        else:
            pos = bisect_right_i8(tdata, utc_val, info.ntrans) - 1
            local_val = utc_val + info.deltas[pos]

        dt64_to_dtstruct(local_val, &dts)
        curr_reso = _reso_stamp(&dts)
        if curr_reso < reso:
            reso = curr_reso

    return Resolution(reso)


# -------------------------------------------------------------------------


@cython.cdivision(False)
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
        int64_t utc_val, local_val
        Py_ssize_t pos, i, n = stamps.shape[0]
        int64_t* tdata = NULL

        int64_t[::1] result = np.empty(n, dtype=np.int64)

    if info.use_dst:
        tdata = <int64_t*>cnp.PyArray_DATA(info.trans)

    for i in range(n):
        utc_val = stamps[i]
        if utc_val == NPY_NAT:
            result[i] = NPY_NAT
            continue

        if info.use_utc:
            local_val = utc_val
        elif info.use_tzlocal:
            local_val = utc_val + localize_tzinfo_api(utc_val, tz)
        elif info.use_fixed:
            local_val = utc_val + info.delta
        else:
            pos = bisect_right_i8(tdata, utc_val, info.ntrans) - 1
            local_val = utc_val + info.deltas[pos]

        result[i] = local_val - (local_val % DAY_NANOS)

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
        int64_t utc_val, local_val
        Py_ssize_t pos, i, n = stamps.shape[0]
        int64_t* tdata = NULL

    if info.use_dst:
        tdata = <int64_t*>cnp.PyArray_DATA(info.trans)

    for i in range(n):
        utc_val = stamps[i]
        if info.use_utc:
            local_val = utc_val
        elif info.use_tzlocal:
            local_val = utc_val + localize_tzinfo_api(utc_val, tz)
        elif info.use_fixed:
            local_val = utc_val + info.delta
        else:
            pos = bisect_right_i8(tdata, utc_val, info.ntrans) - 1
            local_val = utc_val + info.deltas[pos]

        if local_val % DAY_NANOS != 0:
            return False

    return True


# -------------------------------------------------------------------------


@cython.wraparound(False)
@cython.boundscheck(False)
def dt64arr_to_periodarr(ndarray stamps, int freq, tzinfo tz):
    # stamps is int64_t, arbitrary ndim
    cdef:
        Localizer info = Localizer(tz)
        Py_ssize_t pos, i, n = stamps.size
        int64_t utc_val, local_val, res_val
        int64_t* tdata = NULL

        npy_datetimestruct dts
        ndarray result = cnp.PyArray_EMPTY(stamps.ndim, stamps.shape, cnp.NPY_INT64, 0)
        cnp.broadcast mi = cnp.PyArray_MultiIterNew2(result, stamps)

    if info.use_dst:
        tdata = <int64_t*>cnp.PyArray_DATA(info.trans)

    for i in range(n):
        # Analogous to: utc_val = stamps[i]
        utc_val = (<int64_t*>cnp.PyArray_MultiIter_DATA(mi, 1))[0]

        if utc_val == NPY_NAT:
            res_val = NPY_NAT
        else:
            if info.use_utc:
                local_val = utc_val
            elif info.use_tzlocal:
                local_val = utc_val + localize_tzinfo_api(utc_val, tz)
            elif info.use_fixed:
                local_val = utc_val + info.delta
            else:
                pos = bisect_right_i8(tdata, utc_val, info.ntrans) - 1
                local_val = utc_val + info.deltas[pos]

            dt64_to_dtstruct(local_val, &dts)
            res_val = get_period_ordinal(&dts, freq)

        # Analogous to: result[i] = res_val
        (<int64_t*>cnp.PyArray_MultiIter_DATA(mi, 0))[0] = res_val

        cnp.PyArray_MultiIter_NEXT(mi)

    return result
