import cython

import numpy as np
cimport numpy as cnp
from numpy cimport int64_t, int32_t, intp_t, ndarray
cnp.import_array()

import pytz

# stdlib datetime imports
from datetime import time as datetime_time
from cpython.datetime cimport (datetime, tzinfo,
                               PyDateTime_Check, PyDate_Check,
                               PyDateTime_IMPORT)
PyDateTime_IMPORT

from pandas._libs.tslibs.c_timestamp cimport _Timestamp

from pandas._libs.tslibs.np_datetime cimport (
    check_dts_bounds, npy_datetimestruct, pandas_datetime_to_datetimestruct,
    _string_to_dts, npy_datetime, dt64_to_dtstruct, dtstruct_to_dt64,
    get_datetime64_unit, get_datetime64_value, pydatetime_to_dt64,
    NPY_DATETIMEUNIT, NPY_FR_ns)
from pandas._libs.tslibs.np_datetime import OutOfBoundsDatetime

from pandas._libs.tslibs.util cimport (
    is_datetime64_object, is_integer_object, is_float_object)

from pandas._libs.tslibs.timedeltas cimport cast_from_unit
from pandas._libs.tslibs.timezones cimport (
    is_utc, is_tzlocal, is_fixed_offset, get_utcoffset, get_dst_info,
    get_timezone, maybe_get_tz, tz_compare)
from pandas._libs.tslibs.timezones import UTC
from pandas._libs.tslibs.parsing import parse_datetime_string

from pandas._libs.tslibs.nattype import nat_strings
from pandas._libs.tslibs.nattype cimport (
    NPY_NAT, checknull_with_nat, c_NaT as NaT)

from pandas._libs.tslibs.tzconversion import (
    tz_localize_to_utc, tz_convert_single)
from pandas._libs.tslibs.tzconversion cimport (
    _tz_convert_tzlocal_utc, _tz_convert_tzlocal_fromutc)

# ----------------------------------------------------------------------
# Constants

NS_DTYPE = np.dtype('M8[ns]')
TD_DTYPE = np.dtype('m8[ns]')


# ----------------------------------------------------------------------
# Misc Helpers

cdef inline int64_t get_datetime64_nanos(object val) except? -1:
    """
    Extract the value and unit from a np.datetime64 object, then convert the
    value to nanoseconds if necessary.
    """
    cdef:
        npy_datetimestruct dts
        NPY_DATETIMEUNIT unit
        npy_datetime ival

    ival = get_datetime64_value(val)
    if ival == NPY_NAT:
        return NPY_NAT

    unit = get_datetime64_unit(val)

    if unit != NPY_FR_ns:
        pandas_datetime_to_datetimestruct(ival, unit, &dts)
        check_dts_bounds(&dts)
        ival = dtstruct_to_dt64(&dts)

    return ival


@cython.boundscheck(False)
@cython.wraparound(False)
def ensure_datetime64ns(arr: ndarray, copy: bool=True):
    """
    Ensure a np.datetime64 array has dtype specifically 'datetime64[ns]'

    Parameters
    ----------
    arr : ndarray
    copy : bool, default True

    Returns
    -------
    ndarray with dtype datetime64[ns]
    """
    cdef:
        Py_ssize_t i, n = arr.size
        int64_t[:] ivalues, iresult
        NPY_DATETIMEUNIT unit
        npy_datetimestruct dts

    shape = (<object>arr).shape

    if (<object>arr).dtype.byteorder == ">":
        # GH#29684 we incorrectly get OutOfBoundsDatetime if we dont swap
        dtype = arr.dtype
        arr = arr.astype(dtype.newbyteorder("<"))

    ivalues = arr.view(np.int64).ravel()

    result = np.empty(shape, dtype=NS_DTYPE)
    iresult = result.ravel().view(np.int64)

    if len(iresult) == 0:
        result = arr.view(NS_DTYPE)
        if copy:
            result = result.copy()
        return result

    unit = get_datetime64_unit(arr.flat[0])
    if unit == NPY_FR_ns:
        if copy:
            arr = arr.copy()
        result = arr
    else:
        for i in range(n):
            if ivalues[i] != NPY_NAT:
                pandas_datetime_to_datetimestruct(ivalues[i], unit, &dts)
                iresult[i] = dtstruct_to_dt64(&dts)
                check_dts_bounds(&dts)
            else:
                iresult[i] = NPY_NAT

    return result


def ensure_timedelta64ns(arr: ndarray, copy: bool=True):
    """
    Ensure a np.timedelta64 array has dtype specifically 'timedelta64[ns]'

    Parameters
    ----------
    arr : ndarray
    copy : boolean, default True

    Returns
    -------
    result : ndarray with dtype timedelta64[ns]

    """
    return arr.astype(TD_DTYPE, copy=copy)
    # TODO: check for overflows when going from a lower-resolution to nanos


@cython.boundscheck(False)
@cython.wraparound(False)
def datetime_to_datetime64(ndarray[object] values):
    """
    Convert ndarray of datetime-like objects to int64 array representing
    nanosecond timestamps.

    Parameters
    ----------
    values : ndarray[object]

    Returns
    -------
    result : ndarray[int64_t]
    inferred_tz : tzinfo or None
    """
    cdef:
        Py_ssize_t i, n = len(values)
        object val, inferred_tz = None
        int64_t[:] iresult
        npy_datetimestruct dts
        _TSObject _ts
        bint found_naive = False

    result = np.empty(n, dtype='M8[ns]')
    iresult = result.view('i8')
    for i in range(n):
        val = values[i]
        if checknull_with_nat(val):
            iresult[i] = NPY_NAT
        elif PyDateTime_Check(val):
            if val.tzinfo is not None:
                if found_naive:
                    raise ValueError('Cannot mix tz-aware with '
                                     'tz-naive values')
                if inferred_tz is not None:
                    if not tz_compare(val.tzinfo, inferred_tz):
                        raise ValueError('Array must be all same time zone')
                else:
                    inferred_tz = get_timezone(val.tzinfo)

                _ts = convert_datetime_to_tsobject(val, None)
                iresult[i] = _ts.value
                check_dts_bounds(&_ts.dts)
            else:
                found_naive = True
                if inferred_tz is not None:
                    raise ValueError('Cannot mix tz-aware with '
                                     'tz-naive values')
                iresult[i] = pydatetime_to_dt64(val, &dts)
                check_dts_bounds(&dts)
        else:
            raise TypeError(f'Unrecognized value type: {type(val)}')

    return result, inferred_tz


# ----------------------------------------------------------------------
# _TSObject Conversion

# lightweight C object to hold datetime & int64 pair
cdef class _TSObject:
    # cdef:
    #    npy_datetimestruct dts      # npy_datetimestruct
    #    int64_t value               # numpy dt64
    #    object tzinfo
    #    bint fold

    def __cinit__(self):
        # GH 25057. As per PEP 495, set fold to 0 by default
        self.fold = 0

    @property
    def value(self):
        # This is needed in order for `value` to be accessible in lib.pyx
        return self.value


cdef convert_to_tsobject(object ts, object tz, object unit,
                         bint dayfirst, bint yearfirst, int32_t nanos=0):
    """
    Extract datetime and int64 from any of:
        - np.int64 (with unit providing a possible modifier)
        - np.datetime64
        - a float (with unit providing a possible modifier)
        - python int or long object (with unit providing a possible modifier)
        - iso8601 string object
        - python datetime object
        - another timestamp object

    Raises
    ------
    OutOfBoundsDatetime : ts cannot be converted within implementation bounds
    """
    cdef:
        _TSObject obj

    if tz is not None:
        tz = maybe_get_tz(tz)

    obj = _TSObject()

    if isinstance(ts, str):
        return convert_str_to_tsobject(ts, tz, unit, dayfirst, yearfirst)

    if ts is None or ts is NaT:
        obj.value = NPY_NAT
    elif is_datetime64_object(ts):
        obj.value = get_datetime64_nanos(ts)
        if obj.value != NPY_NAT:
            dt64_to_dtstruct(obj.value, &obj.dts)
    elif is_integer_object(ts):
        try:
            ts = <int64_t>ts
        except OverflowError:
            # GH#26651 re-raise as OutOfBoundsDatetime
            raise OutOfBoundsDatetime(ts)
        if ts == NPY_NAT:
            obj.value = NPY_NAT
        else:
            ts = ts * cast_from_unit(None, unit)
            obj.value = ts
            dt64_to_dtstruct(ts, &obj.dts)
    elif is_float_object(ts):
        if ts != ts or ts == NPY_NAT:
            obj.value = NPY_NAT
        else:
            ts = cast_from_unit(ts, unit)
            obj.value = ts
            dt64_to_dtstruct(ts, &obj.dts)
    elif PyDateTime_Check(ts):
        return convert_datetime_to_tsobject(ts, tz, nanos)
    elif PyDate_Check(ts):
        # Keep the converter same as PyDateTime's
        ts = datetime.combine(ts, datetime_time())
        return convert_datetime_to_tsobject(ts, tz)
    elif getattr(ts, '_typ', None) == 'period':
        raise ValueError("Cannot convert Period to Timestamp "
                         "unambiguously. Use to_timestamp")
    else:
        raise TypeError(f'Cannot convert input [{ts}] of type {type(ts)} to '
                        f'Timestamp')

    if tz is not None:
        localize_tso(obj, tz)

    if obj.value != NPY_NAT:
        # check_overflows needs to run after localize_tso
        check_dts_bounds(&obj.dts)
        check_overflows(obj)
    return obj


cdef _TSObject convert_datetime_to_tsobject(datetime ts, object tz,
                                            int32_t nanos=0):
    """
    Convert a datetime (or Timestamp) input `ts`, along with optional timezone
    object `tz` to a _TSObject.

    The optional argument `nanos` allows for cases where datetime input
    needs to be supplemented with higher-precision information.

    Parameters
    ----------
    ts : datetime or Timestamp
        Value to be converted to _TSObject
    tz : tzinfo or None
        timezone for the timezone-aware output
    nanos : int32_t, default is 0
        nanoseconds supplement the precision of the datetime input ts

    Returns
    -------
    obj : _TSObject
    """
    cdef:
        _TSObject obj = _TSObject()

    obj.fold = ts.fold
    if tz is not None:
        tz = maybe_get_tz(tz)

        if ts.tzinfo is not None:
            # Convert the current timezone to the passed timezone
            ts = ts.astimezone(tz)
            obj.value = pydatetime_to_dt64(ts, &obj.dts)
            obj.tzinfo = ts.tzinfo
        elif not is_utc(tz):
            ts = _localize_pydatetime(ts, tz)
            obj.value = pydatetime_to_dt64(ts, &obj.dts)
            obj.tzinfo = ts.tzinfo
        else:
            # UTC
            obj.value = pydatetime_to_dt64(ts, &obj.dts)
            obj.tzinfo = tz
    else:
        obj.value = pydatetime_to_dt64(ts, &obj.dts)
        obj.tzinfo = ts.tzinfo

    if obj.tzinfo is not None and not is_utc(obj.tzinfo):
        offset = get_utcoffset(obj.tzinfo, ts)
        obj.value -= int(offset.total_seconds() * 1e9)

    if isinstance(ts, _Timestamp):
        obj.value += ts.nanosecond
        obj.dts.ps = ts.nanosecond * 1000

    if nanos:
        obj.value += nanos
        obj.dts.ps = nanos * 1000

    check_dts_bounds(&obj.dts)
    check_overflows(obj)
    return obj


cdef _TSObject create_tsobject_tz_using_offset(npy_datetimestruct dts,
                                               int tzoffset, object tz=None):
    """
    Convert a datetimestruct `dts`, along with initial timezone offset
    `tzoffset` to a _TSObject (with timezone object `tz` - optional).

    Parameters
    ----------
    dts: npy_datetimestruct
    tzoffset: int
    tz : tzinfo or None
        timezone for the timezone-aware output.

    Returns
    -------
    obj : _TSObject
    """
    cdef:
        _TSObject obj = _TSObject()
        int64_t value  # numpy dt64
        datetime dt
        ndarray[int64_t] trans
        int64_t[:] deltas

    value = dtstruct_to_dt64(&dts)
    obj.dts = dts
    obj.tzinfo = pytz.FixedOffset(tzoffset)
    obj.value = tz_convert_single(value, obj.tzinfo, UTC)
    if tz is None:
        check_overflows(obj)
        return obj

    # Infer fold from offset-adjusted obj.value
    # see PEP 495 https://www.python.org/dev/peps/pep-0495/#the-fold-attribute
    if is_utc(tz):
        pass
    elif is_tzlocal(tz):
        _tz_convert_tzlocal_fromutc(obj.value, tz, &obj.fold)
    else:
        trans, deltas, typ = get_dst_info(tz)

        if typ == 'dateutil':
            pos = trans.searchsorted(obj.value, side='right') - 1
            obj.fold = _infer_tsobject_fold(obj, trans, deltas, pos)

    # Keep the converter same as PyDateTime's
    dt = datetime(obj.dts.year, obj.dts.month, obj.dts.day,
                  obj.dts.hour, obj.dts.min, obj.dts.sec,
                  obj.dts.us, obj.tzinfo, fold=obj.fold)
    obj = convert_datetime_to_tsobject(
        dt, tz, nanos=obj.dts.ps // 1000)
    return obj


cdef _TSObject convert_str_to_tsobject(object ts, object tz, object unit,
                                       bint dayfirst=False,
                                       bint yearfirst=False):
    """
    Convert a string input `ts`, along with optional timezone object`tz`
    to a _TSObject.

    The optional arguments `dayfirst` and `yearfirst` are passed to the
    dateutil parser.

    Parameters
    ----------
    ts : str
        Value to be converted to _TSObject
    tz : tzinfo or None
        timezone for the timezone-aware output
    dayfirst : bool, default False
        When parsing an ambiguous date string, interpret e.g. "3/4/1975" as
        April 3, as opposed to the standard US interpretation March 4.
    yearfirst : bool, default False
        When parsing an ambiguous date string, interpret e.g. "01/05/09"
        as "May 9, 2001", as opposed to the default "Jan 5, 2009"

    Returns
    -------
    obj : _TSObject
    """
    cdef:
        npy_datetimestruct dts
        int out_local = 0, out_tzoffset = 0
        bint do_parse_datetime_string = False

    if tz is not None:
        tz = maybe_get_tz(tz)

    assert isinstance(ts, str)

    if len(ts) == 0 or ts in nat_strings:
        ts = NaT
    elif ts == 'now':
        # Issue 9000, we short-circuit rather than going
        # into np_datetime_strings which returns utc
        ts = datetime.now(tz)
    elif ts == 'today':
        # Issue 9000, we short-circuit rather than going
        # into np_datetime_strings which returns a normalized datetime
        ts = datetime.now(tz)
        # equiv: datetime.today().replace(tzinfo=tz)
    else:
        string_to_dts_failed = _string_to_dts(
            ts, &dts, &out_local,
            &out_tzoffset, False
        )
        try:
            if not string_to_dts_failed:
                check_dts_bounds(&dts)
                if out_local == 1:
                    return create_tsobject_tz_using_offset(dts,
                                                           out_tzoffset, tz)
                else:
                    ts = dtstruct_to_dt64(&dts)
                    if tz is not None:
                        # shift for localize_tso
                        ts = tz_localize_to_utc(np.array([ts], dtype='i8'), tz,
                                                ambiguous='raise')[0]

        except OutOfBoundsDatetime:
            # GH#19382 for just-barely-OutOfBounds falling back to dateutil
            # parser will return incorrect result because it will ignore
            # nanoseconds
            raise

        except ValueError:
            do_parse_datetime_string = True

        if string_to_dts_failed or do_parse_datetime_string:
            try:
                ts = parse_datetime_string(ts, dayfirst=dayfirst,
                                           yearfirst=yearfirst)
            except (ValueError, OverflowError):
                raise ValueError("could not convert string to Timestamp")

    return convert_to_tsobject(ts, tz, unit, dayfirst, yearfirst)


cdef inline check_overflows(_TSObject obj):
    """
    Check that we haven't silently overflowed in timezone conversion

    Parameters
    ----------
    obj : _TSObject

    Returns
    -------
    None

    Raises
    ------
    OutOfBoundsDatetime
    """
    # GH#12677
    if obj.dts.year == 1677:
        if not (obj.value < 0):
            raise OutOfBoundsDatetime
    elif obj.dts.year == 2262:
        if not (obj.value > 0):
            raise OutOfBoundsDatetime


# ----------------------------------------------------------------------
# Localization

cdef inline void localize_tso(_TSObject obj, tzinfo tz):
    """
    Given the UTC nanosecond timestamp in obj.value, find the wall-clock
    representation of that timestamp in the given timezone.

    Parameters
    ----------
    obj : _TSObject
    tz : tzinfo

    Returns
    -------
    None

    Notes
    -----
    Sets obj.tzinfo inplace, alters obj.dts inplace.
    """
    cdef:
        ndarray[int64_t] trans
        int64_t[:] deltas
        int64_t local_val
        Py_ssize_t pos
        str typ

    assert obj.tzinfo is None

    if is_utc(tz):
        pass
    elif obj.value == NPY_NAT:
        pass
    elif is_tzlocal(tz):
        local_val = _tz_convert_tzlocal_fromutc(obj.value, tz, &obj.fold)
        dt64_to_dtstruct(local_val, &obj.dts)
    else:
        # Adjust datetime64 timestamp, recompute datetimestruct
        trans, deltas, typ = get_dst_info(tz)

        if is_fixed_offset(tz):
            # static/fixed tzinfo; in this case we know len(deltas) == 1
            # This can come back with `typ` of either "fixed" or None
            dt64_to_dtstruct(obj.value + deltas[0], &obj.dts)
        elif typ == 'pytz':
            # i.e. treat_tz_as_pytz(tz)
            pos = trans.searchsorted(obj.value, side='right') - 1
            tz = tz._tzinfos[tz._transition_info[pos]]
            dt64_to_dtstruct(obj.value + deltas[pos], &obj.dts)
        elif typ == 'dateutil':
            # i.e. treat_tz_as_dateutil(tz)
            pos = trans.searchsorted(obj.value, side='right') - 1
            dt64_to_dtstruct(obj.value + deltas[pos], &obj.dts)
            # dateutil supports fold, so we infer fold from value
            obj.fold = _infer_tsobject_fold(obj, trans, deltas, pos)
        else:
            # Note: as of 2018-07-17 all tzinfo objects that are _not_
            # either pytz or dateutil have is_fixed_offset(tz) == True,
            # so this branch will never be reached.
            pass

    obj.tzinfo = tz


cdef inline bint _infer_tsobject_fold(_TSObject obj, ndarray[int64_t] trans,
                                      int64_t[:] deltas, int32_t pos):
    """
    Infer _TSObject fold property from value by assuming 0 and then setting
    to 1 if necessary.

    Parameters
    ----------
    obj : _TSObject
    trans : ndarray[int64_t]
        ndarray of offset transition points in nanoseconds since epoch.
    deltas : int64_t[:]
        array of offsets corresponding to transition points in trans.
    pos : int32_t
        Position of the last transition point before taking fold into account.

    Returns
    -------
    bint
        Due to daylight saving time, one wall clock time can occur twice
        when shifting from summer to winter time; fold describes whether the
        datetime-like corresponds  to the first (0) or the second time (1)
        the wall clock hits the ambiguous time

    References
    ----------
    .. [1] "PEP 495 - Local Time Disambiguation"
           https://www.python.org/dev/peps/pep-0495/#the-fold-attribute
    """
    cdef:
        bint fold = 0

    if pos > 0:
        fold_delta = deltas[pos - 1] - deltas[pos]
        if obj.value - fold_delta < trans[pos]:
            fold = 1

    return fold

cdef inline datetime _localize_pydatetime(datetime dt, tzinfo tz):
    """
    Take a datetime/Timestamp in UTC and localizes to timezone tz.

    NB: Unlike the public version, this treats datetime and Timestamp objects
        identically, i.e. discards nanos from Timestamps.
        It also assumes that the `tz` input is not None.
    """
    try:
        # datetime.replace with pytz may be incorrect result
        return tz.localize(dt)
    except AttributeError:
        return dt.replace(tzinfo=tz)


cpdef inline datetime localize_pydatetime(datetime dt, object tz):
    """
    Take a datetime/Timestamp in UTC and localizes to timezone tz.

    Parameters
    ----------
    dt : datetime or Timestamp
    tz : tzinfo, "UTC", or None

    Returns
    -------
    localized : datetime or Timestamp
    """
    if tz is None:
        return dt
    elif isinstance(dt, _Timestamp):
        return dt.tz_localize(tz)
    elif is_utc(tz):
        return _localize_pydatetime(dt, tz)
    try:
        # datetime.replace with pytz may be incorrect result
        return tz.localize(dt)
    except AttributeError:
        return dt.replace(tzinfo=tz)


# ----------------------------------------------------------------------
# Normalization


def normalize_date(dt: object) -> datetime:
    """
    Normalize datetime.datetime value to midnight. Returns datetime.date as a
    datetime.datetime at midnight

    Parameters
    ----------
    dt : date, datetime, or Timestamp

    Returns
    -------
    normalized : datetime.datetime or Timestamp

    Raises
    ------
    TypeError : if input is not datetime.date, datetime.datetime, or Timestamp
    """
    if PyDateTime_Check(dt):
        if isinstance(dt, _Timestamp):
            return dt.replace(hour=0, minute=0, second=0, microsecond=0,
                              nanosecond=0)
        else:
            # regular datetime object
            return dt.replace(hour=0, minute=0, second=0, microsecond=0)
            # TODO: Make sure DST crossing is handled correctly here
    elif PyDate_Check(dt):
        return datetime(dt.year, dt.month, dt.day)
    else:
        raise TypeError(f'Unrecognized type: {type(dt)}')


@cython.wraparound(False)
@cython.boundscheck(False)
def normalize_i8_timestamps(int64_t[:] stamps, object tz):
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
        Py_ssize_t n = len(stamps)
        int64_t[:] result = np.empty(n, dtype=np.int64)

    result = _normalize_local(stamps, tz)

    return result.base  # .base to access underlying np.ndarray


@cython.wraparound(False)
@cython.boundscheck(False)
cdef int64_t[:] _normalize_local(int64_t[:] stamps, tzinfo tz):
    """
    Normalize each of the (nanosecond) timestamps in the given array by
    rounding down to the beginning of the day (i.e. midnight) for the
    given timezone `tz`.

    Parameters
    ----------
    stamps : int64 ndarray
    tz : tzinfo

    Returns
    -------
    result : int64 ndarray of converted of normalized nanosecond timestamps
    """
    cdef:
        Py_ssize_t i, n = len(stamps)
        int64_t[:] result = np.empty(n, dtype=np.int64)
        ndarray[int64_t] trans
        int64_t[:] deltas
        str typ
        Py_ssize_t[:] pos
        npy_datetimestruct dts
        int64_t delta, local_val

    if is_tzlocal(tz):
        for i in range(n):
            if stamps[i] == NPY_NAT:
                result[i] = NPY_NAT
                continue
            local_val = _tz_convert_tzlocal_utc(stamps[i], tz, to_utc=False)
            dt64_to_dtstruct(local_val, &dts)
            result[i] = _normalized_stamp(&dts)
    else:
        # Adjust datetime64 timestamp, recompute datetimestruct
        trans, deltas, typ = get_dst_info(tz)

        if typ not in ['pytz', 'dateutil']:
            # static/fixed; in this case we know that len(delta) == 1
            delta = deltas[0]
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                dt64_to_dtstruct(stamps[i] + delta, &dts)
                result[i] = _normalized_stamp(&dts)
        else:
            pos = trans.searchsorted(stamps, side='right') - 1
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                dt64_to_dtstruct(stamps[i] + deltas[pos[i]], &dts)
                result[i] = _normalized_stamp(&dts)

    return result


cdef inline int64_t _normalized_stamp(npy_datetimestruct *dts) nogil:
    """
    Normalize the given datetimestruct to midnight, then convert to int64_t.

    Parameters
    ----------
    *dts : pointer to npy_datetimestruct

    Returns
    -------
    stamp : int64
    """
    dts.hour = 0
    dts.min = 0
    dts.sec = 0
    dts.us = 0
    dts.ps = 0
    return dtstruct_to_dt64(dts)


@cython.wraparound(False)
@cython.boundscheck(False)
def is_date_array_normalized(int64_t[:] stamps, object tz=None):
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
        ndarray[int64_t] trans
        int64_t[:] deltas
        intp_t[:] pos
        npy_datetimestruct dts
        int64_t local_val, delta
        str typ

    if tz is None or is_utc(tz):
        for i in range(n):
            dt64_to_dtstruct(stamps[i], &dts)
            if (dts.hour + dts.min + dts.sec + dts.us) > 0:
                return False
    elif is_tzlocal(tz):
        for i in range(n):
            local_val = _tz_convert_tzlocal_utc(stamps[i], tz, to_utc=False)
            dt64_to_dtstruct(local_val, &dts)
            if (dts.hour + dts.min + dts.sec + dts.us) > 0:
                return False
    else:
        trans, deltas, typ = get_dst_info(tz)

        if typ not in ['pytz', 'dateutil']:
            # static/fixed; in this case we know that len(delta) == 1
            delta = deltas[0]
            for i in range(n):
                # Adjust datetime64 timestamp, recompute datetimestruct
                dt64_to_dtstruct(stamps[i] + delta, &dts)
                if (dts.hour + dts.min + dts.sec + dts.us) > 0:
                    return False

        else:
            pos = trans.searchsorted(stamps) - 1
            for i in range(n):
                # Adjust datetime64 timestamp, recompute datetimestruct
                dt64_to_dtstruct(stamps[i] + deltas[pos[i]], &dts)
                if (dts.hour + dts.min + dts.sec + dts.us) > 0:
                    return False

    return True
