# -*- coding: utf-8 -*-
# cython: profile=False

cimport cython
from cython cimport Py_ssize_t

import numpy as np
cimport numpy as cnp
from numpy cimport int64_t, int32_t, ndarray
cnp.import_array()

import pytz

# stdlib datetime imports
from datetime import time as datetime_time
from cpython.datetime cimport (datetime, tzinfo,
                               PyDateTime_Check, PyDate_Check,
                               PyDateTime_CheckExact, PyDateTime_IMPORT)
PyDateTime_IMPORT

from np_datetime cimport (check_dts_bounds,
                          npy_datetimestruct,
                          pandas_datetime_to_datetimestruct, _string_to_dts,
                          npy_datetime,
                          dt64_to_dtstruct, dtstruct_to_dt64,
                          get_datetime64_unit, get_datetime64_value,
                          pydatetime_to_dt64, NPY_DATETIMEUNIT, NPY_FR_ns)
from np_datetime import OutOfBoundsDatetime

from util cimport (is_string_object,
                   is_datetime64_object,
                   is_integer_object, is_float_object, is_array)

from timedeltas cimport cast_from_unit
from timezones cimport (is_utc, is_tzlocal, is_fixed_offset,
                        treat_tz_as_dateutil, treat_tz_as_pytz,
                        get_utcoffset, get_dst_info,
                        get_timezone, maybe_get_tz, tz_compare)
from parsing import parse_datetime_string

from nattype import nat_strings, NaT
from nattype cimport NPY_NAT, checknull_with_nat

# ----------------------------------------------------------------------
# Constants

cdef int64_t DAY_NS = 86400000000000LL
NS_DTYPE = np.dtype('M8[ns]')
TD_DTYPE = np.dtype('m8[ns]')

UTC = pytz.UTC

# ----------------------------------------------------------------------
# Misc Helpers

# TODO: How to declare np.datetime64 as the input type?
cdef inline int64_t get_datetime64_nanos(object val) except? -1:
    """
    Extract the value and unit from a np.datetime64 object, then convert the
    value to nanoseconds if necessary.
    """
    cdef:
        npy_datetimestruct dts
        NPY_DATETIMEUNIT unit
        npy_datetime ival

    unit = get_datetime64_unit(val)
    ival = get_datetime64_value(val)

    if unit != NPY_FR_ns:
        pandas_datetime_to_datetimestruct(ival, unit, &dts)
        check_dts_bounds(&dts)
        ival = dtstruct_to_dt64(&dts)

    return ival


def ensure_datetime64ns(ndarray arr, copy=True):
    """
    Ensure a np.datetime64 array has dtype specifically 'datetime64[ns]'

    Parameters
    ----------
    arr : ndarray
    copy : boolean, default True

    Returns
    -------
    result : ndarray with dtype datetime64[ns]

    """
    cdef:
        Py_ssize_t i, n = arr.size
        int64_t[:] ivalues, iresult
        NPY_DATETIMEUNIT unit
        npy_datetimestruct dts

    shape = (<object> arr).shape

    ivalues = arr.view(np.int64).ravel()

    result = np.empty(shape, dtype='M8[ns]')
    iresult = result.ravel().view(np.int64)

    if len(iresult) == 0:
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


def ensure_timedelta64ns(ndarray arr, copy=True):
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


def datetime_to_datetime64(object[:] values):
    """
    Convert ndarray of datetime-like objects to int64 array representing
    nanosecond timestamps.

    Parameters
    ----------
    values : ndarray

    Returns
    -------
    result : ndarray with dtype int64
    inferred_tz : tzinfo or None
    """
    cdef:
        Py_ssize_t i, n = len(values)
        object val, inferred_tz = None
        int64_t[:] iresult
        npy_datetimestruct dts
        _TSObject _ts

    result = np.empty(n, dtype='M8[ns]')
    iresult = result.view('i8')
    for i in range(n):
        val = values[i]
        if checknull_with_nat(val):
            iresult[i] = NPY_NAT
        elif PyDateTime_Check(val):
            if val.tzinfo is not None:
                if inferred_tz is not None:
                    if not tz_compare(val.tzinfo, inferred_tz):
                        raise ValueError('Array must be all same time zone')
                else:
                    inferred_tz = get_timezone(val.tzinfo)

                _ts = convert_datetime_to_tsobject(val, None)
                iresult[i] = _ts.value
                check_dts_bounds(&_ts.dts)
            else:
                if inferred_tz is not None:
                    raise ValueError('Cannot mix tz-aware with '
                                     'tz-naive values')
                iresult[i] = pydatetime_to_dt64(val, &dts)
                check_dts_bounds(&dts)
        else:
            raise TypeError('Unrecognized value type: %s' % type(val))

    return result, inferred_tz


cdef inline maybe_datetimelike_to_i8(object val):
    """
    Try to convert to a nanosecond timestamp.  Fall back to returning the
    input value.

    Parameters
    ----------
    val : object

    Returns
    -------
    val : int64 timestamp or original input
    """
    cdef:
        npy_datetimestruct dts
    try:
        return val.value
    except AttributeError:
        if is_datetime64_object(val):
            return get_datetime64_value(val)
        elif PyDateTime_Check(val):
            return convert_datetime_to_tsobject(val, None).value
        return val


# ----------------------------------------------------------------------
# _TSObject Conversion

# lightweight C object to hold datetime & int64 pair
cdef class _TSObject:
    # cdef:
    #    npy_datetimestruct dts      # npy_datetimestruct
    #    int64_t value               # numpy dt64
    #    object tzinfo

    @property
    def value(self):
        return self.value


cpdef int64_t pydt_to_i8(object pydt) except? -1:
    """
    Convert to int64 representation compatible with numpy datetime64; converts
    to UTC

    Parameters
    ----------
    pydt : object

    Returns
    -------
    i8value : np.int64
    """
    cdef:
        _TSObject ts

    ts = convert_to_tsobject(pydt, None, None, 0, 0)

    return ts.value


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
    """
    cdef:
        _TSObject obj

    if tz is not None:
        tz = maybe_get_tz(tz)

    obj = _TSObject()

    if is_string_object(ts):
        return convert_str_to_tsobject(ts, tz, unit, dayfirst, yearfirst)

    if ts is None or ts is NaT:
        obj.value = NPY_NAT
    elif is_datetime64_object(ts):
        if ts.view('i8') == NPY_NAT:
            obj.value = NPY_NAT
        else:
            obj.value = get_datetime64_nanos(ts)
            dt64_to_dtstruct(obj.value, &obj.dts)
    elif is_integer_object(ts):
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
        raise TypeError('Cannot convert input [{}] of type {} to '
                        'Timestamp'.format(ts, type(ts)))

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
            obj.tzinfo = pytz.utc
    else:
        obj.value = pydatetime_to_dt64(ts, &obj.dts)
        obj.tzinfo = ts.tzinfo

    if obj.tzinfo is not None and not is_utc(obj.tzinfo):
        offset = get_utcoffset(obj.tzinfo, ts)
        obj.value -= int(offset.total_seconds() * 1e9)

    if not PyDateTime_CheckExact(ts):
        # datetime instance but not datetime type --> Timestamp
        obj.value += ts.nanosecond
        obj.dts.ps = ts.nanosecond * 1000

    if nanos:
        obj.value += nanos
        obj.dts.ps = nanos * 1000

    check_dts_bounds(&obj.dts)
    check_overflows(obj)
    return obj


cdef _TSObject convert_str_to_tsobject(object ts, object tz, object unit,
                                       bint dayfirst=False,
                                       bint yearfirst=False):
    """
    Convert a string-like (bytes or unicode) input `ts`, along with optional
    timezone object `tz` to a _TSObject.

    The optional arguments `dayfirst` and `yearfirst` are passed to the
    dateutil parser.

    Parameters
    ----------
    ts : bytes or unicode
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
        _TSObject obj
        int out_local = 0, out_tzoffset = 0
        datetime dt

    if tz is not None:
        tz = maybe_get_tz(tz)

    obj = _TSObject()

    assert is_string_object(ts)

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
        try:
            _string_to_dts(ts, &obj.dts, &out_local, &out_tzoffset)
            obj.value = dtstruct_to_dt64(&obj.dts)
            check_dts_bounds(&obj.dts)
            if out_local == 1:
                obj.tzinfo = pytz.FixedOffset(out_tzoffset)
                obj.value = tz_convert_single(obj.value, obj.tzinfo, 'UTC')
                if tz is None:
                    check_dts_bounds(&obj.dts)
                    check_overflows(obj)
                    return obj
                else:
                    # Keep the converter same as PyDateTime's
                    obj = convert_to_tsobject(obj.value, obj.tzinfo,
                                              None, 0, 0)
                    dt = datetime(obj.dts.year, obj.dts.month, obj.dts.day,
                                  obj.dts.hour, obj.dts.min, obj.dts.sec,
                                  obj.dts.us, obj.tzinfo)
                    obj = convert_datetime_to_tsobject(dt, tz,
                                                       nanos=obj.dts.ps / 1000)
                    return obj

            else:
                ts = obj.value
                if tz is not None:
                    # shift for localize_tso
                    ts = tz_localize_to_utc(np.array([ts], dtype='i8'), tz,
                                            ambiguous='raise',
                                            errors='raise')[0]

        except OutOfBoundsDatetime:
            # GH#19382 for just-barely-OutOfBounds falling back to dateutil
            # parser will return incorrect result because it will ignore
            # nanoseconds
            raise

        except ValueError:
            try:
                ts = parse_datetime_string(ts, dayfirst=dayfirst,
                                           yearfirst=yearfirst)
            except Exception:
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

    assert obj.tzinfo is None

    if is_utc(tz):
        pass
    elif obj.value == NPY_NAT:
        pass
    elif is_tzlocal(tz):
        local_val = _tz_convert_tzlocal_utc(obj.value, tz, to_utc=False)
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
        else:
            # Note: as of 2018-07-17 all tzinfo objects that are _not_
            # either pytz or dateutil have is_fixed_offset(tz) == True,
            # so this branch will never be reached.
            pass

    obj.tzinfo = tz


cdef inline datetime _localize_pydatetime(datetime dt, tzinfo tz):
    """
    Take a datetime/Timestamp in UTC and localizes to timezone tz.

    NB: Unlike the public version, this treats datetime and Timestamp objects
        identically, i.e. discards nanos from Timestamps.
        It also assumes that the `tz` input is not None.
    """
    if tz == 'UTC' or tz is UTC:
        return UTC.localize(dt)
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
    elif not PyDateTime_CheckExact(dt):
        # i.e. is a Timestamp
        return dt.tz_localize(tz)
    elif tz == 'UTC' or tz is UTC:
        return UTC.localize(dt)
    try:
        # datetime.replace with pytz may be incorrect result
        return tz.localize(dt)
    except AttributeError:
        return dt.replace(tzinfo=tz)


# ----------------------------------------------------------------------
# Timezone Conversion

cdef inline int64_t[:] _tz_convert_dst(int64_t[:] values, tzinfo tz,
                                       bint to_utc=True):
    """
    tz_convert for non-UTC non-tzlocal cases where we have to check
    DST transitions pointwise.

    Parameters
    ----------
    values : ndarray[int64_t]
    tz : tzinfo
    to_utc : bool
        True if converting _to_ UTC, False if converting _from_ utc

    Returns
    -------
    result : ndarray[int64_t]
    """
    cdef:
        Py_ssize_t n = len(values)
        Py_ssize_t i, pos
        int64_t[:] result = np.empty(n, dtype=np.int64)
        ndarray[int64_t] trans
        int64_t[:] deltas
        int64_t v

    trans, deltas, typ = get_dst_info(tz)
    if not to_utc:
        # We add `offset` below instead of subtracting it
        deltas = -1 * np.array(deltas, dtype='i8')

    for i in range(n):
        v = values[i]
        if v == NPY_NAT:
            result[i] = v
        else:
            # TODO: Is it more efficient to call searchsorted pointwise or
            # on `values` outside the loop?  We are not consistent about this.
            # relative effiency of pointwise increases with number of iNaTs
            pos = trans.searchsorted(v, side='right') - 1
            if pos < 0:
                raise ValueError('First time before start of DST info')
            result[i] = v - deltas[pos]

    return result


cdef inline int64_t _tz_convert_tzlocal_utc(int64_t val, tzinfo tz,
                                            bint to_utc=True):
    """
    Convert the i8 representation of a datetime from a tzlocal timezone to
    UTC, or vice-versa.

    Private, not intended for use outside of tslibs.conversion

    Parameters
    ----------
    val : int64_t
    tz : tzinfo
    to_utc : bint
        True if converting tzlocal _to_ UTC, False if going the other direction

    Returns
    -------
    result : int64_t
    """
    cdef:
        npy_datetimestruct dts
        int64_t result, delta
        datetime dt

    dt64_to_dtstruct(val, &dts)
    dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                  dts.min, dts.sec, dts.us, tz)
    delta = int(get_utcoffset(tz, dt).total_seconds()) * 1000000000

    if not to_utc:
        return val + delta
    return val - delta


cdef inline int64_t tz_convert_utc_to_tzlocal(int64_t utc_val, tzinfo tz):
    """
    Parameters
    ----------
    utc_val : int64_t
    tz : tzinfo

    Returns
    -------
    local_val : int64_t
    """
    return _tz_convert_tzlocal_utc(utc_val, tz, to_utc=False)


cpdef int64_t tz_convert_single(int64_t val, object tz1, object tz2):
    """
    Convert the val (in i8) from timezone1 to timezone2

    This is a single timezone version of tz_convert

    Parameters
    ----------
    val : int64
    tz1 : string / timezone object
    tz2 : string / timezone object

    Returns
    -------
    converted: int64
    """
    cdef:
        int64_t[:] deltas
        Py_ssize_t pos
        int64_t v, offset, utc_date
        npy_datetimestruct dts
        int64_t arr[1]

    # See GH#17734 We should always be converting either from UTC or to UTC
    assert (is_utc(tz1) or tz1 == 'UTC') or (is_utc(tz2) or tz2 == 'UTC')

    if val == NPY_NAT:
        return val

    # Convert to UTC
    if is_tzlocal(tz1):
        utc_date = _tz_convert_tzlocal_utc(val, tz1, to_utc=True)
    elif get_timezone(tz1) != 'UTC':
        arr[0] = val
        utc_date = _tz_convert_dst(arr, tz1, to_utc=True)[0]
    else:
        utc_date = val

    if get_timezone(tz2) == 'UTC':
        return utc_date
    elif is_tzlocal(tz2):
        return _tz_convert_tzlocal_utc(utc_date, tz2, to_utc=False)
    else:
        # Convert UTC to other timezone
        arr[0] = utc_date
        # Note: at least with cython 0.28.3, doing a lookup `[0]` in the next
        # line is sensitive to the declared return type of _tz_convert_dst;
        # if it is declared as returning ndarray[int64_t], a compile-time error
        # is raised.
        return _tz_convert_dst(arr, tz2, to_utc=False)[0]


cdef inline int64_t[:] _tz_convert_one_way(int64_t[:] vals, object tz,
                                           bint to_utc):
    """
    Convert the given values (in i8) either to UTC or from UTC.

    Parameters
    ----------
    vals : int64 ndarray
    tz1 : string / timezone object
    to_utc : bint

    Returns
    -------
    converted : ndarray[int64_t]
    """
    cdef:
        int64_t[:] converted, result
        Py_ssize_t i, n = len(vals)
        int64_t val

    if get_timezone(tz) != 'UTC':
        converted = np.empty(n, dtype=np.int64)
        if is_tzlocal(tz):
            for i in range(n):
                val = vals[i]
                if val == NPY_NAT:
                    converted[i] = NPY_NAT
                else:
                    converted[i] = _tz_convert_tzlocal_utc(val, tz, to_utc)
        else:
            converted = _tz_convert_dst(vals, tz, to_utc)
    else:
        converted = vals

    return converted


@cython.boundscheck(False)
@cython.wraparound(False)
def tz_convert(int64_t[:] vals, object tz1, object tz2):
    """
    Convert the values (in i8) from timezone1 to timezone2

    Parameters
    ----------
    vals : int64 ndarray
    tz1 : string / timezone object
    tz2 : string / timezone object

    Returns
    -------
    int64 ndarray of converted
    """
    cdef:
        int64_t[:] utc_dates, converted

    if len(vals) == 0:
        return np.array([], dtype=np.int64)

    # Convert to UTC
    utc_dates = _tz_convert_one_way(vals, tz1, to_utc=True)
    converted = _tz_convert_one_way(utc_dates, tz2, to_utc=False)
    return np.array(converted, dtype=np.int64)


# TODO: cdef scalar version to call from convert_str_to_tsobject
@cython.boundscheck(False)
@cython.wraparound(False)
def tz_localize_to_utc(ndarray[int64_t] vals, object tz, object ambiguous=None,
                       object errors='raise'):
    """
    Localize tzinfo-naive i8 to given time zone (using pytz). If
    there are ambiguities in the values, raise AmbiguousTimeError.

    Parameters
    ----------
    vals : ndarray[int64_t]
    tz : tzinfo or None
    ambiguous : str, bool, or arraylike
        If arraylike, must have the same length as vals
    errors : {"raise", "coerce"}, default "raise"

    Returns
    -------
    localized : ndarray[int64_t]
    """
    cdef:
        ndarray[int64_t] trans
        int64_t[:] deltas, idx_shifted
        ndarray ambiguous_array
        Py_ssize_t i, idx, pos, ntrans, n = len(vals)
        int64_t *tdata
        int64_t v, left, right
        ndarray[int64_t] result, result_a, result_b, dst_hours
        npy_datetimestruct dts
        bint infer_dst = False, is_dst = False, fill = False
        bint is_coerce = errors == 'coerce', is_raise = errors == 'raise'

    # Vectorized version of DstTzInfo.localize

    assert is_coerce or is_raise

    if tz == UTC or tz is None:
        return vals

    result = np.empty(n, dtype=np.int64)

    if is_tzlocal(tz):
        for i in range(n):
            v = vals[i]
            result[i] = _tz_convert_tzlocal_utc(v, tz, to_utc=True)
        return result

    if is_string_object(ambiguous):
        if ambiguous == 'infer':
            infer_dst = True
        elif ambiguous == 'NaT':
            fill = True
    elif isinstance(ambiguous, bool):
        is_dst = True
        if ambiguous:
            ambiguous_array = np.ones(len(vals), dtype=bool)
        else:
            ambiguous_array = np.zeros(len(vals), dtype=bool)
    elif hasattr(ambiguous, '__iter__'):
        is_dst = True
        if len(ambiguous) != len(vals):
            raise ValueError("Length of ambiguous bool-array must be "
                             "the same size as vals")
        ambiguous_array = np.asarray(ambiguous)

    trans, deltas, typ = get_dst_info(tz)

    tdata = <int64_t*> cnp.PyArray_DATA(trans)
    ntrans = len(trans)

    result_a = np.empty(n, dtype=np.int64)
    result_b = np.empty(n, dtype=np.int64)
    result_a.fill(NPY_NAT)
    result_b.fill(NPY_NAT)

    # left side
    idx_shifted = (np.maximum(0, trans.searchsorted(
        vals - DAY_NS, side='right') - 1)).astype(np.int64)

    for i in range(n):
        v = vals[i] - deltas[idx_shifted[i]]
        pos = bisect_right_i8(tdata, v, ntrans) - 1

        # timestamp falls to the left side of the DST transition
        if v + deltas[pos] == vals[i]:
            result_a[i] = v

    # right side
    idx_shifted = (np.maximum(0, trans.searchsorted(
        vals + DAY_NS, side='right') - 1)).astype(np.int64)

    for i in range(n):
        v = vals[i] - deltas[idx_shifted[i]]
        pos = bisect_right_i8(tdata, v, ntrans) - 1

        # timestamp falls to the right side of the DST transition
        if v + deltas[pos] == vals[i]:
            result_b[i] = v

    if infer_dst:
        dst_hours = np.empty(n, dtype=np.int64)
        dst_hours.fill(NPY_NAT)

        # Get the ambiguous hours (given the above, these are the hours
        # where result_a != result_b and neither of them are NAT)
        both_nat = np.logical_and(result_a != NPY_NAT, result_b != NPY_NAT)
        both_eq = result_a == result_b
        trans_idx = np.squeeze(np.nonzero(np.logical_and(both_nat, ~both_eq)))
        if trans_idx.size == 1:
            stamp = _render_tstamp(vals[trans_idx])
            raise pytz.AmbiguousTimeError(
                "Cannot infer dst time from %s as there "
                "are no repeated times" % stamp)
        # Split the array into contiguous chunks (where the difference between
        # indices is 1).  These are effectively dst transitions in different
        # years which is useful for checking that there is not an ambiguous
        # transition in an individual year.
        if trans_idx.size > 0:
            one_diff = np.where(np.diff(trans_idx) != 1)[0] +1
            trans_grp = np.array_split(trans_idx, one_diff)

            # Iterate through each day, if there are no hours where the
            # delta is negative (indicates a repeat of hour) the switch
            # cannot be inferred
            for grp in trans_grp:

                delta = np.diff(result_a[grp])
                if grp.size == 1 or np.all(delta > 0):
                    stamp = _render_tstamp(vals[grp[0]])
                    raise pytz.AmbiguousTimeError(stamp)

                # Find the index for the switch and pull from a for dst and b
                # for standard
                switch_idx = (delta <= 0).nonzero()[0]
                if switch_idx.size > 1:
                    raise pytz.AmbiguousTimeError(
                        "There are %i dst switches when "
                        "there should only be 1." % switch_idx.size)
                switch_idx = switch_idx[0] + 1
                # Pull the only index and adjust
                a_idx = grp[:switch_idx]
                b_idx = grp[switch_idx:]
                dst_hours[grp] = np.hstack((result_a[a_idx], result_b[b_idx]))

    for i in range(n):
        left = result_a[i]
        right = result_b[i]
        if vals[i] == NPY_NAT:
            result[i] = vals[i]
        elif left != NPY_NAT and right != NPY_NAT:
            if left == right:
                result[i] = left
            else:
                if infer_dst and dst_hours[i] != NPY_NAT:
                    result[i] = dst_hours[i]
                elif is_dst:
                    if ambiguous_array[i]:
                        result[i] = left
                    else:
                        result[i] = right
                elif fill:
                    result[i] = NPY_NAT
                else:
                    stamp = _render_tstamp(vals[i])
                    raise pytz.AmbiguousTimeError(
                        "Cannot infer dst time from %r, try using the "
                        "'ambiguous' argument" % stamp)
        elif left != NPY_NAT:
            result[i] = left
        elif right != NPY_NAT:
            result[i] = right
        else:
            if is_coerce:
                result[i] = NPY_NAT
            else:
                stamp = _render_tstamp(vals[i])
                raise pytz.NonExistentTimeError(stamp)

    return result


cdef inline bisect_right_i8(int64_t *data, int64_t val, Py_ssize_t n):
    cdef Py_ssize_t pivot, left = 0, right = n

    assert n >= 1

    # edge cases
    if val > data[n - 1]:
        return n

    if val < data[0]:
        return 0

    while left < right:
        pivot = left + (right - left) // 2

        if data[pivot] <= val:
            left = pivot + 1
        else:
            right = pivot

    return left


cdef inline str _render_tstamp(int64_t val):
    """ Helper function to render exception messages"""
    from timestamps import Timestamp
    return str(Timestamp(val))


# ----------------------------------------------------------------------
# Normalization


def normalize_date(object dt):
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
        if not PyDateTime_CheckExact(dt):
            # i.e. a Timestamp object
            return dt.replace(hour=0, minute=0, second=0, microsecond=0,
                              nanosecond=0)
        else:
            # regular datetime object
            return dt.replace(hour=0, minute=0, second=0, microsecond=0)
            # TODO: Make sure DST crossing is handled correctly here
    elif PyDate_Check(dt):
        return datetime(dt.year, dt.month, dt.day)
    else:
        raise TypeError('Unrecognized type: %s' % type(dt))


@cython.wraparound(False)
@cython.boundscheck(False)
def normalize_i8_timestamps(int64_t[:] stamps, tz=None):
    """
    Normalize each of the (nanosecond) timestamps in the given array by
    rounding down to the beginning of the day (i.e. midnight).  If `tz`
    is not None, then this is midnight for this timezone.

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
        npy_datetimestruct dts
        int64_t[:] result = np.empty(n, dtype=np.int64)

    if tz is not None:
        tz = maybe_get_tz(tz)
        result = _normalize_local(stamps, tz)
    else:
        with nogil:
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                dt64_to_dtstruct(stamps[i], &dts)
                result[i] = _normalized_stamp(&dts)

    return result.base  # .base to access underlying np.ndarray


@cython.wraparound(False)
@cython.boundscheck(False)
cdef int64_t[:] _normalize_local(int64_t[:] stamps, object tz):
    """
    Normalize each of the (nanosecond) timestamps in the given array by
    rounding down to the beginning of the day (i.e. midnight) for the
    given timezone `tz`.

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
        ndarray[int64_t] trans
        int64_t[:] deltas
        Py_ssize_t[:] pos
        npy_datetimestruct dts
        int64_t delta

    if is_utc(tz):
        with nogil:
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                dt64_to_dtstruct(stamps[i], &dts)
                result[i] = _normalized_stamp(&dts)
    elif is_tzlocal(tz):
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


def is_date_array_normalized(int64_t[:] stamps, tz=None):
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
        Py_ssize_t pos, i, n = len(stamps)
        ndarray[int64_t] trans
        int64_t[:] deltas
        npy_datetimestruct dts
        int64_t local_val, delta

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
            for i in range(n):
                # Adjust datetime64 timestamp, recompute datetimestruct
                pos = trans.searchsorted(stamps[i]) - 1

                dt64_to_dtstruct(stamps[i] + deltas[pos], &dts)
                if (dts.hour + dts.min + dts.sec + dts.us) > 0:
                    return False

    return True
