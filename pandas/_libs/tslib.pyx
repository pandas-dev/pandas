# -*- coding: utf-8 -*-
import cython

from cpython.datetime cimport (PyDateTime_Check, PyDate_Check,
                               PyDateTime_CheckExact,
                               PyDateTime_IMPORT,
                               timedelta, datetime, date, time)
# import datetime C API
PyDateTime_IMPORT


cimport numpy as cnp
from numpy cimport int64_t, ndarray, float64_t
import numpy as np
cnp.import_array()

import pytz

from pandas._libs.util cimport (
    is_integer_object, is_float_object, is_string_object, is_datetime64_object)


from pandas._libs.tslibs.np_datetime cimport (
    check_dts_bounds, npy_datetimestruct, _string_to_dts, dt64_to_dtstruct,
    dtstruct_to_dt64, pydatetime_to_dt64, pydate_to_dt64, get_datetime64_value)
from pandas._libs.tslibs.np_datetime import OutOfBoundsDatetime

from pandas._libs.tslibs.parsing import parse_datetime_string

from pandas._libs.tslibs.timedeltas cimport cast_from_unit
from pandas._libs.tslibs.timezones cimport is_utc, is_tzlocal, get_dst_info
from pandas._libs.tslibs.timezones import UTC
from pandas._libs.tslibs.conversion cimport (
    tz_convert_single, _TSObject, convert_datetime_to_tsobject,
    get_datetime64_nanos, tz_convert_utc_to_tzlocal)

# many modules still look for NaT and iNaT here despite them not being needed
from pandas._libs.tslibs.nattype import nat_strings, iNaT  # noqa:F821
from pandas._libs.tslibs.nattype cimport (
    checknull_with_nat, NPY_NAT, c_NaT as NaT)

from pandas._libs.tslibs.offsets cimport to_offset

from pandas._libs.tslibs.timestamps cimport create_timestamp_from_ts
from pandas._libs.tslibs.timestamps import Timestamp


cdef bint PY2 = str == bytes


cdef inline object create_datetime_from_ts(
        int64_t value, npy_datetimestruct dts,
        object tz, object freq):
    """ convenience routine to construct a datetime.datetime from its parts """
    return datetime(dts.year, dts.month, dts.day, dts.hour,
                    dts.min, dts.sec, dts.us, tz)


cdef inline object create_date_from_ts(
        int64_t value, npy_datetimestruct dts,
        object tz, object freq):
    """ convenience routine to construct a datetime.date from its parts """
    return date(dts.year, dts.month, dts.day)


cdef inline object create_time_from_ts(
        int64_t value, npy_datetimestruct dts,
        object tz, object freq):
    """ convenience routine to construct a datetime.time from its parts """
    return time(dts.hour, dts.min, dts.sec, dts.us, tz)


@cython.wraparound(False)
@cython.boundscheck(False)
def ints_to_pydatetime(int64_t[:] arr, object tz=None, object freq=None,
                       str box="datetime"):
    """
    Convert an i8 repr to an ndarray of datetimes, date, time or Timestamp

    Parameters
    ----------
    arr  : array of i8
    tz   : str, default None
         convert to this timezone
    freq : str/Offset, default None
         freq to convert
    box  : {'datetime', 'timestamp', 'date', 'time'}, default 'datetime'
         If datetime, convert to datetime.datetime
         If date, convert to datetime.date
         If time, convert to datetime.time
         If Timestamp, convert to pandas.Timestamp

    Returns
    -------
    result : array of dtype specified by box
    """

    cdef:
        Py_ssize_t i, n = len(arr)
        ndarray[int64_t] trans
        int64_t[:] deltas
        Py_ssize_t pos
        npy_datetimestruct dts
        object dt, new_tz
        str typ
        int64_t value, delta, local_value
        ndarray[object] result = np.empty(n, dtype=object)
        object (*func_create)(int64_t, npy_datetimestruct, object, object)

    if box == "date":
        assert (tz is None), "tz should be None when converting to date"

        func_create = create_date_from_ts
    elif box == "timestamp":
        func_create = create_timestamp_from_ts

        if is_string_object(freq):
            freq = to_offset(freq)
    elif box == "time":
        func_create = create_time_from_ts
    elif box == "datetime":
        func_create = create_datetime_from_ts
    else:
        raise ValueError("box must be one of 'datetime', 'date', 'time' or"
                         " 'timestamp'")

    if is_utc(tz) or tz is None:
        for i in range(n):
            value = arr[i]
            if value == NPY_NAT:
                result[i] = NaT
            else:
                dt64_to_dtstruct(value, &dts)
                result[i] = func_create(value, dts, tz, freq)
    elif is_tzlocal(tz):
        for i in range(n):
            value = arr[i]
            if value == NPY_NAT:
                result[i] = NaT
            else:
                # Python datetime objects do not support nanosecond
                # resolution (yet, PEP 564). Need to compute new value
                # using the i8 representation.
                local_value = tz_convert_utc_to_tzlocal(value, tz)
                dt64_to_dtstruct(local_value, &dts)
                result[i] = func_create(value, dts, tz, freq)
    else:
        trans, deltas, typ = get_dst_info(tz)

        if typ not in ['pytz', 'dateutil']:
            # static/fixed; in this case we know that len(delta) == 1
            delta = deltas[0]
            for i in range(n):
                value = arr[i]
                if value == NPY_NAT:
                    result[i] = NaT
                else:
                    # Adjust datetime64 timestamp, recompute datetimestruct
                    dt64_to_dtstruct(value + delta, &dts)
                    result[i] = func_create(value, dts, tz, freq)

        elif typ == 'dateutil':
            # no zone-name change for dateutil tzs - dst etc
            # represented in single object.
            for i in range(n):
                value = arr[i]
                if value == NPY_NAT:
                    result[i] = NaT
                else:
                    # Adjust datetime64 timestamp, recompute datetimestruct
                    pos = trans.searchsorted(value, side='right') - 1
                    dt64_to_dtstruct(value + deltas[pos], &dts)
                    result[i] = func_create(value, dts, tz, freq)
        else:
            # pytz
            for i in range(n):
                value = arr[i]
                if value == NPY_NAT:
                    result[i] = NaT
                else:
                    # Adjust datetime64 timestamp, recompute datetimestruct
                    pos = trans.searchsorted(value, side='right') - 1
                    # find right representation of dst etc in pytz timezone
                    new_tz = tz._tzinfos[tz._transition_info[pos]]

                    dt64_to_dtstruct(value + deltas[pos], &dts)
                    result[i] = func_create(value, dts, new_tz, freq)

    return result


def _test_parse_iso8601(object ts):
    """
    TESTING ONLY: Parse string into Timestamp using iso8601 parser. Used
    only for testing, actual construction uses `convert_str_to_tsobject`
    """
    cdef:
        _TSObject obj
        int out_local = 0, out_tzoffset = 0

    obj = _TSObject()

    if ts == 'now':
        return Timestamp.utcnow()
    elif ts == 'today':
        return Timestamp.now().normalize()

    _string_to_dts(ts, &obj.dts, &out_local, &out_tzoffset)
    obj.value = dtstruct_to_dt64(&obj.dts)
    check_dts_bounds(&obj.dts)
    if out_local == 1:
        obj.tzinfo = pytz.FixedOffset(out_tzoffset)
        obj.value = tz_convert_single(obj.value, obj.tzinfo, UTC)
        return Timestamp(obj.value, tz=obj.tzinfo)
    else:
        return Timestamp(obj.value)


@cython.wraparound(False)
@cython.boundscheck(False)
def format_array_from_datetime(ndarray[int64_t] values, object tz=None,
                               object format=None, object na_rep=None):
    """
    return a np object array of the string formatted values

    Parameters
    ----------
    values : a 1-d i8 array
    tz : the timezone (or None)
    format : optional, default is None
          a strftime capable string
    na_rep : optional, default is None
          a nat format

    """
    cdef:
        int64_t val, ns, N = len(values)
        ndarray[int64_t] consider_values
        bint show_ms = 0, show_us = 0, show_ns = 0, basic_format = 0
        ndarray[object] result = np.empty(N, dtype=object)
        object ts, res
        npy_datetimestruct dts

    if na_rep is None:
        na_rep = 'NaT'

    # if we don't have a format nor tz, then choose
    # a format based on precision
    basic_format = format is None and tz is None
    if basic_format:
        consider_values = values[values != NPY_NAT]
        show_ns = (consider_values % 1000).any()

        if not show_ns:
            consider_values //= 1000
            show_us = (consider_values % 1000).any()

            if not show_ms:
                consider_values //= 1000
                show_ms = (consider_values % 1000).any()

    for i in range(N):
        val = values[i]

        if val == NPY_NAT:
            result[i] = na_rep
        elif basic_format:

            dt64_to_dtstruct(val, &dts)
            res = '%d-%.2d-%.2d %.2d:%.2d:%.2d' % (dts.year,
                                                   dts.month,
                                                   dts.day,
                                                   dts.hour,
                                                   dts.min,
                                                   dts.sec)

            if show_ns:
                ns = dts.ps / 1000
                res += '.%.9d' % (ns + 1000 * dts.us)
            elif show_us:
                res += '.%.6d' % dts.us
            elif show_ms:
                res += '.%.3d' % (dts.us /1000)

            result[i] = res

        else:

            ts = Timestamp(val, tz=tz)
            if format is None:
                result[i] = str(ts)
            else:

                # invalid format string
                # requires dates > 1900
                try:
                    result[i] = ts.strftime(format)
                except ValueError:
                    result[i] = str(ts)

    return result


def array_with_unit_to_datetime(ndarray values, object unit,
                                str errors='coerce'):
    """
    convert the ndarray according to the unit
    if errors:
      - raise: return converted values or raise OutOfBoundsDatetime
          if out of range on the conversion or
          ValueError for other conversions (e.g. a string)
      - ignore: return non-convertible values as the same unit
      - coerce: NaT for non-convertibles

    """
    cdef:
        Py_ssize_t i, j, n=len(values)
        int64_t m
        ndarray[float64_t] fvalues
        ndarray mask
        bint is_ignore = errors=='ignore'
        bint is_coerce = errors=='coerce'
        bint is_raise = errors=='raise'
        bint need_to_iterate = True
        ndarray[int64_t] iresult
        ndarray[object] oresult

    assert is_ignore or is_coerce or is_raise

    if unit == 'ns':
        if issubclass(values.dtype.type, np.integer):
            return values.astype('M8[ns]')
        return array_to_datetime(values.astype(object), errors=errors)[0]

    m = cast_from_unit(None, unit)

    if is_raise:

        # try a quick conversion to i8
        # if we have nulls that are not type-compat
        # then need to iterate
        try:
            iresult = values.astype('i8', casting='same_kind', copy=False)
            mask = iresult == NPY_NAT
            iresult[mask] = 0
            fvalues = iresult.astype('f8') * m
            need_to_iterate = False
        except:
            pass

        # check the bounds
        if not need_to_iterate:

            if ((fvalues < Timestamp.min.value).any()
                    or (fvalues > Timestamp.max.value).any()):
                raise OutOfBoundsDatetime("cannot convert input with unit "
                                          "'{unit}'".format(unit=unit))
            result = (iresult * m).astype('M8[ns]')
            iresult = result.view('i8')
            iresult[mask] = NPY_NAT
            return result

    result = np.empty(n, dtype='M8[ns]')
    iresult = result.view('i8')

    try:
        for i in range(n):
            val = values[i]

            if checknull_with_nat(val):
                iresult[i] = NPY_NAT

            elif is_integer_object(val) or is_float_object(val):

                if val != val or val == NPY_NAT:
                    iresult[i] = NPY_NAT
                else:
                    try:
                        iresult[i] = cast_from_unit(val, unit)
                    except OverflowError:
                        if is_raise:
                            raise OutOfBoundsDatetime(
                                "cannot convert input {val} with the unit "
                                "'{unit}'".format(val=val, unit=unit))
                        elif is_ignore:
                            raise AssertionError
                        iresult[i] = NPY_NAT

            elif is_string_object(val):
                if len(val) == 0 or val in nat_strings:
                    iresult[i] = NPY_NAT

                else:
                    try:
                        iresult[i] = cast_from_unit(float(val), unit)
                    except ValueError:
                        if is_raise:
                            raise ValueError(
                                "non convertible value {val} with the unit "
                                "'{unit}'".format(val=val, unit=unit))
                        elif is_ignore:
                            raise AssertionError
                        iresult[i] = NPY_NAT
                    except:
                        if is_raise:
                            raise OutOfBoundsDatetime(
                                "cannot convert input {val} with the unit "
                                "'{unit}'".format(val=val, unit=unit))
                        elif is_ignore:
                            raise AssertionError
                        iresult[i] = NPY_NAT

            else:

                if is_raise:
                    raise ValueError("unit='{0}' not valid with non-numerical "
                                     "val='{1}'".format(unit, val))
                if is_ignore:
                    raise AssertionError

                iresult[i] = NPY_NAT

        return result

    except AssertionError:
        pass

    # we have hit an exception
    # and are in ignore mode
    # redo as object

    oresult = np.empty(n, dtype=object)
    for i in range(n):
        val = values[i]

        if checknull_with_nat(val):
            oresult[i] = NaT
        elif is_integer_object(val) or is_float_object(val):

            if val != val or val == NPY_NAT:
                oresult[i] = NaT
            else:
                try:
                    oresult[i] = Timestamp(cast_from_unit(val, unit))
                except:
                    oresult[i] = val

        elif is_string_object(val):
            if len(val) == 0 or val in nat_strings:
                oresult[i] = NaT

            else:
                oresult[i] = val

    return oresult


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef array_to_datetime(ndarray[object] values, str errors='raise',
                        bint dayfirst=False, bint yearfirst=False,
                        object utc=None, bint require_iso8601=False):
    """
    Converts a 1D array of date-like values to a numpy array of either:
        1) datetime64[ns] data
        2) datetime.datetime objects, if OutOfBoundsDatetime or TypeError
           is encountered

    Also returns a pytz.FixedOffset if an array of strings with the same
    timezone offset is passed and utc=True is not passed. Otherwise, None
    is returned

    Handles datetime.date, datetime.datetime, np.datetime64 objects, numeric,
    strings

    Parameters
    ----------
    values : ndarray of object
         date-like objects to convert
    errors : str, default 'raise'
         error behavior when parsing
    dayfirst : bool, default False
         dayfirst parsing behavior when encountering datetime strings
    yearfirst : bool, default False
         yearfirst parsing behavior when encountering datetime strings
    utc : bool, default None
         indicator whether the dates should be UTC
    require_iso8601 : bool, default False
         indicator whether the datetime string should be iso8601

    Returns
    -------
    tuple (ndarray, tzoffset)
    """
    cdef:
        Py_ssize_t i, n = len(values)
        object val, py_dt, tz, tz_out = None
        ndarray[int64_t] iresult
        ndarray[object] oresult
        npy_datetimestruct dts
        bint utc_convert = bool(utc)
        bint seen_integer = 0
        bint seen_string = 0
        bint seen_datetime = 0
        bint seen_datetime_offset = 0
        bint is_raise = errors=='raise'
        bint is_ignore = errors=='ignore'
        bint is_coerce = errors=='coerce'
        bint is_same_offsets
        _TSObject _ts
        int64_t value
        int out_local=0, out_tzoffset=0
        float offset_seconds, tz_offset
        set out_tzoffset_vals = set()

    # specify error conditions
    assert is_raise or is_ignore or is_coerce

    result = np.empty(n, dtype='M8[ns]')
    iresult = result.view('i8')

    try:
        for i in range(n):
            val = values[i]

            try:
                if checknull_with_nat(val):
                    iresult[i] = NPY_NAT

                elif PyDateTime_Check(val):
                    seen_datetime = 1
                    if val.tzinfo is not None:
                        if utc_convert:
                            _ts = convert_datetime_to_tsobject(val, None)
                            iresult[i] = _ts.value
                        else:
                            raise ValueError('Tz-aware datetime.datetime '
                                             'cannot be converted to '
                                             'datetime64 unless utc=True')
                    else:
                        iresult[i] = pydatetime_to_dt64(val, &dts)
                        if not PyDateTime_CheckExact(val):
                            # i.e. a Timestamp object
                            iresult[i] += val.nanosecond
                        check_dts_bounds(&dts)

                elif PyDate_Check(val):
                    seen_datetime = 1
                    iresult[i] = pydate_to_dt64(val, &dts)
                    check_dts_bounds(&dts)

                elif is_datetime64_object(val):
                    seen_datetime = 1
                    iresult[i] = get_datetime64_nanos(val)

                elif is_integer_object(val) or is_float_object(val):
                    # these must be ns unit by-definition
                    seen_integer = 1

                    if val != val or val == NPY_NAT:
                        iresult[i] = NPY_NAT
                    elif is_raise or is_ignore:
                        iresult[i] = val
                    else:
                        # coerce
                        # we now need to parse this as if unit='ns'
                        # we can ONLY accept integers at this point
                        # if we have previously (or in future accept
                        # datetimes/strings, then we must coerce)
                        try:
                            iresult[i] = cast_from_unit(val, 'ns')
                        except:
                            iresult[i] = NPY_NAT

                elif is_string_object(val):
                    # string
                    seen_string = 1

                    if len(val) == 0 or val in nat_strings:
                        iresult[i] = NPY_NAT
                        continue
                    if isinstance(val, unicode) and PY2:
                        val = val.encode('utf-8')

                    try:
                        _string_to_dts(val, &dts, &out_local, &out_tzoffset)
                    except ValueError:
                        # A ValueError at this point is a _parsing_ error
                        # specifically _not_ OutOfBoundsDatetime
                        if _parse_today_now(val, &iresult[i]):
                            continue
                        elif require_iso8601:
                            # if requiring iso8601 strings, skip trying
                            # other formats
                            if is_coerce:
                                iresult[i] = NPY_NAT
                                continue
                            elif is_raise:
                                raise ValueError("time data {val} doesn't "
                                                 "match format specified"
                                                 .format(val=val))
                            return values, tz_out

                        try:
                            py_dt = parse_datetime_string(val,
                                                          dayfirst=dayfirst,
                                                          yearfirst=yearfirst)
                        except Exception:
                            if is_coerce:
                                iresult[i] = NPY_NAT
                                continue
                            raise TypeError("invalid string coercion to "
                                            "datetime")

                        # If the dateutil parser returned tzinfo, capture it
                        # to check if all arguments have the same tzinfo
                        tz = py_dt.utcoffset()
                        if tz is not None:
                            seen_datetime_offset = 1
                            # dateutil timezone objects cannot be hashed, so
                            # store the UTC offsets in seconds instead
                            out_tzoffset_vals.add(tz.total_seconds())
                        else:
                            # Add a marker for naive string, to track if we are
                            # parsing mixed naive and aware strings
                            out_tzoffset_vals.add('naive')

                        _ts = convert_datetime_to_tsobject(py_dt, None)
                        iresult[i] = _ts.value
                    except:
                        # TODO: What exception are we concerned with here?
                        if is_coerce:
                            iresult[i] = NPY_NAT
                            continue
                        raise
                    else:
                        # No error raised by string_to_dts, pick back up
                        # where we left off
                        value = dtstruct_to_dt64(&dts)
                        if out_local == 1:
                            seen_datetime_offset = 1
                            # Store the out_tzoffset in seconds
                            # since we store the total_seconds of
                            # dateutil.tz.tzoffset objects
                            out_tzoffset_vals.add(out_tzoffset * 60.)
                            tz = pytz.FixedOffset(out_tzoffset)
                            value = tz_convert_single(value, tz, UTC)
                        else:
                            # Add a marker for naive string, to track if we are
                            # parsing mixed naive and aware strings
                            out_tzoffset_vals.add('naive')
                        iresult[i] = value
                        check_dts_bounds(&dts)

                else:
                    if is_coerce:
                        iresult[i] = NPY_NAT
                    else:
                        raise TypeError("{typ} is not convertible to datetime"
                                        .format(typ=type(val)))

            except OutOfBoundsDatetime:
                if is_coerce:
                    iresult[i] = NPY_NAT
                    continue
                elif require_iso8601 and is_string_object(val):
                    # GH#19382 for just-barely-OutOfBounds falling back to
                    # dateutil parser will return incorrect result because
                    # it will ignore nanoseconds
                    if is_raise:
                        raise ValueError("time data {val} doesn't "
                                         "match format specified"
                                         .format(val=val))
                    assert is_ignore
                    return values, tz_out
                raise

    except OutOfBoundsDatetime:
        if is_raise:
            raise

        return ignore_errors_out_of_bounds_fallback(values), tz_out

    except TypeError:
        return array_to_datetime_object(values, is_raise, dayfirst, yearfirst)

    if seen_datetime and seen_integer:
        # we have mixed datetimes & integers

        if is_coerce:
            # coerce all of the integers/floats to NaT, preserve
            # the datetimes and other convertibles
            for i in range(n):
                val = values[i]
                if is_integer_object(val) or is_float_object(val):
                    result[i] = NPY_NAT
        elif is_raise:
            raise ValueError("mixed datetimes and integers in passed array")
        else:
            return array_to_datetime_object(values, is_raise,
                                            dayfirst, yearfirst)

    if seen_datetime_offset and not utc_convert:
        # GH#17697
        # 1) If all the offsets are equal, return one offset for
        #    the parsed dates to (maybe) pass to DatetimeIndex
        # 2) If the offsets are different, then force the parsing down the
        #    object path where an array of datetimes
        #    (with individual dateutil.tzoffsets) are returned
        is_same_offsets = len(out_tzoffset_vals) == 1
        if not is_same_offsets:
            return array_to_datetime_object(values, is_raise,
                                            dayfirst, yearfirst)
        else:
            tz_offset = out_tzoffset_vals.pop()
            tz_out = pytz.FixedOffset(tz_offset / 60.)
    return result, tz_out


cdef inline ignore_errors_out_of_bounds_fallback(ndarray[object] values):
    """
    Fallback for array_to_datetime if an OutOfBoundsDatetime is raised
    and errors == "ignore"

    Parameters
    ----------
    values : ndarray[object]

    Returns
    -------
    ndarray[object]
    """
    cdef:
        Py_ssize_t i, n = len(values)
        object val

    oresult = np.empty(n, dtype=object)

    for i in range(n):
        val = values[i]

        # set as nan except if its a NaT
        if checknull_with_nat(val):
            if isinstance(val, float):
                oresult[i] = np.nan
            else:
                oresult[i] = NaT
        elif is_datetime64_object(val):
            if get_datetime64_value(val) == NPY_NAT:
                oresult[i] = NaT
            else:
                oresult[i] = val.item()
        else:
            oresult[i] = val
    return oresult


@cython.wraparound(False)
@cython.boundscheck(False)
cdef array_to_datetime_object(ndarray[object] values, bint is_raise,
                              bint dayfirst=False, bint yearfirst=False):
    """
    Fall back function for array_to_datetime

    Attempts to parse datetime strings with dateutil to return an array
    of datetime objects

    Parameters
    ----------
    values : ndarray of object
         date-like objects to convert
    is_raise : bool
         error behavior when parsing
    dayfirst : bool, default False
         dayfirst parsing behavior when encountering datetime strings
    yearfirst : bool, default False
         yearfirst parsing behavior when encountering datetime strings

    Returns
    -------
    tuple (ndarray, None)
    """
    cdef:
        Py_ssize_t i, n = len(values)
        object val,
        ndarray[object] oresult
        npy_datetimestruct dts

    oresult = np.empty(n, dtype=object)

    # We return an object array and only attempt to parse:
    # 1) NaT or NaT-like values
    # 2) datetime strings, which we return as datetime.datetime
    for i in range(n):
        val = values[i]
        if checknull_with_nat(val):
            oresult[i] = val
        elif is_string_object(val):
            if len(val) == 0 or val in nat_strings:
                oresult[i] = 'NaT'
                continue
            try:
                oresult[i] = parse_datetime_string(val, dayfirst=dayfirst,
                                                   yearfirst=yearfirst)
                pydatetime_to_dt64(oresult[i], &dts)
                check_dts_bounds(&dts)
            except (ValueError, OverflowError):
                if is_raise:
                    raise
                return values, None
        else:
            if is_raise:
                raise
            return values, None
    return oresult, None


cdef inline bint _parse_today_now(str val, int64_t* iresult):
    # We delay this check for as long as possible
    # because it catches relatively rare cases
    if val == 'now':
        # Note: this is *not* the same as Timestamp('now')
        iresult[0] = Timestamp.utcnow().value
        return True
    elif val == 'today':
        iresult[0] = Timestamp.today().value
        return True
    return False
