# -*- coding: utf-8 -*-
# cython: profile=False
# cython: linetrace=False
# distutils: define_macros=CYTHON_TRACE=0
# distutils: define_macros=CYTHON_TRACE_NOGIL=0

cimport numpy as np
from numpy cimport int64_t, import_array, ndarray, float64_t
import numpy as np


from cpython cimport PyTypeObject, PyFloat_Check

cdef extern from "Python.h":
    cdef PyTypeObject *Py_TYPE(object)

from util cimport (is_integer_object, is_float_object, is_string_object,
                   is_datetime64_object)

from cpython.datetime cimport (PyDateTime_Check, PyDate_Check,
                               PyDateTime_IMPORT,
                               timedelta, datetime, date)
# import datetime C API
PyDateTime_IMPORT
# this is our datetime.pxd
from datetime cimport _string_to_dts


from tslibs.np_datetime cimport (check_dts_bounds,
                                 pandas_datetimestruct,
                                 dt64_to_dtstruct, dtstruct_to_dt64,
                                 pydatetime_to_dt64, pydate_to_dt64,
                                 get_datetime64_value,
                                 days_per_month_table,
                                 dayofweek, is_leapyear)
from tslibs.np_datetime import OutOfBoundsDatetime

from .tslibs.parsing import parse_datetime_string

cimport cython

import warnings

import pytz
UTC = pytz.utc

# initialize numpy
import_array()


from tslibs.timedeltas cimport cast_from_unit
from tslibs.timedeltas import Timedelta
from tslibs.timezones cimport (
    is_utc, is_tzlocal, is_fixed_offset,
    treat_tz_as_pytz,
    get_dst_info)
from tslibs.conversion cimport (tz_convert_single, _TSObject,
                                convert_datetime_to_tsobject,
                                get_datetime64_nanos)
from tslibs.conversion import tz_convert_single

from tslibs.nattype import NaT, nat_strings, iNaT
from tslibs.nattype cimport checknull_with_nat, NPY_NAT

from tslibs.timestamps cimport (create_timestamp_from_ts,
                                _NS_UPPER_BOUND, _NS_LOWER_BOUND)
from tslibs.timestamps import Timestamp


cdef inline object create_datetime_from_ts(
        int64_t value, pandas_datetimestruct dts,
        object tz, object freq):
    """ convenience routine to construct a datetime.datetime from its parts """
    return datetime(dts.year, dts.month, dts.day, dts.hour,
                    dts.min, dts.sec, dts.us, tz)

cdef inline object create_date_from_ts(
        int64_t value, pandas_datetimestruct dts,
        object tz, object freq):
    """ convenience routine to construct a datetime.date from its parts """
    return date(dts.year, dts.month, dts.day)


def ints_to_pydatetime(ndarray[int64_t] arr, tz=None, freq=None,
                       box="datetime"):
    """
    Convert an i8 repr to an ndarray of datetimes, date or Timestamp

    Parameters
    ----------
    arr  : array of i8
    tz   : str, default None
         convert to this timezone
    freq : str/Offset, default None
         freq to convert
    box  : {'datetime', 'timestamp', 'date'}, default 'datetime'
         If datetime, convert to datetime.datetime
         If date, convert to datetime.date
         If Timestamp, convert to pandas.Timestamp

    Returns
    -------
    result : array of dtype specified by box
    """

    assert ((box == "datetime") or (box == "date") or (box == "timestamp")), \
        "box must be one of 'datetime', 'date' or 'timestamp'"

    cdef:
        Py_ssize_t i, n = len(arr)
        ndarray[int64_t] trans, deltas
        pandas_datetimestruct dts
        object dt
        int64_t value
        ndarray[object] result = np.empty(n, dtype=object)
        object (*func_create)(int64_t, pandas_datetimestruct, object, object)

    if box == "date":
        assert (tz is None), "tz should be None when converting to date"

        func_create = create_date_from_ts
    elif box == "timestamp":
        func_create = create_timestamp_from_ts

        if is_string_object(freq):
            from pandas.tseries.frequencies import to_offset
            freq = to_offset(freq)
    elif box == "datetime":
        func_create = create_datetime_from_ts

    if tz is not None:
        if is_utc(tz):
            for i in range(n):
                value = arr[i]
                if value == NPY_NAT:
                    result[i] = NaT
                else:
                    dt64_to_dtstruct(value, &dts)
                    result[i] = func_create(value, dts, tz, freq)
        elif is_tzlocal(tz) or is_fixed_offset(tz):
            for i in range(n):
                value = arr[i]
                if value == NPY_NAT:
                    result[i] = NaT
                else:
                    dt64_to_dtstruct(value, &dts)
                    dt = create_datetime_from_ts(value, dts, tz, freq)
                    dt = dt + tz.utcoffset(dt)
                    if box:
                        dt = Timestamp(dt)
                    result[i] = dt
        else:
            trans, deltas, typ = get_dst_info(tz)

            for i in range(n):

                value = arr[i]
                if value == NPY_NAT:
                    result[i] = NaT
                else:

                    # Adjust datetime64 timestamp, recompute datetimestruct
                    pos = trans.searchsorted(value, side='right') - 1
                    if treat_tz_as_pytz(tz):
                        # find right representation of dst etc in pytz timezone
                        new_tz = tz._tzinfos[tz._transition_info[pos]]
                    else:
                        # no zone-name change for dateutil tzs - dst etc
                        # represented in single object.
                        new_tz = tz

                    dt64_to_dtstruct(value + deltas[pos], &dts)
                    result[i] = func_create(value, dts, new_tz, freq)
    else:
        for i in range(n):

            value = arr[i]
            if value == NPY_NAT:
                result[i] = NaT
            else:
                dt64_to_dtstruct(value, &dts)
                result[i] = func_create(value, dts, None, freq)

    return result


def ints_to_pytimedelta(ndarray[int64_t] arr, box=False):
    # convert an i8 repr to an ndarray of timedelta or Timedelta (if box ==
    # True)

    cdef:
        Py_ssize_t i, n = len(arr)
        int64_t value
        ndarray[object] result = np.empty(n, dtype=object)

    for i in range(n):

        value = arr[i]
        if value == NPY_NAT:
            result[i] = NaT
        else:
            if box:
                result[i] = Timedelta(value)
            else:
                result[i] = timedelta(microseconds=int(value) / 1000)

    return result


cdef PyTypeObject* ts_type = <PyTypeObject*> Timestamp


cdef inline bint is_timestamp(object o):
    return Py_TYPE(o) == ts_type  # isinstance(o, Timestamp)


def _test_parse_iso8601(object ts):
    """
    TESTING ONLY: Parse string into Timestamp using iso8601 parser. Used
    only for testing, actual construction uses `convert_str_to_tsobject`
    """
    cdef:
        _TSObject obj
        int out_local = 0, out_tzoffset = 0

    obj = _TSObject()

    _string_to_dts(ts, &obj.dts, &out_local, &out_tzoffset)
    obj.value = dtstruct_to_dt64(&obj.dts)
    check_dts_bounds(&obj.dts)
    if out_local == 1:
        obj.tzinfo = pytz.FixedOffset(out_tzoffset)
        obj.value = tz_convert_single(obj.value, obj.tzinfo, 'UTC')
        return Timestamp(obj.value, tz=obj.tzinfo)
    else:
        return Timestamp(obj.value)


cpdef inline object _localize_pydatetime(object dt, object tz):
    """
    Take a datetime/Timestamp in UTC and localizes to timezone tz.
    """
    if tz is None:
        return dt
    elif isinstance(dt, Timestamp):
        return dt.tz_localize(tz)
    elif tz == 'UTC' or tz is UTC:
        return UTC.localize(dt)
    try:
        # datetime.replace with pytz may be incorrect result
        return tz.localize(dt)
    except AttributeError:
        return dt.replace(tzinfo=tz)


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
        pandas_datetimestruct dts

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


# const for parsers

_MONTHS = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL',
           'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
_MONTH_NUMBERS = {k: i for i, k in enumerate(_MONTHS)}
_MONTH_ALIASES = {(k + 1): v for k, v in enumerate(_MONTHS)}


cpdef array_with_unit_to_datetime(ndarray values, unit, errors='coerce'):
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
        bint need_to_iterate=True
        ndarray[int64_t] iresult
        ndarray[object] oresult

    assert is_ignore or is_coerce or is_raise

    if unit == 'ns':
        if issubclass(values.dtype.type, np.integer):
            return values.astype('M8[ns]')
        return array_to_datetime(values.astype(object), errors=errors)

    m = cast_from_unit(None, unit)

    if is_raise:

        # try a quick conversion to i8
        # if we have nulls that are not type-compat
        # then need to iterate
        try:
            iresult = values.astype('i8', casting='same_kind', copy=False)
            mask = iresult == iNaT
            iresult[mask] = 0
            fvalues = iresult.astype('f8') * m
            need_to_iterate=False
        except:
            pass

        # check the bounds
        if not need_to_iterate:

            if ((fvalues < _NS_LOWER_BOUND).any()
                    or (fvalues > _NS_UPPER_BOUND).any()):
                raise OutOfBoundsDatetime(
                    "cannot convert input with unit '{0}'".format(unit))
            result = (iresult *m).astype('M8[ns]')
            iresult = result.view('i8')
            iresult[mask] = iNaT
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
                                "cannot convert input {0} with the unit "
                                "'{1}'".format(val, unit))
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
                                "non convertible value {0} with the unit "
                                "'{1}'".format(val, unit))
                        elif is_ignore:
                            raise AssertionError
                        iresult[i] = NPY_NAT
                    except:
                        if is_raise:
                            raise OutOfBoundsDatetime(
                                "cannot convert input {0} with the unit "
                                "'{1}'".format(val, unit))
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


cpdef array_to_datetime(ndarray[object] values, errors='raise',
                        dayfirst=False, yearfirst=False,
                        format=None, utc=None,
                        require_iso8601=False):
    cdef:
        Py_ssize_t i, n = len(values)
        object val, py_dt
        ndarray[int64_t] iresult
        ndarray[object] oresult
        pandas_datetimestruct dts
        bint utc_convert = bool(utc)
        bint seen_integer = 0
        bint seen_string = 0
        bint seen_datetime = 0
        bint is_raise = errors=='raise'
        bint is_ignore = errors=='ignore'
        bint is_coerce = errors=='coerce'
        _TSObject _ts
        int out_local=0, out_tzoffset=0

    # specify error conditions
    assert is_raise or is_ignore or is_coerce

    try:
        result = np.empty(n, dtype='M8[ns]')
        iresult = result.view('i8')
        for i in range(n):
            val = values[i]

            if checknull_with_nat(val):
                iresult[i] = NPY_NAT

            elif PyDateTime_Check(val):
                seen_datetime = 1
                if val.tzinfo is not None:
                    if utc_convert:
                        _ts = convert_datetime_to_tsobject(val, None)
                        iresult[i] = _ts.value
                        try:
                            check_dts_bounds(&_ts.dts)
                        except ValueError:
                            if is_coerce:
                                iresult[i] = NPY_NAT
                                continue
                            raise
                    else:
                        raise ValueError('Tz-aware datetime.datetime cannot '
                                         'be converted to datetime64 unless '
                                         'utc=True')
                else:
                    iresult[i] = pydatetime_to_dt64(val, &dts)
                    if is_timestamp(val):
                        iresult[i] += val.nanosecond
                    try:
                        check_dts_bounds(&dts)
                    except ValueError:
                        if is_coerce:
                            iresult[i] = NPY_NAT
                            continue
                        raise

            elif PyDate_Check(val):
                iresult[i] = pydate_to_dt64(val, &dts)
                try:
                    check_dts_bounds(&dts)
                    seen_datetime = 1
                except ValueError:
                    if is_coerce:
                        iresult[i] = NPY_NAT
                        continue
                    raise

            elif is_datetime64_object(val):
                if get_datetime64_value(val) == NPY_NAT:
                    iresult[i] = NPY_NAT
                else:
                    try:
                        iresult[i] = get_datetime64_nanos(val)
                        seen_datetime = 1
                    except ValueError:
                        if is_coerce:
                            iresult[i] = NPY_NAT
                            continue
                        raise

            elif is_integer_object(val) or is_float_object(val):
                # these must be ns unit by-definition

                if val != val or val == NPY_NAT:
                    iresult[i] = NPY_NAT
                elif is_raise or is_ignore:
                    iresult[i] = val
                    seen_integer = 1
                else:
                    # coerce
                    # we now need to parse this as if unit='ns'
                    # we can ONLY accept integers at this point
                    # if we have previously (or in future accept
                    # datetimes/strings, then we must coerce)
                    seen_integer = 1
                    try:
                        iresult[i] = cast_from_unit(val, 'ns')
                    except:
                        iresult[i] = NPY_NAT

            elif is_string_object(val):
                # string

                try:
                    if len(val) == 0 or val in nat_strings:
                        iresult[i] = NPY_NAT
                        continue

                    seen_string = 1
                    _string_to_dts(val, &dts, &out_local, &out_tzoffset)
                    value = dtstruct_to_dt64(&dts)
                    if out_local == 1:
                        tz = pytz.FixedOffset(out_tzoffset)
                        value = tz_convert_single(value, tz, 'UTC')
                    iresult[i] = value
                    check_dts_bounds(&dts)
                except ValueError:
                    # if requiring iso8601 strings, skip trying other formats
                    if require_iso8601:
                        if is_coerce:
                            iresult[i] = NPY_NAT
                            continue
                        elif is_raise:
                            raise ValueError(
                                "time data %r doesn't match format "
                                "specified" % (val,))
                        else:
                            return values

                    try:
                        py_dt = parse_datetime_string(val, dayfirst=dayfirst,
                                                      yearfirst=yearfirst)
                    except Exception:
                        if is_coerce:
                            iresult[i] = NPY_NAT
                            continue
                        raise TypeError("invalid string coercion to datetime")

                    try:
                        _ts = convert_datetime_to_tsobject(py_dt, None)
                        iresult[i] = _ts.value
                    except ValueError:
                        if is_coerce:
                            iresult[i] = NPY_NAT
                            continue
                        raise
                except:
                    if is_coerce:
                        iresult[i] = NPY_NAT
                        continue
                    raise
            else:
                if is_coerce:
                    iresult[i] = NPY_NAT
                else:
                    raise TypeError("{0} is not convertible to datetime"
                                    .format(type(val)))

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
                raise ValueError(
                    "mixed datetimes and integers in passed array")
            else:
                raise TypeError

        return result
    except OutOfBoundsDatetime:
        if is_raise:
            raise

        oresult = np.empty(n, dtype=object)
        for i in range(n):
            val = values[i]

            # set as nan except if its a NaT
            if checknull_with_nat(val):
                if PyFloat_Check(val):
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
    except TypeError:
        oresult = np.empty(n, dtype=object)

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
                except Exception:
                    if is_raise:
                        raise
                    return values
                    # oresult[i] = val
            else:
                if is_raise:
                    raise
                return values

        return oresult


# ----------------------------------------------------------------------
# Some general helper functions


def monthrange(int64_t year, int64_t month):
    cdef:
        int64_t days

    if month < 1 or month > 12:
        raise ValueError("bad month number 0; must be 1-12")

    days = days_per_month_table[is_leapyear(year)][month - 1]

    return (dayofweek(year, month, 1), days)


cpdef normalize_date(object dt):
    """
    Normalize datetime.datetime value to midnight. Returns datetime.date as a
    datetime.datetime at midnight

    Returns
    -------
    normalized : datetime.datetime or Timestamp
    """
    if is_timestamp(dt):
        return dt.replace(hour=0, minute=0, second=0, microsecond=0,
                          nanosecond=0)
    elif PyDateTime_Check(dt):
        return dt.replace(hour=0, minute=0, second=0, microsecond=0)
    elif PyDate_Check(dt):
        return datetime(dt.year, dt.month, dt.day)
    else:
        raise TypeError('Unrecognized type: %s' % type(dt))
