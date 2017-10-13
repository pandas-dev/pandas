# -*- coding: utf-8 -*-
# cython: profile=False
# cython: linetrace=False
# distutils: define_macros=CYTHON_TRACE=0
# distutils: define_macros=CYTHON_TRACE_NOGIL=0

cimport numpy as np
from numpy cimport (int8_t, int32_t, int64_t, import_array, ndarray,
                    float64_t, NPY_DATETIME, NPY_TIMEDELTA)
import numpy as np

import sys
cdef bint PY3 = (sys.version_info[0] >= 3)

from cpython cimport (
    PyTypeObject,
    PyFloat_Check,
    PyComplex_Check,
    PyObject_RichCompareBool,
    PyObject_RichCompare,
    Py_GT, Py_GE, Py_EQ, Py_NE, Py_LT, Py_LE,
    PyUnicode_Check)

cdef extern from "Python.h":
    cdef PyTypeObject *Py_TYPE(object)

from libc.stdlib cimport free

from util cimport (is_integer_object, is_float_object, is_datetime64_object,
                   is_timedelta64_object, INT64_MAX)
cimport util

# this is our datetime.pxd
from datetime cimport (
    pandas_datetimestruct,
    pandas_datetime_to_datetimestruct,
    pandas_datetimestruct_to_datetime,
    days_per_month_table,
    get_datetime64_value,
    get_timedelta64_value,
    get_datetime64_unit,
    PANDAS_DATETIMEUNIT,
    _string_to_dts,
    _pydatetime_to_dts,
    _date_to_datetime64,
    npy_datetime,
    is_leapyear,
    dayofweek,
    check_dts_bounds,
    PANDAS_FR_ns,
    PyDateTime_Check, PyDate_Check,
    PyDateTime_IMPORT,
    timedelta, datetime
    )

# stdlib datetime imports
from datetime import timedelta, datetime
from datetime import time as datetime_time

from khash cimport (
    khiter_t,
    kh_destroy_int64, kh_put_int64,
    kh_init_int64, kh_int64_t,
    kh_resize_int64, kh_get_int64)

from .tslibs.parsing import parse_datetime_string

cimport cython

from pandas.compat import iteritems, callable

import collections
import warnings

# initialize numpy
import_array()
#import_ufunc()

# import datetime C API
PyDateTime_IMPORT

cdef int64_t NPY_NAT = util.get_nat()
iNaT = NPY_NAT

from tslibs.timezones cimport (
    is_utc, is_tzlocal, is_fixed_offset,
    treat_tz_as_dateutil, treat_tz_as_pytz,
    get_timezone, get_utcoffset, maybe_get_tz,
    get_dst_info
    )
from tslibs.fields import (
    get_date_name_field, get_start_end_field, get_date_field,
    build_field_sarray)


cdef inline object create_timestamp_from_ts(
        int64_t value, pandas_datetimestruct dts,
        object tz, object freq):
    """ convenience routine to construct a Timestamp from its parts """
    cdef _Timestamp ts_base
    ts_base = _Timestamp.__new__(Timestamp, dts.year, dts.month,
                                 dts.day, dts.hour, dts.min,
                                 dts.sec, dts.us, tz)
    ts_base.value = value
    ts_base.freq = freq
    ts_base.nanosecond = dts.ps / 1000

    return ts_base


cdef inline object create_datetime_from_ts(
        int64_t value, pandas_datetimestruct dts,
        object tz, object freq):
    """ convenience routine to construct a datetime.datetime from its parts """
    return datetime(dts.year, dts.month, dts.day, dts.hour,
                    dts.min, dts.sec, dts.us, tz)


def ints_to_pydatetime(ndarray[int64_t] arr, tz=None, freq=None, box=False):
    # convert an i8 repr to an ndarray of datetimes or Timestamp (if box ==
    # True)

    cdef:
        Py_ssize_t i, n = len(arr)
        ndarray[int64_t] trans, deltas
        pandas_datetimestruct dts
        object dt
        int64_t value
        ndarray[object] result = np.empty(n, dtype=object)
        object (*func_create)(int64_t, pandas_datetimestruct, object, object)

    if box and util.is_string_object(freq):
        from pandas.tseries.frequencies import to_offset
        freq = to_offset(freq)

    if box:
        func_create = create_timestamp_from_ts
    else:
        func_create = create_datetime_from_ts

    if tz is not None:
        if is_utc(tz):
            for i in range(n):
                value = arr[i]
                if value == NPY_NAT:
                    result[i] = NaT
                else:
                    pandas_datetime_to_datetimestruct(
                        value, PANDAS_FR_ns, &dts)
                    result[i] = func_create(value, dts, tz, freq)
        elif is_tzlocal(tz) or is_fixed_offset(tz):
            for i in range(n):
                value = arr[i]
                if value == NPY_NAT:
                    result[i] = NaT
                else:
                    pandas_datetime_to_datetimestruct(
                        value, PANDAS_FR_ns, &dts)
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

                    pandas_datetime_to_datetimestruct(
                        value + deltas[pos], PANDAS_FR_ns, &dts)
                    result[i] = func_create(value, dts, new_tz, freq)
    else:
        for i in range(n):

            value = arr[i]
            if value == NPY_NAT:
                result[i] = NaT
            else:
                pandas_datetime_to_datetimestruct(value, PANDAS_FR_ns, &dts)
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


_zero_time = datetime_time(0, 0)
_no_input = object()

# Python front end to C extension type _Timestamp
# This serves as the box for datetime64


class Timestamp(_Timestamp):
    """TimeStamp is the pandas equivalent of python's Datetime
    and is interchangable with it in most cases. It's the type used
    for the entries that make up a DatetimeIndex, and other timeseries
    oriented data structures in pandas.

    There are essentially three calling conventions for the constructor. The
    primary form accepts four parameters. They can be passed by position or
    keyword.

    Parameters
    ----------
    ts_input : datetime-like, str, int, float
        Value to be converted to Timestamp
    freq : str, DateOffset
        Offset which Timestamp will have
    tz : string, pytz.timezone, dateutil.tz.tzfile or None
        Time zone for time which Timestamp will have.
    unit : string
        numpy unit used for conversion, if ts_input is int or float
    offset : str, DateOffset
        Deprecated, use freq

    The other two forms mimic the parameters from ``datetime.datetime``. They
    can be passed by either position or keyword, but not both mixed together.

    :func:`datetime.datetime` Parameters
    ------------------------------------

    .. versionadded:: 0.19.0

    year : int
    month : int
    day : int
    hour : int, optional, default is 0
    minute : int, optional, default is 0
    second : int, optional, default is 0
    microsecond : int, optional, default is 0
    tzinfo : datetime.tzinfo, optional, default is None
    """

    @classmethod
    def fromordinal(cls, ordinal, freq=None, tz=None, offset=None):
        """
        passed an ordinal, translate and convert to a ts
        note: by definition there cannot be any tz info on the ordinal itself

        Parameters
        ----------
        ordinal : int
            date corresponding to a proleptic Gregorian ordinal
        freq : str, DateOffset
            Offset which Timestamp will have
        tz : string, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for time which Timestamp will have.
        offset : str, DateOffset
            Deprecated, use freq
        """
        return cls(datetime.fromordinal(ordinal),
                   freq=freq, tz=tz, offset=offset)

    @classmethod
    def now(cls, tz=None):
        """
        Return the current time in the local timezone.  Equivalent
        to datetime.now([tz])

        Parameters
        ----------
        tz : string / timezone object, default None
            Timezone to localize to
        """
        if util.is_string_object(tz):
            tz = maybe_get_tz(tz)
        return cls(datetime.now(tz))

    @classmethod
    def today(cls, tz=None):
        """
        Return the current time in the local timezone.  This differs
        from datetime.today() in that it can be localized to a
        passed timezone.

        Parameters
        ----------
        tz : string / timezone object, default None
            Timezone to localize to
        """
        return cls.now(tz)

    @classmethod
    def utcnow(cls):
        return cls.now('UTC')

    @classmethod
    def utcfromtimestamp(cls, ts):
        return cls(datetime.utcfromtimestamp(ts))

    @classmethod
    def fromtimestamp(cls, ts):
        return cls(datetime.fromtimestamp(ts))

    @classmethod
    def combine(cls, date, time):
        return cls(datetime.combine(date, time))

    def __new__(cls, object ts_input=_no_input,
                object freq=None, tz=None, unit=None,
                year=None, month=None, day=None,
                hour=None, minute=None, second=None, microsecond=None,
                tzinfo=None,
                object offset=None):
        # The parameter list folds together legacy parameter names (the first
        # four) and positional and keyword parameter names from pydatetime.
        #
        # There are three calling forms:
        #
        # - In the legacy form, the first parameter, ts_input, is required
        #   and may be datetime-like, str, int, or float. The second
        #   parameter, offset, is optional and may be str or DateOffset.
        #
        # - ints in the first, second, and third arguments indicate
        #   pydatetime positional arguments. Only the first 8 arguments
        #   (standing in for year, month, day, hour, minute, second,
        #   microsecond, tzinfo) may be non-None. As a shortcut, we just
        #   check that the second argument is an int.
        #
        # - Nones for the first four (legacy) arguments indicate pydatetime
        #   keyword arguments. year, month, and day are required. As a
        #   shortcut, we just check that the first argument was not passed.
        #
        # Mixing pydatetime positional and keyword arguments is forbidden!

        cdef _TSObject ts

        if offset is not None:
            # deprecate offset kwd in 0.19.0, GH13593
            if freq is not None:
                msg = "Can only specify freq or offset, not both"
                raise TypeError(msg)
            warnings.warn("offset is deprecated. Use freq instead",
                          FutureWarning)
            freq = offset

        if ts_input is _no_input:
            # User passed keyword arguments.
            return Timestamp(datetime(year, month, day, hour or 0,
                                      minute or 0, second or 0,
                                      microsecond or 0, tzinfo),
                             tz=tzinfo)
        elif is_integer_object(freq):
            # User passed positional arguments:
            # Timestamp(year, month, day[, hour[, minute[, second[,
            # microsecond[, tzinfo]]]]])
            return Timestamp(datetime(ts_input, freq, tz, unit or 0,
                                      year or 0, month or 0, day or 0,
                                      hour), tz=hour)

        ts = convert_to_tsobject(ts_input, tz, unit, 0, 0)

        if ts.value == NPY_NAT:
            return NaT

        if util.is_string_object(freq):
            from pandas.tseries.frequencies import to_offset
            freq = to_offset(freq)

        return create_timestamp_from_ts(ts.value, ts.dts, ts.tzinfo, freq)

    def _round(self, freq, rounder):

        cdef:
            int64_t unit, r, value,  buff = 1000000
            object result

        from pandas.tseries.frequencies import to_offset
        unit = to_offset(freq).nanos
        if self.tz is not None:
            value = self.tz_localize(None).value
        else:
            value = self.value
        if unit < 1000 and unit % 1000 != 0:
            # for nano rounding, work with the last 6 digits separately
            # due to float precision
            r = (buff * (value // buff) + unit *
                 (rounder((value % buff) / float(unit))).astype('i8'))
        elif unit >= 1000 and unit % 1000 != 0:
            msg = 'Precision will be lost using frequency: {}'
            warnings.warn(msg.format(freq))
            r = (unit * rounder(value / float(unit)).astype('i8'))
        else:
            r = (unit * rounder(value / float(unit)).astype('i8'))
        result = Timestamp(r, unit='ns')
        if self.tz is not None:
            result = result.tz_localize(self.tz)
        return result

    def round(self, freq):
        """
        Round the Timestamp to the specified resolution

        Returns
        -------
        a new Timestamp rounded to the given resolution of `freq`

        Parameters
        ----------
        freq : a freq string indicating the rounding resolution

        Raises
        ------
        ValueError if the freq cannot be converted
        """
        return self._round(freq, np.round)

    def floor(self, freq):
        """
        return a new Timestamp floored to this resolution

        Parameters
        ----------
        freq : a freq string indicating the flooring resolution
        """
        return self._round(freq, np.floor)

    def ceil(self, freq):
        """
        return a new Timestamp ceiled to this resolution

        Parameters
        ----------
        freq : a freq string indicating the ceiling resolution
        """
        return self._round(freq, np.ceil)

    @property
    def tz(self):
        """
        Alias for tzinfo
        """
        return self.tzinfo

    @property
    def offset(self):
        warnings.warn(".offset is deprecated. Use .freq instead",
                      FutureWarning)
        return self.freq

    def __setstate__(self, state):
        self.value = state[0]
        self.freq = state[1]
        self.tzinfo = state[2]

    def __reduce__(self):
        object_state = self.value, self.freq, self.tzinfo
        return (Timestamp, object_state)

    def to_period(self, freq=None):
        """
        Return an period of which this timestamp is an observation.
        """
        from pandas import Period

        if freq is None:
            freq = self.freq

        return Period(self, freq=freq)

    @property
    def dayofweek(self):
        return self.weekday()

    @property
    def weekday_name(self):
        return self._get_named_field('weekday_name')

    @property
    def dayofyear(self):
        return self._get_field('doy')

    @property
    def week(self):
        return self._get_field('woy')

    weekofyear = week

    @property
    def quarter(self):
        return self._get_field('q')

    @property
    def days_in_month(self):
        return self._get_field('dim')

    daysinmonth = days_in_month

    @property
    def freqstr(self):
        return getattr(self.freq, 'freqstr', self.freq)

    @property
    def is_month_start(self):
        return self._get_start_end_field('is_month_start')

    @property
    def is_month_end(self):
        return self._get_start_end_field('is_month_end')

    @property
    def is_quarter_start(self):
        return self._get_start_end_field('is_quarter_start')

    @property
    def is_quarter_end(self):
        return self._get_start_end_field('is_quarter_end')

    @property
    def is_year_start(self):
        return self._get_start_end_field('is_year_start')

    @property
    def is_year_end(self):
        return self._get_start_end_field('is_year_end')

    @property
    def is_leap_year(self):
        return bool(is_leapyear(self.year))

    def tz_localize(self, tz, ambiguous='raise', errors='raise'):
        """
        Convert naive Timestamp to local time zone, or remove
        timezone from tz-aware Timestamp.

        Parameters
        ----------
        tz : string, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for time which Timestamp will be converted to.
            None will remove timezone holding local time.
        ambiguous : bool, 'NaT', default 'raise'
            - bool contains flags to determine if time is dst or not (note
            that this flag is only applicable for ambiguous fall dst dates)
            - 'NaT' will return NaT for an ambiguous time
            - 'raise' will raise an AmbiguousTimeError for an ambiguous time
        errors : 'raise', 'coerce', default 'raise'
            - 'raise' will raise a NonExistentTimeError if a timestamp is not
               valid in the specified timezone (e.g. due to a transition from
               or to DST time)
            - 'coerce' will return NaT if the timestamp can not be converted
              into the specified timezone

              .. versionadded:: 0.19.0

        Returns
        -------
        localized : Timestamp

        Raises
        ------
        TypeError
            If the Timestamp is tz-aware and tz is not None.
        """
        if ambiguous == 'infer':
            raise ValueError('Cannot infer offset with only one time.')

        if self.tzinfo is None:
            # tz naive, localize
            tz = maybe_get_tz(tz)
            if not util.is_string_object(ambiguous):
                ambiguous =   [ambiguous]
            value = tz_localize_to_utc(np.array([self.value], dtype='i8'), tz,
                                       ambiguous=ambiguous, errors=errors)[0]
            return Timestamp(value, tz=tz)
        else:
            if tz is None:
                # reset tz
                value = tz_convert_single(self.value, 'UTC', self.tz)
                return Timestamp(value, tz=None)
            else:
                raise TypeError('Cannot localize tz-aware Timestamp, use '
                                'tz_convert for conversions')

    def tz_convert(self, tz):
        """
        Convert tz-aware Timestamp to another time zone.

        Parameters
        ----------
        tz : string, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for time which Timestamp will be converted to.
            None will remove timezone holding UTC time.

        Returns
        -------
        converted : Timestamp

        Raises
        ------
        TypeError
            If Timestamp is tz-naive.
        """
        if self.tzinfo is None:
            # tz naive, use tz_localize
            raise TypeError('Cannot convert tz-naive Timestamp, use '
                            'tz_localize to localize')
        else:
            # Same UTC timestamp, different time zone
            return Timestamp(self.value, tz=tz)

    astimezone = tz_convert

    def replace(self, year=None, month=None, day=None,
                hour=None, minute=None, second=None, microsecond=None,
                nanosecond=None, tzinfo=object, fold=0):
        """
        implements datetime.replace, handles nanoseconds

        Parameters
        ----------
        year : int, optional
        month : int, optional
        day : int, optional
        hour : int, optional
        minute : int, optional
        second : int, optional
        microsecond : int, optional
        nanosecond: int, optional
        tzinfo : tz-convertible, optional
        fold : int, optional, default is 0
            added in 3.6, NotImplemented

        Returns
        -------
        Timestamp with fields replaced
        """

        cdef:
            pandas_datetimestruct dts
            int64_t value, value_tz, offset
            object _tzinfo, result, k, v
            datetime ts_input

        # set to naive if needed
        _tzinfo = self.tzinfo
        value = self.value
        if _tzinfo is not None:
            value_tz = tz_convert_single(value, _tzinfo, 'UTC')
            value += value - value_tz

        # setup components
        pandas_datetime_to_datetimestruct(value, PANDAS_FR_ns, &dts)
        dts.ps = self.nanosecond * 1000

        # replace
        def validate(k, v):
            """ validate integers """
            if not is_integer_object(v):
                raise ValueError("value must be an integer, received "
                                 "{v} for {k}".format(v=type(v), k=k))
            return v

        if year is not None:
            dts.year = validate('year', year)
        if month is not None:
            dts.month = validate('month', month)
        if day is not None:
            dts.day = validate('day', day)
        if hour is not None:
            dts.hour = validate('hour', hour)
        if minute is not None:
            dts.min = validate('minute', minute)
        if second is not None:
            dts.sec = validate('second', second)
        if microsecond is not None:
            dts.us = validate('microsecond', microsecond)
        if nanosecond is not None:
            dts.ps = validate('nanosecond', nanosecond) * 1000
        if tzinfo is not object:
            _tzinfo = tzinfo

        # reconstruct & check bounds
        ts_input = datetime(dts.year, dts.month, dts.day, dts.hour, dts.min,
                            dts.sec, dts.us, tzinfo=_tzinfo)
        ts = convert_datetime_to_tsobject(ts_input, _tzinfo)
        value = ts.value + (dts.ps // 1000)
        if value != NPY_NAT:
            _check_dts_bounds(&dts)

        return create_timestamp_from_ts(value, dts, _tzinfo, self.freq)

    def isoformat(self, sep='T'):
        base = super(_Timestamp, self).isoformat(sep=sep)
        if self.nanosecond == 0:
            return base

        if self.tzinfo is not None:
            base1, base2 = base[:-6], base[-6:]
        else:
            base1, base2 = base, ""

        if self.microsecond != 0:
            base1 += "%.3d" % self.nanosecond
        else:
            base1 += ".%.9d" % self.nanosecond

        return base1 + base2

    def _has_time_component(self):
        """
        Returns if the Timestamp has a time component
        in addition to the date part
        """
        return (self.time() != _zero_time
                or self.tzinfo is not None
                or self.nanosecond != 0)

    def to_julian_date(self):
        """
        Convert TimeStamp to a Julian Date.
        0 Julian date is noon January 1, 4713 BC.
        """
        year = self.year
        month = self.month
        day = self.day
        if month <= 2:
            year -= 1
            month += 12
        return (day +
                np.fix((153 * month - 457) / 5) +
                365 * year +
                np.floor(year / 4) -
                np.floor(year / 100) +
                np.floor(year / 400) +
                1721118.5 +
                (self.hour +
                 self.minute / 60.0 +
                 self.second / 3600.0 +
                 self.microsecond / 3600.0 / 1e+6 +
                 self.nanosecond / 3600.0 / 1e+9
                ) / 24.0)

    def normalize(self):
        """
        Normalize Timestamp to midnight, preserving
        tz information.
        """
        normalized_value = date_normalize(
            np.array([self.value], dtype='i8'), tz=self.tz)[0]
        return Timestamp(normalized_value).tz_localize(self.tz)

    def __radd__(self, other):
        # __radd__ on cython extension types like _Timestamp is not used, so
        # define it here instead
        return self + other


_nat_strings = set(['NaT', 'nat', 'NAT', 'nan', 'NaN', 'NAN'])


def _make_nat_func(func_name, cls):
    def f(*args, **kwargs):
        return NaT
    f.__name__ = func_name
    f.__doc__ = getattr(cls, func_name).__doc__
    return f


def _make_nan_func(func_name, cls):
    def f(*args, **kwargs):
        return np.nan
    f.__name__ = func_name
    f.__doc__ = getattr(cls, func_name).__doc__
    return f


def _make_error_func(func_name, cls):
    def f(*args, **kwargs):
        raise ValueError("NaTType does not support " + func_name)

    f.__name__ = func_name
    if cls is not None:
        f.__doc__ = getattr(cls, func_name).__doc__
    return f


class NaTType(_NaT):
    """(N)ot-(A)-(T)ime, the time equivalent of NaN"""

    def __new__(cls):
        cdef _NaT base

        base = _NaT.__new__(cls, 1, 1, 1)
        base.value = NPY_NAT

        return base

    def __repr__(self):
        return 'NaT'

    def __str__(self):
        return 'NaT'

    def isoformat(self, sep='T'):
        # This allows Timestamp(ts.isoformat()) to always correctly roundtrip.
        return 'NaT'

    def __hash__(self):
        return NPY_NAT

    def __int__(self):
        return NPY_NAT

    def __long__(self):
        return NPY_NAT

    def __reduce__(self):
        return (__nat_unpickle, (None, ))

    def total_seconds(self):
        """
        Total duration of timedelta in seconds (to ns precision)
        """
        # GH 10939
        return np.nan

    @property
    def is_leap_year(self):
        return False

    @property
    def is_month_start(self):
        return False

    @property
    def is_quarter_start(self):
        return False

    @property
    def is_year_start(self):
        return False

    @property
    def is_month_end(self):
        return False

    @property
    def is_quarter_end(self):
        return False

    @property
    def is_year_end(self):
        return False

    def __rdiv__(self, other):
        return _nat_rdivide_op(self, other)

    def __rtruediv__(self, other):
        return _nat_rdivide_op(self, other)

    def __rfloordiv__(self, other):
        return _nat_rdivide_op(self, other)

    def __rmul__(self, other):
        if is_integer_object(other) or is_float_object(other):
            return NaT
        return NotImplemented

    # ----------------------------------------------------------------------
    # inject the Timestamp field properties
    # these by definition return np.nan

    year = property(fget=lambda self: np.nan)
    quarter = property(fget=lambda self: np.nan)
    month = property(fget=lambda self: np.nan)
    day = property(fget=lambda self: np.nan)
    hour = property(fget=lambda self: np.nan)
    minute = property(fget=lambda self: np.nan)
    second = property(fget=lambda self: np.nan)
    millisecond = property(fget=lambda self: np.nan)
    microsecond = property(fget=lambda self: np.nan)
    nanosecond = property(fget=lambda self: np.nan)

    week = property(fget=lambda self: np.nan)
    dayofyear = property(fget=lambda self: np.nan)
    weekofyear = property(fget=lambda self: np.nan)
    days_in_month = property(fget=lambda self: np.nan)
    daysinmonth = property(fget=lambda self: np.nan)
    dayofweek = property(fget=lambda self: np.nan)
    weekday_name = property(fget=lambda self: np.nan)

    # inject Timedelta properties
    days = property(fget=lambda self: np.nan)
    seconds = property(fget=lambda self: np.nan)
    microseconds = property(fget=lambda self: np.nan)
    nanoseconds = property(fget=lambda self: np.nan)

    # inject pd.Period properties
    qyear = property(fget=lambda self: np.nan)

    # ----------------------------------------------------------------------
    # GH9513 NaT methods (except to_datetime64) to raise, return np.nan, or
    # return NaT create functions that raise, for binding to NaTType
    # These are the ones that can get their docstrings from datetime.

    # nan methods
    weekday = _make_nan_func('weekday', datetime)
    isoweekday = _make_nan_func('isoweekday', datetime)

    # _nat_methods
    date = _make_nat_func('date', datetime)

    utctimetuple = _make_error_func('utctimetuple', datetime)
    timetz = _make_error_func('timetz', datetime)
    timetuple = _make_error_func('timetuple', datetime)
    strptime = _make_error_func('strptime', datetime)
    strftime = _make_error_func('strftime', datetime)
    isocalendar = _make_error_func('isocalendar', datetime)
    dst = _make_error_func('dst', datetime)
    ctime = _make_error_func('ctime', datetime)
    time = _make_error_func('time', datetime)
    toordinal = _make_error_func('toordinal', datetime)
    tzname = _make_error_func('tzname', datetime)
    utcoffset = _make_error_func('utcoffset', datetime)

    # Timestamp has empty docstring for some methods.
    utcfromtimestamp = _make_error_func('utcfromtimestamp', None) 
    fromtimestamp = _make_error_func('fromtimestamp', None)
    combine = _make_error_func('combine', None)
    utcnow = _make_error_func('utcnow', None)

    if PY3:
        timestamp = _make_error_func('timestamp', datetime)

    # GH9513 NaT methods (except to_datetime64) to raise, return np.nan, or
    # return NaT create functions that raise, for binding to NaTType
    astimezone = _make_error_func('astimezone', Timestamp)
    fromordinal = _make_error_func('fromordinal', Timestamp)

    # _nat_methods
    to_pydatetime = _make_nat_func('to_pydatetime', Timestamp)

    now = _make_nat_func('now', Timestamp)
    today = _make_nat_func('today', Timestamp)
    round = _make_nat_func('round', Timestamp)
    floor = _make_nat_func('floor', Timestamp)
    ceil = _make_nat_func('ceil', Timestamp)

    tz_convert = _make_nat_func('tz_convert', Timestamp)
    tz_localize = _make_nat_func('tz_localize', Timestamp)
    replace = _make_nat_func('replace', Timestamp)


def __nat_unpickle(*args):
    # return constant defined in the module
    return NaT

NaT = NaTType()

cdef inline bint _checknull_with_nat(object val):
    """ utility to check if a value is a nat or not """
    return val is None or (
        PyFloat_Check(val) and val != val) or val is NaT

cdef inline bint _check_all_nulls(object val):
    """ utility to check if a value is any type of null """
    cdef bint res
    if PyFloat_Check(val) or PyComplex_Check(val):
        res = val != val
    elif val is NaT:
        res = 1
    elif val is None:
        res = 1
    elif is_datetime64_object(val):
        res = get_datetime64_value(val) == NPY_NAT
    elif is_timedelta64_object(val):
        res = get_timedelta64_value(val) == NPY_NAT
    else:
        res = 0
    return res

cdef inline bint _cmp_nat_dt(_NaT lhs, _Timestamp rhs, int op) except -1:
    return _nat_scalar_rules[op]


cpdef object get_value_box(ndarray arr, object loc):
    cdef:
        Py_ssize_t i, sz
        void* data_ptr

    if util.is_float_object(loc):
        casted = int(loc)
        if casted == loc:
            loc = casted
    i = <Py_ssize_t> loc
    sz = np.PyArray_SIZE(arr)

    if i < 0 and sz > 0:
        i += sz

    if i >= sz or sz == 0 or i < 0:
        raise IndexError('index out of bounds')

    if arr.descr.type_num == NPY_DATETIME:
        return Timestamp(util.get_value_1d(arr, i))
    elif arr.descr.type_num == NPY_TIMEDELTA:
        return Timedelta(util.get_value_1d(arr, i))
    else:
        return util.get_value_1d(arr, i)


# Add the min and max fields at the class level
cdef int64_t _NS_UPPER_BOUND = INT64_MAX
# the smallest value we could actually represent is
#   INT64_MIN + 1 == -9223372036854775807
# but to allow overflow free conversion with a microsecond resolution
# use the smallest value with a 0 nanosecond unit (0s in last 3 digits)
cdef int64_t _NS_LOWER_BOUND = -9223372036854775000

# Resolution is in nanoseconds
Timestamp.min = Timestamp(_NS_LOWER_BOUND)
Timestamp.max = Timestamp(_NS_UPPER_BOUND)


#----------------------------------------------------------------------
# Frequency inference

def unique_deltas(ndarray[int64_t] arr):
    cdef:
        Py_ssize_t i, n = len(arr)
        int64_t val
        khiter_t k
        kh_int64_t *table
        int ret = 0
        list uniques = []

    table = kh_init_int64()
    kh_resize_int64(table, 10)
    for i in range(n - 1):
        val = arr[i + 1] - arr[i]
        k = kh_get_int64(table, val)
        if k == table.n_buckets:
            kh_put_int64(table, val, &ret)
            uniques.append(val)
    kh_destroy_int64(table)

    result = np.array(uniques, dtype=np.int64)
    result.sort()
    return result


cdef inline bint _cmp_scalar(int64_t lhs, int64_t rhs, int op) except -1:
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


cdef int _reverse_ops[6]

_reverse_ops[Py_LT] = Py_GT
_reverse_ops[Py_LE] = Py_GE
_reverse_ops[Py_EQ] = Py_EQ
_reverse_ops[Py_NE] = Py_NE
_reverse_ops[Py_GT] = Py_LT
_reverse_ops[Py_GE] = Py_LE


cdef str _NDIM_STRING = "ndim"

# This is PITA. Because we inherit from datetime, which has very specific
# construction requirements, we need to do object instantiation in python
# (see Timestamp class above). This will serve as a C extension type that
# shadows the python class, where we do any heavy lifting.
cdef class _Timestamp(datetime):

    cdef readonly:
        int64_t value, nanosecond
        object freq       # frequency reference

    def __hash__(_Timestamp self):
        if self.nanosecond:
            return hash(self.value)
        return datetime.__hash__(self)

    def __richcmp__(_Timestamp self, object other, int op):
        cdef:
            _Timestamp ots
            int ndim

        if isinstance(other, _Timestamp):
            if isinstance(other, _NaT):
                return _cmp_nat_dt(other, self, _reverse_ops[op])
            ots = other
        elif isinstance(other, datetime):
            if self.nanosecond == 0:
                val = self.to_pydatetime()
                return PyObject_RichCompareBool(val, other, op)

            try:
                ots = Timestamp(other)
            except ValueError:
                return self._compare_outside_nanorange(other, op)
        else:
            ndim = getattr(other, _NDIM_STRING, -1)

            if ndim != -1:
                if ndim == 0:
                    if isinstance(other, np.datetime64):
                        other = Timestamp(other)
                    else:
                        if op == Py_EQ:
                            return False
                        elif op == Py_NE:
                            return True

                        # only allow ==, != ops
                        raise TypeError('Cannot compare type %r with type %r' %
                                        (type(self).__name__,
                                         type(other).__name__))
                return PyObject_RichCompare(other, self, _reverse_ops[op])
            else:
                if op == Py_EQ:
                    return False
                elif op == Py_NE:
                    return True
                raise TypeError('Cannot compare type %r with type %r' %
                                (type(self).__name__, type(other).__name__))

        self._assert_tzawareness_compat(other)
        return _cmp_scalar(self.value, ots.value, op)

    def __reduce_ex__(self, protocol):
        # python 3.6 compat
        # http://bugs.python.org/issue28730
        # now __reduce_ex__ is defined and higher priority than __reduce__
        return self.__reduce__()

    def __repr__(self):
        stamp = self._repr_base
        zone = None

        try:
            stamp += self.strftime('%z')
            if self.tzinfo:
                zone = get_timezone(self.tzinfo)
        except ValueError:
            year2000 = self.replace(year=2000)
            stamp += year2000.strftime('%z')
            if self.tzinfo:
                zone = get_timezone(self.tzinfo)

        try:
            stamp += zone.strftime(' %%Z')
        except:
            pass

        tz = ", tz='{0}'".format(zone) if zone is not None else ""
        freq = ", freq='{0}'".format(
            self.freq.freqstr) if self.freq is not None else ""

        return "Timestamp('{stamp}'{tz}{freq})".format(
            stamp=stamp, tz=tz, freq=freq)

    cdef bint _compare_outside_nanorange(_Timestamp self, datetime other,
                                         int op) except -1:
        cdef datetime dtval = self.to_pydatetime()

        self._assert_tzawareness_compat(other)

        if self.nanosecond == 0:
            return PyObject_RichCompareBool(dtval, other, op)
        else:
            if op == Py_EQ:
                return False
            elif op == Py_NE:
                return True
            elif op == Py_LT:
                return dtval < other
            elif op == Py_LE:
                return dtval < other
            elif op == Py_GT:
                return dtval >= other
            elif op == Py_GE:
                return dtval >= other

    cdef int _assert_tzawareness_compat(_Timestamp self,
                                        object other) except -1:
        if self.tzinfo is None:
            if other.tzinfo is not None:
                raise TypeError('Cannot compare tz-naive and tz-aware '
                                 'timestamps')
        elif other.tzinfo is None:
            raise TypeError('Cannot compare tz-naive and tz-aware timestamps')

    cpdef datetime to_datetime(_Timestamp self):
        """
        DEPRECATED: use :meth:`to_pydatetime` instead.

        Convert a Timestamp object to a native Python datetime object.
        """
        warnings.warn("to_datetime is deprecated. Use self.to_pydatetime()",
                      FutureWarning, stacklevel=2)
        return self.to_pydatetime(warn=False)

    cpdef datetime to_pydatetime(_Timestamp self, warn=True):
        """
        Convert a Timestamp object to a native Python datetime object.

        If warn=True, issue a warning if nanoseconds is nonzero.
        """
        if self.nanosecond != 0 and warn:
            warnings.warn("Discarding nonzero nanoseconds in conversion",
                          UserWarning, stacklevel=2)

        return datetime(self.year, self.month, self.day,
                        self.hour, self.minute, self.second,
                        self.microsecond, self.tzinfo)

    cpdef to_datetime64(self):
        """ Returns a numpy.datetime64 object with 'ns' precision """
        return np.datetime64(self.value, 'ns')

    def __add__(self, other):
        cdef int64_t other_int, nanos

        if is_timedelta64_object(other):
            other_int = other.astype('timedelta64[ns]').view('i8')
            return Timestamp(self.value + other_int,
                             tz=self.tzinfo, freq=self.freq)

        elif is_integer_object(other):
            if self is NaT:
                # to be compat with Period
                return NaT
            elif self.freq is None:
                raise ValueError("Cannot add integral value to Timestamp "
                                 "without freq.")
            return Timestamp((self.freq * other).apply(self), freq=self.freq)

        elif isinstance(other, timedelta) or hasattr(other, 'delta'):
            nanos = _delta_to_nanoseconds(other)
            result = Timestamp(self.value + nanos,
                               tz=self.tzinfo, freq=self.freq)
            if getattr(other, 'normalize', False):
                result = Timestamp(normalize_date(result))
            return result

        # index/series like
        elif hasattr(other, '_typ'):
            return NotImplemented

        result = datetime.__add__(self, other)
        if isinstance(result, datetime):
            result = Timestamp(result)
            result.nanosecond = self.nanosecond
        return result

    def __sub__(self, other):
        if is_timedelta64_object(other) or is_integer_object(other) \
                or isinstance(other, timedelta) or hasattr(other, 'delta'):
            neg_other = -other
            return self + neg_other

        # a Timestamp-DatetimeIndex -> yields a negative TimedeltaIndex
        elif getattr(other, '_typ', None) == 'datetimeindex':
            # timezone comparison is performed in DatetimeIndex._sub_datelike
            return -other.__sub__(self)

        # a Timestamp-TimedeltaIndex -> yields a negative TimedeltaIndex
        elif getattr(other, '_typ', None) == 'timedeltaindex':
            return (-other).__add__(self)

        elif other is NaT:
            return NaT

        # coerce if necessary if we are a Timestamp-like
        if (isinstance(self, datetime)
            and (isinstance(other, datetime)
                 or is_datetime64_object(other))):
            self = Timestamp(self)
            other = Timestamp(other)

            # validate tz's
            if get_timezone(self.tzinfo) != get_timezone(other.tzinfo):
                raise TypeError(
                    "Timestamp subtraction must have the "
                    "same timezones or no timezones")

            # scalar Timestamp/datetime - Timestamp/datetime -> yields a
            # Timedelta
            try:
                return Timedelta(self.value -other.value)
            except (OverflowError, OutOfBoundsDatetime):
                pass

        # scalar Timestamp/datetime - Timedelta -> yields a Timestamp (with
        # same timezone if specified)
        return datetime.__sub__(self, other)

    cdef int64_t _maybe_convert_value_to_local(self):
        """Convert UTC i8 value to local i8 value if tz exists"""
        cdef:
            int64_t val
        val = self.value
        if self.tz is not None and not is_utc(self.tz):
            val = tz_convert_single(self.value, 'UTC', self.tz)
        return val

    cpdef _get_field(self, field):
        cdef:
            int64_t val
            ndarray[int32_t] out
        val = self._maybe_convert_value_to_local()
        out = get_date_field(np.array([val], dtype=np.int64), field)
        return int(out[0])

    cpdef _get_named_field(self, field):
        cdef:
            int64_t val
            ndarray[object] out
        val = self._maybe_convert_value_to_local()
        out = get_date_name_field(np.array([val], dtype=np.int64), field)
        return out[0]

    cpdef _get_start_end_field(self, field):
        month_kw = self.freq.kwds.get(
            'startingMonth', self.freq.kwds.get(
                'month', 12)) if self.freq else 12
        freqstr = self.freqstr if self.freq else None
        val = self._maybe_convert_value_to_local()
        out = get_start_end_field(
            np.array([val], dtype=np.int64), field, freqstr, month_kw)
        return out[0]

    property _repr_base:
        def __get__(self):
            return '%s %s' % (self._date_repr, self._time_repr)

    property _date_repr:
        def __get__(self):
            # Ideal here would be self.strftime("%Y-%m-%d"), but
            # the datetime strftime() methods require year >= 1900
            return '%d-%.2d-%.2d' % (self.year, self.month, self.day)

    property _time_repr:
        def __get__(self):
            result = '%.2d:%.2d:%.2d' % (self.hour, self.minute, self.second)

            if self.nanosecond != 0:
                result += '.%.9d' % (self.nanosecond + 1000 * self.microsecond)
            elif self.microsecond != 0:
                result += '.%.6d' % self.microsecond

            return result

    property _short_repr:
        def __get__(self):
            # format a Timestamp with only _date_repr if possible
            # otherwise _repr_base
            if (self.hour == 0 and
                self.minute == 0 and
                self.second == 0 and
                self.microsecond == 0 and
                self.nanosecond == 0):
                return self._date_repr
            return self._repr_base

    property asm8:
        def __get__(self):
            return np.datetime64(self.value, 'ns')


cdef PyTypeObject* ts_type = <PyTypeObject*> Timestamp


cdef inline bint is_timestamp(object o):
    return Py_TYPE(o) == ts_type # isinstance(o, Timestamp)


cdef bint _nat_scalar_rules[6]

_nat_scalar_rules[Py_EQ] = False
_nat_scalar_rules[Py_NE] = True
_nat_scalar_rules[Py_LT] = False
_nat_scalar_rules[Py_LE] = False
_nat_scalar_rules[Py_GT] = False
_nat_scalar_rules[Py_GE] = False


cdef _nat_divide_op(self, other):
    if isinstance(other, (Timedelta, np.timedelta64)) or other is NaT:
        return np.nan
    if is_integer_object(other) or is_float_object(other):
        return NaT
    return NotImplemented

cdef _nat_rdivide_op(self, other):
    if isinstance(other, Timedelta):
        return np.nan
    return NotImplemented


cdef class _NaT(_Timestamp):

    def __hash__(_NaT self):
        # py3k needs this defined here
        return hash(self.value)

    def __richcmp__(_NaT self, object other, int op):
        cdef int ndim = getattr(other, 'ndim', -1)

        if ndim == -1:
            return _nat_scalar_rules[op]

        if ndim == 0:
            if isinstance(other, np.datetime64):
                other = Timestamp(other)
            else:
                raise TypeError('Cannot compare type %r with type %r' %
                                (type(self).__name__, type(other).__name__))
        return PyObject_RichCompare(other, self, _reverse_ops[op])

    def __add__(self, other):
        try:
            if isinstance(other, datetime):
                return NaT
            result = _Timestamp.__add__(self, other)
            # Timestamp.__add__ doesn't return DatetimeIndex/TimedeltaIndex
            if result is NotImplemented:
                return result
        except (OverflowError, OutOfBoundsDatetime):
            pass
        return NaT

    def __sub__(self, other):
        if isinstance(other, (datetime, timedelta)):
            return NaT
        try:
            result = _Timestamp.__sub__(self, other)
            # Timestamp.__sub__ may return DatetimeIndex/TimedeltaIndex
            if result is NotImplemented or hasattr(result, '_typ'):
                return result
        except (OverflowError, OutOfBoundsDatetime):
            pass
        return NaT

    def __pos__(self):
        return NaT

    def __neg__(self):
        return NaT

    def __div__(self, other):
        return _nat_divide_op(self, other)

    def __truediv__(self, other):
        return _nat_divide_op(self, other)

    def __floordiv__(self, other):
        return _nat_divide_op(self, other)

    def __mul__(self, other):
        if is_integer_object(other) or is_float_object(other):
            return NaT
        return NotImplemented


# lightweight C object to hold datetime & int64 pair
cdef class _TSObject:
    cdef:
        pandas_datetimestruct dts      # pandas_datetimestruct
        int64_t value               # numpy dt64
        object tzinfo

    property value:
        def __get__(self):
            return self.value


# helper to extract datetime and int64 from several different possibilities
cdef convert_to_tsobject(object ts, object tz, object unit,
                         bint dayfirst, bint yearfirst):
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
        bint utc_convert = 1
        int out_local = 0, out_tzoffset = 0

    if tz is not None:
        tz = maybe_get_tz(tz)

    obj = _TSObject()

    if util.is_string_object(ts):
        return convert_str_to_tsobject(ts, tz, unit, dayfirst, yearfirst)

    if ts is None or ts is NaT:
        obj.value = NPY_NAT
    elif is_datetime64_object(ts):
        if ts.view('i8') == NPY_NAT:
            obj.value = NPY_NAT
        else:
            obj.value = _get_datetime64_nanos(ts)
            pandas_datetime_to_datetimestruct(
                obj.value, PANDAS_FR_ns, &obj.dts)
    elif is_integer_object(ts):
        if ts == NPY_NAT:
            obj.value = NPY_NAT
        else:
            ts = ts * cast_from_unit(None, unit)
            obj.value = ts
            pandas_datetime_to_datetimestruct(ts, PANDAS_FR_ns, &obj.dts)
    elif util.is_float_object(ts):
        if ts != ts or ts == NPY_NAT:
            obj.value = NPY_NAT
        else:
            ts = cast_from_unit(ts, unit)
            obj.value = ts
            pandas_datetime_to_datetimestruct(ts, PANDAS_FR_ns, &obj.dts)
    elif PyDateTime_Check(ts):
        return convert_datetime_to_tsobject(ts, tz)
    elif PyDate_Check(ts):
        # Keep the converter same as PyDateTime's
        ts = datetime.combine(ts, datetime_time())
        return convert_datetime_to_tsobject(ts, tz)
    elif getattr(ts, '_typ', None) == 'period':
        raise ValueError(
            "Cannot convert Period to Timestamp "
            "unambiguously. Use to_timestamp")
    else:
        raise TypeError('Cannot convert input [{}] of type {} to '
                        'Timestamp'.format(ts, type(ts)))

    if obj.value != NPY_NAT:
        _check_dts_bounds(&obj.dts)

    if tz is not None:
        _localize_tso(obj, tz)

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

        # sort of a temporary hack
        if ts.tzinfo is not None:
            if (hasattr(tz, 'normalize') and
                hasattr(ts.tzinfo, '_utcoffset')):
                ts = tz.normalize(ts)
                obj.value = _pydatetime_to_dts(ts, &obj.dts)
                obj.tzinfo = ts.tzinfo
            else:
                # tzoffset
                try:
                    tz = ts.astimezone(tz).tzinfo
                except:
                    pass
                obj.value = _pydatetime_to_dts(ts, &obj.dts)
                ts_offset = get_utcoffset(ts.tzinfo, ts)
                obj.value -= int(ts_offset.total_seconds() * 1e9)
                tz_offset = get_utcoffset(tz, ts)
                obj.value += int(tz_offset.total_seconds() * 1e9)
                pandas_datetime_to_datetimestruct(obj.value,
                                                  PANDAS_FR_ns, &obj.dts)
                obj.tzinfo = tz
        elif not is_utc(tz):
            ts = _localize_pydatetime(ts, tz)
            obj.value = _pydatetime_to_dts(ts, &obj.dts)
            obj.tzinfo = ts.tzinfo
        else:
            # UTC
            obj.value = _pydatetime_to_dts(ts, &obj.dts)
            obj.tzinfo = pytz.utc
    else:
        obj.value = _pydatetime_to_dts(ts, &obj.dts)
        obj.tzinfo = ts.tzinfo

    if obj.tzinfo is not None and not is_utc(obj.tzinfo):
        offset = get_utcoffset(obj.tzinfo, ts)
        obj.value -= int(offset.total_seconds() * 1e9)

    if is_timestamp(ts):
        obj.value += ts.nanosecond
        obj.dts.ps = ts.nanosecond * 1000

    if nanos:
        obj.value += nanos
        obj.dts.ps = nanos * 1000

    _check_dts_bounds(&obj.dts)
    return obj


cpdef convert_str_to_tsobject(object ts, object tz, object unit,
                              dayfirst=False, yearfirst=False):
    """ ts must be a string """

    cdef:
        _TSObject obj
        int out_local = 0, out_tzoffset = 0

    if tz is not None:
        tz = maybe_get_tz(tz)

    obj = _TSObject()

    assert util.is_string_object(ts)

    if len(ts) == 0 or ts in _nat_strings:
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
            obj.value = pandas_datetimestruct_to_datetime(
                PANDAS_FR_ns, &obj.dts)
            _check_dts_bounds(&obj.dts)
            if out_local == 1:
                obj.tzinfo = pytz.FixedOffset(out_tzoffset)
                obj.value = tz_convert_single(obj.value, obj.tzinfo, 'UTC')
                if tz is None:
                    _check_dts_bounds(&obj.dts)
                    return obj
                else:
                    # Keep the converter same as PyDateTime's
                    obj = convert_to_tsobject(obj.value, obj.tzinfo,
                                              None, 0, 0)
                    dtime = datetime(obj.dts.year, obj.dts.month, obj.dts.day,
                                     obj.dts.hour, obj.dts.min, obj.dts.sec,
                                     obj.dts.us, obj.tzinfo)
                    obj = convert_datetime_to_tsobject(dtime, tz,
                                                       nanos=obj.dts.ps / 1000)
                    return obj

            else:
                ts = obj.value
                if tz is not None:
                    # shift for _localize_tso
                    ts = tz_localize_to_utc(np.array([ts], dtype='i8'), tz,
                                            ambiguous='raise',
                                            errors='raise')[0]
        except ValueError:
            try:
                ts = parse_datetime_string(
                    ts, dayfirst=dayfirst, yearfirst=yearfirst)
            except Exception:
                raise ValueError("could not convert string to Timestamp")

    return convert_to_tsobject(ts, tz, unit, dayfirst, yearfirst)


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
    obj.value = pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &obj.dts)
    _check_dts_bounds(&obj.dts)
    if out_local == 1:
        obj.tzinfo = pytz.FixedOffset(out_tzoffset)
        obj.value = tz_convert_single(obj.value, obj.tzinfo, 'UTC')
        return Timestamp(obj.value, tz=obj.tzinfo)
    else:
        return Timestamp(obj.value)

cdef inline void _localize_tso(_TSObject obj, object tz):
    """
    Take a TSObject in UTC and localizes to timezone tz.
    """
    cdef:
        ndarray[int64_t] trans, deltas
        Py_ssize_t delta, posn

    if is_utc(tz):
        obj.tzinfo = tz
    elif is_tzlocal(tz):
        pandas_datetime_to_datetimestruct(obj.value, PANDAS_FR_ns, &obj.dts)
        dt = datetime(obj.dts.year, obj.dts.month, obj.dts.day, obj.dts.hour,
                      obj.dts.min, obj.dts.sec, obj.dts.us, tz)
        delta = int(get_utcoffset(tz, dt).total_seconds()) * 1000000000
        if obj.value != NPY_NAT:
            pandas_datetime_to_datetimestruct(obj.value + delta,
                                              PANDAS_FR_ns, &obj.dts)
        else:
            pandas_datetime_to_datetimestruct(obj.value,
                                              PANDAS_FR_ns, &obj.dts)
        obj.tzinfo = tz
    else:
        # Adjust datetime64 timestamp, recompute datetimestruct
        trans, deltas, typ = get_dst_info(tz)

        pos = trans.searchsorted(obj.value, side='right') - 1

        # static/pytz/dateutil specific code
        if is_fixed_offset(tz):
            # statictzinfo
            if len(deltas) > 0 and obj.value != NPY_NAT:
                pandas_datetime_to_datetimestruct(obj.value + deltas[0],
                                                  PANDAS_FR_ns, &obj.dts)
            else:
                pandas_datetime_to_datetimestruct(
                    obj.value, PANDAS_FR_ns, &obj.dts)
            obj.tzinfo = tz
        elif treat_tz_as_pytz(tz):
            inf = tz._transition_info[pos]
            if obj.value != NPY_NAT:
                pandas_datetime_to_datetimestruct(obj.value + deltas[pos],
                                                  PANDAS_FR_ns, &obj.dts)
            else:
                pandas_datetime_to_datetimestruct(obj.value,
                                                  PANDAS_FR_ns, &obj.dts)
            obj.tzinfo = tz._tzinfos[inf]
        elif treat_tz_as_dateutil(tz):
            if obj.value != NPY_NAT:
                pandas_datetime_to_datetimestruct(obj.value + deltas[pos],
                                                  PANDAS_FR_ns, &obj.dts)
            else:
                pandas_datetime_to_datetimestruct(obj.value,
                                                  PANDAS_FR_ns, &obj.dts)
            obj.tzinfo = tz
        else:
            obj.tzinfo = tz


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


class OutOfBoundsDatetime(ValueError):
    pass

cdef inline _check_dts_bounds(pandas_datetimestruct *dts):
    if check_dts_bounds(dts):
        fmt = '%d-%.2d-%.2d %.2d:%.2d:%.2d' % (dts.year, dts.month,
                                               dts.day, dts.hour,
                                               dts.min, dts.sec)
        raise OutOfBoundsDatetime(
            'Out of bounds nanosecond timestamp: %s' % fmt)


def datetime_to_datetime64(ndarray[object] values):
    cdef:
        Py_ssize_t i, n = len(values)
        object val, inferred_tz = None
        ndarray[int64_t] iresult
        pandas_datetimestruct dts
        _TSObject _ts

    result = np.empty(n, dtype='M8[ns]')
    iresult = result.view('i8')
    for i in range(n):
        val = values[i]
        if _checknull_with_nat(val):
            iresult[i] = NPY_NAT
        elif PyDateTime_Check(val):
            if val.tzinfo is not None:
                if inferred_tz is not None:
                    if get_timezone(val.tzinfo) != inferred_tz:
                        raise ValueError('Array must be all same time zone')
                else:
                    inferred_tz = get_timezone(val.tzinfo)

                _ts = convert_datetime_to_tsobject(val, None)
                iresult[i] = _ts.value
                _check_dts_bounds(&_ts.dts)
            else:
                if inferred_tz is not None:
                    raise ValueError(
                        'Cannot mix tz-aware with tz-naive values')
                iresult[i] = _pydatetime_to_dts(val, &dts)
                _check_dts_bounds(&dts)
        else:
            raise TypeError('Unrecognized value type: %s' % type(val))

    return result, inferred_tz


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
        show_ns = (consider_values%1000).any()

        if not show_ns:
            consider_values //= 1000
            show_us = (consider_values%1000).any()

            if not show_ms:
                consider_values //= 1000
                show_ms = (consider_values%1000).any()

    for i in range(N):
        val = values[i]

        if val == NPY_NAT:
            result[i] = na_rep
        elif basic_format:

            pandas_datetime_to_datetimestruct(val, PANDAS_FR_ns, &dts)
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

_DEFAULT_DATETIME = datetime(1, 1, 1).replace(
    hour=0, minute=0, second=0, microsecond=0)
_MONTHS = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL',
           'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
_MONTH_NUMBERS = {k: i for i, k in enumerate(_MONTHS)}
_MONTH_ALIASES = {(k + 1): v for k, v in enumerate(_MONTHS)}


cpdef object _get_rule_month(object source, object default='DEC'):
    """
    Return starting month of given freq, default is December.

    Example
    -------
    >>> _get_rule_month('D')
    'DEC'

    >>> _get_rule_month('A-JAN')
    'JAN'
    """
    if hasattr(source, 'freqstr'):
        source = source.freqstr
    source = source.upper()
    if '-' not in source:
        return default
    else:
        return source.split('-')[1]


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

            if _checknull_with_nat(val):
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

            elif util.is_string_object(val):
                if len(val) == 0 or val in _nat_strings:
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
                    raise ValueError("non convertible value {0}"
                                     "with the unit '{1}'".format(
                                         val,
                                         unit))
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

        if _checknull_with_nat(val):
            oresult[i] = NaT
        elif is_integer_object(val) or is_float_object(val):

            if val != val or val == NPY_NAT:
                oresult[i] = NaT
            else:
                try:
                    oresult[i] = Timestamp(cast_from_unit(val, unit))
                except:
                    oresult[i] = val

        elif util.is_string_object(val):
            if len(val) == 0 or val in _nat_strings:
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

            if _checknull_with_nat(val):
                iresult[i] = NPY_NAT

            elif PyDateTime_Check(val):
                seen_datetime=1
                if val.tzinfo is not None:
                    if utc_convert:
                        _ts = convert_datetime_to_tsobject(val, None)
                        iresult[i] = _ts.value
                        try:
                            _check_dts_bounds(&_ts.dts)
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
                    iresult[i] = _pydatetime_to_dts(val, &dts)
                    if is_timestamp(val):
                        iresult[i] += (<_Timestamp>val).nanosecond
                    try:
                        _check_dts_bounds(&dts)
                    except ValueError:
                        if is_coerce:
                            iresult[i] = NPY_NAT
                            continue
                        raise

            elif PyDate_Check(val):
                iresult[i] = _date_to_datetime64(val, &dts)
                try:
                    _check_dts_bounds(&dts)
                    seen_datetime=1
                except ValueError:
                    if is_coerce:
                        iresult[i] = NPY_NAT
                        continue
                    raise

            elif util.is_datetime64_object(val):
                if get_datetime64_value(val) == NPY_NAT:
                    iresult[i] = NPY_NAT
                else:
                    try:
                        iresult[i] = _get_datetime64_nanos(val)
                        seen_datetime=1
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
                    seen_integer=1
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

            elif util.is_string_object(val):
                # string

                try:
                    if len(val) == 0 or val in _nat_strings:
                        iresult[i] = NPY_NAT
                        continue

                    seen_string=1
                    _string_to_dts(val, &dts, &out_local, &out_tzoffset)
                    value = pandas_datetimestruct_to_datetime(
                        PANDAS_FR_ns, &dts)
                    if out_local == 1:
                        tz = pytz.FixedOffset(out_tzoffset)
                        value = tz_convert_single(value, tz, 'UTC')
                    iresult[i] = value
                    _check_dts_bounds(&dts)
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

        if  seen_datetime and seen_integer:
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
            if _checknull_with_nat(val):
                if PyFloat_Check(val):
                    oresult[i] = np.nan
                else:
                    oresult[i] = NaT
            elif util.is_datetime64_object(val):
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
            if _checknull_with_nat(val):
                oresult[i] = val
            elif util.is_string_object(val):

                if len(val) == 0 or val in _nat_strings:
                    oresult[i] = 'NaT'
                    continue

                try:
                    oresult[i] = parse_datetime_string(val, dayfirst=dayfirst,
                                                       yearfirst=yearfirst)
                    _pydatetime_to_dts(oresult[i], &dts)
                    _check_dts_bounds(&dts)
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


# Similar to Timestamp/datetime, this is a construction requirement for
# timedeltas that we need to do object instantiation in python. This will
# serve as a C extension type that shadows the Python class, where we do any
# heavy lifting.
cdef class _Timedelta(timedelta):

    cdef readonly:
        int64_t value     # nanoseconds
        object freq       # frequency reference
        bint is_populated # are my components populated
        int64_t _sign, _d, _h, _m, _s, _ms, _us, _ns

    def __hash__(_Timedelta self):
        if self._has_ns():
            return hash(self.value)
        else:
            return timedelta.__hash__(self)

    def __richcmp__(_Timedelta self, object other, int op):
        cdef:
            _Timedelta ots
            int ndim

        if isinstance(other, _Timedelta):
            ots = other
        elif isinstance(other, timedelta):
            ots = Timedelta(other)
        else:
            ndim = getattr(other, _NDIM_STRING, -1)

            if ndim != -1:
                if ndim == 0:
                    if isinstance(other, np.timedelta64):
                        other = Timedelta(other)
                    else:
                        if op == Py_EQ:
                            return False
                        elif op == Py_NE:
                            return True

                        # only allow ==, != ops
                        raise TypeError('Cannot compare type %r with type %r' %
                                        (type(self).__name__,
                                         type(other).__name__))
                if isinstance(other, np.ndarray):
                    return PyObject_RichCompare(np.array([self]), other, op)
                return PyObject_RichCompare(other, self, _reverse_ops[op])
            else:
                if op == Py_EQ:
                    return False
                elif op == Py_NE:
                    return True
                raise TypeError('Cannot compare type %r with type %r' %
                                (type(self).__name__, type(other).__name__))

        return _cmp_scalar(self.value, ots.value, op)

    def _ensure_components(_Timedelta self):
        """
        compute the components
        """
        cdef int64_t sfrac, ifrac, frac, ivalue = self.value

        if self.is_populated:
            return

        # put frac in seconds
        frac = ivalue /(1000 *1000 *1000)
        if frac < 0:
            self._sign = -1

            # even fraction
            if (-frac % 86400) != 0:
                self._d = -frac /86400 + 1
                frac += 86400 *self._d
            else:
                frac = -frac
        else:
            self._sign = 1
            self._d = 0

        if frac >= 86400:
            self._d += frac / 86400
            frac -= self._d * 86400

        if frac >= 3600:
            self._h = frac / 3600
            frac -= self._h * 3600
        else:
            self._h = 0

        if frac >= 60:
            self._m = frac / 60
            frac -= self._m * 60
        else:
            self._m = 0

        if frac >= 0:
            self._s = frac
            frac -= self._s
        else:
            self._s = 0

        sfrac = (self._h * 3600 + self._m * 60
                 + self._s) * (1000 * 1000 * 1000)
        if self._sign < 0:
            ifrac = ivalue + self._d *DAY_NS - sfrac
        else:
            ifrac = ivalue - (self._d *DAY_NS + sfrac)

        if ifrac != 0:
            self._ms = ifrac /(1000 *1000)
            ifrac -= self._ms *1000 *1000
            self._us = ifrac /1000
            ifrac -= self._us *1000
            self._ns = ifrac
        else:
            self._ms = 0
            self._us = 0
            self._ns = 0

        self.is_populated = 1

    cpdef timedelta to_pytimedelta(_Timedelta self):
        """
        return an actual datetime.timedelta object
        note: we lose nanosecond resolution if any
        """
        return timedelta(microseconds=int(self.value) /1000)

    cpdef bint _has_ns(self):
        return self.value % 1000 != 0

# components named tuple
Components = collections.namedtuple('Components', [
    'days', 'hours', 'minutes', 'seconds',
    'milliseconds', 'microseconds', 'nanoseconds'])

# Python front end to C extension type _Timedelta
# This serves as the box for timedelta64


class Timedelta(_Timedelta):
    """
    Represents a duration, the difference between two dates or times.

    Timedelta is the pandas equivalent of python's ``datetime.timedelta``
    and is interchangable with it in most cases.

    Parameters
    ----------
    value : Timedelta, timedelta, np.timedelta64, string, or integer
    unit : string, [D,h,m,s,ms,us,ns]
        Denote the unit of the input, if input is an integer. Default 'ns'.
    days, seconds, microseconds,
    milliseconds, minutes, hours, weeks : numeric, optional
        Values for construction in compat with datetime.timedelta.
        np ints and floats will be coereced to python ints and floats.

    Notes
    -----
    The ``.value`` attribute is always in ns.

    """

    def __new__(cls, object value=_no_input, unit=None, **kwargs):
        cdef _Timedelta td_base

        if value is _no_input:
            if not len(kwargs):
                raise ValueError(
                    "cannot construct a Timedelta without a value/unit or "
                    "descriptive keywords (days,seconds....)")

            def _to_py_int_float(v):
                if is_integer_object(v):
                    return int(v)
                elif is_float_object(v):
                    return float(v)
                raise TypeError(
                    "Invalid type {0}. Must be int or float.".format(type(v)))

            kwargs = dict([(k, _to_py_int_float(v))
                            for k, v in iteritems(kwargs)])

            try:
                nano = kwargs.pop('nanoseconds', 0)
                value = convert_to_timedelta64(
                    timedelta(**kwargs), 'ns') + nano
            except TypeError as e:
                raise ValueError("cannot construct a Timedelta from the "
                                 "passed arguments, allowed keywords are "
                                 "[weeks, days, hours, minutes, seconds, "
                                 "milliseconds, microseconds, nanoseconds]")

        if isinstance(value, Timedelta):
            value = value.value
        elif util.is_string_object(value):
            value = np.timedelta64(parse_timedelta_string(value))
        elif isinstance(value, timedelta):
            value = convert_to_timedelta64(value, 'ns')
        elif isinstance(value, np.timedelta64):
            if unit is not None:
                value = value.astype('timedelta64[{0}]'.format(unit))
            value = value.astype('timedelta64[ns]')
        elif hasattr(value, 'delta'):
            value = np.timedelta64(_delta_to_nanoseconds(value.delta), 'ns')
        elif is_integer_object(value) or util.is_float_object(value):
            # unit=None is de-facto 'ns'
            value = convert_to_timedelta64(value, unit)
        elif _checknull_with_nat(value):
            return NaT
        else:
            raise ValueError(
                "Value must be Timedelta, string, integer, "
                "float, timedelta or convertible")

        if isinstance(value, np.timedelta64):
            value = value.view('i8')

        # nat
        if value == NPY_NAT:
            return NaT

        # make timedelta happy
        td_base = _Timedelta.__new__(cls, microseconds=int(value) /1000)
        td_base.value = value
        td_base.is_populated = 0
        return td_base

    @property
    def delta(self):
        """ return out delta in ns (for internal compat) """
        return self.value

    @property
    def asm8(self):
        """ return a numpy timedelta64 array view of myself """
        return np.int64(self.value).view('m8[ns]')

    @property
    def resolution(self):
        """ return a string representing the lowest resolution that we have """

        self._ensure_components()
        if self._ns:
            return "N"
        elif self._us:
            return "U"
        elif self._ms:
            return "L"
        elif self._s:
            return "S"
        elif self._m:
            return "T"
        elif self._h:
            return "H"
        else:
            return "D"

    def _round(self, freq, rounder):

        cdef int64_t result, unit

        from pandas.tseries.frequencies import to_offset
        unit = to_offset(freq).nanos
        result = unit *rounder(self.value /float(unit))
        return Timedelta(result, unit='ns')

    def round(self, freq):
        """
        Round the Timedelta to the specified resolution

        Returns
        -------
        a new Timedelta rounded to the given resolution of `freq`

        Parameters
        ----------
        freq : a freq string indicating the rounding resolution

        Raises
        ------
        ValueError if the freq cannot be converted
        """
        return self._round(freq, np.round)

    def floor(self, freq):
        """
        return a new Timedelta floored to this resolution

        Parameters
        ----------
        freq : a freq string indicating the flooring resolution
        """
        return self._round(freq, np.floor)

    def ceil(self, freq):
        """
        return a new Timedelta ceiled to this resolution

        Parameters
        ----------
        freq : a freq string indicating the ceiling resolution
        """
        return self._round(freq, np.ceil)

    def _repr_base(self, format=None):
        """

        Parameters
        ----------
        format : None|all|even_day|sub_day|long

        Returns
        -------
        converted : string of a Timedelta

        """
        cdef object sign_pretty, sign2_pretty, seconds_pretty, subs

        self._ensure_components()

        if self._sign < 0:
            sign_pretty = "-"
            sign2_pretty = " +"
        else:
            sign_pretty = ""
            sign2_pretty = " "

        # show everything
        if format == 'all':
            seconds_pretty = "%02d.%03d%03d%03d" % (
                self._s, self._ms, self._us, self._ns)
            return "%s%d days%s%02d:%02d:%s" % (sign_pretty, self._d,
                                                sign2_pretty, self._h,
                                                self._m, seconds_pretty)

        # by default not showing nano
        if self._ms or self._us or self._ns:
            seconds_pretty = "%02d.%03d%03d" % (self._s, self._ms, self._us)
        else:
            seconds_pretty = "%02d" % self._s

        # if we have a partial day
        subs = (self._h or self._m or self._s or
                self._ms or self._us or self._ns)

        if format == 'even_day':
            if not subs:
                return "%s%d days" % (sign_pretty, self._d)

        elif format == 'sub_day':
            if not self._d:

                # degenerate, don't need the extra space
                if self._sign > 0:
                    sign2_pretty = ""
                return "%s%s%02d:%02d:%s" % (sign_pretty, sign2_pretty,
                                             self._h, self._m, seconds_pretty)

        if subs or format=='long':
            return "%s%d days%s%02d:%02d:%s" % (sign_pretty, self._d,
                                                sign2_pretty, self._h,
                                                self._m, seconds_pretty)
        return "%s%d days" % (sign_pretty, self._d)

    def __repr__(self):
        return "Timedelta('{0}')".format(self._repr_base(format='long'))
    def __str__(self):
        return self._repr_base(format='long')

    @property
    def components(self):
        """ Return a Components NamedTuple-like """
        self._ensure_components()
        if self._sign < 0:
            return Components(-self._d, self._h, self._m, self._s,
                              self._ms, self._us, self._ns)

        # return the named tuple
        return Components(self._d, self._h, self._m, self._s,
                          self._ms, self._us, self._ns)

    @property
    def days(self):
        """
        Number of Days

        .components will return the shown components
        """
        self._ensure_components()
        if self._sign < 0:
            return -1 *self._d
        return self._d

    @property
    def seconds(self):
        """
        Number of seconds (>= 0 and less than 1 day).

        .components will return the shown components
        """
        self._ensure_components()
        return self._h *3600 + self._m *60 + self._s

    @property
    def microseconds(self):
        """
        Number of microseconds (>= 0 and less than 1 second).

        .components will return the shown components
        """
        self._ensure_components()
        return self._ms *1000 + self._us

    @property
    def nanoseconds(self):
        """
        Number of nanoseconds (>= 0 and less than 1 microsecond).

        .components will return the shown components
        """
        self._ensure_components()
        return self._ns

    def total_seconds(self):
        """
        Total duration of timedelta in seconds (to ns precision)
        """
        return 1e-9 *self.value

    def isoformat(self):
        """
        Format Timedelta as ISO 8601 Duration like
        `P[n]Y[n]M[n]DT[n]H[n]M[n]S`, where the `[n]`s are replaced by the
        values. See https://en.wikipedia.org/wiki/ISO_8601#Durations

        .. versionadded:: 0.20.0

        Returns
        -------
        formatted : str

        Notes
        -----
        The longest component is days, whose value may be larger than
        365.
        Every component is always included, even if its value is 0.
        Pandas uses nanosecond precision, so up to 9 decimal places may
        be included in the seconds component.
        Trailing 0's are removed from the seconds component after the decimal.
        We do not 0 pad components, so it's `...T5H...`, not `...T05H...`

        Examples
        --------
        >>> td = pd.Timedelta(days=6, minutes=50, seconds=3,
        ...                   milliseconds=10, microseconds=10, nanoseconds=12)
        >>> td.isoformat()
        'P6DT0H50M3.010010012S'
        >>> pd.Timedelta(hours=1, seconds=10).isoformat()
        'P0DT0H0M10S'
        >>> pd.Timedelta(hours=1, seconds=10).isoformat()
        'P0DT0H0M10S'
        >>> pd.Timedelta(days=500.5).isoformat()
        'P500DT12H0MS'

        See Also
        --------
        Timestamp.isoformat
        """
        components = self.components
        seconds = '{}.{:0>3}{:0>3}{:0>3}'.format(components.seconds,
                                                 components.milliseconds,
                                                 components.microseconds,
                                                 components.nanoseconds)
        # Trim unnecessary 0s, 1.000000000 -> 1
        seconds = seconds.rstrip('0').rstrip('.')
        tpl = 'P{td.days}DT{td.hours}H{td.minutes}M{seconds}S'.format(
            td=components, seconds=seconds)
        return tpl

    def __setstate__(self, state):
        (value) = state
        self.value = value

    def __reduce__(self):
        object_state = self.value,
        return (Timedelta, object_state)

    def view(self, dtype):
        """ array view compat """
        return np.timedelta64(self.value).view(dtype)

    def to_timedelta64(self):
        """ Returns a numpy.timedelta64 object with 'ns' precision """
        return np.timedelta64(self.value, 'ns')

    def _validate_ops_compat(self, other):

        # return True if we are compat with operating
        if _checknull_with_nat(other):
            return True
        elif isinstance(other, (Timedelta, timedelta, np.timedelta64)):
            return True
        elif util.is_string_object(other):
            return True
        elif hasattr(other, 'delta'):
            return True
        return False

    # higher than np.ndarray and np.matrix
    __array_priority__ = 100

    def _binary_op_method_timedeltalike(op, name):
        # define a binary operation that only works if the other argument is
        # timedelta like or an array of timedeltalike
        def f(self, other):
            # an offset
            if hasattr(other, 'delta') and not isinstance(other, Timedelta):
                return op(self, other.delta)

            # a datetimelike
            if (isinstance(other, (datetime, np.datetime64))
                    and not isinstance(other, (Timestamp, NaTType))):
                return op(self, Timestamp(other))

            # nd-array like
            if hasattr(other, 'dtype'):
                if other.dtype.kind not in ['m', 'M']:
                    # raise rathering than letting numpy return wrong answer
                    return NotImplemented
                return op(self.to_timedelta64(), other)

            if not self._validate_ops_compat(other):
                return NotImplemented

            if other is NaT:
                return NaT

            try:
                other = Timedelta(other)
            except ValueError:
                # failed to parse as timedelta
                return NotImplemented

            return Timedelta(op(self.value, other.value), unit='ns')

        f.__name__ = name
        return f

    __add__ = _binary_op_method_timedeltalike(lambda x, y: x + y, '__add__')
    __radd__ = _binary_op_method_timedeltalike(lambda x, y: x + y, '__radd__')
    __sub__ = _binary_op_method_timedeltalike(lambda x, y: x - y, '__sub__')
    __rsub__ = _binary_op_method_timedeltalike(lambda x, y: y - x, '__rsub__')

    def __mul__(self, other):

        # nd-array like
        if hasattr(other, 'dtype'):
            return other * self.to_timedelta64()

        if other is NaT:
            return NaT

        # only integers and floats allowed
        if not (is_integer_object(other) or is_float_object(other)):
            return NotImplemented

        return Timedelta(other * self.value, unit='ns')

    __rmul__ = __mul__

    def __truediv__(self, other):

        if hasattr(other, 'dtype'):
            return self.to_timedelta64() / other

        # integers or floats
        if is_integer_object(other) or is_float_object(other):
            return Timedelta(self.value /other, unit='ns')

        if not self._validate_ops_compat(other):
            return NotImplemented

        other = Timedelta(other)
        if other is NaT:
            return np.nan
        return self.value /float(other.value)

    def __rtruediv__(self, other):
        if hasattr(other, 'dtype'):
            return other / self.to_timedelta64()

        if not self._validate_ops_compat(other):
            return NotImplemented

        other = Timedelta(other)
        if other is NaT:
            return NaT
        return float(other.value) / self.value

    if not PY3:
        __div__ = __truediv__
        __rdiv__ = __rtruediv__

    def __floordiv__(self, other):

        if hasattr(other, 'dtype'):

            # work with i8
            other = other.astype('m8[ns]').astype('i8')

            return self.value // other

        # integers only
        if is_integer_object(other):
            return Timedelta(self.value // other, unit='ns')

        if not self._validate_ops_compat(other):
            return NotImplemented

        other = Timedelta(other)
        if other is NaT:
            return np.nan
        return self.value // other.value

    def __rfloordiv__(self, other):
        if hasattr(other, 'dtype'):

            # work with i8
            other = other.astype('m8[ns]').astype('i8')
            return other // self.value

        if not self._validate_ops_compat(other):
            return NotImplemented

        other = Timedelta(other)
        if other is NaT:
            return NaT
        return other.value // self.value

    def _op_unary_method(func, name):

        def f(self):
            return Timedelta(func(self.value), unit='ns')
        f.__name__ = name
        return f

    __inv__ = _op_unary_method(lambda x: -x, '__inv__')
    __neg__ = _op_unary_method(lambda x: -x, '__neg__')
    __pos__ = _op_unary_method(lambda x: x, '__pos__')
    __abs__ = _op_unary_method(lambda x: abs(x), '__abs__')

# resolution in ns
Timedelta.min = Timedelta(np.iinfo(np.int64).min +1)
Timedelta.max = Timedelta(np.iinfo(np.int64).max)

cdef PyTypeObject* td_type = <PyTypeObject*> Timedelta


cdef inline bint is_timedelta(object o):
    return Py_TYPE(o) == td_type # isinstance(o, Timedelta)


cpdef array_to_timedelta64(ndarray[object] values, unit='ns', errors='raise'):
    """
    Convert an ndarray to an array of timedeltas. If errors == 'coerce',
    coerce non-convertible objects to NaT. Otherwise, raise.
    """

    cdef:
        Py_ssize_t i, n
        ndarray[int64_t] iresult

    if errors not in ('ignore', 'raise', 'coerce'):
        raise ValueError("errors must be one of 'ignore', "
                         "'raise', or 'coerce'}")

    n = values.shape[0]
    result = np.empty(n, dtype='m8[ns]')
    iresult = result.view('i8')

    # Usually, we have all strings. If so, we hit the fast path.
    # If this path fails, we try conversion a different way, and
    # this is where all of the error handling will take place.
    try:
        for i in range(n):
            result[i] = parse_timedelta_string(values[i])
    except:
        for i in range(n):
            try:
                result[i] = convert_to_timedelta64(values[i], unit)
            except ValueError:
                if errors == 'coerce':
                    result[i] = NPY_NAT
                else:
                    raise

    return iresult

cdef dict timedelta_abbrevs = { 'D': 'd',
                                'd': 'd',
                                'days': 'd',
                                'day': 'd',
                                'hours': 'h',
                                'hour': 'h',
                                'hr': 'h',
                                'h': 'h',
                                'm': 'm',
                                'minute': 'm',
                                'min': 'm',
                                'minutes': 'm',
                                's': 's',
                                'seconds': 's',
                                'sec': 's',
                                'second': 's',
                                'ms': 'ms',
                                'milliseconds': 'ms',
                                'millisecond': 'ms',
                                'milli': 'ms',
                                'millis': 'ms',
                                'us': 'us',
                                'microseconds': 'us',
                                'microsecond': 'us',
                                'micro': 'us',
                                'micros': 'us',
                                'ns': 'ns',
                                'nanoseconds': 'ns',
                                'nano': 'ns',
                                'nanos': 'ns',
                                'nanosecond': 'ns',
                                }
timedelta_abbrevs_map = timedelta_abbrevs

cdef inline int64_t timedelta_as_neg(int64_t value, bint neg):
    """

    Parameters
    ----------
    value : int64_t of the timedelta value
    neg : boolean if the a negative value
    """
    if neg:
        return -value
    return value

cdef inline timedelta_from_spec(object number, object frac, object unit):
    """

    Parameters
    ----------
    number : a list of number digits
    frac : a list of frac digits
    unit : a list of unit characters
    """
    cdef object n

    try:
        unit = ''.join(unit)
        unit = timedelta_abbrevs[unit.lower()]
    except KeyError:
        raise ValueError("invalid abbreviation: {0}".format(unit))

    n = ''.join(number) + '.' + ''.join(frac)
    return cast_from_unit(float(n), unit)

cdef inline parse_timedelta_string(object ts):
    """
    Parse a regular format timedelta string. Return an int64_t (in ns)
    or raise a ValueError on an invalid parse.
    """

    cdef:
        unicode c
        bint neg=0, have_dot=0, have_value=0, have_hhmmss=0
        object current_unit=None
        int64_t result=0, m=0, r
        list number=[], frac=[], unit=[]

    # neg : tracks if we have a leading negative for the value
    # have_dot : tracks if we are processing a dot (either post hhmmss or
    #            inside an expression)
    # have_value : track if we have at least 1 leading unit
    # have_hhmmss : tracks if we have a regular format hh:mm:ss

    if len(ts) == 0 or ts in _nat_strings:
        return NPY_NAT

    # decode ts if necessary
    if not PyUnicode_Check(ts) and not PY3:
        ts = str(ts).decode('utf-8')

    for c in ts:

        # skip whitespace / commas
        if c == ' ' or c == ',':
            pass

        # positive signs are ignored
        elif c == '+':
            pass

        # neg
        elif c == '-':

            if neg or have_value or have_hhmmss:
                raise ValueError("only leading negative signs are allowed")

            neg = 1

        # number (ascii codes)
        elif ord(c) >= 48 and ord(c) <= 57:

            if have_dot:

                # we found a dot, but now its just a fraction
                if len(unit):
                    number.append(c)
                    have_dot = 0
                else:
                    frac.append(c)

            elif not len(unit):
                number.append(c)

            else:
                r = timedelta_from_spec(number, frac, unit)
                unit, number, frac = [], [c], []

                result += timedelta_as_neg(r, neg)

        # hh:mm:ss.
        elif c == ':':

            # we flip this off if we have a leading value
            if have_value:
                neg = 0

            # we are in the pattern hh:mm:ss pattern
            if len(number):
                if current_unit is None:
                    current_unit = 'h'
                    m = 1000000000L * 3600
                elif current_unit == 'h':
                    current_unit = 'm'
                    m = 1000000000L * 60
                elif current_unit == 'm':
                    current_unit = 's'
                    m = 1000000000L
                r = <int64_t> int(''.join(number)) * m
                result += timedelta_as_neg(r, neg)
                have_hhmmss = 1
            else:
                raise ValueError("expecting hh:mm:ss format, "
                                 "received: {0}".format(ts))

            unit, number = [], []

        # after the decimal point
        elif c == '.':

            if len(number) and current_unit is not None:

                # by definition we had something like
                # so we need to evaluate the final field from a
                # hh:mm:ss (so current_unit is 'm')
                if current_unit != 'm':
                    raise ValueError("expected hh:mm:ss format before .")
                m = 1000000000L
                r = <int64_t> int(''.join(number)) * m
                result += timedelta_as_neg(r, neg)
                have_value = 1
                unit, number, frac = [], [], []

            have_dot = 1

        # unit
        else:
            unit.append(c)
            have_value = 1
            have_dot = 0

    # we had a dot, but we have a fractional
    # value since we have an unit
    if have_dot and len(unit):
        r = timedelta_from_spec(number, frac, unit)
        result += timedelta_as_neg(r, neg)

    # we have a dot as part of a regular format
    # e.g. hh:mm:ss.fffffff
    elif have_dot:

        if ((len(number) or len(frac)) and not len(unit)
            and current_unit is None):
            raise ValueError("no units specified")

        if len(frac) > 0 and len(frac) <= 3:
            m = 10**(3 -len(frac)) * 1000L * 1000L
        elif len(frac) > 3 and len(frac) <= 6:
            m = 10**(6 -len(frac)) * 1000L
        else:
            m = 10**(9 -len(frac))

        r = <int64_t> int(''.join(frac)) * m
        result += timedelta_as_neg(r, neg)

    # we have a regular format
    # we must have seconds at this point (hence the unit is still 'm')
    elif current_unit is not None:
        if current_unit != 'm':
            raise ValueError("expected hh:mm:ss format")
        m = 1000000000L
        r = <int64_t> int(''.join(number)) * m
        result += timedelta_as_neg(r, neg)

    # we have a last abbreviation
    elif len(unit):
        if len(number):
            r = timedelta_from_spec(number, frac, unit)
            result += timedelta_as_neg(r, neg)
        else:
            raise ValueError("unit abbreviation w/o a number")

    # treat as nanoseconds
    # but only if we don't have anything else
    else:
        if have_value:
            raise ValueError("have leftover units")
        if len(number):
            r = timedelta_from_spec(number, frac, 'ns')
            result += timedelta_as_neg(r, neg)

    return result

cpdef convert_to_timedelta64(object ts, object unit):
    """
    Convert an incoming object to a timedelta64 if possible

    Handle these types of objects:
        - timedelta/Timedelta
        - timedelta64
        - an offset
        - np.int64 (with unit providing a possible modifier)
        - None/NaT

    Return an ns based int64

    # kludgy here until we have a timedelta scalar
    # handle the numpy < 1.7 case
    """
    if _checknull_with_nat(ts):
        return np.timedelta64(NPY_NAT)
    elif isinstance(ts, Timedelta):
        # already in the proper format
        ts = np.timedelta64(ts.value)
    elif util.is_datetime64_object(ts):
        # only accept a NaT here
        if ts.astype('int64') == NPY_NAT:
            return np.timedelta64(NPY_NAT)
    elif isinstance(ts, np.timedelta64):
        ts = ts.astype("m8[{0}]".format(unit.lower()))
    elif is_integer_object(ts):
        if ts == NPY_NAT:
            return np.timedelta64(NPY_NAT)
        else:
            if util.is_array(ts):
                ts = ts.astype('int64').item()
            if unit in ['Y', 'M', 'W']:
                ts = np.timedelta64(ts, unit)
            else:
                ts = cast_from_unit(ts, unit)
                ts = np.timedelta64(ts)
    elif is_float_object(ts):
        if util.is_array(ts):
            ts = ts.astype('int64').item()
        if unit in ['Y', 'M', 'W']:
            ts = np.timedelta64(int(ts), unit)
        else:
            ts = cast_from_unit(ts, unit)
            ts = np.timedelta64(ts)
    elif util.is_string_object(ts):
        ts = np.timedelta64(parse_timedelta_string(ts))
    elif hasattr(ts, 'delta'):
        ts = np.timedelta64(_delta_to_nanoseconds(ts), 'ns')

    if isinstance(ts, timedelta):
        ts = np.timedelta64(_delta_to_nanoseconds(ts), 'ns')
    elif not isinstance(ts, np.timedelta64):
        raise ValueError("Invalid type for timedelta "
                         "scalar: %s" % type(ts))
    return ts.astype('timedelta64[ns]')


#----------------------------------------------------------------------
# Conversion routines

cpdef int64_t _delta_to_nanoseconds(delta) except? -1:
    if isinstance(delta, np.ndarray):
        return delta.astype('m8[ns]').astype('int64')
    if hasattr(delta, 'nanos'):
        return delta.nanos
    if hasattr(delta, 'delta'):
        delta = delta.delta
    if is_timedelta64_object(delta):
        return delta.astype("timedelta64[ns]").item()
    if is_integer_object(delta):
        return delta

    return (delta.days * 24 * 60 * 60 * 1000000
            + delta.seconds * 1000000
            + delta.microseconds) * 1000


cdef inline _get_datetime64_nanos(object val):
    cdef:
        pandas_datetimestruct dts
        PANDAS_DATETIMEUNIT unit
        npy_datetime ival

    unit = get_datetime64_unit(val)
    ival = get_datetime64_value(val)

    if unit != PANDAS_FR_ns:
        pandas_datetime_to_datetimestruct(ival, unit, &dts)
        _check_dts_bounds(&dts)
        return pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &dts)
    else:
        return ival

cpdef inline int64_t cast_from_unit(object ts, object unit) except? -1:
    """ return a casting of the unit represented to nanoseconds
        round the fractional part of a float to our precision, p """
    cdef:
        int64_t m
        int p

    if unit == 'D' or unit == 'd':
        m = 1000000000L * 86400
        p = 6
    elif unit == 'h':
        m = 1000000000L * 3600
        p = 6
    elif unit == 'm':
        m = 1000000000L * 60
        p = 6
    elif unit == 's':
        m = 1000000000L
        p = 6
    elif unit == 'ms':
        m = 1000000L
        p = 3
    elif unit == 'us':
        m = 1000L
        p = 0
    elif unit == 'ns' or unit is None:
        m = 1L
        p = 0
    else:
        raise ValueError("cannot cast unit {0}".format(unit))

    # just give me the unit back
    if ts is None:
        return m

    # cast the unit, multiply base/frace separately
    # to avoid precision issues from float -> int
    base = <int64_t> ts
    frac = ts -base
    if p:
        frac = round(frac, p)
    return <int64_t> (base *m) + <int64_t> (frac *m)


def cast_to_nanoseconds(ndarray arr):
    cdef:
        Py_ssize_t i, n = arr.size
        ndarray[int64_t] ivalues, iresult
        PANDAS_DATETIMEUNIT unit
        pandas_datetimestruct dts

    shape = (<object> arr).shape

    ivalues = arr.view(np.int64).ravel()

    result = np.empty(shape, dtype='M8[ns]')
    iresult = result.ravel().view(np.int64)

    if len(iresult) == 0:
        return result

    unit = get_datetime64_unit(arr.flat[0])
    for i in range(n):
        if ivalues[i] != NPY_NAT:
            pandas_datetime_to_datetimestruct(ivalues[i], unit, &dts)
            iresult[i] = pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &dts)
            _check_dts_bounds(&dts)
        else:
            iresult[i] = NPY_NAT

    return result


cdef inline _to_i8(object val):
    cdef pandas_datetimestruct dts
    try:
        return val.value
    except AttributeError:
        if is_datetime64_object(val):
            return get_datetime64_value(val)
        elif PyDateTime_Check(val):
            return Timestamp(val).value
        return val

cpdef pydt_to_i8(object pydt):
    """
    Convert to int64 representation compatible with numpy datetime64; converts
    to UTC
    """
    cdef:
        _TSObject ts

    ts = convert_to_tsobject(pydt, None, None, 0, 0)

    return ts.value


def i8_to_pydt(int64_t i8, object tzinfo = None):
    """
    Inverse of pydt_to_i8
    """
    return Timestamp(i8)

#----------------------------------------------------------------------
# time zone conversion helpers

import pytz
UTC = pytz.utc


@cython.boundscheck(False)
@cython.wraparound(False)
def tz_convert(ndarray[int64_t] vals, object tz1, object tz2):
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
        ndarray[int64_t] utc_dates, tt, result, trans, deltas
        Py_ssize_t i, j, pos, n = len(vals)
        ndarray[Py_ssize_t] posn
        int64_t v, offset, delta
        pandas_datetimestruct dts

    if len(vals) == 0:
        return np.array([], dtype=np.int64)

    # Convert to UTC
    if get_timezone(tz1) != 'UTC':
        utc_dates = np.empty(n, dtype=np.int64)
        if is_tzlocal(tz1):
            for i in range(n):
                v = vals[i]
                if v == NPY_NAT:
                    utc_dates[i] = NPY_NAT
                else:
                    pandas_datetime_to_datetimestruct(v, PANDAS_FR_ns, &dts)
                    dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                                  dts.min, dts.sec, dts.us, tz1)
                    delta = (int(get_utcoffset(tz1, dt).total_seconds())
                             * 1000000000)
                    utc_dates[i] = v - delta
        else:
            trans, deltas, typ = get_dst_info(tz1)

            # all-NaT
            tt = vals[vals!=NPY_NAT]
            if not len(tt):
                return vals

            posn = trans.searchsorted(tt, side='right')
            j = 0
            for i in range(n):
                v = vals[i]
                if v == NPY_NAT:
                    utc_dates[i] = NPY_NAT
                else:
                    pos = posn[j] - 1
                    j = j + 1
                    if pos < 0:
                        raise ValueError('First time before start of DST info')
                    offset = deltas[pos]
                    utc_dates[i] = v - offset
    else:
        utc_dates = vals

    if get_timezone(tz2) == 'UTC':
        return utc_dates

    result = np.zeros(n, dtype=np.int64)
    if is_tzlocal(tz2):
        for i in range(n):
            v = utc_dates[i]
            if v == NPY_NAT:
                result[i] = NPY_NAT
            else:
                pandas_datetime_to_datetimestruct(v, PANDAS_FR_ns, &dts)
                dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                              dts.min, dts.sec, dts.us, tz2)
                delta = (int(get_utcoffset(tz2, dt).total_seconds())
                             * 1000000000)
                result[i] = v + delta
        return result

    # Convert UTC to other timezone
    trans, deltas, typ = get_dst_info(tz2)

    # use first non-NaT element
    # if all-NaT, return all-NaT
    if (result==NPY_NAT).all():
        return result

    # if all NaT, return all NaT
    tt = utc_dates[utc_dates!=NPY_NAT]
    if not len(tt):
        return utc_dates

    posn = trans.searchsorted(tt, side='right')

    j = 0
    for i in range(n):
        v = utc_dates[i]
        if vals[i] == NPY_NAT:
            result[i] = vals[i]
        else:
            pos = posn[j] - 1
            j = j + 1
            if pos < 0:
                raise ValueError('First time before start of DST info')
            offset = deltas[pos]
            result[i] = v + offset
    return result


cpdef int64_t tz_convert_single(int64_t val, object tz1, object tz2):
    """
    Convert the val (in i8) from timezone1 to timezone2

    This is a single timezone versoin of tz_convert

    Parameters
    ----------
    val : int64
    tz1 : string / timezone object
    tz2 : string / timezone object

    Returns
    -------
    int64 converted

    """

    cdef:
        ndarray[int64_t] trans, deltas
        Py_ssize_t pos
        int64_t v, offset, utc_date
        pandas_datetimestruct dts

    if val == NPY_NAT:
        return val

    # Convert to UTC
    if is_tzlocal(tz1):
        pandas_datetime_to_datetimestruct(val, PANDAS_FR_ns, &dts)
        dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                      dts.min, dts.sec, dts.us, tz1)
        delta = int(get_utcoffset(tz1, dt).total_seconds()) * 1000000000
        utc_date = val - delta
    elif get_timezone(tz1) != 'UTC':
        trans, deltas, typ = get_dst_info(tz1)
        pos = trans.searchsorted(val, side='right') - 1
        if pos < 0:
            raise ValueError('First time before start of DST info')
        offset = deltas[pos]
        utc_date = val - offset
    else:
        utc_date = val

    if get_timezone(tz2) == 'UTC':
        return utc_date
    if is_tzlocal(tz2):
        pandas_datetime_to_datetimestruct(val, PANDAS_FR_ns, &dts)
        dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                      dts.min, dts.sec, dts.us, tz2)
        delta = int(get_utcoffset(tz2, dt).total_seconds()) * 1000000000
        return utc_date + delta

    # Convert UTC to other timezone
    trans, deltas, typ = get_dst_info(tz2)

    pos = trans.searchsorted(utc_date, side='right') - 1
    if pos < 0:
        raise ValueError('First time before start of DST info')

    offset = deltas[pos]
    return utc_date + offset


@cython.boundscheck(False)
@cython.wraparound(False)
def tz_localize_to_utc(ndarray[int64_t] vals, object tz, object ambiguous=None,
                       object errors='raise'):
    """
    Localize tzinfo-naive i8 to given time zone (using pytz). If
    there are ambiguities in the values, raise AmbiguousTimeError.

    Returns
    -------
    localized : DatetimeIndex
    """
    cdef:
        ndarray[int64_t] trans, deltas, idx_shifted
        ndarray ambiguous_array
        Py_ssize_t i, idx, pos, ntrans, n = len(vals)
        int64_t *tdata
        int64_t v, left, right
        ndarray[int64_t] result, result_a, result_b, dst_hours
        pandas_datetimestruct dts
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
            pandas_datetime_to_datetimestruct(v, PANDAS_FR_ns, &dts)
            dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                          dts.min, dts.sec, dts.us, tz)
            delta = int(get_utcoffset(tz, dt).total_seconds()) * 1000000000
            result[i] = v - delta
        return result

    if util.is_string_object(ambiguous):
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
            raise ValueError(
                "Length of ambiguous bool-array must be the same size as vals")
        ambiguous_array = np.asarray(ambiguous)

    trans, deltas, typ = get_dst_info(tz)

    tdata = <int64_t*> trans.data
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
            stamp = Timestamp(vals[trans_idx])
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
                    stamp = Timestamp(vals[grp[0]])
                    raise pytz.AmbiguousTimeError(stamp)

                # Find the index for the switch and pull from a for dst and b
                # for standard
                switch_idx = (delta <= 0).nonzero()[0]
                if switch_idx.size > 1:
                    raise pytz.AmbiguousTimeError(
                        "There are %i dst switches when "
                        "there should only be 1." % switch_idx.size)
                switch_idx = switch_idx[0] + 1 # Pull the only index and adjust
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
                    stamp = Timestamp(vals[i])
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
                stamp = Timestamp(vals[i])
                raise pytz.NonExistentTimeError(stamp)

    return result


cdef inline bisect_right_i8(int64_t *data, int64_t val, Py_ssize_t n):
    cdef Py_ssize_t pivot, left = 0, right = n

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


# Accessors
#----------------------------------------------------------------------

def get_time_micros(ndarray[int64_t] dtindex):
    """
    Datetime as int64 representation to a structured array of fields
    """
    cdef:
        Py_ssize_t i, n = len(dtindex)
        pandas_datetimestruct dts
        ndarray[int64_t] micros

    micros = np.empty(n, dtype=np.int64)

    for i in range(n):
        pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
        micros[i] = 1000000LL * (dts.hour * 60 * 60 +
                                 60 * dts.min + dts.sec) + dts.us

    return micros


cdef int64_t DAY_NS = 86400000000000LL


@cython.wraparound(False)
@cython.boundscheck(False)
def date_normalize(ndarray[int64_t] stamps, tz=None):
    cdef:
        Py_ssize_t i, n = len(stamps)
        pandas_datetimestruct dts
        ndarray[int64_t] result = np.empty(n, dtype=np.int64)

    if tz is not None:
        tz = maybe_get_tz(tz)
        result = _normalize_local(stamps, tz)
    else:
        with nogil:
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                pandas_datetime_to_datetimestruct(
                    stamps[i], PANDAS_FR_ns, &dts)
                result[i] = _normalized_stamp(&dts)

    return result


@cython.wraparound(False)
@cython.boundscheck(False)
cdef _normalize_local(ndarray[int64_t] stamps, object tz):
    cdef:
        Py_ssize_t n = len(stamps)
        ndarray[int64_t] result = np.empty(n, dtype=np.int64)
        ndarray[int64_t] trans, deltas, pos
        pandas_datetimestruct dts

    if is_utc(tz):
        with nogil:
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                pandas_datetime_to_datetimestruct(
                    stamps[i], PANDAS_FR_ns, &dts)
                result[i] = _normalized_stamp(&dts)
    elif is_tzlocal(tz):
        for i in range(n):
            if stamps[i] == NPY_NAT:
                result[i] = NPY_NAT
                continue
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns, &dts)
            dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                          dts.min, dts.sec, dts.us, tz)
            delta = int(get_utcoffset(tz, dt).total_seconds()) * 1000000000
            pandas_datetime_to_datetimestruct(stamps[i] + delta,
                                              PANDAS_FR_ns, &dts)
            result[i] = _normalized_stamp(&dts)
    else:
        # Adjust datetime64 timestamp, recompute datetimestruct
        trans, deltas, typ = get_dst_info(tz)

        _pos = trans.searchsorted(stamps, side='right') - 1
        if _pos.dtype != np.int64:
            _pos = _pos.astype(np.int64)
        pos = _pos

        # statictzinfo
        if typ not in ['pytz', 'dateutil']:
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                pandas_datetime_to_datetimestruct(stamps[i] + deltas[0],
                                                  PANDAS_FR_ns, &dts)
                result[i] = _normalized_stamp(&dts)
        else:
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                pandas_datetime_to_datetimestruct(stamps[i] + deltas[pos[i]],
                                                  PANDAS_FR_ns, &dts)
                result[i] = _normalized_stamp(&dts)

    return result

cdef inline int64_t _normalized_stamp(pandas_datetimestruct *dts) nogil:
    dts.hour = 0
    dts.min = 0
    dts.sec = 0
    dts.us = 0
    dts.ps = 0
    return pandas_datetimestruct_to_datetime(PANDAS_FR_ns, dts)


def dates_normalized(ndarray[int64_t] stamps, tz=None):
    cdef:
        Py_ssize_t i, n = len(stamps)
        ndarray[int64_t] trans, deltas
        pandas_datetimestruct dts

    if tz is None or is_utc(tz):
        for i in range(n):
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns, &dts)
            if (dts.hour + dts.min + dts.sec + dts.us) > 0:
                return False
    elif is_tzlocal(tz):
        for i in range(n):
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns, &dts)
            dt = datetime(dts.year, dts.month, dts.day, dts.hour, dts.min,
                          dts.sec, dts.us, tz)
            dt = dt + tz.utcoffset(dt)
            if (dt.hour + dt.minute + dt.second + dt.microsecond) > 0:
                return False
    else:
        trans, deltas, typ = get_dst_info(tz)

        for i in range(n):
            # Adjust datetime64 timestamp, recompute datetimestruct
            pos = trans.searchsorted(stamps[i]) - 1
            inf = tz._transition_info[pos]

            pandas_datetime_to_datetimestruct(stamps[i] + deltas[pos],
                                              PANDAS_FR_ns, &dts)
            if (dts.hour + dts.min + dts.sec + dts.us) > 0:
                return False

    return True

# Some general helper functions
#----------------------------------------------------------------------


def monthrange(int64_t year, int64_t month):
    cdef:
        int64_t days

    if month < 1 or month > 12:
        raise ValueError("bad month number 0; must be 1-12")

    days = days_per_month_table[is_leapyear(year)][month -1]

    return (dayofweek(year, month, 1), days)


cdef inline int days_in_month(pandas_datetimestruct dts) nogil:
    return days_per_month_table[is_leapyear(dts.year)][dts.month -1]

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


cdef inline int _year_add_months(pandas_datetimestruct dts,
                                 int months) nogil:
    """new year number after shifting pandas_datetimestruct number of months"""
    return dts.year + (dts.month + months - 1) / 12

cdef inline int _month_add_months(pandas_datetimestruct dts,
                                  int months) nogil:
    """
    New month number after shifting pandas_datetimestruct
    number of months.
    """
    cdef int new_month = (dts.month + months) % 12
    return 12 if new_month == 0 else new_month


@cython.wraparound(False)
@cython.boundscheck(False)
def shift_months(int64_t[:] dtindex, int months, object day=None):
    """
    Given an int64-based datetime index, shift all elements
    specified number of months using DateOffset semantics

    day: {None, 'start', 'end'}
       * None: day of month
       * 'start' 1st day of month
       * 'end' last day of month
    """
    cdef:
        Py_ssize_t i
        pandas_datetimestruct dts
        int count = len(dtindex)
        int months_to_roll
        bint roll_check
        int64_t[:] out = np.empty(count, dtype='int64')

    if day is None:
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = NPY_NAT; continue
                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                dts.year = _year_add_months(dts, months)
                dts.month = _month_add_months(dts, months)

                dts.day = min(dts.day, days_in_month(dts))
                out[i] = pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &dts)
    elif day == 'start':
        roll_check = False
        if months <= 0:
            months += 1
            roll_check = True
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = NPY_NAT; continue
                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                months_to_roll = months

                # offset semantics - if on the anchor point and going backwards
                # shift to next
                if roll_check and dts.day == 1:
                    months_to_roll -= 1

                dts.year = _year_add_months(dts, months_to_roll)
                dts.month = _month_add_months(dts, months_to_roll)
                dts.day = 1

                out[i] = pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &dts)
    elif day == 'end':
        roll_check = False
        if months > 0:
            months -= 1
            roll_check = True
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = NPY_NAT; continue
                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                months_to_roll = months

                # similar semantics - when adding shift forward by one
                # month if already at an end of month
                if roll_check and dts.day == days_in_month(dts):
                    months_to_roll += 1

                dts.year = _year_add_months(dts, months_to_roll)
                dts.month = _month_add_months(dts, months_to_roll)

                dts.day = days_in_month(dts)
                out[i] = pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &dts)
    else:
        raise ValueError("day must be None, 'start' or 'end'")

    return np.asarray(out)
