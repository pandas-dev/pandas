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

from util cimport (is_integer_object, is_float_object, is_string_object,
                   is_datetime64_object, is_timedelta64_object,
                   INT64_MAX)
cimport util

from cpython.datetime cimport (PyDelta_Check, PyTZInfo_Check,
                               PyDateTime_Check, PyDate_Check,
                               PyDateTime_IMPORT,
                               timedelta, datetime)
# import datetime C API
PyDateTime_IMPORT
# this is our datetime.pxd
from datetime cimport pandas_datetime_to_datetimestruct, _string_to_dts

# stdlib datetime imports
from datetime import time as datetime_time


from tslibs.np_datetime cimport (check_dts_bounds,
                                 reverse_ops,
                                 cmp_scalar,
                                 pandas_datetimestruct,
                                 PANDAS_DATETIMEUNIT, PANDAS_FR_ns,
                                 dt64_to_dtstruct, dtstruct_to_dt64,
                                 pydatetime_to_dt64, pydate_to_dt64,
                                 npy_datetime,
                                 get_datetime64_unit, get_datetime64_value,
                                 get_timedelta64_value,
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


from tslibs.timedeltas cimport cast_from_unit, delta_to_nanoseconds
from tslibs.timedeltas import Timedelta
from tslibs.timezones cimport (
    is_utc, is_tzlocal, is_fixed_offset,
    treat_tz_as_dateutil, treat_tz_as_pytz,
    get_timezone, get_utcoffset, maybe_get_tz,
    get_dst_info)
from tslibs.fields import (
    get_date_name_field, get_start_end_field, get_date_field,
    build_field_sarray)
from tslibs.conversion cimport (tz_convert_single, _TSObject,
                                convert_to_tsobject,
                                convert_datetime_to_tsobject,
                                get_datetime64_nanos)
from tslibs.conversion import (tz_localize_to_utc,
                               tz_convert_single, date_normalize)

from tslibs.nattype import NaT, nat_strings, iNaT
from tslibs.nattype cimport _checknull_with_nat, NPY_NAT


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

    if box and is_string_object(freq):
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


_zero_time = datetime_time(0, 0)
_no_input = object()

# Python front end to C extension type _Timestamp
# This serves as the box for datetime64


class Timestamp(_Timestamp):
    """Pandas replacement for datetime.datetime

    TimeStamp is the pandas equivalent of python's Datetime
    and is interchangable with it in most cases. It's the type used
    for the entries that make up a DatetimeIndex, and other timeseries
    oriented data structures in pandas.

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

    year, month, day : int
        .. versionadded:: 0.19.0
    hour, minute, second, microsecond : int, optional, default 0
        .. versionadded:: 0.19.0
    tzinfo : datetime.tzinfo, optional, default None
        .. versionadded:: 0.19.0

    Notes
    -----
    There are essentially three calling conventions for the constructor. The
    primary form accepts four parameters. They can be passed by position or
    keyword.

    The other two forms mimic the parameters from ``datetime.datetime``. They
    can be passed by either position or keyword, but not both mixed together.

    Examples
    --------
    >>> pd.Timestamp('2017-01-01T12')
    Timestamp('2017-01-01 12:00:00')

    >>> pd.Timestamp(2017, 1, 1, 12)
    Timestamp('2017-01-01 12:00:00')

    >>> pd.Timestamp(year=2017, month=1, day=1, hour=12)
    Timestamp('2017-01-01 12:00:00')
    """

    @classmethod
    def fromordinal(cls, ordinal, freq=None, tz=None, offset=None):
        """
        Timestamp.fromordinal(ordinal, freq=None, tz=None, offset=None)

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
        Timestamp.now(tz=None)

        Returns new Timestamp object representing current time local to
        tz.

        Parameters
        ----------
        tz : string / timezone object, default None
            Timezone to localize to
        """
        if is_string_object(tz):
            tz = maybe_get_tz(tz)
        return cls(datetime.now(tz))

    @classmethod
    def today(cls, tz=None):
        """
        Timestamp.today(cls, tz=None)

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
        """
        Timestamp.utcnow()

        Return a new Timestamp representing UTC day and time.
        """
        return cls.now('UTC')

    @classmethod
    def utcfromtimestamp(cls, ts):
        """
        Timestamp.utcfromtimestamp(ts)

        Construct a naive UTC datetime from a POSIX timestamp.
        """
        return cls(datetime.utcfromtimestamp(ts))

    @classmethod
    def fromtimestamp(cls, ts):
        """
        Timestamp.fromtimestamp(ts)

        timestamp[, tz] -> tz's local time from POSIX timestamp.
        """
        return cls(datetime.fromtimestamp(ts))

    @classmethod
    def combine(cls, date, time):
        """
        Timsetamp.combine(date, time)

        date, time -> datetime with same date and time fields
        """
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

        if tzinfo is not None:
            if not PyTZInfo_Check(tzinfo):
                # tzinfo must be a datetime.tzinfo object, GH#17690
                raise TypeError('tzinfo must be a datetime.tzinfo object, '
                                'not %s' % type(tzinfo))
            elif tz is not None:
                raise ValueError('Can provide at most one of tz, tzinfo')

        if ts_input is _no_input:
            # User passed keyword arguments.
            if tz is None:
                # Handle the case where the user passes `tz` and not `tzinfo`
                tz = tzinfo
            return Timestamp(datetime(year, month, day, hour or 0,
                                      minute or 0, second or 0,
                                      microsecond or 0, tzinfo),
                             tz=tz)
        elif is_integer_object(freq):
            # User passed positional arguments:
            # Timestamp(year, month, day[, hour[, minute[, second[,
            # microsecond[, tzinfo]]]]])
            return Timestamp(datetime(ts_input, freq, tz, unit or 0,
                                      year or 0, month or 0, day or 0,
                                      hour), tz=hour)

        if tzinfo is not None:
            # User passed tzinfo instead of tz; avoid silently ignoring
            tz, tzinfo = tzinfo, None

        ts = convert_to_tsobject(ts_input, tz, unit, 0, 0)

        if ts.value == NPY_NAT:
            return NaT

        if is_string_object(freq):
            from pandas.tseries.frequencies import to_offset
            freq = to_offset(freq)

        return create_timestamp_from_ts(ts.value, ts.dts, ts.tzinfo, freq)

    def _round(self, freq, rounder):

        cdef:
            int64_t unit, r, value, buff = 1000000
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
        cdef dict wdays = {0: 'Monday', 1: 'Tuesday', 2: 'Wednesday',
                           3: 'Thursday', 4: 'Friday', 5: 'Saturday',
                           6: 'Sunday'}
        return wdays[self.weekday()]

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
            if not is_string_object(ambiguous):
                ambiguous = [ambiguous]
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
        dt64_to_dtstruct(value, &dts)
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
            check_dts_bounds(&dts)

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


# ----------------------------------------------------------------------


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


cpdef object get_value_box(ndarray arr, object loc):
    cdef:
        Py_ssize_t i, sz

    if is_float_object(loc):
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
            ots = other
        elif other is NaT:
            return op == Py_NE
        elif PyDateTime_Check(other):
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
                    if is_datetime64_object(other):
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
                return PyObject_RichCompare(other, self, reverse_ops[op])
            else:
                if op == Py_EQ:
                    return False
                elif op == Py_NE:
                    return True
                raise TypeError('Cannot compare type %r with type %r' %
                                (type(self).__name__, type(other).__name__))

        self._assert_tzawareness_compat(other)
        return cmp_scalar(self.value, ots.value, op)

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

        elif PyDelta_Check(other) or hasattr(other, 'delta'):
            nanos = delta_to_nanoseconds(other)
            result = Timestamp(self.value + nanos,
                               tz=self.tzinfo, freq=self.freq)
            if getattr(other, 'normalize', False):
                result = Timestamp(normalize_date(result))
            return result

        # index/series like
        elif hasattr(other, '_typ'):
            return NotImplemented

        result = datetime.__add__(self, other)
        if PyDateTime_Check(result):
            result = Timestamp(result)
            result.nanosecond = self.nanosecond
        return result

    def __sub__(self, other):
        if (is_timedelta64_object(other) or is_integer_object(other) or
                PyDelta_Check(other) or hasattr(other, 'delta')):
            # `delta` attribute is for offsets.Tick or offsets.Week obj
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
        if (PyDateTime_Check(self)
                and (PyDateTime_Check(other) or is_datetime64_object(other))):
            self = Timestamp(self)
            other = Timestamp(other)

            # validate tz's
            if get_timezone(self.tzinfo) != get_timezone(other.tzinfo):
                raise TypeError("Timestamp subtraction must have the "
                                "same timezones or no timezones")

            # scalar Timestamp/datetime - Timestamp/datetime -> yields a
            # Timedelta
            try:
                return Timedelta(self.value - other.value)
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

    cpdef int _get_field(self, field):
        cdef:
            int64_t val
            ndarray[int32_t] out
        val = self._maybe_convert_value_to_local()
        out = get_date_field(np.array([val], dtype=np.int64), field)
        return int(out[0])

    cpdef _get_start_end_field(self, field):
        cdef:
            int64_t val
            dict kwds

        freq = self.freq
        if freq:
            kwds = freq.kwds
            month_kw = kwds.get('startingMonth', kwds.get('month', 12))
            freqstr = self.freqstr
        else:
            month_kw = 12
            freqstr = None

        val = self._maybe_convert_value_to_local()
        out = get_start_end_field(np.array([val], dtype=np.int64),
                                  field, freqstr, month_kw)
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

    def timestamp(self):
        """Return POSIX timestamp as float."""
        # py27 compat, see GH#17329
        return round(self.value / 1e9, 6)


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

            if _checknull_with_nat(val):
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
            if _checknull_with_nat(val):
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
            if _checknull_with_nat(val):
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


cdef PyTypeObject* td_type = <PyTypeObject*> Timedelta


cdef inline bint is_timedelta(object o):
    return Py_TYPE(o) == td_type  # isinstance(o, Timedelta)


# ----------------------------------------------------------------------
# Conversion routines

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
            iresult[i] = dtstruct_to_dt64(&dts)
            check_dts_bounds(&dts)
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


# ----------------------------------------------------------------------
# Accessors


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
        dt64_to_dtstruct(dtindex[i], &dts)
        micros[i] = 1000000LL * (dts.hour * 60 * 60 +
                                 60 * dts.min + dts.sec) + dts.us

    return micros


# ----------------------------------------------------------------------
# Some general helper functions


def monthrange(int64_t year, int64_t month):
    cdef:
        int64_t days

    if month < 1 or month > 12:
        raise ValueError("bad month number 0; must be 1-12")

    days = days_per_month_table[is_leapyear(year)][month - 1]

    return (dayofweek(year, month, 1), days)


cdef inline int days_in_month(pandas_datetimestruct dts) nogil:
    return days_per_month_table[is_leapyear(dts.year)][dts.month - 1]


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


cdef inline int _year_add_months(pandas_datetimestruct dts, int months) nogil:
    """new year number after shifting pandas_datetimestruct number of months"""
    return dts.year + (dts.month + months - 1) / 12


cdef inline int _month_add_months(pandas_datetimestruct dts, int months) nogil:
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
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                dts.year = _year_add_months(dts, months)
                dts.month = _month_add_months(dts, months)

                dts.day = min(dts.day, days_in_month(dts))
                out[i] = dtstruct_to_dt64(&dts)
    elif day == 'start':
        roll_check = False
        if months <= 0:
            months += 1
            roll_check = True
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                months_to_roll = months

                # offset semantics - if on the anchor point and going backwards
                # shift to next
                if roll_check and dts.day == 1:
                    months_to_roll -= 1

                dts.year = _year_add_months(dts, months_to_roll)
                dts.month = _month_add_months(dts, months_to_roll)
                dts.day = 1

                out[i] = dtstruct_to_dt64(&dts)
    elif day == 'end':
        roll_check = False
        if months > 0:
            months -= 1
            roll_check = True
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT:
                    out[i] = NPY_NAT
                    continue

                dt64_to_dtstruct(dtindex[i], &dts)
                months_to_roll = months

                # similar semantics - when adding shift forward by one
                # month if already at an end of month
                if roll_check and dts.day == days_in_month(dts):
                    months_to_roll += 1

                dts.year = _year_add_months(dts, months_to_roll)
                dts.month = _month_add_months(dts, months_to_roll)

                dts.day = days_in_month(dts)
                out[i] = dtstruct_to_dt64(&dts)
    else:
        raise ValueError("day must be None, 'start' or 'end'")

    return np.asarray(out)
