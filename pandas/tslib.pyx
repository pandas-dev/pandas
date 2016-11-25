# cython: profile=False

import warnings

cimport numpy as np
from numpy cimport (int8_t, int32_t, int64_t, import_array, ndarray,
                    NPY_INT64, NPY_DATETIME, NPY_TIMEDELTA)
from datetime cimport get_datetime64_value, get_timedelta64_value
import numpy as np

# GH3363
from sys import version_info
cdef bint PY2 = version_info[0] == 2
cdef bint PY3 = not PY2

from cpython cimport (
    PyTypeObject,
    PyFloat_Check,
    PyLong_Check,
    PyObject_RichCompareBool,
    PyObject_RichCompare,
    Py_GT, Py_GE, Py_EQ, Py_NE, Py_LT, Py_LE,
    PyUnicode_Check,
    PyUnicode_AsUTF8String,
)


# Cython < 0.17 doesn't have this in cpython
cdef extern from "Python.h":
    cdef PyTypeObject *Py_TYPE(object)
    int PySlice_Check(object)

cdef extern from "datetime_helper.h":
    double total_seconds(object)

# this is our datetime.pxd
from datetime cimport cmp_pandas_datetimestruct
from libc.stdlib cimport free

from util cimport (is_integer_object, is_float_object, is_datetime64_object,
                   is_timedelta64_object, INT64_MAX)
cimport util

from datetime cimport *
from khash cimport *
cimport cython

from datetime import timedelta, datetime
from datetime import time as datetime_time

import re

# dateutil compat
from dateutil.tz import (tzoffset, tzlocal as _dateutil_tzlocal,
                         tzfile as _dateutil_tzfile,
                         tzutc as _dateutil_tzutc,
                         tzstr as _dateutil_tzstr)

from pandas.compat import is_platform_windows
if is_platform_windows():
    from dateutil.zoneinfo import gettz as _dateutil_gettz
else:
    from dateutil.tz import gettz as _dateutil_gettz
from dateutil.relativedelta import relativedelta
from dateutil.parser import DEFAULTPARSER

from pytz.tzinfo import BaseTzInfo as _pytz_BaseTzInfo
from pandas.compat import (parse_date, string_types, iteritems,
                           StringIO, callable)

import operator
import collections
import warnings

# initialize numpy
import_array()
#import_ufunc()

# import datetime C API
PyDateTime_IMPORT

# in numpy 1.7, will prob need the following:
# numpy_pydatetime_import

cdef int64_t NPY_NAT = util.get_nat()
iNaT = NPY_NAT

# < numpy 1.7 compat for NaT
compat_NaT = np.array([NPY_NAT]).astype('m8[ns]').item()


try:
    basestring
except NameError: # py3
    basestring = str


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
        if _is_utc(tz):
            for i in range(n):
                value = arr[i]
                if value == NPY_NAT:
                    result[i] = NaT
                else:
                    pandas_datetime_to_datetimestruct(
                        value, PANDAS_FR_ns, &dts)
                    result[i] = func_create(value, dts, tz, freq)
        elif _is_tzlocal(tz) or _is_fixed_offset(tz):
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
            trans, deltas, typ = _get_dst_info(tz)

            for i in range(n):

                value = arr[i]
                if value == NPY_NAT:
                    result[i] = NaT
                else:

                    # Adjust datetime64 timestamp, recompute datetimestruct
                    pos = trans.searchsorted(value, side='right') - 1
                    if _treat_tz_as_pytz(tz):
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


cdef inline bint _is_tzlocal(object tz):
    return isinstance(tz, _dateutil_tzlocal)


cdef inline bint _is_fixed_offset(object tz):
    if _treat_tz_as_dateutil(tz):
        if len(tz._trans_idx) == 0 and len(tz._trans_list) == 0:
            return 1
        else:
            return 0
    elif _treat_tz_as_pytz(tz):
        if (len(tz._transition_info) == 0
            and len(tz._utc_transition_times) == 0):
            return 1
        else:
            return 0
    return 1

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
        if isinstance(tz, basestring):
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

        cdef int64_t unit
        cdef object result, value

        from pandas.tseries.frequencies import to_offset
        unit = to_offset(freq).nanos
        if self.tz is not None:
            value = self.tz_localize(None).value
        else:
            value = self.value
        result = Timestamp(unit * rounder(value / float(unit)), unit='ns')
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
        from pandas.tseries.period import Period

        if freq is None:
            freq = self.freq

        return Period(self, freq=freq)

    @property
    def dayofweek(self):
        return self.weekday()

    @property
    def weekday_name(self):
        out = get_date_name_field(
            np.array([self.value], dtype=np.int64), 'weekday_name')
        return out[0]

    @property
    def dayofyear(self):
        return self._get_field('doy')

    @property
    def week(self):
        return self._get_field('woy')

    weekofyear = week

    @property
    def microsecond(self):
        return self._get_field('us')

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
            if not isinstance(ambiguous, basestring):
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

    def replace(self, **kwds):
        """
        implements datetime.replace, handles nanoseconds

        Parameters
        ----------
        kwargs: key-value dict

        accepted keywords are:
        year, month, day, hour, minute, second, microsecond, nanosecond, tzinfo

        values must be integer, or for tzinfo, a tz-convertible

        Returns
        -------
        Timestamp with fields replaced
        """

        cdef:
            pandas_datetimestruct dts
            int64_t value
            object tzinfo, result, k, v
            _TSObject ts

        # set to naive if needed
        tzinfo = self.tzinfo
        value = self.value
        if tzinfo is not None:
            value = tz_convert_single(value, 'UTC', tzinfo)

        # setup components
        pandas_datetime_to_datetimestruct(value, PANDAS_FR_ns, &dts)
        dts.ps = self.nanosecond * 1000

        # replace
        def validate(k, v):
            """ validate integers """
            if not isinstance(v, int):
                raise ValueError("value must be an integer, received "
                                 "{v} for {k}".format(v=type(v), k=k))
            return v

        for k, v in kwds.items():
            if k == 'year':
                dts.year = validate(k, v)
            elif k == 'month':
                dts.month = validate(k, v)
            elif k == 'day':
                dts.day = validate(k, v)
            elif k == 'hour':
                dts.hour = validate(k, v)
            elif k == 'minute':
                dts.min = validate(k, v)
            elif k == 'second':
                dts.sec = validate(k, v)
            elif k == 'microsecond':
                dts.us = validate(k, v)
            elif k == 'nanosecond':
                dts.ps = validate(k, v) * 1000
            elif k == 'tzinfo':
                tzinfo = v
            else:
                raise ValueError("invalid name {} passed".format(k))

        # reconstruct & check bounds
        value = pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &dts)
        if value != NPY_NAT:
            _check_dts_bounds(&dts)

        # set tz if needed
        if tzinfo is not None:
            value = tz_convert_single(value, tzinfo, 'UTC')

        result = create_timestamp_from_ts(value, dts, tzinfo, self.freq)

        return result

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


class NaTType(_NaT):
    """(N)ot-(A)-(T)ime, the time equivalent of NaN"""

    def __new__(cls):
        cdef _NaT base

        base = _NaT.__new__(cls, 1, 1, 1)
        base._day = -1
        base._month = -1
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
        # GH 10939
        return np.nan

    @property
    def is_leap_year(self):
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


fields = ['year', 'quarter', 'month', 'day', 'hour',
          'minute', 'second', 'millisecond', 'microsecond', 'nanosecond',
          'week', 'dayofyear', 'days_in_month', 'daysinmonth', 'dayofweek',
          'weekday_name']
for field in fields:
    prop = property(fget=lambda self: np.nan)
    setattr(NaTType, field, prop)


# GH9513 NaT methods (except to_datetime64) to raise, return np.nan, or
# return NaT create functions that raise, for binding to NaTType
def _make_error_func(func_name):
    def f(*args, **kwargs):
        raise ValueError("NaTType does not support " + func_name)
    f.__name__ = func_name
    return f


def _make_nat_func(func_name):
    def f(*args, **kwargs):
        return NaT
    f.__name__ = func_name
    return f


def _make_nan_func(func_name):
    def f(*args, **kwargs):
        return np.nan
    f.__name__ = func_name
    return f

_nat_methods = ['date', 'now', 'replace', 'to_pydatetime', 'today']

_nan_methods = ['weekday', 'isoweekday', 'total_seconds']

_implemented_methods = ['to_datetime', 'to_datetime64', 'isoformat']
_implemented_methods.extend(_nat_methods)
_implemented_methods.extend(_nan_methods)

for _method_name in _nat_methods:
    # not all methods exist in all versions of Python
    if hasattr(NaTType, _method_name):
        setattr(NaTType, _method_name, _make_nat_func(_method_name))

for _method_name in _nan_methods:
    if hasattr(NaTType, _method_name):
        setattr(NaTType, _method_name, _make_nan_func(_method_name))

for _maybe_method_name in dir(NaTType):
    _maybe_method = getattr(NaTType, _maybe_method_name)
    if (callable(_maybe_method)
        and not _maybe_method_name.startswith("_")
        and _maybe_method_name not in _implemented_methods):
        setattr(NaTType, _maybe_method_name,
                _make_error_func(_maybe_method_name))


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
    if PyFloat_Check(val):
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

cdef pandas_datetimestruct _NS_MIN_DTS, _NS_MAX_DTS
pandas_datetime_to_datetimestruct(_NS_LOWER_BOUND, PANDAS_FR_ns, &_NS_MIN_DTS)
pandas_datetime_to_datetimestruct(_NS_UPPER_BOUND, PANDAS_FR_ns, &_NS_MAX_DTS)

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


cdef inline bint _is_multiple(int64_t us, int64_t mult):
    return us % mult == 0


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
                zone = _get_zone(self.tzinfo)
        except ValueError:
            year2000 = self.replace(year=2000)
            stamp += year2000.strftime('%z')
            if self.tzinfo:
                zone = _get_zone(self.tzinfo)

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
        cdef:
            pandas_datetimestruct dts
            _TSObject ts

        if self.nanosecond != 0 and warn:
            warnings.warn("Discarding nonzero nanoseconds in conversion",
                          UserWarning, stacklevel=2)
        ts = convert_to_tsobject(self, self.tzinfo, None, 0, 0)
        dts = ts.dts
        return datetime(dts.year, dts.month, dts.day,
                        dts.hour, dts.min, dts.sec,
                        dts.us, ts.tzinfo)

    cpdef to_datetime64(self):
        """ Returns a numpy.datetime64 object with 'ns' precision """
        return np.datetime64(self.value, 'ns')

    def __add__(self, other):
        cdef int64_t other_int

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

    cpdef _get_field(self, field):
        out = get_date_field(np.array([self.value], dtype=np.int64), field)
        return int(out[0])

    cpdef _get_start_end_field(self, field):
        month_kw = self.freq.kwds.get(
            'startingMonth', self.freq.kwds.get(
                'month', 12)) if self.freq else 12
        freqstr = self.freqstr if self.freq else None
        out = get_start_end_field(
            np.array([self.value], dtype=np.int64), field, freqstr, month_kw)
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


def _delta_to_nanoseconds(delta):
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


# lightweight C object to hold datetime & int64 pair
cdef class _TSObject:
    cdef:
        pandas_datetimestruct dts      # pandas_datetimestruct
        int64_t value               # numpy dt64
        object tzinfo

    property value:
        def __get__(self):
            return self.value

cpdef _get_utcoffset(tzinfo, obj):
    try:
        return tzinfo._utcoffset
    except AttributeError:
        return tzinfo.utcoffset(obj)

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
        if tz is not None:
            # sort of a temporary hack
            if ts.tzinfo is not None:
                if (hasattr(tz, 'normalize') and
                    hasattr(ts.tzinfo, '_utcoffset')):
                    ts = tz.normalize(ts)
                    obj.value = _pydatetime_to_dts(ts, &obj.dts)
                    obj.tzinfo = ts.tzinfo
                else: #tzoffset
                    try:
                        tz = ts.astimezone(tz).tzinfo
                    except:
                        pass
                    obj.value = _pydatetime_to_dts(ts, &obj.dts)
                    ts_offset = _get_utcoffset(ts.tzinfo, ts)
                    obj.value -= _delta_to_nanoseconds(ts_offset)
                    tz_offset = _get_utcoffset(tz, ts)
                    obj.value += _delta_to_nanoseconds(tz_offset)
                    pandas_datetime_to_datetimestruct(obj.value,
                                                      PANDAS_FR_ns, &obj.dts)
                    obj.tzinfo = tz
            elif not _is_utc(tz):
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

        if obj.tzinfo is not None and not _is_utc(obj.tzinfo):
            offset = _get_utcoffset(obj.tzinfo, ts)
            obj.value -= _delta_to_nanoseconds(offset)

        if is_timestamp(ts):
            obj.value += ts.nanosecond
            obj.dts.ps = ts.nanosecond * 1000
        _check_dts_bounds(&obj.dts)
        return obj
    elif PyDate_Check(ts):
        # Keep the converter same as PyDateTime's
        ts = datetime.combine(ts, datetime_time())
        return convert_to_tsobject(ts, tz, None, 0, 0)
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
        ts = Timestamp.now(tz)
    elif ts == 'today':
        # Issue 9000, we short-circuit rather than going
        # into np_datetime_strings which returns a normalized datetime
        ts = Timestamp.today(tz)
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
                    ts = Timestamp(obj.value, tz=obj.tzinfo)
            else:
                ts = obj.value
                if tz is not None:
                    # shift for _localize_tso
                    ts = tz_convert_single(ts, tz, 'UTC')
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
    if _is_utc(tz):
        obj.tzinfo = tz
    elif _is_tzlocal(tz):
        pandas_datetime_to_datetimestruct(obj.value, PANDAS_FR_ns, &obj.dts)
        dt = datetime(obj.dts.year, obj.dts.month, obj.dts.day, obj.dts.hour,
                      obj.dts.min, obj.dts.sec, obj.dts.us, tz)
        delta = int(total_seconds(_get_utcoffset(tz, dt))) * 1000000000
        if obj.value != NPY_NAT:
            pandas_datetime_to_datetimestruct(obj.value + delta,
                                              PANDAS_FR_ns, &obj.dts)
        else:
            pandas_datetime_to_datetimestruct(obj.value,
                                              PANDAS_FR_ns, &obj.dts)
        obj.tzinfo = tz
    else:
        # Adjust datetime64 timestamp, recompute datetimestruct
        trans, deltas, typ = _get_dst_info(tz)

        pos = trans.searchsorted(obj.value, side='right') - 1

        # static/pytz/dateutil specific code
        if _is_fixed_offset(tz):
            # statictzinfo
            if len(deltas) > 0 and obj.value != NPY_NAT:
                pandas_datetime_to_datetimestruct(obj.value + deltas[0],
                                                  PANDAS_FR_ns, &obj.dts)
            else:
                pandas_datetime_to_datetimestruct(
                    obj.value, PANDAS_FR_ns, &obj.dts)
            obj.tzinfo = tz
        elif _treat_tz_as_pytz(tz):
            inf = tz._transition_info[pos]
            if obj.value != NPY_NAT:
                pandas_datetime_to_datetimestruct(obj.value + deltas[pos],
                                                  PANDAS_FR_ns, &obj.dts)
            else:
                pandas_datetime_to_datetimestruct(obj.value,
                                                  PANDAS_FR_ns, &obj.dts)
            obj.tzinfo = tz._tzinfos[inf]
        elif _treat_tz_as_dateutil(tz):
            if obj.value != NPY_NAT:
                pandas_datetime_to_datetimestruct(obj.value + deltas[pos],
                                                  PANDAS_FR_ns, &obj.dts)
            else:
                pandas_datetime_to_datetimestruct(obj.value,
                                                  PANDAS_FR_ns, &obj.dts)
            obj.tzinfo = tz
        else:
            obj.tzinfo = tz


def _localize_pydatetime(object dt, object tz):
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


def get_timezone(tz):
    return _get_zone(tz)

cdef inline bint _is_utc(object tz):
    return tz is UTC or isinstance(tz, _dateutil_tzutc)

cdef inline object _get_zone(object tz):
    """
    We need to do several things here:
    1) Distinguish between pytz and dateutil timezones
    2) Not be over-specific (e.g. US/Eastern with/without DST is same *zone*
       but a different tz object)
    3) Provide something to serialize when we're storing a datetime object
       in pytables.

    We return a string prefaced with dateutil if it's a dateutil tz, else just
    the tz name. It needs to be a string so that we can serialize it with
    UJSON/pytables. maybe_get_tz (below) is the inverse of this process.
    """
    if _is_utc(tz):
        return 'UTC'
    else:
        if _treat_tz_as_dateutil(tz):
            if '.tar.gz' in tz._filename:
                raise ValueError(
                    'Bad tz filename. Dateutil on python 3 on windows has a '
                    'bug which causes tzfile._filename to be the same for all '
                    'timezone files. Please construct dateutil timezones '
                    'implicitly by passing a string like "dateutil/Europe'
                    '/London" when you construct your pandas objects instead '
                    'of passing a timezone object. See '
                    'https://github.com/pydata/pandas/pull/7362')
            return 'dateutil/' + tz._filename
        else:
            # tz is a pytz timezone or unknown.
            try:
                zone = tz.zone
                if zone is None:
                    return tz
                return zone
            except AttributeError:
                return tz


cpdef inline object maybe_get_tz(object tz):
    """
    (Maybe) Construct a timezone object from a string. If tz is a string, use
    it to construct a timezone object. Otherwise, just return tz.
    """
    if isinstance(tz, string_types):
        if tz == 'tzlocal()':
            tz = _dateutil_tzlocal()
        elif tz.startswith('dateutil/'):
            zone = tz[9:]
            tz = _dateutil_gettz(zone)
            # On Python 3 on Windows, the filename is not always set correctly.
            if isinstance(tz, _dateutil_tzfile) and '.tar.gz' in tz._filename:
                tz._filename = zone
        else:
            tz = pytz.timezone(tz)
    elif is_integer_object(tz):
        tz = pytz.FixedOffset(tz / 60)
    return tz


class OutOfBoundsDatetime(ValueError):
    pass

cdef inline _check_dts_bounds(pandas_datetimestruct *dts):
    cdef:
        bint error = False

    if dts.year <= 1677 and cmp_pandas_datetimestruct(dts, &_NS_MIN_DTS) == -1:
        error = True
    elif (
            dts.year >= 2262 and
            cmp_pandas_datetimestruct(dts, &_NS_MAX_DTS) == 1):
        error = True

    if error:
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
                    if _get_zone(val.tzinfo) != inferred_tz:
                        raise ValueError('Array must be all same time zone')
                else:
                    inferred_tz = _get_zone(val.tzinfo)

                _ts = convert_to_tsobject(val, None, None, 0, 0)
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

cdef:
    set _not_datelike_strings = set(['a', 'A', 'm', 'M', 'p', 'P', 't', 'T'])

cpdef bint _does_string_look_like_datetime(object date_string):
    if date_string.startswith('0'):
        # Strings starting with 0 are more consistent with a
        # date-like string than a number
        return True

    try:
        if float(date_string) < 1000:
            return False
    except ValueError:
        pass

    if date_string in _not_datelike_strings:
        return False

    return True


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


class DateParseError(ValueError):
    pass


cdef object _TIMEPAT = re.compile(r'^([01]?[0-9]|2[0-3]):([0-5][0-9])')


def parse_datetime_string(object date_string, object freq=None,
                          dayfirst=False, yearfirst=False, **kwargs):
    """parse datetime string, only returns datetime.
    Also cares special handling matching time patterns.

    Returns
    -------
    datetime
    """

    cdef:
        object dt

    if not _does_string_look_like_datetime(date_string):
        raise ValueError('Given date string not likely a datetime.')

    if _TIMEPAT.match(date_string):
        # use current datetime as default, not pass _DEFAULT_DATETIME
        dt = parse_date(date_string, dayfirst=dayfirst,
                        yearfirst=yearfirst, **kwargs)
        return dt
    try:
        dt, _, _ = _parse_dateabbr_string(date_string, _DEFAULT_DATETIME, freq)
        return dt
    except DateParseError:
        raise
    except ValueError:
        pass

    try:
        dt = parse_date(date_string, default=_DEFAULT_DATETIME,
                        dayfirst=dayfirst, yearfirst=yearfirst, **kwargs)
    except TypeError:
        # following may be raised from dateutil
        # TypeError: 'NoneType' object is not iterable
        raise ValueError('Given date string not likely a datetime.')

    return dt


def parse_datetime_string_with_reso(object date_string, object freq=None,
                                    dayfirst=False, yearfirst=False, **kwargs):
    """parse datetime string, only returns datetime

    Returns
    -------
    datetime
    """

    cdef:
        object parsed, reso

    if not _does_string_look_like_datetime(date_string):
        raise ValueError('Given date string not likely a datetime.')

    try:
        return _parse_dateabbr_string(date_string, _DEFAULT_DATETIME, freq)
    except DateParseError:
        raise
    except ValueError:
        pass

    try:
        parsed, reso = dateutil_parse(date_string, _DEFAULT_DATETIME,
                                      dayfirst=dayfirst, yearfirst=yearfirst)
    except Exception as e:
        # TODO: allow raise of errors within instead
        raise DateParseError(e)
    if parsed is None:
        raise DateParseError("Could not parse %s" % date_string)
    return parsed, parsed, reso


cdef inline object _parse_dateabbr_string(object date_string, object default,
                                          object freq):
    cdef:
        object ret
        int year, quarter, month, mnum, date_len

    # special handling for possibilities eg, 2Q2005, 2Q05, 2005Q1, 05Q1
    assert util.is_string_object(date_string)

    # len(date_string) == 0
    # should be NaT???

    if date_string in _nat_strings:
        return NaT, NaT, ''

    date_string = date_string.upper()
    date_len = len(date_string)

    if date_len == 4:
        # parse year only like 2000
        try:
            ret = default.replace(year=int(date_string))
            return ret, ret, 'year'
        except ValueError:
            pass

    try:
        if 4 <= date_len <= 7:
            i = date_string.index('Q', 1, 6)
            if i == 1:
                quarter = int(date_string[0])
                if date_len == 4 or (date_len == 5
                                     and date_string[i + 1] == '-'):
                    # r'(\d)Q-?(\d\d)')
                    year = 2000 + int(date_string[-2:])
                elif date_len == 6 or (date_len == 7
                                       and date_string[i + 1] == '-'):
                    # r'(\d)Q-?(\d\d\d\d)')
                    year = int(date_string[-4:])
                else:
                    raise ValueError
            elif i == 2 or i == 3:
                # r'(\d\d)-?Q(\d)'
                if date_len == 4 or (date_len == 5
                                     and date_string[i - 1] == '-'):
                    quarter = int(date_string[-1])
                    year = 2000 + int(date_string[:2])
                else:
                    raise ValueError
            elif i == 4 or i == 5:
                if date_len == 6 or (date_len == 7
                                     and date_string[i - 1] == '-'):
                    # r'(\d\d\d\d)-?Q(\d)'
                    quarter = int(date_string[-1])
                    year = int(date_string[:4])
                else:
                    raise ValueError

            if not (1 <= quarter <= 4):
                msg = ('Incorrect quarterly string is given, quarter must be '
                       'between 1 and 4: {0}')
                raise DateParseError(msg.format(date_string))

            if freq is not None:
                # hack attack, #1228
                try:
                    mnum = _MONTH_NUMBERS[_get_rule_month(freq)] + 1
                except (KeyError, ValueError):
                    msg = ('Unable to retrieve month information from given '
                           'freq: {0}').format(freq)
                    raise DateParseError(msg)

                month = (mnum + (quarter - 1) * 3) % 12 + 1
                if month > mnum:
                    year -= 1
            else:
                month = (quarter - 1) * 3 + 1

            ret = default.replace(year=year, month=month)
            return ret, ret, 'quarter'

    except DateParseError:
        raise
    except ValueError:
        pass

    if date_len == 6 and (freq == 'M' or getattr(
            freq, 'rule_code', None) == 'M'):
        year = int(date_string[:4])
        month = int(date_string[4:6])
        try:
            ret = default.replace(year=year, month=month)
            return ret, ret, 'month'
        except ValueError:
            pass

    for pat in ['%Y-%m', '%m-%Y', '%b %Y', '%b-%Y']:
        try:
            ret = datetime.strptime(date_string, pat)
            return ret, ret, 'month'
        except ValueError:
            pass

    raise ValueError('Unable to parse {0}'.format(date_string))


def dateutil_parse(object timestr, object default, ignoretz=False,
                   tzinfos=None, **kwargs):
    """ lifted from dateutil to get resolution"""

    cdef:
        object fobj, res, attr, ret, tzdata
        object reso = None
        dict repl = {}

    fobj = StringIO(str(timestr))
    res = DEFAULTPARSER._parse(fobj, **kwargs)

    # dateutil 2.2 compat
    if isinstance(res, tuple):
        res, _ = res

    if res is None:
        msg = "Unknown datetime string format, unable to parse: {0}"
        raise ValueError(msg.format(timestr))

    for attr in ["year", "month", "day", "hour",
                 "minute", "second", "microsecond"]:
        value = getattr(res, attr)
        if value is not None:
            repl[attr] = value
            reso = attr

    if reso is None:
        msg = "Unable to parse datetime string: {0}"
        raise ValueError(msg.format(timestr))

    if reso == 'microsecond':
        if repl['microsecond'] == 0:
            reso = 'second'
        elif repl['microsecond'] % 1000 == 0:
            reso = 'millisecond'

    ret = default.replace(**repl)
    if res.weekday is not None and not res.day:
        ret = ret + relativedelta.relativedelta(weekday=res.weekday)
    if not ignoretz:
        if callable(tzinfos) or tzinfos and res.tzname in tzinfos:
            if callable(tzinfos):
                tzdata = tzinfos(res.tzname, res.tzoffset)
            else:
                tzdata = tzinfos.get(res.tzname)
            if isinstance(tzdata, datetime.tzinfo):
                tzinfo = tzdata
            elif isinstance(tzdata, string_types):
                tzinfo = _dateutil_tzstr(tzdata)
            elif isinstance(tzdata, int):
                tzinfo = tzoffset(res.tzname, tzdata)
            else:
                raise ValueError("offset must be tzinfo subclass, "
                                 "tz string, or int offset")
            ret = ret.replace(tzinfo=tzinfo)
        elif res.tzname and res.tzname in time.tzname:
            ret = ret.replace(tzinfo=_dateutil_tzlocal())
        elif res.tzoffset == 0:
            ret = ret.replace(tzinfo=_dateutil_tzutc())
        elif res.tzoffset:
            ret = ret.replace(tzinfo=tzoffset(res.tzname, res.tzoffset))
    return ret, reso


# const for parsers

_DEFAULT_DATETIME = datetime(1, 1, 1).replace(
    hour=0, minute=0, second=0, microsecond=0)
_MONTHS = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL',
           'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
_MONTH_NUMBERS = dict((k, i) for i, k in enumerate(_MONTHS))
_MONTH_ALIASES = dict((k + 1, v) for k, v in enumerate(_MONTHS))


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
                        _ts = convert_to_tsobject(val, None, 'ns', 0, 0)
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
                        _ts = convert_to_tsobject(py_dt, None, None, 0, 0)
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
            if isinstance(other, _NaT):
                return _cmp_nat_dt(other, self, _reverse_ops[op])
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

            kwargs = dict([ (k, _to_py_int_float(v))
                            for k, v in iteritems(kwargs) ])

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

        return Timedelta(other *self.value, unit='ns')

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

    def _not_implemented(self, *args, **kwargs):
        return NotImplemented

    __floordiv__ = _not_implemented
    __rfloordiv__ = _not_implemented

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
        ts = np.timedelta64(ts)
    elif not isinstance(ts, np.timedelta64):
        raise ValueError("Invalid type for timedelta "
                         "scalar: %s" % type(ts))
    return ts.astype('timedelta64[ns]')


def array_strptime(ndarray[object] values, object fmt,
                   bint exact=True, errors='raise'):
    """
    Parameters
    ----------
    values : ndarray of string-like objects
    fmt : string-like regex
    exact : matches must be exact if True, search if False
    coerce : if invalid values found, coerce to NaT
    """

    cdef:
        Py_ssize_t i, n = len(values)
        pandas_datetimestruct dts
        ndarray[int64_t] iresult
        int year, month, day, minute, hour, second, weekday, julian, tz
        int week_of_year, week_of_year_start
        int64_t us, ns
        object val, group_key, ampm, found
        dict found_key
        bint is_raise = errors=='raise'
        bint is_ignore = errors=='ignore'
        bint is_coerce = errors=='coerce'

    assert is_raise or is_ignore or is_coerce

    global _TimeRE_cache, _regex_cache
    with _cache_lock:
        if _getlang() != _TimeRE_cache.locale_time.lang:
            _TimeRE_cache = TimeRE()
            _regex_cache.clear()
        if len(_regex_cache) > _CACHE_MAX_SIZE:
            _regex_cache.clear()
        locale_time = _TimeRE_cache.locale_time
        format_regex = _regex_cache.get(fmt)
        if not format_regex:
            try:
                format_regex = _TimeRE_cache.compile(fmt)
            # KeyError raised when a bad format is found; can be specified as
            # \\, in which case it was a stray % but with a space after it
            except KeyError, err:
                bad_directive = err.args[0]
                if bad_directive == "\\":
                    bad_directive = "%"
                del err
                raise ValueError("'%s' is a bad directive in format '%s'" %
                                    (bad_directive, fmt))
            # IndexError only occurs when the format string is "%"
            except IndexError:
                raise ValueError("stray %% in format '%s'" % fmt)
            _regex_cache[fmt] = format_regex

    result = np.empty(n, dtype='M8[ns]')
    iresult = result.view('i8')

    dts.us = dts.ps = dts.as = 0

    cdef dict _parse_code_table = {
        'y': 0,
        'Y': 1,
        'm': 2,
        'B': 3,
        'b': 4,
        'd': 5,
        'H': 6,
        'I': 7,
        'M': 8,
        'S': 9,
        'f': 10,
        'A': 11,
        'a': 12,
        'w': 13,
        'j': 14,
        'U': 15,
        'W': 16,
        'Z': 17,
        'p': 18   # just an additional key, works only with I
    }
    cdef int parse_code

    for i in range(n):
        val = values[i]
        if util.is_string_object(val):
            if val in _nat_strings:
                iresult[i] = NPY_NAT
                continue
        else:
            if _checknull_with_nat(val):
                iresult[i] = NPY_NAT
                continue
            else:
                val = str(val)

        # exact matching
        if exact:
            found = format_regex.match(val)
            if not found:
                if is_coerce:
                    iresult[i] = NPY_NAT
                    continue
                raise ValueError("time data %r does not match "
                                 "format %r (match)" % (values[i], fmt))
            if len(val) != found.end():
                if is_coerce:
                    iresult[i] = NPY_NAT
                    continue
                raise ValueError("unconverted data remains: %s" %
                                  values[i][found.end():])

        # search
        else:
            found = format_regex.search(val)
            if not found:
                if is_coerce:
                    iresult[i] = NPY_NAT
                    continue
                raise ValueError("time data %r does not match format "
                                 "%r (search)" % (values[i], fmt))

        year = 1900
        month = day = 1
        hour = minute = second = ns = us = 0
        tz = -1
        # Default to -1 to signify that values not known; not critical to have,
        # though
        week_of_year = -1
        week_of_year_start = -1
        # weekday and julian defaulted to -1 so as to signal need to calculate
        # values
        weekday = julian = -1
        found_dict = found.groupdict()
        for group_key in found_dict.iterkeys():
            # Directives not explicitly handled below:
            #   c, x, X
            #      handled by making out of other directives
            #   U, W
            #      worthless without day of the week
            parse_code = _parse_code_table[group_key]

            if parse_code == 0:
                year = int(found_dict['y'])
                # Open Group specification for strptime() states that a %y
                #value in the range of [00, 68] is in the century 2000, while
                #[69,99] is in the century 1900
                if year <= 68:
                    year += 2000
                else:
                    year += 1900
            elif parse_code == 1:
                year = int(found_dict['Y'])
            elif parse_code == 2:
                month = int(found_dict['m'])
            elif parse_code == 3:
            # elif group_key == 'B':
                month = locale_time.f_month.index(found_dict['B'].lower())
            elif parse_code == 4:
            # elif group_key == 'b':
                month = locale_time.a_month.index(found_dict['b'].lower())
            elif parse_code == 5:
            # elif group_key == 'd':
                day = int(found_dict['d'])
            elif parse_code == 6:
            # elif group_key == 'H':
                hour = int(found_dict['H'])
            elif parse_code == 7:
                hour = int(found_dict['I'])
                ampm = found_dict.get('p', '').lower()
                # If there was no AM/PM indicator, we'll treat this like AM
                if ampm in ('', locale_time.am_pm[0]):
                    # We're in AM so the hour is correct unless we're
                    # looking at 12 midnight.
                    # 12 midnight == 12 AM == hour 0
                    if hour == 12:
                        hour = 0
                elif ampm == locale_time.am_pm[1]:
                    # We're in PM so we need to add 12 to the hour unless
                    # we're looking at 12 noon.
                    # 12 noon == 12 PM == hour 12
                    if hour != 12:
                        hour += 12
            elif parse_code == 8:
                minute = int(found_dict['M'])
            elif parse_code == 9:
                second = int(found_dict['S'])
            elif parse_code == 10:
                s = found_dict['f']
                # Pad to always return nanoseconds
                s += "0" * (9 - len(s))
                us = long(s)
                ns = us % 1000
                us = us / 1000
            elif parse_code == 11:
                weekday = locale_time.f_weekday.index(found_dict['A'].lower())
            elif parse_code == 12:
                weekday = locale_time.a_weekday.index(found_dict['a'].lower())
            elif parse_code == 13:
                weekday = int(found_dict['w'])
                if weekday == 0:
                    weekday = 6
                else:
                    weekday -= 1
            elif parse_code == 14:
                julian = int(found_dict['j'])
            elif parse_code == 15 or parse_code == 16:
                week_of_year = int(found_dict[group_key])
                if group_key == 'U':
                    # U starts week on Sunday.
                    week_of_year_start = 6
                else:
                    # W starts week on Monday.
                    week_of_year_start = 0
            elif parse_code == 17:
                # Since -1 is default value only need to worry about setting tz
                # if it can be something other than -1.
                found_zone = found_dict['Z'].lower()
                for value, tz_values in enumerate(locale_time.timezone):
                    if found_zone in tz_values:
                        # Deal w/ bad locale setup where timezone names are the
                        # same and yet time.daylight is true; too ambiguous to
                        # be able to tell what timezone has daylight savings
                        if (time.tzname[0] == time.tzname[1] and
                            time.daylight and found_zone not in (
                                "utc", "gmt")):
                            break
                        else:
                            tz = value
                            break
        # If we know the wk of the year and what day of that wk, we can figure
        # out the Julian day of the year.
        if julian == -1 and week_of_year != -1 and weekday != -1:
            week_starts_Mon = True if week_of_year_start == 0 else False
            julian = _calc_julian_from_U_or_W(year, week_of_year, weekday,
                                                week_starts_Mon)
        # Cannot pre-calculate datetime_date() since can change in Julian
        # calculation and thus could have different value for the day of the wk
        # calculation.
        try:
            if julian == -1:
                # Need to add 1 to result since first day of the year is 1, not
                # 0.
                julian = datetime_date(year, month, day).toordinal() - \
                    datetime_date(year, 1, 1).toordinal() + 1
            else: # Assume that if they bothered to include Julian day it will
                # be accurate.
                datetime_result = datetime_date.fromordinal(
                    (julian - 1) + datetime_date(year, 1, 1).toordinal())
                year = datetime_result.year
                month = datetime_result.month
                day = datetime_result.day
        except ValueError:
            if is_coerce:
                iresult[i] = NPY_NAT
                continue
            raise
        if weekday == -1:
            weekday = datetime_date(year, month, day).weekday()

        dts.year = year
        dts.month = month
        dts.day = day
        dts.hour = hour
        dts.min = minute
        dts.sec = second
        dts.us = us
        dts.ps = ns * 1000

        iresult[i] = pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &dts)
        try:
            _check_dts_bounds(&dts)
        except ValueError:
            if is_coerce:
                iresult[i] = NPY_NAT
                continue
            raise

    return result


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

#----------------------------------------------------------------------
# Conversion routines


def pydt_to_i8(object pydt):
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

try:
    import pytz
    UTC = pytz.utc
    have_pytz = True
except:
    have_pytz = False


def tz_convert(ndarray[int64_t] vals, object tz1, object tz2):
    cdef:
        ndarray[int64_t] utc_dates, tt, result, trans, deltas
        Py_ssize_t i, j, pos, n = len(vals)
        ndarray[Py_ssize_t] posn
        int64_t v, offset, delta
        pandas_datetimestruct dts

    if not have_pytz:
        import pytz

    if len(vals) == 0:
        return np.array([], dtype=np.int64)

    # Convert to UTC
    if _get_zone(tz1) != 'UTC':
        utc_dates = np.empty(n, dtype=np.int64)
        if _is_tzlocal(tz1):
            for i in range(n):
                v = vals[i]
                if v == NPY_NAT:
                    utc_dates[i] = NPY_NAT
                else:
                    pandas_datetime_to_datetimestruct(v, PANDAS_FR_ns, &dts)
                    dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                                  dts.min, dts.sec, dts.us, tz1)
                    delta = (int(total_seconds(_get_utcoffset(tz1, dt)))
                             * 1000000000)
                    utc_dates[i] = v - delta
        else:
            trans, deltas, typ = _get_dst_info(tz1)

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

    if _get_zone(tz2) == 'UTC':
        return utc_dates

    result = np.zeros(n, dtype=np.int64)
    if _is_tzlocal(tz2):
        for i in range(n):
            v = utc_dates[i]
            if v == NPY_NAT:
                result[i] = NPY_NAT
            else:
                pandas_datetime_to_datetimestruct(v, PANDAS_FR_ns, &dts)
                dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                              dts.min, dts.sec, dts.us, tz2)
                delta = int(total_seconds(
                    _get_utcoffset(tz2, dt))) * 1000000000
                result[i] = v + delta
        return result

    # Convert UTC to other timezone
    trans, deltas, typ = _get_dst_info(tz2)

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


def tz_convert_single(int64_t val, object tz1, object tz2):
    cdef:
        ndarray[int64_t] trans, deltas
        Py_ssize_t pos
        int64_t v, offset, utc_date
        pandas_datetimestruct dts

    if not have_pytz:
        import pytz

    if val == NPY_NAT:
        return val

    # Convert to UTC
    if _is_tzlocal(tz1):
        pandas_datetime_to_datetimestruct(val, PANDAS_FR_ns, &dts)
        dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                      dts.min, dts.sec, dts.us, tz1)
        delta = int(total_seconds(_get_utcoffset(tz1, dt))) * 1000000000
        utc_date = val - delta
    elif _get_zone(tz1) != 'UTC':
        trans, deltas, typ = _get_dst_info(tz1)
        pos = trans.searchsorted(val, side='right') - 1
        if pos < 0:
            raise ValueError('First time before start of DST info')
        offset = deltas[pos]
        utc_date = val - offset
    else:
        utc_date = val

    if _get_zone(tz2) == 'UTC':
        return utc_date
    if _is_tzlocal(tz2):
        pandas_datetime_to_datetimestruct(val, PANDAS_FR_ns, &dts)
        dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                      dts.min, dts.sec, dts.us, tz2)
        delta = int(total_seconds(_get_utcoffset(tz2, dt))) * 1000000000
        return utc_date + delta

    # Convert UTC to other timezone
    trans, deltas, typ = _get_dst_info(tz2)

    pos = trans.searchsorted(utc_date, side='right') - 1
    if pos < 0:
        raise ValueError('First time before start of DST info')

    offset = deltas[pos]
    return utc_date + offset

# Timezone data caches, key is the pytz string or dateutil file name.
dst_cache = {}

cdef inline bint _treat_tz_as_pytz(object tz):
    return hasattr(tz, '_utc_transition_times') and hasattr(
        tz, '_transition_info')

cdef inline bint _treat_tz_as_dateutil(object tz):
    return hasattr(tz, '_trans_list') and hasattr(tz, '_trans_idx')


def _p_tz_cache_key(tz):
    """ Python interface for cache function to facilitate testing."""
    return _tz_cache_key(tz)


cdef inline object _tz_cache_key(object tz):
    """
    Return the key in the cache for the timezone info object or None
    if unknown.

    The key is currently the tz string for pytz timezones, the filename for
    dateutil timezones.

    Notes
    =====
    This cannot just be the hash of a timezone object. Unfortunately, the
    hashes of two dateutil tz objects which represent the same timezone are
    not equal (even though the tz objects will compare equal and represent
    the same tz file). Also, pytz objects are not always hashable so we use
    str(tz) instead.
    """
    if isinstance(tz, _pytz_BaseTzInfo):
        return tz.zone
    elif isinstance(tz, _dateutil_tzfile):
        if '.tar.gz' in tz._filename:
            raise ValueError('Bad tz filename. Dateutil on python 3 on '
                             'windows has a bug which causes tzfile._filename '
                             'to be the same for all timezone files. Please '
                             'construct dateutil timezones implicitly by '
                             'passing a string like "dateutil/Europe/London" '
                             'when you construct your pandas objects instead '
                             'of passing a timezone object. See '
                             'https://github.com/pydata/pandas/pull/7362')
        return 'dateutil' + tz._filename
    else:
        return None


cdef object _get_dst_info(object tz):
    """
    return a tuple of :
      (UTC times of DST transitions,
       UTC offsets in microseconds corresponding to DST transitions,
       string of type of transitions)

    """
    cache_key = _tz_cache_key(tz)
    if cache_key is None:
        num = int(total_seconds(_get_utcoffset(tz, None))) * 1000000000
        return (np.array([NPY_NAT + 1], dtype=np.int64),
                np.array([num], dtype=np.int64),
                None)

    if cache_key not in dst_cache:
        if _treat_tz_as_pytz(tz):
            trans = np.array(tz._utc_transition_times, dtype='M8[ns]')
            trans = trans.view('i8')
            try:
                if tz._utc_transition_times[0].year == 1:
                    trans[0] = NPY_NAT + 1
            except Exception:
                pass
            deltas = _unbox_utcoffsets(tz._transition_info)
            typ = 'pytz'

        elif _treat_tz_as_dateutil(tz):
            if len(tz._trans_list):
                # get utc trans times
                trans_list = _get_utc_trans_times_from_dateutil_tz(tz)
                trans = np.hstack([
                    np.array([0], dtype='M8[s]'), # place holder for first item
                    np.array(trans_list, dtype='M8[s]')]).astype(
                    'M8[ns]')  # all trans listed
                trans = trans.view('i8')
                trans[0] = NPY_NAT + 1

                # deltas
                deltas = np.array([v.offset for v in (
                    tz._ttinfo_before,) + tz._trans_idx], dtype='i8')
                deltas *= 1000000000
                typ = 'dateutil'

            elif _is_fixed_offset(tz):
                trans = np.array([NPY_NAT + 1], dtype=np.int64)
                deltas = np.array([tz._ttinfo_std.offset],
                                  dtype='i8') * 1000000000
                typ = 'fixed'
            else:
                trans = np.array([], dtype='M8[ns]')
                deltas = np.array([], dtype='i8')
                typ = None

        else:
            # static tzinfo
            trans = np.array([NPY_NAT + 1], dtype=np.int64)
            num = int(total_seconds(_get_utcoffset(tz, None))) * 1000000000
            deltas = np.array([num], dtype=np.int64)
            typ = 'static'

        dst_cache[cache_key] = (trans, deltas, typ)

    return dst_cache[cache_key]

cdef object _get_utc_trans_times_from_dateutil_tz(object tz):
    """
    Transition times in dateutil timezones are stored in local non-dst
    time.  This code converts them to UTC. It's the reverse of the code
    in dateutil.tz.tzfile.__init__.
    """
    new_trans = list(tz._trans_list)
    last_std_offset = 0
    for i, (trans, tti) in enumerate(zip(tz._trans_list, tz._trans_idx)):
        if not tti.isdst:
            last_std_offset = tti.offset
        new_trans[i] = trans - last_std_offset
    return new_trans


def tot_seconds(td):
    return total_seconds(td)

cpdef ndarray _unbox_utcoffsets(object transinfo):
    cdef:
        Py_ssize_t i, sz
        ndarray[int64_t] arr

    sz = len(transinfo)
    arr = np.empty(sz, dtype='i8')

    for i in range(sz):
        arr[i] = int(total_seconds(transinfo[i][0])) * 1000000000

    return arr


@cython.boundscheck(False)
@cython.wraparound(False)
def tz_localize_to_utc(ndarray[int64_t] vals, object tz, object ambiguous=None,
                       object errors='raise'):
    """
    Localize tzinfo-naive DateRange to given time zone (using pytz). If
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

    if not have_pytz:
        raise Exception("Could not find pytz module")

    if tz == UTC or tz is None:
        return vals

    result = np.empty(n, dtype=np.int64)

    if _is_tzlocal(tz):
        for i in range(n):
            v = vals[i]
            pandas_datetime_to_datetimestruct(v, PANDAS_FR_ns, &dts)
            dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                          dts.min, dts.sec, dts.us, tz)
            delta = int(total_seconds(_get_utcoffset(tz, dt))) * 1000000000
            result[i] = v - delta
        return result

    if isinstance(ambiguous, string_types):
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

    trans, deltas, typ = _get_dst_info(tz)

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

def build_field_sarray(ndarray[int64_t] dtindex):
    """
    Datetime as int64 representation to a structured array of fields
    """
    cdef:
        Py_ssize_t i, count = 0
        int isleap
        pandas_datetimestruct dts
        ndarray[int32_t] years, months, days, hours, minutes, seconds, mus

    count = len(dtindex)

    sa_dtype = [('Y', 'i4'), # year
                ('M', 'i4'), # month
                ('D', 'i4'), # day
                ('h', 'i4'), # hour
                ('m', 'i4'), # min
                ('s', 'i4'), # second
                ('u', 'i4')] # microsecond

    out = np.empty(count, dtype=sa_dtype)

    years = out['Y']
    months = out['M']
    days = out['D']
    hours = out['h']
    minutes = out['m']
    seconds = out['s']
    mus = out['u']

    for i in range(count):
        pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
        years[i] = dts.year
        months[i] = dts.month
        days[i] = dts.day
        hours[i] = dts.hour
        minutes[i] = dts.min
        seconds[i] = dts.sec
        mus[i] = dts.us

    return out


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


@cython.wraparound(False)
@cython.boundscheck(False)
def get_date_field(ndarray[int64_t] dtindex, object field):
    """
    Given a int64-based datetime index, extract the year, month, etc.,
    field and return an array of these values.
    """
    cdef:
        _TSObject ts
        Py_ssize_t i, count = 0
        ndarray[int32_t] out
        ndarray[int32_t, ndim=2] _month_offset
        int isleap, isleap_prev
        pandas_datetimestruct dts
        int mo_off, doy, dow, woy

    _month_offset = np.array(
        [[ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 ],
         [ 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 ]],
        dtype=np.int32 )

    count = len(dtindex)
    out = np.empty(count, dtype='i4')

    if field == 'Y':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                out[i] = dts.year
        return out

    elif field == 'M':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                out[i] = dts.month
        return out

    elif field == 'D':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                out[i] = dts.day
        return out

    elif field == 'h':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                out[i] = dts.hour
        return out

    elif field == 'm':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                out[i] = dts.min
        return out

    elif field == 's':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                out[i] = dts.sec
        return out

    elif field == 'us':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                out[i] = dts.us
        return out

    elif field == 'ns':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                out[i] = dts.ps / 1000
        return out
    elif field == 'doy':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                isleap = is_leapyear(dts.year)
                out[i] = _month_offset[isleap, dts.month -1] + dts.day
        return out

    elif field == 'dow':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                out[i] = dayofweek(dts.year, dts.month, dts.day)
        return out

    elif field == 'woy':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                isleap = is_leapyear(dts.year)
                isleap_prev = is_leapyear(dts.year - 1)
                mo_off = _month_offset[isleap, dts.month - 1]
                doy = mo_off + dts.day
                dow = dayofweek(dts.year, dts.month, dts.day)

                #estimate
                woy = (doy - 1) - dow + 3
                if woy >= 0:
                    woy = woy / 7 + 1

                # verify
                if woy < 0:
                    if (woy > -2) or (woy == -2 and isleap_prev):
                        woy = 53
                    else:
                        woy = 52
                elif woy == 53:
                    if 31 - dts.day + dow < 3:
                        woy = 1

                out[i] = woy
        return out

    elif field == 'q':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                out[i] = dts.month
                out[i] = ((out[i] - 1) / 3) + 1
        return out

    elif field == 'dim':
        with nogil:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                out[i] = days_in_month(dts)
        return out
    elif field == 'is_leap_year':
        return _isleapyear_arr(get_date_field(dtindex, 'Y'))

    raise ValueError("Field %s not supported" % field)


@cython.wraparound(False)
def get_start_end_field(ndarray[int64_t] dtindex, object field,
                        object freqstr=None, int month_kw=12):
    """
    Given an int64-based datetime index return array of indicators
    of whether timestamps are at the start/end of the month/quarter/year
    (defined by frequency).
    """
    cdef:
        _TSObject ts
        Py_ssize_t i
        int count = 0
        bint is_business = 0
        int end_month = 12
        int start_month = 1
        ndarray[int8_t] out
        ndarray[int32_t, ndim=2] _month_offset
        bint isleap
        pandas_datetimestruct dts
        int mo_off, dom, doy, dow, ldom

    _month_offset = np.array(
        [[ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 ],
         [ 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 ]],
        dtype=np.int32 )

    count = len(dtindex)
    out = np.zeros(count, dtype='int8')

    if freqstr:
        if freqstr == 'C':
            raise ValueError(
                "Custom business days is not supported by %s" % field)
        is_business = freqstr[0] == 'B'

        # YearBegin(), BYearBegin() use month = starting month of year.
        # QuarterBegin(), BQuarterBegin() use startingMonth = starting
        # month of year. Other offests use month, startingMonth as ending
        # month of year.

        if (freqstr[0:2] in ['MS', 'QS', 'AS']) or (
                freqstr[1:3] in ['MS', 'QS', 'AS']):
            end_month = 12 if month_kw == 1 else month_kw - 1
            start_month = month_kw
        else:
            end_month = month_kw
            start_month = (end_month % 12) + 1
    else:
        end_month = 12
        start_month = 1

    if field == 'is_month_start':
        if is_business:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                ts = convert_to_tsobject(dtindex[i], None, None, 0, 0)
                dom = dts.day
                dow = ts_dayofweek(ts)

                if (dom == 1 and dow < 5) or (dom <= 3 and dow == 0):
                    out[i] = 1
            return out.view(bool)
        else:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                dom = dts.day

                if dom == 1:
                    out[i] = 1
            return out.view(bool)

    elif field == 'is_month_end':
        if is_business:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                ts = convert_to_tsobject(dtindex[i], None, None, 0, 0)
                isleap = is_leapyear(dts.year)
                mo_off = _month_offset[isleap, dts.month - 1]
                dom = dts.day
                doy = mo_off + dom
                ldom = _month_offset[isleap, dts.month]
                dow = ts_dayofweek(ts)

                if (ldom == doy and dow < 5) or (
                        dow == 4 and (ldom - doy <= 2)):
                    out[i] = 1
            return out.view(bool)
        else:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                isleap = is_leapyear(dts.year)
                mo_off = _month_offset[isleap, dts.month - 1]
                dom = dts.day
                doy = mo_off + dom
                ldom = _month_offset[isleap, dts.month]

                if ldom == doy:
                    out[i] = 1
            return out.view(bool)

    elif field == 'is_quarter_start':
        if is_business:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                ts = convert_to_tsobject(dtindex[i], None, None, 0, 0)
                dom = dts.day
                dow = ts_dayofweek(ts)

                if ((dts.month - start_month) % 3 == 0) and (
                        (dom == 1 and dow < 5) or (dom <= 3 and dow == 0)):
                    out[i] = 1
            return out.view(bool)
        else:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                dom = dts.day

                if ((dts.month - start_month) % 3 == 0) and dom == 1:
                    out[i] = 1
            return out.view(bool)

    elif field == 'is_quarter_end':
        if is_business:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                ts = convert_to_tsobject(dtindex[i], None, None, 0, 0)
                isleap = is_leapyear(dts.year)
                mo_off = _month_offset[isleap, dts.month - 1]
                dom = dts.day
                doy = mo_off + dom
                ldom = _month_offset[isleap, dts.month]
                dow = ts_dayofweek(ts)

                if ((dts.month - end_month) % 3 == 0) and (
                        (ldom == doy and dow < 5) or (
                            dow == 4 and (ldom - doy <= 2))):
                    out[i] = 1
            return out.view(bool)
        else:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                isleap = is_leapyear(dts.year)
                mo_off = _month_offset[isleap, dts.month - 1]
                dom = dts.day
                doy = mo_off + dom
                ldom = _month_offset[isleap, dts.month]

                if ((dts.month - end_month) % 3 == 0) and (ldom == doy):
                    out[i] = 1
            return out.view(bool)

    elif field == 'is_year_start':
        if is_business:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                ts = convert_to_tsobject(dtindex[i], None, None, 0, 0)
                dom = dts.day
                dow = ts_dayofweek(ts)

                if (dts.month == start_month) and (
                        (dom == 1 and dow < 5) or (dom <= 3 and dow == 0)):
                    out[i] = 1
            return out.view(bool)
        else:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                dom = dts.day

                if (dts.month == start_month) and dom == 1:
                    out[i] = 1
            return out.view(bool)

    elif field == 'is_year_end':
        if is_business:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                ts = convert_to_tsobject(dtindex[i], None, None, 0, 0)
                isleap = is_leapyear(dts.year)
                dom = dts.day
                mo_off = _month_offset[isleap, dts.month - 1]
                doy = mo_off + dom
                dow = ts_dayofweek(ts)
                ldom = _month_offset[isleap, dts.month]

                if (dts.month == end_month) and (
                        (ldom == doy and dow < 5) or (
                            dow == 4 and (ldom - doy <= 2))):
                    out[i] = 1
            return out.view(bool)
        else:
            for i in range(count):
                if dtindex[i] == NPY_NAT: out[i] = -1; continue

                pandas_datetime_to_datetimestruct(
                    dtindex[i], PANDAS_FR_ns, &dts)
                ts = convert_to_tsobject(dtindex[i], None, None, 0, 0)
                isleap = is_leapyear(dts.year)
                mo_off = _month_offset[isleap, dts.month - 1]
                dom = dts.day
                doy = mo_off + dom
                ldom = _month_offset[isleap, dts.month]

                if (dts.month == end_month) and (ldom == doy):
                    out[i] = 1
            return out.view(bool)

    raise ValueError("Field %s not supported" % field)


@cython.wraparound(False)
@cython.boundscheck(False)
def get_date_name_field(ndarray[int64_t] dtindex, object field):
    """
    Given a int64-based datetime index, return array of strings of date
    name based on requested field (e.g. weekday_name)
    """
    cdef:
        _TSObject ts
        Py_ssize_t i, count = 0
        ndarray[object] out
        pandas_datetimestruct dts
        int dow

    _dayname = np.array(
        ['Monday', 'Tuesday', 'Wednesday', 'Thursday',
            'Friday', 'Saturday', 'Sunday'],
        dtype=np.object_ )

    count = len(dtindex)
    out = np.empty(count, dtype=object)

    if field == 'weekday_name':
        for i in range(count):
            if dtindex[i] == NPY_NAT: out[i] = np.nan; continue

            pandas_datetime_to_datetimestruct(dtindex[i], PANDAS_FR_ns, &dts)
            dow = dayofweek(dts.year, dts.month, dts.day)
            out[i] = _dayname[dow]
        return out

    raise ValueError("Field %s not supported" % field)


cdef inline int m8_weekday(int64_t val):
    ts = convert_to_tsobject(val, None, None, 0, 0)
    return ts_dayofweek(ts)

cdef int64_t DAY_NS = 86400000000000LL


@cython.wraparound(False)
@cython.boundscheck(False)
def date_normalize(ndarray[int64_t] stamps, tz=None):
    cdef:
        Py_ssize_t i, n = len(stamps)
        pandas_datetimestruct dts
        _TSObject tso
        ndarray[int64_t] result = np.empty(n, dtype=np.int64)

    if tz is not None:
        tso = _TSObject()
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

    if _is_utc(tz):
        with nogil:
            for i in range(n):
                if stamps[i] == NPY_NAT:
                    result[i] = NPY_NAT
                    continue
                pandas_datetime_to_datetimestruct(
                    stamps[i], PANDAS_FR_ns, &dts)
                result[i] = _normalized_stamp(&dts)
    elif _is_tzlocal(tz):
        for i in range(n):
            if stamps[i] == NPY_NAT:
                result[i] = NPY_NAT
                continue
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns, &dts)
            dt = datetime(dts.year, dts.month, dts.day, dts.hour,
                          dts.min, dts.sec, dts.us, tz)
            delta = int(total_seconds(_get_utcoffset(tz, dt))) * 1000000000
            pandas_datetime_to_datetimestruct(stamps[i] + delta,
                                              PANDAS_FR_ns, &dts)
            result[i] = _normalized_stamp(&dts)
    else:
        # Adjust datetime64 timestamp, recompute datetimestruct
        trans, deltas, typ = _get_dst_info(tz)

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
        pandas_datetimestruct dts

    if tz is None or _is_utc(tz):
        for i in range(n):
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns, &dts)
            if (dts.hour + dts.min + dts.sec + dts.us) > 0:
                return False
    elif _is_tzlocal(tz):
        for i in range(n):
            pandas_datetime_to_datetimestruct(stamps[i], PANDAS_FR_ns, &dts)
            dt = datetime(dts.year, dts.month, dts.day, dts.hour, dts.min,
                          dts.sec, dts.us, tz)
            dt = dt + tz.utcoffset(dt)
            if (dt.hour + dt.minute + dt.second + dt.microsecond) > 0:
                return False
    else:
        trans, deltas, typ = _get_dst_info(tz)

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


cpdef _isleapyear_arr(ndarray years):
    cdef:
        ndarray[int8_t] out

    # to make NaT result as False
    out = np.zeros(len(years), dtype='int8')
    out[np.logical_or(years % 400 == 0,
                      np.logical_and(years % 4 == 0,
                                     years % 100 > 0))] = 1
    return out.view(bool)


def monthrange(int64_t year, int64_t month):
    cdef:
        int64_t days
        int64_t day_of_week

    if month < 1 or month > 12:
        raise ValueError("bad month number 0; must be 1-12")

    days = days_per_month_table[is_leapyear(year)][month -1]

    return (dayofweek(year, month, 1), days)

cdef inline int64_t ts_dayofweek(_TSObject ts):
    return dayofweek(ts.dts.year, ts.dts.month, ts.dts.day)

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

#----------------------------------------------------------------------
# Don't even ask

"""Strptime-related classes and functions.

CLASSES:
    LocaleTime -- Discovers and stores locale-specific time information
    TimeRE -- Creates regexes for pattern matching a string of text containing
                time information

FUNCTIONS:
    _getlang -- Figure out what language is being used for the locale
    strptime -- Calculates the time struct represented by the passed-in string

"""
import time
import locale
import calendar
from re import compile as re_compile
from re import IGNORECASE
from re import escape as re_escape
from datetime import date as datetime_date

# Python 2 vs Python 3
try:
    from thread import allocate_lock as _thread_allocate_lock
except:
    try:
        from _thread import allocate_lock as _thread_allocate_lock
    except:
        try:
            from dummy_thread import allocate_lock as _thread_allocate_lock
        except:
            from _dummy_thread import allocate_lock as _thread_allocate_lock

__all__ = []


def _getlang():
    # Figure out what the current language is set to.
    return locale.getlocale(locale.LC_TIME)


class LocaleTime(object):
    """Stores and handles locale-specific information related to time.

    ATTRIBUTES:
        f_weekday -- full weekday names (7-item list)
        a_weekday -- abbreviated weekday names (7-item list)
        f_month -- full month names (13-item list; dummy value in [0], which
                    is added by code)
        a_month -- abbreviated month names (13-item list, dummy value in
                    [0], which is added by code)
        am_pm -- AM/PM representation (2-item list)
        LC_date_time -- format string for date/time representation (string)
        LC_date -- format string for date representation (string)
        LC_time -- format string for time representation (string)
        timezone -- daylight- and non-daylight-savings timezone representation
                    (2-item list of sets)
        lang -- Language used by instance (2-item tuple)
    """

    def __init__(self):
        """Set all attributes.

        Order of methods called matters for dependency reasons.

        The locale language is set at the offset and then checked again before
        exiting.  This is to make sure that the attributes were not set with a
        mix of information from more than one locale.  This would most likely
        happen when using threads where one thread calls a locale-dependent
        function while another thread changes the locale while the function in
        the other thread is still running.  Proper coding would call for
        locks to prevent changing the locale while locale-dependent code is
        running.  The check here is done in case someone does not think about
        doing this.

        Only other possible issue is if someone changed the timezone and did
        not call tz.tzset .  That is an issue for the programmer, though,
        since changing the timezone is worthless without that call.

        """
        self.lang = _getlang()
        self.__calc_weekday()
        self.__calc_month()
        self.__calc_am_pm()
        self.__calc_timezone()
        self.__calc_date_time()
        if _getlang() != self.lang:
            raise ValueError("locale changed during initialization")

    def __pad(self, seq, front):
        # Add '' to seq to either the front (is True), else the back.
        seq = list(seq)
        if front:
            seq.insert(0, '')
        else:
            seq.append('')
        return seq

    def __calc_weekday(self):
        # Set self.a_weekday and self.f_weekday using the calendar
        # module.
        a_weekday = [calendar.day_abbr[i].lower() for i in range(7)]
        f_weekday = [calendar.day_name[i].lower() for i in range(7)]
        self.a_weekday = a_weekday
        self.f_weekday = f_weekday

    def __calc_month(self):
        # Set self.f_month and self.a_month using the calendar module.
        a_month = [calendar.month_abbr[i].lower() for i in range(13)]
        f_month = [calendar.month_name[i].lower() for i in range(13)]
        self.a_month = a_month
        self.f_month = f_month

    def __calc_am_pm(self):
        # Set self.am_pm by using time.strftime().

        # The magic date (1999,3,17,hour,44,55,2,76,0) is not really that
        # magical; just happened to have used it everywhere else where a
        # static date was needed.
        am_pm = []
        for hour in (01, 22):
            time_tuple = time.struct_time(
                (1999, 3, 17, hour, 44, 55, 2, 76, 0))
            am_pm.append(time.strftime("%p", time_tuple).lower())
        self.am_pm = am_pm

    def __calc_date_time(self):
        # Set self.date_time, self.date, & self.time by using
        # time.strftime().

        # Use (1999,3,17,22,44,55,2,76,0) for magic date because the amount of
        # overloaded numbers is minimized.  The order in which searches for
        # values within the format string is very important; it eliminates
        # possible ambiguity for what something represents.
        time_tuple = time.struct_time((1999, 3, 17, 22, 44, 55, 2, 76, 0))
        date_time = [None, None, None]
        date_time[0] = time.strftime("%c", time_tuple).lower()
        date_time[1] = time.strftime("%x", time_tuple).lower()
        date_time[2] = time.strftime("%X", time_tuple).lower()
        replacement_pairs = [('%', '%%'), (self.f_weekday[2], '%A'),
                             (self.f_month[3],
                              '%B'), (self.a_weekday[2], '%a'),
                             (self.a_month[3], '%b'), (self.am_pm[1], '%p'),
                             ('1999', '%Y'), ('99', '%y'), ('22', '%H'),
                             ('44', '%M'), ('55', '%S'), ('76', '%j'),
                             ('17', '%d'), ('03', '%m'), ('3', '%m'),
                             # '3' needed for when no leading zero.
                             ('2', '%w'), ('10', '%I')]
        replacement_pairs.extend([(tz, "%Z") for tz_values in self.timezone
                                                for tz in tz_values])
        for offset, directive in ((0, '%c'), (1, '%x'), (2, '%X')):
            current_format = date_time[offset]
            for old, new in replacement_pairs:
                # Must deal with possible lack of locale info
                # manifesting itself as the empty string (e.g., Swedish's
                # lack of AM/PM info) or a platform returning a tuple of empty
                # strings (e.g., MacOS 9 having timezone as ('','')).
                if old:
                    current_format = current_format.replace(old, new)
            # If %W is used, then Sunday, 2005-01-03 will fall on week 0 since
            # 2005-01-03 occurs before the first Monday of the year.  Otherwise
            # %U is used.
            time_tuple = time.struct_time((1999, 1, 3, 1, 1, 1, 6, 3, 0))
            if '00' in time.strftime(directive, time_tuple):
                U_W = '%W'
            else:
                U_W = '%U'
            date_time[offset] = current_format.replace('11', U_W)
        self.LC_date_time = date_time[0]
        self.LC_date = date_time[1]
        self.LC_time = date_time[2]

    def __calc_timezone(self):
        # Set self.timezone by using time.tzname.
        # Do not worry about possibility of time.tzname[0] == timetzname[1]
        # and time.daylight; handle that in strptime .
        try:
            time.tzset()
        except AttributeError:
            pass
        no_saving = frozenset(["utc", "gmt", time.tzname[0].lower()])
        if time.daylight:
            has_saving = frozenset([time.tzname[1].lower()])
        else:
            has_saving = frozenset()
        self.timezone = (no_saving, has_saving)


class TimeRE(dict):
    """Handle conversion from format directives to regexes."""

    def __init__(self, locale_time=None):
        """Create keys/values.

        Order of execution is important for dependency reasons.

        """
        if locale_time:
            self.locale_time = locale_time
        else:
            self.locale_time = LocaleTime()
        base = super(TimeRE, self)
        base.__init__({
            # The " \d" part of the regex is to make %c from ANSI C work
            'd': r"(?P<d>3[0-1]|[1-2]\d|0[1-9]|[1-9]| [1-9])",
            'f': r"(?P<f>[0-9]{1,9})",
            'H': r"(?P<H>2[0-3]|[0-1]\d|\d)",
            'I': r"(?P<I>1[0-2]|0[1-9]|[1-9])",
            'j': (r"(?P<j>36[0-6]|3[0-5]\d|[1-2]\d\d|0[1-9]\d|00[1-9]|"
                  r"[1-9]\d|0[1-9]|[1-9])"),
            'm': r"(?P<m>1[0-2]|0[1-9]|[1-9])",
            'M': r"(?P<M>[0-5]\d|\d)",
            'S': r"(?P<S>6[0-1]|[0-5]\d|\d)",
            'U': r"(?P<U>5[0-3]|[0-4]\d|\d)",
            'w': r"(?P<w>[0-6])",
            # W is set below by using 'U'
            'y': r"(?P<y>\d\d)",
            #XXX: Does 'Y' need to worry about having less or more than
            #     4 digits?
            'Y': r"(?P<Y>\d\d\d\d)",
            'A': self.__seqToRE(self.locale_time.f_weekday, 'A'),
            'a': self.__seqToRE(self.locale_time.a_weekday, 'a'),
            'B': self.__seqToRE(self.locale_time.f_month[1:], 'B'),
            'b': self.__seqToRE(self.locale_time.a_month[1:], 'b'),
            'p': self.__seqToRE(self.locale_time.am_pm, 'p'),
            'Z': self.__seqToRE((tz for tz_names in self.locale_time.timezone
                                        for tz in tz_names),
                                'Z'),
            '%': '%'})
        base.__setitem__('W', base.__getitem__('U').replace('U', 'W'))
        base.__setitem__('c', self.pattern(self.locale_time.LC_date_time))
        base.__setitem__('x', self.pattern(self.locale_time.LC_date))
        base.__setitem__('X', self.pattern(self.locale_time.LC_time))

    def __seqToRE(self, to_convert, directive):
        """Convert a list to a regex string for matching a directive.

        Want possible matching values to be from longest to shortest.  This
        prevents the possibility of a match occuring for a value that also
        a substring of a larger value that should have matched (e.g., 'abc'
        matching when 'abcdef' should have been the match).

        """
        to_convert = sorted(to_convert, key=len, reverse=True)
        for value in to_convert:
            if value != '':
                break
        else:
            return ''
        regex = '|'.join(re_escape(stuff) for stuff in to_convert)
        regex = '(?P<%s>%s' % (directive, regex)
        return '%s)' % regex

    def pattern(self, format):
        """Return regex pattern for the format string.

        Need to make sure that any characters that might be interpreted as
        regex syntax are escaped.

        """
        processed_format = ''
        # The sub() call escapes all characters that might be misconstrued
        # as regex syntax.  Cannot use re.escape since we have to deal with
        # format directives (%m, etc.).
        regex_chars = re_compile(r"([\\.^$*+?\(\){}\[\]|])")
        format = regex_chars.sub(r"\\\1", format)
        whitespace_replacement = re_compile(r'\s+')
        format = whitespace_replacement.sub(r'\\s+', format)
        while '%' in format:
            directive_index = format.index('%') +1
            processed_format = "%s%s%s" % (processed_format,
                                           format[:directive_index -1],
                                           self[format[directive_index]])
            format = format[directive_index +1:]
        return "%s%s" % (processed_format, format)

    def compile(self, format):
        """Return a compiled re object for the format string."""
        return re_compile(self.pattern(format), IGNORECASE)

_cache_lock = _thread_allocate_lock()
# DO NOT modify _TimeRE_cache or _regex_cache without acquiring the cache lock
# first!
_TimeRE_cache = TimeRE()
_CACHE_MAX_SIZE = 5 # Max number of regexes stored in _regex_cache
_regex_cache = {}

cdef _calc_julian_from_U_or_W(int year, int week_of_year,
                              int day_of_week, int week_starts_Mon):
    """Calculate the Julian day based on the year, week of the year, and day of
    the week, with week_start_day representing whether the week of the year
    assumes the week starts on Sunday or Monday (6 or 0)."""

    cdef:
        int first_weekday,  week_0_length, days_to_week

    first_weekday = datetime_date(year, 1, 1).weekday()
    # If we are dealing with the %U directive (week starts on Sunday), it's
    # easier to just shift the view to Sunday being the first day of the
    # week.
    if not week_starts_Mon:
        first_weekday = (first_weekday + 1) % 7
        day_of_week = (day_of_week + 1) % 7
    # Need to watch out for a week 0 (when the first day of the year is not
    # the same as that specified by %U or %W).
    week_0_length = (7 - first_weekday) % 7
    if week_of_year == 0:
        return 1 + day_of_week - first_weekday
    else:
        days_to_week = week_0_length + (7 * (week_of_year - 1))
        return 1 + days_to_week + day_of_week

# def _strptime_time(data_string, format="%a %b %d %H:%M:%S %Y"):
#     return _strptime(data_string, format)[0]
