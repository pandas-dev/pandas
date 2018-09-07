# -*- coding: utf-8 -*-
import warnings

from cpython cimport (PyObject_RichCompareBool, PyObject_RichCompare,
                      Py_GT, Py_GE, Py_EQ, Py_NE, Py_LT, Py_LE)

import numpy as np
cimport numpy as cnp
from numpy cimport int64_t, int32_t, int8_t
cnp.import_array()

from datetime import time as datetime_time
from cpython.datetime cimport (datetime,
                               PyDateTime_Check, PyDelta_Check, PyTZInfo_Check,
                               PyDateTime_IMPORT)
PyDateTime_IMPORT

from util cimport (is_datetime64_object, is_timedelta64_object,
                   is_integer_object, is_string_object, is_array)

cimport ccalendar
from conversion import tz_localize_to_utc, normalize_i8_timestamps
from conversion cimport (tz_convert_single, _TSObject,
                         convert_to_tsobject, convert_datetime_to_tsobject)
from fields import get_start_end_field, get_date_name_field
from nattype import NaT
from nattype cimport NPY_NAT
from np_datetime import OutOfBoundsDatetime
from np_datetime cimport (reverse_ops, cmp_scalar, check_dts_bounds,
                          npy_datetimestruct, dt64_to_dtstruct)
from offsets cimport to_offset
from timedeltas import Timedelta
from timedeltas cimport delta_to_nanoseconds
from timezones cimport (
    get_timezone, is_utc, maybe_get_tz, treat_tz_as_pytz, tz_compare)

# ----------------------------------------------------------------------
# Constants
_zero_time = datetime_time(0, 0)
_no_input = object()

# ----------------------------------------------------------------------


cdef inline object create_timestamp_from_ts(int64_t value,
                                            npy_datetimestruct dts,
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


def round_ns(values, rounder, freq):
    """
    Applies rounding function at given frequency

    Parameters
    ----------
    values : :obj:`ndarray`
    rounder : function, eg. 'ceil', 'floor', 'round'
    freq : str, obj

    Returns
    -------
    :obj:`ndarray`
    """
    unit = to_offset(freq).nanos

    # GH21262 If the Timestamp is multiple of the freq str
    # don't apply any rounding
    mask = values % unit == 0
    if mask.all():
        return values
    r = values.copy()

    if unit < 1000:
        # for nano rounding, work with the last 6 digits separately
        # due to float precision
        buff = 1000000
        r[~mask] = (buff * (values[~mask] // buff) +
                    unit * (rounder((values[~mask] % buff) *
                            (1 / float(unit)))).astype('i8'))
    else:
        if unit % 1000 != 0:
            msg = 'Precision will be lost using frequency: {}'
            warnings.warn(msg.format(freq))
        # GH19206
        # to deal with round-off when unit is large
        if unit >= 1e9:
            divisor = 10 ** int(np.log10(unit / 1e7))
        else:
            divisor = 10
        r[~mask] = (unit * rounder((values[~mask] *
                    (divisor / float(unit))) / divisor)
                    .astype('i8'))
    return r


# This is PITA. Because we inherit from datetime, which has very specific
# construction requirements, we need to do object instantiation in python
# (see Timestamp class above). This will serve as a C extension type that
# shadows the python class, where we do any heavy lifting.
cdef class _Timestamp(datetime):

    cdef readonly:
        int64_t value, nanosecond
        object freq       # frequency reference
        list _date_attributes

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
            ndim = getattr(other, "ndim", -1)

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
                elif is_array(other):
                    # avoid recursion error GH#15183
                    return PyObject_RichCompare(np.array([self]), other, op)
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
        freq = "" if self.freq is None else ", freq='{0}'".format(self.freqstr)

        return "Timestamp('{stamp}'{tz}{freq})".format(stamp=stamp,
                                                       tz=tz, freq=freq)

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
            # delta --> offsets.Tick
            nanos = delta_to_nanoseconds(other)
            result = Timestamp(self.value + nanos,
                               tz=self.tzinfo, freq=self.freq)
            if getattr(other, 'normalize', False):
                # DateOffset
                result = result.normalize()
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
            if not tz_compare(self.tzinfo, other.tzinfo):
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

    cpdef bint _get_start_end_field(self, str field):
        cdef:
            int64_t val
            dict kwds
            int8_t out[1]
            int month_kw

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

    cpdef _get_date_name_field(self, object field, object locale):
        cdef:
            int64_t val
            object[:] out

        val = self._maybe_convert_value_to_local()
        out = get_date_name_field(np.array([val], dtype=np.int64),
                                  field, locale=locale)
        return out[0]

    @property
    def _repr_base(self):
        return '{date} {time}'.format(date=self._date_repr,
                                      time=self._time_repr)

    @property
    def _date_repr(self):
        # Ideal here would be self.strftime("%Y-%m-%d"), but
        # the datetime strftime() methods require year >= 1900
        return '%d-%.2d-%.2d' % (self.year, self.month, self.day)

    @property
    def _time_repr(self):
        result = '%.2d:%.2d:%.2d' % (self.hour, self.minute, self.second)

        if self.nanosecond != 0:
            result += '.%.9d' % (self.nanosecond + 1000 * self.microsecond)
        elif self.microsecond != 0:
            result += '.%.6d' % self.microsecond

        return result

    @property
    def _short_repr(self):
        # format a Timestamp with only _date_repr if possible
        # otherwise _repr_base
        if (self.hour == 0 and
                self.minute == 0 and
                self.second == 0 and
                self.microsecond == 0 and
                self.nanosecond == 0):
            return self._date_repr
        return self._repr_base

    @property
    def asm8(self):
        return np.datetime64(self.value, 'ns')

    @property
    def resolution(self):
        """
        Return resolution describing the smallest difference between two
        times that can be represented by Timestamp object_state
        """
        # GH#21336, GH#21365
        return Timedelta(nanoseconds=1)

    def timestamp(self):
        """Return POSIX timestamp as float."""
        # py27 compat, see GH#17329
        return round(self.value / 1e9, 6)


# ----------------------------------------------------------------------

# Python front end to C extension type _Timestamp
# This serves as the box for datetime64


class Timestamp(_Timestamp):
    """Pandas replacement for datetime.datetime

    Timestamp is the pandas equivalent of python's Datetime
    and is interchangeable with it in most cases. It's the type used
    for the entries that make up a DatetimeIndex, and other timeseries
    oriented data structures in pandas.

    Parameters
    ----------
    ts_input : datetime-like, str, int, float
        Value to be converted to Timestamp
    freq : str, DateOffset
        Offset which Timestamp will have
    tz : str, pytz.timezone, dateutil.tz.tzfile or None
        Time zone for time which Timestamp will have.
    unit : str
        Unit used for conversion if ts_input is of type int or float. The
        valid values are 'D', 'h', 'm', 's', 'ms', 'us', and 'ns'. For
        example, 's' means seconds and 'ms' means milliseconds.
    year, month, day : int
        .. versionadded:: 0.19.0
    hour, minute, second, microsecond : int, optional, default 0
        .. versionadded:: 0.19.0
    nanosecond : int, optional, default 0
        .. versionadded:: 0.23.0
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
    Using the primary calling convention:

    This converts a datetime-like string
    >>> pd.Timestamp('2017-01-01T12')
    Timestamp('2017-01-01 12:00:00')

    This converts a float representing a Unix epoch in units of seconds
    >>> pd.Timestamp(1513393355.5, unit='s')
    Timestamp('2017-12-16 03:02:35.500000')

    This converts an int representing a Unix-epoch in units of seconds
    and for a particular timezone
    >>> pd.Timestamp(1513393355, unit='s', tz='US/Pacific')
    Timestamp('2017-12-15 19:02:35-0800', tz='US/Pacific')

    Using the other two forms that mimic the API for ``datetime.datetime``:

    >>> pd.Timestamp(2017, 1, 1, 12)
    Timestamp('2017-01-01 12:00:00')

    >>> pd.Timestamp(year=2017, month=1, day=1, hour=12)
    Timestamp('2017-01-01 12:00:00')
    """

    @classmethod
    def fromordinal(cls, ordinal, freq=None, tz=None):
        """
        Timestamp.fromordinal(ordinal, freq=None, tz=None)

        passed an ordinal, translate and convert to a ts
        note: by definition there cannot be any tz info on the ordinal itself

        Parameters
        ----------
        ordinal : int
            date corresponding to a proleptic Gregorian ordinal
        freq : str, DateOffset
            Offset which Timestamp will have
        tz : str, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for time which Timestamp will have.
        """
        return cls(datetime.fromordinal(ordinal),
                   freq=freq, tz=tz)

    @classmethod
    def now(cls, tz=None):
        """
        Timestamp.now(tz=None)

        Returns new Timestamp object representing current time local to
        tz.

        Parameters
        ----------
        tz : str or timezone object, default None
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
        tz : str or timezone object, default None
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
                nanosecond=None, tzinfo=None):
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

        _date_attributes = [year, month, day, hour, minute, second,
                            microsecond, nanosecond]

        if tzinfo is not None:
            if not PyTZInfo_Check(tzinfo):
                # tzinfo must be a datetime.tzinfo object, GH#17690
                raise TypeError('tzinfo must be a datetime.tzinfo object, '
                                'not %s' % type(tzinfo))
            elif tz is not None:
                raise ValueError('Can provide at most one of tz, tzinfo')

        if is_string_object(ts_input):
            # User passed a date string to parse.
            # Check that the user didn't also pass a date attribute kwarg.
            if any(arg is not None for arg in _date_attributes):
                raise ValueError('Cannot pass a date attribute keyword '
                                 'argument when passing a date string')

        elif ts_input is _no_input:
            # User passed keyword arguments.
            if tz is None:
                # Handle the case where the user passes `tz` and not `tzinfo`
                tz = tzinfo
            return Timestamp(datetime(year, month, day, hour or 0,
                                      minute or 0, second or 0,
                                      microsecond or 0, tzinfo),
                             nanosecond=nanosecond, tz=tz)
        elif is_integer_object(freq):
            # User passed positional arguments:
            # Timestamp(year, month, day[, hour[, minute[, second[,
            # microsecond[, nanosecond[, tzinfo]]]]]])
            return Timestamp(datetime(ts_input, freq, tz, unit or 0,
                                      year or 0, month or 0, day or 0,
                                      minute), nanosecond=hour, tz=minute)

        if tzinfo is not None:
            # User passed tzinfo instead of tz; avoid silently ignoring
            tz, tzinfo = tzinfo, None

        ts = convert_to_tsobject(ts_input, tz, unit, 0, 0, nanosecond or 0)

        if ts.value == NPY_NAT:
            return NaT

        if is_string_object(freq):
            freq = to_offset(freq)

        return create_timestamp_from_ts(ts.value, ts.dts, ts.tzinfo, freq)

    def _round(self, freq, rounder):
        if self.tz is not None:
            value = self.tz_localize(None).value
        else:
            value = self.value

        value = np.array([value], dtype=np.int64)

        # Will only ever contain 1 element for timestamp
        r = round_ns(value, rounder, freq)[0]
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

    @tz.setter
    def tz(self, value):
        # GH 3746: Prevent localizing or converting the index by setting tz
        raise AttributeError("Cannot directly set timezone. Use tz_localize() "
                             "or tz_convert() as appropriate")

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

    def day_name(self, locale=None):
        """
        Return the day name of the Timestamp with specified locale.

        Parameters
        ----------
        locale : string, default None (English locale)
            locale determining the language in which to return the day name

        Returns
        -------
        day_name : string

        .. versionadded:: 0.23.0
        """
        return self._get_date_name_field('day_name', locale)

    def month_name(self, locale=None):
        """
        Return the month name of the Timestamp with specified locale.

        Parameters
        ----------
        locale : string, default None (English locale)
            locale determining the language in which to return the month name

        Returns
        -------
        month_name : string

        .. versionadded:: 0.23.0
        """
        return self._get_date_name_field('month_name', locale)

    @property
    def weekday_name(self):
        """
        .. deprecated:: 0.23.0
            Use ``Timestamp.day_name()`` instead
        """
        warnings.warn("`weekday_name` is deprecated and will be removed in a "
                      "future version. Use `day_name` instead",
                      FutureWarning)
        return self.day_name()

    @property
    def dayofyear(self):
        return ccalendar.get_day_of_year(self.year, self.month, self.day)

    @property
    def week(self):
        return ccalendar.get_week_of_year(self.year, self.month, self.day)

    weekofyear = week

    @property
    def quarter(self):
        return ((self.month - 1) // 3) + 1

    @property
    def days_in_month(self):
        return ccalendar.get_days_in_month(self.year, self.month)

    daysinmonth = days_in_month

    @property
    def freqstr(self):
        return getattr(self.freq, 'freqstr', self.freq)

    @property
    def is_month_start(self):
        if self.freq is None:
            # fast-path for non-business frequencies
            return self.day == 1
        return self._get_start_end_field('is_month_start')

    @property
    def is_month_end(self):
        if self.freq is None:
            # fast-path for non-business frequencies
            return self.day == self.days_in_month
        return self._get_start_end_field('is_month_end')

    @property
    def is_quarter_start(self):
        if self.freq is None:
            # fast-path for non-business frequencies
            return self.day == 1 and self.month % 3 == 1
        return self._get_start_end_field('is_quarter_start')

    @property
    def is_quarter_end(self):
        if self.freq is None:
            # fast-path for non-business frequencies
            return (self.month % 3) == 0 and self.day == self.days_in_month
        return self._get_start_end_field('is_quarter_end')

    @property
    def is_year_start(self):
        if self.freq is None:
            # fast-path for non-business frequencies
            return self.day == self.month == 1
        return self._get_start_end_field('is_year_start')

    @property
    def is_year_end(self):
        if self.freq is None:
            # fast-path for non-business frequencies
            return self.month == 12 and self.day == 31
        return self._get_start_end_field('is_year_end')

    @property
    def is_leap_year(self):
        return bool(ccalendar.is_leapyear(self.year))

    def tz_localize(self, tz, ambiguous='raise', nonexistent='raise',
                    errors='raise'):
        """
        Convert naive Timestamp to local time zone, or remove
        timezone from tz-aware Timestamp.

        Parameters
        ----------
        tz : str, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for time which Timestamp will be converted to.
            None will remove timezone holding local time.

        ambiguous : bool, 'NaT', default 'raise'
            - bool contains flags to determine if time is dst or not (note
              that this flag is only applicable for ambiguous fall dst dates)
            - 'NaT' will return NaT for an ambiguous time
            - 'raise' will raise an AmbiguousTimeError for an ambiguous time

        nonexistent : str {'NaT', 'raise'}

            - 'infer' will shift the non-existent time to a real local time
            - 'NaT' will return NaT where there are ambiguous times
            - 'raise' will raise an NonExistentTimeError if there are ambiguous
              times

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
                                       ambiguous=ambiguous,
                                       nonexistent=nonexistent,
                                       errors=errors)[0]
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
        tz : str, pytz.timezone, dateutil.tz.tzfile or None
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
            npy_datetimestruct dts
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
        if _tzinfo is not None and treat_tz_as_pytz(_tzinfo):
            # replacing across a DST boundary may induce a new tzinfo object
            # see GH#18319
            ts_input = _tzinfo.localize(datetime(dts.year, dts.month, dts.day,
                                                 dts.hour, dts.min, dts.sec,
                                                 dts.us))
            _tzinfo = ts_input.tzinfo
        else:
            ts_input = datetime(dts.year, dts.month, dts.day,
                                dts.hour, dts.min, dts.sec, dts.us,
                                tzinfo=_tzinfo)

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
        normalized_value = normalize_i8_timestamps(
            np.array([self.value], dtype='i8'), tz=self.tz)[0]
        return Timestamp(normalized_value).tz_localize(self.tz)

    def __radd__(self, other):
        # __radd__ on cython extension types like _Timestamp is not used, so
        # define it here instead
        return self + other


# Add the min and max fields at the class level
cdef int64_t _NS_UPPER_BOUND = np.iinfo(np.int64).max
# the smallest value we could actually represent is
#   INT64_MIN + 1 == -9223372036854775807
# but to allow overflow free conversion with a microsecond resolution
# use the smallest value with a 0 nanosecond unit (0s in last 3 digits)
cdef int64_t _NS_LOWER_BOUND = -9223372036854775000

# Resolution is in nanoseconds
Timestamp.min = Timestamp(_NS_LOWER_BOUND)
Timestamp.max = Timestamp(_NS_UPPER_BOUND)
