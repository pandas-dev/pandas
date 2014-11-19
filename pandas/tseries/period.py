# pylint: disable=E1101,E1103,W0232
import operator

from datetime import datetime, date, timedelta
import numpy as np
from pandas.core.base import PandasObject

import pandas.tseries.frequencies as frequencies
from pandas.tseries.frequencies import get_freq_code as _gfc
from pandas.tseries.index import DatetimeIndex, Int64Index, Index
from pandas.tseries.base import DatetimeIndexOpsMixin
from pandas.tseries.tools import parse_time_string
import pandas.tseries.offsets as offsets

import pandas.core.common as com
from pandas.core.common import (isnull, _INT64_DTYPE, _maybe_box,
                                _values_from_object, ABCSeries)
from pandas import compat
from pandas.lib import Timestamp, Timedelta
import pandas.lib as lib
import pandas.tslib as tslib
import pandas.algos as _algos
from pandas.compat import zip, u


#---------------
# Period logic

def _period_field_accessor(name, alias):
    def f(self):
        base, mult = _gfc(self.freq)
        return tslib.get_period_field(alias, self.ordinal, base)
    f.__name__ = name
    return property(f)


def _field_accessor(name, alias, docstring=None):
    def f(self):
        base, mult = _gfc(self.freq)
        return tslib.get_period_field_arr(alias, self.values, base)
    f.__name__ = name
    f.__doc__ = docstring
    return property(f)

class Period(PandasObject):
    """
    Represents an period of time

    Parameters
    ----------
    value : Period or compat.string_types, default None
        The time period represented (e.g., '4Q2005')
    freq : str, default None
        e.g., 'B' for businessday. Must be a singular rule-code (e.g. 5T is not
        allowed).
    year : int, default None
    month : int, default 1
    quarter : int, default None
    day : int, default 1
    hour : int, default 0
    minute : int, default 0
    second : int, default 0
    """
    __slots__ = ['freq', 'ordinal']
    _comparables = ['name','freqstr']
    _typ = 'period'

    @classmethod
    def _from_ordinal(cls, ordinal, freq):
        """ fast creation from an ordinal and freq that are already validated! """
        self = object.__new__(cls)
        self.ordinal = ordinal
        self.freq = freq
        return self

    def __init__(self, value=None, freq=None, ordinal=None,
                 year=None, month=1, quarter=None, day=1,
                 hour=0, minute=0, second=0):

        # freq points to a tuple (base, mult);  base is one of the defined
        # periods such as A, Q, etc. Every five minutes would be, e.g.,
        # ('T', 5) but may be passed in as a string like '5T'

        self.freq = None

        # ordinal is the period offset from the gregorian proleptic epoch
        self.ordinal = None

        if ordinal is not None and value is not None:
            raise ValueError(("Only value or ordinal but not both should be "
                              "given but not both"))
        elif ordinal is not None:
            if not com.is_integer(ordinal):
                raise ValueError("Ordinal must be an integer")
            if freq is None:
                raise ValueError('Must supply freq for ordinal value')
            self.ordinal = ordinal

        elif value is None:
            if freq is None:
                raise ValueError("If value is None, freq cannot be None")

            self.ordinal = _ordinal_from_fields(year, month, quarter, day,
                                                hour, minute, second, freq)

        elif isinstance(value, Period):
            other = value
            if freq is None or _gfc(freq) == _gfc(other.freq):
                self.ordinal = other.ordinal
                freq = other.freq
            else:
                converted = other.asfreq(freq)
                self.ordinal = converted.ordinal

        elif com._is_null_datelike_scalar(value) or value in tslib._nat_strings:
            self.ordinal = tslib.iNaT
            if freq is None:
                raise ValueError("If value is NaT, freq cannot be None "
                                 "because it cannot be inferred")

        elif isinstance(value, compat.string_types) or com.is_integer(value):
            if com.is_integer(value):
                value = str(value)
            value = value.upper()

            dt, _, reso = parse_time_string(value, freq)
            if freq is None:
                try:
                    freq = frequencies.Resolution.get_freq(reso)
                except KeyError:
                    raise ValueError("Invalid frequency or could not infer: %s" % reso)

        elif isinstance(value, datetime):
            dt = value
            if freq is None:
                raise ValueError('Must supply freq for datetime value')
        elif isinstance(value, date):
            dt = datetime(year=value.year, month=value.month, day=value.day)
            if freq is None:
                raise ValueError('Must supply freq for datetime value')
        else:
            msg = "Value must be Period, string, integer, or datetime"
            raise ValueError(msg)

        base, mult = _gfc(freq)
        if mult != 1:
            # TODO: Better error message - this is slightly confusing
            raise ValueError('Only mult == 1 supported')

        if self.ordinal is None:
            self.ordinal = tslib.period_ordinal(dt.year, dt.month, dt.day,
                                                dt.hour, dt.minute, dt.second, dt.microsecond, 0,
                                                base)

        self.freq = frequencies._get_freq_str(base)

    def __eq__(self, other):
        if isinstance(other, Period):
            if other.freq != self.freq:
                raise ValueError("Cannot compare non-conforming periods")
            if self.ordinal == tslib.iNaT or other.ordinal == tslib.iNaT:
                return False
            return (self.ordinal == other.ordinal
                    and _gfc(self.freq) == _gfc(other.freq))
        return NotImplemented

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((self.ordinal, self.freq))

    def _add_delta(self, other):
        if isinstance(other, (timedelta, np.timedelta64, offsets.Tick, Timedelta)):
            offset = frequencies.to_offset(self.freq)
            if isinstance(offset, offsets.Tick):
                nanos = tslib._delta_to_nanoseconds(other)
                offset_nanos = tslib._delta_to_nanoseconds(offset)

                if nanos % offset_nanos == 0:
                    if self.ordinal == tslib.iNaT:
                        ordinal = self.ordinal
                    else:
                        ordinal = self.ordinal + (nanos // offset_nanos)
                    return Period(ordinal=ordinal, freq=self.freq)
        elif isinstance(other, offsets.DateOffset):
            freqstr = frequencies.get_standard_freq(other)
            base = frequencies.get_base_alias(freqstr)

            if base == self.freq:
                if self.ordinal == tslib.iNaT:
                    ordinal = self.ordinal
                else:
                    ordinal = self.ordinal + other.n
                return Period(ordinal=ordinal, freq=self.freq)

        raise ValueError("Input has different freq from Period(freq={0})".format(self.freq))

    def __add__(self, other):
        if isinstance(other, (timedelta, np.timedelta64,
                              offsets.Tick, offsets.DateOffset, Timedelta)):
            return self._add_delta(other)
        elif com.is_integer(other):
            if self.ordinal == tslib.iNaT:
                ordinal = self.ordinal
            else:
                ordinal = self.ordinal + other
            return Period(ordinal=ordinal, freq=self.freq)
        else:  # pragma: no cover
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, (timedelta, np.timedelta64,
                              offsets.Tick, offsets.DateOffset, Timedelta)):
            neg_other = -other
            return self + neg_other
        elif com.is_integer(other):
            if self.ordinal == tslib.iNaT:
                ordinal = self.ordinal
            else:
                ordinal = self.ordinal - other
            return Period(ordinal=ordinal, freq=self.freq)
        elif isinstance(other, Period):
            if other.freq != self.freq:
                raise ValueError("Cannot do arithmetic with "
                                 "non-conforming periods")
            if self.ordinal == tslib.iNaT or other.ordinal == tslib.iNaT:
                return Period(ordinal=tslib.iNaT, freq=self.freq)
            return self.ordinal - other.ordinal
        else:  # pragma: no cover
            return NotImplemented

    def _comp_method(func, name):
        def f(self, other):
            if isinstance(other, Period):
                if other.freq != self.freq:
                    raise ValueError("Cannot compare non-conforming periods")
                if self.ordinal == tslib.iNaT or other.ordinal == tslib.iNaT:
                    return False
                return func(self.ordinal, other.ordinal)
            else:
                raise TypeError(other)

        f.__name__ = name
        return f

    __lt__ = _comp_method(operator.lt, '__lt__')
    __le__ = _comp_method(operator.le, '__le__')
    __gt__ = _comp_method(operator.gt, '__gt__')
    __ge__ = _comp_method(operator.ge, '__ge__')

    def asfreq(self, freq, how='E'):
        """
        Convert Period to desired frequency, either at the start or end of the
        interval

        Parameters
        ----------
        freq : string
        how : {'E', 'S', 'end', 'start'}, default 'end'
            Start or end of the timespan

        Returns
        -------
        resampled : Period
        """
        how = _validate_end_alias(how)
        base1, mult1 = _gfc(self.freq)
        base2, mult2 = _gfc(freq)

        if mult2 != 1:
            raise ValueError('Only mult == 1 supported')

        end = how == 'E'
        new_ordinal = tslib.period_asfreq(self.ordinal, base1, base2, end)

        return Period(ordinal=new_ordinal, freq=base2)

    @property
    def start_time(self):
        return self.to_timestamp(how='S')

    @property
    def end_time(self):
        if self.ordinal == tslib.iNaT:
            ordinal = self.ordinal
        else:
            ordinal = (self + 1).start_time.value - 1
        return Timestamp(ordinal)

    def to_timestamp(self, freq=None, how='start', tz=None):
        """
        Return the Timestamp representation of the Period at the target
        frequency at the specified end (how) of the Period

        Parameters
        ----------
        freq : string or DateOffset, default is 'D' if self.freq is week or
               longer and 'S' otherwise
            Target frequency
        how: str, default 'S' (start)
            'S', 'E'. Can be aliased as case insensitive
            'Start', 'Finish', 'Begin', 'End'

        Returns
        -------
        Timestamp
        """
        how = _validate_end_alias(how)

        if freq is None:
            base, mult = _gfc(self.freq)
            freq = frequencies.get_to_timestamp_base(base)

        base, mult = _gfc(freq)
        val = self.asfreq(freq, how)

        dt64 = tslib.period_ordinal_to_dt64(val.ordinal, base)
        return Timestamp(dt64, tz=tz)

    year = _period_field_accessor('year', 0)
    month = _period_field_accessor('month', 3)
    day = _period_field_accessor('day', 4)
    hour = _period_field_accessor('hour', 5)
    minute = _period_field_accessor('minute', 6)
    second = _period_field_accessor('second', 7)
    weekofyear = _period_field_accessor('week', 8)
    week = weekofyear
    dayofweek = _period_field_accessor('dayofweek', 10)
    weekday = dayofweek
    dayofyear = _period_field_accessor('dayofyear', 9)
    quarter = _period_field_accessor('quarter', 2)
    qyear = _period_field_accessor('qyear', 1)

    @classmethod
    def now(cls, freq=None):
        return Period(datetime.now(), freq=freq)

    def __repr__(self):
        base, mult = _gfc(self.freq)
        formatted = tslib.period_format(self.ordinal, base)
        freqstr = frequencies._reverse_period_code_map[base]

        if not compat.PY3:
            encoding = com.get_option("display.encoding")
            formatted = formatted.encode(encoding)

        return "Period('%s', '%s')" % (formatted, freqstr)

    def __unicode__(self):
        """
        Return a string representation for a particular DataFrame

        Invoked by unicode(df) in py2 only. Yields a Unicode String in both
        py2/py3.
        """
        base, mult = _gfc(self.freq)
        formatted = tslib.period_format(self.ordinal, base)
        value = ("%s" % formatted)
        return value

    def strftime(self, fmt):
        """
        Returns the string representation of the :class:`Period`, depending
        on the selected :keyword:`format`. :keyword:`format` must be a string
        containing one or several directives.  The method recognizes the same
        directives as the :func:`time.strftime` function of the standard Python
        distribution, as well as the specific additional directives ``%f``,
        ``%F``, ``%q``. (formatting & docs originally from scikits.timeries)

        +-----------+--------------------------------+-------+
        | Directive | Meaning                        | Notes |
        +===========+================================+=======+
        | ``%a``    | Locale's abbreviated weekday   |       |
        |           | name.                          |       |
        +-----------+--------------------------------+-------+
        | ``%A``    | Locale's full weekday name.    |       |
        +-----------+--------------------------------+-------+
        | ``%b``    | Locale's abbreviated month     |       |
        |           | name.                          |       |
        +-----------+--------------------------------+-------+
        | ``%B``    | Locale's full month name.      |       |
        +-----------+--------------------------------+-------+
        | ``%c``    | Locale's appropriate date and  |       |
        |           | time representation.           |       |
        +-----------+--------------------------------+-------+
        | ``%d``    | Day of the month as a decimal  |       |
        |           | number [01,31].                |       |
        +-----------+--------------------------------+-------+
        | ``%f``    | 'Fiscal' year without a        | \(1)  |
        |           | century  as a decimal number   |       |
        |           | [00,99]                        |       |
        +-----------+--------------------------------+-------+
        | ``%F``    | 'Fiscal' year with a century   | \(2)  |
        |           | as a decimal number            |       |
        +-----------+--------------------------------+-------+
        | ``%H``    | Hour (24-hour clock) as a      |       |
        |           | decimal number [00,23].        |       |
        +-----------+--------------------------------+-------+
        | ``%I``    | Hour (12-hour clock) as a      |       |
        |           | decimal number [01,12].        |       |
        +-----------+--------------------------------+-------+
        | ``%j``    | Day of the year as a decimal   |       |
        |           | number [001,366].              |       |
        +-----------+--------------------------------+-------+
        | ``%m``    | Month as a decimal number      |       |
        |           | [01,12].                       |       |
        +-----------+--------------------------------+-------+
        | ``%M``    | Minute as a decimal number     |       |
        |           | [00,59].                       |       |
        +-----------+--------------------------------+-------+
        | ``%p``    | Locale's equivalent of either  | \(3)  |
        |           | AM or PM.                      |       |
        +-----------+--------------------------------+-------+
        | ``%q``    | Quarter as a decimal number    |       |
        |           | [01,04]                        |       |
        +-----------+--------------------------------+-------+
        | ``%S``    | Second as a decimal number     | \(4)  |
        |           | [00,61].                       |       |
        +-----------+--------------------------------+-------+
        | ``%U``    | Week number of the year        | \(5)  |
        |           | (Sunday as the first day of    |       |
        |           | the week) as a decimal number  |       |
        |           | [00,53].  All days in a new    |       |
        |           | year preceding the first       |       |
        |           | Sunday are considered to be in |       |
        |           | week 0.                        |       |
        +-----------+--------------------------------+-------+
        | ``%w``    | Weekday as a decimal number    |       |
        |           | [0(Sunday),6].                 |       |
        +-----------+--------------------------------+-------+
        | ``%W``    | Week number of the year        | \(5)  |
        |           | (Monday as the first day of    |       |
        |           | the week) as a decimal number  |       |
        |           | [00,53].  All days in a new    |       |
        |           | year preceding the first       |       |
        |           | Monday are considered to be in |       |
        |           | week 0.                        |       |
        +-----------+--------------------------------+-------+
        | ``%x``    | Locale's appropriate date      |       |
        |           | representation.                |       |
        +-----------+--------------------------------+-------+
        | ``%X``    | Locale's appropriate time      |       |
        |           | representation.                |       |
        +-----------+--------------------------------+-------+
        | ``%y``    | Year without century as a      |       |
        |           | decimal number [00,99].        |       |
        +-----------+--------------------------------+-------+
        | ``%Y``    | Year with century as a decimal |       |
        |           | number.                        |       |
        +-----------+--------------------------------+-------+
        | ``%Z``    | Time zone name (no characters  |       |
        |           | if no time zone exists).       |       |
        +-----------+--------------------------------+-------+
        | ``%%``    | A literal ``'%'`` character.   |       |
        +-----------+--------------------------------+-------+

        .. note::

            (1)
                The ``%f`` directive is the same as ``%y`` if the frequency is
                not quarterly.
                Otherwise, it corresponds to the 'fiscal' year, as defined by
                the :attr:`qyear` attribute.

            (2)
                The ``%F`` directive is the same as ``%Y`` if the frequency is
                not quarterly.
                Otherwise, it corresponds to the 'fiscal' year, as defined by
                the :attr:`qyear` attribute.

            (3)
                The ``%p`` directive only affects the output hour field
                if the ``%I`` directive is used to parse the hour.

            (4)
                The range really is ``0`` to ``61``; this accounts for leap
                seconds and the (very rare) double leap seconds.

            (5)
                The ``%U`` and ``%W`` directives are only used in calculations
                when the day of the week and the year are specified.

        .. rubric::  Examples

            >>> a = Period(freq='Q@JUL', year=2006, quarter=1)
            >>> a.strftime('%F-Q%q')
            '2006-Q1'
            >>> # Output the last month in the quarter of this date
            >>> a.strftime('%b-%Y')
            'Oct-2005'
            >>>
            >>> a = Period(freq='D', year=2001, month=1, day=1)
            >>> a.strftime('%d-%b-%Y')
            '01-Jan-2006'
            >>> a.strftime('%b. %d, %Y was a %A')
            'Jan. 01, 2001 was a Monday'
        """
        base, mult = _gfc(self.freq)
        return tslib.period_format(self.ordinal, base, fmt)

def _get_ordinals(data, freq):
    f = lambda x: Period(x, freq=freq).ordinal
    if isinstance(data[0], Period):
        return tslib.extract_ordinals(data, freq)
    else:
        return lib.map_infer(data, f)


def dt64arr_to_periodarr(data, freq, tz):
    if data.dtype != np.dtype('M8[ns]'):
        raise ValueError('Wrong dtype: %s' % data.dtype)

    base, mult = _gfc(freq)
    return tslib.dt64arr_to_periodarr(data.view('i8'), base, tz)

# --- Period index sketch

def _period_index_cmp(opname, nat_result=False):
    """
    Wrap comparison operations to convert datetime-like to datetime64
    """
    def wrapper(self, other):
        if isinstance(other, Period):
            func = getattr(self.values, opname)
            if other.freq != self.freq:
                raise AssertionError("Frequencies must be equal")

            result = func(other.ordinal)
        elif isinstance(other, PeriodIndex):
            if other.freq != self.freq:
                raise AssertionError("Frequencies must be equal")

            result = getattr(self.values, opname)(other.values)

            mask = (com.mask_missing(self.values, tslib.iNaT) |
                    com.mask_missing(other.values, tslib.iNaT))
            if mask.any():
                result[mask] = nat_result

            return result
        else:
            other = Period(other, freq=self.freq)
            func = getattr(self.values, opname)
            result = func(other.ordinal)

        if other.ordinal == tslib.iNaT:
            result.fill(nat_result)
        mask = self.values == tslib.iNaT
        if mask.any():
            result[mask] = nat_result

        return result
    return wrapper


class PeriodIndex(DatetimeIndexOpsMixin, Int64Index):
    """
    Immutable ndarray holding ordinal values indicating regular periods in
    time such as particular years, quarters, months, etc. A value of 1 is the
    period containing the Gregorian proleptic datetime Jan 1, 0001 00:00:00.
    This ordinal representation is from the scikits.timeseries project.

    For instance,
        # construct period for day 1/1/1 and get the first second
        i = Period(year=1,month=1,day=1,freq='D').asfreq('S', 'S')
        i.ordinal
        ===> 1

    Index keys are boxed to Period objects which carries the metadata (eg,
    frequency information).

    Parameters
    ----------
    data  : array-like (1-dimensional), optional
        Optional period-like data to construct index with
    dtype : NumPy dtype (default: i8)
    copy  : bool
        Make a copy of input ndarray
    freq : string or period object, optional
        One of pandas period strings or corresponding objects
    start : starting value, period-like, optional
        If data is None, used as the start point in generating regular
        period data.
    periods  : int, optional, > 0
        Number of periods to generate, if generating index. Takes precedence
        over end argument
    end   : end value, period-like, optional
        If periods is none, generated index will extend to first conforming
        period on or just past end argument
    year : int, array, or Series, default None
    month : int, array, or Series, default None
    quarter : int, array, or Series, default None
    day : int, array, or Series, default None
    hour : int, array, or Series, default None
    minute : int, array, or Series, default None
    second : int, array, or Series, default None
    tz : object, default None
        Timezone for converting datetime64 data to Periods

    Examples
    --------
    >>> idx = PeriodIndex(year=year_arr, quarter=q_arr)

    >>> idx2 = PeriodIndex(start='2000', end='2010', freq='A')
    """
    _box_scalars = True
    _typ = 'periodindex'
    _attributes = ['name','freq']
    _datetimelike_ops = ['year','month','day','hour','minute','second',
                         'weekofyear','week','dayofweek','weekday','dayofyear','quarter', 'qyear', 'freq']
    _is_numeric_dtype = False
    freq = None

    __eq__ = _period_index_cmp('__eq__')
    __ne__ = _period_index_cmp('__ne__', nat_result=True)
    __lt__ = _period_index_cmp('__lt__')
    __gt__ = _period_index_cmp('__gt__')
    __le__ = _period_index_cmp('__le__')
    __ge__ = _period_index_cmp('__ge__')

    def __new__(cls, data=None, ordinal=None, freq=None, start=None, end=None,
                periods=None, copy=False, name=None, tz=None, **kwargs):

        freq = frequencies.get_standard_freq(freq)

        if periods is not None:
            if com.is_float(periods):
                periods = int(periods)
            elif not com.is_integer(periods):
                raise ValueError('Periods must be a number, got %s' %
                                 str(periods))

        if data is None:
            if ordinal is not None:
                data = np.asarray(ordinal, dtype=np.int64)
            else:
                data, freq = cls._generate_range(start, end, periods,
                                                 freq, kwargs)
        else:
            ordinal, freq = cls._from_arraylike(data, freq, tz)
            data = np.array(ordinal, dtype=np.int64, copy=False)

        return cls._simple_new(data, name=name, freq=freq)

    @classmethod
    def _generate_range(cls, start, end, periods, freq, fields):
        field_count = len(fields)
        if com._count_not_none(start, end) > 0:
            if field_count > 0:
                raise ValueError('Can either instantiate from fields '
                                 'or endpoints, but not both')
            subarr, freq = _get_ordinal_range(start, end, periods, freq)
        elif field_count > 0:
            subarr, freq = _range_from_fields(freq=freq, **fields)
        else:
            raise ValueError('Not enough parameters to construct '
                             'Period range')

        return subarr, freq

    @classmethod
    def _from_arraylike(cls, data, freq, tz):

        if not isinstance(data, (np.ndarray, PeriodIndex, DatetimeIndex, Int64Index)):
            if np.isscalar(data) or isinstance(data, Period):
                raise ValueError('PeriodIndex() must be called with a '
                                 'collection of some kind, %s was passed'
                                 % repr(data))

            # other iterable of some kind
            if not isinstance(data, (list, tuple)):
                data = list(data)

            try:
                data = com._ensure_int64(data)
                if freq is None:
                    raise ValueError('freq not specified')
                data = np.array([Period(x, freq=freq).ordinal for x in data],
                                dtype=np.int64)
            except (TypeError, ValueError):
                data = com._ensure_object(data)

                if freq is None and len(data) > 0:
                    freq = getattr(data[0], 'freq', None)

                if freq is None:
                    raise ValueError('freq not specified and cannot be '
                                     'inferred from first element')

                data = _get_ordinals(data, freq)
        else:
            if isinstance(data, PeriodIndex):
                if freq is None or freq == data.freq:
                    freq = data.freq
                    data = data.values
                else:
                    base1, _ = _gfc(data.freq)
                    base2, _ = _gfc(freq)
                    data = tslib.period_asfreq_arr(data.values, base1,
                                                   base2, 1)
            else:
                if freq is None and len(data) > 0:
                    freq = getattr(data[0], 'freq', None)

                if freq is None:
                    raise ValueError('freq not specified and cannot be '
                                     'inferred from first element')

                if data.dtype != np.int64:
                    if np.issubdtype(data.dtype, np.datetime64):
                        data = dt64arr_to_periodarr(data, freq, tz)
                    else:
                        try:
                            data = com._ensure_int64(data)
                        except (TypeError, ValueError):
                            data = com._ensure_object(data)
                            data = _get_ordinals(data, freq)

        return data, freq

    @classmethod
    def _simple_new(cls, values, name=None, freq=None, **kwargs):
        result = object.__new__(cls)
        result._data = values
        result.name = name
        result.freq = freq
        result._reset_identity()
        return result

    @property
    def _na_value(self):
        return self._box_func(tslib.iNaT)

    def __contains__(self, key):
        if not isinstance(key, Period) or key.freq != self.freq:
            if isinstance(key, compat.string_types):
                try:
                    self.get_loc(key)
                    return True
                except Exception:
                    return False
            return False
        return key.ordinal in self._engine

    @property
    def _box_func(self):
        return lambda x: Period._from_ordinal(ordinal=x, freq=self.freq)

    def _to_embed(self, keep_tz=False):
        """ return an array repr of this object, potentially casting to object """
        return self.asobject.values

    def asof_locs(self, where, mask):
        """
        where : array of timestamps
        mask : array of booleans where data is not NA

        """
        where_idx = where
        if isinstance(where_idx, DatetimeIndex):
            where_idx = PeriodIndex(where_idx.values, freq=self.freq)

        locs = self.values[mask].searchsorted(where_idx.values, side='right')

        locs = np.where(locs > 0, locs - 1, 0)
        result = np.arange(len(self))[mask].take(locs)

        first = mask.argmax()
        result[(locs == 0) & (where_idx.values < self.values[first])] = -1

        return result

    def _array_values(self):
        return self.asobject

    def astype(self, dtype):
        dtype = np.dtype(dtype)
        if dtype == np.object_:
            return Index(np.array(list(self), dtype), dtype)
        elif dtype == _INT64_DTYPE:
            return Index(self.values, dtype)
        raise ValueError('Cannot cast PeriodIndex to dtype %s' % dtype)

    def searchsorted(self, key, side='left'):
        if isinstance(key, Period):
            if key.freq != self.freq:
                raise ValueError("Different period frequency: %s" % key.freq)
            key = key.ordinal
        elif isinstance(key, compat.string_types):
            key = Period(key, freq=self.freq).ordinal

        return self.values.searchsorted(key, side=side)

    @property
    def is_all_dates(self):
        return True

    @property
    def is_full(self):
        """
        Returns True if there are any missing periods from start to end
        """
        if len(self) == 0:
            return True
        if not self.is_monotonic:
            raise ValueError('Index is not monotonic')
        values = self.values
        return ((values[1:] - values[:-1]) < 2).all()

    @property
    def freqstr(self):
        return self.freq

    def asfreq(self, freq=None, how='E'):
        how = _validate_end_alias(how)

        freq = frequencies.get_standard_freq(freq)

        base1, mult1 = _gfc(self.freq)
        base2, mult2 = _gfc(freq)

        if mult2 != 1:
            raise ValueError('Only mult == 1 supported')

        end = how == 'E'
        new_data = tslib.period_asfreq_arr(self.values, base1, base2, end)
        return self._simple_new(new_data, self.name, freq=freq)

    def to_datetime(self, dayfirst=False):
        return self.to_timestamp()

    year = _field_accessor('year', 0, "The year of the period")
    month = _field_accessor('month', 3, "The month as January=1, December=12")
    day = _field_accessor('day', 4, "The days of the period")
    hour = _field_accessor('hour', 5, "The hour of the period")
    minute = _field_accessor('minute', 6, "The minute of the period")
    second = _field_accessor('second', 7, "The second of the period")
    weekofyear = _field_accessor('week', 8, "The week ordinal of the year")
    week = weekofyear
    dayofweek = _field_accessor('dayofweek', 10, "The day of the week with Monday=0, Sunday=6")
    weekday = dayofweek
    dayofyear = day_of_year = _field_accessor('dayofyear', 9, "The ordinal day of the year")
    quarter = _field_accessor('quarter', 2, "The quarter of the date")
    qyear = _field_accessor('qyear', 1)

    def _get_object_array(self):
        freq = self.freq
        return np.array([ Period._from_ordinal(ordinal=x, freq=freq) for x in self.values], copy=False)

    def _mpl_repr(self):
        # how to represent ourselves to matplotlib
        return self._get_object_array()

    def equals(self, other):
        """
        Determines if two Index objects contain the same elements.
        """
        if self.is_(other):
            return True

        if (not hasattr(other, 'inferred_type') or
                other.inferred_type != 'int64'):
            try:
                other = PeriodIndex(other)
            except:
                return False

        return np.array_equal(self.asi8, other.asi8)

    def to_timestamp(self, freq=None, how='start'):
        """
        Cast to DatetimeIndex

        Parameters
        ----------
        freq : string or DateOffset, default 'D' for week or longer, 'S'
               otherwise
            Target frequency
        how : {'s', 'e', 'start', 'end'}

        Returns
        -------
        DatetimeIndex
        """
        how = _validate_end_alias(how)

        if freq is None:
            base, mult = _gfc(self.freq)
            freq = frequencies.get_to_timestamp_base(base)

        base, mult = _gfc(freq)
        new_data = self.asfreq(freq, how)

        new_data = tslib.periodarr_to_dt64arr(new_data.values, base)
        return DatetimeIndex(new_data, freq='infer', name=self.name)

    def _add_delta(self, other):
        if isinstance(other, (timedelta, np.timedelta64, offsets.Tick, Timedelta)):
            offset = frequencies.to_offset(self.freq)
            if isinstance(offset, offsets.Tick):
                nanos = tslib._delta_to_nanoseconds(other)
                offset_nanos = tslib._delta_to_nanoseconds(offset)
                if nanos % offset_nanos == 0:
                    return self.shift(nanos // offset_nanos)
        elif isinstance(other, offsets.DateOffset):
            freqstr = frequencies.get_standard_freq(other)
            base = frequencies.get_base_alias(freqstr)

            if base == self.freq:
                return self.shift(other.n)
        raise ValueError("Input has different freq from PeriodIndex(freq={0})".format(self.freq))

    def shift(self, n):
        """
        Specialized shift which produces an PeriodIndex

        Parameters
        ----------
        n : int
            Periods to shift by
        freq : freq string

        Returns
        -------
        shifted : PeriodIndex
        """
        mask = self.values == tslib.iNaT
        values = self.values + n
        values[mask] = tslib.iNaT
        return PeriodIndex(data=values, name=self.name, freq=self.freq)

    @property
    def inferred_type(self):
        # b/c data is represented as ints make sure we can't have ambiguous
        # indexing
        return 'period'

    def get_value(self, series, key):
        """
        Fast lookup of value from 1-dimensional ndarray. Only use this if you
        know what you're doing
        """
        s = _values_from_object(series)
        try:
            return _maybe_box(self, super(PeriodIndex, self).get_value(s, key), series, key)
        except (KeyError, IndexError):
            try:
                asdt, parsed, reso = parse_time_string(key, self.freq)
                grp = frequencies._infer_period_group(reso)
                freqn = frequencies._period_group(self.freq)

                vals = self.values

                # if our data is higher resolution than requested key, slice
                if grp < freqn:
                    iv = Period(asdt, freq=(grp, 1))
                    ord1 = iv.asfreq(self.freq, how='S').ordinal
                    ord2 = iv.asfreq(self.freq, how='E').ordinal

                    if ord2 < vals[0] or ord1 > vals[-1]:
                        raise KeyError(key)

                    pos = np.searchsorted(self.values, [ord1, ord2])
                    key = slice(pos[0], pos[1] + 1)
                    return series[key]
                elif grp == freqn:
                    key = Period(asdt, freq=self.freq).ordinal
                    return _maybe_box(self, self._engine.get_value(s, key), series, key)
                else:
                    raise KeyError(key)
            except TypeError:
                pass

            key = Period(key, self.freq).ordinal
            return _maybe_box(self, self._engine.get_value(s, key), series, key)

    def get_loc(self, key):
        """
        Get integer location for requested label

        Returns
        -------
        loc : int
        """
        try:
            return self._engine.get_loc(key)
        except KeyError:
            if com.is_integer(key):
                raise

            try:
                asdt, parsed, reso = parse_time_string(key, self.freq)
                key = asdt
            except TypeError:
                pass

            key = Period(key, self.freq)
            try:
                return self._engine.get_loc(key.ordinal)
            except KeyError:
                raise KeyError(key)

    def _maybe_cast_slice_bound(self, label, side):
        """
        If label is a string or a datetime, cast it to Period.ordinal according to
        resolution.

        Parameters
        ----------
        label : object
        side : {'left', 'right'}

        Returns
        -------
        bound : Period or object

        Notes
        -----
        Value of `side` parameter should be validated in caller.

        """
        if isinstance(label, datetime):
            return Period(label, freq=self.freq)
        elif isinstance(label, compat.string_types):
            try:
                _, parsed, reso = parse_time_string(label, self.freq)
                bounds = self._parsed_string_to_bounds(reso, parsed)
                return bounds[0 if side == 'left' else 1]
            except Exception:
                raise KeyError(label)

        return label

    def _parsed_string_to_bounds(self, reso, parsed):
        if reso == 'year':
            t1 = Period(year=parsed.year, freq='A')
        elif reso == 'month':
            t1 = Period(year=parsed.year, month=parsed.month, freq='M')
        elif reso == 'quarter':
            q = (parsed.month - 1) // 3 + 1
            t1 = Period(year=parsed.year, quarter=q, freq='Q-DEC')
        elif reso == 'day':
            t1 = Period(year=parsed.year, month=parsed.month, day=parsed.day,
                        freq='D')
        elif reso == 'hour':
            t1 = Period(year=parsed.year, month=parsed.month, day=parsed.day,
                        hour=parsed.hour, freq='H')
        elif reso == 'minute':
            t1 = Period(year=parsed.year, month=parsed.month, day=parsed.day,
                        hour=parsed.hour, minute=parsed.minute, freq='T')
        elif reso == 'second':
            t1 = Period(year=parsed.year, month=parsed.month, day=parsed.day,
                        hour=parsed.hour, minute=parsed.minute, second=parsed.second,
                        freq='S')
        else:
            raise KeyError(key)
        return (t1.asfreq(self.freq, how='start'),
                t1.asfreq(self.freq, how='end'))

    def _get_string_slice(self, key):
        if not self.is_monotonic:
            raise ValueError('Partial indexing only valid for '
                             'ordered time series')

        key, parsed, reso = parse_time_string(key, self.freq)

        grp = frequencies._infer_period_group(reso)
        freqn = frequencies._period_group(self.freq)
        if reso in ['day', 'hour', 'minute', 'second'] and not grp < freqn:
            raise KeyError(key)

        t1, t2 = self._parsed_string_to_bounds(reso, parsed)
        return slice(self.searchsorted(t1.ordinal, side='left'),
                     self.searchsorted(t2.ordinal, side='right'))

    def join(self, other, how='left', level=None, return_indexers=False):
        """
        See Index.join
        """
        self._assert_can_do_setop(other)

        result = Int64Index.join(self, other, how=how, level=level,
                                 return_indexers=return_indexers)

        if return_indexers:
            result, lidx, ridx = result
            return self._apply_meta(result), lidx, ridx
        return self._apply_meta(result)

    def _assert_can_do_setop(self, other):
        if not isinstance(other, PeriodIndex):
            raise ValueError('can only call with other PeriodIndex-ed objects')

        if self.freq != other.freq:
            raise ValueError('Only like-indexed PeriodIndexes compatible '
                             'for join (for now)')

    def _wrap_union_result(self, other, result):
        name = self.name if self.name == other.name else None
        result = self._apply_meta(result)
        result.name = name
        return result

    def _apply_meta(self, rawarr):
        if not isinstance(rawarr, PeriodIndex):
            rawarr = PeriodIndex(rawarr, freq=self.freq)
        return rawarr

    def __getitem__(self, key):
        getitem = self._data.__getitem__
        if np.isscalar(key):
            val = getitem(key)
            return Period(ordinal=val, freq=self.freq)
        else:
            if com._is_bool_indexer(key):
                key = np.asarray(key)

            result = getitem(key)
            if result.ndim > 1:
                # MPL kludge
                # values = np.asarray(list(values), dtype=object)
                # return values.reshape(result.shape)

                return PeriodIndex(result, name=self.name, freq=self.freq)

            return PeriodIndex(result, name=self.name, freq=self.freq)

    def _format_native_types(self, na_rep=u('NaT'), **kwargs):

        values = np.array(list(self), dtype=object)
        mask = isnull(self.values)
        values[mask] = na_rep

        imask = ~mask
        values[imask] = np.array([u('%s') % dt for dt in values[imask]])
        return values.tolist()

    def __array_finalize__(self, obj):
        if not self.ndim:  # pragma: no cover
            return self.item()

        self.freq = getattr(obj, 'freq', None)
        self.name = getattr(obj, 'name', None)
        self._reset_identity()

    def _format_footer(self):
        tagline = 'Length: %d, Freq: %s'
        return tagline % (len(self), self.freqstr)

    def take(self, indices, axis=None):
        """
        Analogous to ndarray.take
        """
        indices = com._ensure_platform_int(indices)
        taken = self.values.take(indices, axis=axis)
        return self._simple_new(taken, self.name, freq=self.freq)

    def append(self, other):
        """
        Append a collection of Index options together

        Parameters
        ----------
        other : Index or list/tuple of indices

        Returns
        -------
        appended : Index
        """
        name = self.name
        to_concat = [self]

        if isinstance(other, (list, tuple)):
            to_concat = to_concat + list(other)
        else:
            to_concat.append(other)

        for obj in to_concat:
            if isinstance(obj, Index) and obj.name != name:
                name = None
                break

        to_concat = self._ensure_compat_concat(to_concat)

        if isinstance(to_concat[0], PeriodIndex):
            if len(set([x.freq for x in to_concat])) > 1:
                # box
                to_concat = [x.asobject.values for x in to_concat]
            else:
                cat_values = np.concatenate([x.values for x in to_concat])
                return PeriodIndex(cat_values, freq=self.freq, name=name)

        to_concat = [x.values if isinstance(x, Index) else x
                     for x in to_concat]
        return Index(com._concat_compat(to_concat), name=name)

    def __setstate__(self, state):
        """Necessary for making this object picklable"""

        if isinstance(state, dict):
            super(PeriodIndex, self).__setstate__(state)

        elif isinstance(state, tuple):

            # < 0.15 compat
            if len(state) == 2:
                nd_state, own_state = state
                data = np.empty(nd_state[1], dtype=nd_state[2])
                np.ndarray.__setstate__(data, nd_state)

                try:  # backcompat
                    self.freq = own_state[1]
                except:
                    pass

            else:  # pragma: no cover
                data = np.empty(state)
                np.ndarray.__setstate__(self, state)

            self._data = data

        else:
            raise Exception("invalid pickle state")
    _unpickle_compat = __setstate__

    def tz_convert(self, tz):
        """
        Convert tz-aware DatetimeIndex from one time zone to another (using pytz/dateutil)

        Parameters
        ----------
        tz : string, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for time. Corresponding timestamps would be converted to
            time zone of the TimeSeries.
            None will remove timezone holding UTC time.

        Returns
        -------
        normalized : DatetimeIndex

        Note
        ----
        Not currently implemented for PeriodIndex
        """
        raise NotImplementedError("Not yet implemented for PeriodIndex")

    def tz_localize(self, tz, infer_dst=False):
        """
        Localize tz-naive DatetimeIndex to given time zone (using pytz/dateutil),
        or remove timezone from tz-aware DatetimeIndex

        Parameters
        ----------
        tz : string, pytz.timezone, dateutil.tz.tzfile or None
            Time zone for time. Corresponding timestamps would be converted to
            time zone of the TimeSeries.
            None will remove timezone holding local time.
        infer_dst : boolean, default False
            Attempt to infer fall dst-transition hours based on order

        Returns
        -------
        localized : DatetimeIndex

        Note
        ----
        Not currently implemented for PeriodIndex
        """
        raise NotImplementedError("Not yet implemented for PeriodIndex")


PeriodIndex._add_numeric_methods_disabled()
PeriodIndex._add_logical_methods_disabled()
PeriodIndex._add_datetimelike_methods()


def _get_ordinal_range(start, end, periods, freq):
    if com._count_not_none(start, end, periods) < 2:
        raise ValueError('Must specify 2 of start, end, periods')

    if start is not None:
        start = Period(start, freq)
    if end is not None:
        end = Period(end, freq)

    is_start_per = isinstance(start, Period)
    is_end_per = isinstance(end, Period)

    if is_start_per and is_end_per and start.freq != end.freq:
        raise ValueError('Start and end must have same freq')
    if ((is_start_per and start.ordinal == tslib.iNaT) or
        (is_end_per and end.ordinal == tslib.iNaT)):
        raise ValueError('Start and end must not be NaT')

    if freq is None:
        if is_start_per:
            freq = start.freq
        elif is_end_per:
            freq = end.freq
        else:  # pragma: no cover
            raise ValueError('Could not infer freq from start/end')

    if periods is not None:
        if start is None:
            data = np.arange(end.ordinal - periods + 1,
                             end.ordinal + 1,
                             dtype=np.int64)
        else:
            data = np.arange(start.ordinal, start.ordinal + periods,
                             dtype=np.int64)
    else:
        data = np.arange(start.ordinal, end.ordinal + 1, dtype=np.int64)

    return data, freq


def _range_from_fields(year=None, month=None, quarter=None, day=None,
                       hour=None, minute=None, second=None, freq=None):
    if hour is None:
        hour = 0
    if minute is None:
        minute = 0
    if second is None:
        second = 0
    if day is None:
        day = 1

    ordinals = []

    if quarter is not None:
        if freq is None:
            freq = 'Q'
            base = frequencies.FreqGroup.FR_QTR
        else:
            base, mult = _gfc(freq)
            if mult != 1:
                raise ValueError('Only mult == 1 supported')
            if base != frequencies.FreqGroup.FR_QTR:
                raise AssertionError("base must equal FR_QTR")

        year, quarter = _make_field_arrays(year, quarter)
        for y, q in zip(year, quarter):
            y, m = _quarter_to_myear(y, q, freq)
            val = tslib.period_ordinal(y, m, 1, 1, 1, 1, 0, 0, base)
            ordinals.append(val)
    else:
        base, mult = _gfc(freq)
        if mult != 1:
            raise ValueError('Only mult == 1 supported')

        arrays = _make_field_arrays(year, month, day, hour, minute, second)
        for y, mth, d, h, mn, s in zip(*arrays):
            ordinals.append(tslib.period_ordinal(y, mth, d, h, mn, s, 0, 0, base))

    return np.array(ordinals, dtype=np.int64), freq


def _make_field_arrays(*fields):
    length = None
    for x in fields:
        if isinstance(x, (list, np.ndarray, ABCSeries)):
            if length is not None and len(x) != length:
                raise ValueError('Mismatched Period array lengths')
            elif length is None:
                length = len(x)

    arrays = [np.asarray(x) if isinstance(x, (np.ndarray, list, ABCSeries))
              else np.repeat(x, length) for x in fields]

    return arrays


def _ordinal_from_fields(year, month, quarter, day, hour, minute,
                         second, freq):
    base, mult = _gfc(freq)
    if mult != 1:
        raise ValueError('Only mult == 1 supported')

    if quarter is not None:
        year, month = _quarter_to_myear(year, quarter, freq)

    return tslib.period_ordinal(year, month, day, hour, minute, second, 0, 0, base)


def _quarter_to_myear(year, quarter, freq):
    if quarter is not None:
        if quarter <= 0 or quarter > 4:
            raise ValueError('Quarter must be 1 <= q <= 4')

        mnum = frequencies._month_numbers[frequencies._get_rule_month(freq)] + 1
        month = (mnum + (quarter - 1) * 3) % 12 + 1
        if month > mnum:
            year -= 1

    return year, month


def _validate_end_alias(how):
    how_dict = {'S': 'S', 'E': 'E',
                'START': 'S', 'FINISH': 'E',
                'BEGIN': 'S', 'END': 'E'}
    how = how_dict.get(str(how).upper())
    if how not in set(['S', 'E']):
        raise ValueError('How must be one of S or E')
    return how


def pnow(freq=None):
    return Period(datetime.now(), freq=freq)


def period_range(start=None, end=None, periods=None, freq='D', name=None):
    """
    Return a fixed frequency datetime index, with day (calendar) as the default
    frequency


    Parameters
    ----------
    start :
    end :
    periods : int, default None
        Number of periods in the index
    freq : str/DateOffset, default 'D'
        Frequency alias
    name : str, default None
        Name for the resulting PeriodIndex

    Returns
    -------
    prng : PeriodIndex
    """
    return PeriodIndex(start=start, end=end, periods=periods,
                       freq=freq, name=name)
