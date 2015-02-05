from datetime import datetime, date, timedelta
import operator
import numpy as np

from pandas import compat
from pandas.core import common as com
from pandas.core.base import PandasObject

from pandas.tseries import frequencies
from pandas.tseries.frequencies import get_freq_code as _gfc
from pandas.tseries import offsets
from pandas.tseries.tools import parse_time_string

from pandas import tslib
from tslib import Timedelta, Timestamp


#---------------
# Period logic

def _period_field_accessor(name, alias):
    def f(self):
        base, mult = _gfc(self.freq)
        return tslib.get_period_field(alias, self.ordinal, base)
    f.__name__ = name
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

        elif com.is_null_datelike_scalar(value) or value in tslib._nat_strings:
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
