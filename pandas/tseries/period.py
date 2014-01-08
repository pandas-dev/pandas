# pylint: disable=E1101,E1103,W0232
import operator

from datetime import datetime, date
import numpy as np
from pandas.core.base import PandasObject

from pandas.tseries.frequencies import (get_freq_code as _gfc,
                                        _month_numbers, FreqGroup)
from pandas.tseries.index import DatetimeIndex, Int64Index, Index
from pandas.tseries.tools import parse_time_string
import pandas.tseries.frequencies as _freq_mod

import pandas.core.common as com
from pandas.core.common import (isnull, _INT64_DTYPE, _maybe_box,
                                _values_from_object)
from pandas import compat
from pandas.lib import Timestamp
import pandas.lib as lib
import pandas.tslib as tslib
import pandas.algos as _algos
from pandas.compat import map, zip, u


#---------------
# Period logic

def _period_field_accessor(name, alias):
    def f(self):
        base, mult = _gfc(self.freq)
        return tslib.get_period_field(alias, self.ordinal, base)
    f.__name__ = name
    return property(f)


def _field_accessor(name, alias):
    def f(self):
        base, mult = _gfc(self.freq)
        return tslib.get_period_field_arr(alias, self.values, base)
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

        elif isinstance(value, compat.string_types) or com.is_integer(value):
            if com.is_integer(value):
                value = str(value)

            dt, freq = _get_date_and_freq(value, freq)

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

        self.freq = _freq_mod._get_freq_str(base)

    def __eq__(self, other):
        if isinstance(other, Period):
            if other.freq != self.freq:
                raise ValueError("Cannot compare non-conforming periods")
            return (self.ordinal == other.ordinal
                    and _gfc(self.freq) == _gfc(other.freq))
        else:
            raise TypeError(other)
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.ordinal, self.freq))

    def __add__(self, other):
        if com.is_integer(other):
            return Period(ordinal=self.ordinal + other, freq=self.freq)
        else:  # pragma: no cover
            raise TypeError(other)

    def __sub__(self, other):
        if com.is_integer(other):
            return Period(ordinal=self.ordinal - other, freq=self.freq)
        if isinstance(other, Period):
            if other.freq != self.freq:
                raise ValueError("Cannot do arithmetic with "
                                 "non-conforming periods")
            return self.ordinal - other.ordinal
        else:  # pragma: no cover
            raise TypeError(other)

    def _comp_method(func, name):
        def f(self, other):
            if isinstance(other, Period):
                if other.freq != self.freq:
                    raise ValueError("Cannot compare non-conforming periods")
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
            freq = _freq_mod.get_to_timestamp_base(base)

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
        freqstr = _freq_mod._reverse_period_code_map[base]

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


def _get_date_and_freq(value, freq):
    value = value.upper()
    dt, _, reso = parse_time_string(value, freq)

    if freq is None:
        if reso == 'year':
            freq = 'A'
        elif reso == 'quarter':
            freq = 'Q'
        elif reso == 'month':
            freq = 'M'
        elif reso == 'day':
            freq = 'D'
        elif reso == 'hour':
            freq = 'H'
        elif reso == 'minute':
            freq = 'T'
        elif reso == 'second':
            freq = 'S'
        elif reso == 'microsecond':
            if dt.microsecond % 1000 == 0:
                freq = 'L'
            else:
                freq = 'U'
        else:
            raise ValueError("Invalid frequency or could not infer: %s" % reso)

    return dt, freq


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

def _period_index_cmp(opname):
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
            return getattr(self.values, opname)(other.values)
        else:
            other = Period(other, freq=self.freq)
            func = getattr(self.values, opname)
            result = func(other.ordinal)

        return result
    return wrapper


class PeriodIndex(Int64Index):
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
    year : int or array, default None
    month : int or array, default None
    quarter : int or array, default None
    day : int or array, default None
    hour : int or array, default None
    minute : int or array, default None
    second : int or array, default None
    tz : object, default None
        Timezone for converting datetime64 data to Periods

    Examples
    --------
    >>> idx = PeriodIndex(year=year_arr, quarter=q_arr)

    >>> idx2 = PeriodIndex(start='2000', end='2010', freq='A')
    """
    _box_scalars = True

    __eq__ = _period_index_cmp('__eq__')
    __ne__ = _period_index_cmp('__ne__')
    __lt__ = _period_index_cmp('__lt__')
    __gt__ = _period_index_cmp('__gt__')
    __le__ = _period_index_cmp('__le__')
    __ge__ = _period_index_cmp('__ge__')

    def __new__(cls, data=None, ordinal=None, freq=None, start=None, end=None,
                periods=None, copy=False, name=None, year=None, month=None,
                quarter=None, day=None, hour=None, minute=None, second=None,
                tz=None):

        freq = _freq_mod.get_standard_freq(freq)

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
                fields = [year, month, quarter, day, hour, minute, second]
                data, freq = cls._generate_range(start, end, periods,
                                                 freq, fields)
        else:
            ordinal, freq = cls._from_arraylike(data, freq, tz)
            data = np.array(ordinal, dtype=np.int64, copy=False)

        subarr = data.view(cls)
        subarr.name = name
        subarr.freq = freq

        return subarr

    @classmethod
    def _generate_range(cls, start, end, periods, freq, fields):
        field_count = com._count_not_none(*fields)
        if com._count_not_none(start, end) > 0:
            if field_count > 0:
                raise ValueError('Can either instantiate from fields '
                                 'or endpoints, but not both')
            subarr, freq = _get_ordinal_range(start, end, periods, freq)
        elif field_count > 0:
            y, mth, q, d, h, minute, s = fields
            subarr, freq = _range_from_fields(year=y, month=mth, quarter=q,
                                              day=d, hour=h, minute=minute,
                                              second=s, freq=freq)
        else:
            raise ValueError('Not enough parameters to construct '
                             'Period range')

        return subarr, freq

    @classmethod
    def _from_arraylike(cls, data, freq, tz):
        if not isinstance(data, np.ndarray):
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

    def _box_values(self, values):
        f = lambda x: Period(ordinal=x, freq=self.freq)
        return lib.map_infer(values, f)

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

    @property
    def asobject(self):
        return Index(self._box_values(self.values), dtype=object)

    def _array_values(self):
        return self.asobject

    def astype(self, dtype):
        dtype = np.dtype(dtype)
        if dtype == np.object_:
            return Index(np.array(list(self), dtype), dtype)
        elif dtype == _INT64_DTYPE:
            return Index(self.values, dtype)
        raise ValueError('Cannot cast PeriodIndex to dtype %s' % dtype)

    def __iter__(self):
        for val in self.values:
            yield Period(ordinal=val, freq=self.freq)

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

    def factorize(self):
        """
        Specialized factorize that boxes uniques
        """
        from pandas.core.algorithms import factorize
        labels, uniques = factorize(self.values)
        uniques = PeriodIndex(ordinal=uniques, freq=self.freq)
        return labels, uniques

    @property
    def freqstr(self):
        return self.freq

    def asfreq(self, freq=None, how='E'):
        how = _validate_end_alias(how)

        freq = _freq_mod.get_standard_freq(freq)

        base1, mult1 = _gfc(self.freq)
        base2, mult2 = _gfc(freq)

        if mult2 != 1:
            raise ValueError('Only mult == 1 supported')

        end = how == 'E'
        new_data = tslib.period_asfreq_arr(self.values, base1, base2, end)

        result = new_data.view(PeriodIndex)
        result.name = self.name
        result.freq = freq
        return result

    def to_datetime(self, dayfirst=False):
        return self.to_timestamp()

    year = _field_accessor('year', 0)
    month = _field_accessor('month', 3)
    day = _field_accessor('day', 4)
    hour = _field_accessor('hour', 5)
    minute = _field_accessor('minute', 6)
    second = _field_accessor('second', 7)
    weekofyear = _field_accessor('week', 8)
    week = weekofyear
    dayofweek = _field_accessor('dayofweek', 10)
    weekday = dayofweek
    dayofyear = day_of_year = _field_accessor('dayofyear', 9)
    quarter = _field_accessor('quarter', 2)
    qyear = _field_accessor('qyear', 1)

    # Try to run function on index first, and then on elements of index
    # Especially important for group-by functionality
    def map(self, f):
        try:
            result = f(self)
            if not isinstance(result, np.ndarray):
                raise TypeError
            return result
        except Exception:
            return _algos.arrmap_object(self.asobject, f)

    def _get_object_array(self):
        freq = self.freq
        boxfunc = lambda x: Period(ordinal=x, freq=freq)
        boxer = np.frompyfunc(boxfunc, 1, 1)
        return boxer(self.values)

    def _mpl_repr(self):
        # how to represent ourselves to matplotlib
        return self._get_object_array()

    def equals(self, other):
        """
        Determines if two Index objects contain the same elements.
        """
        if self.is_(other):
            return True

        return np.array_equal(self.asi8, other.asi8)

    def tolist(self):
        """
        Return a list of Period objects
        """
        return self._get_object_array().tolist()

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
            freq = _freq_mod.get_to_timestamp_base(base)

        base, mult = _gfc(freq)
        new_data = self.asfreq(freq, how)

        new_data = tslib.periodarr_to_dt64arr(new_data.values, base)
        return DatetimeIndex(new_data, freq='infer', name=self.name)

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
        if n == 0:
            return self

        return PeriodIndex(data=self.values + n, freq=self.freq)

    def __add__(self, other):
        return PeriodIndex(ordinal=self.values + other, freq=self.freq)

    def __sub__(self, other):
        return PeriodIndex(ordinal=self.values - other, freq=self.freq)

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
                grp = _freq_mod._infer_period_group(reso)
                freqn = _freq_mod._period_group(self.freq)

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
                else:
                    key = Period(asdt, freq=self.freq).ordinal
                    return _maybe_box(self, self._engine.get_value(s, key), series, key)
            except TypeError:
                pass
            except KeyError:
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

    def slice_locs(self, start=None, end=None):
        """
        Index.slice_locs, customized to handle partial ISO-8601 string slicing
        """
        if isinstance(start, compat.string_types) or isinstance(end, compat.string_types):
            try:
                if start:
                    start_loc = self._get_string_slice(start).start
                else:
                    start_loc = 0

                if end:
                    end_loc = self._get_string_slice(end).stop
                else:
                    end_loc = len(self)

                return start_loc, end_loc
            except KeyError:
                pass

        if isinstance(start, datetime) and isinstance(end, datetime):
            ordinals = self.values
            t1 = Period(start, freq=self.freq)
            t2 = Period(end, freq=self.freq)

            left = ordinals.searchsorted(t1.ordinal, side='left')
            right = ordinals.searchsorted(t2.ordinal, side='right')
            return left, right

        return Int64Index.slice_locs(self, start, end)

    def _get_string_slice(self, key):
        if not self.is_monotonic:
            raise ValueError('Partial indexing only valid for '
                             'ordered time series')

        asdt, parsed, reso = parse_time_string(key, self.freq)
        key = asdt

        if reso == 'year':
            t1 = Period(year=parsed.year, freq='A')
        elif reso == 'month':
            t1 = Period(year=parsed.year, month=parsed.month, freq='M')
        elif reso == 'quarter':
            q = (parsed.month - 1) // 3 + 1
            t1 = Period(year=parsed.year, quarter=q, freq='Q-DEC')
        else:
            raise KeyError(key)

        ordinals = self.values

        t2 = t1.asfreq(self.freq, how='end')
        t1 = t1.asfreq(self.freq, how='start')

        left = ordinals.searchsorted(t1.ordinal, side='left')
        right = ordinals.searchsorted(t2.ordinal, side='right')
        return slice(left, right)

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
            rawarr = rawarr.view(PeriodIndex)
        rawarr.freq = self.freq
        return rawarr

    def __getitem__(self, key):
        """Override numpy.ndarray's __getitem__ method to work as desired"""
        arr_idx = self.view(np.ndarray)
        if np.isscalar(key):
            val = arr_idx[key]
            return Period(ordinal=val, freq=self.freq)
        else:
            if com._is_bool_indexer(key):
                key = np.asarray(key)

            result = arr_idx[key]
            if result.ndim > 1:
                # MPL kludge
                # values = np.asarray(list(values), dtype=object)
                # return values.reshape(result.shape)

                return PeriodIndex(result, name=self.name, freq=self.freq)

            return PeriodIndex(result, name=self.name, freq=self.freq)

    _getitem_slice = __getitem__

    def _format_with_header(self, header, **kwargs):
        return header + self._format_native_types(**kwargs)

    def _format_native_types(self, na_rep=u('NaT'), **kwargs):

        values = np.array(list(self), dtype=object)
        mask = isnull(self.values)
        values[mask] = na_rep

        imask = -mask
        values[imask] = np.array([u('%s') % dt for dt in values[imask]])
        return values.tolist()

    def __array_finalize__(self, obj):
        if not self.ndim:  # pragma: no cover
            return self.item()

        self.freq = getattr(obj, 'freq', None)
        self.name = getattr(obj, 'name', None)
        self._reset_identity()

    def __repr__(self):
        output = com.pprint_thing(self.__class__) + '\n'
        output += 'freq: %s\n' % self.freq
        n = len(self)
        if n == 1:
            output += '[%s]\n' % (self[0])
        elif n == 2:
            output += '[%s, %s]\n' % (self[0], self[-1])
        elif n:
            output += '[%s, ..., %s]\n' % (self[0], self[-1])
        output += 'length: %d' % n
        return output

    def __unicode__(self):
        output = self.__class__.__name__
        output += u('(')
        prefix = '' if compat.PY3 else 'u'
        mapper = "{0}'{{0}}'".format(prefix)
        output += '[{0}]'.format(', '.join(map(mapper.format, self)))
        output += ", freq='{0}'".format(self.freq)
        output += ')'
        return output

    def __bytes__(self):
        encoding = com.get_option('display.encoding')
        return self.__unicode__().encode(encoding, 'replace')

    def __str__(self):
        if compat.PY3:
            return self.__unicode__()
        return self.__bytes__()

    def take(self, indices, axis=None):
        """
        Analogous to ndarray.take
        """
        indices = com._ensure_platform_int(indices)
        taken = self.values.take(indices, axis=axis)
        taken = taken.view(PeriodIndex)
        taken.freq = self.freq
        taken.name = self.name
        return taken

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
                to_concat = [x.asobject for x in to_concat]
            else:
                cat_values = np.concatenate([x.values for x in to_concat])
                return PeriodIndex(cat_values, freq=self.freq, name=name)

        to_concat = [x.values if isinstance(x, Index) else x
                     for x in to_concat]
        return Index(com._concat_compat(to_concat), name=name)

    def __reduce__(self):
        """Necessary for making this object picklable"""
        object_state = list(np.ndarray.__reduce__(self))
        subclass_state = (self.name, self.freq)
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        """Necessary for making this object picklable"""
        if len(state) == 2:
            nd_state, own_state = state
            np.ndarray.__setstate__(self, nd_state)
            self.name = own_state[0]
            try:  # backcompat
                self.freq = own_state[1]
            except:
                pass
        else:  # pragma: no cover
            np.ndarray.__setstate__(self, state)


def _get_ordinal_range(start, end, periods, freq):
    if com._count_not_none(start, end, periods) < 2:
        raise ValueError('Must specify 2 of start, end, periods')

    if start is not None:
        start = Period(start, freq)
    if end is not None:
        end = Period(end, freq)

    is_start_per = isinstance(start, Period)
    is_end_per = isinstance(end, Period)

    if is_start_per and is_end_per and (start.freq != end.freq):
        raise ValueError('Start and end must have same freq')

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
            base = FreqGroup.FR_QTR
        else:
            base, mult = _gfc(freq)
            if mult != 1:
                raise ValueError('Only mult == 1 supported')
            if base != FreqGroup.FR_QTR:
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
        if isinstance(x, (list, np.ndarray)):
            if length is not None and len(x) != length:
                raise ValueError('Mismatched Period array lengths')
            elif length is None:
                length = len(x)

    arrays = [np.asarray(x) if isinstance(x, (np.ndarray, list))
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

        mnum = _month_numbers[_freq_mod._get_rule_month(freq)] + 1
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
