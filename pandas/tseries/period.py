from datetime import datetime
import numpy as np

from pandas.tseries.frequencies import get_freq_code as _gfc, to_offset
from pandas.tseries.index import DatetimeIndex, Int64Index
from pandas.tseries.tools import parse_time_string
import pandas.tseries.frequencies as _freq_mod

import pandas.core.common as com
from pandas.util import py3compat

from pandas._tseries import Timestamp
import pandas._tseries as lib


#---------------
# Period logic


def _period_field_accessor(name, alias=None):
    if alias is None:
        alias = name
    def f(self):
        base, mult = _gfc(self.freq)
        g = getattr(lib, 'get_period_%s' % alias)
        return g(self.ordinal, base, mult)
    f.__name__ = name
    return property(f)

def _field_accessor(name, alias=None):
    if alias is None:
        alias = name
    def f(self):
        base, mult = _gfc(self.freq)
        g = getattr(lib, 'get_period_%s_arr' % alias)
        return g(self.ordinal, base, mult)
    f.__name__ = name
    return property(f)

def to_period(arg, freq=None):
    """ Attempts to convert arg to timestamp """
    if arg is None:
        return arg

    if type(arg) == float:
        raise TypeError("Cannot convert a float to period")

    return Period(arg, freq=freq)

class Period(object):

    def __init__(self, value=None, freq=None,
                 year=None, month=1, quarter=None, day=1,
                 hour=0, minute=0, second=0):
        """
        Represents an period of time

        Parameters
        ----------
        value : Period or basestring, default None
            The time period represented (e.g., '4Q2005')
        freq : str, default None
            e.g., 'B' for businessday, ('T', 5) or '5T' for 5 minutes
        year : int, default None
        month : int, default 1
        quarter : int, default None
        day : int, default 1
        hour : int, default 0
        minute : int, default 0
        second : int, default 0
        """
        # freq points to a tuple (base, mult);  base is one of the defined
        # periods such as A, Q, etc. Every five minutes would be, e.g.,
        # ('T', 5) but may be passed in as a string like '5T'

        self.freq = None

        # ordinal is the period offset from the gregorian proleptic epoch

        self.ordinal = None

        if value is None:
            if freq is None:
                raise ValueError("If value is None, freq cannot be None")

            if year is None:
                raise ValueError("If value is None, year cannot be None")

            if quarter is not None:
                month = (quarter - 1) * 3 + 1

            base, mult = _gfc(freq)

            self.ordinal = lib.period_ordinal(year, month, day, hour, minute,
                                            second, base, mult)

        elif isinstance(value, Period):
            other = value
            if freq is None or _gfc(freq) == _gfc(other.freq):
                self.ordinal = other.ordinal
                freq = other.freq
            else:
                converted = other.asfreq(freq)
                self.ordinal = converted.ordinal

        elif isinstance(value, basestring):
            value = value.upper()
            dt, parsed, reso = parse_time_string(value)

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
                else:
                    raise ValueError("Could not infer frequency for period")

        elif isinstance(value, datetime):
            dt = value
            if freq is None:
                raise ValueError('Must supply freq for datetime value')
        elif com.is_integer(value):
            if value <= 0:
                raise ValueError("Value must be positive")
            self.ordinal = value
            if freq is None:
                raise ValueError('Must supply freq for ordinal value')
        else:
            msg = "Value must be Period, string, integer, or datetime"
            raise ValueError(msg)

        base, mult = _gfc(freq)

        if self.ordinal is None:
            self.ordinal = lib.period_ordinal(dt.year, dt.month, dt.day, dt.hour,
                                            dt.minute, dt.second, base, mult)

        self.freq = _freq_mod._get_freq_str(base, mult)

    def __eq__(self, other):
        if isinstance(other, Period):
            return (self.ordinal == other.ordinal
                    and _gfc(self.freq) == _gfc(other.freq))
        return False

    def __add__(self, other):
        if isinstance(other, (int, long)):
            return Period(self.ordinal + other, self.freq)
        raise ValueError("Cannot add with non-integer value")

    def __sub__(self, other):
        if isinstance(other, (int, long)):
            return Period(self.ordinal - other, self.freq)
        if isinstance(other, Period):
            if other.freq != self.freq:
                raise ValueError("Cannot do arithmetic with "
                                 "non-conforming periods")
            return self.ordinal - other.ordinal
        raise ValueError("Cannot sub with non-integer value")

    def asfreq(self, freq=None, how='E'):
        """

        Parameters
        ----------
        freq :
        how :

        Returns
        -------
        resampled : Period
        """
        how = _validate_end_alias(how)
        base1, mult1 = _gfc(self.freq)
        base2, mult2 = _gfc(freq)

        new_ordinal = lib.period_asfreq(self.ordinal, base1, mult1,
                                        base2, mult2, py3compat.str_to_bytes(how))

        return Period(new_ordinal, (base2, mult2))

    @property
    def start_time(self):
        return self.to_timestamp(how='S')

    @property
    def end_time(self):
        return self.to_timestamp(how='E')

    def to_timestamp(self, freq='D', how='S'):
        """
        Return the Timestamp at the start/end of the period

        Parameters
        ----------
        freq : string or DateOffset, default 'D'
            Target frequency
        how: str, default 'S' (start)
            'S', 'E'. Can be aliased as case insensitive
            'Start', 'Finish', 'Begin', 'End'

        Returns
        -------
        Timestamp
        """
        # how = _validate_end_alias(how)

        base, mult = _gfc(freq)
        new_val = self.asfreq(freq, how)
        new_val = lib.period_ordinal_to_dt64(new_val.ordinal, base, mult)
        ts_freq = _period_rule_to_timestamp_rule(self.freq, how=how)
        return Timestamp(new_val, offset=to_offset(ts_freq))

    year = _period_field_accessor('year')
    month = _period_field_accessor('month')
    day = _period_field_accessor('day')
    hour = _period_field_accessor('hour')
    minute = _period_field_accessor('minute')
    second = _period_field_accessor('second')
    weekofyear = _period_field_accessor('week')
    week = weekofyear
    dayofweek = _period_field_accessor('dayofweek', 'dow')
    weekday = dayofweek
    dayofyear = day_of_year = _period_field_accessor('dayofyear', 'doy')
    quarter = _period_field_accessor('quarter')
    qyear = _period_field_accessor('qyear')

    @classmethod
    def now(cls, freq=None):
        return Period(datetime.now(), freq=freq)

    def __repr__(self):
        base, mult = _gfc(self.freq)
        formatted = lib.period_ordinal_to_string(self.ordinal, base, mult)
        freqstr = _freq_mod._reverse_period_code_map[base]
        if mult == 1:
            return "Period('%s', '%s')" % (formatted, freqstr)
        return ("Period('%s', '%d%s')" % (formatted, mult, freqstr))

    def __str__(self):
        base, mult = _gfc(self.freq)
        formatted = lib.period_ordinal_to_string(self.ordinal, base, mult)
        return ("%s" % formatted)

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
        if fmt is not None:
            return lib.period_strftime(self.ordinal, base, mult, fmt)
        else:
            return lib.period_ordinal_to_string(self.ordinal, base, mult)

def _period_unbox(key, check=None):
    '''
    Period-like => int64
    '''
    if not isinstance(key, Period):
        key = Period(key, freq=check)
    elif check is not None:
        if key.freq != check:
            raise ValueError("%s is wrong freq" % key)
    return np.int64(key.ordinal)

def _period_unbox_array(arr, check=None):
    if arr is None:
        return arr
    unboxer = np.frompyfunc(lambda x: _period_unbox(x, check=check), 1, 1)
    return unboxer(arr)

def _period_box_array(arr, freq):
    if arr is None:
        return arr

    if not isinstance(arr, np.ndarray):
        return arr

    boxfunc = lambda x: Period(x, freq)
    boxer = np.frompyfunc(boxfunc, 1, 1)
    return boxer(arr)

def dt64arr_to_periodarr(data, freq):
    if data is None:
        return data

    if isinstance(freq, basestring):
        base, mult = _gfc(freq)
    else:
        base, mult = freq

    return lib.dt64arr_to_periodarr(data.view('i8'), base, mult)

# --- Period index sketch

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
    """

    def __new__(cls, data=None,
                freq=None, start=None, end=None, periods=None,
                copy=False, name=None):

        if isinstance(freq, Period):
            freq = freq.freq
        else:
            freq = _freq_mod.get_standard_freq(freq)

        if data is None:
            subarr, freq = _get_ordinal_range(start, end, periods, freq)
            subarr = subarr.view(cls)
            subarr.name = name
            subarr.freq = freq

            return subarr

        if not isinstance(data, np.ndarray):
            if np.isscalar(data):
                raise ValueError('PeriodIndex() must be called with a '
                                 'collection of some kind, %s was passed'
                                 % repr(data))

            elif isinstance(data, Period):
                raise ValueError('Data must be array of dates, strings, '
                                 'or Period objects')

            # other iterable of some kind
            if not isinstance(data, (list, tuple)):
                data = list(data)

            try:
                data = np.array(data, dtype='i8')
            except:
                data = np.array(data, dtype='O')

            if freq is None and len(data) > 0:
                freq = getattr(data[0], 'freq')

            if freq is None:
                raise ValueError(('freq not specified and cannot be inferred '
                                  'from first element'))

            data = _period_unbox_array(data, check=freq)
        else:
            if isinstance(data, PeriodIndex):
                if freq is None or freq == data.freq:
                    freq = data.freq
                    data = data.values
                else:
                    base1, mult1 = _gfc(data.freq)
                    base2, mult2 = _gfc(freq)
                    data = lib.period_asfreq_arr(data.values, base1, mult1,
                                                 base2, mult2, b'E')
            else:
                if freq is None and len(data) > 0:
                    freq = getattr(data[0], 'freq')

                if freq is None:
                    raise ValueError(('freq not specified and cannot be '
                                      'inferred from first element'))

                if data.dtype == np.datetime64:
                    data = dt64arr_to_periodarr(data, freq)
                elif data.dtype == np.int64:
                    pass
                else:
                    try:
                        data = data.astype('i8')
                    except:
                        data = data.astype('O')
                        data = _period_unbox_array(data, check=freq)

        data = np.array(data, dtype=np.int64, copy=False)

        if (data <= 0).any():
            raise ValueError("Found illegal (<= 0) values in data")

        subarr = data.view(cls)
        subarr.name = name
        subarr.freq = freq

        return subarr

    def __iter__(self):
        for val in self.values:
            yield Period(val, freq=self.freq)

    @property
    def is_all_dates(self):
        return True

    @property
    def freqstr(self):
        return self.freq

    def asfreq(self, freq=None, how='E'):
        how = _validate_end_alias(how)

        base1, mult1 = _gfc(self.freq)

        if isinstance(freq, basestring):
            base2, mult2 = _gfc(freq)
        else:
            base2, mult2 = freq


        new_data = lib.period_asfreq_arr(self.values,
                                         base1, mult1,
                                         base2, mult2, py3compat.str_to_bytes(how))

        return PeriodIndex(new_data, freq=freq)

    year = _period_field_accessor('year')
    month = _period_field_accessor('month')
    day = _period_field_accessor('day')
    hour = _period_field_accessor('hour')
    minute = _period_field_accessor('minute')
    second = _period_field_accessor('second')
    weekofyear = _period_field_accessor('week')
    week = weekofyear
    dayofweek = _period_field_accessor('dayofweek', 'dow')
    weekday = dayofweek
    dayofyear = day_of_year = _period_field_accessor('dayofyear', 'doy')
    quarter = _period_field_accessor('quarter')
    qyear = _period_field_accessor('qyear')

    # Try to run function on index first, and then on elements of index
    # Especially important for group-by functionality
    def map(self, func_to_map):
        try:
            return func_to_map(self)
        except:
            return super(DatetimeIndex, self).map(func_to_map)

    def _mpl_repr(self):
        # how to represent ourselves to matplotlib
        return _period_box_array(self, self.freq)

    def to_timestamp(self, freq='D', how='start'):
        """
        Cast to DatetimeIndex

        Parameters
        ----------
        freq : string or DateOffset, default 'D'
            Target frequency
        how : {'s', 'e', 'start', 'end'}

        Returns
        -------
        DatetimeIndex
        """
        base, mult = _gfc(freq)
        new_data = self.asfreq(freq, how)
        new_data = lib.periodarr_to_dt64arr(new_data.values, base, mult)
        return DatetimeIndex(new_data, freq='infer')

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
        if isinstance(other, (int, long)):
            return PeriodIndex(self.values + other, self.freq)
        return super(PeriodIndex, self).__add__(other)

    def __sub__(self, other):
        if isinstance(other, (int, long)):
            return PeriodIndex(self.values - other, self.freq)
        if isinstance(other, Period):
            if other.freq != self.freq:
                raise ValueError("Cannot do arithmetic with "
                                 "non-conforming periods")
            return PeriodIndex(self.values - other.ordinal)
        return super(PeriodIndex, self).__sub__(other)

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
        try:
            return super(PeriodIndex, self).get_value(series, key)
        except KeyError:
            try:
                asdt, parsed, reso = parse_time_string(key)
                grp = _freq_mod._infer_period_group(reso)
                freqn = _freq_mod._period_group(self.freq)

                # if our data is higher resolution than requested key, slice
                if grp < freqn:
                    iv = Period(asdt, freq=(grp,1))
                    ord1 = iv.asfreq(self.freq, how='S').ordinal
                    ord2 = iv.asfreq(self.freq, how='E').ordinal
                    pos = np.searchsorted(self.values, [ord1, ord2])
                    key = slice(pos[0], pos[1]+1)
                    return series[key]
                else:
                    key = to_period(asdt, freq=self.freq).ordinal
                    return self._engine.get_value(series, key)
            except TypeError:
                pass
            except KeyError:
                pass
            except IndexError:
                ival = Period(key, freq=self.freq)
                raise IndexError("%s is out of bounds" % ival)

            key = to_period(key, self.freq).ordinal
            return self._engine.get_value(series, key)

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
                asdt, parsed, reso = parse_time_string(key)
                key = asdt
            except TypeError:
                pass
            except KeyError:
                pass

            key = to_period(key, self.freq).ordinal
            return self._engine.get_loc(key)

    def __getitem__(self, key):
        """Override numpy.ndarray's __getitem__ method to work as desired"""
        arr_idx = self.view(np.ndarray)
        if np.isscalar(key):
            val = arr_idx[key]
            return Period(val, freq=self.freq)
        else:
            if com._is_bool_indexer(key):
                key = np.asarray(key)

            result = arr_idx[key]
            if result.ndim > 1:
                return PeriodIndex(result, name=self.name, freq=self.freq)

            return PeriodIndex(result, name=self.name, freq=self.freq)

    def format(self, name=False):
        """
        Render a string representation of the Index
        """
        header = []

        if name:
            header.append(str(self.name) if self.name is not None else '')

        return header + ['%s' % Period(x, freq=self.freq) for x in self]

    def _view_like(self, ndarray):
        result = ndarray.view(type(self))
        result.freq = self.freq
        result.name = self.name
        return result

    def __array_finalize__(self, obj):
        if self.ndim == 0: # pragma: no cover
            return self.item()

        self.freq = getattr(obj, 'freq', None)

    def __repr__(self):
        output = str(self.__class__) + '\n'
        output += 'freq: ''%s''\n' % self.freq
        if len(self) > 0:
            output += '[%s, ..., %s]\n' % (self[0], self[-1])
        output += 'length: %d' % len(self)
        return output

    def take(self, indices, axis=None):
        """
        Analogous to ndarray.take
        """
        taken = self.values.take(indices, axis=axis)
        taken = taken.view(PeriodIndex)
        taken.freq = self.freq
        taken.name = self.name
        return taken

def _get_ordinal_range(start, end, periods, freq):
    if com._count_not_none(start, end, periods) < 2:
        raise ValueError('Must specify 2 of start, end, periods')

    start = to_period(start, freq)
    end = to_period(end, freq)

    is_start_per = isinstance(start, Period)
    is_end_per = isinstance(end, Period)
    if (start is not None and not is_start_per):
        raise ValueError('Failed to convert %s to period' % start)

    if (end is not None and not is_end_per):
        raise ValueError('Failed to convert %s to period' % end)

    if is_start_per and is_end_per and (start.freq != end.freq):
        raise ValueError('Start and end must have same freq')

    if freq is None:
        if is_start_per:
            freq = start.freq
        elif is_end_per:
            freq = end.freq
        else:
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
        if start is None or end is None:
            msg = 'Must specify both start and end if periods is None'
            raise ValueError(msg)
        data = np.arange(start.ordinal, end.ordinal+1, dtype=np.int64)

    return data, freq


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

def period_range(start=None, end=None, periods=None, freq='D'):
    """
    Return a fixed frequency datetime index, with day (calendar) as the default
    frequency


    Parameters
    ----------
    start :
    end :
    normalize : bool, default False
        Normalize start/end dates to midnight before generating date range

    Returns
    -------

    """
    return PeriodIndex(start=start, end=end, periods=periods,
                       freq=freq)

def _period_rule_to_timestamp_rule(freq, how='end'):
    how = how.lower()
    if how in ('end', 'e'):
        return freq
    else:
        if freq.startswith('A-') or freq.startswith('BA-'):
            base, color = freq.split('-')
            return '%sS-%s' % (base, color)
        return freq
