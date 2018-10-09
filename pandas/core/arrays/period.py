# -*- coding: utf-8 -*-
from datetime import timedelta
import operator
import warnings

import numpy as np

from pandas._libs import lib
from pandas._libs.tslib import NaT, iNaT
from pandas._libs.tslibs.period import (
    Period, IncompatibleFrequency, DIFFERENT_FREQ_INDEX,
    get_period_field_arr, period_asfreq_arr)
from pandas._libs.tslibs import period as libperiod
from pandas._libs.tslibs.timedeltas import delta_to_nanoseconds, Timedelta
from pandas._libs.tslibs.fields import isleapyear_arr

from pandas import compat
from pandas.util._decorators import (cache_readonly, deprecate_kwarg)

from pandas.core.dtypes.common import (
    is_integer_dtype, is_float_dtype, is_period_dtype, is_timedelta64_dtype,
    is_datetime64_dtype, _TD_DTYPE)
from pandas.core.dtypes.dtypes import PeriodDtype
from pandas.core.dtypes.generic import ABCSeries

import pandas.core.common as com

from pandas.tseries import frequencies
from pandas.tseries.offsets import Tick, DateOffset

from pandas.core.arrays.datetimelike import DatetimeLikeArrayMixin


def _field_accessor(name, alias, docstring=None):
    def f(self):
        base, mult = frequencies.get_freq_code(self.freq)
        result = get_period_field_arr(alias, self._ndarray_values, base)
        return result

    f.__name__ = name
    f.__doc__ = docstring
    return property(f)


def _period_array_cmp(cls, op):
    """
    Wrap comparison operations to convert Period-like to PeriodDtype
    """
    opname = '__{name}__'.format(name=op.__name__)
    nat_result = True if opname == '__ne__' else False

    def wrapper(self, other):
        op = getattr(self._ndarray_values, opname)
        if isinstance(other, Period):
            if other.freq != self.freq:
                msg = DIFFERENT_FREQ_INDEX.format(self.freqstr, other.freqstr)
                raise IncompatibleFrequency(msg)

            result = op(other.ordinal)
        elif isinstance(other, PeriodArrayMixin):
            if other.freq != self.freq:
                msg = DIFFERENT_FREQ_INDEX.format(self.freqstr, other.freqstr)
                raise IncompatibleFrequency(msg)

            result = op(other._ndarray_values)

            mask = self._isnan | other._isnan
            if mask.any():
                result[mask] = nat_result

            return result
        elif other is NaT:
            result = np.empty(len(self._ndarray_values), dtype=bool)
            result.fill(nat_result)
        else:
            other = Period(other, freq=self.freq)
            result = op(other.ordinal)

        if self.hasnans:
            result[self._isnan] = nat_result

        return result

    return compat.set_function_name(wrapper, opname, cls)


class PeriodArrayMixin(DatetimeLikeArrayMixin):
    @property
    def _box_func(self):
        return lambda x: Period._from_ordinal(ordinal=x, freq=self.freq)

    @cache_readonly
    def dtype(self):
        return PeriodDtype.construct_from_string(self.freq)

    @property
    def _ndarray_values(self):
        # Ordinals
        return self._data

    @property
    def asi8(self):
        return self._ndarray_values.view('i8')

    @property
    def freq(self):
        """Return the frequency object if it is set, otherwise None"""
        return self._freq

    @freq.setter
    def freq(self, value):
        msg = ('Setting {cls}.freq has been deprecated and will be '
               'removed in a future version; use {cls}.asfreq instead. '
               'The {cls}.freq setter is not guaranteed to work.')
        warnings.warn(msg.format(cls=type(self).__name__),
                      FutureWarning, stacklevel=2)
        self._freq = value

    # --------------------------------------------------------------------
    # Constructors

    _attributes = ["freq"]

    def __new__(cls, values, freq=None, **kwargs):
        if is_period_dtype(values):
            # PeriodArray, PeriodIndex
            if freq is not None and values.freq != freq:
                raise IncompatibleFrequency(freq, values.freq)
            freq = values.freq
            values = values.asi8

        elif is_datetime64_dtype(values):
            # TODO: what if it has tz?
            values = dt64arr_to_periodarr(values, freq)

        return cls._simple_new(values, freq, **kwargs)

    @classmethod
    def _simple_new(cls, values, freq=None, **kwargs):
        """
        Values can be any type that can be coerced to Periods.
        Ordinals in an ndarray are fastpath-ed to `_from_ordinals`
        """

        if not is_integer_dtype(values):
            values = np.array(values, copy=False)
            if len(values) > 0 and is_float_dtype(values):
                raise TypeError("{cls} can't take floats"
                                .format(cls=cls.__name__))
            return cls(values, freq=freq)

        return cls._from_ordinals(values, freq)

    @classmethod
    def _from_ordinals(cls, values, freq=None):
        """
        Values should be int ordinals
        `__new__` & `_simple_new` cooerce to ordinals and call this method
        """

        values = np.array(values, dtype='int64', copy=False)

        result = object.__new__(cls)
        result._data = values
        if freq is None:
            raise ValueError('freq is not specified and cannot be inferred')
        result._freq = Period._maybe_convert_freq(freq)
        return result

    @classmethod
    def _generate_range(cls, start, end, periods, freq, fields):
        if freq is not None:
            freq = Period._maybe_convert_freq(freq)

        field_count = len(fields)
        if com.count_not_none(start, end) > 0:
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

    # --------------------------------------------------------------------
    # Vectorized analogues of Period properties

    year = _field_accessor('year', 0, "The year of the period")
    month = _field_accessor('month', 3, "The month as January=1, December=12")
    day = _field_accessor('day', 4, "The days of the period")
    hour = _field_accessor('hour', 5, "The hour of the period")
    minute = _field_accessor('minute', 6, "The minute of the period")
    second = _field_accessor('second', 7, "The second of the period")
    weekofyear = _field_accessor('week', 8, "The week ordinal of the year")
    week = weekofyear
    dayofweek = _field_accessor('dayofweek', 10,
                                "The day of the week with Monday=0, Sunday=6")
    weekday = dayofweek
    dayofyear = day_of_year = _field_accessor('dayofyear', 9,
                                              "The ordinal day of the year")
    quarter = _field_accessor('quarter', 2, "The quarter of the date")
    qyear = _field_accessor('qyear', 1)
    days_in_month = _field_accessor('days_in_month', 11,
                                    "The number of days in the month")
    daysinmonth = days_in_month

    @property
    def is_leap_year(self):
        """ Logical indicating if the date belongs to a leap year """
        return isleapyear_arr(np.asarray(self.year))

    @property
    def start_time(self):
        return self.to_timestamp(how='start')

    @property
    def end_time(self):
        return self.to_timestamp(how='end')

    def asfreq(self, freq=None, how='E'):
        """
        Convert the Period Array/Index to the specified frequency `freq`.

        Parameters
        ----------
        freq : str
            a frequency
        how : str {'E', 'S'}
            'E', 'END', or 'FINISH' for end,
            'S', 'START', or 'BEGIN' for start.
            Whether the elements should be aligned to the end
            or start within pa period. January 31st ('END') vs.
            January 1st ('START') for example.

        Returns
        -------
        new : Period Array/Index with the new frequency

        Examples
        --------
        >>> pidx = pd.period_range('2010-01-01', '2015-01-01', freq='A')
        >>> pidx
        <class 'pandas.core.indexes.period.PeriodIndex'>
        [2010, ..., 2015]
        Length: 6, Freq: A-DEC

        >>> pidx.asfreq('M')
        <class 'pandas.core.indexes.period.PeriodIndex'>
        [2010-12, ..., 2015-12]
        Length: 6, Freq: M

        >>> pidx.asfreq('M', how='S')
        <class 'pandas.core.indexes.period.PeriodIndex'>
        [2010-01, ..., 2015-01]
        Length: 6, Freq: M
        """
        how = libperiod._validate_end_alias(how)

        freq = Period._maybe_convert_freq(freq)

        base1, mult1 = frequencies.get_freq_code(self.freq)
        base2, mult2 = frequencies.get_freq_code(freq)

        asi8 = self.asi8
        # mult1 can't be negative or 0
        end = how == 'E'
        if end:
            ordinal = asi8 + mult1 - 1
        else:
            ordinal = asi8

        new_data = period_asfreq_arr(ordinal, base1, base2, end)

        if self.hasnans:
            new_data[self._isnan] = iNaT

        return self._shallow_copy(new_data, freq=freq)

    def to_timestamp(self, freq=None, how='start'):
        """
        Cast to DatetimeArray/Index

        Parameters
        ----------
        freq : string or DateOffset, optional
            Target frequency. The default is 'D' for week or longer,
            'S' otherwise
        how : {'s', 'e', 'start', 'end'}

        Returns
        -------
        DatetimeArray/Index
        """
        from pandas.core.arrays.datetimes import DatetimeArrayMixin

        how = libperiod._validate_end_alias(how)

        end = how == 'E'
        if end:
            if freq == 'B':
                # roll forward to ensure we land on B date
                adjust = Timedelta(1, 'D') - Timedelta(1, 'ns')
                return self.to_timestamp(how='start') + adjust
            else:
                adjust = Timedelta(1, 'ns')
                return (self + 1).to_timestamp(how='start') - adjust

        if freq is None:
            base, mult = frequencies.get_freq_code(self.freq)
            freq = frequencies.get_to_timestamp_base(base)
        else:
            freq = Period._maybe_convert_freq(freq)

        base, mult = frequencies.get_freq_code(freq)
        new_data = self.asfreq(freq, how=how)

        new_data = libperiod.periodarr_to_dt64arr(new_data._ndarray_values,
                                                  base)
        return DatetimeArrayMixin(new_data, freq='infer')

    # ------------------------------------------------------------------
    # Arithmetic Methods

    _create_comparison_method = classmethod(_period_array_cmp)

    def _sub_datelike(self, other):
        assert other is not NaT
        return NotImplemented

    def _sub_period(self, other):
        # If the operation is well-defined, we return an object-Index
        # of DateOffsets.  Null entries are filled with pd.NaT
        if self.freq != other.freq:
            msg = DIFFERENT_FREQ_INDEX.format(self.freqstr, other.freqstr)
            raise IncompatibleFrequency(msg)

        asi8 = self.asi8
        new_data = asi8 - other.ordinal
        new_data = np.array([self.freq * x for x in new_data])

        if self.hasnans:
            new_data[self._isnan] = NaT

        return new_data

    def _add_offset(self, other):
        assert not isinstance(other, Tick)
        base = frequencies.get_base_alias(other.rule_code)
        if base != self.freq.rule_code:
            msg = DIFFERENT_FREQ_INDEX.format(self.freqstr, other.freqstr)
            raise IncompatibleFrequency(msg)
        return self._time_shift(other.n)

    def _add_delta_td(self, other):
        assert isinstance(self.freq, Tick)  # checked by calling function
        assert isinstance(other, (timedelta, np.timedelta64, Tick))

        delta = self._check_timedeltalike_freq_compat(other)

        # Note: when calling parent class's _add_delta_td, it will call
        #  delta_to_nanoseconds(delta).  Because delta here is an integer,
        #  delta_to_nanoseconds will return it unchanged.
        return DatetimeLikeArrayMixin._add_delta_td(self, delta)

    def _add_delta_tdi(self, other):
        assert isinstance(self.freq, Tick)  # checked by calling function

        delta = self._check_timedeltalike_freq_compat(other)
        return self._addsub_int_array(delta, operator.add)

    def _add_delta(self, other):
        """
        Add a timedelta-like, Tick, or TimedeltaIndex-like object
        to self.

        Parameters
        ----------
        other : {timedelta, np.timedelta64, Tick,
                 TimedeltaIndex, ndarray[timedelta64]}

        Returns
        -------
        result : same type as self
        """
        if not isinstance(self.freq, Tick):
            # We cannot add timedelta-like to non-tick PeriodArray
            raise IncompatibleFrequency("Input has different freq from "
                                        "{cls}(freq={freqstr})"
                                        .format(cls=type(self).__name__,
                                                freqstr=self.freqstr))

        # TODO: standardize across datetimelike subclasses whether to return
        #  i8 view or _shallow_copy
        if isinstance(other, (Tick, timedelta, np.timedelta64)):
            new_values = self._add_delta_td(other)
            return self._shallow_copy(new_values)
        elif is_timedelta64_dtype(other):
            # ndarray[timedelta64] or TimedeltaArray/index
            new_values = self._add_delta_tdi(other)
            return self._shallow_copy(new_values)
        else:  # pragma: no cover
            raise TypeError(type(other).__name__)

    @deprecate_kwarg(old_arg_name='n', new_arg_name='periods')
    def shift(self, periods):
        """
        Shift index by desired number of increments.

        This method is for shifting the values of period indexes
        by a specified time increment.

        Parameters
        ----------
        periods : int
            Number of periods (or increments) to shift by,
            can be positive or negative.

            .. versionchanged:: 0.24.0

        Returns
        -------
        pandas.PeriodIndex
            Shifted index.

        See Also
        --------
        DatetimeIndex.shift : Shift values of DatetimeIndex.
        """
        return self._time_shift(periods)

    def _time_shift(self, n):
        values = self._ndarray_values + n * self.freq.n
        if self.hasnans:
            values[self._isnan] = iNaT
        return self._shallow_copy(values=values)

    def _maybe_convert_timedelta(self, other):
        """
        Convert timedelta-like input to an integer multiple of self.freq

        Parameters
        ----------
        other : timedelta, np.timedelta64, DateOffset, int,
                np.ndarray[timedelta64], TimedeltaIndex

        Returns
        -------
        converted : int, np.ndarray[int64]

        Raises
        ------
        IncompatibleFrequency : if the input cannot be written as a multiple
            of self.freq.  Note IncompatibleFrequency subclasses ValueError.
        """
        if isinstance(
                other, (timedelta, np.timedelta64, Tick, np.ndarray)):
            # TODO: is the np.ndarray case still relevant now that Arithmetic
            #  methods don't call this method?
            offset = frequencies.to_offset(self.freq.rule_code)
            if isinstance(offset, Tick):
                # _check_timedeltalike_freq_compat will raise if incompatible
                delta = self._check_timedeltalike_freq_compat(other)
                return delta
        elif isinstance(other, DateOffset):
            freqstr = other.rule_code
            base = frequencies.get_base_alias(freqstr)
            if base == self.freq.rule_code:
                return other.n
            msg = DIFFERENT_FREQ_INDEX.format(self.freqstr, other.freqstr)
            raise IncompatibleFrequency(msg)
        elif lib.is_integer(other):
            # integer is passed to .shift via
            # _add_datetimelike_methods basically
            # but ufunc may pass integer to _add_delta
            return other

        # raise when input doesn't have freq
        msg = "Input has different freq from {cls}(freq={freqstr})"
        raise IncompatibleFrequency(msg.format(cls=type(self).__name__,
                                               freqstr=self.freqstr))

    def _check_timedeltalike_freq_compat(self, other):
        assert isinstance(self.freq, Tick)  # checked by calling function
        own_offset = frequencies.to_offset(self.freq.rule_code)
        base_nanos = delta_to_nanoseconds(own_offset)

        if isinstance(other, (timedelta, np.timedelta64, Tick)):
            nanos = delta_to_nanoseconds(other)

        elif isinstance(other, np.ndarray):
            assert other.dtype.kind == 'm'
            if other.dtype != _TD_DTYPE:
                # i.e. non-nano unit
                # TODO: disallow unit-less timedelta64
                other = other.astype(_TD_DTYPE)
            nanos = other.view('i8')
        else:
            # TimedeltaArray/Index
            nanos = other.asi8

        if np.all(nanos % base_nanos == 0):
            # nanos being added is an integer multiple of the
            #  base-frequency to self.freq
            delta = nanos // base_nanos
            # delta is the integer (or integer-array) number of periods
            # by which will be added to self.
            return delta

        raise IncompatibleFrequency("Input has different freq from "
                                    "{cls}(freq={freqstr})"
                                    .format(cls=type(self).__name__,
                                            freqstr=self.freqstr))


PeriodArrayMixin._add_comparison_ops()
PeriodArrayMixin._add_datetimelike_methods()


# -------------------------------------------------------------------
# Constructor Helpers

def dt64arr_to_periodarr(data, freq, tz=None):
    if data.dtype != np.dtype('M8[ns]'):
        raise ValueError('Wrong dtype: %s' % data.dtype)

    freq = Period._maybe_convert_freq(freq)
    base, mult = frequencies.get_freq_code(freq)
    return libperiod.dt64arr_to_periodarr(data.view('i8'), base, tz)


def _get_ordinal_range(start, end, periods, freq, mult=1):
    if com.count_not_none(start, end, periods) != 2:
        raise ValueError('Of the three parameters: start, end, and periods, '
                         'exactly two must be specified')

    if freq is not None:
        _, mult = frequencies.get_freq_code(freq)

    if start is not None:
        start = Period(start, freq)
    if end is not None:
        end = Period(end, freq)

    is_start_per = isinstance(start, Period)
    is_end_per = isinstance(end, Period)

    if is_start_per and is_end_per and start.freq != end.freq:
        raise ValueError('start and end must have same freq')
    if (start is NaT or end is NaT):
        raise ValueError('start and end must not be NaT')

    if freq is None:
        if is_start_per:
            freq = start.freq
        elif is_end_per:
            freq = end.freq
        else:  # pragma: no cover
            raise ValueError('Could not infer freq from start/end')

    if periods is not None:
        periods = periods * mult
        if start is None:
            data = np.arange(end.ordinal - periods + mult,
                             end.ordinal + 1, mult,
                             dtype=np.int64)
        else:
            data = np.arange(start.ordinal, start.ordinal + periods, mult,
                             dtype=np.int64)
    else:
        data = np.arange(start.ordinal, end.ordinal + 1, mult, dtype=np.int64)

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
            base, mult = frequencies.get_freq_code(freq)
            if base != frequencies.FreqGroup.FR_QTR:
                raise AssertionError("base must equal FR_QTR")

        year, quarter = _make_field_arrays(year, quarter)
        for y, q in compat.zip(year, quarter):
            y, m = libperiod.quarter_to_myear(y, q, freq)
            val = libperiod.period_ordinal(y, m, 1, 1, 1, 1, 0, 0, base)
            ordinals.append(val)
    else:
        base, mult = frequencies.get_freq_code(freq)
        arrays = _make_field_arrays(year, month, day, hour, minute, second)
        for y, mth, d, h, mn, s in compat.zip(*arrays):
            ordinals.append(libperiod.period_ordinal(
                y, mth, d, h, mn, s, 0, 0, base))

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
