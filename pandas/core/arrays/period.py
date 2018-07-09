# -*- coding: utf-8 -*-
from datetime import timedelta
import warnings

import numpy as np

from pandas._libs import lib
from pandas._libs.tslib import NaT, iNaT
from pandas._libs.tslibs.period import (
    Period, IncompatibleFrequency, DIFFERENT_FREQ_INDEX,
    get_period_field_arr, period_asfreq_arr)
from pandas._libs.tslibs import period as libperiod
from pandas._libs.tslibs.timedeltas import delta_to_nanoseconds
from pandas._libs.tslibs.fields import isleapyear_arr

from pandas import compat
from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.common import (
    is_integer_dtype, is_float_dtype, is_period_dtype)
from pandas.core.dtypes.dtypes import PeriodDtype

from pandas.tseries import frequencies
from pandas.tseries.offsets import Tick, DateOffset

from .datetimelike import DatetimeLikeArrayMixin


def _field_accessor(name, alias, docstring=None):
    def f(self):
        base, mult = frequencies.get_freq_code(self.freq)
        result = get_period_field_arr(alias, self._ndarray_values, base)
        return result

    f.__name__ = name
    f.__doc__ = docstring
    return property(f)


def _period_array_cmp(opname, cls):
    """
    Wrap comparison operations to convert Period-like to PeriodDtype
    """
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

        return self._simple_new(new_data, self.name, freq=freq)

    # ------------------------------------------------------------------
    # Arithmetic Methods

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
        return self.shift(other.n)

    def _add_delta_td(self, other):
        assert isinstance(other, (timedelta, np.timedelta64, Tick))
        nanos = delta_to_nanoseconds(other)
        own_offset = frequencies.to_offset(self.freq.rule_code)

        if isinstance(own_offset, Tick):
            offset_nanos = delta_to_nanoseconds(own_offset)
            if np.all(nanos % offset_nanos == 0):
                return self.shift(nanos // offset_nanos)

        # raise when input doesn't have freq
        raise IncompatibleFrequency("Input has different freq from "
                                    "{cls}(freq={freqstr})"
                                    .format(cls=type(self).__name__,
                                            freqstr=self.freqstr))

    def _add_delta(self, other):
        ordinal_delta = self._maybe_convert_timedelta(other)
        return self.shift(ordinal_delta)

    def shift(self, n):
        """
        Specialized shift which produces an Period Array/Index

        Parameters
        ----------
        n : int
            Periods to shift by

        Returns
        -------
        shifted : Period Array/Index
        """
        values = self._ndarray_values + n * self.freq.n
        if self.hasnans:
            values[self._isnan] = iNaT
        return self._shallow_copy(values=values)

    def _maybe_convert_timedelta(self, other):
        """
        Convert timedelta-like input to an integer multiple of self.freq

        Parameters
        ----------
        other : timedelta, np.timedelta64, DateOffset, int, np.ndarray

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
            offset = frequencies.to_offset(self.freq.rule_code)
            if isinstance(offset, Tick):
                if isinstance(other, np.ndarray):
                    nanos = np.vectorize(delta_to_nanoseconds)(other)
                else:
                    nanos = delta_to_nanoseconds(other)
                offset_nanos = delta_to_nanoseconds(offset)
                check = np.all(nanos % offset_nanos == 0)
                if check:
                    return nanos // offset_nanos
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

    @classmethod
    def _add_comparison_methods(cls):
        """ add in comparison methods """
        cls.__eq__ = _period_array_cmp('__eq__', cls)
        cls.__ne__ = _period_array_cmp('__ne__', cls)
        cls.__lt__ = _period_array_cmp('__lt__', cls)
        cls.__gt__ = _period_array_cmp('__gt__', cls)
        cls.__le__ = _period_array_cmp('__le__', cls)
        cls.__ge__ = _period_array_cmp('__ge__', cls)


PeriodArrayMixin._add_comparison_methods()
