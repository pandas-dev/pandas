# -*- coding: utf-8 -*-
from datetime import timedelta
import warnings

import numpy as np

from pandas._libs import lib
from pandas._libs.tslib import NaT
from pandas._libs.tslibs.period import (
    Period, IncompatibleFrequency, DIFFERENT_FREQ_INDEX,
    get_period_field_arr)
from pandas._libs.tslibs.timedeltas import delta_to_nanoseconds
from pandas._libs.tslibs.fields import isleapyear_arr

from pandas.util._decorators import cache_readonly

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
               'removed in a future version; use PeriodIndex.asfreq instead. '
               'The {cls}.freq setter is not guaranteed to work.')
        warnings.warn(msg.format(cls=type(self).__name__),
                      FutureWarning, stacklevel=2)
        self._freq = value

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
