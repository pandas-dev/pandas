# -*- coding: utf-8 -*-
from datetime import timedelta
import warnings

import numpy as np

from pandas._libs import lib
from pandas._libs.tslib import NaT
from pandas._libs.tslibs.period import (
    Period, IncompatibleFrequency, DIFFERENT_FREQ_INDEX)
from pandas._libs.tslibs.timedeltas import delta_to_nanoseconds

from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.dtypes import PeriodDtype

from pandas.tseries import frequencies
from pandas.tseries.offsets import Tick, DateOffset

from .datetimelike import DatetimeLikeArrayMixin


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
        msg = ('Setting PeriodIndex.freq has been deprecated and will be '
               'removed in a future version; use PeriodIndex.asfreq instead. '
               'The PeriodIndex.freq setter is not guaranteed to work.')
        warnings.warn(msg, FutureWarning, stacklevel=2)
        self._freq = value

    # ------------------------------------------------------------------
    # Arithmetic Methods

    def _sub_datelike(self, other):
        assert other is not NaT
        return NotImplemented

    def _maybe_convert_timedelta(self, other):
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
        msg = "Input has different freq from PeriodIndex(freq={0})"
        raise IncompatibleFrequency(msg.format(self.freqstr))
