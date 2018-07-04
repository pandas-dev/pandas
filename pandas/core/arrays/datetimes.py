# -*- coding: utf-8 -*-
import warnings

import numpy as np
from pytz import utc

from pandas._libs.tslib import Timestamp, NaT, iNaT
from pandas._libs.tslibs import conversion, timezones

from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.common import _NS_DTYPE, is_datetime64tz_dtype
from pandas.core.dtypes.dtypes import DatetimeTZDtype

from .datetimelike import DatetimeLikeArrayMixin


class DatetimeArrayMixin(DatetimeLikeArrayMixin):
    """
    Assumes that subclass __new__/__init__ defines:
        tz
        _freq
        _data
    """

    # -----------------------------------------------------------------
    # Descriptive Properties

    @property
    def _box_func(self):
        return lambda x: Timestamp(x, freq=self.freq, tz=self.tz)

    @cache_readonly
    def dtype(self):
        if self.tz is None:
            return _NS_DTYPE
        return DatetimeTZDtype('ns', self.tz)

    @property
    def tzinfo(self):
        """
        Alias for tz attribute
        """
        return self.tz

    @property  # NB: override with cache_readonly in immutable subclasses
    def _timezone(self):
        """ Comparable timezone both for pytz / dateutil"""
        return timezones.get_timezone(self.tzinfo)

    @property
    def offset(self):
        """get/set the frequency of the instance"""
        msg = ('DatetimeIndex.offset has been deprecated and will be removed '
               'in a future version; use DatetimeIndex.freq instead.')
        warnings.warn(msg, FutureWarning, stacklevel=2)
        return self.freq

    @offset.setter
    def offset(self, value):
        """get/set the frequency of the instance"""
        msg = ('DatetimeIndex.offset has been deprecated and will be removed '
               'in a future version; use DatetimeIndex.freq instead.')
        warnings.warn(msg, FutureWarning, stacklevel=2)
        self.freq = value

    # -----------------------------------------------------------------
    # Comparison Methods

    def _has_same_tz(self, other):
        zzone = self._timezone

        # vzone sholdn't be None if value is non-datetime like
        if isinstance(other, np.datetime64):
            # convert to Timestamp as np.datetime64 doesn't have tz attr
            other = Timestamp(other)
        vzone = timezones.get_timezone(getattr(other, 'tzinfo', '__no_tz__'))
        return zzone == vzone

    def _assert_tzawareness_compat(self, other):
        # adapted from _Timestamp._assert_tzawareness_compat
        other_tz = getattr(other, 'tzinfo', None)
        if is_datetime64tz_dtype(other):
            # Get tzinfo from Series dtype
            other_tz = other.dtype.tz
        if other is NaT:
            # pd.NaT quacks both aware and naive
            pass
        elif self.tz is None:
            if other_tz is not None:
                raise TypeError('Cannot compare tz-naive and tz-aware '
                                'datetime-like objects.')
        elif other_tz is None:
            raise TypeError('Cannot compare tz-naive and tz-aware '
                            'datetime-like objects')

    # -----------------------------------------------------------------
    # Arithmetic Methods

    def _sub_datelike_dti(self, other):
        """subtraction of two DatetimeIndexes"""
        if not len(self) == len(other):
            raise ValueError("cannot add indices of unequal length")

        self_i8 = self.asi8
        other_i8 = other.asi8
        new_values = self_i8 - other_i8
        if self.hasnans or other.hasnans:
            mask = (self._isnan) | (other._isnan)
            new_values[mask] = iNaT
        return new_values.view('timedelta64[ns]')

    # -----------------------------------------------------------------
    # Timezone Conversion and Localization Methods

    def _local_timestamps(self):
        values = self.asi8
        indexer = values.argsort()
        result = conversion.tz_convert(values.take(indexer), utc, self.tz)

        n = len(indexer)
        reverse = np.empty(n, dtype=np.int_)
        reverse.put(indexer, np.arange(n))
        return result.take(reverse)
