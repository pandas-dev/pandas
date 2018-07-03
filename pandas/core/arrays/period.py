# -*- coding: utf-8 -*-

from pandas._libs.tslibs.period import Period

from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.dtypes import PeriodDtype

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
