# -*- coding: utf-8 -*-

from pandas._libs.tslib import Timedelta

from pandas.core.dtypes.common import _TD_DTYPE

from .datetimelike import DatetimeLikeArrayMixin


class TimedeltaArray(DatetimeLikeArrayMixin):
    @property
    def _box_func(self):
        return lambda x: Timedelta(x, unit='ns')

    @property
    def dtype(self):
        return _TD_DTYPE
