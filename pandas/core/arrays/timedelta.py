# -*- coding: utf-8 -*-

from pandas._libs import tslib
from pandas._libs.tslib import Timedelta, NaT

from pandas.core.dtypes.common import _TD_DTYPE

from pandas.tseries.offsets import Tick

from .datetimelike import DatetimeLikeArrayMixin


class TimedeltaArrayMixin(DatetimeLikeArrayMixin):
    @property
    def _box_func(self):
        return lambda x: Timedelta(x, unit='ns')

    @property
    def dtype(self):
        return _TD_DTYPE

    # ----------------------------------------------------------------
    # Arithmetic Methods

    def _add_offset(self, other):
        assert not isinstance(other, Tick)
        raise TypeError("cannot add the type {typ} to a {cls}"
                        .format(typ=type(other).__name__,
                                cls=type(self).__name__))

    def _sub_datelike(self, other):
        assert other is not NaT
        raise TypeError("cannot subtract a datelike from a {cls}"
                        .format(cls=type(self).__name__))

    # ----------------------------------------------------------------
    # Conversion Methods - Vectorized analogues of Timedelta methods

    def to_pytimedelta(self):
        """
        Return Timedelta Array/Index as object ndarray of datetime.timedelta
        objects

        Returns
        -------
        datetimes : ndarray
        """
        return tslib.ints_to_pytimedelta(self.asi8)
