# -*- coding: utf-8 -*-
import numpy as np

import pandas as pd

from pandas.core.arrays.datetimes import DatetimeArrayMixin
from pandas.core.arrays.timedelta import TimedeltaArrayMixin
from pandas.core.arrays.period import PeriodArrayMixin


class TestDatetimeArray(object):

    def test_from_dti(self, tz_naive_fixture):
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        arr = DatetimeArrayMixin(dti)
        assert list(dti) == list(arr)

    def test_astype_object(self, tz_naive_fixture):
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        arr = DatetimeArrayMixin(dti)
        asobj = arr.astype('O')
        assert isinstance(asobj, np.ndarray)
        assert asobj.dtype == 'O'
        assert list(asobj) == list(dti)


class TestTimedeltaArray(object):
    def test_from_tdi(self):
        tdi = pd.TimedeltaIndex(['1 Day', '3 Hours'])
        arr = TimedeltaArrayMixin(tdi)
        assert list(arr) == list(tdi)

    def test_astype_object(self):
        tdi = pd.TimedeltaIndex(['1 Day', '3 Hours'])
        arr = TimedeltaArrayMixin(tdi)
        asobj = arr.astype('O')
        assert isinstance(asobj, np.ndarray)
        assert asobj.dtype == 'O'
        assert list(asobj) == list(tdi)


class TestPeriodArray(object):

    def test_from_pi(self):
        pi = pd.period_range('2016', freq='Q', periods=3)
        arr = PeriodArrayMixin(pi)
        assert list(arr) == list(pi)

    def test_astype_object(self):
        pi = pd.period_range('2016', freq='Q', periods=3)
        arr = PeriodArrayMixin(pi)
        asobj = arr.astype('O')
        assert isinstance(asobj, np.ndarray)
        assert asobj.dtype == 'O'
        assert list(asobj) == list(pi)
