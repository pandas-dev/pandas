# -*- coding: utf-8 -*-
from datetime import timedelta

import numpy as np

from pandas._libs.tslib import array_to_datetime
from pandas import Timestamp

import pandas.util.testing as tm


class TestArrayToDatetime(object):
    def test_coerce_out_of_bounds_utc(self):
        ts = Timestamp('1900-01-01', tz='US/Pacific')
        dt = ts.to_pydatetime() - timedelta(days=365 * 300)  # ~1600AD
        arr = np.array([dt])
        result = array_to_datetime(arr, utc=True, errors='coerce')
        expected = np.array(['NaT'], dtype='datetime64[ns]')
        tm.assert_numpy_array_equal(result, expected)
