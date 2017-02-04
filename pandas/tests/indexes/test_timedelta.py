import numpy as np
from datetime import timedelta

import pandas as pd
import pandas.util.testing as tm
from pandas import (timedelta_range, date_range, Series, Timedelta,
                    DatetimeIndex)


class TestSlicing(tm.TestCase):

    def test_timedelta(self):
        # this is valid too
        index = date_range('1/1/2000', periods=50, freq='B')
        shifted = index + timedelta(1)
        back = shifted + timedelta(-1)
        self.assertTrue(tm.equalContents(index, back))
        self.assertEqual(shifted.freq, index.freq)
        self.assertEqual(shifted.freq, back.freq)

        result = index - timedelta(1)
        expected = index + timedelta(-1)
        tm.assert_index_equal(result, expected)

        # GH4134, buggy with timedeltas
        rng = date_range('2013', '2014')
        s = Series(rng)
        result1 = rng - pd.offsets.Hour(1)
        result2 = DatetimeIndex(s - np.timedelta64(100000000))
        result3 = rng - np.timedelta64(100000000)
        result4 = DatetimeIndex(s - pd.offsets.Hour(1))
        tm.assert_index_equal(result1, result4)
        tm.assert_index_equal(result2, result3)


class TestTimeSeries(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_series_box_timedelta(self):
        rng = timedelta_range('1 day 1 s', periods=5, freq='h')
        s = Series(rng)
        tm.assertIsInstance(s[1], Timedelta)
        tm.assertIsInstance(s.iat[2], Timedelta)
