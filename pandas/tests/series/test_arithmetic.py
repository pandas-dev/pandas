# -*- coding: utf-8 -*-
from datetime import timedelta

import pandas as pd
import pandas.util.testing as tm


class TestTimedeltaSeriesComparisons(object):
    def test_compare_timedelta_series(self):
        # regresssion test for GH5963
        s = pd.Series([timedelta(days=1), timedelta(days=2)])
        actual = s > timedelta(days=1)
        expected = pd.Series([False, True])
        tm.assert_series_equal(actual, expected)


class TestPeriodSeriesArithmetic(object):
    def test_ops_series_timedelta(self):
        # GH 13043
        ser = pd.Series([pd.Period('2015-01-01', freq='D'),
                         pd.Period('2015-01-02', freq='D')], name='xxx')
        assert ser.dtype == object

        expected = pd.Series([pd.Period('2015-01-02', freq='D'),
                              pd.Period('2015-01-03', freq='D')], name='xxx')

        result = ser + pd.Timedelta('1 days')
        tm.assert_series_equal(result, expected)

        result = pd.Timedelta('1 days') + ser
        tm.assert_series_equal(result, expected)

        result = ser + pd.tseries.offsets.Day()
        tm.assert_series_equal(result, expected)

        result = pd.tseries.offsets.Day() + ser
        tm.assert_series_equal(result, expected)

    def test_ops_series_period(self):
        # GH 13043
        ser = pd.Series([pd.Period('2015-01-01', freq='D'),
                         pd.Period('2015-01-02', freq='D')], name='xxx')
        assert ser.dtype == object

        per = pd.Period('2015-01-10', freq='D')
        # dtype will be object because of original dtype
        expected = pd.Series([9, 8], name='xxx', dtype=object)
        tm.assert_series_equal(per - ser, expected)
        tm.assert_series_equal(ser - per, -expected)

        s2 = pd.Series([pd.Period('2015-01-05', freq='D'),
                        pd.Period('2015-01-04', freq='D')], name='xxx')
        assert s2.dtype == object

        expected = pd.Series([4, 2], name='xxx', dtype=object)
        tm.assert_series_equal(s2 - ser, expected)
        tm.assert_series_equal(ser - s2, -expected)
