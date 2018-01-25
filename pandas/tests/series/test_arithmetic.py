# -*- coding: utf-8 -*-
from datetime import timedelta
import operator

import numpy as np

import pandas as pd
import pandas.util.testing as tm


class TestSeriesComparison(object):
    def test_compare_invalid(self):
        # GH#8058
        # ops testing
        a = pd.Series(np.random.randn(5), name=0)
        b = pd.Series(np.random.randn(5))
        b.name = pd.Timestamp('2000-01-01')
        tm.assert_series_equal(a / b, 1 / (b / a))


class TestTimestampSeriesComparison(object):
    def test_timestamp_compare_series(self):
        # make sure we can compare Timestamps on the right AND left hand side
        # GH#4982
        ser = pd.Series(pd.date_range('20010101', periods=10), name='dates')
        s_nat = ser.copy(deep=True)

        ser[0] = pd.Timestamp('nat')
        ser[3] = pd.Timestamp('nat')

        ops = {'lt': 'gt', 'le': 'ge', 'eq': 'eq', 'ne': 'ne'}

        for left, right in ops.items():
            left_f = getattr(operator, left)
            right_f = getattr(operator, right)

            # no nats
            expected = left_f(ser, pd.Timestamp('20010109'))
            result = right_f(pd.Timestamp('20010109'), ser)
            tm.assert_series_equal(result, expected)

            # nats
            expected = left_f(ser, pd.Timestamp('nat'))
            result = right_f(pd.Timestamp('nat'), ser)
            tm.assert_series_equal(result, expected)

            # compare to timestamp with series containing nats
            expected = left_f(s_nat, pd.Timestamp('20010109'))
            result = right_f(pd.Timestamp('20010109'), s_nat)
            tm.assert_series_equal(result, expected)

            # compare to nat with series containing nats
            expected = left_f(s_nat, pd.Timestamp('nat'))
            result = right_f(pd.Timestamp('nat'), s_nat)
            tm.assert_series_equal(result, expected)

    def test_timestamp_equality(self):
        # GH#11034
        ser = pd.Series([pd.Timestamp('2000-01-29 01:59:00'), 'NaT'])
        result = ser != ser
        tm.assert_series_equal(result, pd.Series([False, True]))
        result = ser != ser[0]
        tm.assert_series_equal(result, pd.Series([False, True]))
        result = ser != ser[1]
        tm.assert_series_equal(result, pd.Series([True, True]))

        result = ser == ser
        tm.assert_series_equal(result, pd.Series([True, False]))
        result = ser == ser[0]
        tm.assert_series_equal(result, pd.Series([True, False]))
        result = ser == ser[1]
        tm.assert_series_equal(result, pd.Series([False, False]))


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


class TestTimestampSeriesArithmetic(object):
    def test_timestamp_sub_series(self):
        ser = pd.Series(pd.date_range('2014-03-17', periods=2, freq='D',
                                      tz='US/Eastern'))
        ts = ser[0]

        delta_series = pd.Series([np.timedelta64(0, 'D'),
                                  np.timedelta64(1, 'D')])
        tm.assert_series_equal(ser - ts, delta_series)
        tm.assert_series_equal(ts - ser, -delta_series)
