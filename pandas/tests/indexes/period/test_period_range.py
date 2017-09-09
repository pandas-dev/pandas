import pytest
import pandas.util.testing as tm
from pandas import date_range, period_range, PeriodIndex


class TestPeriodRange(object):

    @pytest.mark.parametrize('freq', ['D', 'W', 'M', 'Q', 'A'])
    def test_construction(self, freq):
        # non-empty
        expected = date_range(start='2017-01-01', periods=5,
                              freq=freq, name='foo').to_period()
        start, end = str(expected[0]), str(expected[-1])

        result = period_range(start=start, end=end, freq=freq, name='foo')
        tm.assert_index_equal(result, expected)

        result = period_range(start=start, periods=5, freq=freq, name='foo')
        tm.assert_index_equal(result, expected)

        result = period_range(end=end, periods=5, freq=freq, name='foo')
        tm.assert_index_equal(result, expected)

        # empty
        expected = PeriodIndex([], freq=freq, name='foo')

        result = period_range(start=start, periods=0, freq=freq, name='foo')
        tm.assert_index_equal(result, expected)

        result = period_range(end=end, periods=0, freq=freq, name='foo')
        tm.assert_index_equal(result, expected)

        result = period_range(start=end, end=start, freq=freq, name='foo')
        tm.assert_index_equal(result, expected)

    def test_errors(self):
        # not enough params
        msg = ('Of the three parameters, start, end, and periods, '
               'exactly two must be specified')
        with tm.assert_raises_regex(ValueError, msg):
            period_range(start='2017Q1')

        with tm.assert_raises_regex(ValueError, msg):
            period_range(end='2017Q1')

        with tm.assert_raises_regex(ValueError, msg):
            period_range(periods=5)

        with tm.assert_raises_regex(ValueError, msg):
            period_range()

        # too many params
        with tm.assert_raises_regex(ValueError, msg):
            period_range(start='2017Q1', end='2018Q1', periods=8, freq='Q')
