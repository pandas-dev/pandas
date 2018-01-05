from __future__ import division

import pytest
import numpy as np
from datetime import timedelta
from pandas import (
    Interval, IntervalIndex, Timestamp, Timedelta, DateOffset,
    interval_range, date_range, timedelta_range)
from pandas.tseries.offsets import Day
import pandas.util.testing as tm
import pandas as pd


@pytest.fixture(scope='class', params=['left', 'right', 'both', 'neither'])
def closed(request):
    return request.param


@pytest.fixture(scope='class', params=[None, 'foo'])
def name(request):
    return request.param


class TestIntervalRange(object):

    def test_construction_from_numeric(self, closed, name):
        # combinations of start/end/periods without freq
        expected = IntervalIndex.from_breaks(
            np.arange(0, 6), name=name, closed=closed)

        result = interval_range(start=0, end=5, name=name, closed=closed)
        tm.assert_index_equal(result, expected)

        result = interval_range(start=0, periods=5, name=name, closed=closed)
        tm.assert_index_equal(result, expected)

        result = interval_range(end=5, periods=5, name=name, closed=closed)
        tm.assert_index_equal(result, expected)

        # combinations of start/end/periods with freq
        expected = IntervalIndex.from_tuples([(0, 2), (2, 4), (4, 6)],
                                             name=name, closed=closed)

        result = interval_range(start=0, end=6, freq=2, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        result = interval_range(start=0, periods=3, freq=2, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        result = interval_range(end=6, periods=3, freq=2, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        # output truncates early if freq causes end to be skipped.
        expected = IntervalIndex.from_tuples([(0.0, 1.5), (1.5, 3.0)],
                                             name=name, closed=closed)
        result = interval_range(start=0, end=4, freq=1.5, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('tz', [None, 'US/Eastern'])
    def test_construction_from_timestamp(self, closed, name, tz):
        # combinations of start/end/periods without freq
        start = Timestamp('2017-01-01', tz=tz)
        end = Timestamp('2017-01-06', tz=tz)
        breaks = date_range(start=start, end=end)
        expected = IntervalIndex.from_breaks(breaks, name=name, closed=closed)

        result = interval_range(start=start, end=end, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        result = interval_range(start=start, periods=5, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        result = interval_range(end=end, periods=5, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        # combinations of start/end/periods with fixed freq
        freq = '2D'
        start = Timestamp('2017-01-01', tz=tz)
        end = Timestamp('2017-01-07', tz=tz)
        breaks = date_range(start=start, end=end, freq=freq)
        expected = IntervalIndex.from_breaks(breaks, name=name, closed=closed)

        result = interval_range(start=start, end=end, freq=freq, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        result = interval_range(start=start, periods=3, freq=freq, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        result = interval_range(end=end, periods=3, freq=freq, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        # output truncates early if freq causes end to be skipped.
        end = Timestamp('2017-01-08', tz=tz)
        result = interval_range(start=start, end=end, freq=freq, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        # combinations of start/end/periods with non-fixed freq
        freq = 'M'
        start = Timestamp('2017-01-01', tz=tz)
        end = Timestamp('2017-12-31', tz=tz)
        breaks = date_range(start=start, end=end, freq=freq)
        expected = IntervalIndex.from_breaks(breaks, name=name, closed=closed)

        result = interval_range(start=start, end=end, freq=freq, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        result = interval_range(start=start, periods=11, freq=freq, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        result = interval_range(end=end, periods=11, freq=freq, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        # output truncates early if freq causes end to be skipped.
        end = Timestamp('2018-01-15', tz=tz)
        result = interval_range(start=start, end=end, freq=freq, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

    def test_construction_from_timedelta(self, closed, name):
        # combinations of start/end/periods without freq
        start, end = Timedelta('1 day'), Timedelta('6 days')
        breaks = timedelta_range(start=start, end=end)
        expected = IntervalIndex.from_breaks(breaks, name=name, closed=closed)

        result = interval_range(start=start, end=end, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        result = interval_range(start=start, periods=5, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        result = interval_range(end=end, periods=5, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        # combinations of start/end/periods with fixed freq
        freq = '2D'
        start, end = Timedelta('1 day'), Timedelta('7 days')
        breaks = timedelta_range(start=start, end=end, freq=freq)
        expected = IntervalIndex.from_breaks(breaks, name=name, closed=closed)

        result = interval_range(start=start, end=end, freq=freq, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        result = interval_range(start=start, periods=3, freq=freq, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        result = interval_range(end=end, periods=3, freq=freq, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

        # output truncates early if freq causes end to be skipped.
        end = Timedelta('7 days 1 hour')
        result = interval_range(start=start, end=end, freq=freq, name=name,
                                closed=closed)
        tm.assert_index_equal(result, expected)

    def test_constructor_coverage(self):
        # float value for periods
        expected = pd.interval_range(start=0, periods=10)
        result = pd.interval_range(start=0, periods=10.5)
        tm.assert_index_equal(result, expected)

        # equivalent timestamp-like start/end
        start, end = Timestamp('2017-01-01'), Timestamp('2017-01-15')
        expected = pd.interval_range(start=start, end=end)

        result = pd.interval_range(start=start.to_pydatetime(),
                                   end=end.to_pydatetime())
        tm.assert_index_equal(result, expected)

        result = pd.interval_range(start=start.asm8, end=end.asm8)
        tm.assert_index_equal(result, expected)

        # equivalent freq with timestamp
        equiv_freq = ['D', Day(), Timedelta(days=1), timedelta(days=1),
                      DateOffset(days=1)]
        for freq in equiv_freq:
            result = pd.interval_range(start=start, end=end, freq=freq)
            tm.assert_index_equal(result, expected)

        # equivalent timedelta-like start/end
        start, end = Timedelta(days=1), Timedelta(days=10)
        expected = pd.interval_range(start=start, end=end)

        result = pd.interval_range(start=start.to_pytimedelta(),
                                   end=end.to_pytimedelta())
        tm.assert_index_equal(result, expected)

        result = pd.interval_range(start=start.asm8, end=end.asm8)
        tm.assert_index_equal(result, expected)

        # equivalent freq with timedelta
        equiv_freq = ['D', Day(), Timedelta(days=1), timedelta(days=1)]
        for freq in equiv_freq:
            result = pd.interval_range(start=start, end=end, freq=freq)
            tm.assert_index_equal(result, expected)

    def test_errors(self):
        # not enough params
        msg = ('Of the three parameters: start, end, and periods, '
               'exactly two must be specified')

        with tm.assert_raises_regex(ValueError, msg):
            interval_range(start=0)

        with tm.assert_raises_regex(ValueError, msg):
            interval_range(end=5)

        with tm.assert_raises_regex(ValueError, msg):
            interval_range(periods=2)

        with tm.assert_raises_regex(ValueError, msg):
            interval_range()

        # too many params
        with tm.assert_raises_regex(ValueError, msg):
            interval_range(start=0, end=5, periods=6)

        # mixed units
        msg = 'start, end, freq need to be type compatible'
        with tm.assert_raises_regex(TypeError, msg):
            interval_range(start=0, end=Timestamp('20130101'), freq=2)

        with tm.assert_raises_regex(TypeError, msg):
            interval_range(start=0, end=Timedelta('1 day'), freq=2)

        with tm.assert_raises_regex(TypeError, msg):
            interval_range(start=0, end=10, freq='D')

        with tm.assert_raises_regex(TypeError, msg):
            interval_range(start=Timestamp('20130101'), end=10, freq='D')

        with tm.assert_raises_regex(TypeError, msg):
            interval_range(start=Timestamp('20130101'),
                           end=Timedelta('1 day'), freq='D')

        with tm.assert_raises_regex(TypeError, msg):
            interval_range(start=Timestamp('20130101'),
                           end=Timestamp('20130110'), freq=2)

        with tm.assert_raises_regex(TypeError, msg):
            interval_range(start=Timedelta('1 day'), end=10, freq='D')

        with tm.assert_raises_regex(TypeError, msg):
            interval_range(start=Timedelta('1 day'),
                           end=Timestamp('20130110'), freq='D')

        with tm.assert_raises_regex(TypeError, msg):
            interval_range(start=Timedelta('1 day'),
                           end=Timedelta('10 days'), freq=2)

        # invalid periods
        msg = 'periods must be a number, got foo'
        with tm.assert_raises_regex(TypeError, msg):
            interval_range(start=0, periods='foo')

        # invalid start
        msg = 'start must be numeric or datetime-like, got foo'
        with tm.assert_raises_regex(ValueError, msg):
            interval_range(start='foo', periods=10)

        # invalid end
        msg = r'end must be numeric or datetime-like, got \(0, 1\]'
        with tm.assert_raises_regex(ValueError, msg):
            interval_range(end=Interval(0, 1), periods=10)

        # invalid freq for datetime-like
        msg = 'freq must be numeric or convertible to DateOffset, got foo'
        with tm.assert_raises_regex(ValueError, msg):
            interval_range(start=0, end=10, freq='foo')

        with tm.assert_raises_regex(ValueError, msg):
            interval_range(start=Timestamp('20130101'), periods=10, freq='foo')

        with tm.assert_raises_regex(ValueError, msg):
            interval_range(end=Timedelta('1 day'), periods=10, freq='foo')

        # mixed tz
        start = Timestamp('2017-01-01', tz='US/Eastern')
        end = Timestamp('2017-01-07', tz='US/Pacific')
        msg = 'Start and end cannot both be tz-aware with different timezones'
        with tm.assert_raises_regex(TypeError, msg):
            interval_range(start=start, end=end)
