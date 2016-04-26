# coding=utf-8
# pylint: disable-msg=E1101,W0612

from datetime import datetime

import numpy as np
import pandas as pd

from pandas import (Index, Series, DataFrame, bdate_range,
                    date_range, period_range, timedelta_range)
from pandas.tseries.period import PeriodIndex
from pandas.tseries.index import Timestamp, DatetimeIndex
from pandas.tseries.tdi import TimedeltaIndex
import pandas.core.common as com

from pandas.util.testing import assert_series_equal
import pandas.util.testing as tm

from .common import TestData


class TestSeriesDatetimeValues(TestData, tm.TestCase):

    _multiprocess_can_split_ = True

    def test_dt_namespace_accessor(self):

        # GH 7207, 11128
        # test .dt namespace accessor

        ok_for_base = ['year', 'month', 'day', 'hour', 'minute', 'second',
                       'weekofyear', 'week', 'dayofweek', 'weekday',
                       'dayofyear', 'quarter', 'freq', 'days_in_month',
                       'daysinmonth']
        ok_for_period = ok_for_base + ['qyear', 'start_time', 'end_time']
        ok_for_period_methods = ['strftime', 'to_timestamp', 'asfreq']
        ok_for_dt = ok_for_base + ['date', 'time', 'microsecond', 'nanosecond',
                                   'is_month_start', 'is_month_end',
                                   'is_quarter_start', 'is_quarter_end',
                                   'is_year_start', 'is_year_end', 'tz',
                                   'weekday_name']
        ok_for_dt_methods = ['to_period', 'to_pydatetime', 'tz_localize',
                             'tz_convert', 'normalize', 'strftime', 'round',
                             'floor', 'ceil', 'weekday_name']
        ok_for_td = ['days', 'seconds', 'microseconds', 'nanoseconds']
        ok_for_td_methods = ['components', 'to_pytimedelta', 'total_seconds',
                             'round', 'floor', 'ceil']

        def get_expected(s, name):
            result = getattr(Index(s._values), prop)
            if isinstance(result, np.ndarray):
                if com.is_integer_dtype(result):
                    result = result.astype('int64')
            elif not com.is_list_like(result):
                return result
            return Series(result, index=s.index, name=s.name)

        def compare(s, name):
            a = getattr(s.dt, prop)
            b = get_expected(s, prop)
            if not (com.is_list_like(a) and com.is_list_like(b)):
                self.assertEqual(a, b)
            else:
                tm.assert_series_equal(a, b)

        # datetimeindex
        cases = [Series(date_range('20130101', periods=5), name='xxx'),
                 Series(date_range('20130101', periods=5, freq='s'),
                        name='xxx'),
                 Series(date_range('20130101 00:00:00', periods=5, freq='ms'),
                        name='xxx')]
        for s in cases:
            for prop in ok_for_dt:
                # we test freq below
                if prop != 'freq':
                    compare(s, prop)

            for prop in ok_for_dt_methods:
                getattr(s.dt, prop)

            result = s.dt.to_pydatetime()
            self.assertIsInstance(result, np.ndarray)
            self.assertTrue(result.dtype == object)

            result = s.dt.tz_localize('US/Eastern')
            exp_values = DatetimeIndex(s.values).tz_localize('US/Eastern')
            expected = Series(exp_values, index=s.index, name='xxx')
            tm.assert_series_equal(result, expected)

            tz_result = result.dt.tz
            self.assertEqual(str(tz_result), 'US/Eastern')
            freq_result = s.dt.freq
            self.assertEqual(freq_result, DatetimeIndex(s.values,
                                                        freq='infer').freq)

            # let's localize, then convert
            result = s.dt.tz_localize('UTC').dt.tz_convert('US/Eastern')
            exp_values = (DatetimeIndex(s.values).tz_localize('UTC')
                                                 .tz_convert('US/Eastern'))
            expected = Series(exp_values, index=s.index, name='xxx')
            tm.assert_series_equal(result, expected)

        # round
        s = Series(pd.to_datetime(['2012-01-01 13:00:00',
                                   '2012-01-01 12:01:00',
                                   '2012-01-01 08:00:00']), name='xxx')
        result = s.dt.round('D')
        expected = Series(pd.to_datetime(['2012-01-02', '2012-01-02',
                                          '2012-01-01']), name='xxx')
        tm.assert_series_equal(result, expected)

        # round with tz
        result = (s.dt.tz_localize('UTC')
                   .dt.tz_convert('US/Eastern')
                   .dt.round('D'))
        exp_values = pd.to_datetime(['2012-01-01', '2012-01-01',
                                     '2012-01-01']).tz_localize('US/Eastern')
        expected = Series(exp_values, name='xxx')
        tm.assert_series_equal(result, expected)

        # floor
        s = Series(pd.to_datetime(['2012-01-01 13:00:00',
                                   '2012-01-01 12:01:00',
                                   '2012-01-01 08:00:00']), name='xxx')
        result = s.dt.floor('D')
        expected = Series(pd.to_datetime(['2012-01-01', '2012-01-01',
                                          '2012-01-01']), name='xxx')
        tm.assert_series_equal(result, expected)

        # ceil
        s = Series(pd.to_datetime(['2012-01-01 13:00:00',
                                   '2012-01-01 12:01:00',
                                   '2012-01-01 08:00:00']), name='xxx')
        result = s.dt.ceil('D')
        expected = Series(pd.to_datetime(['2012-01-02', '2012-01-02',
                                          '2012-01-02']), name='xxx')
        tm.assert_series_equal(result, expected)

        # datetimeindex with tz
        s = Series(date_range('20130101', periods=5, tz='US/Eastern'),
                   name='xxx')
        for prop in ok_for_dt:

            # we test freq below
            if prop != 'freq':
                compare(s, prop)

        for prop in ok_for_dt_methods:
            getattr(s.dt, prop)

        result = s.dt.to_pydatetime()
        self.assertIsInstance(result, np.ndarray)
        self.assertTrue(result.dtype == object)

        result = s.dt.tz_convert('CET')
        expected = Series(s._values.tz_convert('CET'),
                          index=s.index, name='xxx')
        tm.assert_series_equal(result, expected)

        tz_result = result.dt.tz
        self.assertEqual(str(tz_result), 'CET')
        freq_result = s.dt.freq
        self.assertEqual(freq_result, DatetimeIndex(s.values,
                                                    freq='infer').freq)

        # timedeltaindex
        cases = [Series(timedelta_range('1 day', periods=5),
                        index=list('abcde'), name='xxx'),
                 Series(timedelta_range('1 day 01:23:45', periods=5,
                        freq='s'), name='xxx'),
                 Series(timedelta_range('2 days 01:23:45.012345', periods=5,
                        freq='ms'), name='xxx')]
        for s in cases:
            for prop in ok_for_td:
                # we test freq below
                if prop != 'freq':
                    compare(s, prop)

            for prop in ok_for_td_methods:
                getattr(s.dt, prop)

            result = s.dt.components
            self.assertIsInstance(result, DataFrame)
            tm.assert_index_equal(result.index, s.index)

            result = s.dt.to_pytimedelta()
            self.assertIsInstance(result, np.ndarray)
            self.assertTrue(result.dtype == object)

            result = s.dt.total_seconds()
            self.assertIsInstance(result, pd.Series)
            self.assertTrue(result.dtype == 'float64')

            freq_result = s.dt.freq
            self.assertEqual(freq_result, TimedeltaIndex(s.values,
                                                         freq='infer').freq)

        # both
        index = date_range('20130101', periods=3, freq='D')
        s = Series(date_range('20140204', periods=3, freq='s'),
                   index=index, name='xxx')
        exp = Series(np.array([2014, 2014, 2014], dtype='int64'),
                     index=index, name='xxx')
        tm.assert_series_equal(s.dt.year, exp)

        exp = Series(np.array([2, 2, 2], dtype='int64'),
                     index=index, name='xxx')
        tm.assert_series_equal(s.dt.month, exp)

        exp = Series(np.array([0, 1, 2], dtype='int64'),
                     index=index, name='xxx')
        tm.assert_series_equal(s.dt.second, exp)

        exp = pd.Series([s[0]] * 3, index=index, name='xxx')
        tm.assert_series_equal(s.dt.normalize(), exp)

        # periodindex
        cases = [Series(period_range('20130101', periods=5, freq='D'),
                        name='xxx')]
        for s in cases:
            for prop in ok_for_period:
                # we test freq below
                if prop != 'freq':
                    compare(s, prop)

            for prop in ok_for_period_methods:
                getattr(s.dt, prop)

            freq_result = s.dt.freq
            self.assertEqual(freq_result, PeriodIndex(s.values).freq)

        # test limited display api
        def get_dir(s):
            results = [r for r in s.dt.__dir__() if not r.startswith('_')]
            return list(sorted(set(results)))

        s = Series(date_range('20130101', periods=5, freq='D'), name='xxx')
        results = get_dir(s)
        tm.assert_almost_equal(
            results, list(sorted(set(ok_for_dt + ok_for_dt_methods))))

        s = Series(period_range('20130101', periods=5,
                                freq='D', name='xxx').asobject)
        results = get_dir(s)
        tm.assert_almost_equal(
            results, list(sorted(set(ok_for_period + ok_for_period_methods))))

        # 11295
        # ambiguous time error on the conversions
        s = Series(pd.date_range('2015-01-01', '2016-01-01',
                                 freq='T'), name='xxx')
        s = s.dt.tz_localize('UTC').dt.tz_convert('America/Chicago')
        results = get_dir(s)
        tm.assert_almost_equal(
            results, list(sorted(set(ok_for_dt + ok_for_dt_methods))))
        exp_values = pd.date_range('2015-01-01', '2016-01-01', freq='T',
                                   tz='UTC').tz_convert('America/Chicago')
        expected = Series(exp_values, name='xxx')
        tm.assert_series_equal(s, expected)

        # no setting allowed
        s = Series(date_range('20130101', periods=5, freq='D'), name='xxx')
        with tm.assertRaisesRegexp(ValueError, "modifications"):
            s.dt.hour = 5

        # trying to set a copy
        with pd.option_context('chained_assignment', 'raise'):

            def f():
                s.dt.hour[0] = 5

            self.assertRaises(com.SettingWithCopyError, f)

    def test_dt_accessor_no_new_attributes(self):
        # https://github.com/pydata/pandas/issues/10673
        s = Series(date_range('20130101', periods=5, freq='D'))
        with tm.assertRaisesRegexp(AttributeError,
                                   "You cannot add any new attribute"):
            s.dt.xlabel = "a"

    def test_strftime(self):
        # GH 10086
        s = Series(date_range('20130101', periods=5))
        result = s.dt.strftime('%Y/%m/%d')
        expected = Series(['2013/01/01', '2013/01/02', '2013/01/03',
                           '2013/01/04', '2013/01/05'])
        tm.assert_series_equal(result, expected)

        s = Series(date_range('2015-02-03 11:22:33.4567', periods=5))
        result = s.dt.strftime('%Y/%m/%d %H-%M-%S')
        expected = Series(['2015/02/03 11-22-33', '2015/02/04 11-22-33',
                           '2015/02/05 11-22-33', '2015/02/06 11-22-33',
                           '2015/02/07 11-22-33'])
        tm.assert_series_equal(result, expected)

        s = Series(period_range('20130101', periods=5))
        result = s.dt.strftime('%Y/%m/%d')
        expected = Series(['2013/01/01', '2013/01/02', '2013/01/03',
                           '2013/01/04', '2013/01/05'])
        tm.assert_series_equal(result, expected)

        s = Series(period_range(
            '2015-02-03 11:22:33.4567', periods=5, freq='s'))
        result = s.dt.strftime('%Y/%m/%d %H-%M-%S')
        expected = Series(['2015/02/03 11-22-33', '2015/02/03 11-22-34',
                           '2015/02/03 11-22-35', '2015/02/03 11-22-36',
                           '2015/02/03 11-22-37'])
        tm.assert_series_equal(result, expected)

        s = Series(date_range('20130101', periods=5))
        s.iloc[0] = pd.NaT
        result = s.dt.strftime('%Y/%m/%d')
        expected = Series(['NaT', '2013/01/02', '2013/01/03', '2013/01/04',
                           '2013/01/05'])
        tm.assert_series_equal(result, expected)

        datetime_index = date_range('20150301', periods=5)
        result = datetime_index.strftime("%Y/%m/%d")
        expected = np.array(
            ['2015/03/01', '2015/03/02', '2015/03/03', '2015/03/04',
             '2015/03/05'], dtype=object)
        self.assert_numpy_array_equal(result, expected)

        period_index = period_range('20150301', periods=5)
        result = period_index.strftime("%Y/%m/%d")
        expected = np.array(
            ['2015/03/01', '2015/03/02', '2015/03/03', '2015/03/04',
             '2015/03/05'], dtype=object)
        self.assert_numpy_array_equal(result, expected)

        s = Series([datetime(2013, 1, 1, 2, 32, 59), datetime(2013, 1, 2, 14,
                                                              32, 1)])
        result = s.dt.strftime('%Y-%m-%d %H:%M:%S')
        expected = Series(["2013-01-01 02:32:59", "2013-01-02 14:32:01"])
        tm.assert_series_equal(result, expected)

        s = Series(period_range('20130101', periods=4, freq='H'))
        result = s.dt.strftime('%Y/%m/%d %H:%M:%S')
        expected = Series(["2013/01/01 00:00:00", "2013/01/01 01:00:00",
                           "2013/01/01 02:00:00", "2013/01/01 03:00:00"])

        s = Series(period_range('20130101', periods=4, freq='L'))
        result = s.dt.strftime('%Y/%m/%d %H:%M:%S.%l')
        expected = Series(
            ["2013/01/01 00:00:00.000", "2013/01/01 00:00:00.001",
             "2013/01/01 00:00:00.002", "2013/01/01 00:00:00.003"])
        tm.assert_series_equal(result, expected)

    def test_valid_dt_with_missing_values(self):

        from datetime import date, time

        # GH 8689
        s = Series(date_range('20130101', periods=5, freq='D'))
        s.iloc[2] = pd.NaT

        for attr in ['microsecond', 'nanosecond', 'second', 'minute', 'hour',
                     'day']:
            expected = getattr(s.dt, attr).copy()
            expected.iloc[2] = np.nan
            result = getattr(s.dt, attr)
            tm.assert_series_equal(result, expected)

        result = s.dt.date
        expected = Series(
            [date(2013, 1, 1), date(2013, 1, 2), np.nan, date(2013, 1, 4),
             date(2013, 1, 5)], dtype='object')
        tm.assert_series_equal(result, expected)

        result = s.dt.time
        expected = Series(
            [time(0), time(0), np.nan, time(0), time(0)], dtype='object')
        tm.assert_series_equal(result, expected)

    def test_dt_accessor_api(self):
        # GH 9322
        from pandas.tseries.common import (CombinedDatetimelikeProperties,
                                           DatetimeProperties)
        self.assertIs(Series.dt, CombinedDatetimelikeProperties)

        s = Series(date_range('2000-01-01', periods=3))
        self.assertIsInstance(s.dt, DatetimeProperties)

        for s in [Series(np.arange(5)), Series(list('abcde')),
                  Series(np.random.randn(5))]:
            with tm.assertRaisesRegexp(AttributeError,
                                       "only use .dt accessor"):
                s.dt
            self.assertFalse(hasattr(s, 'dt'))

    def test_sub_of_datetime_from_TimeSeries(self):
        from pandas.tseries.timedeltas import to_timedelta
        from datetime import datetime
        a = Timestamp(datetime(1993, 0o1, 0o7, 13, 30, 00))
        b = datetime(1993, 6, 22, 13, 30)
        a = Series([a])
        result = to_timedelta(np.abs(a - b))
        self.assertEqual(result.dtype, 'timedelta64[ns]')

    def test_between(self):
        s = Series(bdate_range('1/1/2000', periods=20).asobject)
        s[::2] = np.nan

        result = s[s.between(s[3], s[17])]
        expected = s[3:18].dropna()
        assert_series_equal(result, expected)

        result = s[s.between(s[3], s[17], inclusive=False)]
        expected = s[5:16].dropna()
        assert_series_equal(result, expected)
