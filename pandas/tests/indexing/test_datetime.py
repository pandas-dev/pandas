import numpy as np
import pandas as pd
from pandas import date_range, Index, DataFrame, Series, Timestamp
from pandas.util import testing as tm


class TestDatetimeIndex(tm.TestCase):

    def test_indexing_with_datetime_tz(self):

        # 8260
        # support datetime64 with tz

        idx = Index(date_range('20130101', periods=3, tz='US/Eastern'),
                    name='foo')
        dr = date_range('20130110', periods=3)
        df = DataFrame({'A': idx, 'B': dr})
        df['C'] = idx
        df.iloc[1, 1] = pd.NaT
        df.iloc[1, 2] = pd.NaT

        # indexing
        result = df.iloc[1]
        expected = Series([Timestamp('2013-01-02 00:00:00-0500',
                                     tz='US/Eastern'), np.nan, np.nan],
                          index=list('ABC'), dtype='object', name=1)
        tm.assert_series_equal(result, expected)
        result = df.loc[1]
        expected = Series([Timestamp('2013-01-02 00:00:00-0500',
                                     tz='US/Eastern'), np.nan, np.nan],
                          index=list('ABC'), dtype='object', name=1)
        tm.assert_series_equal(result, expected)

        # indexing - fast_xs
        df = DataFrame({'a': date_range('2014-01-01', periods=10, tz='UTC')})
        result = df.iloc[5]
        expected = Timestamp('2014-01-06 00:00:00+0000', tz='UTC', freq='D')
        self.assertEqual(result, expected)

        result = df.loc[5]
        self.assertEqual(result, expected)

        # indexing - boolean
        result = df[df.a > df.a[3]]
        expected = df.iloc[4:]
        tm.assert_frame_equal(result, expected)

        # indexing - setting an element
        df = DataFrame(data=pd.to_datetime(
            ['2015-03-30 20:12:32', '2015-03-12 00:11:11']), columns=['time'])
        df['new_col'] = ['new', 'old']
        df.time = df.set_index('time').index.tz_localize('UTC')
        v = df[df.new_col == 'new'].set_index('time').index.tz_convert(
            'US/Pacific')

        # trying to set a single element on a part of a different timezone
        def f():
            df.loc[df.new_col == 'new', 'time'] = v

        self.assertRaises(ValueError, f)

        v = df.loc[df.new_col == 'new', 'time'] + pd.Timedelta('1s')
        df.loc[df.new_col == 'new', 'time'] = v
        tm.assert_series_equal(df.loc[df.new_col == 'new', 'time'], v)

    def test_indexing_with_datetimeindex_tz(self):

        # GH 12050
        # indexing on a series with a datetimeindex with tz
        index = pd.date_range('2015-01-01', periods=2, tz='utc')

        ser = pd.Series(range(2), index=index,
                        dtype='int64')

        # list-like indexing

        for sel in (index, list(index)):
            # getitem
            tm.assert_series_equal(ser[sel], ser)

            # setitem
            result = ser.copy()
            result[sel] = 1
            expected = pd.Series(1, index=index)
            tm.assert_series_equal(result, expected)

            # .loc getitem
            tm.assert_series_equal(ser.loc[sel], ser)

            # .loc setitem
            result = ser.copy()
            result.loc[sel] = 1
            expected = pd.Series(1, index=index)
            tm.assert_series_equal(result, expected)

        # single element indexing

        # getitem
        self.assertEqual(ser[index[1]], 1)

        # setitem
        result = ser.copy()
        result[index[1]] = 5
        expected = pd.Series([0, 5], index=index)
        tm.assert_series_equal(result, expected)

        # .loc getitem
        self.assertEqual(ser.loc[index[1]], 1)

        # .loc setitem
        result = ser.copy()
        result.loc[index[1]] = 5
        expected = pd.Series([0, 5], index=index)
        tm.assert_series_equal(result, expected)

    def test_partial_setting_with_datetimelike_dtype(self):

        # GH9478
        # a datetimeindex alignment issue with partial setting
        df = pd.DataFrame(np.arange(6.).reshape(3, 2), columns=list('AB'),
                          index=pd.date_range('1/1/2000', periods=3,
                                              freq='1H'))
        expected = df.copy()
        expected['C'] = [expected.index[0]] + [pd.NaT, pd.NaT]

        mask = df.A < 1
        df.loc[mask, 'C'] = df.loc[mask].index
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_datetime(self):

        # GH 9516
        dt1 = Timestamp('20130101 09:00:00')
        dt2 = Timestamp('20130101 10:00:00')

        for conv in [lambda x: x, lambda x: x.to_datetime64(),
                     lambda x: x.to_pydatetime(), lambda x: np.datetime64(x)]:

            df = pd.DataFrame()
            df.loc[conv(dt1), 'one'] = 100
            df.loc[conv(dt2), 'one'] = 200

            expected = DataFrame({'one': [100.0, 200.0]}, index=[dt1, dt2])
            tm.assert_frame_equal(df, expected)

    def test_series_partial_set_datetime(self):
        # GH 11497

        idx = date_range('2011-01-01', '2011-01-02', freq='D', name='idx')
        ser = Series([0.1, 0.2], index=idx, name='s')

        result = ser.loc[[Timestamp('2011-01-01'), Timestamp('2011-01-02')]]
        exp = Series([0.1, 0.2], index=idx, name='s')
        tm.assert_series_equal(result, exp, check_index_type=True)

        keys = [Timestamp('2011-01-02'), Timestamp('2011-01-02'),
                Timestamp('2011-01-01')]
        exp = Series([0.2, 0.2, 0.1], index=pd.DatetimeIndex(keys, name='idx'),
                     name='s')
        tm.assert_series_equal(ser.loc[keys], exp, check_index_type=True)

        keys = [Timestamp('2011-01-03'), Timestamp('2011-01-02'),
                Timestamp('2011-01-03')]
        exp = Series([np.nan, 0.2, np.nan],
                     index=pd.DatetimeIndex(keys, name='idx'), name='s')
        tm.assert_series_equal(ser.loc[keys], exp, check_index_type=True)

    def test_series_partial_set_period(self):
        # GH 11497

        idx = pd.period_range('2011-01-01', '2011-01-02', freq='D', name='idx')
        ser = Series([0.1, 0.2], index=idx, name='s')

        result = ser.loc[[pd.Period('2011-01-01', freq='D'),
                          pd.Period('2011-01-02', freq='D')]]
        exp = Series([0.1, 0.2], index=idx, name='s')
        tm.assert_series_equal(result, exp, check_index_type=True)

        keys = [pd.Period('2011-01-02', freq='D'),
                pd.Period('2011-01-02', freq='D'),
                pd.Period('2011-01-01', freq='D')]
        exp = Series([0.2, 0.2, 0.1], index=pd.PeriodIndex(keys, name='idx'),
                     name='s')
        tm.assert_series_equal(ser.loc[keys], exp, check_index_type=True)

        keys = [pd.Period('2011-01-03', freq='D'),
                pd.Period('2011-01-02', freq='D'),
                pd.Period('2011-01-03', freq='D')]
        exp = Series([np.nan, 0.2, np.nan],
                     index=pd.PeriodIndex(keys, name='idx'), name='s')
        result = ser.loc[keys]
        tm.assert_series_equal(result, exp)
