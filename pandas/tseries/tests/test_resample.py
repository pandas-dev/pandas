# pylint: disable=E1101

from datetime import datetime, timedelta
from functools import partial

from pandas.compat import range, lrange, zip, product
import numpy as np

from pandas import (Series, TimeSeries, DataFrame, Panel, Index,
                    isnull, notnull, Timestamp)

from pandas.core.groupby import DataError
from pandas.tseries.index import date_range
from pandas.tseries.tdi import timedelta_range
from pandas.tseries.offsets import Minute, BDay
from pandas.tseries.period import period_range, PeriodIndex, Period
from pandas.tseries.resample import DatetimeIndex, TimeGrouper
from pandas.tseries.frequencies import MONTHS, DAYS

import pandas.tseries.offsets as offsets
import pandas as pd

import nose

from pandas.util.testing import (assert_series_equal, assert_almost_equal,
                                 assert_frame_equal)
import pandas.util.testing as tm

bday = BDay()


class TestResample(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        dti = DatetimeIndex(start=datetime(2005, 1, 1),
                            end=datetime(2005, 1, 10), freq='Min')

        self.series = Series(np.random.rand(len(dti)), dti)

    def test_custom_grouper(self):

        dti = DatetimeIndex(freq='Min', start=datetime(2005, 1, 1),
                            end=datetime(2005, 1, 10))

        s = Series(np.array([1] * len(dti)), index=dti, dtype='int64')

        b = TimeGrouper(Minute(5))
        g = s.groupby(b)

        # check all cython functions work
        funcs = ['add', 'mean', 'prod', 'ohlc', 'min', 'max', 'var']
        for f in funcs:
            g._cython_agg_general(f)

        b = TimeGrouper(Minute(5), closed='right', label='right')
        g = s.groupby(b)
        # check all cython functions work
        funcs = ['add', 'mean', 'prod', 'ohlc', 'min', 'max', 'var']
        for f in funcs:
            g._cython_agg_general(f)

        self.assertEqual(g.ngroups, 2593)
        self.assertTrue(notnull(g.mean()).all())

        # construct expected val
        arr = [1] + [5] * 2592
        idx = dti[0:-1:5]
        idx = idx.append(dti[-1:])
        expect = Series(arr, index=idx)

        # GH2763 - return in put dtype if we can
        result = g.agg(np.sum)
        assert_series_equal(result, expect)

        df = DataFrame(np.random.rand(len(dti), 10), index=dti, dtype='float64')
        r = df.groupby(b).agg(np.sum)

        self.assertEqual(len(r.columns), 10)
        self.assertEqual(len(r.index), 2593)

    def test_resample_basic(self):
        rng = date_range('1/1/2000 00:00:00', '1/1/2000 00:13:00', freq='min',
                         name='index')
        s = Series(np.random.randn(14), index=rng)
        result = s.resample('5min', how='mean', closed='right', label='right')

        exp_idx = date_range('1/1/2000', periods=4, freq='5min', name='index')
        expected = Series([s[0], s[1:6].mean(), s[6:11].mean(), s[11:].mean()],
                          index=exp_idx)
        assert_series_equal(result, expected)
        self.assertEqual(result.index.name, 'index')

        result = s.resample('5min', how='mean', closed='left', label='right')

        exp_idx = date_range('1/1/2000 00:05', periods=3, freq='5min', name='index')
        expected = Series([s[:5].mean(), s[5:10].mean(), s[10:].mean()], index=exp_idx)
        assert_series_equal(result, expected)

        s = self.series
        result = s.resample('5Min', how='last')
        grouper = TimeGrouper(Minute(5), closed='left', label='left')
        expect = s.groupby(grouper).agg(lambda x: x[-1])
        assert_series_equal(result, expect)

    def test_resample_how(self):
        rng = date_range('1/1/2000 00:00:00', '1/1/2000 00:13:00',
                            freq='min', name='index')
        s = Series(np.random.randn(14), index=rng)
        grouplist = np.ones_like(s)
        grouplist[0] = 0
        grouplist[1:6] = 1
        grouplist[6:11] = 2
        grouplist[11:] = 3
        args = ['sum', 'mean', 'std', 'sem', 'max', 'min',
                'median', 'first', 'last', 'ohlc']

        def _ohlc(group):
            if isnull(group).all():
                return np.repeat(np.nan, 4)
            return [group[0], group.max(), group.min(), group[-1]]
        inds = date_range('1/1/2000', periods=4, freq='5min', name='index')

        for arg in args:
            if arg == 'ohlc':
                func = _ohlc
            else:
                func = arg
            try:
                result = s.resample('5min', how=arg,
                                    closed='right', label='right')

                expected = s.groupby(grouplist).agg(func)
                self.assertEqual(result.index.name, 'index')
                if arg == 'ohlc':
                    expected = DataFrame(expected.values.tolist())
                    expected.columns = ['open', 'high', 'low', 'close']
                    expected.index = Index(inds, name='index')
                    assert_frame_equal(result, expected)
                else:
                    expected.index = inds
                    assert_series_equal(result, expected)
            except BaseException as exc:

                exc.args += ('how=%s' % arg,)
                raise

    def test_resample_how_callables(self):
        # GH 7929
        data = np.arange(5, dtype=np.int64)
        ind = pd.DatetimeIndex(start='2014-01-01', periods=len(data), freq='d')
        df = pd.DataFrame({"A": data, "B": data}, index=ind)

        def fn(x, a=1):
            return str(type(x))

        class fn_class:
            def __call__(self, x):
                return str(type(x))

        df_standard = df.resample("M", how=fn)
        df_lambda = df.resample("M", how=lambda x: str(type(x)))
        df_partial = df.resample("M", how=partial(fn))
        df_partial2 = df.resample("M", how=partial(fn, a=2))
        df_class = df.resample("M", how=fn_class())

        assert_frame_equal(df_standard, df_lambda)
        assert_frame_equal(df_standard, df_partial)
        assert_frame_equal(df_standard, df_partial2)
        assert_frame_equal(df_standard, df_class)

    def test_resample_with_timedeltas(self):

        expected = DataFrame({'A' : np.arange(1480)})
        expected = expected.groupby(expected.index // 30).sum()
        expected.index = pd.timedelta_range('0 days',freq='30T',periods=50)

        df = DataFrame({'A' : np.arange(1480)},index=pd.to_timedelta(np.arange(1480),unit='T'))
        result = df.resample('30T',how='sum')

        assert_frame_equal(result, expected)

    def test_resample_rounding(self):
        # GH 8371
        # odd results when rounding is needed

        data = """date,time,value
11-08-2014,00:00:01.093,1
11-08-2014,00:00:02.159,1
11-08-2014,00:00:02.667,1
11-08-2014,00:00:03.175,1
11-08-2014,00:00:07.058,1
11-08-2014,00:00:07.362,1
11-08-2014,00:00:08.324,1
11-08-2014,00:00:08.830,1
11-08-2014,00:00:08.982,1
11-08-2014,00:00:09.815,1
11-08-2014,00:00:10.540,1
11-08-2014,00:00:11.061,1
11-08-2014,00:00:11.617,1
11-08-2014,00:00:13.607,1
11-08-2014,00:00:14.535,1
11-08-2014,00:00:15.525,1
11-08-2014,00:00:17.960,1
11-08-2014,00:00:20.674,1
11-08-2014,00:00:21.191,1"""

        from pandas.compat import StringIO
        df = pd.read_csv(StringIO(data), parse_dates={'timestamp': ['date', 'time']}, index_col='timestamp')
        df.index.name = None
        result = df.resample('6s', how='sum')
        expected = DataFrame({'value' : [4,9,4,2]},index=date_range('2014-11-08',freq='6s',periods=4))
        assert_frame_equal(result,expected)

        result = df.resample('7s', how='sum')
        expected = DataFrame({'value' : [4,10,4,1]},index=date_range('2014-11-08',freq='7s',periods=4))
        assert_frame_equal(result,expected)

        result = df.resample('11s', how='sum')
        expected = DataFrame({'value' : [11,8]},index=date_range('2014-11-08',freq='11s',periods=2))
        assert_frame_equal(result,expected)

        result = df.resample('13s', how='sum')
        expected = DataFrame({'value' : [13,6]},index=date_range('2014-11-08',freq='13s',periods=2))
        assert_frame_equal(result,expected)

        result = df.resample('17s', how='sum')
        expected = DataFrame({'value' : [16,3]},index=date_range('2014-11-08',freq='17s',periods=2))
        assert_frame_equal(result,expected)

    def test_resample_basic_from_daily(self):
        # from daily
        dti = DatetimeIndex(
            start=datetime(2005, 1, 1), end=datetime(2005, 1, 10),
            freq='D', name='index')

        s = Series(np.random.rand(len(dti)), dti)

        # to weekly
        result = s.resample('w-sun', how='last')

        self.assertEqual(len(result), 3)
        self.assertTrue((result.index.dayofweek == [6, 6, 6]).all())
        self.assertEqual(result.iloc[0], s['1/2/2005'])
        self.assertEqual(result.iloc[1], s['1/9/2005'])
        self.assertEqual(result.iloc[2], s.iloc[-1])

        result = s.resample('W-MON', how='last')
        self.assertEqual(len(result), 2)
        self.assertTrue((result.index.dayofweek == [0, 0]).all())
        self.assertEqual(result.iloc[0], s['1/3/2005'])
        self.assertEqual(result.iloc[1], s['1/10/2005'])

        result = s.resample('W-TUE', how='last')
        self.assertEqual(len(result), 2)
        self.assertTrue((result.index.dayofweek == [1, 1]).all())
        self.assertEqual(result.iloc[0], s['1/4/2005'])
        self.assertEqual(result.iloc[1], s['1/10/2005'])

        result = s.resample('W-WED', how='last')
        self.assertEqual(len(result), 2)
        self.assertTrue((result.index.dayofweek == [2, 2]).all())
        self.assertEqual(result.iloc[0], s['1/5/2005'])
        self.assertEqual(result.iloc[1], s['1/10/2005'])

        result = s.resample('W-THU', how='last')
        self.assertEqual(len(result), 2)
        self.assertTrue((result.index.dayofweek == [3, 3]).all())
        self.assertEqual(result.iloc[0], s['1/6/2005'])
        self.assertEqual(result.iloc[1], s['1/10/2005'])

        result = s.resample('W-FRI', how='last')
        self.assertEqual(len(result), 2)
        self.assertTrue((result.index.dayofweek == [4, 4]).all())
        self.assertEqual(result.iloc[0], s['1/7/2005'])
        self.assertEqual(result.iloc[1], s['1/10/2005'])

        # to biz day
        result = s.resample('B', how='last')
        self.assertEqual(len(result), 7)
        self.assertTrue((result.index.dayofweek == [4, 0, 1, 2, 3, 4, 0]).all())
        self.assertEqual(result.iloc[0], s['1/2/2005'])
        self.assertEqual(result.iloc[1], s['1/3/2005'])
        self.assertEqual(result.iloc[5], s['1/9/2005'])
        self.assertEqual(result.index.name, 'index')

    def test_resample_upsampling_picked_but_not_correct(self):

        # Test for issue #3020
        dates = date_range('01-Jan-2014','05-Jan-2014', freq='D')
        series = Series(1, index=dates)

        result = series.resample('D')
        self.assertEqual(result.index[0], dates[0])

        # GH 5955
        # incorrect deciding to upsample when the axis frequency matches the resample frequency

        import datetime
        s = Series(np.arange(1.,6),index=[datetime.datetime(1975, 1, i, 12, 0) for i in range(1, 6)])
        expected = Series(np.arange(1.,6),index=date_range('19750101',periods=5,freq='D'))

        result = s.resample('D',how='count')
        assert_series_equal(result,Series(1,index=expected.index))

        result1 = s.resample('D',how='sum')
        result2 = s.resample('D',how='mean')
        result3 = s.resample('D')
        assert_series_equal(result1,expected)
        assert_series_equal(result2,expected)
        assert_series_equal(result3,expected)

    def test_resample_frame_basic(self):
        df = tm.makeTimeDataFrame()

        b = TimeGrouper('M')
        g = df.groupby(b)

        # check all cython functions work
        funcs = ['add', 'mean', 'prod', 'min', 'max', 'var']
        for f in funcs:
            g._cython_agg_general(f)

        result = df.resample('A')
        assert_series_equal(result['A'], df['A'].resample('A'))

        result = df.resample('M')
        assert_series_equal(result['A'], df['A'].resample('M'))

        df.resample('M', kind='period')
        df.resample('W-WED', kind='period')

    def test_resample_loffset(self):
        rng = date_range('1/1/2000 00:00:00', '1/1/2000 00:13:00', freq='min')
        s = Series(np.random.randn(14), index=rng)

        result = s.resample('5min', how='mean', closed='right', label='right',
                            loffset=timedelta(minutes=1))
        idx = date_range('1/1/2000', periods=4, freq='5min')
        expected = Series([s[0], s[1:6].mean(), s[6:11].mean(), s[11:].mean()],
                          index=idx + timedelta(minutes=1))
        assert_series_equal(result, expected)

        expected = s.resample(
            '5min', how='mean', closed='right', label='right',
            loffset='1min')
        assert_series_equal(result, expected)

        expected = s.resample(
            '5min', how='mean', closed='right', label='right',
            loffset=Minute(1))
        assert_series_equal(result, expected)

        self.assertEqual(result.index.freq, Minute(5))

                # from daily
        dti = DatetimeIndex(
            start=datetime(2005, 1, 1), end=datetime(2005, 1, 10),
            freq='D')
        ser = Series(np.random.rand(len(dti)), dti)

        # to weekly
        result = ser.resample('w-sun', how='last')
        expected = ser.resample('w-sun', how='last', loffset=-bday)
        self.assertEqual(result.index[0] - bday, expected.index[0])

    def test_resample_upsample(self):
        # from daily
        dti = DatetimeIndex(
            start=datetime(2005, 1, 1), end=datetime(2005, 1, 10),
            freq='D', name='index')

        s = Series(np.random.rand(len(dti)), dti)

        # to minutely, by padding
        result = s.resample('Min', fill_method='pad')
        self.assertEqual(len(result), 12961)
        self.assertEqual(result[0], s[0])
        self.assertEqual(result[-1], s[-1])

        self.assertEqual(result.index.name, 'index')

    def test_resample_extra_index_point(self):
        # GH 9756
        index = DatetimeIndex(start='20150101', end='20150331', freq='BM')
        expected = DataFrame({'A' : Series([21,41,63], index=index)})

        index = DatetimeIndex(start='20150101', end='20150331', freq='B')
        df = DataFrame({'A' : Series(range(len(index)),index=index)},dtype='int64')
        result = df.resample('BM', how='last')
        assert_frame_equal(result, expected)

    def test_upsample_with_limit(self):
        rng = date_range('1/1/2000', periods=3, freq='5t')
        ts = Series(np.random.randn(len(rng)), rng)

        result = ts.resample('t', fill_method='ffill', limit=2)
        expected = ts.reindex(result.index, method='ffill', limit=2)
        assert_series_equal(result, expected)

    def test_resample_ohlc(self):
        s = self.series

        grouper = TimeGrouper(Minute(5))
        expect = s.groupby(grouper).agg(lambda x: x[-1])
        result = s.resample('5Min', how='ohlc')

        self.assertEqual(len(result), len(expect))
        self.assertEqual(len(result.columns), 4)

        xs = result.iloc[-2]
        self.assertEqual(xs['open'], s[-6])
        self.assertEqual(xs['high'], s[-6:-1].max())
        self.assertEqual(xs['low'], s[-6:-1].min())
        self.assertEqual(xs['close'], s[-2])

        xs = result.iloc[0]
        self.assertEqual(xs['open'], s[0])
        self.assertEqual(xs['high'], s[:5].max())
        self.assertEqual(xs['low'], s[:5].min())
        self.assertEqual(xs['close'], s[4])

    def test_resample_ohlc_dataframe(self):
        df = (pd.DataFrame({'PRICE': {Timestamp('2011-01-06 10:59:05', tz=None): 24990,
                                     Timestamp('2011-01-06 12:43:33', tz=None): 25499,
                                     Timestamp('2011-01-06 12:54:09', tz=None): 25499},
                           'VOLUME': {Timestamp('2011-01-06 10:59:05', tz=None): 1500000000,
                                      Timestamp('2011-01-06 12:43:33', tz=None): 5000000000,
                                      Timestamp('2011-01-06 12:54:09', tz=None): 100000000}})
            ).reindex_axis(['VOLUME', 'PRICE'], axis=1)
        res = df.resample('H', how='ohlc')
        exp = pd.concat([df['VOLUME'].resample('H', how='ohlc'),
                          df['PRICE'].resample('H', how='ohlc')],
                        axis=1,
                        keys=['VOLUME', 'PRICE'])
        assert_frame_equal(exp, res)

        df.columns = [['a', 'b'], ['c', 'd']]
        res = df.resample('H', how='ohlc')
        exp.columns = pd.MultiIndex.from_tuples([('a', 'c', 'open'), ('a', 'c', 'high'),
            ('a', 'c', 'low'), ('a', 'c', 'close'), ('b', 'd', 'open'),
            ('b', 'd', 'high'), ('b', 'd', 'low'), ('b', 'd', 'close')])
        assert_frame_equal(exp, res)

        # dupe columns fail atm
        # df.columns = ['PRICE', 'PRICE']

    def test_resample_dup_index(self):

        # GH 4812
        # dup columns with resample raising
        df = DataFrame(np.random.randn(4,12),index=[2000,2000,2000,2000],columns=[ Period(year=2000,month=i+1,freq='M') for i in range(12) ])
        df.iloc[3,:] = np.nan
        result = df.resample('Q',axis=1)
        expected = df.groupby(lambda x: int((x.month-1)/3),axis=1).mean()
        expected.columns = [ Period(year=2000,quarter=i+1,freq='Q') for i in range(4) ]
        assert_frame_equal(result, expected)

    def test_resample_reresample(self):
        dti = DatetimeIndex(
            start=datetime(2005, 1, 1), end=datetime(2005, 1, 10),
            freq='D')
        s = Series(np.random.rand(len(dti)), dti)
        bs = s.resample('B', closed='right', label='right')
        result = bs.resample('8H')
        self.assertEqual(len(result), 22)
        tm.assertIsInstance(result.index.freq, offsets.DateOffset)
        self.assertEqual(result.index.freq, offsets.Hour(8))

    def test_resample_timestamp_to_period(self):
        ts = _simple_ts('1/1/1990', '1/1/2000')

        result = ts.resample('A-DEC', kind='period')
        expected = ts.resample('A-DEC')
        expected.index = period_range('1990', '2000', freq='a-dec')
        assert_series_equal(result, expected)

        result = ts.resample('A-JUN', kind='period')
        expected = ts.resample('A-JUN')
        expected.index = period_range('1990', '2000', freq='a-jun')
        assert_series_equal(result, expected)

        result = ts.resample('M', kind='period')
        expected = ts.resample('M')
        expected.index = period_range('1990-01', '2000-01', freq='M')
        assert_series_equal(result, expected)

        result = ts.resample('M', kind='period')
        expected = ts.resample('M')
        expected.index = period_range('1990-01', '2000-01', freq='M')
        assert_series_equal(result, expected)

    def test_ohlc_5min(self):
        def _ohlc(group):
            if isnull(group).all():
                return np.repeat(np.nan, 4)
            return [group[0], group.max(), group.min(), group[-1]]

        rng = date_range('1/1/2000 00:00:00', '1/1/2000 5:59:50',
                         freq='10s')
        ts = Series(np.random.randn(len(rng)), index=rng)

        resampled = ts.resample('5min', how='ohlc', closed='right',
                                label='right')

        self.assertTrue((resampled.ix['1/1/2000 00:00'] == ts[0]).all())

        exp = _ohlc(ts[1:31])
        self.assertTrue((resampled.ix['1/1/2000 00:05'] == exp).all())

        exp = _ohlc(ts['1/1/2000 5:55:01':])
        self.assertTrue((resampled.ix['1/1/2000 6:00:00'] == exp).all())

    def test_downsample_non_unique(self):
        rng = date_range('1/1/2000', '2/29/2000')
        rng2 = rng.repeat(5).values
        ts = Series(np.random.randn(len(rng2)), index=rng2)

        result = ts.resample('M', how='mean')

        expected = ts.groupby(lambda x: x.month).mean()
        self.assertEqual(len(result), 2)
        assert_almost_equal(result[0], expected[1])
        assert_almost_equal(result[1], expected[2])

    def test_asfreq_non_unique(self):
        # GH #1077
        rng = date_range('1/1/2000', '2/29/2000')
        rng2 = rng.repeat(2).values
        ts = Series(np.random.randn(len(rng2)), index=rng2)

        self.assertRaises(Exception, ts.asfreq, 'B')

    def test_resample_axis1(self):
        rng = date_range('1/1/2000', '2/29/2000')
        df = DataFrame(np.random.randn(3, len(rng)), columns=rng,
                       index=['a', 'b', 'c'])

        result = df.resample('M', axis=1)
        expected = df.T.resample('M').T
        tm.assert_frame_equal(result, expected)

    def test_resample_panel(self):
        rng = date_range('1/1/2000', '6/30/2000')
        n = len(rng)

        panel = Panel(np.random.randn(3, n, 5),
                      items=['one', 'two', 'three'],
                      major_axis=rng,
                      minor_axis=['a', 'b', 'c', 'd', 'e'])

        result = panel.resample('M', axis=1)

        def p_apply(panel, f):
            result = {}
            for item in panel.items:
                result[item] = f(panel[item])
            return Panel(result, items=panel.items)

        expected = p_apply(panel, lambda x: x.resample('M'))
        tm.assert_panel_equal(result, expected)

        panel2 = panel.swapaxes(1, 2)
        result = panel2.resample('M', axis=2)
        expected = p_apply(panel2, lambda x: x.resample('M', axis=1))
        tm.assert_panel_equal(result, expected)

    def test_resample_panel_numpy(self):
        rng = date_range('1/1/2000', '6/30/2000')
        n = len(rng)

        panel = Panel(np.random.randn(3, n, 5),
                      items=['one', 'two', 'three'],
                      major_axis=rng,
                      minor_axis=['a', 'b', 'c', 'd', 'e'])

        result = panel.resample('M', how=lambda x: x.mean(1), axis=1)
        expected = panel.resample('M', how='mean', axis=1)
        tm.assert_panel_equal(result, expected)

        panel = panel.swapaxes(1, 2)
        result = panel.resample('M', how=lambda x: x.mean(2), axis=2)
        expected = panel.resample('M', how='mean', axis=2)
        tm.assert_panel_equal(result, expected)

    def test_resample_anchored_ticks(self):
        # If a fixed delta (5 minute, 4 hour) evenly divides a day, we should
        # "anchor" the origin at midnight so we get regular intervals rather
        # than starting from the first timestamp which might start in the middle
        # of a desired interval

        rng = date_range('1/1/2000 04:00:00', periods=86400, freq='s')
        ts = Series(np.random.randn(len(rng)), index=rng)
        ts[:2] = np.nan  # so results are the same

        freqs = ['t', '5t', '15t', '30t', '4h', '12h']
        for freq in freqs:
            result = ts[2:].resample(freq, closed='left', label='left')
            expected = ts.resample(freq, closed='left', label='left')
            assert_series_equal(result, expected)

    def test_resample_single_group(self):
        mysum = lambda x: x.sum()

        rng = date_range('2000-1-1', '2000-2-10', freq='D')
        ts = Series(np.random.randn(len(rng)), index=rng)
        assert_series_equal(ts.resample('M', how='sum'),
                            ts.resample('M', how=mysum))

        rng = date_range('2000-1-1', '2000-1-10', freq='D')
        ts = Series(np.random.randn(len(rng)), index=rng)
        assert_series_equal(ts.resample('M', how='sum'),
                            ts.resample('M', how=mysum))

        # GH 3849
        s = Series([30.1, 31.6], index=[Timestamp('20070915 15:30:00'),
                                        Timestamp('20070915 15:40:00')])
        expected = Series([0.75], index=[Timestamp('20070915')])
        result = s.resample('D', how=lambda x: np.std(x))
        assert_series_equal(result, expected)

    def test_resample_base(self):
        rng = date_range('1/1/2000 00:00:00', '1/1/2000 02:00', freq='s')
        ts = Series(np.random.randn(len(rng)), index=rng)

        resampled = ts.resample('5min', base=2)
        exp_rng = date_range('12/31/1999 23:57:00', '1/1/2000 01:57',
                             freq='5min')
        self.assertTrue(resampled.index.equals(exp_rng))

    def test_resample_base_with_timedeltaindex(self):

        # GH 10530
        rng = timedelta_range(start = '0s', periods = 25, freq = 's')
        ts = Series(np.random.randn(len(rng)), index = rng)

        with_base = ts.resample('2s', base = 5)
        without_base = ts.resample('2s')

        exp_without_base = timedelta_range(start = '0s', end = '25s', freq = '2s')
        exp_with_base = timedelta_range(start = '5s', end = '29s', freq = '2s')

        self.assertTrue(without_base.index.equals(exp_without_base))
        self.assertTrue(with_base.index.equals(exp_with_base))

    def test_resample_daily_anchored(self):
        rng = date_range('1/1/2000 0:00:00', periods=10000, freq='T')
        ts = Series(np.random.randn(len(rng)), index=rng)
        ts[:2] = np.nan  # so results are the same

        result = ts[2:].resample('D', closed='left', label='left')
        expected = ts.resample('D', closed='left', label='left')
        assert_series_equal(result, expected)

    def test_resample_to_period_monthly_buglet(self):
        # GH #1259

        rng = date_range('1/1/2000', '12/31/2000')
        ts = Series(np.random.randn(len(rng)), index=rng)

        result = ts.resample('M', kind='period')
        exp_index = period_range('Jan-2000', 'Dec-2000', freq='M')
        self.assertTrue(result.index.equals(exp_index))

    def test_resample_empty(self):
        ts = _simple_ts('1/1/2000', '2/1/2000')[:0]

        result = ts.resample('A')
        self.assertEqual(len(result), 0)
        self.assertEqual(result.index.freqstr, 'A-DEC')

        result = ts.resample('A', kind='period')
        self.assertEqual(len(result), 0)
        self.assertEqual(result.index.freqstr, 'A-DEC')

        xp = DataFrame()
        rs = xp.resample('A')
        assert_frame_equal(xp, rs)

        # Empty series were sometimes causing a segfault (for the functions
        # with Cython bounds-checking disabled) or an IndexError.  We just run
        # them to ensure they no longer do.  (GH #10228)
        for index in tm.all_timeseries_index_generator(0):
            for dtype in (np.float, np.int, np.object, 'datetime64[ns]'):
                for how in ('count', 'mean', 'min', 'ohlc', 'last', 'prod'):
                    empty_series = pd.Series([], index, dtype)
                    try:
                        empty_series.resample('d', how)
                    except DataError:
                        # Ignore these since some combinations are invalid
                        # (ex: doing mean with dtype of np.object)
                        pass

    def test_weekly_resample_buglet(self):
        # #1327
        rng = date_range('1/1/2000', freq='B', periods=20)
        ts = Series(np.random.randn(len(rng)), index=rng)

        resampled = ts.resample('W')
        expected = ts.resample('W-SUN')
        assert_series_equal(resampled, expected)

    def test_monthly_resample_error(self):
        # #1451
        dates = date_range('4/16/2012 20:00', periods=5000, freq='h')
        ts = Series(np.random.randn(len(dates)), index=dates)
        # it works!
        result = ts.resample('M')

    def test_resample_anchored_intraday(self):
        # #1471, #1458

        rng = date_range('1/1/2012', '4/1/2012', freq='100min')
        df = DataFrame(rng.month, index=rng)

        result = df.resample('M')
        expected = df.resample('M', kind='period').to_timestamp(how='end')
        tm.assert_frame_equal(result, expected)

        result = df.resample('M', closed='left')
        exp = df.tshift(1, freq='D').resample('M', kind='period')
        exp = exp.to_timestamp(how='end')

        tm.assert_frame_equal(result, exp)

        rng = date_range('1/1/2012', '4/1/2012', freq='100min')
        df = DataFrame(rng.month, index=rng)

        result = df.resample('Q')
        expected = df.resample('Q', kind='period').to_timestamp(how='end')
        tm.assert_frame_equal(result, expected)

        result = df.resample('Q', closed='left')
        expected = df.tshift(1, freq='D').resample('Q', kind='period',
                                                   closed='left')
        expected = expected.to_timestamp(how='end')
        tm.assert_frame_equal(result, expected)

        ts = _simple_ts('2012-04-29 23:00', '2012-04-30 5:00', freq='h')
        resampled = ts.resample('M')
        self.assertEqual(len(resampled), 1)

    def test_resample_anchored_monthstart(self):
        ts = _simple_ts('1/1/2000', '12/31/2002')

        freqs = ['MS', 'BMS', 'QS-MAR', 'AS-DEC', 'AS-JUN']

        for freq in freqs:
            result = ts.resample(freq, how='mean')

    def test_resample_anchored_multiday(self):
        # When resampling a range spanning multiple days, ensure that the
        # start date gets used to determine the offset.  Fixes issue where
        # a one day period is not a multiple of the frequency.
        #
        # See: https://github.com/pydata/pandas/issues/8683

        s = pd.Series(np.random.randn(5),
                      index=pd.date_range('2014-10-14 23:06:23.206',
                                          periods=3, freq='400L')
                      | pd.date_range('2014-10-15 23:00:00',
                                      periods=2, freq='2200L'))

        # Ensure left closing works
        result = s.resample('2200L', 'mean')
        self.assertEqual(result.index[-1],
                         pd.Timestamp('2014-10-15 23:00:02.000'))

        # Ensure right closing works
        result = s.resample('2200L', 'mean', label='right')
        self.assertEqual(result.index[-1],
                         pd.Timestamp('2014-10-15 23:00:04.200'))

    def test_corner_cases(self):
        # miscellaneous test coverage

        rng = date_range('1/1/2000', periods=12, freq='t')
        ts = Series(np.random.randn(len(rng)), index=rng)

        result = ts.resample('5t', closed='right', label='left')
        ex_index = date_range('1999-12-31 23:55', periods=4, freq='5t')
        self.assertTrue(result.index.equals(ex_index))

        len0pts = _simple_pts('2007-01', '2010-05', freq='M')[:0]
        # it works
        result = len0pts.resample('A-DEC')
        self.assertEqual(len(result), 0)

        # resample to periods
        ts = _simple_ts('2000-04-28', '2000-04-30 11:00', freq='h')
        result = ts.resample('M', kind='period')
        self.assertEqual(len(result), 1)
        self.assertEqual(result.index[0], Period('2000-04', freq='M'))

    def test_anchored_lowercase_buglet(self):
        dates = date_range('4/16/2012 20:00', periods=50000, freq='s')
        ts = Series(np.random.randn(len(dates)), index=dates)
        # it works!
        ts.resample('d')

    def test_upsample_apply_functions(self):
        # #1596
        rng = pd.date_range('2012-06-12', periods=4, freq='h')

        ts = Series(np.random.randn(len(rng)), index=rng)

        result = ts.resample('20min', how=['mean', 'sum'])
        tm.assertIsInstance(result, DataFrame)

    def test_resample_not_monotonic(self):
        rng = pd.date_range('2012-06-12', periods=200, freq='h')
        ts = Series(np.random.randn(len(rng)), index=rng)

        ts = ts.take(np.random.permutation(len(ts)))

        result = ts.resample('D', how='sum')
        exp = ts.sort_index().resample('D', how='sum')
        assert_series_equal(result, exp)

    def test_resample_median_bug_1688(self):

        for dtype in ['int64','int32','float64','float32']:
            df = DataFrame([1, 2], index=[datetime(2012, 1, 1, 0, 0, 0),
                                          datetime(2012, 1, 1, 0, 5, 0)],
                           dtype = dtype)

            result = df.resample("T", how=lambda x: x.mean())
            exp = df.asfreq('T')
            tm.assert_frame_equal(result, exp)

            result = df.resample("T", how="median")
            exp = df.asfreq('T')
            tm.assert_frame_equal(result, exp)

    def test_how_lambda_functions(self):

        ts = _simple_ts('1/1/2000', '4/1/2000')

        result = ts.resample('M', how=lambda x: x.mean())
        exp = ts.resample('M', how='mean')
        tm.assert_series_equal(result, exp)

        self.assertRaises(Exception, ts.resample, 'M',
                          how=[lambda x: x.mean(), lambda x: x.std(ddof=1)])

        result = ts.resample('M', how={'foo': lambda x: x.mean(),
                                       'bar': lambda x: x.std(ddof=1)})
        foo_exp = ts.resample('M', how='mean')
        foo_exp.name = 'foo'
        bar_exp = ts.resample('M', how='std')
        bar_exp.name = 'bar'

        tm.assert_series_equal(result['foo'], foo_exp)
        tm.assert_series_equal(result['bar'], bar_exp)

    def test_resample_unequal_times(self):
        # #1772
        start = datetime(1999, 3, 1, 5)
        # end hour is less than start
        end = datetime(2012, 7, 31, 4)
        bad_ind = date_range(start, end, freq="30min")
        df = DataFrame({'close': 1}, index=bad_ind)

        # it works!
        df.resample('AS', 'sum')

    def test_resample_consistency(self):

        # GH 6418
        # resample with bfill / limit / reindex consistency

        i30 = index=pd.date_range('2002-02-02', periods=4, freq='30T')
        s=pd.Series(np.arange(4.), index=i30)
        s[2] = np.NaN

        # Upsample by factor 3 with reindex() and resample() methods:
        i10 = pd.date_range(i30[0], i30[-1], freq='10T')

        s10 = s.reindex(index=i10, method='bfill')
        s10_2 = s.reindex(index=i10, method='bfill', limit=2)
        rl = s.reindex_like(s10, method='bfill', limit=2)
        r10_2 = s.resample('10Min', fill_method='bfill', limit=2)
        r10 = s.resample('10Min', fill_method='bfill')

        # s10_2, r10, r10_2, rl should all be equal
        assert_series_equal(s10_2, r10)
        assert_series_equal(s10_2, r10_2)
        assert_series_equal(s10_2, rl)

    def test_resample_timegrouper(self):
        # GH 7227
        dates1 = [datetime(2014, 10, 1), datetime(2014, 9, 3),
                  datetime(2014, 11, 5), datetime(2014, 9, 5),
                  datetime(2014, 10, 8), datetime(2014, 7, 15)]

        dates2 = dates1[:2] + [pd.NaT] + dates1[2:4] + [pd.NaT] + dates1[4:]
        dates3 = [pd.NaT] + dates1 + [pd.NaT]

        for dates in [dates1, dates2, dates3]:
            df = DataFrame(dict(A=dates, B=np.arange(len(dates))))
            result = df.set_index('A').resample('M', how='count')
            exp_idx = pd.DatetimeIndex(['2014-07-31', '2014-08-31', '2014-09-30',
                                        '2014-10-31', '2014-11-30'], freq='M', name='A')
            expected = DataFrame({'B': [1, 0, 2, 2, 1]}, index=exp_idx)
            assert_frame_equal(result, expected)

            result = df.groupby(pd.Grouper(freq='M', key='A')).count()
            assert_frame_equal(result, expected)

            df = DataFrame(dict(A=dates, B=np.arange(len(dates)), C=np.arange(len(dates))))
            result = df.set_index('A').resample('M', how='count')
            expected = DataFrame({'B': [1, 0, 2, 2, 1], 'C': [1, 0, 2, 2, 1]},
                                 index=exp_idx, columns=['B', 'C'])
            assert_frame_equal(result, expected)

            result = df.groupby(pd.Grouper(freq='M', key='A')).count()
            assert_frame_equal(result, expected)

    def test_resample_group_info(self):  # GH10914
        for n, k in product((10000, 100000), (10, 100, 1000)):
            dr = date_range(start='2015-08-27', periods=n // 10, freq='T')
            ts = Series(np.random.randint(0, n // k, n).astype('int64'),
                        index=np.random.choice(dr, n))

            left = ts.resample('30T', how='nunique')
            ix = date_range(start=ts.index.min(),
                            end=ts.index.max(),
                            freq='30T')

            vals = ts.values
            bins = np.searchsorted(ix.values, ts.index, side='right')

            sorter = np.lexsort((vals, bins))
            vals, bins = vals[sorter], bins[sorter]

            mask = np.r_[True, vals[1:] != vals[:-1]]
            mask |= np.r_[True, bins[1:] != bins[:-1]]

            arr = np.bincount(bins[mask] - 1, minlength=len(ix)).astype('int64',copy=False)
            right = Series(arr, index=ix)

            assert_series_equal(left, right)

    def test_resample_size(self):
        n = 10000
        dr = date_range('2015-09-19', periods=n, freq='T')
        ts = Series(np.random.randn(n), index=np.random.choice(dr, n))

        left = ts.resample('7T', how='size')
        ix = date_range(start=left.index.min(), end=ts.index.max(), freq='7T')

        bins = np.searchsorted(ix.values, ts.index.values, side='right')
        val = np.bincount(bins, minlength=len(ix) + 1)[1:].astype('int64',copy=False)

        right = Series(val, index=ix)
        assert_series_equal(left, right)

    def test_resmaple_dst_anchor(self):
        # 5172
        dti = DatetimeIndex([datetime(2012, 11, 4, 23)], tz='US/Eastern')
        df = DataFrame([5], index=dti)
        assert_frame_equal(df.resample(rule='D', how='sum'),
                           DataFrame([5], index=df.index.normalize()))
        df.resample(rule='MS', how='sum')
        assert_frame_equal(df.resample(rule='MS', how='sum'),
                           DataFrame([5], index=DatetimeIndex([datetime(2012, 11, 1)],
                                                              tz='US/Eastern')))

        dti = date_range('2013-09-30', '2013-11-02', freq='30Min', tz='Europe/Paris')
        values = range(dti.size)
        df = DataFrame({"a": values, "b": values, "c": values}, index=dti, dtype='int64')
        how = {"a": "min", "b": "max", "c": "count"}

        assert_frame_equal(df.resample("W-MON", how=how)[["a", "b", "c"]],
                           DataFrame({"a": [0, 48, 384, 720, 1056, 1394],
                                      "b": [47, 383, 719, 1055, 1393, 1586],
                                      "c": [48, 336, 336, 336, 338, 193]},
                                     index=date_range('9/30/2013', '11/4/2013',
                                                      freq='W-MON', tz='Europe/Paris')),
                           'W-MON Frequency')

        assert_frame_equal(df.resample("2W-MON", how=how)[["a", "b", "c"]],
                           DataFrame({"a": [0, 48, 720, 1394],
                                      "b": [47, 719, 1393, 1586],
                                      "c": [48, 672, 674, 193]},
                                     index=date_range('9/30/2013', '11/11/2013',
                                                      freq='2W-MON', tz='Europe/Paris')),
                           '2W-MON Frequency')

        assert_frame_equal(df.resample("MS", how=how)[["a", "b", "c"]],
                           DataFrame({"a": [0, 48, 1538],
                                      "b": [47, 1537, 1586],
                                      "c": [48, 1490, 49]},
                                     index=date_range('9/1/2013', '11/1/2013',
                                                      freq='MS', tz='Europe/Paris')),
                           'MS Frequency')

        assert_frame_equal(df.resample("2MS", how=how)[["a", "b", "c"]],
                           DataFrame({"a": [0, 1538],
                                      "b": [1537, 1586],
                                      "c": [1538, 49]},
                                     index=date_range('9/1/2013', '11/1/2013',
                                                      freq='2MS', tz='Europe/Paris')),
                           '2MS Frequency')

        df_daily = df['10/26/2013':'10/29/2013']
        assert_frame_equal(df_daily.resample("D", how={"a": "min", "b": "max", "c": "count"})[["a", "b", "c"]],
                           DataFrame({"a": [1248, 1296, 1346, 1394],
                                      "b": [1295, 1345, 1393, 1441],
                                      "c": [48, 50, 48, 48]},
                                     index=date_range('10/26/2013', '10/29/2013',
                                                      freq='D', tz='Europe/Paris')),
                           'D Frequency')

def _simple_ts(start, end, freq='D'):
    rng = date_range(start, end, freq=freq)
    return Series(np.random.randn(len(rng)), index=rng)


def _simple_pts(start, end, freq='D'):
    rng = period_range(start, end, freq=freq)
    return Series(np.random.randn(len(rng)), index=rng)


class TestResamplePeriodIndex(tm.TestCase):

    _multiprocess_can_split_ = True

    def test_annual_upsample_D_s_f(self):
        self._check_annual_upsample_cases('D', 'start', 'ffill')

    def test_annual_upsample_D_e_f(self):
        self._check_annual_upsample_cases('D', 'end', 'ffill')

    def test_annual_upsample_D_s_b(self):
        self._check_annual_upsample_cases('D', 'start', 'bfill')

    def test_annual_upsample_D_e_b(self):
        self._check_annual_upsample_cases('D', 'end', 'bfill')

    def test_annual_upsample_B_s_f(self):
        self._check_annual_upsample_cases('B', 'start', 'ffill')

    def test_annual_upsample_B_e_f(self):
        self._check_annual_upsample_cases('B', 'end', 'ffill')

    def test_annual_upsample_B_s_b(self):
        self._check_annual_upsample_cases('B', 'start', 'bfill')

    def test_annual_upsample_B_e_b(self):
        self._check_annual_upsample_cases('B', 'end', 'bfill')

    def test_annual_upsample_M_s_f(self):
        self._check_annual_upsample_cases('M', 'start', 'ffill')

    def test_annual_upsample_M_e_f(self):
        self._check_annual_upsample_cases('M', 'end', 'ffill')

    def test_annual_upsample_M_s_b(self):
        self._check_annual_upsample_cases('M', 'start', 'bfill')

    def test_annual_upsample_M_e_b(self):
        self._check_annual_upsample_cases('M', 'end', 'bfill')

    def _check_annual_upsample_cases(self, targ, conv, meth, end='12/31/1991'):
        for month in MONTHS:
            ts = _simple_pts('1/1/1990', end, freq='A-%s' % month)

            result = ts.resample(targ, fill_method=meth,
                                 convention=conv)
            expected = result.to_timestamp(targ, how=conv)
            expected = expected.asfreq(targ, meth).to_period()
            assert_series_equal(result, expected)

    def test_basic_downsample(self):
        ts = _simple_pts('1/1/1990', '6/30/1995', freq='M')
        result = ts.resample('a-dec')

        expected = ts.groupby(ts.index.year).mean()
        expected.index = period_range('1/1/1990', '6/30/1995',
                                      freq='a-dec')
        assert_series_equal(result, expected)

        # this is ok
        assert_series_equal(ts.resample('a-dec'), result)
        assert_series_equal(ts.resample('a'), result)

    def test_not_subperiod(self):
        # These are incompatible period rules for resampling
        ts = _simple_pts('1/1/1990', '6/30/1995', freq='w-wed')
        self.assertRaises(ValueError, ts.resample, 'a-dec')
        self.assertRaises(ValueError, ts.resample, 'q-mar')
        self.assertRaises(ValueError, ts.resample, 'M')
        self.assertRaises(ValueError, ts.resample, 'w-thu')

    def test_basic_upsample(self):
        ts = _simple_pts('1/1/1990', '6/30/1995', freq='M')
        result = ts.resample('a-dec')

        resampled = result.resample('D', fill_method='ffill', convention='end')

        expected = result.to_timestamp('D', how='end')
        expected = expected.asfreq('D', 'ffill').to_period()

        assert_series_equal(resampled, expected)

    def test_upsample_with_limit(self):
        rng = period_range('1/1/2000', periods=5, freq='A')
        ts = Series(np.random.randn(len(rng)), rng)

        result = ts.resample('M', fill_method='ffill', limit=2,
                             convention='end')
        expected = ts.asfreq('M').reindex(result.index, method='ffill',
                                          limit=2)
        assert_series_equal(result, expected)

    def test_annual_upsample(self):
        ts = _simple_pts('1/1/1990', '12/31/1995', freq='A-DEC')
        df = DataFrame({'a': ts})
        rdf = df.resample('D', fill_method='ffill')
        exp = df['a'].resample('D', fill_method='ffill')
        assert_series_equal(rdf['a'], exp)

        rng = period_range('2000', '2003', freq='A-DEC')
        ts = Series([1, 2, 3, 4], index=rng)

        result = ts.resample('M', fill_method='ffill')
        ex_index = period_range('2000-01', '2003-12', freq='M')

        expected = ts.asfreq('M', how='start').reindex(ex_index,
                                                       method='ffill')
        assert_series_equal(result, expected)

    def test_quarterly_upsample(self):
        targets = ['D', 'B', 'M']

        for month in MONTHS:
            ts = _simple_pts('1/1/1990', '12/31/1995', freq='Q-%s' % month)

            for targ, conv in product(targets, ['start', 'end']):
                result = ts.resample(targ, fill_method='ffill',
                                     convention=conv)
                expected = result.to_timestamp(targ, how=conv)
                expected = expected.asfreq(targ, 'ffill').to_period()
                assert_series_equal(result, expected)

    def test_monthly_upsample(self):
        targets = ['D', 'B']

        ts = _simple_pts('1/1/1990', '12/31/1995', freq='M')

        for targ, conv in product(targets, ['start', 'end']):
            result = ts.resample(targ, fill_method='ffill',
                                 convention=conv)
            expected = result.to_timestamp(targ, how=conv)
            expected = expected.asfreq(targ, 'ffill').to_period()
            assert_series_equal(result, expected)

    def test_fill_method_and_how_upsample(self):
        # GH2073
        s = Series(np.arange(9,dtype='int64'),
                   index=date_range('2010-01-01', periods=9, freq='Q'))
        last = s.resample('M', fill_method='ffill')
        both = s.resample('M', how='last', fill_method='ffill').astype('int64')
        assert_series_equal(last, both)

    def test_weekly_upsample(self):
        targets = ['D', 'B']

        for day in DAYS:
            ts = _simple_pts('1/1/1990', '12/31/1995', freq='W-%s' % day)

            for targ, conv in product(targets, ['start', 'end']):
                result = ts.resample(targ, fill_method='ffill',
                                     convention=conv)
                expected = result.to_timestamp(targ, how=conv)
                expected = expected.asfreq(targ, 'ffill').to_period()
                assert_series_equal(result, expected)

    def test_resample_to_timestamps(self):
        ts = _simple_pts('1/1/1990', '12/31/1995', freq='M')

        result = ts.resample('A-DEC', kind='timestamp')
        expected = ts.to_timestamp(how='end').resample('A-DEC')
        assert_series_equal(result, expected)

    def test_resample_to_quarterly(self):
        for month in MONTHS:
            ts = _simple_pts('1990', '1992', freq='A-%s' % month)
            quar_ts = ts.resample('Q-%s' % month, fill_method='ffill')

            stamps = ts.to_timestamp('D', how='start')
            qdates = period_range(ts.index[0].asfreq('D', 'start'),
                                  ts.index[-1].asfreq('D', 'end'),
                                  freq='Q-%s' % month)

            expected = stamps.reindex(qdates.to_timestamp('D', 's'),
                                      method='ffill')
            expected.index = qdates

            assert_series_equal(quar_ts, expected)

        # conforms, but different month
        ts = _simple_pts('1990', '1992', freq='A-JUN')

        for how in ['start', 'end']:
            result = ts.resample('Q-MAR', convention=how, fill_method='ffill')
            expected = ts.asfreq('Q-MAR', how=how)
            expected = expected.reindex(result.index, method='ffill')

            # .to_timestamp('D')
            # expected = expected.resample('Q-MAR', fill_method='ffill')

            assert_series_equal(result, expected)

    def test_resample_fill_missing(self):
        rng = PeriodIndex([2000, 2005, 2007, 2009], freq='A')

        s = Series(np.random.randn(4), index=rng)

        stamps = s.to_timestamp()

        filled = s.resample('A')
        expected = stamps.resample('A').to_period('A')
        assert_series_equal(filled, expected)

        filled = s.resample('A', fill_method='ffill')
        expected = stamps.resample('A', fill_method='ffill').to_period('A')
        assert_series_equal(filled, expected)

    def test_cant_fill_missing_dups(self):
        rng = PeriodIndex([2000, 2005, 2005, 2007, 2007], freq='A')
        s = Series(np.random.randn(5), index=rng)
        self.assertRaises(Exception, s.resample, 'A')

    def test_resample_5minute(self):
        rng = period_range('1/1/2000', '1/5/2000', freq='T')
        ts = Series(np.random.randn(len(rng)), index=rng)

        result = ts.resample('5min')
        expected = ts.to_timestamp().resample('5min')
        assert_series_equal(result, expected)

    def test_upsample_daily_business_daily(self):
        ts = _simple_pts('1/1/2000', '2/1/2000', freq='B')

        result = ts.resample('D')
        expected = ts.asfreq('D').reindex(period_range('1/3/2000', '2/1/2000'))
        assert_series_equal(result, expected)

        ts = _simple_pts('1/1/2000', '2/1/2000')
        result = ts.resample('H', convention='s')
        exp_rng = period_range('1/1/2000', '2/1/2000 23:00', freq='H')
        expected = ts.asfreq('H', how='s').reindex(exp_rng)
        assert_series_equal(result, expected)

    def test_resample_empty(self):
        ts = _simple_pts('1/1/2000', '2/1/2000')[:0]

        result = ts.resample('A')
        self.assertEqual(len(result), 0)

    def test_resample_irregular_sparse(self):
        dr = date_range(start='1/1/2012', freq='5min', periods=1000)
        s = Series(np.array(100), index=dr)
        # subset the data.
        subset = s[:'2012-01-04 06:55']

        result = subset.resample('10min', how=len)
        expected = s.resample('10min', how=len).ix[result.index]
        assert_series_equal(result, expected)

    def test_resample_weekly_all_na(self):
        rng = date_range('1/1/2000', periods=10, freq='W-WED')
        ts = Series(np.random.randn(len(rng)), index=rng)

        result = ts.resample('W-THU')

        self.assertTrue(result.isnull().all())

        result = ts.resample('W-THU', fill_method='ffill')[:-1]
        expected = ts.asfreq('W-THU', method='ffill')
        assert_series_equal(result, expected)

    def test_resample_tz_localized(self):
        dr = date_range(start='2012-4-13', end='2012-5-1')
        ts = Series(lrange(len(dr)), dr)

        ts_utc = ts.tz_localize('UTC')
        ts_local = ts_utc.tz_convert('America/Los_Angeles')

        result = ts_local.resample('W')

        ts_local_naive = ts_local.copy()
        ts_local_naive.index = [x.replace(tzinfo=None)
                                for x in ts_local_naive.index.to_pydatetime()]

        exp = ts_local_naive.resample('W').tz_localize('America/Los_Angeles')

        assert_series_equal(result, exp)

        # it works
        result = ts_local.resample('D')

        # #2245
        idx = date_range('2001-09-20 15:59', '2001-09-20 16:00', freq='T',
                         tz='Australia/Sydney')
        s = Series([1, 2], index=idx)

        result = s.resample('D', closed='right', label='right')
        ex_index = date_range('2001-09-21', periods=1, freq='D',
                              tz='Australia/Sydney')
        expected = Series([1.5], index=ex_index)

        assert_series_equal(result, expected)

        # for good measure
        result = s.resample('D', kind='period')
        ex_index = period_range('2001-09-20', periods=1, freq='D')
        expected = Series([1.5], index=ex_index)
        assert_series_equal(result, expected)

        # GH 6397
        # comparing an offset that doesn't propogate tz's
        rng = date_range('1/1/2011', periods=20000, freq='H')
        rng = rng.tz_localize('EST')
        ts = DataFrame(index=rng)
        ts['first']=np.random.randn(len(rng))
        ts['second']=np.cumsum(np.random.randn(len(rng)))
        expected = DataFrame({ 'first' : ts.resample('A',how=np.sum)['first'],
                               'second' : ts.resample('A',how=np.mean)['second'] },columns=['first','second'])
        result = ts.resample('A', how={'first':np.sum, 'second':np.mean}).reindex(columns=['first','second'])
        assert_frame_equal(result,expected)

    def test_closed_left_corner(self):
        # #1465
        s = Series(np.random.randn(21),
                   index=date_range(start='1/1/2012 9:30',
                                    freq='1min', periods=21))
        s[0] = np.nan

        result = s.resample('10min', how='mean', closed='left', label='right')
        exp = s[1:].resample('10min', how='mean', closed='left', label='right')
        assert_series_equal(result, exp)

        result = s.resample('10min', how='mean', closed='left', label='left')
        exp = s[1:].resample('10min', how='mean', closed='left', label='left')

        ex_index = date_range(start='1/1/2012 9:30', freq='10min', periods=3)

        self.assertTrue(result.index.equals(ex_index))
        assert_series_equal(result, exp)

    def test_quarterly_resampling(self):
        rng = period_range('2000Q1', periods=10, freq='Q-DEC')
        ts = Series(np.arange(10), index=rng)

        result = ts.resample('A')
        exp = ts.to_timestamp().resample('A').to_period()
        assert_series_equal(result, exp)

    def test_resample_weekly_bug_1726(self):
        # 8/6/12 is a Monday
        ind = DatetimeIndex(start="8/6/2012", end="8/26/2012", freq="D")
        n = len(ind)
        data = [[x] * 5 for x in range(n)]
        df = DataFrame(data, columns=['open', 'high', 'low', 'close', 'vol'],
                       index=ind)

        # it works!
        df.resample('W-MON', how='first', closed='left', label='left')

    def test_resample_bms_2752(self):
        # GH2753
        foo = pd.Series(index=pd.bdate_range('20000101','20000201'))
        res1 = foo.resample("BMS")
        res2 = foo.resample("BMS").resample("B")
        self.assertEqual(res1.index[0], Timestamp('20000103'))
        self.assertEqual(res1.index[0], res2.index[0])

    # def test_monthly_convention_span(self):
    #     rng = period_range('2000-01', periods=3, freq='M')
    #     ts = Series(np.arange(3), index=rng)

    #     # hacky way to get same thing
    #     exp_index = period_range('2000-01-01', '2000-03-31', freq='D')
    #     expected = ts.asfreq('D', how='end').reindex(exp_index)
    #     expected = expected.fillna(method='bfill')

    #     result = ts.resample('D', convention='span')

    #     assert_series_equal(result, expected)

    def test_default_right_closed_label(self):
        end_freq = ['D', 'Q', 'M', 'D']
        end_types = ['M', 'A', 'Q', 'W']

        for from_freq, to_freq in zip(end_freq, end_types):
            idx = DatetimeIndex(start='8/15/2012', periods=100,
                                freq=from_freq)
            df = DataFrame(np.random.randn(len(idx), 2), idx)

            resampled = df.resample(to_freq)
            assert_frame_equal(resampled, df.resample(to_freq, closed='right',
                                                      label='right'))

    def test_default_left_closed_label(self):
        others = ['MS', 'AS', 'QS', 'D', 'H']
        others_freq = ['D', 'Q', 'M', 'H', 'T']

        for from_freq, to_freq in zip(others_freq, others):
            idx = DatetimeIndex(start='8/15/2012', periods=100,
                                freq=from_freq)
            df = DataFrame(np.random.randn(len(idx), 2), idx)

            resampled = df.resample(to_freq)
            assert_frame_equal(resampled, df.resample(to_freq, closed='left',
                                                      label='left'))

    def test_all_values_single_bin(self):
        # 2070
        index = period_range(start="2012-01-01", end="2012-12-31", freq="M")
        s = Series(np.random.randn(len(index)), index=index)

        result = s.resample("A", how='mean')
        tm.assert_almost_equal(result[0], s.mean())

    def test_evenly_divisible_with_no_extra_bins(self):
        # 4076
        # when the frequency is evenly divisible, sometimes extra bins

        df = DataFrame(np.random.randn(9, 3), index=date_range('2000-1-1', periods=9))
        result = df.resample('5D')
        expected = pd.concat([df.iloc[0:5].mean(),df.iloc[5:].mean()],axis=1).T
        expected.index = [Timestamp('2000-1-1'),Timestamp('2000-1-6')]
        assert_frame_equal(result,expected)

        index = date_range(start='2001-5-4', periods=28)
        df = DataFrame(
            [{'REST_KEY': 1, 'DLY_TRN_QT': 80, 'DLY_SLS_AMT': 90,
              'COOP_DLY_TRN_QT': 30, 'COOP_DLY_SLS_AMT': 20}] * 28 +
            [{'REST_KEY': 2, 'DLY_TRN_QT': 70, 'DLY_SLS_AMT': 10,
              'COOP_DLY_TRN_QT': 50, 'COOP_DLY_SLS_AMT': 20}] * 28,
            index=index.append(index)).sort_index()

        index = date_range('2001-5-4',periods=4,freq='7D')
        expected = DataFrame(
            [{'REST_KEY': 14, 'DLY_TRN_QT': 14, 'DLY_SLS_AMT': 14,
              'COOP_DLY_TRN_QT': 14, 'COOP_DLY_SLS_AMT': 14}] * 4,
            index=index)
        result = df.resample('7D', how='count')
        assert_frame_equal(result,expected)

        expected = DataFrame(
            [{'REST_KEY': 21, 'DLY_TRN_QT': 1050, 'DLY_SLS_AMT': 700,
              'COOP_DLY_TRN_QT': 560, 'COOP_DLY_SLS_AMT': 280}] * 4,
            index=index)
        result = df.resample('7D', how='sum')
        assert_frame_equal(result,expected)

class TestTimeGrouper(tm.TestCase):

    def setUp(self):
        self.ts = Series(np.random.randn(1000),
                         index=date_range('1/1/2000', periods=1000))

    def test_apply(self):
        grouper = TimeGrouper('A', label='right', closed='right')

        grouped = self.ts.groupby(grouper)

        f = lambda x: x.sort_values()[-3:]

        applied = grouped.apply(f)
        expected = self.ts.groupby(lambda x: x.year).apply(f)

        applied.index = applied.index.droplevel(0)
        expected.index = expected.index.droplevel(0)
        assert_series_equal(applied, expected)

    def test_count(self):
        self.ts[::3] = np.nan

        grouper = TimeGrouper('A', label='right', closed='right')
        result = self.ts.resample('A', how='count')

        expected = self.ts.groupby(lambda x: x.year).count()
        expected.index = result.index

        assert_series_equal(result, expected)

    def test_numpy_reduction(self):
        result = self.ts.resample('A', how='prod', closed='right')

        expected = self.ts.groupby(lambda x: x.year).agg(np.prod)
        expected.index = result.index

        assert_series_equal(result, expected)

    def test_apply_iteration(self):
        # #2300
        N = 1000
        ind = pd.date_range(start="2000-01-01", freq="D", periods=N)
        df = DataFrame({'open': 1, 'close': 2}, index=ind)
        tg = TimeGrouper('M')

        _, grouper, _ = tg._get_grouper(df)

        # Errors
        grouped = df.groupby(grouper, group_keys=False)
        f = lambda df: df['close'] / df['open']

        # it works!
        result = grouped.apply(f)
        self.assertTrue(result.index.equals(df.index))

    def test_panel_aggregation(self):
        ind = pd.date_range('1/1/2000', periods=100)
        data = np.random.randn(2, len(ind), 4)
        wp = pd.Panel(data, items=['Item1', 'Item2'], major_axis=ind,
                      minor_axis=['A', 'B', 'C', 'D'])

        tg = TimeGrouper('M', axis=1)
        _, grouper, _ = tg._get_grouper(wp)
        bingrouped = wp.groupby(grouper)
        binagg = bingrouped.mean()

        def f(x):
            assert(isinstance(x, Panel))
            return x.mean(1)
        result = bingrouped.agg(f)
        tm.assert_panel_equal(result, binagg)

    def test_fails_on_no_datetime_index(self):
        index_names = ('Int64Index', 'PeriodIndex', 'Index', 'Float64Index',
                       'MultiIndex')
        index_funcs = (tm.makeIntIndex, tm.makePeriodIndex,
                       tm.makeUnicodeIndex, tm.makeFloatIndex,
                       lambda m: tm.makeCustomIndex(m, 2))
        n = 2
        for name, func in zip(index_names, index_funcs):
            index = func(n)
            df = DataFrame({'a': np.random.randn(n)}, index=index)
            with tm.assertRaisesRegexp(TypeError,
                                       "axis must be a DatetimeIndex, "
                                       "but got an instance of %r" % name):
                df.groupby(TimeGrouper('D'))

    def test_aggregate_normal(self):
        # check TimeGrouper's aggregation is identical as normal groupby

        n = 20
        data = np.random.randn(n, 4)
        normal_df = DataFrame(data, columns=['A', 'B', 'C', 'D'])
        normal_df['key'] = [1, 2, 3, 4, 5] * 4

        dt_df = DataFrame(data, columns=['A', 'B', 'C', 'D'])
        dt_df['key'] = [datetime(2013, 1, 1), datetime(2013, 1, 2), datetime(2013, 1, 3),
                        datetime(2013, 1, 4), datetime(2013, 1, 5)] * 4

        normal_grouped = normal_df.groupby('key')
        dt_grouped = dt_df.groupby(TimeGrouper(key='key', freq='D'))

        for func in ['min', 'max', 'prod', 'var', 'std', 'mean']:
            expected = getattr(normal_grouped, func)()
            dt_result = getattr(dt_grouped, func)()
            expected.index = date_range(start='2013-01-01', freq='D', periods=5, name='key')
            assert_frame_equal(expected, dt_result)

        for func in ['count', 'sum']:
            expected = getattr(normal_grouped, func)()
            expected.index = date_range(start='2013-01-01', freq='D', periods=5, name='key')
            dt_result = getattr(dt_grouped, func)()
            assert_frame_equal(expected, dt_result)

        # GH 7453
        for func in ['size']:
            expected = getattr(normal_grouped, func)()
            expected.index = date_range(start='2013-01-01', freq='D', periods=5, name='key')
            dt_result = getattr(dt_grouped, func)()
            assert_series_equal(expected, dt_result)

        """
        for func in ['first', 'last']:
            expected = getattr(normal_grouped, func)()
            expected.index = date_range(start='2013-01-01', freq='D', periods=5, name='key')
            dt_result = getattr(dt_grouped, func)()
            assert_frame_equal(expected, dt_result)

        for func in ['nth']:
            expected = getattr(normal_grouped, func)(3)
            expected.index = date_range(start='2013-01-01', freq='D', periods=5, name='key')
            dt_result = getattr(dt_grouped, func)(3)
            assert_frame_equal(expected, dt_result)
        """
        # if TimeGrouper is used included, 'first','last' and 'nth' doesn't work yet

    def test_aggregate_with_nat(self):
        # check TimeGrouper's aggregation is identical as normal groupby

        n = 20
        data = np.random.randn(n, 4).astype('int64')
        normal_df = DataFrame(data, columns=['A', 'B', 'C', 'D'])
        normal_df['key'] = [1, 2, np.nan, 4, 5] * 4

        dt_df = DataFrame(data, columns=['A', 'B', 'C', 'D'])
        dt_df['key'] = [datetime(2013, 1, 1), datetime(2013, 1, 2), pd.NaT,
                        datetime(2013, 1, 4), datetime(2013, 1, 5)] * 4

        normal_grouped = normal_df.groupby('key')
        dt_grouped = dt_df.groupby(TimeGrouper(key='key', freq='D'))

        for func in ['min', 'max', 'sum', 'prod']:
            normal_result = getattr(normal_grouped, func)()
            dt_result = getattr(dt_grouped, func)()
            pad = DataFrame([[np.nan, np.nan, np.nan, np.nan]],
                            index=[3], columns=['A', 'B', 'C', 'D'])
            expected = normal_result.append(pad)
            expected = expected.sort_index()
            expected.index = date_range(start='2013-01-01', freq='D', periods=5, name='key')
            assert_frame_equal(expected, dt_result)

        for func in ['count']:
            normal_result = getattr(normal_grouped, func)()
            pad = DataFrame([[0, 0, 0, 0]], index=[3], columns=['A', 'B', 'C', 'D'])
            expected = normal_result.append(pad)
            expected = expected.sort_index()
            expected.index = date_range(start='2013-01-01', freq='D', periods=5, name='key')
            dt_result = getattr(dt_grouped, func)()
            assert_frame_equal(expected, dt_result)

        for func in ['size']:
            normal_result = getattr(normal_grouped, func)()
            pad = Series([0], index=[3])
            expected = normal_result.append(pad)
            expected = expected.sort_index()
            expected.index = date_range(start='2013-01-01', freq='D', periods=5, name='key')
            dt_result = getattr(dt_grouped, func)()
            assert_series_equal(expected, dt_result)
            # GH 9925
            self.assertEqual(dt_result.index.name, 'key')

        # if NaT is included, 'var', 'std', 'mean', 'first','last' and 'nth' doesn't work yet


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
