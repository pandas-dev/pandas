# pylint: disable=E1101

from datetime import datetime, timedelta

import numpy as np

from pandas import Series, TimeSeries, DataFrame, Panel, isnull, notnull, Timestamp

from pandas.tseries.index import date_range
from pandas.tseries.offsets import Minute, BDay
from pandas.tseries.period import period_range, PeriodIndex, Period
from pandas.tseries.resample import DatetimeIndex, TimeGrouper
import pandas.tseries.offsets as offsets
import pandas as pd

import unittest
import nose

from pandas.util.testing import (assert_series_equal, assert_almost_equal,
                                 assert_frame_equal)
import pandas.util.testing as tm

bday = BDay()


def _skip_if_no_pytz():
    try:
        import pytz
    except ImportError:
        raise nose.SkipTest


class TestResample(unittest.TestCase):
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

        self.assertEquals(g.ngroups, 2593)
        self.assert_(notnull(g.mean()).all())

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

        self.assertEquals(len(r.columns), 10)
        self.assertEquals(len(r.index), 2593)

    def test_resample_basic(self):
        rng = date_range('1/1/2000 00:00:00', '1/1/2000 00:13:00', freq='min',
                         name='index')
        s = Series(np.random.randn(14), index=rng)
        result = s.resample('5min', how='mean', closed='right', label='right')
        expected = Series([s[0], s[1:6].mean(), s[6:11].mean(), s[11:].mean()],
                          index=date_range('1/1/2000', periods=4, freq='5min'))
        assert_series_equal(result, expected)
        self.assert_(result.index.name == 'index')

        result = s.resample('5min', how='mean', closed='left', label='right')
        expected = Series([s[:5].mean(), s[5:10].mean(), s[10:].mean()],
                          index=date_range('1/1/2000 00:05', periods=3,
                                           freq='5min'))
        assert_series_equal(result, expected)

        s = self.series
        result = s.resample('5Min', how='last')
        grouper = TimeGrouper(Minute(5), closed='left', label='left')
        expect = s.groupby(grouper).agg(lambda x: x[-1])
        assert_series_equal(result, expect)

    def test_resample_basic_from_daily(self):
        # from daily
        dti = DatetimeIndex(
            start=datetime(2005, 1, 1), end=datetime(2005, 1, 10),
            freq='D', name='index')

        s = Series(np.random.rand(len(dti)), dti)

        # to weekly
        result = s.resample('w-sun', how='last')

        self.assertEquals(len(result), 3)
        self.assert_((result.index.dayofweek == [6, 6, 6]).all())
        self.assertEquals(result.irow(0), s['1/2/2005'])
        self.assertEquals(result.irow(1), s['1/9/2005'])
        self.assertEquals(result.irow(2), s.irow(-1))

        result = s.resample('W-MON', how='last')
        self.assertEquals(len(result), 2)
        self.assert_((result.index.dayofweek == [0, 0]).all())
        self.assertEquals(result.irow(0), s['1/3/2005'])
        self.assertEquals(result.irow(1), s['1/10/2005'])

        result = s.resample('W-TUE', how='last')
        self.assertEquals(len(result), 2)
        self.assert_((result.index.dayofweek == [1, 1]).all())
        self.assertEquals(result.irow(0), s['1/4/2005'])
        self.assertEquals(result.irow(1), s['1/10/2005'])

        result = s.resample('W-WED', how='last')
        self.assertEquals(len(result), 2)
        self.assert_((result.index.dayofweek == [2, 2]).all())
        self.assertEquals(result.irow(0), s['1/5/2005'])
        self.assertEquals(result.irow(1), s['1/10/2005'])

        result = s.resample('W-THU', how='last')
        self.assertEquals(len(result), 2)
        self.assert_((result.index.dayofweek == [3, 3]).all())
        self.assertEquals(result.irow(0), s['1/6/2005'])
        self.assertEquals(result.irow(1), s['1/10/2005'])

        result = s.resample('W-FRI', how='last')
        self.assertEquals(len(result), 2)
        self.assert_((result.index.dayofweek == [4, 4]).all())
        self.assertEquals(result.irow(0), s['1/7/2005'])
        self.assertEquals(result.irow(1), s['1/10/2005'])

        # to biz day
        result = s.resample('B', how='last')
        self.assertEquals(len(result), 7)
        self.assert_((result.index.dayofweek == [4, 0, 1, 2, 3, 4, 0]).all())
        self.assertEquals(result.irow(0), s['1/2/2005'])
        self.assertEquals(result.irow(1), s['1/3/2005'])
        self.assertEquals(result.irow(5), s['1/9/2005'])
        self.assert_(result.index.name == 'index')

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

        self.assert_(result.index.freq == Minute(5))

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
        self.assertEquals(len(result), 12961)
        self.assertEquals(result[0], s[0])
        self.assertEquals(result[-1], s[-1])

        self.assert_(result.index.name == 'index')

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

        self.assertEquals(len(result), len(expect))
        self.assertEquals(len(result.columns), 4)

        xs = result.irow(-2)
        self.assertEquals(xs['open'], s[-6])
        self.assertEquals(xs['high'], s[-6:-1].max())
        self.assertEquals(xs['low'], s[-6:-1].min())
        self.assertEquals(xs['close'], s[-2])

        xs = result.irow(0)
        self.assertEquals(xs['open'], s[0])
        self.assertEquals(xs['high'], s[:5].max())
        self.assertEquals(xs['low'], s[:5].min())
        self.assertEquals(xs['close'], s[4])

    def test_resample_reresample(self):
        dti = DatetimeIndex(
            start=datetime(2005, 1, 1), end=datetime(2005, 1, 10),
            freq='D')
        s = Series(np.random.rand(len(dti)), dti)
        bs = s.resample('B', closed='right', label='right')
        result = bs.resample('8H')
        self.assertEquals(len(result), 22)
        self.assert_(isinstance(result.index.freq, offsets.DateOffset))
        self.assert_(result.index.freq == offsets.Hour(8))

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

        self.assert_((resampled.ix['1/1/2000 00:00'] == ts[0]).all())

        exp = _ohlc(ts[1:31])
        self.assert_((resampled.ix['1/1/2000 00:05'] == exp).all())

        exp = _ohlc(ts['1/1/2000 5:55:01':])
        self.assert_((resampled.ix['1/1/2000 6:00:00'] == exp).all())

    def test_downsample_non_unique(self):
        rng = date_range('1/1/2000', '2/29/2000')
        rng2 = rng.repeat(5).values
        ts = Series(np.random.randn(len(rng2)), index=rng2)

        result = ts.resample('M', how='mean')

        expected = ts.groupby(lambda x: x.month).mean()
        self.assertEquals(len(result), 2)
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

    def test_resample_base(self):
        rng = date_range('1/1/2000 00:00:00', '1/1/2000 02:00', freq='s')
        ts = Series(np.random.randn(len(rng)), index=rng)

        resampled = ts.resample('5min', base=2)
        exp_rng = date_range('12/31/1999 23:57:00', '1/1/2000 01:57',
                             freq='5min')
        self.assert_(resampled.index.equals(exp_rng))

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
        self.assert_(result.index.equals(exp_index))

    def test_resample_empty(self):
        ts = _simple_ts('1/1/2000', '2/1/2000')[:0]

        result = ts.resample('A')
        self.assert_(len(result) == 0)
        self.assert_(result.index.freqstr == 'A-DEC')

        result = ts.resample('A', kind='period')
        self.assert_(len(result) == 0)
        self.assert_(result.index.freqstr == 'A-DEC')

        xp = DataFrame()
        rs = xp.resample('A')
        assert_frame_equal(xp, rs)

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
        self.assert_(len(resampled) == 1)

    def test_resample_anchored_monthstart(self):
        ts = _simple_ts('1/1/2000', '12/31/2002')

        freqs = ['MS', 'BMS', 'QS-MAR', 'AS-DEC', 'AS-JUN']

        for freq in freqs:
            result = ts.resample(freq, how='mean')

    def test_corner_cases(self):
        # miscellaneous test coverage

        rng = date_range('1/1/2000', periods=12, freq='t')
        ts = Series(np.random.randn(len(rng)), index=rng)

        result = ts.resample('5t', closed='right', label='left')
        ex_index = date_range('1999-12-31 23:55', periods=4, freq='5t')
        self.assert_(result.index.equals(ex_index))

        len0pts = _simple_pts('2007-01', '2010-05', freq='M')[:0]
        # it works
        result = len0pts.resample('A-DEC')
        self.assert_(len(result) == 0)

        # resample to periods
        ts = _simple_ts('2000-04-28', '2000-04-30 11:00', freq='h')
        result = ts.resample('M', kind='period')
        self.assert_(len(result) == 1)
        self.assert_(result.index[0] == Period('2000-04', freq='M'))

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
        self.assert_(isinstance(result, DataFrame))

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
                          how=[lambda x: x.mean(), lambda x: x.std()])

        result = ts.resample('M', how={'foo': lambda x: x.mean(),
                                       'bar': lambda x: x.std()})
        foo_exp = ts.resample('M', how='mean')
        bar_exp = ts.resample('M', how='std')

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


def _simple_ts(start, end, freq='D'):
    rng = date_range(start, end, freq=freq)
    return Series(np.random.randn(len(rng)), index=rng)


def _simple_pts(start, end, freq='D'):
    rng = period_range(start, end, freq=freq)
    return TimeSeries(np.random.randn(len(rng)), index=rng)


from pandas.tseries.frequencies import MONTHS, DAYS
from pandas.util.compat import product


class TestResamplePeriodIndex(unittest.TestCase):

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

            expected = stamps.reindex(qdates.to_timestamp('D', 'e'),
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

        s = TimeSeries(np.random.randn(4), index=rng)

        stamps = s.to_timestamp()

        filled = s.resample('A')
        expected = stamps.resample('A').to_period('A')
        assert_series_equal(filled, expected)

        filled = s.resample('A', fill_method='ffill')
        expected = stamps.resample('A', fill_method='ffill').to_period('A')
        assert_series_equal(filled, expected)

    def test_cant_fill_missing_dups(self):
        rng = PeriodIndex([2000, 2005, 2005, 2007, 2007], freq='A')
        s = TimeSeries(np.random.randn(5), index=rng)
        self.assertRaises(Exception, s.resample, 'A')

    def test_resample_5minute(self):
        rng = period_range('1/1/2000', '1/5/2000', freq='T')
        ts = TimeSeries(np.random.randn(len(rng)), index=rng)

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
        self.assert_(len(result) == 0)

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

        self.assert_(result.isnull().all())

        result = ts.resample('W-THU', fill_method='ffill')[:-1]
        expected = ts.asfreq('W-THU', method='ffill')
        assert_series_equal(result, expected)

    def test_resample_tz_localized(self):
        dr = date_range(start='2012-4-13', end='2012-5-1')
        ts = Series(range(len(dr)), dr)

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

        self.assert_(result.index.equals(ex_index))
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

    def test_resample_doesnt_truncate(self):
        """Test for issue #3020"""
        import pandas as pd
        dates = pd.date_range('01-Jan-2014','05-Jan-2014', freq='D')
        series = Series(1, index=dates)

        result = series.resample('D')
        self.assertEquals(result.index[0], dates[0])

class TestTimeGrouper(unittest.TestCase):

    def setUp(self):
        self.ts = Series(np.random.randn(1000),
                         index=date_range('1/1/2000', periods=1000))

    def test_apply(self):
        grouper = TimeGrouper('A', label='right', closed='right')

        grouped = self.ts.groupby(grouper)

        f = lambda x: x.order()[-3:]

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

        grouper = tg.get_grouper(df)

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
        grouper = tg.get_grouper(wp)
        bingrouped = wp.groupby(grouper)
        binagg = bingrouped.mean()

        def f(x):
            assert(isinstance(x, Panel))
            return x.mean(1)
        result = bingrouped.agg(f)
        tm.assert_panel_equal(result, binagg)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
