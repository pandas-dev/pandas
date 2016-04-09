# coding=utf-8
# pylint: disable-msg=E1101,W0612

from datetime import datetime

from numpy import nan
import numpy as np

from pandas import Index, Series, notnull, date_range
from pandas.tseries.index import DatetimeIndex
from pandas.tseries.tdi import TimedeltaIndex

import pandas.core.datetools as datetools

from pandas.util.testing import assert_series_equal, assert_almost_equal
import pandas.util.testing as tm

from .common import TestData


class TestSeriesTimeSeries(TestData, tm.TestCase):
    _multiprocess_can_split_ = True

    def test_shift(self):
        shifted = self.ts.shift(1)
        unshifted = shifted.shift(-1)

        tm.assert_dict_equal(unshifted.valid(), self.ts, compare_keys=False)

        offset = datetools.bday
        shifted = self.ts.shift(1, freq=offset)
        unshifted = shifted.shift(-1, freq=offset)

        assert_series_equal(unshifted, self.ts)

        unshifted = self.ts.shift(0, freq=offset)
        assert_series_equal(unshifted, self.ts)

        shifted = self.ts.shift(1, freq='B')
        unshifted = shifted.shift(-1, freq='B')

        assert_series_equal(unshifted, self.ts)

        # corner case
        unshifted = self.ts.shift(0)
        assert_series_equal(unshifted, self.ts)

        # Shifting with PeriodIndex
        ps = tm.makePeriodSeries()
        shifted = ps.shift(1)
        unshifted = shifted.shift(-1)
        tm.assert_dict_equal(unshifted.valid(), ps, compare_keys=False)

        shifted2 = ps.shift(1, 'B')
        shifted3 = ps.shift(1, datetools.bday)
        assert_series_equal(shifted2, shifted3)
        assert_series_equal(ps, shifted2.shift(-1, 'B'))

        self.assertRaises(ValueError, ps.shift, freq='D')

        # legacy support
        shifted4 = ps.shift(1, freq='B')
        assert_series_equal(shifted2, shifted4)

        shifted5 = ps.shift(1, freq=datetools.bday)
        assert_series_equal(shifted5, shifted4)

        # 32-bit taking
        # GH 8129
        index = date_range('2000-01-01', periods=5)
        for dtype in ['int32', 'int64']:
            s1 = Series(np.arange(5, dtype=dtype), index=index)
            p = s1.iloc[1]
            result = s1.shift(periods=p)
            expected = Series([np.nan, 0, 1, 2, 3], index=index)
            assert_series_equal(result, expected)

        # xref 8260
        # with tz
        s = Series(
            date_range('2000-01-01 09:00:00', periods=5,
                       tz='US/Eastern'), name='foo')
        result = s - s.shift()
        assert_series_equal(result, Series(
            TimedeltaIndex(['NaT'] + ['1 days'] * 4), name='foo'))

        # incompat tz
        s2 = Series(
            date_range('2000-01-01 09:00:00', periods=5, tz='CET'), name='foo')
        self.assertRaises(ValueError, lambda: s - s2)

    def test_tshift(self):
        # PeriodIndex
        ps = tm.makePeriodSeries()
        shifted = ps.tshift(1)
        unshifted = shifted.tshift(-1)

        assert_series_equal(unshifted, ps)

        shifted2 = ps.tshift(freq='B')
        assert_series_equal(shifted, shifted2)

        shifted3 = ps.tshift(freq=datetools.bday)
        assert_series_equal(shifted, shifted3)

        self.assertRaises(ValueError, ps.tshift, freq='M')

        # DatetimeIndex
        shifted = self.ts.tshift(1)
        unshifted = shifted.tshift(-1)

        assert_series_equal(self.ts, unshifted)

        shifted2 = self.ts.tshift(freq=self.ts.index.freq)
        assert_series_equal(shifted, shifted2)

        inferred_ts = Series(self.ts.values, Index(np.asarray(self.ts.index)),
                             name='ts')
        shifted = inferred_ts.tshift(1)
        unshifted = shifted.tshift(-1)
        assert_series_equal(shifted, self.ts.tshift(1))
        assert_series_equal(unshifted, inferred_ts)

        no_freq = self.ts[[0, 5, 7]]
        self.assertRaises(ValueError, no_freq.tshift)

    def test_truncate(self):
        offset = datetools.bday

        ts = self.ts[::3]

        start, end = self.ts.index[3], self.ts.index[6]
        start_missing, end_missing = self.ts.index[2], self.ts.index[7]

        # neither specified
        truncated = ts.truncate()
        assert_series_equal(truncated, ts)

        # both specified
        expected = ts[1:3]

        truncated = ts.truncate(start, end)
        assert_series_equal(truncated, expected)

        truncated = ts.truncate(start_missing, end_missing)
        assert_series_equal(truncated, expected)

        # start specified
        expected = ts[1:]

        truncated = ts.truncate(before=start)
        assert_series_equal(truncated, expected)

        truncated = ts.truncate(before=start_missing)
        assert_series_equal(truncated, expected)

        # end specified
        expected = ts[:3]

        truncated = ts.truncate(after=end)
        assert_series_equal(truncated, expected)

        truncated = ts.truncate(after=end_missing)
        assert_series_equal(truncated, expected)

        # corner case, empty series returned
        truncated = ts.truncate(after=self.ts.index[0] - offset)
        assert (len(truncated) == 0)

        truncated = ts.truncate(before=self.ts.index[-1] + offset)
        assert (len(truncated) == 0)

        self.assertRaises(ValueError, ts.truncate,
                          before=self.ts.index[-1] + offset,
                          after=self.ts.index[0] - offset)

    def test_asof(self):
        # array or list or dates
        N = 50
        rng = date_range('1/1/1990', periods=N, freq='53s')
        ts = Series(np.random.randn(N), index=rng)
        ts[15:30] = np.nan
        dates = date_range('1/1/1990', periods=N * 3, freq='25s')

        result = ts.asof(dates)
        self.assertTrue(notnull(result).all())
        lb = ts.index[14]
        ub = ts.index[30]

        result = ts.asof(list(dates))
        self.assertTrue(notnull(result).all())
        lb = ts.index[14]
        ub = ts.index[30]

        mask = (result.index >= lb) & (result.index < ub)
        rs = result[mask]
        self.assertTrue((rs == ts[lb]).all())

        val = result[result.index[result.index >= ub][0]]
        self.assertEqual(ts[ub], val)

        self.ts[5:10] = np.NaN
        self.ts[15:20] = np.NaN

        val1 = self.ts.asof(self.ts.index[7])
        val2 = self.ts.asof(self.ts.index[19])

        self.assertEqual(val1, self.ts[4])
        self.assertEqual(val2, self.ts[14])

        # accepts strings
        val1 = self.ts.asof(str(self.ts.index[7]))
        self.assertEqual(val1, self.ts[4])

        # in there
        self.assertEqual(self.ts.asof(self.ts.index[3]), self.ts[3])

        # no as of value
        d = self.ts.index[0] - datetools.bday
        self.assertTrue(np.isnan(self.ts.asof(d)))

    def test_getitem_setitem_datetimeindex(self):
        from pandas import date_range

        N = 50
        # testing with timezone, GH #2785
        rng = date_range('1/1/1990', periods=N, freq='H', tz='US/Eastern')
        ts = Series(np.random.randn(N), index=rng)

        result = ts["1990-01-01 04:00:00"]
        expected = ts[4]
        self.assertEqual(result, expected)

        result = ts.copy()
        result["1990-01-01 04:00:00"] = 0
        result["1990-01-01 04:00:00"] = ts[4]
        assert_series_equal(result, ts)

        result = ts["1990-01-01 04:00:00":"1990-01-01 07:00:00"]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        result = ts.copy()
        result["1990-01-01 04:00:00":"1990-01-01 07:00:00"] = 0
        result["1990-01-01 04:00:00":"1990-01-01 07:00:00"] = ts[4:8]
        assert_series_equal(result, ts)

        lb = "1990-01-01 04:00:00"
        rb = "1990-01-01 07:00:00"
        result = ts[(ts.index >= lb) & (ts.index <= rb)]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        # repeat all the above with naive datetimes
        result = ts[datetime(1990, 1, 1, 4)]
        expected = ts[4]
        self.assertEqual(result, expected)

        result = ts.copy()
        result[datetime(1990, 1, 1, 4)] = 0
        result[datetime(1990, 1, 1, 4)] = ts[4]
        assert_series_equal(result, ts)

        result = ts[datetime(1990, 1, 1, 4):datetime(1990, 1, 1, 7)]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        result = ts.copy()
        result[datetime(1990, 1, 1, 4):datetime(1990, 1, 1, 7)] = 0
        result[datetime(1990, 1, 1, 4):datetime(1990, 1, 1, 7)] = ts[4:8]
        assert_series_equal(result, ts)

        lb = datetime(1990, 1, 1, 4)
        rb = datetime(1990, 1, 1, 7)
        result = ts[(ts.index >= lb) & (ts.index <= rb)]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        result = ts[ts.index[4]]
        expected = ts[4]
        self.assertEqual(result, expected)

        result = ts[ts.index[4:8]]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        result = ts.copy()
        result[ts.index[4:8]] = 0
        result[4:8] = ts[4:8]
        assert_series_equal(result, ts)

        # also test partial date slicing
        result = ts["1990-01-02"]
        expected = ts[24:48]
        assert_series_equal(result, expected)

        result = ts.copy()
        result["1990-01-02"] = 0
        result["1990-01-02"] = ts[24:48]
        assert_series_equal(result, ts)

    def test_getitem_setitem_datetime_tz_pytz(self):
        tm._skip_if_no_pytz()
        from pytz import timezone as tz

        from pandas import date_range

        N = 50
        # testing with timezone, GH #2785
        rng = date_range('1/1/1990', periods=N, freq='H', tz='US/Eastern')
        ts = Series(np.random.randn(N), index=rng)

        # also test Timestamp tz handling, GH #2789
        result = ts.copy()
        result["1990-01-01 09:00:00+00:00"] = 0
        result["1990-01-01 09:00:00+00:00"] = ts[4]
        assert_series_equal(result, ts)

        result = ts.copy()
        result["1990-01-01 03:00:00-06:00"] = 0
        result["1990-01-01 03:00:00-06:00"] = ts[4]
        assert_series_equal(result, ts)

        # repeat with datetimes
        result = ts.copy()
        result[datetime(1990, 1, 1, 9, tzinfo=tz('UTC'))] = 0
        result[datetime(1990, 1, 1, 9, tzinfo=tz('UTC'))] = ts[4]
        assert_series_equal(result, ts)

        result = ts.copy()

        # comparison dates with datetime MUST be localized!
        date = tz('US/Central').localize(datetime(1990, 1, 1, 3))
        result[date] = 0
        result[date] = ts[4]
        assert_series_equal(result, ts)

    def test_getitem_setitem_datetime_tz_dateutil(self):
        tm._skip_if_no_dateutil()
        from dateutil.tz import tzutc
        from pandas.tslib import _dateutil_gettz as gettz

        tz = lambda x: tzutc() if x == 'UTC' else gettz(
            x)  # handle special case for utc in dateutil

        from pandas import date_range

        N = 50
        # testing with timezone, GH #2785
        rng = date_range('1/1/1990', periods=N, freq='H', tz='US/Eastern')
        ts = Series(np.random.randn(N), index=rng)

        # also test Timestamp tz handling, GH #2789
        result = ts.copy()
        result["1990-01-01 09:00:00+00:00"] = 0
        result["1990-01-01 09:00:00+00:00"] = ts[4]
        assert_series_equal(result, ts)

        result = ts.copy()
        result["1990-01-01 03:00:00-06:00"] = 0
        result["1990-01-01 03:00:00-06:00"] = ts[4]
        assert_series_equal(result, ts)

        # repeat with datetimes
        result = ts.copy()
        result[datetime(1990, 1, 1, 9, tzinfo=tz('UTC'))] = 0
        result[datetime(1990, 1, 1, 9, tzinfo=tz('UTC'))] = ts[4]
        assert_series_equal(result, ts)

        result = ts.copy()
        result[datetime(1990, 1, 1, 3, tzinfo=tz('US/Central'))] = 0
        result[datetime(1990, 1, 1, 3, tzinfo=tz('US/Central'))] = ts[4]
        assert_series_equal(result, ts)

    def test_getitem_setitem_periodindex(self):
        from pandas import period_range

        N = 50
        rng = period_range('1/1/1990', periods=N, freq='H')
        ts = Series(np.random.randn(N), index=rng)

        result = ts["1990-01-01 04"]
        expected = ts[4]
        self.assertEqual(result, expected)

        result = ts.copy()
        result["1990-01-01 04"] = 0
        result["1990-01-01 04"] = ts[4]
        assert_series_equal(result, ts)

        result = ts["1990-01-01 04":"1990-01-01 07"]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        result = ts.copy()
        result["1990-01-01 04":"1990-01-01 07"] = 0
        result["1990-01-01 04":"1990-01-01 07"] = ts[4:8]
        assert_series_equal(result, ts)

        lb = "1990-01-01 04"
        rb = "1990-01-01 07"
        result = ts[(ts.index >= lb) & (ts.index <= rb)]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        # GH 2782
        result = ts[ts.index[4]]
        expected = ts[4]
        self.assertEqual(result, expected)

        result = ts[ts.index[4:8]]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        result = ts.copy()
        result[ts.index[4:8]] = 0
        result[4:8] = ts[4:8]
        assert_series_equal(result, ts)

    def test_asof_periodindex(self):
        from pandas import period_range, PeriodIndex
        # array or list or dates
        N = 50
        rng = period_range('1/1/1990', periods=N, freq='H')
        ts = Series(np.random.randn(N), index=rng)
        ts[15:30] = np.nan
        dates = date_range('1/1/1990', periods=N * 3, freq='37min')

        result = ts.asof(dates)
        self.assertTrue(notnull(result).all())
        lb = ts.index[14]
        ub = ts.index[30]

        result = ts.asof(list(dates))
        self.assertTrue(notnull(result).all())
        lb = ts.index[14]
        ub = ts.index[30]

        pix = PeriodIndex(result.index.values, freq='H')
        mask = (pix >= lb) & (pix < ub)
        rs = result[mask]
        self.assertTrue((rs == ts[lb]).all())

        ts[5:10] = np.NaN
        ts[15:20] = np.NaN

        val1 = ts.asof(ts.index[7])
        val2 = ts.asof(ts.index[19])

        self.assertEqual(val1, ts[4])
        self.assertEqual(val2, ts[14])

        # accepts strings
        val1 = ts.asof(str(ts.index[7]))
        self.assertEqual(val1, ts[4])

        # in there
        self.assertEqual(ts.asof(ts.index[3]), ts[3])

        # no as of value
        d = ts.index[0].to_timestamp() - datetools.bday
        self.assertTrue(np.isnan(ts.asof(d)))

    def test_asof_more(self):
        from pandas import date_range

        s = Series([nan, nan, 1, 2, nan, nan, 3, 4, 5],
                   index=date_range('1/1/2000', periods=9))

        dates = s.index[[4, 5, 6, 2, 1]]

        result = s.asof(dates)
        expected = Series([2, 2, 3, 1, np.nan], index=dates)

        assert_series_equal(result, expected)

        s = Series([1.5, 2.5, 1, 2, nan, nan, 3, 4, 5],
                   index=date_range('1/1/2000', periods=9))
        result = s.asof(s.index[0])
        self.assertEqual(result, s[0])

    def test_asfreq(self):
        ts = Series([0., 1., 2.], index=[datetime(2009, 10, 30), datetime(
            2009, 11, 30), datetime(2009, 12, 31)])

        daily_ts = ts.asfreq('B')
        monthly_ts = daily_ts.asfreq('BM')
        self.assert_numpy_array_equal(monthly_ts, ts)

        daily_ts = ts.asfreq('B', method='pad')
        monthly_ts = daily_ts.asfreq('BM')
        self.assert_numpy_array_equal(monthly_ts, ts)

        daily_ts = ts.asfreq(datetools.bday)
        monthly_ts = daily_ts.asfreq(datetools.bmonthEnd)
        self.assert_numpy_array_equal(monthly_ts, ts)

        result = ts[:0].asfreq('M')
        self.assertEqual(len(result), 0)
        self.assertIsNot(result, ts)

    def test_diff(self):
        # Just run the function
        self.ts.diff()

        # int dtype
        a = 10000000000000000
        b = a + 1
        s = Series([a, b])

        rs = s.diff()
        self.assertEqual(rs[1], 1)

        # neg n
        rs = self.ts.diff(-1)
        xp = self.ts - self.ts.shift(-1)
        assert_series_equal(rs, xp)

        # 0
        rs = self.ts.diff(0)
        xp = self.ts - self.ts
        assert_series_equal(rs, xp)

        # datetime diff (GH3100)
        s = Series(date_range('20130102', periods=5))
        rs = s - s.shift(1)
        xp = s.diff()
        assert_series_equal(rs, xp)

        # timedelta diff
        nrs = rs - rs.shift(1)
        nxp = xp.diff()
        assert_series_equal(nrs, nxp)

        # with tz
        s = Series(
            date_range('2000-01-01 09:00:00', periods=5,
                       tz='US/Eastern'), name='foo')
        result = s.diff()
        assert_series_equal(result, Series(
            TimedeltaIndex(['NaT'] + ['1 days'] * 4), name='foo'))

    def test_pct_change(self):
        rs = self.ts.pct_change(fill_method=None)
        assert_series_equal(rs, self.ts / self.ts.shift(1) - 1)

        rs = self.ts.pct_change(2)
        filled = self.ts.fillna(method='pad')
        assert_series_equal(rs, filled / filled.shift(2) - 1)

        rs = self.ts.pct_change(fill_method='bfill', limit=1)
        filled = self.ts.fillna(method='bfill', limit=1)
        assert_series_equal(rs, filled / filled.shift(1) - 1)

        rs = self.ts.pct_change(freq='5D')
        filled = self.ts.fillna(method='pad')
        assert_series_equal(rs, filled / filled.shift(freq='5D') - 1)

    def test_pct_change_shift_over_nas(self):
        s = Series([1., 1.5, np.nan, 2.5, 3.])

        chg = s.pct_change()
        expected = Series([np.nan, 0.5, np.nan, 2.5 / 1.5 - 1, .2])
        assert_series_equal(chg, expected)

    def test_autocorr(self):
        # Just run the function
        corr1 = self.ts.autocorr()

        # Now run it with the lag parameter
        corr2 = self.ts.autocorr(lag=1)

        # corr() with lag needs Series of at least length 2
        if len(self.ts) <= 2:
            self.assertTrue(np.isnan(corr1))
            self.assertTrue(np.isnan(corr2))
        else:
            self.assertEqual(corr1, corr2)

        # Choose a random lag between 1 and length of Series - 2
        # and compare the result with the Series corr() function
        n = 1 + np.random.randint(max(1, len(self.ts) - 2))
        corr1 = self.ts.corr(self.ts.shift(n))
        corr2 = self.ts.autocorr(lag=n)

        # corr() with lag needs Series of at least length 2
        if len(self.ts) <= 2:
            self.assertTrue(np.isnan(corr1))
            self.assertTrue(np.isnan(corr2))
        else:
            self.assertEqual(corr1, corr2)

    def test_first_last_valid(self):
        ts = self.ts.copy()
        ts[:5] = np.NaN

        index = ts.first_valid_index()
        self.assertEqual(index, ts.index[5])

        ts[-5:] = np.NaN
        index = ts.last_valid_index()
        self.assertEqual(index, ts.index[-6])

        ts[:] = np.nan
        self.assertIsNone(ts.last_valid_index())
        self.assertIsNone(ts.first_valid_index())

        ser = Series([], index=[])
        self.assertIsNone(ser.last_valid_index())
        self.assertIsNone(ser.first_valid_index())

        # GH12800
        empty = Series()
        self.assertIsNone(empty.last_valid_index())
        self.assertIsNone(empty.first_valid_index())

    def test_mpl_compat_hack(self):
        result = self.ts[:, np.newaxis]
        expected = self.ts.values[:, np.newaxis]
        assert_almost_equal(result, expected)

    def test_timeseries_coercion(self):
        idx = tm.makeDateIndex(10000)
        ser = Series(np.random.randn(len(idx)), idx.astype(object))
        with tm.assert_produces_warning(FutureWarning):
            self.assertTrue(ser.is_time_series)
        self.assertTrue(ser.index.is_all_dates)
        self.assertIsInstance(ser.index, DatetimeIndex)
