# coding=utf-8

import numpy as np
from pandas import (offsets, Series, notnull,
                    isnull, date_range, Timestamp)

import pandas.util.testing as tm

from .common import TestData


class TestSeriesAsof(TestData, tm.TestCase):

    def test_basic(self):

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

    def test_scalar(self):

        N = 30
        rng = date_range('1/1/1990', periods=N, freq='53s')
        ts = Series(np.arange(N), index=rng)
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
        result = ts.asof(ts.index[3])
        self.assertEqual(result, ts[3])

        # no as of value
        d = ts.index[0] - offsets.BDay()
        self.assertTrue(np.isnan(ts.asof(d)))

    def test_with_nan(self):
        # basic asof test
        rng = date_range('1/1/2000', '1/2/2000', freq='4h')
        s = Series(np.arange(len(rng)), index=rng)
        r = s.resample('2h').mean()

        result = r.asof(r.index)
        expected = Series([0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6.],
                          index=date_range('1/1/2000', '1/2/2000', freq='2h'))
        tm.assert_series_equal(result, expected)

        r.iloc[3:5] = np.nan
        result = r.asof(r.index)
        expected = Series([0, 0, 1, 1, 1, 1, 3, 3, 4, 4, 5, 5, 6.],
                          index=date_range('1/1/2000', '1/2/2000', freq='2h'))
        tm.assert_series_equal(result, expected)

        r.iloc[-3:] = np.nan
        result = r.asof(r.index)
        expected = Series([0, 0, 1, 1, 1, 1, 3, 3, 4, 4, 4, 4, 4.],
                          index=date_range('1/1/2000', '1/2/2000', freq='2h'))
        tm.assert_series_equal(result, expected)

    def test_periodindex(self):
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

        ts[5:10] = np.nan
        ts[15:20] = np.nan

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
        d = ts.index[0].to_timestamp() - offsets.BDay()
        self.assertTrue(isnull(ts.asof(d)))

    def test_errors(self):

        s = Series([1, 2, 3],
                   index=[Timestamp('20130101'),
                          Timestamp('20130103'),
                          Timestamp('20130102')])

        # non-monotonic
        self.assertFalse(s.index.is_monotonic)
        with self.assertRaises(ValueError):
            s.asof(s.index[0])

        # subset with Series
        N = 10
        rng = date_range('1/1/1990', periods=N, freq='53s')
        s = Series(np.random.randn(N), index=rng)
        with self.assertRaises(ValueError):
            s.asof(s.index[0], subset='foo')
