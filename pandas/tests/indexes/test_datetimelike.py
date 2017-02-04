# -*- coding: utf-8 -*-

import numpy as np
from datetime import timedelta

import pandas as pd
from pandas.util import testing as tm
from pandas import (DatetimeIndex, Float64Index, Index, Int64Index,
                    NaT, Period, PeriodIndex, Series, Timedelta,
                    TimedeltaIndex, period_range,
                    timedelta_range, notnull)


from .datetimelike import DatetimeLike


class TestPeriodIndex(DatetimeLike, tm.TestCase):
    _holder = PeriodIndex
    _multiprocess_can_split_ = True

    def setUp(self):
        self.indices = dict(index=tm.makePeriodIndex(10))
        self.setup_indices()

    def create_index(self):
        return period_range('20130101', periods=5, freq='D')

    def test_construction_base_constructor(self):
        # GH 13664
        arr = [pd.Period('2011-01', freq='M'), pd.NaT,
               pd.Period('2011-03', freq='M')]
        tm.assert_index_equal(pd.Index(arr), pd.PeriodIndex(arr))
        tm.assert_index_equal(pd.Index(np.array(arr)),
                              pd.PeriodIndex(np.array(arr)))

        arr = [np.nan, pd.NaT, pd.Period('2011-03', freq='M')]
        tm.assert_index_equal(pd.Index(arr), pd.PeriodIndex(arr))
        tm.assert_index_equal(pd.Index(np.array(arr)),
                              pd.PeriodIndex(np.array(arr)))

        arr = [pd.Period('2011-01', freq='M'), pd.NaT,
               pd.Period('2011-03', freq='D')]
        tm.assert_index_equal(pd.Index(arr), pd.Index(arr, dtype=object))

        tm.assert_index_equal(pd.Index(np.array(arr)),
                              pd.Index(np.array(arr), dtype=object))

    def test_astype(self):
        # GH 13149, GH 13209
        idx = PeriodIndex(['2016-05-16', 'NaT', NaT, np.NaN], freq='D')

        result = idx.astype(object)
        expected = Index([Period('2016-05-16', freq='D')] +
                         [Period(NaT, freq='D')] * 3, dtype='object')
        tm.assert_index_equal(result, expected)

        result = idx.astype(int)
        expected = Int64Index([16937] + [-9223372036854775808] * 3,
                              dtype=np.int64)
        tm.assert_index_equal(result, expected)

        idx = period_range('1990', '2009', freq='A')
        result = idx.astype('i8')
        self.assert_index_equal(result, Index(idx.asi8))
        self.assert_numpy_array_equal(result.values, idx.asi8)

    def test_astype_raises(self):
        # GH 13149, GH 13209
        idx = PeriodIndex(['2016-05-16', 'NaT', NaT, np.NaN], freq='D')

        self.assertRaises(ValueError, idx.astype, str)
        self.assertRaises(ValueError, idx.astype, float)
        self.assertRaises(ValueError, idx.astype, 'timedelta64')
        self.assertRaises(ValueError, idx.astype, 'timedelta64[ns]')

    def test_shift(self):

        # test shift for PeriodIndex
        # GH8083
        drange = self.create_index()
        result = drange.shift(1)
        expected = PeriodIndex(['2013-01-02', '2013-01-03', '2013-01-04',
                                '2013-01-05', '2013-01-06'], freq='D')
        self.assert_index_equal(result, expected)

    def test_pickle_compat_construction(self):
        pass

    def test_get_loc(self):
        idx = pd.period_range('2000-01-01', periods=3)

        for method in [None, 'pad', 'backfill', 'nearest']:
            self.assertEqual(idx.get_loc(idx[1], method), 1)
            self.assertEqual(
                idx.get_loc(idx[1].asfreq('H', how='start'), method), 1)
            self.assertEqual(idx.get_loc(idx[1].to_timestamp(), method), 1)
            self.assertEqual(
                idx.get_loc(idx[1].to_timestamp().to_pydatetime(), method), 1)
            self.assertEqual(idx.get_loc(str(idx[1]), method), 1)

        idx = pd.period_range('2000-01-01', periods=5)[::2]
        self.assertEqual(idx.get_loc('2000-01-02T12', method='nearest',
                                     tolerance='1 day'), 1)
        self.assertEqual(idx.get_loc('2000-01-02T12', method='nearest',
                                     tolerance=pd.Timedelta('1D')), 1)
        self.assertEqual(idx.get_loc('2000-01-02T12', method='nearest',
                                     tolerance=np.timedelta64(1, 'D')), 1)
        self.assertEqual(idx.get_loc('2000-01-02T12', method='nearest',
                                     tolerance=timedelta(1)), 1)
        with tm.assertRaisesRegexp(ValueError, 'must be convertible'):
            idx.get_loc('2000-01-10', method='nearest', tolerance='foo')

        msg = 'Input has different freq from PeriodIndex\\(freq=D\\)'
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.get_loc('2000-01-10', method='nearest', tolerance='1 hour')
        with tm.assertRaises(KeyError):
            idx.get_loc('2000-01-10', method='nearest', tolerance='1 day')

    def test_where(self):
        i = self.create_index()
        result = i.where(notnull(i))
        expected = i
        tm.assert_index_equal(result, expected)

        i2 = i.copy()
        i2 = pd.PeriodIndex([pd.NaT, pd.NaT] + i[2:].tolist(),
                            freq='D')
        result = i.where(notnull(i2))
        expected = i2
        tm.assert_index_equal(result, expected)

    def test_where_other(self):

        i = self.create_index()
        for arr in [np.nan, pd.NaT]:
            result = i.where(notnull(i), other=np.nan)
            expected = i
            tm.assert_index_equal(result, expected)

        i2 = i.copy()
        i2 = pd.PeriodIndex([pd.NaT, pd.NaT] + i[2:].tolist(),
                            freq='D')
        result = i.where(notnull(i2), i2)
        tm.assert_index_equal(result, i2)

        i2 = i.copy()
        i2 = pd.PeriodIndex([pd.NaT, pd.NaT] + i[2:].tolist(),
                            freq='D')
        result = i.where(notnull(i2), i2.values)
        tm.assert_index_equal(result, i2)

    def test_get_indexer(self):
        idx = pd.period_range('2000-01-01', periods=3).asfreq('H', how='start')
        tm.assert_numpy_array_equal(idx.get_indexer(idx),
                                    np.array([0, 1, 2], dtype=np.intp))

        target = pd.PeriodIndex(['1999-12-31T23', '2000-01-01T12',
                                 '2000-01-02T01'], freq='H')
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'pad'),
                                    np.array([-1, 0, 1], dtype=np.intp))
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'backfill'),
                                    np.array([0, 1, 2], dtype=np.intp))
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'nearest'),
                                    np.array([0, 1, 1], dtype=np.intp))
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'nearest',
                                                    tolerance='1 hour'),
                                    np.array([0, -1, 1], dtype=np.intp))

        msg = 'Input has different freq from PeriodIndex\\(freq=H\\)'
        with self.assertRaisesRegexp(ValueError, msg):
            idx.get_indexer(target, 'nearest', tolerance='1 minute')

        tm.assert_numpy_array_equal(idx.get_indexer(target, 'nearest',
                                                    tolerance='1 day'),
                                    np.array([0, 1, 1], dtype=np.intp))

    def test_repeat(self):
        # GH10183
        idx = pd.period_range('2000-01-01', periods=3, freq='D')
        res = idx.repeat(3)
        exp = PeriodIndex(idx.values.repeat(3), freq='D')
        self.assert_index_equal(res, exp)
        self.assertEqual(res.freqstr, 'D')

    def test_period_index_indexer(self):
        # GH4125
        idx = pd.period_range('2002-01', '2003-12', freq='M')
        df = pd.DataFrame(pd.np.random.randn(24, 10), index=idx)
        self.assert_frame_equal(df, df.loc[idx])
        self.assert_frame_equal(df, df.loc[list(idx)])
        self.assert_frame_equal(df, df.loc[list(idx)])
        self.assert_frame_equal(df.iloc[0:5], df.loc[idx[0:5]])
        self.assert_frame_equal(df, df.loc[list(idx)])

    def test_fillna_period(self):
        # GH 11343
        idx = pd.PeriodIndex(['2011-01-01 09:00', pd.NaT,
                              '2011-01-01 11:00'], freq='H')

        exp = pd.PeriodIndex(['2011-01-01 09:00', '2011-01-01 10:00',
                              '2011-01-01 11:00'], freq='H')
        self.assert_index_equal(
            idx.fillna(pd.Period('2011-01-01 10:00', freq='H')), exp)

        exp = pd.Index([pd.Period('2011-01-01 09:00', freq='H'), 'x',
                        pd.Period('2011-01-01 11:00', freq='H')], dtype=object)
        self.assert_index_equal(idx.fillna('x'), exp)

        exp = pd.Index([pd.Period('2011-01-01 09:00', freq='H'),
                        pd.Period('2011-01-01', freq='D'),
                        pd.Period('2011-01-01 11:00', freq='H')], dtype=object)
        self.assert_index_equal(idx.fillna(pd.Period('2011-01-01', freq='D')),
                                exp)

    def test_no_millisecond_field(self):
        with self.assertRaises(AttributeError):
            DatetimeIndex.millisecond

        with self.assertRaises(AttributeError):
            DatetimeIndex([]).millisecond

    def test_difference_freq(self):
        # GH14323: difference of Period MUST preserve frequency
        # but the ability to union results must be preserved

        index = period_range("20160920", "20160925", freq="D")

        other = period_range("20160921", "20160924", freq="D")
        expected = PeriodIndex(["20160920", "20160925"], freq='D')
        idx_diff = index.difference(other)
        tm.assert_index_equal(idx_diff, expected)
        tm.assert_attr_equal('freq', idx_diff, expected)

        other = period_range("20160922", "20160925", freq="D")
        idx_diff = index.difference(other)
        expected = PeriodIndex(["20160920", "20160921"], freq='D')
        tm.assert_index_equal(idx_diff, expected)
        tm.assert_attr_equal('freq', idx_diff, expected)


class TestTimedeltaIndex(DatetimeLike, tm.TestCase):
    _holder = TimedeltaIndex
    _multiprocess_can_split_ = True

    def setUp(self):
        self.indices = dict(index=tm.makeTimedeltaIndex(10))
        self.setup_indices()

    def create_index(self):
        return pd.to_timedelta(range(5), unit='d') + pd.offsets.Hour(1)

    def test_construction_base_constructor(self):
        arr = [pd.Timedelta('1 days'), pd.NaT, pd.Timedelta('3 days')]
        tm.assert_index_equal(pd.Index(arr), pd.TimedeltaIndex(arr))
        tm.assert_index_equal(pd.Index(np.array(arr)),
                              pd.TimedeltaIndex(np.array(arr)))

        arr = [np.nan, pd.NaT, pd.Timedelta('1 days')]
        tm.assert_index_equal(pd.Index(arr), pd.TimedeltaIndex(arr))
        tm.assert_index_equal(pd.Index(np.array(arr)),
                              pd.TimedeltaIndex(np.array(arr)))

    def test_shift(self):
        # test shift for TimedeltaIndex
        # err8083

        drange = self.create_index()
        result = drange.shift(1)
        expected = TimedeltaIndex(['1 days 01:00:00', '2 days 01:00:00',
                                   '3 days 01:00:00',
                                   '4 days 01:00:00', '5 days 01:00:00'],
                                  freq='D')
        self.assert_index_equal(result, expected)

        result = drange.shift(3, freq='2D 1s')
        expected = TimedeltaIndex(['6 days 01:00:03', '7 days 01:00:03',
                                   '8 days 01:00:03', '9 days 01:00:03',
                                   '10 days 01:00:03'], freq='D')
        self.assert_index_equal(result, expected)

    def test_astype(self):
        # GH 13149, GH 13209
        idx = TimedeltaIndex([1e14, 'NaT', pd.NaT, np.NaN])

        result = idx.astype(object)
        expected = Index([Timedelta('1 days 03:46:40')] + [pd.NaT] * 3,
                         dtype=object)
        tm.assert_index_equal(result, expected)

        result = idx.astype(int)
        expected = Int64Index([100000000000000] + [-9223372036854775808] * 3,
                              dtype=np.int64)
        tm.assert_index_equal(result, expected)

        rng = timedelta_range('1 days', periods=10)

        result = rng.astype('i8')
        self.assert_index_equal(result, Index(rng.asi8))
        self.assert_numpy_array_equal(rng.asi8, result.values)

    def test_astype_timedelta64(self):
        # GH 13149, GH 13209
        idx = TimedeltaIndex([1e14, 'NaT', pd.NaT, np.NaN])

        result = idx.astype('timedelta64')
        expected = Float64Index([1e+14] + [np.NaN] * 3, dtype='float64')
        tm.assert_index_equal(result, expected)

        result = idx.astype('timedelta64[ns]')
        tm.assert_index_equal(result, idx)
        self.assertFalse(result is idx)

        result = idx.astype('timedelta64[ns]', copy=False)
        tm.assert_index_equal(result, idx)
        self.assertTrue(result is idx)

    def test_astype_raises(self):
        # GH 13149, GH 13209
        idx = TimedeltaIndex([1e14, 'NaT', pd.NaT, np.NaN])

        self.assertRaises(ValueError, idx.astype, float)
        self.assertRaises(ValueError, idx.astype, str)
        self.assertRaises(ValueError, idx.astype, 'datetime64')
        self.assertRaises(ValueError, idx.astype, 'datetime64[ns]')

    def test_get_loc(self):
        idx = pd.to_timedelta(['0 days', '1 days', '2 days'])

        for method in [None, 'pad', 'backfill', 'nearest']:
            self.assertEqual(idx.get_loc(idx[1], method), 1)
            self.assertEqual(idx.get_loc(idx[1].to_pytimedelta(), method), 1)
            self.assertEqual(idx.get_loc(str(idx[1]), method), 1)

        self.assertEqual(
            idx.get_loc(idx[1], 'pad', tolerance=pd.Timedelta(0)), 1)
        self.assertEqual(
            idx.get_loc(idx[1], 'pad', tolerance=np.timedelta64(0, 's')), 1)
        self.assertEqual(idx.get_loc(idx[1], 'pad', tolerance=timedelta(0)), 1)

        with tm.assertRaisesRegexp(ValueError, 'must be convertible'):
            idx.get_loc(idx[1], method='nearest', tolerance='foo')

        for method, loc in [('pad', 1), ('backfill', 2), ('nearest', 1)]:
            self.assertEqual(idx.get_loc('1 day 1 hour', method), loc)

    def test_get_indexer(self):
        idx = pd.to_timedelta(['0 days', '1 days', '2 days'])
        tm.assert_numpy_array_equal(idx.get_indexer(idx),
                                    np.array([0, 1, 2], dtype=np.intp))

        target = pd.to_timedelta(['-1 hour', '12 hours', '1 day 1 hour'])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'pad'),
                                    np.array([-1, 0, 1], dtype=np.intp))
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'backfill'),
                                    np.array([0, 1, 2], dtype=np.intp))
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'nearest'),
                                    np.array([0, 1, 1], dtype=np.intp))

        res = idx.get_indexer(target, 'nearest',
                              tolerance=pd.Timedelta('1 hour'))
        tm.assert_numpy_array_equal(res, np.array([0, -1, 1], dtype=np.intp))

    def test_numeric_compat(self):

        idx = self._holder(np.arange(5, dtype='int64'))
        didx = self._holder(np.arange(5, dtype='int64') ** 2)
        result = idx * 1
        tm.assert_index_equal(result, idx)

        result = 1 * idx
        tm.assert_index_equal(result, idx)

        result = idx / 1
        tm.assert_index_equal(result, idx)

        result = idx // 1
        tm.assert_index_equal(result, idx)

        result = idx * np.array(5, dtype='int64')
        tm.assert_index_equal(result,
                              self._holder(np.arange(5, dtype='int64') * 5))

        result = idx * np.arange(5, dtype='int64')
        tm.assert_index_equal(result, didx)

        result = idx * Series(np.arange(5, dtype='int64'))
        tm.assert_index_equal(result, didx)

        result = idx * Series(np.arange(5, dtype='float64') + 0.1)
        tm.assert_index_equal(result, self._holder(np.arange(
            5, dtype='float64') * (np.arange(5, dtype='float64') + 0.1)))

        # invalid
        self.assertRaises(TypeError, lambda: idx * idx)
        self.assertRaises(ValueError, lambda: idx * self._holder(np.arange(3)))
        self.assertRaises(ValueError, lambda: idx * np.array([1, 2]))

    def test_pickle_compat_construction(self):
        pass

    def test_ufunc_coercions(self):
        # normal ops are also tested in tseries/test_timedeltas.py
        idx = TimedeltaIndex(['2H', '4H', '6H', '8H', '10H'],
                             freq='2H', name='x')

        for result in [idx * 2, np.multiply(idx, 2)]:
            tm.assertIsInstance(result, TimedeltaIndex)
            exp = TimedeltaIndex(['4H', '8H', '12H', '16H', '20H'],
                                 freq='4H', name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, '4H')

        for result in [idx / 2, np.divide(idx, 2)]:
            tm.assertIsInstance(result, TimedeltaIndex)
            exp = TimedeltaIndex(['1H', '2H', '3H', '4H', '5H'],
                                 freq='H', name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, 'H')

        idx = TimedeltaIndex(['2H', '4H', '6H', '8H', '10H'],
                             freq='2H', name='x')
        for result in [-idx, np.negative(idx)]:
            tm.assertIsInstance(result, TimedeltaIndex)
            exp = TimedeltaIndex(['-2H', '-4H', '-6H', '-8H', '-10H'],
                                 freq='-2H', name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, '-2H')

        idx = TimedeltaIndex(['-2H', '-1H', '0H', '1H', '2H'],
                             freq='H', name='x')
        for result in [abs(idx), np.absolute(idx)]:
            tm.assertIsInstance(result, TimedeltaIndex)
            exp = TimedeltaIndex(['2H', '1H', '0H', '1H', '2H'],
                                 freq=None, name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, None)

    def test_fillna_timedelta(self):
        # GH 11343
        idx = pd.TimedeltaIndex(['1 day', pd.NaT, '3 day'])

        exp = pd.TimedeltaIndex(['1 day', '2 day', '3 day'])
        self.assert_index_equal(idx.fillna(pd.Timedelta('2 day')), exp)

        exp = pd.TimedeltaIndex(['1 day', '3 hour', '3 day'])
        idx.fillna(pd.Timedelta('3 hour'))

        exp = pd.Index(
            [pd.Timedelta('1 day'), 'x', pd.Timedelta('3 day')], dtype=object)
        self.assert_index_equal(idx.fillna('x'), exp)

    def test_difference_freq(self):
        # GH14323: Difference of TimedeltaIndex should not preserve frequency

        index = timedelta_range("0 days", "5 days", freq="D")

        other = timedelta_range("1 days", "4 days", freq="D")
        expected = TimedeltaIndex(["0 days", "5 days"], freq=None)
        idx_diff = index.difference(other)
        tm.assert_index_equal(idx_diff, expected)
        tm.assert_attr_equal('freq', idx_diff, expected)

        other = timedelta_range("2 days", "5 days", freq="D")
        idx_diff = index.difference(other)
        expected = TimedeltaIndex(["0 days", "1 days"], freq=None)
        tm.assert_index_equal(idx_diff, expected)
        tm.assert_attr_equal('freq', idx_diff, expected)
