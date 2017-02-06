import numpy as np
from datetime import timedelta

import pandas as pd
from pandas.util import testing as tm
from pandas import (PeriodIndex, period_range, notnull, DatetimeIndex, NaT,
                    Index, Period, Int64Index)

from ..datetimelike import DatetimeLike


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
