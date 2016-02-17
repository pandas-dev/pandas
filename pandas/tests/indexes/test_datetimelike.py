# -*- coding: utf-8 -*-

from datetime import timedelta, time

import numpy as np

from pandas import (date_range, period_range,
                    Series, Index, DatetimeIndex,
                    TimedeltaIndex, PeriodIndex)

import pandas.util.testing as tm

import pandas as pd
from pandas.lib import Timestamp

from .common import Base


class DatetimeLike(Base):

    def test_shift_identity(self):

        idx = self.create_index()
        self.assert_index_equal(idx, idx.shift(0))

    def test_str(self):

        # test the string repr
        idx = self.create_index()
        idx.name = 'foo'
        self.assertFalse("length=%s" % len(idx) in str(idx))
        self.assertTrue("'foo'" in str(idx))
        self.assertTrue(idx.__class__.__name__ in str(idx))

        if hasattr(idx, 'tz'):
            if idx.tz is not None:
                self.assertTrue(idx.tz in str(idx))
        if hasattr(idx, 'freq'):
            self.assertTrue("freq='%s'" % idx.freqstr in str(idx))

    def test_view(self):
        super(DatetimeLike, self).test_view()

        i = self.create_index()

        i_view = i.view('i8')
        result = self._holder(i)
        tm.assert_index_equal(result, i)

        i_view = i.view(self._holder)
        result = self._holder(i)
        tm.assert_index_equal(result, i_view)


class TestDatetimeIndex(DatetimeLike, tm.TestCase):
    _holder = DatetimeIndex
    _multiprocess_can_split_ = True

    def setUp(self):
        self.indices = dict(index=tm.makeDateIndex(10))
        self.setup_indices()

    def create_index(self):
        return date_range('20130101', periods=5)

    def test_shift(self):

        # test shift for datetimeIndex and non datetimeIndex
        # GH8083

        drange = self.create_index()
        result = drange.shift(1)
        expected = DatetimeIndex(['2013-01-02', '2013-01-03', '2013-01-04',
                                  '2013-01-05',
                                  '2013-01-06'], freq='D')
        self.assert_index_equal(result, expected)

        result = drange.shift(-1)
        expected = DatetimeIndex(['2012-12-31', '2013-01-01', '2013-01-02',
                                  '2013-01-03', '2013-01-04'],
                                 freq='D')
        self.assert_index_equal(result, expected)

        result = drange.shift(3, freq='2D')
        expected = DatetimeIndex(['2013-01-07', '2013-01-08', '2013-01-09',
                                  '2013-01-10',
                                  '2013-01-11'], freq='D')
        self.assert_index_equal(result, expected)

    def test_construction_with_alt(self):

        i = pd.date_range('20130101', periods=5, freq='H', tz='US/Eastern')
        i2 = DatetimeIndex(i, dtype=i.dtype)
        self.assert_index_equal(i, i2)

        i2 = DatetimeIndex(i.tz_localize(None).asi8, tz=i.dtype.tz)
        self.assert_index_equal(i, i2)

        i2 = DatetimeIndex(i.tz_localize(None).asi8, dtype=i.dtype)
        self.assert_index_equal(i, i2)

        i2 = DatetimeIndex(
            i.tz_localize(None).asi8, dtype=i.dtype, tz=i.dtype.tz)
        self.assert_index_equal(i, i2)

        # localize into the provided tz
        i2 = DatetimeIndex(i.tz_localize(None).asi8, tz='UTC')
        expected = i.tz_localize(None).tz_localize('UTC')
        self.assert_index_equal(i2, expected)

        # incompat tz/dtype
        self.assertRaises(ValueError, lambda: DatetimeIndex(
            i.tz_localize(None).asi8, dtype=i.dtype, tz='US/Pacific'))

    def test_pickle_compat_construction(self):
        pass

    def test_construction_index_with_mixed_timezones(self):
        # GH 11488
        # no tz results in DatetimeIndex
        result = Index(
            [Timestamp('2011-01-01'), Timestamp('2011-01-02')], name='idx')
        exp = DatetimeIndex(
            [Timestamp('2011-01-01'), Timestamp('2011-01-02')], name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))
        self.assertIsNone(result.tz)

        # same tz results in DatetimeIndex
        result = Index([Timestamp('2011-01-01 10:00', tz='Asia/Tokyo'),
                        Timestamp('2011-01-02 10:00', tz='Asia/Tokyo')],
                       name='idx')
        exp = DatetimeIndex(
            [Timestamp('2011-01-01 10:00'), Timestamp('2011-01-02 10:00')
             ], tz='Asia/Tokyo', name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))
        self.assertIsNotNone(result.tz)
        self.assertEqual(result.tz, exp.tz)

        # same tz results in DatetimeIndex (DST)
        result = Index([Timestamp('2011-01-01 10:00', tz='US/Eastern'),
                        Timestamp('2011-08-01 10:00', tz='US/Eastern')],
                       name='idx')
        exp = DatetimeIndex([Timestamp('2011-01-01 10:00'),
                             Timestamp('2011-08-01 10:00')],
                            tz='US/Eastern', name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))
        self.assertIsNotNone(result.tz)
        self.assertEqual(result.tz, exp.tz)

        # different tz results in Index(dtype=object)
        result = Index([Timestamp('2011-01-01 10:00'),
                        Timestamp('2011-01-02 10:00', tz='US/Eastern')],
                       name='idx')
        exp = Index([Timestamp('2011-01-01 10:00'),
                     Timestamp('2011-01-02 10:00', tz='US/Eastern')],
                    dtype='object', name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertFalse(isinstance(result, DatetimeIndex))

        result = Index([Timestamp('2011-01-01 10:00', tz='Asia/Tokyo'),
                        Timestamp('2011-01-02 10:00', tz='US/Eastern')],
                       name='idx')
        exp = Index([Timestamp('2011-01-01 10:00', tz='Asia/Tokyo'),
                     Timestamp('2011-01-02 10:00', tz='US/Eastern')],
                    dtype='object', name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertFalse(isinstance(result, DatetimeIndex))

        # passing tz results in DatetimeIndex
        result = Index([Timestamp('2011-01-01 10:00'),
                        Timestamp('2011-01-02 10:00', tz='US/Eastern')],
                       tz='Asia/Tokyo', name='idx')
        exp = DatetimeIndex([Timestamp('2011-01-01 19:00'),
                             Timestamp('2011-01-03 00:00')],
                            tz='Asia/Tokyo', name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))

        # length = 1
        result = Index([Timestamp('2011-01-01')], name='idx')
        exp = DatetimeIndex([Timestamp('2011-01-01')], name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))
        self.assertIsNone(result.tz)

        # length = 1 with tz
        result = Index(
            [Timestamp('2011-01-01 10:00', tz='Asia/Tokyo')], name='idx')
        exp = DatetimeIndex([Timestamp('2011-01-01 10:00')], tz='Asia/Tokyo',
                            name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))
        self.assertIsNotNone(result.tz)
        self.assertEqual(result.tz, exp.tz)

    def test_construction_index_with_mixed_timezones_with_NaT(self):
        # GH 11488
        result = Index([pd.NaT, Timestamp('2011-01-01'),
                        pd.NaT, Timestamp('2011-01-02')], name='idx')
        exp = DatetimeIndex([pd.NaT, Timestamp('2011-01-01'),
                             pd.NaT, Timestamp('2011-01-02')], name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))
        self.assertIsNone(result.tz)

        # same tz results in DatetimeIndex
        result = Index([pd.NaT, Timestamp('2011-01-01 10:00', tz='Asia/Tokyo'),
                        pd.NaT, Timestamp('2011-01-02 10:00',
                                          tz='Asia/Tokyo')],
                       name='idx')
        exp = DatetimeIndex([pd.NaT, Timestamp('2011-01-01 10:00'),
                             pd.NaT, Timestamp('2011-01-02 10:00')],
                            tz='Asia/Tokyo', name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))
        self.assertIsNotNone(result.tz)
        self.assertEqual(result.tz, exp.tz)

        # same tz results in DatetimeIndex (DST)
        result = Index([Timestamp('2011-01-01 10:00', tz='US/Eastern'),
                        pd.NaT,
                        Timestamp('2011-08-01 10:00', tz='US/Eastern')],
                       name='idx')
        exp = DatetimeIndex([Timestamp('2011-01-01 10:00'), pd.NaT,
                             Timestamp('2011-08-01 10:00')],
                            tz='US/Eastern', name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))
        self.assertIsNotNone(result.tz)
        self.assertEqual(result.tz, exp.tz)

        # different tz results in Index(dtype=object)
        result = Index([pd.NaT, Timestamp('2011-01-01 10:00'),
                        pd.NaT, Timestamp('2011-01-02 10:00',
                                          tz='US/Eastern')],
                       name='idx')
        exp = Index([pd.NaT, Timestamp('2011-01-01 10:00'),
                     pd.NaT, Timestamp('2011-01-02 10:00', tz='US/Eastern')],
                    dtype='object', name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertFalse(isinstance(result, DatetimeIndex))

        result = Index([pd.NaT, Timestamp('2011-01-01 10:00', tz='Asia/Tokyo'),
                        pd.NaT, Timestamp('2011-01-02 10:00',
                                          tz='US/Eastern')], name='idx')
        exp = Index([pd.NaT, Timestamp('2011-01-01 10:00', tz='Asia/Tokyo'),
                     pd.NaT, Timestamp('2011-01-02 10:00', tz='US/Eastern')],
                    dtype='object', name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertFalse(isinstance(result, DatetimeIndex))

        # passing tz results in DatetimeIndex
        result = Index([pd.NaT, Timestamp('2011-01-01 10:00'),
                        pd.NaT, Timestamp('2011-01-02 10:00',
                                          tz='US/Eastern')],
                       tz='Asia/Tokyo', name='idx')
        exp = DatetimeIndex([pd.NaT, Timestamp('2011-01-01 19:00'),
                             pd.NaT, Timestamp('2011-01-03 00:00')],
                            tz='Asia/Tokyo', name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))

        # all NaT
        result = Index([pd.NaT, pd.NaT], name='idx')
        exp = DatetimeIndex([pd.NaT, pd.NaT], name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))
        self.assertIsNone(result.tz)

        # all NaT with tz
        result = Index([pd.NaT, pd.NaT], tz='Asia/Tokyo', name='idx')
        exp = DatetimeIndex([pd.NaT, pd.NaT], tz='Asia/Tokyo', name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))
        self.assertIsNotNone(result.tz)
        self.assertEqual(result.tz, exp.tz)

    def test_construction_dti_with_mixed_timezones(self):
        # GH 11488 (not changed, added explicit tests)

        # no tz results in DatetimeIndex
        result = DatetimeIndex(
            [Timestamp('2011-01-01'), Timestamp('2011-01-02')], name='idx')
        exp = DatetimeIndex(
            [Timestamp('2011-01-01'), Timestamp('2011-01-02')], name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))

        # same tz results in DatetimeIndex
        result = DatetimeIndex([Timestamp('2011-01-01 10:00', tz='Asia/Tokyo'),
                                Timestamp('2011-01-02 10:00',
                                          tz='Asia/Tokyo')],
                               name='idx')
        exp = DatetimeIndex(
            [Timestamp('2011-01-01 10:00'), Timestamp('2011-01-02 10:00')
             ], tz='Asia/Tokyo', name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))

        # same tz results in DatetimeIndex (DST)
        result = DatetimeIndex([Timestamp('2011-01-01 10:00', tz='US/Eastern'),
                                Timestamp('2011-08-01 10:00',
                                          tz='US/Eastern')],
                               name='idx')
        exp = DatetimeIndex([Timestamp('2011-01-01 10:00'),
                             Timestamp('2011-08-01 10:00')],
                            tz='US/Eastern', name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))

        # different tz coerces tz-naive to tz-awareIndex(dtype=object)
        result = DatetimeIndex([Timestamp('2011-01-01 10:00'),
                                Timestamp('2011-01-02 10:00',
                                          tz='US/Eastern')], name='idx')
        exp = DatetimeIndex([Timestamp('2011-01-01 05:00'),
                             Timestamp('2011-01-02 10:00')],
                            tz='US/Eastern', name='idx')
        self.assert_index_equal(result, exp, exact=True)
        self.assertTrue(isinstance(result, DatetimeIndex))

        # tz mismatch affecting to tz-aware raises TypeError/ValueError
        with tm.assertRaises(ValueError):
            DatetimeIndex([Timestamp('2011-01-01 10:00', tz='Asia/Tokyo'),
                           Timestamp('2011-01-02 10:00', tz='US/Eastern')],
                          name='idx')

        with tm.assertRaises(TypeError):
            DatetimeIndex([Timestamp('2011-01-01 10:00'),
                           Timestamp('2011-01-02 10:00', tz='US/Eastern')],
                          tz='Asia/Tokyo', name='idx')

        with tm.assertRaises(ValueError):
            DatetimeIndex([Timestamp('2011-01-01 10:00', tz='Asia/Tokyo'),
                           Timestamp('2011-01-02 10:00', tz='US/Eastern')],
                          tz='US/Eastern', name='idx')

    def test_get_loc(self):
        idx = pd.date_range('2000-01-01', periods=3)

        for method in [None, 'pad', 'backfill', 'nearest']:
            self.assertEqual(idx.get_loc(idx[1], method), 1)
            self.assertEqual(idx.get_loc(idx[1].to_pydatetime(), method), 1)
            self.assertEqual(idx.get_loc(str(idx[1]), method), 1)
            if method is not None:
                self.assertEqual(idx.get_loc(idx[1], method,
                                             tolerance=pd.Timedelta('0 days')),
                                 1)

        self.assertEqual(idx.get_loc('2000-01-01', method='nearest'), 0)
        self.assertEqual(idx.get_loc('2000-01-01T12', method='nearest'), 1)

        self.assertEqual(idx.get_loc('2000-01-01T12', method='nearest',
                                     tolerance='1 day'), 1)
        self.assertEqual(idx.get_loc('2000-01-01T12', method='nearest',
                                     tolerance=pd.Timedelta('1D')), 1)
        self.assertEqual(idx.get_loc('2000-01-01T12', method='nearest',
                                     tolerance=np.timedelta64(1, 'D')), 1)
        self.assertEqual(idx.get_loc('2000-01-01T12', method='nearest',
                                     tolerance=timedelta(1)), 1)
        with tm.assertRaisesRegexp(ValueError, 'must be convertible'):
            idx.get_loc('2000-01-01T12', method='nearest', tolerance='foo')
        with tm.assertRaises(KeyError):
            idx.get_loc('2000-01-01T03', method='nearest', tolerance='2 hours')

        self.assertEqual(idx.get_loc('2000', method='nearest'), slice(0, 3))
        self.assertEqual(idx.get_loc('2000-01', method='nearest'), slice(0, 3))

        self.assertEqual(idx.get_loc('1999', method='nearest'), 0)
        self.assertEqual(idx.get_loc('2001', method='nearest'), 2)

        with tm.assertRaises(KeyError):
            idx.get_loc('1999', method='pad')
        with tm.assertRaises(KeyError):
            idx.get_loc('2001', method='backfill')

        with tm.assertRaises(KeyError):
            idx.get_loc('foobar')
        with tm.assertRaises(TypeError):
            idx.get_loc(slice(2))

        idx = pd.to_datetime(['2000-01-01', '2000-01-04'])
        self.assertEqual(idx.get_loc('2000-01-02', method='nearest'), 0)
        self.assertEqual(idx.get_loc('2000-01-03', method='nearest'), 1)
        self.assertEqual(idx.get_loc('2000-01', method='nearest'), slice(0, 2))

        # time indexing
        idx = pd.date_range('2000-01-01', periods=24, freq='H')
        tm.assert_numpy_array_equal(idx.get_loc(time(12)), [12])
        tm.assert_numpy_array_equal(idx.get_loc(time(12, 30)), [])
        with tm.assertRaises(NotImplementedError):
            idx.get_loc(time(12, 30), method='pad')

    def test_get_indexer(self):
        idx = pd.date_range('2000-01-01', periods=3)
        tm.assert_numpy_array_equal(idx.get_indexer(idx), [0, 1, 2])

        target = idx[0] + pd.to_timedelta(['-1 hour', '12 hours',
                                           '1 day 1 hour'])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'pad'), [-1, 0, 1])
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'backfill'), [0, 1, 2])
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'nearest'), [0, 1, 1])
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'nearest',
                            tolerance=pd.Timedelta('1 hour')),
            [0, -1, 1])
        with tm.assertRaises(ValueError):
            idx.get_indexer(idx[[0]], method='nearest', tolerance='foo')

    def test_roundtrip_pickle_with_tz(self):

        # GH 8367
        # round-trip of timezone
        index = date_range('20130101', periods=3, tz='US/Eastern', name='foo')
        unpickled = self.round_trip_pickle(index)
        self.assertTrue(index.equals(unpickled))

    def test_reindex_preserves_tz_if_target_is_empty_list_or_array(self):
        # GH7774
        index = date_range('20130101', periods=3, tz='US/Eastern')
        self.assertEqual(str(index.reindex([])[0].tz), 'US/Eastern')
        self.assertEqual(str(index.reindex(np.array([]))[0].tz), 'US/Eastern')

    def test_time_loc(self):  # GH8667
        from datetime import time
        from pandas.index import _SIZE_CUTOFF

        ns = _SIZE_CUTOFF + np.array([-100, 100], dtype=np.int64)
        key = time(15, 11, 30)
        start = key.hour * 3600 + key.minute * 60 + key.second
        step = 24 * 3600

        for n in ns:
            idx = pd.date_range('2014-11-26', periods=n, freq='S')
            ts = pd.Series(np.random.randn(n), index=idx)
            i = np.arange(start, n, step)

            tm.assert_numpy_array_equal(ts.index.get_loc(key), i)
            tm.assert_series_equal(ts[key], ts.iloc[i])

            left, right = ts.copy(), ts.copy()
            left[key] *= -10
            right.iloc[i] *= -10
            tm.assert_series_equal(left, right)

    def test_time_overflow_for_32bit_machines(self):
        # GH8943.  On some machines NumPy defaults to np.int32 (for example,
        # 32-bit Linux machines).  In the function _generate_regular_range
        # found in tseries/index.py, `periods` gets multiplied by `strides`
        # (which has value 1e9) and since the max value for np.int32 is ~2e9,
        # and since those machines won't promote np.int32 to np.int64, we get
        # overflow.
        periods = np.int_(1000)

        idx1 = pd.date_range(start='2000', periods=periods, freq='S')
        self.assertEqual(len(idx1), periods)

        idx2 = pd.date_range(end='2000', periods=periods, freq='S')
        self.assertEqual(len(idx2), periods)

    def test_intersection(self):
        first = self.index
        second = self.index[5:]
        intersect = first.intersection(second)
        self.assertTrue(tm.equalContents(intersect, second))

        # GH 10149
        cases = [klass(second.values) for klass in [np.array, Series, list]]
        for case in cases:
            result = first.intersection(case)
            self.assertTrue(tm.equalContents(result, second))

        third = Index(['a', 'b', 'c'])
        result = first.intersection(third)
        expected = pd.Index([], dtype=object)
        self.assert_index_equal(result, expected)

    def test_union(self):
        first = self.index[:5]
        second = self.index[5:]
        everything = self.index
        union = first.union(second)
        self.assertTrue(tm.equalContents(union, everything))

        # GH 10149
        cases = [klass(second.values) for klass in [np.array, Series, list]]
        for case in cases:
            result = first.union(case)
            self.assertTrue(tm.equalContents(result, everything))

    def test_nat(self):
        self.assertIs(DatetimeIndex([np.nan])[0], pd.NaT)

    def test_ufunc_coercions(self):
        idx = date_range('2011-01-01', periods=3, freq='2D', name='x')

        delta = np.timedelta64(1, 'D')
        for result in [idx + delta, np.add(idx, delta)]:
            tm.assertIsInstance(result, DatetimeIndex)
            exp = date_range('2011-01-02', periods=3, freq='2D', name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, '2D')

        for result in [idx - delta, np.subtract(idx, delta)]:
            tm.assertIsInstance(result, DatetimeIndex)
            exp = date_range('2010-12-31', periods=3, freq='2D', name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, '2D')

        delta = np.array([np.timedelta64(1, 'D'), np.timedelta64(2, 'D'),
                          np.timedelta64(3, 'D')])
        for result in [idx + delta, np.add(idx, delta)]:
            tm.assertIsInstance(result, DatetimeIndex)
            exp = DatetimeIndex(['2011-01-02', '2011-01-05', '2011-01-08'],
                                freq='3D', name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, '3D')

        for result in [idx - delta, np.subtract(idx, delta)]:
            tm.assertIsInstance(result, DatetimeIndex)
            exp = DatetimeIndex(['2010-12-31', '2011-01-01', '2011-01-02'],
                                freq='D', name='x')
            tm.assert_index_equal(result, exp)
            self.assertEqual(result.freq, 'D')

    def test_fillna_datetime64(self):
        # GH 11343
        for tz in ['US/Eastern', 'Asia/Tokyo']:
            idx = pd.DatetimeIndex(['2011-01-01 09:00', pd.NaT,
                                    '2011-01-01 11:00'])

            exp = pd.DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00',
                                    '2011-01-01 11:00'])
            self.assert_index_equal(
                idx.fillna(pd.Timestamp('2011-01-01 10:00')), exp)

            # tz mismatch
            exp = pd.Index([pd.Timestamp('2011-01-01 09:00'),
                            pd.Timestamp('2011-01-01 10:00', tz=tz),
                            pd.Timestamp('2011-01-01 11:00')], dtype=object)
            self.assert_index_equal(
                idx.fillna(pd.Timestamp('2011-01-01 10:00', tz=tz)), exp)

            # object
            exp = pd.Index([pd.Timestamp('2011-01-01 09:00'), 'x',
                            pd.Timestamp('2011-01-01 11:00')], dtype=object)
            self.assert_index_equal(idx.fillna('x'), exp)

            idx = pd.DatetimeIndex(
                ['2011-01-01 09:00', pd.NaT, '2011-01-01 11:00'], tz=tz)

            exp = pd.DatetimeIndex(
                ['2011-01-01 09:00', '2011-01-01 10:00', '2011-01-01 11:00'
                 ], tz=tz)
            self.assert_index_equal(
                idx.fillna(pd.Timestamp('2011-01-01 10:00', tz=tz)), exp)

            exp = pd.Index([pd.Timestamp('2011-01-01 09:00', tz=tz),
                            pd.Timestamp('2011-01-01 10:00'),
                            pd.Timestamp('2011-01-01 11:00', tz=tz)],
                           dtype=object)
            self.assert_index_equal(
                idx.fillna(pd.Timestamp('2011-01-01 10:00')), exp)

            # object
            exp = pd.Index([pd.Timestamp('2011-01-01 09:00', tz=tz),
                            'x',
                            pd.Timestamp('2011-01-01 11:00', tz=tz)],
                           dtype=object)
            self.assert_index_equal(idx.fillna('x'), exp)


class TestPeriodIndex(DatetimeLike, tm.TestCase):
    _holder = PeriodIndex
    _multiprocess_can_split_ = True

    def setUp(self):
        self.indices = dict(index=tm.makePeriodIndex(10))
        self.setup_indices()

    def create_index(self):
        return period_range('20130101', periods=5, freq='D')

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

    def test_get_indexer(self):
        idx = pd.period_range('2000-01-01', periods=3).asfreq('H', how='start')
        tm.assert_numpy_array_equal(idx.get_indexer(idx), [0, 1, 2])

        target = pd.PeriodIndex(['1999-12-31T23', '2000-01-01T12',
                                 '2000-01-02T01'], freq='H')
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'pad'), [-1, 0, 1])
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'backfill'), [0, 1, 2])
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'nearest'), [0, 1, 1])
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'nearest', tolerance='1 hour'),
            [0, -1, 1])

        msg = 'Input has different freq from PeriodIndex\\(freq=H\\)'
        with self.assertRaisesRegexp(ValueError, msg):
            idx.get_indexer(target, 'nearest', tolerance='1 minute')

        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'nearest', tolerance='1 day'), [0, 1, 1])

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
        self.assert_frame_equal(df, df.ix[idx])
        self.assert_frame_equal(df, df.ix[list(idx)])
        self.assert_frame_equal(df, df.loc[list(idx)])
        self.assert_frame_equal(df.iloc[0:5], df.loc[idx[0:5]])
        self.assert_frame_equal(df, df.loc[list(idx)])

    def test_fillna_period(self):
        # GH 11343
        idx = pd.PeriodIndex(
            ['2011-01-01 09:00', pd.NaT, '2011-01-01 11:00'], freq='H')

        exp = pd.PeriodIndex(
            ['2011-01-01 09:00', '2011-01-01 10:00', '2011-01-01 11:00'
             ], freq='H')
        self.assert_index_equal(
            idx.fillna(pd.Period('2011-01-01 10:00', freq='H')), exp)

        exp = pd.Index([pd.Period('2011-01-01 09:00', freq='H'), 'x',
                        pd.Period('2011-01-01 11:00', freq='H')], dtype=object)
        self.assert_index_equal(idx.fillna('x'), exp)

        with tm.assertRaisesRegexp(
                ValueError,
                'Input has different freq=D from PeriodIndex\\(freq=H\\)'):
            idx.fillna(pd.Period('2011-01-01', freq='D'))

    def test_no_millisecond_field(self):
        with self.assertRaises(AttributeError):
            DatetimeIndex.millisecond

        with self.assertRaises(AttributeError):
            DatetimeIndex([]).millisecond


class TestTimedeltaIndex(DatetimeLike, tm.TestCase):
    _holder = TimedeltaIndex
    _multiprocess_can_split_ = True

    def setUp(self):
        self.indices = dict(index=tm.makeTimedeltaIndex(10))
        self.setup_indices()

    def create_index(self):
        return pd.to_timedelta(range(5), unit='d') + pd.offsets.Hour(1)

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
        tm.assert_numpy_array_equal(idx.get_indexer(idx), [0, 1, 2])

        target = pd.to_timedelta(['-1 hour', '12 hours', '1 day 1 hour'])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'pad'), [-1, 0, 1])
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'backfill'), [0, 1, 2])
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'nearest'), [0, 1, 1])
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'nearest',
                            tolerance=pd.Timedelta('1 hour')),
            [0, -1, 1])

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
