# -*- coding: utf-8 -*-

from datetime import datetime, timedelta, time

import numpy as np

from pandas import (DatetimeIndex, Float64Index, Index, Int64Index,
                    NaT, Period, PeriodIndex, Series, Timedelta,
                    TimedeltaIndex, date_range, period_range,
                    timedelta_range, notnull)

import pandas.util.testing as tm

import pandas as pd
from pandas.tslib import Timestamp, OutOfBoundsDatetime

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
        result = Index([Timestamp('2011-01-01'),
                        Timestamp('2011-01-02')], name='idx')
        exp = DatetimeIndex([Timestamp('2011-01-01'),
                             Timestamp('2011-01-02')], name='idx')
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
        exp = DatetimeIndex([Timestamp('2011-01-01 10:00'),
                             Timestamp('2011-01-02 10:00')],
                            tz='Asia/Tokyo', name='idx')
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

        with tm.assertRaisesRegexp(TypeError, 'data is already tz-aware'):
            DatetimeIndex([Timestamp('2011-01-01 10:00'),
                           Timestamp('2011-01-02 10:00', tz='US/Eastern')],
                          tz='Asia/Tokyo', name='idx')

        with tm.assertRaises(ValueError):
            DatetimeIndex([Timestamp('2011-01-01 10:00', tz='Asia/Tokyo'),
                           Timestamp('2011-01-02 10:00', tz='US/Eastern')],
                          tz='US/Eastern', name='idx')

        with tm.assertRaisesRegexp(TypeError, 'data is already tz-aware'):
            # passing tz should results in DatetimeIndex, then mismatch raises
            # TypeError
            Index([pd.NaT, Timestamp('2011-01-01 10:00'),
                   pd.NaT, Timestamp('2011-01-02 10:00', tz='US/Eastern')],
                  tz='Asia/Tokyo', name='idx')

    def test_construction_base_constructor(self):
        arr = [pd.Timestamp('2011-01-01'), pd.NaT, pd.Timestamp('2011-01-03')]
        tm.assert_index_equal(pd.Index(arr), pd.DatetimeIndex(arr))
        tm.assert_index_equal(pd.Index(np.array(arr)),
                              pd.DatetimeIndex(np.array(arr)))

        arr = [np.nan, pd.NaT, pd.Timestamp('2011-01-03')]
        tm.assert_index_equal(pd.Index(arr), pd.DatetimeIndex(arr))
        tm.assert_index_equal(pd.Index(np.array(arr)),
                              pd.DatetimeIndex(np.array(arr)))

    def test_construction_outofbounds(self):
        # GH 13663
        dates = [datetime(3000, 1, 1), datetime(4000, 1, 1),
                 datetime(5000, 1, 1), datetime(6000, 1, 1)]
        exp = Index(dates, dtype=object)
        # coerces to object
        tm.assert_index_equal(Index(dates), exp)

        with tm.assertRaises(OutOfBoundsDatetime):
            # can't create DatetimeIndex
            DatetimeIndex(dates)

    def test_astype(self):
        # GH 13149, GH 13209
        idx = DatetimeIndex(['2016-05-16', 'NaT', NaT, np.NaN])

        result = idx.astype(object)
        expected = Index([Timestamp('2016-05-16')] + [NaT] * 3, dtype=object)
        tm.assert_index_equal(result, expected)

        result = idx.astype(int)
        expected = Int64Index([1463356800000000000] +
                              [-9223372036854775808] * 3, dtype=np.int64)
        tm.assert_index_equal(result, expected)

        rng = date_range('1/1/2000', periods=10)
        result = rng.astype('i8')
        self.assert_index_equal(result, Index(rng.asi8))
        self.assert_numpy_array_equal(result.values, rng.asi8)

    def test_astype_with_tz(self):

        # with tz
        rng = date_range('1/1/2000', periods=10, tz='US/Eastern')
        result = rng.astype('datetime64[ns]')
        expected = (date_range('1/1/2000', periods=10,
                               tz='US/Eastern')
                    .tz_convert('UTC').tz_localize(None))
        tm.assert_index_equal(result, expected)

        # BUG#10442 : testing astype(str) is correct for Series/DatetimeIndex
        result = pd.Series(pd.date_range('2012-01-01', periods=3)).astype(str)
        expected = pd.Series(
            ['2012-01-01', '2012-01-02', '2012-01-03'], dtype=object)
        tm.assert_series_equal(result, expected)

        result = Series(pd.date_range('2012-01-01', periods=3,
                                      tz='US/Eastern')).astype(str)
        expected = Series(['2012-01-01 00:00:00-05:00',
                           '2012-01-02 00:00:00-05:00',
                           '2012-01-03 00:00:00-05:00'],
                          dtype=object)
        tm.assert_series_equal(result, expected)

    def test_astype_str_compat(self):
        # GH 13149, GH 13209
        # verify that we are returing NaT as a string (and not unicode)

        idx = DatetimeIndex(['2016-05-16', 'NaT', NaT, np.NaN])
        result = idx.astype(str)
        expected = Index(['2016-05-16', 'NaT', 'NaT', 'NaT'], dtype=object)
        tm.assert_index_equal(result, expected)

    def test_astype_str(self):
        # test astype string - #10442
        result = date_range('2012-01-01', periods=4,
                            name='test_name').astype(str)
        expected = Index(['2012-01-01', '2012-01-02', '2012-01-03',
                          '2012-01-04'], name='test_name', dtype=object)
        tm.assert_index_equal(result, expected)

        # test astype string with tz and name
        result = date_range('2012-01-01', periods=3, name='test_name',
                            tz='US/Eastern').astype(str)
        expected = Index(['2012-01-01 00:00:00-05:00',
                          '2012-01-02 00:00:00-05:00',
                          '2012-01-03 00:00:00-05:00'],
                         name='test_name', dtype=object)
        tm.assert_index_equal(result, expected)

        # test astype string with freqH and name
        result = date_range('1/1/2011', periods=3, freq='H',
                            name='test_name').astype(str)
        expected = Index(['2011-01-01 00:00:00', '2011-01-01 01:00:00',
                          '2011-01-01 02:00:00'],
                         name='test_name', dtype=object)
        tm.assert_index_equal(result, expected)

        # test astype string with freqH and timezone
        result = date_range('3/6/2012 00:00', periods=2, freq='H',
                            tz='Europe/London', name='test_name').astype(str)
        expected = Index(['2012-03-06 00:00:00+00:00',
                          '2012-03-06 01:00:00+00:00'],
                         dtype=object, name='test_name')
        tm.assert_index_equal(result, expected)

    def test_astype_datetime64(self):
        # GH 13149, GH 13209
        idx = DatetimeIndex(['2016-05-16', 'NaT', NaT, np.NaN])

        result = idx.astype('datetime64[ns]')
        tm.assert_index_equal(result, idx)
        self.assertFalse(result is idx)

        result = idx.astype('datetime64[ns]', copy=False)
        tm.assert_index_equal(result, idx)
        self.assertTrue(result is idx)

        idx_tz = DatetimeIndex(['2016-05-16', 'NaT', NaT, np.NaN], tz='EST')
        result = idx_tz.astype('datetime64[ns]')
        expected = DatetimeIndex(['2016-05-16 05:00:00', 'NaT', 'NaT', 'NaT'],
                                 dtype='datetime64[ns]')
        tm.assert_index_equal(result, expected)

    def test_astype_raises(self):
        # GH 13149, GH 13209
        idx = DatetimeIndex(['2016-05-16', 'NaT', NaT, np.NaN])

        self.assertRaises(ValueError, idx.astype, float)
        self.assertRaises(ValueError, idx.astype, 'timedelta64')
        self.assertRaises(ValueError, idx.astype, 'timedelta64[ns]')
        self.assertRaises(ValueError, idx.astype, 'datetime64')
        self.assertRaises(ValueError, idx.astype, 'datetime64[D]')

    def test_where_other(self):

        # other is ndarray or Index
        i = pd.date_range('20130101', periods=3, tz='US/Eastern')

        for arr in [np.nan, pd.NaT]:
            result = i.where(notnull(i), other=np.nan)
            expected = i
            tm.assert_index_equal(result, expected)

        i2 = i.copy()
        i2 = Index([pd.NaT, pd.NaT] + i[2:].tolist())
        result = i.where(notnull(i2), i2)
        tm.assert_index_equal(result, i2)

        i2 = i.copy()
        i2 = Index([pd.NaT, pd.NaT] + i[2:].tolist())
        result = i.where(notnull(i2), i2.values)
        tm.assert_index_equal(result, i2)

    def test_where_tz(self):
        i = pd.date_range('20130101', periods=3, tz='US/Eastern')
        result = i.where(notnull(i))
        expected = i
        tm.assert_index_equal(result, expected)

        i2 = i.copy()
        i2 = Index([pd.NaT, pd.NaT] + i[2:].tolist())
        result = i.where(notnull(i2))
        expected = i2
        tm.assert_index_equal(result, expected)

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
        tm.assert_numpy_array_equal(idx.get_loc(time(12)),
                                    np.array([12]), check_dtype=False)
        tm.assert_numpy_array_equal(idx.get_loc(time(12, 30)),
                                    np.array([]), check_dtype=False)
        with tm.assertRaises(NotImplementedError):
            idx.get_loc(time(12, 30), method='pad')

    def test_get_indexer(self):
        idx = pd.date_range('2000-01-01', periods=3)
        exp = np.array([0, 1, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(idx.get_indexer(idx), exp)

        target = idx[0] + pd.to_timedelta(['-1 hour', '12 hours',
                                           '1 day 1 hour'])
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'pad'),
                                    np.array([-1, 0, 1], dtype=np.intp))
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'backfill'),
                                    np.array([0, 1, 2], dtype=np.intp))
        tm.assert_numpy_array_equal(idx.get_indexer(target, 'nearest'),
                                    np.array([0, 1, 1], dtype=np.intp))
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, 'nearest',
                            tolerance=pd.Timedelta('1 hour')),
            np.array([0, -1, 1], dtype=np.intp))
        with tm.assertRaises(ValueError):
            idx.get_indexer(idx[[0]], method='nearest', tolerance='foo')

    def test_roundtrip_pickle_with_tz(self):

        # GH 8367
        # round-trip of timezone
        index = date_range('20130101', periods=3, tz='US/Eastern', name='foo')
        unpickled = self.round_trip_pickle(index)
        self.assert_index_equal(index, unpickled)

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

            tm.assert_numpy_array_equal(ts.index.get_loc(key), i,
                                        check_dtype=False)
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

            idx = pd.DatetimeIndex(['2011-01-01 09:00', pd.NaT,
                                    '2011-01-01 11:00'], tz=tz)

            exp = pd.DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00',
                                    '2011-01-01 11:00'], tz=tz)
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

    def test_difference_of_union(self):
        # GH14323: Test taking the union of differences of an Index.
        # Difference of DatetimeIndex does not preserve frequency,
        # so a differencing operation should not retain the freq field of the
        # original index.
        i = pd.date_range("20160920", "20160925", freq="D")

        a = pd.date_range("20160921", "20160924", freq="D")
        expected = pd.DatetimeIndex(["20160920", "20160925"], freq=None)
        a_diff = i.difference(a)
        tm.assert_index_equal(a_diff, expected)
        tm.assert_attr_equal('freq', a_diff, expected)

        b = pd.date_range("20160922", "20160925", freq="D")
        b_diff = i.difference(b)
        expected = pd.DatetimeIndex(["20160920", "20160921"], freq=None)
        tm.assert_index_equal(b_diff, expected)
        tm.assert_attr_equal('freq', b_diff, expected)

        union_of_diff = a_diff.union(b_diff)
        expected = pd.DatetimeIndex(["20160920", "20160921", "20160925"],
                                    freq=None)
        tm.assert_index_equal(union_of_diff, expected)
        tm.assert_attr_equal('freq', union_of_diff, expected)


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
        self.assert_frame_equal(df, df.ix[idx])
        self.assert_frame_equal(df, df.ix[list(idx)])
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

    def test_difference_of_union(self):
        # GH14323: Test taking the union of differences of an Index.
        # Difference of Period MUST preserve frequency, but the ability
        # to union results must be preserved
        i = pd.period_range("20160920", "20160925", freq="D")

        a = pd.period_range("20160921", "20160924", freq="D")
        expected = pd.PeriodIndex(["20160920", "20160925"], freq='D')
        a_diff = i.difference(a)
        tm.assert_index_equal(a_diff, expected)
        tm.assert_attr_equal('freq', a_diff, expected)

        b = pd.period_range("20160922", "20160925", freq="D")
        b_diff = i.difference(b)
        expected = pd.PeriodIndex(["20160920", "20160921"], freq='D')
        tm.assert_index_equal(b_diff, expected)
        tm.assert_attr_equal('freq', b_diff, expected)

        union_of_diff = a_diff.union(b_diff)
        expected = pd.PeriodIndex(["20160920", "20160921", "20160925"],
                                  freq='D')
        tm.assert_index_equal(union_of_diff, expected)
        tm.assert_attr_equal('freq', union_of_diff, expected)


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

    def test_difference_of_union(self):
        # GH14323: Test taking the union of differences of an Index.
        # Difference of TimedeltaIndex does not preserve frequency,
        # so a differencing operation should not retain the freq field of the
        # original index.
        i = pd.timedelta_range("0 days", "5 days", freq="D")

        a = pd.timedelta_range("1 days", "4 days", freq="D")
        expected = pd.TimedeltaIndex(["0 days", "5 days"], freq=None)
        a_diff = i.difference(a)
        tm.assert_index_equal(a_diff, expected)
        tm.assert_attr_equal('freq', a_diff, expected)

        b = pd.timedelta_range("2 days", "5 days", freq="D")
        b_diff = i.difference(b)
        expected = pd.TimedeltaIndex(["0 days", "1 days"], freq=None)
        tm.assert_index_equal(b_diff, expected)
        tm.assert_attr_equal('freq', b_diff, expected)

        union_of_difference = a_diff.union(b_diff)
        expected = pd.TimedeltaIndex(["0 days", "1 days", "5 days"],
                                     freq=None)
        tm.assert_index_equal(union_of_difference, expected)
        tm.assert_attr_equal('freq', union_of_difference, expected)
