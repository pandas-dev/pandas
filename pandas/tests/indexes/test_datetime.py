import warnings
import numpy as np
from datetime import time, timedelta, date

import pandas as pd
import pandas.tslib as tslib
import pandas.util.testing as tm
import pandas.compat as compat
from pandas.compat import lrange
from pandas.core.common import PerformanceWarning
from pandas.compat.numpy import np_datetime64_compat
from pandas.tslib import OutOfBoundsDatetime
from pandas import (
    Series, DataFrame, Timestamp, Index, Int64Index, Float64Index, date_range,
    bdate_range, offsets, datetime, DatetimeIndex, _np_version_under1p8, NaT,
    notnull, PeriodIndex, _np_version_under1p10, Timedelta, TimedeltaIndex)

from pandas.util.testing import assert_series_equal, assert_almost_equal

from pandas.tests.test_base import Ops

randn = np.random.randn


class TestDatetimeIndex(tm.TestCase):
    _multiprocess_can_split_ = True

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

    def test_construction_with_ndarray(self):
        # GH 5152
        dates = [datetime(2013, 10, 7),
                 datetime(2013, 10, 8),
                 datetime(2013, 10, 9)]
        data = DatetimeIndex(dates, freq=pd.tseries.frequencies.BDay()).values
        result = DatetimeIndex(data, freq=pd.tseries.frequencies.BDay())
        expected = DatetimeIndex(['2013-10-07',
                                  '2013-10-08',
                                  '2013-10-09'],
                                 freq='B')
        tm.assert_index_equal(result, expected)

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

    def test_difference_freq(self):
        # GH14323: difference of DatetimeIndex should not preserve frequency

        index = date_range("20160920", "20160925", freq="D")
        other = date_range("20160921", "20160924", freq="D")
        expected = DatetimeIndex(["20160920", "20160925"], freq=None)
        idx_diff = index.difference(other)
        tm.assert_index_equal(idx_diff, expected)
        tm.assert_attr_equal('freq', idx_diff, expected)

        other = date_range("20160922", "20160925", freq="D")
        idx_diff = index.difference(other)
        expected = DatetimeIndex(["20160920", "20160921"], freq=None)
        tm.assert_index_equal(idx_diff, expected)
        tm.assert_attr_equal('freq', idx_diff, expected)

    def test_week_of_month_frequency(self):
        # GH 5348: "ValueError: Could not evaluate WOM-1SUN" shouldn't raise
        d1 = date(2002, 9, 1)
        d2 = date(2013, 10, 27)
        d3 = date(2012, 9, 30)
        idx1 = DatetimeIndex([d1, d2])
        idx2 = DatetimeIndex([d3])
        result_append = idx1.append(idx2)
        expected = DatetimeIndex([d1, d2, d3])
        tm.assert_index_equal(result_append, expected)
        result_union = idx1.union(idx2)
        expected = DatetimeIndex([d1, d3, d2])
        tm.assert_index_equal(result_union, expected)

        # GH 5115
        result = date_range("2013-1-1", periods=4, freq='WOM-1SAT')
        dates = ['2013-01-05', '2013-02-02', '2013-03-02', '2013-04-06']
        expected = DatetimeIndex(dates, freq='WOM-1SAT')
        tm.assert_index_equal(result, expected)

    def test_hash_error(self):
        index = date_range('20010101', periods=10)
        with tm.assertRaisesRegexp(TypeError, "unhashable type: %r" %
                                   type(index).__name__):
            hash(index)

    def test_stringified_slice_with_tz(self):
        # GH2658
        import datetime
        start = datetime.datetime.now()
        idx = DatetimeIndex(start=start, freq="1d", periods=10)
        df = DataFrame(lrange(10), index=idx)
        df["2013-01-14 23:44:34.437768-05:00":]  # no exception here

    def test_append_join_nondatetimeindex(self):
        rng = date_range('1/1/2000', periods=10)
        idx = Index(['a', 'b', 'c', 'd'])

        result = rng.append(idx)
        tm.assertIsInstance(result[0], Timestamp)

        # it works
        rng.join(idx, how='outer')

    def test_to_period_nofreq(self):
        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-04'])
        self.assertRaises(ValueError, idx.to_period)

        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-03'],
                            freq='infer')
        self.assertEqual(idx.freqstr, 'D')
        expected = pd.PeriodIndex(['2000-01-01', '2000-01-02',
                                   '2000-01-03'], freq='D')
        tm.assert_index_equal(idx.to_period(), expected)

        # GH 7606
        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-03'])
        self.assertEqual(idx.freqstr, None)
        tm.assert_index_equal(idx.to_period(), expected)

    def test_000constructor_resolution(self):
        # 2252
        t1 = Timestamp((1352934390 * 1000000000) + 1000000 + 1000 + 1)
        idx = DatetimeIndex([t1])

        self.assertEqual(idx.nanosecond[0], t1.nanosecond)

    def test_constructor_coverage(self):
        rng = date_range('1/1/2000', periods=10.5)
        exp = date_range('1/1/2000', periods=10)
        tm.assert_index_equal(rng, exp)

        self.assertRaises(ValueError, DatetimeIndex, start='1/1/2000',
                          periods='foo', freq='D')

        self.assertRaises(ValueError, DatetimeIndex, start='1/1/2000',
                          end='1/10/2000')

        self.assertRaises(ValueError, DatetimeIndex, '1/1/2000')

        # generator expression
        gen = (datetime(2000, 1, 1) + timedelta(i) for i in range(10))
        result = DatetimeIndex(gen)
        expected = DatetimeIndex([datetime(2000, 1, 1) + timedelta(i)
                                  for i in range(10)])
        tm.assert_index_equal(result, expected)

        # NumPy string array
        strings = np.array(['2000-01-01', '2000-01-02', '2000-01-03'])
        result = DatetimeIndex(strings)
        expected = DatetimeIndex(strings.astype('O'))
        tm.assert_index_equal(result, expected)

        from_ints = DatetimeIndex(expected.asi8)
        tm.assert_index_equal(from_ints, expected)

        # string with NaT
        strings = np.array(['2000-01-01', '2000-01-02', 'NaT'])
        result = DatetimeIndex(strings)
        expected = DatetimeIndex(strings.astype('O'))
        tm.assert_index_equal(result, expected)

        from_ints = DatetimeIndex(expected.asi8)
        tm.assert_index_equal(from_ints, expected)

        # non-conforming
        self.assertRaises(ValueError, DatetimeIndex,
                          ['2000-01-01', '2000-01-02', '2000-01-04'], freq='D')

        self.assertRaises(ValueError, DatetimeIndex, start='2011-01-01',
                          freq='b')
        self.assertRaises(ValueError, DatetimeIndex, end='2011-01-01',
                          freq='B')
        self.assertRaises(ValueError, DatetimeIndex, periods=10, freq='D')

    def test_constructor_datetime64_tzformat(self):
        # GH 6572
        tm._skip_if_no_pytz()
        import pytz
        # ISO 8601 format results in pytz.FixedOffset
        for freq in ['AS', 'W-SUN']:
            idx = date_range('2013-01-01T00:00:00-05:00',
                             '2016-01-01T23:59:59-05:00', freq=freq)
            expected = date_range('2013-01-01T00:00:00', '2016-01-01T23:59:59',
                                  freq=freq, tz=pytz.FixedOffset(-300))
            tm.assert_index_equal(idx, expected)
            # Unable to use `US/Eastern` because of DST
            expected_i8 = date_range('2013-01-01T00:00:00',
                                     '2016-01-01T23:59:59', freq=freq,
                                     tz='America/Lima')
            self.assert_numpy_array_equal(idx.asi8, expected_i8.asi8)

            idx = date_range('2013-01-01T00:00:00+09:00',
                             '2016-01-01T23:59:59+09:00', freq=freq)
            expected = date_range('2013-01-01T00:00:00', '2016-01-01T23:59:59',
                                  freq=freq, tz=pytz.FixedOffset(540))
            tm.assert_index_equal(idx, expected)
            expected_i8 = date_range('2013-01-01T00:00:00',
                                     '2016-01-01T23:59:59', freq=freq,
                                     tz='Asia/Tokyo')
            self.assert_numpy_array_equal(idx.asi8, expected_i8.asi8)

        tm._skip_if_no_dateutil()

        # Non ISO 8601 format results in dateutil.tz.tzoffset
        for freq in ['AS', 'W-SUN']:
            idx = date_range('2013/1/1 0:00:00-5:00', '2016/1/1 23:59:59-5:00',
                             freq=freq)
            expected = date_range('2013-01-01T00:00:00', '2016-01-01T23:59:59',
                                  freq=freq, tz=pytz.FixedOffset(-300))
            tm.assert_index_equal(idx, expected)
            # Unable to use `US/Eastern` because of DST
            expected_i8 = date_range('2013-01-01T00:00:00',
                                     '2016-01-01T23:59:59', freq=freq,
                                     tz='America/Lima')
            self.assert_numpy_array_equal(idx.asi8, expected_i8.asi8)

            idx = date_range('2013/1/1 0:00:00+9:00',
                             '2016/1/1 23:59:59+09:00', freq=freq)
            expected = date_range('2013-01-01T00:00:00', '2016-01-01T23:59:59',
                                  freq=freq, tz=pytz.FixedOffset(540))
            tm.assert_index_equal(idx, expected)
            expected_i8 = date_range('2013-01-01T00:00:00',
                                     '2016-01-01T23:59:59', freq=freq,
                                     tz='Asia/Tokyo')
            self.assert_numpy_array_equal(idx.asi8, expected_i8.asi8)

    def test_constructor_dtype(self):

        # passing a dtype with a tz should localize
        idx = DatetimeIndex(['2013-01-01', '2013-01-02'],
                            dtype='datetime64[ns, US/Eastern]')
        expected = DatetimeIndex(['2013-01-01', '2013-01-02']
                                 ).tz_localize('US/Eastern')
        tm.assert_index_equal(idx, expected)

        idx = DatetimeIndex(['2013-01-01', '2013-01-02'],
                            tz='US/Eastern')
        tm.assert_index_equal(idx, expected)

        # if we already have a tz and its not the same, then raise
        idx = DatetimeIndex(['2013-01-01', '2013-01-02'],
                            dtype='datetime64[ns, US/Eastern]')

        self.assertRaises(ValueError,
                          lambda: DatetimeIndex(idx,
                                                dtype='datetime64[ns]'))

        # this is effectively trying to convert tz's
        self.assertRaises(TypeError,
                          lambda: DatetimeIndex(idx,
                                                dtype='datetime64[ns, CET]'))
        self.assertRaises(ValueError,
                          lambda: DatetimeIndex(
                              idx, tz='CET',
                              dtype='datetime64[ns, US/Eastern]'))
        result = DatetimeIndex(idx, dtype='datetime64[ns, US/Eastern]')
        tm.assert_index_equal(idx, result)

    def test_constructor_name(self):
        idx = DatetimeIndex(start='2000-01-01', periods=1, freq='A',
                            name='TEST')
        self.assertEqual(idx.name, 'TEST')

    def test_comparisons_coverage(self):
        rng = date_range('1/1/2000', periods=10)

        # raise TypeError for now
        self.assertRaises(TypeError, rng.__lt__, rng[3].value)

        result = rng == list(rng)
        exp = rng == rng
        self.assert_numpy_array_equal(result, exp)

    def test_comparisons_nat(self):

        fidx1 = pd.Index([1.0, np.nan, 3.0, np.nan, 5.0, 7.0])
        fidx2 = pd.Index([2.0, 3.0, np.nan, np.nan, 6.0, 7.0])

        didx1 = pd.DatetimeIndex(['2014-01-01', pd.NaT, '2014-03-01', pd.NaT,
                                  '2014-05-01', '2014-07-01'])
        didx2 = pd.DatetimeIndex(['2014-02-01', '2014-03-01', pd.NaT, pd.NaT,
                                  '2014-06-01', '2014-07-01'])
        darr = np.array([np_datetime64_compat('2014-02-01 00:00Z'),
                         np_datetime64_compat('2014-03-01 00:00Z'),
                         np_datetime64_compat('nat'), np.datetime64('nat'),
                         np_datetime64_compat('2014-06-01 00:00Z'),
                         np_datetime64_compat('2014-07-01 00:00Z')])

        if _np_version_under1p8:
            # cannot test array because np.datetime('nat') returns today's date
            cases = [(fidx1, fidx2), (didx1, didx2)]
        else:
            cases = [(fidx1, fidx2), (didx1, didx2), (didx1, darr)]

        # Check pd.NaT is handles as the same as np.nan
        with tm.assert_produces_warning(None):
            for idx1, idx2 in cases:

                result = idx1 < idx2
                expected = np.array([True, False, False, False, True, False])
                self.assert_numpy_array_equal(result, expected)

                result = idx2 > idx1
                expected = np.array([True, False, False, False, True, False])
                self.assert_numpy_array_equal(result, expected)

                result = idx1 <= idx2
                expected = np.array([True, False, False, False, True, True])
                self.assert_numpy_array_equal(result, expected)

                result = idx2 >= idx1
                expected = np.array([True, False, False, False, True, True])
                self.assert_numpy_array_equal(result, expected)

                result = idx1 == idx2
                expected = np.array([False, False, False, False, False, True])
                self.assert_numpy_array_equal(result, expected)

                result = idx1 != idx2
                expected = np.array([True, True, True, True, True, False])
                self.assert_numpy_array_equal(result, expected)

        with tm.assert_produces_warning(None):
            for idx1, val in [(fidx1, np.nan), (didx1, pd.NaT)]:
                result = idx1 < val
                expected = np.array([False, False, False, False, False, False])
                self.assert_numpy_array_equal(result, expected)
                result = idx1 > val
                self.assert_numpy_array_equal(result, expected)

                result = idx1 <= val
                self.assert_numpy_array_equal(result, expected)
                result = idx1 >= val
                self.assert_numpy_array_equal(result, expected)

                result = idx1 == val
                self.assert_numpy_array_equal(result, expected)

                result = idx1 != val
                expected = np.array([True, True, True, True, True, True])
                self.assert_numpy_array_equal(result, expected)

        # Check pd.NaT is handles as the same as np.nan
        with tm.assert_produces_warning(None):
            for idx1, val in [(fidx1, 3), (didx1, datetime(2014, 3, 1))]:
                result = idx1 < val
                expected = np.array([True, False, False, False, False, False])
                self.assert_numpy_array_equal(result, expected)
                result = idx1 > val
                expected = np.array([False, False, False, False, True, True])
                self.assert_numpy_array_equal(result, expected)

                result = idx1 <= val
                expected = np.array([True, False, True, False, False, False])
                self.assert_numpy_array_equal(result, expected)
                result = idx1 >= val
                expected = np.array([False, False, True, False, True, True])
                self.assert_numpy_array_equal(result, expected)

                result = idx1 == val
                expected = np.array([False, False, True, False, False, False])
                self.assert_numpy_array_equal(result, expected)

                result = idx1 != val
                expected = np.array([True, True, False, True, True, True])
                self.assert_numpy_array_equal(result, expected)

    def test_map(self):
        rng = date_range('1/1/2000', periods=10)

        f = lambda x: x.strftime('%Y%m%d')
        result = rng.map(f)
        exp = Index([f(x) for x in rng], dtype='<U8')
        tm.assert_index_equal(result, exp)

    def test_iteration_preserves_tz(self):

        tm._skip_if_no_dateutil()

        # GH 8890
        import dateutil
        index = date_range("2012-01-01", periods=3, freq='H', tz='US/Eastern')

        for i, ts in enumerate(index):
            result = ts
            expected = index[i]
            self.assertEqual(result, expected)

        index = date_range("2012-01-01", periods=3, freq='H',
                           tz=dateutil.tz.tzoffset(None, -28800))

        for i, ts in enumerate(index):
            result = ts
            expected = index[i]
            self.assertEqual(result._repr_base, expected._repr_base)
            self.assertEqual(result, expected)

        # 9100
        index = pd.DatetimeIndex(['2014-12-01 03:32:39.987000-08:00',
                                  '2014-12-01 04:12:34.987000-08:00'])
        for i, ts in enumerate(index):
            result = ts
            expected = index[i]
            self.assertEqual(result._repr_base, expected._repr_base)
            self.assertEqual(result, expected)

    def test_misc_coverage(self):
        rng = date_range('1/1/2000', periods=5)
        result = rng.groupby(rng.day)
        tm.assertIsInstance(list(result.values())[0][0], Timestamp)

        idx = DatetimeIndex(['2000-01-03', '2000-01-01', '2000-01-02'])
        self.assertFalse(idx.equals(list(idx)))

        non_datetime = Index(list('abc'))
        self.assertFalse(idx.equals(list(non_datetime)))

    def test_union_coverage(self):
        idx = DatetimeIndex(['2000-01-03', '2000-01-01', '2000-01-02'])
        ordered = DatetimeIndex(idx.sort_values(), freq='infer')
        result = ordered.union(idx)
        tm.assert_index_equal(result, ordered)

        result = ordered[:0].union(ordered)
        tm.assert_index_equal(result, ordered)
        self.assertEqual(result.freq, ordered.freq)

    def test_union_bug_1730(self):
        rng_a = date_range('1/1/2012', periods=4, freq='3H')
        rng_b = date_range('1/1/2012', periods=4, freq='4H')

        result = rng_a.union(rng_b)
        exp = DatetimeIndex(sorted(set(list(rng_a)) | set(list(rng_b))))
        tm.assert_index_equal(result, exp)

    def test_union_bug_1745(self):
        left = DatetimeIndex(['2012-05-11 15:19:49.695000'])
        right = DatetimeIndex(['2012-05-29 13:04:21.322000',
                               '2012-05-11 15:27:24.873000',
                               '2012-05-11 15:31:05.350000'])

        result = left.union(right)
        exp = DatetimeIndex(sorted(set(list(left)) | set(list(right))))
        tm.assert_index_equal(result, exp)

    def test_union_bug_4564(self):
        from pandas import DateOffset
        left = date_range("2013-01-01", "2013-02-01")
        right = left + DateOffset(minutes=15)

        result = left.union(right)
        exp = DatetimeIndex(sorted(set(list(left)) | set(list(right))))
        tm.assert_index_equal(result, exp)

    def test_union_freq_both_none(self):
        # GH11086
        expected = bdate_range('20150101', periods=10)
        expected.freq = None

        result = expected.union(expected)
        tm.assert_index_equal(result, expected)
        self.assertIsNone(result.freq)

    def test_union_dataframe_index(self):
        rng1 = date_range('1/1/1999', '1/1/2012', freq='MS')
        s1 = Series(np.random.randn(len(rng1)), rng1)

        rng2 = date_range('1/1/1980', '12/1/2001', freq='MS')
        s2 = Series(np.random.randn(len(rng2)), rng2)
        df = DataFrame({'s1': s1, 's2': s2})

        exp = pd.date_range('1/1/1980', '1/1/2012', freq='MS')
        tm.assert_index_equal(df.index, exp)

    def test_intersection_bug_1708(self):
        from pandas import DateOffset
        index_1 = date_range('1/1/2012', periods=4, freq='12H')
        index_2 = index_1 + DateOffset(hours=1)

        result = index_1 & index_2
        self.assertEqual(len(result), 0)

    def test_intersection(self):
        # GH 4690 (with tz)
        for tz in [None, 'Asia/Tokyo', 'US/Eastern', 'dateutil/US/Pacific']:
            base = date_range('6/1/2000', '6/30/2000', freq='D', name='idx')

            # if target has the same name, it is preserved
            rng2 = date_range('5/15/2000', '6/20/2000', freq='D', name='idx')
            expected2 = date_range('6/1/2000', '6/20/2000', freq='D',
                                   name='idx')

            # if target name is different, it will be reset
            rng3 = date_range('5/15/2000', '6/20/2000', freq='D', name='other')
            expected3 = date_range('6/1/2000', '6/20/2000', freq='D',
                                   name=None)

            rng4 = date_range('7/1/2000', '7/31/2000', freq='D', name='idx')
            expected4 = DatetimeIndex([], name='idx')

            for (rng, expected) in [(rng2, expected2), (rng3, expected3),
                                    (rng4, expected4)]:
                result = base.intersection(rng)
                tm.assert_index_equal(result, expected)
                self.assertEqual(result.name, expected.name)
                self.assertEqual(result.freq, expected.freq)
                self.assertEqual(result.tz, expected.tz)

            # non-monotonic
            base = DatetimeIndex(['2011-01-05', '2011-01-04',
                                  '2011-01-02', '2011-01-03'],
                                 tz=tz, name='idx')

            rng2 = DatetimeIndex(['2011-01-04', '2011-01-02',
                                  '2011-02-02', '2011-02-03'],
                                 tz=tz, name='idx')
            expected2 = DatetimeIndex(
                ['2011-01-04', '2011-01-02'], tz=tz, name='idx')

            rng3 = DatetimeIndex(['2011-01-04', '2011-01-02',
                                  '2011-02-02', '2011-02-03'],
                                 tz=tz, name='other')
            expected3 = DatetimeIndex(
                ['2011-01-04', '2011-01-02'], tz=tz, name=None)

            # GH 7880
            rng4 = date_range('7/1/2000', '7/31/2000', freq='D', tz=tz,
                              name='idx')
            expected4 = DatetimeIndex([], tz=tz, name='idx')

            for (rng, expected) in [(rng2, expected2), (rng3, expected3),
                                    (rng4, expected4)]:
                result = base.intersection(rng)
                tm.assert_index_equal(result, expected)
                self.assertEqual(result.name, expected.name)
                self.assertIsNone(result.freq)
                self.assertEqual(result.tz, expected.tz)

        # empty same freq GH2129
        rng = date_range('6/1/2000', '6/15/2000', freq='T')
        result = rng[0:0].intersection(rng)
        self.assertEqual(len(result), 0)

        result = rng.intersection(rng[0:0])
        self.assertEqual(len(result), 0)

    def test_string_index_series_name_converted(self):
        # #1644
        df = DataFrame(np.random.randn(10, 4),
                       index=date_range('1/1/2000', periods=10))

        result = df.loc['1/3/2000']
        self.assertEqual(result.name, df.index[2])

        result = df.T['1/3/2000']
        self.assertEqual(result.name, df.index[2])

    # GH 10699
    def test_datetime64_with_DateOffset(self):
        for klass, assert_func in zip([Series, DatetimeIndex],
                                      [self.assert_series_equal,
                                       tm.assert_index_equal]):
            s = klass(date_range('2000-01-01', '2000-01-31'), name='a')
            result = s + pd.DateOffset(years=1)
            result2 = pd.DateOffset(years=1) + s
            exp = klass(date_range('2001-01-01', '2001-01-31'), name='a')
            assert_func(result, exp)
            assert_func(result2, exp)

            result = s - pd.DateOffset(years=1)
            exp = klass(date_range('1999-01-01', '1999-01-31'), name='a')
            assert_func(result, exp)

            s = klass([Timestamp('2000-01-15 00:15:00', tz='US/Central'),
                       pd.Timestamp('2000-02-15', tz='US/Central')], name='a')
            result = s + pd.offsets.Day()
            result2 = pd.offsets.Day() + s
            exp = klass([Timestamp('2000-01-16 00:15:00', tz='US/Central'),
                         Timestamp('2000-02-16', tz='US/Central')], name='a')
            assert_func(result, exp)
            assert_func(result2, exp)

            s = klass([Timestamp('2000-01-15 00:15:00', tz='US/Central'),
                       pd.Timestamp('2000-02-15', tz='US/Central')], name='a')
            result = s + pd.offsets.MonthEnd()
            result2 = pd.offsets.MonthEnd() + s
            exp = klass([Timestamp('2000-01-31 00:15:00', tz='US/Central'),
                         Timestamp('2000-02-29', tz='US/Central')], name='a')
            assert_func(result, exp)
            assert_func(result2, exp)

            # array of offsets - valid for Series only
            if klass is Series:
                with tm.assert_produces_warning(PerformanceWarning):
                    s = klass([Timestamp('2000-1-1'), Timestamp('2000-2-1')])
                    result = s + Series([pd.offsets.DateOffset(years=1),
                                         pd.offsets.MonthEnd()])
                    exp = klass([Timestamp('2001-1-1'), Timestamp('2000-2-29')
                                 ])
                    assert_func(result, exp)

                    # same offset
                    result = s + Series([pd.offsets.DateOffset(years=1),
                                         pd.offsets.DateOffset(years=1)])
                    exp = klass([Timestamp('2001-1-1'), Timestamp('2001-2-1')])
                    assert_func(result, exp)

            s = klass([Timestamp('2000-01-05 00:15:00'),
                       Timestamp('2000-01-31 00:23:00'),
                       Timestamp('2000-01-01'),
                       Timestamp('2000-03-31'),
                       Timestamp('2000-02-29'),
                       Timestamp('2000-12-31'),
                       Timestamp('2000-05-15'),
                       Timestamp('2001-06-15')])

            # DateOffset relativedelta fastpath
            relative_kwargs = [('years', 2), ('months', 5), ('days', 3),
                               ('hours', 5), ('minutes', 10), ('seconds', 2),
                               ('microseconds', 5)]
            for i, kwd in enumerate(relative_kwargs):
                op = pd.DateOffset(**dict([kwd]))
                assert_func(klass([x + op for x in s]), s + op)
                assert_func(klass([x - op for x in s]), s - op)
                op = pd.DateOffset(**dict(relative_kwargs[:i + 1]))
                assert_func(klass([x + op for x in s]), s + op)
                assert_func(klass([x - op for x in s]), s - op)

            # assert these are equal on a piecewise basis
            offsets = ['YearBegin', ('YearBegin', {'month': 5}), 'YearEnd',
                       ('YearEnd', {'month': 5}), 'MonthBegin', 'MonthEnd',
                       'SemiMonthEnd', 'SemiMonthBegin',
                       'Week', ('Week', {
                           'weekday': 3
                       }), 'BusinessDay', 'BDay', 'QuarterEnd', 'QuarterBegin',
                       'CustomBusinessDay', 'CDay', 'CBMonthEnd',
                       'CBMonthBegin', 'BMonthBegin', 'BMonthEnd',
                       'BusinessHour', 'BYearBegin', 'BYearEnd',
                       'BQuarterBegin', ('LastWeekOfMonth', {
                           'weekday': 2
                       }), ('FY5253Quarter', {'qtr_with_extra_week': 1,
                                              'startingMonth': 1,
                                              'weekday': 2,
                                              'variation': 'nearest'}),
                       ('FY5253', {'weekday': 0,
                                   'startingMonth': 2,
                                   'variation':
                                   'nearest'}), ('WeekOfMonth', {'weekday': 2,
                                                                 'week': 2}),
                       'Easter', ('DateOffset', {'day': 4}),
                       ('DateOffset', {'month': 5})]

            with warnings.catch_warnings(record=True):
                for normalize in (True, False):
                    for do in offsets:
                        if isinstance(do, tuple):
                            do, kwargs = do
                        else:
                            do = do
                            kwargs = {}

                        for n in [0, 5]:
                            if (do in ['WeekOfMonth', 'LastWeekOfMonth',
                                       'FY5253Quarter', 'FY5253'] and n == 0):
                                continue
                            op = getattr(pd.offsets, do)(n,
                                                         normalize=normalize,
                                                         **kwargs)
                            assert_func(klass([x + op for x in s]), s + op)
                            assert_func(klass([x - op for x in s]), s - op)
                            assert_func(klass([op + x for x in s]), op + s)

    def test_overflow_offset(self):
        # xref https://github.com/statsmodels/statsmodels/issues/3374
        # ends up multiplying really large numbers which overflow

        t = Timestamp('2017-01-13 00:00:00', freq='D')
        offset = 20169940 * pd.offsets.Day(1)

        def f():
            t + offset
        self.assertRaises(OverflowError, f)

        def f():
            offset + t
        self.assertRaises(OverflowError, f)

        def f():
            t - offset
        self.assertRaises(OverflowError, f)

    def test_get_duplicates(self):
        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-02',
                             '2000-01-03', '2000-01-03', '2000-01-04'])

        result = idx.get_duplicates()
        ex = DatetimeIndex(['2000-01-02', '2000-01-03'])
        tm.assert_index_equal(result, ex)

    def test_argmin_argmax(self):
        idx = DatetimeIndex(['2000-01-04', '2000-01-01', '2000-01-02'])
        self.assertEqual(idx.argmin(), 1)
        self.assertEqual(idx.argmax(), 0)

    def test_sort_values(self):
        idx = DatetimeIndex(['2000-01-04', '2000-01-01', '2000-01-02'])

        ordered = idx.sort_values()
        self.assertTrue(ordered.is_monotonic)

        ordered = idx.sort_values(ascending=False)
        self.assertTrue(ordered[::-1].is_monotonic)

        ordered, dexer = idx.sort_values(return_indexer=True)
        self.assertTrue(ordered.is_monotonic)
        self.assert_numpy_array_equal(dexer,
                                      np.array([1, 2, 0], dtype=np.intp))

        ordered, dexer = idx.sort_values(return_indexer=True, ascending=False)
        self.assertTrue(ordered[::-1].is_monotonic)
        self.assert_numpy_array_equal(dexer,
                                      np.array([0, 2, 1], dtype=np.intp))

    def test_round(self):

        # round
        dt = Timestamp('20130101 09:10:11')
        result = dt.round('D')
        expected = Timestamp('20130101')
        self.assertEqual(result, expected)

        dt = Timestamp('20130101 19:10:11')
        result = dt.round('D')
        expected = Timestamp('20130102')
        self.assertEqual(result, expected)

        dt = Timestamp('20130201 12:00:00')
        result = dt.round('D')
        expected = Timestamp('20130202')
        self.assertEqual(result, expected)

        dt = Timestamp('20130104 12:00:00')
        result = dt.round('D')
        expected = Timestamp('20130105')
        self.assertEqual(result, expected)

        dt = Timestamp('20130104 12:32:00')
        result = dt.round('30Min')
        expected = Timestamp('20130104 12:30:00')
        self.assertEqual(result, expected)

        dti = date_range('20130101 09:10:11', periods=5)
        result = dti.round('D')
        expected = date_range('20130101', periods=5)
        tm.assert_index_equal(result, expected)

        # floor
        dt = Timestamp('20130101 09:10:11')
        result = dt.floor('D')
        expected = Timestamp('20130101')
        self.assertEqual(result, expected)

        # ceil
        dt = Timestamp('20130101 09:10:11')
        result = dt.ceil('D')
        expected = Timestamp('20130102')
        self.assertEqual(result, expected)

        # round with tz
        dt = Timestamp('20130101 09:10:11', tz='US/Eastern')
        result = dt.round('D')
        expected = Timestamp('20130101', tz='US/Eastern')
        self.assertEqual(result, expected)

        dt = Timestamp('20130101 09:10:11', tz='US/Eastern')
        result = dt.round('s')
        self.assertEqual(result, dt)

        dti = date_range('20130101 09:10:11',
                         periods=5).tz_localize('UTC').tz_convert('US/Eastern')
        result = dti.round('D')
        expected = date_range('20130101', periods=5).tz_localize('US/Eastern')
        tm.assert_index_equal(result, expected)

        result = dti.round('s')
        tm.assert_index_equal(result, dti)

        # invalid
        for freq in ['Y', 'M', 'foobar']:
            self.assertRaises(ValueError, lambda: dti.round(freq))

    def test_insert(self):
        idx = DatetimeIndex(
            ['2000-01-04', '2000-01-01', '2000-01-02'], name='idx')

        result = idx.insert(2, datetime(2000, 1, 5))
        exp = DatetimeIndex(['2000-01-04', '2000-01-01', '2000-01-05',
                             '2000-01-02'], name='idx')
        tm.assert_index_equal(result, exp)

        # insertion of non-datetime should coerce to object index
        result = idx.insert(1, 'inserted')
        expected = Index([datetime(2000, 1, 4), 'inserted',
                          datetime(2000, 1, 1),
                          datetime(2000, 1, 2)], name='idx')
        self.assertNotIsInstance(result, DatetimeIndex)
        tm.assert_index_equal(result, expected)
        self.assertEqual(result.name, expected.name)

        idx = date_range('1/1/2000', periods=3, freq='M', name='idx')

        # preserve freq
        expected_0 = DatetimeIndex(['1999-12-31', '2000-01-31', '2000-02-29',
                                    '2000-03-31'], name='idx', freq='M')
        expected_3 = DatetimeIndex(['2000-01-31', '2000-02-29', '2000-03-31',
                                    '2000-04-30'], name='idx', freq='M')

        # reset freq to None
        expected_1_nofreq = DatetimeIndex(['2000-01-31', '2000-01-31',
                                           '2000-02-29',
                                           '2000-03-31'], name='idx',
                                          freq=None)
        expected_3_nofreq = DatetimeIndex(['2000-01-31', '2000-02-29',
                                           '2000-03-31',
                                           '2000-01-02'], name='idx',
                                          freq=None)

        cases = [(0, datetime(1999, 12, 31), expected_0),
                 (-3, datetime(1999, 12, 31), expected_0),
                 (3, datetime(2000, 4, 30), expected_3),
                 (1, datetime(2000, 1, 31), expected_1_nofreq),
                 (3, datetime(2000, 1, 2), expected_3_nofreq)]

        for n, d, expected in cases:
            result = idx.insert(n, d)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)

        # reset freq to None
        result = idx.insert(3, datetime(2000, 1, 2))
        expected = DatetimeIndex(['2000-01-31', '2000-02-29', '2000-03-31',
                                  '2000-01-02'], name='idx', freq=None)
        tm.assert_index_equal(result, expected)
        self.assertEqual(result.name, expected.name)
        self.assertTrue(result.freq is None)

        # GH 7299
        tm._skip_if_no_pytz()
        import pytz

        idx = date_range('1/1/2000', periods=3, freq='D', tz='Asia/Tokyo',
                         name='idx')
        with tm.assertRaises(ValueError):
            result = idx.insert(3, pd.Timestamp('2000-01-04'))
        with tm.assertRaises(ValueError):
            result = idx.insert(3, datetime(2000, 1, 4))
        with tm.assertRaises(ValueError):
            result = idx.insert(3, pd.Timestamp('2000-01-04', tz='US/Eastern'))
        with tm.assertRaises(ValueError):
            result = idx.insert(3,
                                datetime(2000, 1, 4,
                                         tzinfo=pytz.timezone('US/Eastern')))

        for tz in ['US/Pacific', 'Asia/Singapore']:
            idx = date_range('1/1/2000 09:00', periods=6, freq='H', tz=tz,
                             name='idx')
            # preserve freq
            expected = date_range('1/1/2000 09:00', periods=7, freq='H', tz=tz,
                                  name='idx')
            for d in [pd.Timestamp('2000-01-01 15:00', tz=tz),
                      pytz.timezone(tz).localize(datetime(2000, 1, 1, 15))]:

                result = idx.insert(6, d)
                tm.assert_index_equal(result, expected)
                self.assertEqual(result.name, expected.name)
                self.assertEqual(result.freq, expected.freq)
                self.assertEqual(result.tz, expected.tz)

            expected = DatetimeIndex(['2000-01-01 09:00', '2000-01-01 10:00',
                                      '2000-01-01 11:00',
                                      '2000-01-01 12:00', '2000-01-01 13:00',
                                      '2000-01-01 14:00',
                                      '2000-01-01 10:00'], name='idx',
                                     tz=tz, freq=None)
            # reset freq to None
            for d in [pd.Timestamp('2000-01-01 10:00', tz=tz),
                      pytz.timezone(tz).localize(datetime(2000, 1, 1, 10))]:
                result = idx.insert(6, d)
                tm.assert_index_equal(result, expected)
                self.assertEqual(result.name, expected.name)
                self.assertTrue(result.freq is None)
                self.assertEqual(result.tz, expected.tz)

    def test_delete(self):
        idx = date_range(start='2000-01-01', periods=5, freq='M', name='idx')

        # prserve freq
        expected_0 = date_range(start='2000-02-01', periods=4, freq='M',
                                name='idx')
        expected_4 = date_range(start='2000-01-01', periods=4, freq='M',
                                name='idx')

        # reset freq to None
        expected_1 = DatetimeIndex(['2000-01-31', '2000-03-31', '2000-04-30',
                                    '2000-05-31'], freq=None, name='idx')

        cases = {0: expected_0,
                 -5: expected_0,
                 -1: expected_4,
                 4: expected_4,
                 1: expected_1}
        for n, expected in compat.iteritems(cases):
            result = idx.delete(n)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)

        with tm.assertRaises((IndexError, ValueError)):
            # either depeidnig on numpy version
            result = idx.delete(5)

        for tz in [None, 'Asia/Tokyo', 'US/Pacific']:
            idx = date_range(start='2000-01-01 09:00', periods=10, freq='H',
                             name='idx', tz=tz)

            expected = date_range(start='2000-01-01 10:00', periods=9,
                                  freq='H', name='idx', tz=tz)
            result = idx.delete(0)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freqstr, 'H')
            self.assertEqual(result.tz, expected.tz)

            expected = date_range(start='2000-01-01 09:00', periods=9,
                                  freq='H', name='idx', tz=tz)
            result = idx.delete(-1)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freqstr, 'H')
            self.assertEqual(result.tz, expected.tz)

    def test_delete_slice(self):
        idx = date_range(start='2000-01-01', periods=10, freq='D', name='idx')

        # prserve freq
        expected_0_2 = date_range(start='2000-01-04', periods=7, freq='D',
                                  name='idx')
        expected_7_9 = date_range(start='2000-01-01', periods=7, freq='D',
                                  name='idx')

        # reset freq to None
        expected_3_5 = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-03',
                                      '2000-01-07', '2000-01-08', '2000-01-09',
                                      '2000-01-10'], freq=None, name='idx')

        cases = {(0, 1, 2): expected_0_2,
                 (7, 8, 9): expected_7_9,
                 (3, 4, 5): expected_3_5}
        for n, expected in compat.iteritems(cases):
            result = idx.delete(n)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)

            result = idx.delete(slice(n[0], n[-1] + 1))
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)

        for tz in [None, 'Asia/Tokyo', 'US/Pacific']:
            ts = pd.Series(1, index=pd.date_range(
                '2000-01-01 09:00', periods=10, freq='H', name='idx', tz=tz))
            # preserve freq
            result = ts.drop(ts.index[:5]).index
            expected = pd.date_range('2000-01-01 14:00', periods=5, freq='H',
                                     name='idx', tz=tz)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)
            self.assertEqual(result.tz, expected.tz)

            # reset freq to None
            result = ts.drop(ts.index[[1, 3, 5, 7, 9]]).index
            expected = DatetimeIndex(['2000-01-01 09:00', '2000-01-01 11:00',
                                      '2000-01-01 13:00',
                                      '2000-01-01 15:00', '2000-01-01 17:00'],
                                     freq=None, name='idx', tz=tz)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)
            self.assertEqual(result.tz, expected.tz)

    def test_take(self):
        dates = [datetime(2010, 1, 1, 14), datetime(2010, 1, 1, 15),
                 datetime(2010, 1, 1, 17), datetime(2010, 1, 1, 21)]

        for tz in [None, 'US/Eastern', 'Asia/Tokyo']:
            idx = DatetimeIndex(start='2010-01-01 09:00',
                                end='2010-02-01 09:00', freq='H', tz=tz,
                                name='idx')
            expected = DatetimeIndex(dates, freq=None, name='idx', tz=tz)

            taken1 = idx.take([5, 6, 8, 12])
            taken2 = idx[[5, 6, 8, 12]]

            for taken in [taken1, taken2]:
                tm.assert_index_equal(taken, expected)
                tm.assertIsInstance(taken, DatetimeIndex)
                self.assertIsNone(taken.freq)
                self.assertEqual(taken.tz, expected.tz)
                self.assertEqual(taken.name, expected.name)

    def test_take_fill_value(self):
        # GH 12631
        idx = pd.DatetimeIndex(['2011-01-01', '2011-02-01', '2011-03-01'],
                               name='xxx')
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.DatetimeIndex(['2011-02-01', '2011-01-01', '2011-03-01'],
                                    name='xxx')
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.DatetimeIndex(['2011-02-01', '2011-01-01', 'NaT'],
                                    name='xxx')
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        expected = pd.DatetimeIndex(['2011-02-01', '2011-01-01', '2011-03-01'],
                                    name='xxx')
        tm.assert_index_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with tm.assertRaises(IndexError):
            idx.take(np.array([1, -5]))

    def test_take_fill_value_with_timezone(self):
        idx = pd.DatetimeIndex(['2011-01-01', '2011-02-01', '2011-03-01'],
                               name='xxx', tz='US/Eastern')
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.DatetimeIndex(['2011-02-01', '2011-01-01', '2011-03-01'],
                                    name='xxx', tz='US/Eastern')
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.DatetimeIndex(['2011-02-01', '2011-01-01', 'NaT'],
                                    name='xxx', tz='US/Eastern')
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        expected = pd.DatetimeIndex(['2011-02-01', '2011-01-01', '2011-03-01'],
                                    name='xxx', tz='US/Eastern')
        tm.assert_index_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with tm.assertRaises(IndexError):
            idx.take(np.array([1, -5]))

    def test_map_bug_1677(self):
        index = DatetimeIndex(['2012-04-25 09:30:00.393000'])
        f = index.asof

        result = index.map(f)
        expected = Index([f(index[0])])
        tm.assert_index_equal(result, expected)

    def test_groupby_function_tuple_1677(self):
        df = DataFrame(np.random.rand(100),
                       index=date_range("1/1/2000", periods=100))
        monthly_group = df.groupby(lambda x: (x.year, x.month))

        result = monthly_group.mean()
        tm.assertIsInstance(result.index[0], tuple)

    def test_append_numpy_bug_1681(self):
        # another datetime64 bug
        dr = date_range('2011/1/1', '2012/1/1', freq='W-FRI')
        a = DataFrame()
        c = DataFrame({'A': 'foo', 'B': dr}, index=dr)

        result = a.append(c)
        self.assertTrue((result['B'] == dr).all())

    def test_isin(self):
        index = tm.makeDateIndex(4)
        result = index.isin(index)
        self.assertTrue(result.all())

        result = index.isin(list(index))
        self.assertTrue(result.all())

        assert_almost_equal(index.isin([index[2], 5]),
                            np.array([False, False, True, False]))

    def test_union(self):
        i1 = Int64Index(np.arange(0, 20, 2))
        i2 = Int64Index(np.arange(10, 30, 2))
        result = i1.union(i2)
        expected = Int64Index(np.arange(0, 30, 2))
        tm.assert_index_equal(result, expected)

    def test_union_with_DatetimeIndex(self):
        i1 = Int64Index(np.arange(0, 20, 2))
        i2 = DatetimeIndex(start='2012-01-03 00:00:00', periods=10, freq='D')
        i1.union(i2)  # Works
        i2.union(i1)  # Fails with "AttributeError: can't set attribute"

    def test_time(self):
        rng = pd.date_range('1/1/2000', freq='12min', periods=10)
        result = pd.Index(rng).time
        expected = [t.time() for t in rng]
        self.assertTrue((result == expected).all())

    def test_date(self):
        rng = pd.date_range('1/1/2000', freq='12H', periods=10)
        result = pd.Index(rng).date
        expected = [t.date() for t in rng]
        self.assertTrue((result == expected).all())

    def test_does_not_convert_mixed_integer(self):
        df = tm.makeCustomDataframe(10, 10,
                                    data_gen_f=lambda *args, **kwargs: randn(),
                                    r_idx_type='i', c_idx_type='dt')
        cols = df.columns.join(df.index, how='outer')
        joined = cols.join(df.columns)
        self.assertEqual(cols.dtype, np.dtype('O'))
        self.assertEqual(cols.dtype, joined.dtype)
        tm.assert_numpy_array_equal(cols.values, joined.values)

    def test_slice_keeps_name(self):
        # GH4226
        st = pd.Timestamp('2013-07-01 00:00:00', tz='America/Los_Angeles')
        et = pd.Timestamp('2013-07-02 00:00:00', tz='America/Los_Angeles')
        dr = pd.date_range(st, et, freq='H', name='timebucket')
        self.assertEqual(dr[1:].name, dr.name)

    def test_join_self(self):
        index = date_range('1/1/2000', periods=10)
        kinds = 'outer', 'inner', 'left', 'right'
        for kind in kinds:
            joined = index.join(index, how=kind)
            self.assertIs(index, joined)

    def assert_index_parameters(self, index):
        assert index.freq == '40960N'
        assert index.inferred_freq == '40960N'

    def test_ns_index(self):
        nsamples = 400
        ns = int(1e9 / 24414)
        dtstart = np.datetime64('2012-09-20T00:00:00')

        dt = dtstart + np.arange(nsamples) * np.timedelta64(ns, 'ns')
        freq = ns * offsets.Nano()
        index = pd.DatetimeIndex(dt, freq=freq, name='time')
        self.assert_index_parameters(index)

        new_index = pd.DatetimeIndex(start=index[0], end=index[-1],
                                     freq=index.freq)
        self.assert_index_parameters(new_index)

    def test_join_with_period_index(self):
        df = tm.makeCustomDataframe(
            10, 10, data_gen_f=lambda *args: np.random.randint(2),
            c_idx_type='p', r_idx_type='dt')
        s = df.iloc[:5, 0]
        joins = 'left', 'right', 'inner', 'outer'

        for join in joins:
            with tm.assertRaisesRegexp(ValueError, 'can only call with other '
                                       'PeriodIndex-ed objects'):
                df.columns.join(s.index, how=join)

    def test_factorize(self):
        idx1 = DatetimeIndex(['2014-01', '2014-01', '2014-02', '2014-02',
                              '2014-03', '2014-03'])

        exp_arr = np.array([0, 0, 1, 1, 2, 2], dtype=np.intp)
        exp_idx = DatetimeIndex(['2014-01', '2014-02', '2014-03'])

        arr, idx = idx1.factorize()
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        arr, idx = idx1.factorize(sort=True)
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        # tz must be preserved
        idx1 = idx1.tz_localize('Asia/Tokyo')
        exp_idx = exp_idx.tz_localize('Asia/Tokyo')

        arr, idx = idx1.factorize()
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        idx2 = pd.DatetimeIndex(['2014-03', '2014-03', '2014-02', '2014-01',
                                 '2014-03', '2014-01'])

        exp_arr = np.array([2, 2, 1, 0, 2, 0], dtype=np.intp)
        exp_idx = DatetimeIndex(['2014-01', '2014-02', '2014-03'])
        arr, idx = idx2.factorize(sort=True)
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        exp_arr = np.array([0, 0, 1, 2, 0, 2], dtype=np.intp)
        exp_idx = DatetimeIndex(['2014-03', '2014-02', '2014-01'])
        arr, idx = idx2.factorize()
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        # freq must be preserved
        idx3 = date_range('2000-01', periods=4, freq='M', tz='Asia/Tokyo')
        exp_arr = np.array([0, 1, 2, 3], dtype=np.intp)
        arr, idx = idx3.factorize()
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, idx3)

    def test_factorize_tz(self):
        # GH 13750
        for tz in [None, 'UTC', 'US/Eastern', 'Asia/Tokyo']:
            base = pd.date_range('2016-11-05', freq='H', periods=100, tz=tz)
            idx = base.repeat(5)

            exp_arr = np.arange(100, dtype=np.intp).repeat(5)

            for obj in [idx, pd.Series(idx)]:
                arr, res = obj.factorize()
                self.assert_numpy_array_equal(arr, exp_arr)
                tm.assert_index_equal(res, base)

    def test_factorize_dst(self):
        # GH 13750
        idx = pd.date_range('2016-11-06', freq='H', periods=12,
                            tz='US/Eastern')

        for obj in [idx, pd.Series(idx)]:
            arr, res = obj.factorize()
            self.assert_numpy_array_equal(arr, np.arange(12, dtype=np.intp))
            tm.assert_index_equal(res, idx)

        idx = pd.date_range('2016-06-13', freq='H', periods=12,
                            tz='US/Eastern')

        for obj in [idx, pd.Series(idx)]:
            arr, res = obj.factorize()
            self.assert_numpy_array_equal(arr, np.arange(12, dtype=np.intp))
            tm.assert_index_equal(res, idx)

    def test_slice_with_negative_step(self):
        ts = Series(np.arange(20),
                    date_range('2014-01-01', periods=20, freq='MS'))
        SLC = pd.IndexSlice

        def assert_slices_equivalent(l_slc, i_slc):
            assert_series_equal(ts[l_slc], ts.iloc[i_slc])
            assert_series_equal(ts.loc[l_slc], ts.iloc[i_slc])
            assert_series_equal(ts.loc[l_slc], ts.iloc[i_slc])

        assert_slices_equivalent(SLC[Timestamp('2014-10-01')::-1], SLC[9::-1])
        assert_slices_equivalent(SLC['2014-10-01'::-1], SLC[9::-1])

        assert_slices_equivalent(SLC[:Timestamp('2014-10-01'):-1], SLC[:8:-1])
        assert_slices_equivalent(SLC[:'2014-10-01':-1], SLC[:8:-1])

        assert_slices_equivalent(SLC['2015-02-01':'2014-10-01':-1],
                                 SLC[13:8:-1])
        assert_slices_equivalent(SLC[Timestamp('2015-02-01'):Timestamp(
            '2014-10-01'):-1], SLC[13:8:-1])
        assert_slices_equivalent(SLC['2015-02-01':Timestamp('2014-10-01'):-1],
                                 SLC[13:8:-1])
        assert_slices_equivalent(SLC[Timestamp('2015-02-01'):'2014-10-01':-1],
                                 SLC[13:8:-1])

        assert_slices_equivalent(SLC['2014-10-01':'2015-02-01':-1], SLC[:0])

    def test_slice_with_zero_step_raises(self):
        ts = Series(np.arange(20),
                    date_range('2014-01-01', periods=20, freq='MS'))
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: ts[::0])
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: ts.loc[::0])
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: ts.loc[::0])

    def test_slice_bounds_empty(self):
        # GH 14354
        empty_idx = DatetimeIndex(freq='1H', periods=0, end='2015')

        right = empty_idx._maybe_cast_slice_bound('2015-01-02', 'right', 'loc')
        exp = Timestamp('2015-01-02 23:59:59.999999999')
        self.assertEqual(right, exp)

        left = empty_idx._maybe_cast_slice_bound('2015-01-02', 'left', 'loc')
        exp = Timestamp('2015-01-02 00:00:00')
        self.assertEqual(left, exp)


class TestDatetimeIndexOps(Ops):
    tz = [None, 'UTC', 'Asia/Tokyo', 'US/Eastern', 'dateutil/Asia/Singapore',
          'dateutil/US/Pacific']

    def setUp(self):
        super(TestDatetimeIndexOps, self).setUp()
        mask = lambda x: (isinstance(x, DatetimeIndex) or
                          isinstance(x, PeriodIndex))
        self.is_valid_objs = [o for o in self.objs if mask(o)]
        self.not_valid_objs = [o for o in self.objs if not mask(o)]

    def test_ops_properties(self):
        self.check_ops_properties(
            ['year', 'month', 'day', 'hour', 'minute', 'second', 'weekofyear',
             'week', 'dayofweek', 'dayofyear', 'quarter'])
        self.check_ops_properties(['date', 'time', 'microsecond', 'nanosecond',
                                   'is_month_start', 'is_month_end',
                                   'is_quarter_start',
                                   'is_quarter_end', 'is_year_start',
                                   'is_year_end', 'weekday_name'],
                                  lambda x: isinstance(x, DatetimeIndex))

    def test_ops_properties_basic(self):

        # sanity check that the behavior didn't change
        # GH7206
        for op in ['year', 'day', 'second', 'weekday']:
            self.assertRaises(TypeError, lambda x: getattr(self.dt_series, op))

        # attribute access should still work!
        s = Series(dict(year=2000, month=1, day=10))
        self.assertEqual(s.year, 2000)
        self.assertEqual(s.month, 1)
        self.assertEqual(s.day, 10)
        self.assertRaises(AttributeError, lambda: s.weekday)

    def test_asobject_tolist(self):
        idx = pd.date_range(start='2013-01-01', periods=4, freq='M',
                            name='idx')
        expected_list = [Timestamp('2013-01-31'),
                         Timestamp('2013-02-28'),
                         Timestamp('2013-03-31'),
                         Timestamp('2013-04-30')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.asobject
        self.assertTrue(isinstance(result, Index))

        self.assertEqual(result.dtype, object)
        self.assert_index_equal(result, expected)
        self.assertEqual(result.name, expected.name)
        self.assertEqual(idx.tolist(), expected_list)

        idx = pd.date_range(start='2013-01-01', periods=4, freq='M',
                            name='idx', tz='Asia/Tokyo')
        expected_list = [Timestamp('2013-01-31', tz='Asia/Tokyo'),
                         Timestamp('2013-02-28', tz='Asia/Tokyo'),
                         Timestamp('2013-03-31', tz='Asia/Tokyo'),
                         Timestamp('2013-04-30', tz='Asia/Tokyo')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.asobject
        self.assertTrue(isinstance(result, Index))
        self.assertEqual(result.dtype, object)
        self.assert_index_equal(result, expected)
        self.assertEqual(result.name, expected.name)
        self.assertEqual(idx.tolist(), expected_list)

        idx = DatetimeIndex([datetime(2013, 1, 1), datetime(2013, 1, 2),
                             pd.NaT, datetime(2013, 1, 4)], name='idx')
        expected_list = [Timestamp('2013-01-01'),
                         Timestamp('2013-01-02'), pd.NaT,
                         Timestamp('2013-01-04')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.asobject
        self.assertTrue(isinstance(result, Index))
        self.assertEqual(result.dtype, object)
        self.assert_index_equal(result, expected)
        self.assertEqual(result.name, expected.name)
        self.assertEqual(idx.tolist(), expected_list)

    def test_minmax(self):
        for tz in self.tz:
            # monotonic
            idx1 = pd.DatetimeIndex(['2011-01-01', '2011-01-02',
                                     '2011-01-03'], tz=tz)
            self.assertTrue(idx1.is_monotonic)

            # non-monotonic
            idx2 = pd.DatetimeIndex(['2011-01-01', pd.NaT, '2011-01-03',
                                     '2011-01-02', pd.NaT], tz=tz)
            self.assertFalse(idx2.is_monotonic)

            for idx in [idx1, idx2]:
                self.assertEqual(idx.min(), Timestamp('2011-01-01', tz=tz))
                self.assertEqual(idx.max(), Timestamp('2011-01-03', tz=tz))
                self.assertEqual(idx.argmin(), 0)
                self.assertEqual(idx.argmax(), 2)

        for op in ['min', 'max']:
            # Return NaT
            obj = DatetimeIndex([])
            self.assertTrue(pd.isnull(getattr(obj, op)()))

            obj = DatetimeIndex([pd.NaT])
            self.assertTrue(pd.isnull(getattr(obj, op)()))

            obj = DatetimeIndex([pd.NaT, pd.NaT, pd.NaT])
            self.assertTrue(pd.isnull(getattr(obj, op)()))

    def test_numpy_minmax(self):
        dr = pd.date_range(start='2016-01-15', end='2016-01-20')

        self.assertEqual(np.min(dr),
                         Timestamp('2016-01-15 00:00:00', freq='D'))
        self.assertEqual(np.max(dr),
                         Timestamp('2016-01-20 00:00:00', freq='D'))

        errmsg = "the 'out' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, errmsg, np.min, dr, out=0)
        tm.assertRaisesRegexp(ValueError, errmsg, np.max, dr, out=0)

        self.assertEqual(np.argmin(dr), 0)
        self.assertEqual(np.argmax(dr), 5)

        if not _np_version_under1p10:
            errmsg = "the 'out' parameter is not supported"
            tm.assertRaisesRegexp(ValueError, errmsg, np.argmin, dr, out=0)
            tm.assertRaisesRegexp(ValueError, errmsg, np.argmax, dr, out=0)

    def test_round(self):
        for tz in self.tz:
            rng = pd.date_range(start='2016-01-01', periods=5,
                                freq='30Min', tz=tz)
            elt = rng[1]

            expected_rng = DatetimeIndex([
                Timestamp('2016-01-01 00:00:00', tz=tz, freq='30T'),
                Timestamp('2016-01-01 00:00:00', tz=tz, freq='30T'),
                Timestamp('2016-01-01 01:00:00', tz=tz, freq='30T'),
                Timestamp('2016-01-01 02:00:00', tz=tz, freq='30T'),
                Timestamp('2016-01-01 02:00:00', tz=tz, freq='30T'),
            ])
            expected_elt = expected_rng[1]

            tm.assert_index_equal(rng.round(freq='H'), expected_rng)
            self.assertEqual(elt.round(freq='H'), expected_elt)

            msg = pd.tseries.frequencies._INVALID_FREQ_ERROR
            with tm.assertRaisesRegexp(ValueError, msg):
                rng.round(freq='foo')
            with tm.assertRaisesRegexp(ValueError, msg):
                elt.round(freq='foo')

            msg = "<MonthEnd> is a non-fixed frequency"
            tm.assertRaisesRegexp(ValueError, msg, rng.round, freq='M')
            tm.assertRaisesRegexp(ValueError, msg, elt.round, freq='M')

    def test_repeat_range(self):
        rng = date_range('1/1/2000', '1/1/2001')

        result = rng.repeat(5)
        self.assertIsNone(result.freq)
        self.assertEqual(len(result), 5 * len(rng))

        for tz in self.tz:
            index = pd.date_range('2001-01-01', periods=2, freq='D', tz=tz)
            exp = pd.DatetimeIndex(['2001-01-01', '2001-01-01',
                                    '2001-01-02', '2001-01-02'], tz=tz)
            for res in [index.repeat(2), np.repeat(index, 2)]:
                tm.assert_index_equal(res, exp)
                self.assertIsNone(res.freq)

            index = pd.date_range('2001-01-01', periods=2, freq='2D', tz=tz)
            exp = pd.DatetimeIndex(['2001-01-01', '2001-01-01',
                                    '2001-01-03', '2001-01-03'], tz=tz)
            for res in [index.repeat(2), np.repeat(index, 2)]:
                tm.assert_index_equal(res, exp)
                self.assertIsNone(res.freq)

            index = pd.DatetimeIndex(['2001-01-01', 'NaT', '2003-01-01'],
                                     tz=tz)
            exp = pd.DatetimeIndex(['2001-01-01', '2001-01-01', '2001-01-01',
                                    'NaT', 'NaT', 'NaT',
                                    '2003-01-01', '2003-01-01', '2003-01-01'],
                                   tz=tz)
            for res in [index.repeat(3), np.repeat(index, 3)]:
                tm.assert_index_equal(res, exp)
                self.assertIsNone(res.freq)

    def test_repeat(self):
        reps = 2
        msg = "the 'axis' parameter is not supported"

        for tz in self.tz:
            rng = pd.date_range(start='2016-01-01', periods=2,
                                freq='30Min', tz=tz)

            expected_rng = DatetimeIndex([
                Timestamp('2016-01-01 00:00:00', tz=tz, freq='30T'),
                Timestamp('2016-01-01 00:00:00', tz=tz, freq='30T'),
                Timestamp('2016-01-01 00:30:00', tz=tz, freq='30T'),
                Timestamp('2016-01-01 00:30:00', tz=tz, freq='30T'),
            ])

            res = rng.repeat(reps)
            tm.assert_index_equal(res, expected_rng)
            self.assertIsNone(res.freq)

            tm.assert_index_equal(np.repeat(rng, reps), expected_rng)
            tm.assertRaisesRegexp(ValueError, msg, np.repeat,
                                  rng, reps, axis=1)

    def test_representation(self):

        idx = []
        idx.append(DatetimeIndex([], freq='D'))
        idx.append(DatetimeIndex(['2011-01-01'], freq='D'))
        idx.append(DatetimeIndex(['2011-01-01', '2011-01-02'], freq='D'))
        idx.append(DatetimeIndex(
            ['2011-01-01', '2011-01-02', '2011-01-03'], freq='D'))
        idx.append(DatetimeIndex(
            ['2011-01-01 09:00', '2011-01-01 10:00', '2011-01-01 11:00'
             ], freq='H', tz='Asia/Tokyo'))
        idx.append(DatetimeIndex(
            ['2011-01-01 09:00', '2011-01-01 10:00', pd.NaT], tz='US/Eastern'))
        idx.append(DatetimeIndex(
            ['2011-01-01 09:00', '2011-01-01 10:00', pd.NaT], tz='UTC'))

        exp = []
        exp.append("""DatetimeIndex([], dtype='datetime64[ns]', freq='D')""")
        exp.append("DatetimeIndex(['2011-01-01'], dtype='datetime64[ns]', "
                   "freq='D')")
        exp.append("DatetimeIndex(['2011-01-01', '2011-01-02'], "
                   "dtype='datetime64[ns]', freq='D')")
        exp.append("DatetimeIndex(['2011-01-01', '2011-01-02', '2011-01-03'], "
                   "dtype='datetime64[ns]', freq='D')")
        exp.append("DatetimeIndex(['2011-01-01 09:00:00+09:00', "
                   "'2011-01-01 10:00:00+09:00', '2011-01-01 11:00:00+09:00']"
                   ", dtype='datetime64[ns, Asia/Tokyo]', freq='H')")
        exp.append("DatetimeIndex(['2011-01-01 09:00:00-05:00', "
                   "'2011-01-01 10:00:00-05:00', 'NaT'], "
                   "dtype='datetime64[ns, US/Eastern]', freq=None)")
        exp.append("DatetimeIndex(['2011-01-01 09:00:00+00:00', "
                   "'2011-01-01 10:00:00+00:00', 'NaT'], "
                   "dtype='datetime64[ns, UTC]', freq=None)""")

        with pd.option_context('display.width', 300):
            for indx, expected in zip(idx, exp):
                for func in ['__repr__', '__unicode__', '__str__']:
                    result = getattr(indx, func)()
                    self.assertEqual(result, expected)

    def test_representation_to_series(self):
        idx1 = DatetimeIndex([], freq='D')
        idx2 = DatetimeIndex(['2011-01-01'], freq='D')
        idx3 = DatetimeIndex(['2011-01-01', '2011-01-02'], freq='D')
        idx4 = DatetimeIndex(
            ['2011-01-01', '2011-01-02', '2011-01-03'], freq='D')
        idx5 = DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00',
                              '2011-01-01 11:00'], freq='H', tz='Asia/Tokyo')
        idx6 = DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00', pd.NaT],
                             tz='US/Eastern')
        idx7 = DatetimeIndex(['2011-01-01 09:00', '2011-01-02 10:15'])

        exp1 = """Series([], dtype: datetime64[ns])"""

        exp2 = """0   2011-01-01
dtype: datetime64[ns]"""

        exp3 = """0   2011-01-01
1   2011-01-02
dtype: datetime64[ns]"""

        exp4 = """0   2011-01-01
1   2011-01-02
2   2011-01-03
dtype: datetime64[ns]"""

        exp5 = """0   2011-01-01 09:00:00+09:00
1   2011-01-01 10:00:00+09:00
2   2011-01-01 11:00:00+09:00
dtype: datetime64[ns, Asia/Tokyo]"""

        exp6 = """0   2011-01-01 09:00:00-05:00
1   2011-01-01 10:00:00-05:00
2                         NaT
dtype: datetime64[ns, US/Eastern]"""

        exp7 = """0   2011-01-01 09:00:00
1   2011-01-02 10:15:00
dtype: datetime64[ns]"""

        with pd.option_context('display.width', 300):
            for idx, expected in zip([idx1, idx2, idx3, idx4,
                                      idx5, idx6, idx7],
                                     [exp1, exp2, exp3, exp4,
                                      exp5, exp6, exp7]):
                result = repr(Series(idx))
                self.assertEqual(result, expected)

    def test_summary(self):
        # GH9116
        idx1 = DatetimeIndex([], freq='D')
        idx2 = DatetimeIndex(['2011-01-01'], freq='D')
        idx3 = DatetimeIndex(['2011-01-01', '2011-01-02'], freq='D')
        idx4 = DatetimeIndex(
            ['2011-01-01', '2011-01-02', '2011-01-03'], freq='D')
        idx5 = DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00',
                              '2011-01-01 11:00'],
                             freq='H', tz='Asia/Tokyo')
        idx6 = DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00', pd.NaT],
                             tz='US/Eastern')

        exp1 = """DatetimeIndex: 0 entries
Freq: D"""

        exp2 = """DatetimeIndex: 1 entries, 2011-01-01 to 2011-01-01
Freq: D"""

        exp3 = """DatetimeIndex: 2 entries, 2011-01-01 to 2011-01-02
Freq: D"""

        exp4 = """DatetimeIndex: 3 entries, 2011-01-01 to 2011-01-03
Freq: D"""

        exp5 = ("DatetimeIndex: 3 entries, 2011-01-01 09:00:00+09:00 "
                "to 2011-01-01 11:00:00+09:00\n"
                "Freq: H")

        exp6 = """DatetimeIndex: 3 entries, 2011-01-01 09:00:00-05:00 to NaT"""

        for idx, expected in zip([idx1, idx2, idx3, idx4, idx5, idx6],
                                 [exp1, exp2, exp3, exp4, exp5, exp6]):
            result = idx.summary()
            self.assertEqual(result, expected)

    def test_resolution(self):
        for freq, expected in zip(['A', 'Q', 'M', 'D', 'H', 'T',
                                   'S', 'L', 'U'],
                                  ['day', 'day', 'day', 'day', 'hour',
                                   'minute', 'second', 'millisecond',
                                   'microsecond']):
            for tz in self.tz:
                idx = pd.date_range(start='2013-04-01', periods=30, freq=freq,
                                    tz=tz)
                self.assertEqual(idx.resolution, expected)

    def test_union(self):
        for tz in self.tz:
            # union
            rng1 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)
            other1 = pd.date_range('1/6/2000', freq='D', periods=5, tz=tz)
            expected1 = pd.date_range('1/1/2000', freq='D', periods=10, tz=tz)

            rng2 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)
            other2 = pd.date_range('1/4/2000', freq='D', periods=5, tz=tz)
            expected2 = pd.date_range('1/1/2000', freq='D', periods=8, tz=tz)

            rng3 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)
            other3 = pd.DatetimeIndex([], tz=tz)
            expected3 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)

            for rng, other, expected in [(rng1, other1, expected1),
                                         (rng2, other2, expected2),
                                         (rng3, other3, expected3)]:

                result_union = rng.union(other)
                tm.assert_index_equal(result_union, expected)

    def test_add_iadd(self):
        for tz in self.tz:

            # offset
            offsets = [pd.offsets.Hour(2), timedelta(hours=2),
                       np.timedelta64(2, 'h'), Timedelta(hours=2)]

            for delta in offsets:
                rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
                result = rng + delta
                expected = pd.date_range('2000-01-01 02:00',
                                         '2000-02-01 02:00', tz=tz)
                tm.assert_index_equal(result, expected)
                rng += delta
                tm.assert_index_equal(rng, expected)

            # int
            rng = pd.date_range('2000-01-01 09:00', freq='H', periods=10,
                                tz=tz)
            result = rng + 1
            expected = pd.date_range('2000-01-01 10:00', freq='H', periods=10,
                                     tz=tz)
            tm.assert_index_equal(result, expected)
            rng += 1
            tm.assert_index_equal(rng, expected)

        idx = DatetimeIndex(['2011-01-01', '2011-01-02'])
        msg = "cannot add a datelike to a DatetimeIndex"
        with tm.assertRaisesRegexp(TypeError, msg):
            idx + Timestamp('2011-01-01')

        with tm.assertRaisesRegexp(TypeError, msg):
            Timestamp('2011-01-01') + idx

    def test_add_dti_dti(self):
        # previously performed setop (deprecated in 0.16.0), now raises
        # TypeError (GH14164)

        dti = date_range('20130101', periods=3)
        dti_tz = date_range('20130101', periods=3).tz_localize('US/Eastern')

        with tm.assertRaises(TypeError):
            dti + dti

        with tm.assertRaises(TypeError):
            dti_tz + dti_tz

        with tm.assertRaises(TypeError):
            dti_tz + dti

        with tm.assertRaises(TypeError):
            dti + dti_tz

    def test_difference(self):
        for tz in self.tz:
            # diff
            rng1 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)
            other1 = pd.date_range('1/6/2000', freq='D', periods=5, tz=tz)
            expected1 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)

            rng2 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)
            other2 = pd.date_range('1/4/2000', freq='D', periods=5, tz=tz)
            expected2 = pd.date_range('1/1/2000', freq='D', periods=3, tz=tz)

            rng3 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)
            other3 = pd.DatetimeIndex([], tz=tz)
            expected3 = pd.date_range('1/1/2000', freq='D', periods=5, tz=tz)

            for rng, other, expected in [(rng1, other1, expected1),
                                         (rng2, other2, expected2),
                                         (rng3, other3, expected3)]:
                result_diff = rng.difference(other)
                tm.assert_index_equal(result_diff, expected)

    def test_sub_isub(self):
        for tz in self.tz:

            # offset
            offsets = [pd.offsets.Hour(2), timedelta(hours=2),
                       np.timedelta64(2, 'h'), Timedelta(hours=2)]

            for delta in offsets:
                rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
                expected = pd.date_range('1999-12-31 22:00',
                                         '2000-01-31 22:00', tz=tz)

                result = rng - delta
                tm.assert_index_equal(result, expected)
                rng -= delta
                tm.assert_index_equal(rng, expected)

            # int
            rng = pd.date_range('2000-01-01 09:00', freq='H', periods=10,
                                tz=tz)
            result = rng - 1
            expected = pd.date_range('2000-01-01 08:00', freq='H', periods=10,
                                     tz=tz)
            tm.assert_index_equal(result, expected)
            rng -= 1
            tm.assert_index_equal(rng, expected)

    def test_sub_dti_dti(self):
        # previously performed setop (deprecated in 0.16.0), now changed to
        # return subtraction -> TimeDeltaIndex (GH ...)

        dti = date_range('20130101', periods=3)
        dti_tz = date_range('20130101', periods=3).tz_localize('US/Eastern')
        dti_tz2 = date_range('20130101', periods=3).tz_localize('UTC')
        expected = TimedeltaIndex([0, 0, 0])

        result = dti - dti
        tm.assert_index_equal(result, expected)

        result = dti_tz - dti_tz
        tm.assert_index_equal(result, expected)

        with tm.assertRaises(TypeError):
            dti_tz - dti

        with tm.assertRaises(TypeError):
            dti - dti_tz

        with tm.assertRaises(TypeError):
            dti_tz - dti_tz2

        # isub
        dti -= dti
        tm.assert_index_equal(dti, expected)

        # different length raises ValueError
        dti1 = date_range('20130101', periods=3)
        dti2 = date_range('20130101', periods=4)
        with tm.assertRaises(ValueError):
            dti1 - dti2

        # NaN propagation
        dti1 = DatetimeIndex(['2012-01-01', np.nan, '2012-01-03'])
        dti2 = DatetimeIndex(['2012-01-02', '2012-01-03', np.nan])
        expected = TimedeltaIndex(['1 days', np.nan, np.nan])
        result = dti2 - dti1
        tm.assert_index_equal(result, expected)

    def test_sub_period(self):
        # GH 13078
        # not supported, check TypeError
        p = pd.Period('2011-01-01', freq='D')

        for freq in [None, 'D']:
            idx = pd.DatetimeIndex(['2011-01-01', '2011-01-02'], freq=freq)

            with tm.assertRaises(TypeError):
                idx - p

            with tm.assertRaises(TypeError):
                p - idx

    def test_comp_nat(self):
        left = pd.DatetimeIndex([pd.Timestamp('2011-01-01'), pd.NaT,
                                 pd.Timestamp('2011-01-03')])
        right = pd.DatetimeIndex([pd.NaT, pd.NaT, pd.Timestamp('2011-01-03')])

        for l, r in [(left, right), (left.asobject, right.asobject)]:
            result = l == r
            expected = np.array([False, False, True])
            tm.assert_numpy_array_equal(result, expected)

            result = l != r
            expected = np.array([True, True, False])
            tm.assert_numpy_array_equal(result, expected)

            expected = np.array([False, False, False])
            tm.assert_numpy_array_equal(l == pd.NaT, expected)
            tm.assert_numpy_array_equal(pd.NaT == r, expected)

            expected = np.array([True, True, True])
            tm.assert_numpy_array_equal(l != pd.NaT, expected)
            tm.assert_numpy_array_equal(pd.NaT != l, expected)

            expected = np.array([False, False, False])
            tm.assert_numpy_array_equal(l < pd.NaT, expected)
            tm.assert_numpy_array_equal(pd.NaT > l, expected)

    def test_value_counts_unique(self):
        # GH 7735
        for tz in self.tz:
            idx = pd.date_range('2011-01-01 09:00', freq='H', periods=10)
            # create repeated values, 'n'th element is repeated by n+1 times
            idx = DatetimeIndex(np.repeat(idx.values, range(1, len(idx) + 1)),
                                tz=tz)

            exp_idx = pd.date_range('2011-01-01 18:00', freq='-1H', periods=10,
                                    tz=tz)
            expected = Series(range(10, 0, -1), index=exp_idx, dtype='int64')

            for obj in [idx, Series(idx)]:
                tm.assert_series_equal(obj.value_counts(), expected)

            expected = pd.date_range('2011-01-01 09:00', freq='H', periods=10,
                                     tz=tz)
            tm.assert_index_equal(idx.unique(), expected)

            idx = DatetimeIndex(['2013-01-01 09:00', '2013-01-01 09:00',
                                 '2013-01-01 09:00', '2013-01-01 08:00',
                                 '2013-01-01 08:00', pd.NaT], tz=tz)

            exp_idx = DatetimeIndex(['2013-01-01 09:00', '2013-01-01 08:00'],
                                    tz=tz)
            expected = Series([3, 2], index=exp_idx)

            for obj in [idx, Series(idx)]:
                tm.assert_series_equal(obj.value_counts(), expected)

            exp_idx = DatetimeIndex(['2013-01-01 09:00', '2013-01-01 08:00',
                                     pd.NaT], tz=tz)
            expected = Series([3, 2, 1], index=exp_idx)

            for obj in [idx, Series(idx)]:
                tm.assert_series_equal(obj.value_counts(dropna=False),
                                       expected)

            tm.assert_index_equal(idx.unique(), exp_idx)

    def test_nonunique_contains(self):
        # GH 9512
        for idx in map(DatetimeIndex,
                       ([0, 1, 0], [0, 0, -1], [0, -1, -1],
                        ['2015', '2015', '2016'], ['2015', '2015', '2014'])):
            tm.assertIn(idx[0], idx)

    def test_order(self):
        # with freq
        idx1 = DatetimeIndex(['2011-01-01', '2011-01-02',
                              '2011-01-03'], freq='D', name='idx')
        idx2 = DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00',
                              '2011-01-01 11:00'], freq='H',
                             tz='Asia/Tokyo', name='tzidx')

        for idx in [idx1, idx2]:
            ordered = idx.sort_values()
            self.assert_index_equal(ordered, idx)
            self.assertEqual(ordered.freq, idx.freq)

            ordered = idx.sort_values(ascending=False)
            expected = idx[::-1]
            self.assert_index_equal(ordered, expected)
            self.assertEqual(ordered.freq, expected.freq)
            self.assertEqual(ordered.freq.n, -1)

            ordered, indexer = idx.sort_values(return_indexer=True)
            self.assert_index_equal(ordered, idx)
            self.assert_numpy_array_equal(indexer,
                                          np.array([0, 1, 2]),
                                          check_dtype=False)
            self.assertEqual(ordered.freq, idx.freq)

            ordered, indexer = idx.sort_values(return_indexer=True,
                                               ascending=False)
            expected = idx[::-1]
            self.assert_index_equal(ordered, expected)
            self.assert_numpy_array_equal(indexer,
                                          np.array([2, 1, 0]),
                                          check_dtype=False)
            self.assertEqual(ordered.freq, expected.freq)
            self.assertEqual(ordered.freq.n, -1)

        # without freq
        for tz in self.tz:
            idx1 = DatetimeIndex(['2011-01-01', '2011-01-03', '2011-01-05',
                                  '2011-01-02', '2011-01-01'],
                                 tz=tz, name='idx1')
            exp1 = DatetimeIndex(['2011-01-01', '2011-01-01', '2011-01-02',
                                  '2011-01-03', '2011-01-05'],
                                 tz=tz, name='idx1')

            idx2 = DatetimeIndex(['2011-01-01', '2011-01-03', '2011-01-05',
                                  '2011-01-02', '2011-01-01'],
                                 tz=tz, name='idx2')

            exp2 = DatetimeIndex(['2011-01-01', '2011-01-01', '2011-01-02',
                                  '2011-01-03', '2011-01-05'],
                                 tz=tz, name='idx2')

            idx3 = DatetimeIndex([pd.NaT, '2011-01-03', '2011-01-05',
                                  '2011-01-02', pd.NaT], tz=tz, name='idx3')
            exp3 = DatetimeIndex([pd.NaT, pd.NaT, '2011-01-02', '2011-01-03',
                                  '2011-01-05'], tz=tz, name='idx3')

            for idx, expected in [(idx1, exp1), (idx2, exp2), (idx3, exp3)]:
                ordered = idx.sort_values()
                self.assert_index_equal(ordered, expected)
                self.assertIsNone(ordered.freq)

                ordered = idx.sort_values(ascending=False)
                self.assert_index_equal(ordered, expected[::-1])
                self.assertIsNone(ordered.freq)

                ordered, indexer = idx.sort_values(return_indexer=True)
                self.assert_index_equal(ordered, expected)

                exp = np.array([0, 4, 3, 1, 2])
                self.assert_numpy_array_equal(indexer, exp, check_dtype=False)
                self.assertIsNone(ordered.freq)

                ordered, indexer = idx.sort_values(return_indexer=True,
                                                   ascending=False)
                self.assert_index_equal(ordered, expected[::-1])

                exp = np.array([2, 1, 3, 4, 0])
                self.assert_numpy_array_equal(indexer, exp, check_dtype=False)
                self.assertIsNone(ordered.freq)

    def test_getitem(self):
        idx1 = pd.date_range('2011-01-01', '2011-01-31', freq='D', name='idx')
        idx2 = pd.date_range('2011-01-01', '2011-01-31', freq='D',
                             tz='Asia/Tokyo', name='idx')

        for idx in [idx1, idx2]:
            result = idx[0]
            self.assertEqual(result, Timestamp('2011-01-01', tz=idx.tz))

            result = idx[0:5]
            expected = pd.date_range('2011-01-01', '2011-01-05', freq='D',
                                     tz=idx.tz, name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)

            result = idx[0:10:2]
            expected = pd.date_range('2011-01-01', '2011-01-09', freq='2D',
                                     tz=idx.tz, name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)

            result = idx[-20:-5:3]
            expected = pd.date_range('2011-01-12', '2011-01-24', freq='3D',
                                     tz=idx.tz, name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)

            result = idx[4::-1]
            expected = DatetimeIndex(['2011-01-05', '2011-01-04', '2011-01-03',
                                      '2011-01-02', '2011-01-01'],
                                     freq='-1D', tz=idx.tz, name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)

    def test_drop_duplicates_metadata(self):
        # GH 10115
        idx = pd.date_range('2011-01-01', '2011-01-31', freq='D', name='idx')
        result = idx.drop_duplicates()
        self.assert_index_equal(idx, result)
        self.assertEqual(idx.freq, result.freq)

        idx_dup = idx.append(idx)
        self.assertIsNone(idx_dup.freq)  # freq is reset
        result = idx_dup.drop_duplicates()
        self.assert_index_equal(idx, result)
        self.assertIsNone(result.freq)

    def test_drop_duplicates(self):
        # to check Index/Series compat
        base = pd.date_range('2011-01-01', '2011-01-31', freq='D', name='idx')
        idx = base.append(base[:5])

        res = idx.drop_duplicates()
        tm.assert_index_equal(res, base)
        res = Series(idx).drop_duplicates()
        tm.assert_series_equal(res, Series(base))

        res = idx.drop_duplicates(keep='last')
        exp = base[5:].append(base[:5])
        tm.assert_index_equal(res, exp)
        res = Series(idx).drop_duplicates(keep='last')
        tm.assert_series_equal(res, Series(exp, index=np.arange(5, 36)))

        res = idx.drop_duplicates(keep=False)
        tm.assert_index_equal(res, base[5:])
        res = Series(idx).drop_duplicates(keep=False)
        tm.assert_series_equal(res, Series(base[5:], index=np.arange(5, 31)))

    def test_take(self):
        # GH 10295
        idx1 = pd.date_range('2011-01-01', '2011-01-31', freq='D', name='idx')
        idx2 = pd.date_range('2011-01-01', '2011-01-31', freq='D',
                             tz='Asia/Tokyo', name='idx')

        for idx in [idx1, idx2]:
            result = idx.take([0])
            self.assertEqual(result, Timestamp('2011-01-01', tz=idx.tz))

            result = idx.take([0, 1, 2])
            expected = pd.date_range('2011-01-01', '2011-01-03', freq='D',
                                     tz=idx.tz, name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)

            result = idx.take([0, 2, 4])
            expected = pd.date_range('2011-01-01', '2011-01-05', freq='2D',
                                     tz=idx.tz, name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)

            result = idx.take([7, 4, 1])
            expected = pd.date_range('2011-01-08', '2011-01-02', freq='-3D',
                                     tz=idx.tz, name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)

            result = idx.take([3, 2, 5])
            expected = DatetimeIndex(['2011-01-04', '2011-01-03',
                                      '2011-01-06'],
                                     freq=None, tz=idx.tz, name='idx')
            self.assert_index_equal(result, expected)
            self.assertIsNone(result.freq)

            result = idx.take([-3, 2, 5])
            expected = DatetimeIndex(['2011-01-29', '2011-01-03',
                                      '2011-01-06'],
                                     freq=None, tz=idx.tz, name='idx')
            self.assert_index_equal(result, expected)
            self.assertIsNone(result.freq)

    def test_take_invalid_kwargs(self):
        idx = pd.date_range('2011-01-01', '2011-01-31', freq='D', name='idx')
        indices = [1, 6, 5, 9, 10, 13, 15, 3]

        msg = r"take\(\) got an unexpected keyword argument 'foo'"
        tm.assertRaisesRegexp(TypeError, msg, idx.take,
                              indices, foo=2)

        msg = "the 'out' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, idx.take,
                              indices, out=indices)

        msg = "the 'mode' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, idx.take,
                              indices, mode='clip')

    def test_infer_freq(self):
        # GH 11018
        for freq in ['A', '2A', '-2A', 'Q', '-1Q', 'M', '-1M', 'D', '3D',
                     '-3D', 'W', '-1W', 'H', '2H', '-2H', 'T', '2T', 'S',
                     '-3S']:
            idx = pd.date_range('2011-01-01 09:00:00', freq=freq, periods=10)
            result = pd.DatetimeIndex(idx.asi8, freq='infer')
            tm.assert_index_equal(idx, result)
            self.assertEqual(result.freq, freq)

    def test_nat_new(self):
        idx = pd.date_range('2011-01-01', freq='D', periods=5, name='x')
        result = idx._nat_new()
        exp = pd.DatetimeIndex([pd.NaT] * 5, name='x')
        tm.assert_index_equal(result, exp)

        result = idx._nat_new(box=False)
        exp = np.array([tslib.iNaT] * 5, dtype=np.int64)
        tm.assert_numpy_array_equal(result, exp)

    def test_shift(self):
        # GH 9903
        for tz in self.tz:
            idx = pd.DatetimeIndex([], name='xxx', tz=tz)
            tm.assert_index_equal(idx.shift(0, freq='H'), idx)
            tm.assert_index_equal(idx.shift(3, freq='H'), idx)

            idx = pd.DatetimeIndex(['2011-01-01 10:00', '2011-01-01 11:00'
                                    '2011-01-01 12:00'], name='xxx', tz=tz)
            tm.assert_index_equal(idx.shift(0, freq='H'), idx)
            exp = pd.DatetimeIndex(['2011-01-01 13:00', '2011-01-01 14:00'
                                    '2011-01-01 15:00'], name='xxx', tz=tz)
            tm.assert_index_equal(idx.shift(3, freq='H'), exp)
            exp = pd.DatetimeIndex(['2011-01-01 07:00', '2011-01-01 08:00'
                                    '2011-01-01 09:00'], name='xxx', tz=tz)
            tm.assert_index_equal(idx.shift(-3, freq='H'), exp)

    def test_nat(self):
        self.assertIs(pd.DatetimeIndex._na_value, pd.NaT)
        self.assertIs(pd.DatetimeIndex([])._na_value, pd.NaT)

        for tz in [None, 'US/Eastern', 'UTC']:
            idx = pd.DatetimeIndex(['2011-01-01', '2011-01-02'], tz=tz)
            self.assertTrue(idx._can_hold_na)

            tm.assert_numpy_array_equal(idx._isnan, np.array([False, False]))
            self.assertFalse(idx.hasnans)
            tm.assert_numpy_array_equal(idx._nan_idxs,
                                        np.array([], dtype=np.intp))

            idx = pd.DatetimeIndex(['2011-01-01', 'NaT'], tz=tz)
            self.assertTrue(idx._can_hold_na)

            tm.assert_numpy_array_equal(idx._isnan, np.array([False, True]))
            self.assertTrue(idx.hasnans)
            tm.assert_numpy_array_equal(idx._nan_idxs,
                                        np.array([1], dtype=np.intp))

    def test_equals(self):
        # GH 13107
        for tz in [None, 'UTC', 'US/Eastern', 'Asia/Tokyo']:
            idx = pd.DatetimeIndex(['2011-01-01', '2011-01-02', 'NaT'])
            self.assertTrue(idx.equals(idx))
            self.assertTrue(idx.equals(idx.copy()))
            self.assertTrue(idx.equals(idx.asobject))
            self.assertTrue(idx.asobject.equals(idx))
            self.assertTrue(idx.asobject.equals(idx.asobject))
            self.assertFalse(idx.equals(list(idx)))
            self.assertFalse(idx.equals(pd.Series(idx)))

            idx2 = pd.DatetimeIndex(['2011-01-01', '2011-01-02', 'NaT'],
                                    tz='US/Pacific')
            self.assertFalse(idx.equals(idx2))
            self.assertFalse(idx.equals(idx2.copy()))
            self.assertFalse(idx.equals(idx2.asobject))
            self.assertFalse(idx.asobject.equals(idx2))
            self.assertFalse(idx.equals(list(idx2)))
            self.assertFalse(idx.equals(pd.Series(idx2)))

            # same internal, different tz
            idx3 = pd.DatetimeIndex._simple_new(idx.asi8, tz='US/Pacific')
            tm.assert_numpy_array_equal(idx.asi8, idx3.asi8)
            self.assertFalse(idx.equals(idx3))
            self.assertFalse(idx.equals(idx3.copy()))
            self.assertFalse(idx.equals(idx3.asobject))
            self.assertFalse(idx.asobject.equals(idx3))
            self.assertFalse(idx.equals(list(idx3)))
            self.assertFalse(idx.equals(pd.Series(idx3)))


class TestDateTimeIndexToJulianDate(tm.TestCase):

    def test_1700(self):
        r1 = Float64Index([2345897.5, 2345898.5, 2345899.5, 2345900.5,
                           2345901.5])
        r2 = date_range(start=Timestamp('1710-10-01'), periods=5,
                        freq='D').to_julian_date()
        self.assertIsInstance(r2, Float64Index)
        tm.assert_index_equal(r1, r2)

    def test_2000(self):
        r1 = Float64Index([2451601.5, 2451602.5, 2451603.5, 2451604.5,
                           2451605.5])
        r2 = date_range(start=Timestamp('2000-02-27'), periods=5,
                        freq='D').to_julian_date()
        self.assertIsInstance(r2, Float64Index)
        tm.assert_index_equal(r1, r2)

    def test_hour(self):
        r1 = Float64Index(
            [2451601.5, 2451601.5416666666666666, 2451601.5833333333333333,
             2451601.625, 2451601.6666666666666666])
        r2 = date_range(start=Timestamp('2000-02-27'), periods=5,
                        freq='H').to_julian_date()
        self.assertIsInstance(r2, Float64Index)
        tm.assert_index_equal(r1, r2)

    def test_minute(self):
        r1 = Float64Index(
            [2451601.5, 2451601.5006944444444444, 2451601.5013888888888888,
             2451601.5020833333333333, 2451601.5027777777777777])
        r2 = date_range(start=Timestamp('2000-02-27'), periods=5,
                        freq='T').to_julian_date()
        self.assertIsInstance(r2, Float64Index)
        tm.assert_index_equal(r1, r2)

    def test_second(self):
        r1 = Float64Index(
            [2451601.5, 2451601.500011574074074, 2451601.5000231481481481,
             2451601.5000347222222222, 2451601.5000462962962962])
        r2 = date_range(start=Timestamp('2000-02-27'), periods=5,
                        freq='S').to_julian_date()
        self.assertIsInstance(r2, Float64Index)
        tm.assert_index_equal(r1, r2)
