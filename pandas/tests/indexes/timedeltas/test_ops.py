import numpy as np
from datetime import timedelta
from distutils.version import LooseVersion

import pandas as pd
import pandas.util.testing as tm
from pandas import to_timedelta
from pandas.util.testing import assert_series_equal, assert_frame_equal
from pandas import (Series, Timedelta, DataFrame, Timestamp, TimedeltaIndex,
                    timedelta_range, date_range, DatetimeIndex, Int64Index,
                    _np_version_under1p10, Float64Index, Index, tslib)

from pandas.tests.test_base import Ops


class TestTimedeltaIndexOps(Ops):
    def setUp(self):
        super(TestTimedeltaIndexOps, self).setUp()
        mask = lambda x: isinstance(x, TimedeltaIndex)
        self.is_valid_objs = [o for o in self.objs if mask(o)]
        self.not_valid_objs = []

    def test_ops_properties(self):
        self.check_ops_properties(['days', 'hours', 'minutes', 'seconds',
                                   'milliseconds'])
        self.check_ops_properties(['microseconds', 'nanoseconds'])

    def test_asobject_tolist(self):
        idx = timedelta_range(start='1 days', periods=4, freq='D', name='idx')
        expected_list = [Timedelta('1 days'), Timedelta('2 days'),
                         Timedelta('3 days'), Timedelta('4 days')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.asobject
        self.assertTrue(isinstance(result, Index))

        self.assertEqual(result.dtype, object)
        self.assert_index_equal(result, expected)
        self.assertEqual(result.name, expected.name)
        self.assertEqual(idx.tolist(), expected_list)

        idx = TimedeltaIndex([timedelta(days=1), timedelta(days=2), pd.NaT,
                              timedelta(days=4)], name='idx')
        expected_list = [Timedelta('1 days'), Timedelta('2 days'), pd.NaT,
                         Timedelta('4 days')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.asobject
        self.assertTrue(isinstance(result, Index))
        self.assertEqual(result.dtype, object)
        self.assert_index_equal(result, expected)
        self.assertEqual(result.name, expected.name)
        self.assertEqual(idx.tolist(), expected_list)

    def test_minmax(self):

        # monotonic
        idx1 = TimedeltaIndex(['1 days', '2 days', '3 days'])
        self.assertTrue(idx1.is_monotonic)

        # non-monotonic
        idx2 = TimedeltaIndex(['1 days', np.nan, '3 days', 'NaT'])
        self.assertFalse(idx2.is_monotonic)

        for idx in [idx1, idx2]:
            self.assertEqual(idx.min(), Timedelta('1 days')),
            self.assertEqual(idx.max(), Timedelta('3 days')),
            self.assertEqual(idx.argmin(), 0)
            self.assertEqual(idx.argmax(), 2)

        for op in ['min', 'max']:
            # Return NaT
            obj = TimedeltaIndex([])
            self.assertTrue(pd.isnull(getattr(obj, op)()))

            obj = TimedeltaIndex([pd.NaT])
            self.assertTrue(pd.isnull(getattr(obj, op)()))

            obj = TimedeltaIndex([pd.NaT, pd.NaT, pd.NaT])
            self.assertTrue(pd.isnull(getattr(obj, op)()))

    def test_numpy_minmax(self):
        dr = pd.date_range(start='2016-01-15', end='2016-01-20')
        td = TimedeltaIndex(np.asarray(dr))

        self.assertEqual(np.min(td), Timedelta('16815 days'))
        self.assertEqual(np.max(td), Timedelta('16820 days'))

        errmsg = "the 'out' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, errmsg, np.min, td, out=0)
        tm.assertRaisesRegexp(ValueError, errmsg, np.max, td, out=0)

        self.assertEqual(np.argmin(td), 0)
        self.assertEqual(np.argmax(td), 5)

        if not _np_version_under1p10:
            errmsg = "the 'out' parameter is not supported"
            tm.assertRaisesRegexp(ValueError, errmsg, np.argmin, td, out=0)
            tm.assertRaisesRegexp(ValueError, errmsg, np.argmax, td, out=0)

    def test_round(self):
        td = pd.timedelta_range(start='16801 days', periods=5, freq='30Min')
        elt = td[1]

        expected_rng = TimedeltaIndex([
            Timedelta('16801 days 00:00:00'),
            Timedelta('16801 days 00:00:00'),
            Timedelta('16801 days 01:00:00'),
            Timedelta('16801 days 02:00:00'),
            Timedelta('16801 days 02:00:00'),
        ])
        expected_elt = expected_rng[1]

        tm.assert_index_equal(td.round(freq='H'), expected_rng)
        self.assertEqual(elt.round(freq='H'), expected_elt)

        msg = pd.tseries.frequencies._INVALID_FREQ_ERROR
        with self.assertRaisesRegexp(ValueError, msg):
            td.round(freq='foo')
        with tm.assertRaisesRegexp(ValueError, msg):
            elt.round(freq='foo')

        msg = "<MonthEnd> is a non-fixed frequency"
        tm.assertRaisesRegexp(ValueError, msg, td.round, freq='M')
        tm.assertRaisesRegexp(ValueError, msg, elt.round, freq='M')

    def test_representation(self):
        idx1 = TimedeltaIndex([], freq='D')
        idx2 = TimedeltaIndex(['1 days'], freq='D')
        idx3 = TimedeltaIndex(['1 days', '2 days'], freq='D')
        idx4 = TimedeltaIndex(['1 days', '2 days', '3 days'], freq='D')
        idx5 = TimedeltaIndex(['1 days 00:00:01', '2 days', '3 days'])

        exp1 = """TimedeltaIndex([], dtype='timedelta64[ns]', freq='D')"""

        exp2 = ("TimedeltaIndex(['1 days'], dtype='timedelta64[ns]', "
                "freq='D')")

        exp3 = ("TimedeltaIndex(['1 days', '2 days'], "
                "dtype='timedelta64[ns]', freq='D')")

        exp4 = ("TimedeltaIndex(['1 days', '2 days', '3 days'], "
                "dtype='timedelta64[ns]', freq='D')")

        exp5 = ("TimedeltaIndex(['1 days 00:00:01', '2 days 00:00:00', "
                "'3 days 00:00:00'], dtype='timedelta64[ns]', freq=None)")

        with pd.option_context('display.width', 300):
            for idx, expected in zip([idx1, idx2, idx3, idx4, idx5],
                                     [exp1, exp2, exp3, exp4, exp5]):
                for func in ['__repr__', '__unicode__', '__str__']:
                    result = getattr(idx, func)()
                    self.assertEqual(result, expected)

    def test_representation_to_series(self):
        idx1 = TimedeltaIndex([], freq='D')
        idx2 = TimedeltaIndex(['1 days'], freq='D')
        idx3 = TimedeltaIndex(['1 days', '2 days'], freq='D')
        idx4 = TimedeltaIndex(['1 days', '2 days', '3 days'], freq='D')
        idx5 = TimedeltaIndex(['1 days 00:00:01', '2 days', '3 days'])

        exp1 = """Series([], dtype: timedelta64[ns])"""

        exp2 = """0   1 days
dtype: timedelta64[ns]"""

        exp3 = """0   1 days
1   2 days
dtype: timedelta64[ns]"""

        exp4 = """0   1 days
1   2 days
2   3 days
dtype: timedelta64[ns]"""

        exp5 = """0   1 days 00:00:01
1   2 days 00:00:00
2   3 days 00:00:00
dtype: timedelta64[ns]"""

        with pd.option_context('display.width', 300):
            for idx, expected in zip([idx1, idx2, idx3, idx4, idx5],
                                     [exp1, exp2, exp3, exp4, exp5]):
                result = repr(pd.Series(idx))
                self.assertEqual(result, expected)

    def test_summary(self):
        # GH9116
        idx1 = TimedeltaIndex([], freq='D')
        idx2 = TimedeltaIndex(['1 days'], freq='D')
        idx3 = TimedeltaIndex(['1 days', '2 days'], freq='D')
        idx4 = TimedeltaIndex(['1 days', '2 days', '3 days'], freq='D')
        idx5 = TimedeltaIndex(['1 days 00:00:01', '2 days', '3 days'])

        exp1 = """TimedeltaIndex: 0 entries
Freq: D"""

        exp2 = """TimedeltaIndex: 1 entries, 1 days to 1 days
Freq: D"""

        exp3 = """TimedeltaIndex: 2 entries, 1 days to 2 days
Freq: D"""

        exp4 = """TimedeltaIndex: 3 entries, 1 days to 3 days
Freq: D"""

        exp5 = ("TimedeltaIndex: 3 entries, 1 days 00:00:01 to 3 days "
                "00:00:00")

        for idx, expected in zip([idx1, idx2, idx3, idx4, idx5],
                                 [exp1, exp2, exp3, exp4, exp5]):
            result = idx.summary()
            self.assertEqual(result, expected)

    def test_add_iadd(self):

        # only test adding/sub offsets as + is now numeric

        # offset
        offsets = [pd.offsets.Hour(2), timedelta(hours=2),
                   np.timedelta64(2, 'h'), Timedelta(hours=2)]

        for delta in offsets:
            rng = timedelta_range('1 days', '10 days')
            result = rng + delta
            expected = timedelta_range('1 days 02:00:00', '10 days 02:00:00',
                                       freq='D')
            tm.assert_index_equal(result, expected)
            rng += delta
            tm.assert_index_equal(rng, expected)

        # int
        rng = timedelta_range('1 days 09:00:00', freq='H', periods=10)
        result = rng + 1
        expected = timedelta_range('1 days 10:00:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)
        rng += 1
        tm.assert_index_equal(rng, expected)

    def test_sub_isub(self):
        # only test adding/sub offsets as - is now numeric

        # offset
        offsets = [pd.offsets.Hour(2), timedelta(hours=2),
                   np.timedelta64(2, 'h'), Timedelta(hours=2)]

        for delta in offsets:
            rng = timedelta_range('1 days', '10 days')
            result = rng - delta
            expected = timedelta_range('0 days 22:00:00', '9 days 22:00:00')
            tm.assert_index_equal(result, expected)
            rng -= delta
            tm.assert_index_equal(rng, expected)

        # int
        rng = timedelta_range('1 days 09:00:00', freq='H', periods=10)
        result = rng - 1
        expected = timedelta_range('1 days 08:00:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)
        rng -= 1
        tm.assert_index_equal(rng, expected)

        idx = TimedeltaIndex(['1 day', '2 day'])
        msg = "cannot subtract a datelike from a TimedeltaIndex"
        with tm.assertRaisesRegexp(TypeError, msg):
            idx - Timestamp('2011-01-01')

        result = Timestamp('2011-01-01') + idx
        expected = DatetimeIndex(['2011-01-02', '2011-01-03'])
        tm.assert_index_equal(result, expected)

    def test_ops_compat(self):

        offsets = [pd.offsets.Hour(2), timedelta(hours=2),
                   np.timedelta64(2, 'h'), Timedelta(hours=2)]

        rng = timedelta_range('1 days', '10 days', name='foo')

        # multiply
        for offset in offsets:
            self.assertRaises(TypeError, lambda: rng * offset)

        # divide
        expected = Int64Index((np.arange(10) + 1) * 12, name='foo')
        for offset in offsets:
            result = rng / offset
            tm.assert_index_equal(result, expected, exact=False)

        # divide with nats
        rng = TimedeltaIndex(['1 days', pd.NaT, '2 days'], name='foo')
        expected = Float64Index([12, np.nan, 24], name='foo')
        for offset in offsets:
            result = rng / offset
            tm.assert_index_equal(result, expected)

        # don't allow division by NaT (make could in the future)
        self.assertRaises(TypeError, lambda: rng / pd.NaT)

    def test_subtraction_ops(self):

        # with datetimes/timedelta and tdi/dti
        tdi = TimedeltaIndex(['1 days', pd.NaT, '2 days'], name='foo')
        dti = date_range('20130101', periods=3, name='bar')
        td = Timedelta('1 days')
        dt = Timestamp('20130101')

        self.assertRaises(TypeError, lambda: tdi - dt)
        self.assertRaises(TypeError, lambda: tdi - dti)
        self.assertRaises(TypeError, lambda: td - dt)
        self.assertRaises(TypeError, lambda: td - dti)

        result = dt - dti
        expected = TimedeltaIndex(['0 days', '-1 days', '-2 days'], name='bar')
        tm.assert_index_equal(result, expected)

        result = dti - dt
        expected = TimedeltaIndex(['0 days', '1 days', '2 days'], name='bar')
        tm.assert_index_equal(result, expected)

        result = tdi - td
        expected = TimedeltaIndex(['0 days', pd.NaT, '1 days'], name='foo')
        tm.assert_index_equal(result, expected, check_names=False)

        result = td - tdi
        expected = TimedeltaIndex(['0 days', pd.NaT, '-1 days'], name='foo')
        tm.assert_index_equal(result, expected, check_names=False)

        result = dti - td
        expected = DatetimeIndex(
            ['20121231', '20130101', '20130102'], name='bar')
        tm.assert_index_equal(result, expected, check_names=False)

        result = dt - tdi
        expected = DatetimeIndex(['20121231', pd.NaT, '20121230'], name='foo')
        tm.assert_index_equal(result, expected)

    def test_subtraction_ops_with_tz(self):

        # check that dt/dti subtraction ops with tz are validated
        dti = date_range('20130101', periods=3)
        ts = Timestamp('20130101')
        dt = ts.to_pydatetime()
        dti_tz = date_range('20130101', periods=3).tz_localize('US/Eastern')
        ts_tz = Timestamp('20130101').tz_localize('US/Eastern')
        ts_tz2 = Timestamp('20130101').tz_localize('CET')
        dt_tz = ts_tz.to_pydatetime()
        td = Timedelta('1 days')

        def _check(result, expected):
            self.assertEqual(result, expected)
            self.assertIsInstance(result, Timedelta)

        # scalars
        result = ts - ts
        expected = Timedelta('0 days')
        _check(result, expected)

        result = dt_tz - ts_tz
        expected = Timedelta('0 days')
        _check(result, expected)

        result = ts_tz - dt_tz
        expected = Timedelta('0 days')
        _check(result, expected)

        # tz mismatches
        self.assertRaises(TypeError, lambda: dt_tz - ts)
        self.assertRaises(TypeError, lambda: dt_tz - dt)
        self.assertRaises(TypeError, lambda: dt_tz - ts_tz2)
        self.assertRaises(TypeError, lambda: dt - dt_tz)
        self.assertRaises(TypeError, lambda: ts - dt_tz)
        self.assertRaises(TypeError, lambda: ts_tz2 - ts)
        self.assertRaises(TypeError, lambda: ts_tz2 - dt)
        self.assertRaises(TypeError, lambda: ts_tz - ts_tz2)

        # with dti
        self.assertRaises(TypeError, lambda: dti - ts_tz)
        self.assertRaises(TypeError, lambda: dti_tz - ts)
        self.assertRaises(TypeError, lambda: dti_tz - ts_tz2)

        result = dti_tz - dt_tz
        expected = TimedeltaIndex(['0 days', '1 days', '2 days'])
        tm.assert_index_equal(result, expected)

        result = dt_tz - dti_tz
        expected = TimedeltaIndex(['0 days', '-1 days', '-2 days'])
        tm.assert_index_equal(result, expected)

        result = dti_tz - ts_tz
        expected = TimedeltaIndex(['0 days', '1 days', '2 days'])
        tm.assert_index_equal(result, expected)

        result = ts_tz - dti_tz
        expected = TimedeltaIndex(['0 days', '-1 days', '-2 days'])
        tm.assert_index_equal(result, expected)

        result = td - td
        expected = Timedelta('0 days')
        _check(result, expected)

        result = dti_tz - td
        expected = DatetimeIndex(
            ['20121231', '20130101', '20130102'], tz='US/Eastern')
        tm.assert_index_equal(result, expected)

    def test_dti_tdi_numeric_ops(self):

        # These are normally union/diff set-like ops
        tdi = TimedeltaIndex(['1 days', pd.NaT, '2 days'], name='foo')
        dti = date_range('20130101', periods=3, name='bar')

        # TODO(wesm): unused?
        # td = Timedelta('1 days')
        # dt = Timestamp('20130101')

        result = tdi - tdi
        expected = TimedeltaIndex(['0 days', pd.NaT, '0 days'], name='foo')
        tm.assert_index_equal(result, expected)

        result = tdi + tdi
        expected = TimedeltaIndex(['2 days', pd.NaT, '4 days'], name='foo')
        tm.assert_index_equal(result, expected)

        result = dti - tdi  # name will be reset
        expected = DatetimeIndex(['20121231', pd.NaT, '20130101'])
        tm.assert_index_equal(result, expected)

    def test_sub_period(self):
        # GH 13078
        # not supported, check TypeError
        p = pd.Period('2011-01-01', freq='D')

        for freq in [None, 'H']:
            idx = pd.TimedeltaIndex(['1 hours', '2 hours'], freq=freq)

            with tm.assertRaises(TypeError):
                idx - p

            with tm.assertRaises(TypeError):
                p - idx

    def test_addition_ops(self):

        # with datetimes/timedelta and tdi/dti
        tdi = TimedeltaIndex(['1 days', pd.NaT, '2 days'], name='foo')
        dti = date_range('20130101', periods=3, name='bar')
        td = Timedelta('1 days')
        dt = Timestamp('20130101')

        result = tdi + dt
        expected = DatetimeIndex(['20130102', pd.NaT, '20130103'], name='foo')
        tm.assert_index_equal(result, expected)

        result = dt + tdi
        expected = DatetimeIndex(['20130102', pd.NaT, '20130103'], name='foo')
        tm.assert_index_equal(result, expected)

        result = td + tdi
        expected = TimedeltaIndex(['2 days', pd.NaT, '3 days'], name='foo')
        tm.assert_index_equal(result, expected)

        result = tdi + td
        expected = TimedeltaIndex(['2 days', pd.NaT, '3 days'], name='foo')
        tm.assert_index_equal(result, expected)

        # unequal length
        self.assertRaises(ValueError, lambda: tdi + dti[0:1])
        self.assertRaises(ValueError, lambda: tdi[0:1] + dti)

        # random indexes
        self.assertRaises(TypeError, lambda: tdi + Int64Index([1, 2, 3]))

        # this is a union!
        # self.assertRaises(TypeError, lambda : Int64Index([1,2,3]) + tdi)

        result = tdi + dti  # name will be reset
        expected = DatetimeIndex(['20130102', pd.NaT, '20130105'])
        tm.assert_index_equal(result, expected)

        result = dti + tdi  # name will be reset
        expected = DatetimeIndex(['20130102', pd.NaT, '20130105'])
        tm.assert_index_equal(result, expected)

        result = dt + td
        expected = Timestamp('20130102')
        self.assertEqual(result, expected)

        result = td + dt
        expected = Timestamp('20130102')
        self.assertEqual(result, expected)

    def test_comp_nat(self):
        left = pd.TimedeltaIndex([pd.Timedelta('1 days'), pd.NaT,
                                 pd.Timedelta('3 days')])
        right = pd.TimedeltaIndex([pd.NaT, pd.NaT, pd.Timedelta('3 days')])

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

        idx = timedelta_range('1 days 09:00:00', freq='H', periods=10)
        # create repeated values, 'n'th element is repeated by n+1 times
        idx = TimedeltaIndex(np.repeat(idx.values, range(1, len(idx) + 1)))

        exp_idx = timedelta_range('1 days 18:00:00', freq='-1H', periods=10)
        expected = Series(range(10, 0, -1), index=exp_idx, dtype='int64')

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(), expected)

        expected = timedelta_range('1 days 09:00:00', freq='H', periods=10)
        tm.assert_index_equal(idx.unique(), expected)

        idx = TimedeltaIndex(['1 days 09:00:00', '1 days 09:00:00',
                              '1 days 09:00:00', '1 days 08:00:00',
                              '1 days 08:00:00', pd.NaT])

        exp_idx = TimedeltaIndex(['1 days 09:00:00', '1 days 08:00:00'])
        expected = Series([3, 2], index=exp_idx)

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(), expected)

        exp_idx = TimedeltaIndex(['1 days 09:00:00', '1 days 08:00:00',
                                  pd.NaT])
        expected = Series([3, 2, 1], index=exp_idx)

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(dropna=False), expected)

        tm.assert_index_equal(idx.unique(), exp_idx)

    def test_nonunique_contains(self):
        # GH 9512
        for idx in map(TimedeltaIndex, ([0, 1, 0], [0, 0, -1], [0, -1, -1],
                                        ['00:01:00', '00:01:00', '00:02:00'],
                                        ['00:01:00', '00:01:00', '00:00:01'])):
            tm.assertIn(idx[0], idx)

    def test_unknown_attribute(self):
        # GH 9680
        tdi = pd.timedelta_range(start=0, periods=10, freq='1s')
        ts = pd.Series(np.random.normal(size=10), index=tdi)
        self.assertNotIn('foo', ts.__dict__.keys())
        self.assertRaises(AttributeError, lambda: ts.foo)

    def test_order(self):
        # GH 10295
        idx1 = TimedeltaIndex(['1 day', '2 day', '3 day'], freq='D',
                              name='idx')
        idx2 = TimedeltaIndex(
            ['1 hour', '2 hour', '3 hour'], freq='H', name='idx')

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
            self.assert_index_equal(ordered, idx[::-1])
            self.assertEqual(ordered.freq, expected.freq)
            self.assertEqual(ordered.freq.n, -1)

        idx1 = TimedeltaIndex(['1 hour', '3 hour', '5 hour',
                               '2 hour ', '1 hour'], name='idx1')
        exp1 = TimedeltaIndex(['1 hour', '1 hour', '2 hour',
                               '3 hour', '5 hour'], name='idx1')

        idx2 = TimedeltaIndex(['1 day', '3 day', '5 day',
                               '2 day', '1 day'], name='idx2')

        # TODO(wesm): unused?
        # exp2 = TimedeltaIndex(['1 day', '1 day', '2 day',
        #                        '3 day', '5 day'], name='idx2')

        # idx3 = TimedeltaIndex([pd.NaT, '3 minute', '5 minute',
        #                        '2 minute', pd.NaT], name='idx3')
        # exp3 = TimedeltaIndex([pd.NaT, pd.NaT, '2 minute', '3 minute',
        #                        '5 minute'], name='idx3')

        for idx, expected in [(idx1, exp1), (idx1, exp1), (idx1, exp1)]:
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
        idx1 = pd.timedelta_range('1 day', '31 day', freq='D', name='idx')

        for idx in [idx1]:
            result = idx[0]
            self.assertEqual(result, pd.Timedelta('1 day'))

            result = idx[0:5]
            expected = pd.timedelta_range('1 day', '5 day', freq='D',
                                          name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)

            result = idx[0:10:2]
            expected = pd.timedelta_range('1 day', '9 day', freq='2D',
                                          name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)

            result = idx[-20:-5:3]
            expected = pd.timedelta_range('12 day', '24 day', freq='3D',
                                          name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)

            result = idx[4::-1]
            expected = TimedeltaIndex(['5 day', '4 day', '3 day',
                                       '2 day', '1 day'],
                                      freq='-1D', name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)

    def test_drop_duplicates_metadata(self):
        # GH 10115
        idx = pd.timedelta_range('1 day', '31 day', freq='D', name='idx')
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
        base = pd.timedelta_range('1 day', '31 day', freq='D', name='idx')
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
        idx1 = pd.timedelta_range('1 day', '31 day', freq='D', name='idx')

        for idx in [idx1]:
            result = idx.take([0])
            self.assertEqual(result, pd.Timedelta('1 day'))

            result = idx.take([-1])
            self.assertEqual(result, pd.Timedelta('31 day'))

            result = idx.take([0, 1, 2])
            expected = pd.timedelta_range('1 day', '3 day', freq='D',
                                          name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)

            result = idx.take([0, 2, 4])
            expected = pd.timedelta_range('1 day', '5 day', freq='2D',
                                          name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)

            result = idx.take([7, 4, 1])
            expected = pd.timedelta_range('8 day', '2 day', freq='-3D',
                                          name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)

            result = idx.take([3, 2, 5])
            expected = TimedeltaIndex(['4 day', '3 day', '6 day'], name='idx')
            self.assert_index_equal(result, expected)
            self.assertIsNone(result.freq)

            result = idx.take([-3, 2, 5])
            expected = TimedeltaIndex(['29 day', '3 day', '6 day'], name='idx')
            self.assert_index_equal(result, expected)
            self.assertIsNone(result.freq)

    def test_take_invalid_kwargs(self):
        idx = pd.timedelta_range('1 day', '31 day', freq='D', name='idx')
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
        for freq in ['D', '3D', '-3D', 'H', '2H', '-2H', 'T', '2T', 'S', '-3S'
                     ]:
            idx = pd.timedelta_range('1', freq=freq, periods=10)
            result = pd.TimedeltaIndex(idx.asi8, freq='infer')
            tm.assert_index_equal(idx, result)
            self.assertEqual(result.freq, freq)

    def test_nat_new(self):

        idx = pd.timedelta_range('1', freq='D', periods=5, name='x')
        result = idx._nat_new()
        exp = pd.TimedeltaIndex([pd.NaT] * 5, name='x')
        tm.assert_index_equal(result, exp)

        result = idx._nat_new(box=False)
        exp = np.array([tslib.iNaT] * 5, dtype=np.int64)
        tm.assert_numpy_array_equal(result, exp)

    def test_shift(self):
        # GH 9903
        idx = pd.TimedeltaIndex([], name='xxx')
        tm.assert_index_equal(idx.shift(0, freq='H'), idx)
        tm.assert_index_equal(idx.shift(3, freq='H'), idx)

        idx = pd.TimedeltaIndex(['5 hours', '6 hours', '9 hours'], name='xxx')
        tm.assert_index_equal(idx.shift(0, freq='H'), idx)
        exp = pd.TimedeltaIndex(['8 hours', '9 hours', '12 hours'], name='xxx')
        tm.assert_index_equal(idx.shift(3, freq='H'), exp)
        exp = pd.TimedeltaIndex(['2 hours', '3 hours', '6 hours'], name='xxx')
        tm.assert_index_equal(idx.shift(-3, freq='H'), exp)

        tm.assert_index_equal(idx.shift(0, freq='T'), idx)
        exp = pd.TimedeltaIndex(['05:03:00', '06:03:00', '9:03:00'],
                                name='xxx')
        tm.assert_index_equal(idx.shift(3, freq='T'), exp)
        exp = pd.TimedeltaIndex(['04:57:00', '05:57:00', '8:57:00'],
                                name='xxx')
        tm.assert_index_equal(idx.shift(-3, freq='T'), exp)

    def test_repeat(self):
        index = pd.timedelta_range('1 days', periods=2, freq='D')
        exp = pd.TimedeltaIndex(['1 days', '1 days', '2 days', '2 days'])
        for res in [index.repeat(2), np.repeat(index, 2)]:
            tm.assert_index_equal(res, exp)
            self.assertIsNone(res.freq)

        index = TimedeltaIndex(['1 days', 'NaT', '3 days'])
        exp = TimedeltaIndex(['1 days', '1 days', '1 days',
                              'NaT', 'NaT', 'NaT',
                              '3 days', '3 days', '3 days'])
        for res in [index.repeat(3), np.repeat(index, 3)]:
            tm.assert_index_equal(res, exp)
            self.assertIsNone(res.freq)

    def test_nat(self):
        self.assertIs(pd.TimedeltaIndex._na_value, pd.NaT)
        self.assertIs(pd.TimedeltaIndex([])._na_value, pd.NaT)

        idx = pd.TimedeltaIndex(['1 days', '2 days'])
        self.assertTrue(idx._can_hold_na)

        tm.assert_numpy_array_equal(idx._isnan, np.array([False, False]))
        self.assertFalse(idx.hasnans)
        tm.assert_numpy_array_equal(idx._nan_idxs,
                                    np.array([], dtype=np.intp))

        idx = pd.TimedeltaIndex(['1 days', 'NaT'])
        self.assertTrue(idx._can_hold_na)

        tm.assert_numpy_array_equal(idx._isnan, np.array([False, True]))
        self.assertTrue(idx.hasnans)
        tm.assert_numpy_array_equal(idx._nan_idxs,
                                    np.array([1], dtype=np.intp))

    def test_equals(self):
        # GH 13107
        idx = pd.TimedeltaIndex(['1 days', '2 days', 'NaT'])
        self.assertTrue(idx.equals(idx))
        self.assertTrue(idx.equals(idx.copy()))
        self.assertTrue(idx.equals(idx.asobject))
        self.assertTrue(idx.asobject.equals(idx))
        self.assertTrue(idx.asobject.equals(idx.asobject))
        self.assertFalse(idx.equals(list(idx)))
        self.assertFalse(idx.equals(pd.Series(idx)))

        idx2 = pd.TimedeltaIndex(['2 days', '1 days', 'NaT'])
        self.assertFalse(idx.equals(idx2))
        self.assertFalse(idx.equals(idx2.copy()))
        self.assertFalse(idx.equals(idx2.asobject))
        self.assertFalse(idx.asobject.equals(idx2))
        self.assertFalse(idx.asobject.equals(idx2.asobject))
        self.assertFalse(idx.equals(list(idx2)))
        self.assertFalse(idx.equals(pd.Series(idx2)))


class TestTimedeltas(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_ops(self):

        td = Timedelta(10, unit='d')
        self.assertEqual(-td, Timedelta(-10, unit='d'))
        self.assertEqual(+td, Timedelta(10, unit='d'))
        self.assertEqual(td - td, Timedelta(0, unit='ns'))
        self.assertTrue((td - pd.NaT) is pd.NaT)
        self.assertEqual(td + td, Timedelta(20, unit='d'))
        self.assertTrue((td + pd.NaT) is pd.NaT)
        self.assertEqual(td * 2, Timedelta(20, unit='d'))
        self.assertTrue((td * pd.NaT) is pd.NaT)
        self.assertEqual(td / 2, Timedelta(5, unit='d'))
        self.assertEqual(abs(td), td)
        self.assertEqual(abs(-td), td)
        self.assertEqual(td / td, 1)
        self.assertTrue((td / pd.NaT) is np.nan)

        # invert
        self.assertEqual(-td, Timedelta('-10d'))
        self.assertEqual(td * -1, Timedelta('-10d'))
        self.assertEqual(-1 * td, Timedelta('-10d'))
        self.assertEqual(abs(-td), Timedelta('10d'))

        # invalid
        self.assertRaises(TypeError, lambda: Timedelta(11, unit='d') // 2)

        # invalid multiply with another timedelta
        self.assertRaises(TypeError, lambda: td * td)

        # can't operate with integers
        self.assertRaises(TypeError, lambda: td + 2)
        self.assertRaises(TypeError, lambda: td - 2)

    def test_ops_offsets(self):
        td = Timedelta(10, unit='d')
        self.assertEqual(Timedelta(241, unit='h'), td + pd.offsets.Hour(1))
        self.assertEqual(Timedelta(241, unit='h'), pd.offsets.Hour(1) + td)
        self.assertEqual(240, td / pd.offsets.Hour(1))
        self.assertEqual(1 / 240.0, pd.offsets.Hour(1) / td)
        self.assertEqual(Timedelta(239, unit='h'), td - pd.offsets.Hour(1))
        self.assertEqual(Timedelta(-239, unit='h'), pd.offsets.Hour(1) - td)

    def test_ops_ndarray(self):
        td = Timedelta('1 day')

        # timedelta, timedelta
        other = pd.to_timedelta(['1 day']).values
        expected = pd.to_timedelta(['2 days']).values
        self.assert_numpy_array_equal(td + other, expected)
        if LooseVersion(np.__version__) >= '1.8':
            self.assert_numpy_array_equal(other + td, expected)
        self.assertRaises(TypeError, lambda: td + np.array([1]))
        self.assertRaises(TypeError, lambda: np.array([1]) + td)

        expected = pd.to_timedelta(['0 days']).values
        self.assert_numpy_array_equal(td - other, expected)
        if LooseVersion(np.__version__) >= '1.8':
            self.assert_numpy_array_equal(-other + td, expected)
        self.assertRaises(TypeError, lambda: td - np.array([1]))
        self.assertRaises(TypeError, lambda: np.array([1]) - td)

        expected = pd.to_timedelta(['2 days']).values
        self.assert_numpy_array_equal(td * np.array([2]), expected)
        self.assert_numpy_array_equal(np.array([2]) * td, expected)
        self.assertRaises(TypeError, lambda: td * other)
        self.assertRaises(TypeError, lambda: other * td)

        self.assert_numpy_array_equal(td / other,
                                      np.array([1], dtype=np.float64))
        if LooseVersion(np.__version__) >= '1.8':
            self.assert_numpy_array_equal(other / td,
                                          np.array([1], dtype=np.float64))

        # timedelta, datetime
        other = pd.to_datetime(['2000-01-01']).values
        expected = pd.to_datetime(['2000-01-02']).values
        self.assert_numpy_array_equal(td + other, expected)
        if LooseVersion(np.__version__) >= '1.8':
            self.assert_numpy_array_equal(other + td, expected)

        expected = pd.to_datetime(['1999-12-31']).values
        self.assert_numpy_array_equal(-td + other, expected)
        if LooseVersion(np.__version__) >= '1.8':
            self.assert_numpy_array_equal(other - td, expected)

    def test_ops_series(self):
        # regression test for GH8813
        td = Timedelta('1 day')
        other = pd.Series([1, 2])
        expected = pd.Series(pd.to_timedelta(['1 day', '2 days']))
        tm.assert_series_equal(expected, td * other)
        tm.assert_series_equal(expected, other * td)

    def test_ops_series_object(self):
        # GH 13043
        s = pd.Series([pd.Timestamp('2015-01-01', tz='US/Eastern'),
                       pd.Timestamp('2015-01-01', tz='Asia/Tokyo')],
                      name='xxx')
        self.assertEqual(s.dtype, object)

        exp = pd.Series([pd.Timestamp('2015-01-02', tz='US/Eastern'),
                         pd.Timestamp('2015-01-02', tz='Asia/Tokyo')],
                        name='xxx')
        tm.assert_series_equal(s + pd.Timedelta('1 days'), exp)
        tm.assert_series_equal(pd.Timedelta('1 days') + s, exp)

        # object series & object series
        s2 = pd.Series([pd.Timestamp('2015-01-03', tz='US/Eastern'),
                        pd.Timestamp('2015-01-05', tz='Asia/Tokyo')],
                       name='xxx')
        self.assertEqual(s2.dtype, object)
        exp = pd.Series([pd.Timedelta('2 days'), pd.Timedelta('4 days')],
                        name='xxx')
        tm.assert_series_equal(s2 - s, exp)
        tm.assert_series_equal(s - s2, -exp)

        s = pd.Series([pd.Timedelta('01:00:00'), pd.Timedelta('02:00:00')],
                      name='xxx', dtype=object)
        self.assertEqual(s.dtype, object)

        exp = pd.Series([pd.Timedelta('01:30:00'), pd.Timedelta('02:30:00')],
                        name='xxx')
        tm.assert_series_equal(s + pd.Timedelta('00:30:00'), exp)
        tm.assert_series_equal(pd.Timedelta('00:30:00') + s, exp)

    def test_ops_notimplemented(self):
        class Other:
            pass

        other = Other()

        td = Timedelta('1 day')
        self.assertTrue(td.__add__(other) is NotImplemented)
        self.assertTrue(td.__sub__(other) is NotImplemented)
        self.assertTrue(td.__truediv__(other) is NotImplemented)
        self.assertTrue(td.__mul__(other) is NotImplemented)
        self.assertTrue(td.__floordiv__(td) is NotImplemented)

    def test_ops_error_str(self):
        # GH 13624
        tdi = TimedeltaIndex(['1 day', '2 days'])

        for l, r in [(tdi, 'a'), ('a', tdi)]:
            with tm.assertRaises(TypeError):
                l + r

            with tm.assertRaises(TypeError):
                l > r

            with tm.assertRaises(TypeError):
                l == r

            with tm.assertRaises(TypeError):
                l != r

    def test_timedelta_ops(self):
        # GH4984
        # make sure ops return Timedelta
        s = Series([Timestamp('20130101') + timedelta(seconds=i * i)
                    for i in range(10)])
        td = s.diff()

        result = td.mean()
        expected = to_timedelta(timedelta(seconds=9))
        self.assertEqual(result, expected)

        result = td.to_frame().mean()
        self.assertEqual(result[0], expected)

        result = td.quantile(.1)
        expected = Timedelta(np.timedelta64(2600, 'ms'))
        self.assertEqual(result, expected)

        result = td.median()
        expected = to_timedelta('00:00:09')
        self.assertEqual(result, expected)

        result = td.to_frame().median()
        self.assertEqual(result[0], expected)

        # GH 6462
        # consistency in returned values for sum
        result = td.sum()
        expected = to_timedelta('00:01:21')
        self.assertEqual(result, expected)

        result = td.to_frame().sum()
        self.assertEqual(result[0], expected)

        # std
        result = td.std()
        expected = to_timedelta(Series(td.dropna().values).std())
        self.assertEqual(result, expected)

        result = td.to_frame().std()
        self.assertEqual(result[0], expected)

        # invalid ops
        for op in ['skew', 'kurt', 'sem', 'prod']:
            self.assertRaises(TypeError, getattr(td, op))

        # GH 10040
        # make sure NaT is properly handled by median()
        s = Series([Timestamp('2015-02-03'), Timestamp('2015-02-07')])
        self.assertEqual(s.diff().median(), timedelta(days=4))

        s = Series([Timestamp('2015-02-03'), Timestamp('2015-02-07'),
                    Timestamp('2015-02-15')])
        self.assertEqual(s.diff().median(), timedelta(days=6))

    def test_timedelta_ops_scalar(self):
        # GH 6808
        base = pd.to_datetime('20130101 09:01:12.123456')
        expected_add = pd.to_datetime('20130101 09:01:22.123456')
        expected_sub = pd.to_datetime('20130101 09:01:02.123456')

        for offset in [pd.to_timedelta(10, unit='s'), timedelta(seconds=10),
                       np.timedelta64(10, 's'),
                       np.timedelta64(10000000000, 'ns'),
                       pd.offsets.Second(10)]:
            result = base + offset
            self.assertEqual(result, expected_add)

            result = base - offset
            self.assertEqual(result, expected_sub)

        base = pd.to_datetime('20130102 09:01:12.123456')
        expected_add = pd.to_datetime('20130103 09:01:22.123456')
        expected_sub = pd.to_datetime('20130101 09:01:02.123456')

        for offset in [pd.to_timedelta('1 day, 00:00:10'),
                       pd.to_timedelta('1 days, 00:00:10'),
                       timedelta(days=1, seconds=10),
                       np.timedelta64(1, 'D') + np.timedelta64(10, 's'),
                       pd.offsets.Day() + pd.offsets.Second(10)]:
            result = base + offset
            self.assertEqual(result, expected_add)

            result = base - offset
            self.assertEqual(result, expected_sub)

    def test_timedelta_ops_with_missing_values(self):
        # setup
        s1 = pd.to_timedelta(Series(['00:00:01']))
        s2 = pd.to_timedelta(Series(['00:00:02']))
        sn = pd.to_timedelta(Series([pd.NaT]))
        df1 = DataFrame(['00:00:01']).apply(pd.to_timedelta)
        df2 = DataFrame(['00:00:02']).apply(pd.to_timedelta)
        dfn = DataFrame([pd.NaT]).apply(pd.to_timedelta)
        scalar1 = pd.to_timedelta('00:00:01')
        scalar2 = pd.to_timedelta('00:00:02')
        timedelta_NaT = pd.to_timedelta('NaT')
        NA = np.nan

        actual = scalar1 + scalar1
        self.assertEqual(actual, scalar2)
        actual = scalar2 - scalar1
        self.assertEqual(actual, scalar1)

        actual = s1 + s1
        assert_series_equal(actual, s2)
        actual = s2 - s1
        assert_series_equal(actual, s1)

        actual = s1 + scalar1
        assert_series_equal(actual, s2)
        actual = scalar1 + s1
        assert_series_equal(actual, s2)
        actual = s2 - scalar1
        assert_series_equal(actual, s1)
        actual = -scalar1 + s2
        assert_series_equal(actual, s1)

        actual = s1 + timedelta_NaT
        assert_series_equal(actual, sn)
        actual = timedelta_NaT + s1
        assert_series_equal(actual, sn)
        actual = s1 - timedelta_NaT
        assert_series_equal(actual, sn)
        actual = -timedelta_NaT + s1
        assert_series_equal(actual, sn)

        actual = s1 + NA
        assert_series_equal(actual, sn)
        actual = NA + s1
        assert_series_equal(actual, sn)
        actual = s1 - NA
        assert_series_equal(actual, sn)
        actual = -NA + s1
        assert_series_equal(actual, sn)

        actual = s1 + pd.NaT
        assert_series_equal(actual, sn)
        actual = s2 - pd.NaT
        assert_series_equal(actual, sn)

        actual = s1 + df1
        assert_frame_equal(actual, df2)
        actual = s2 - df1
        assert_frame_equal(actual, df1)
        actual = df1 + s1
        assert_frame_equal(actual, df2)
        actual = df2 - s1
        assert_frame_equal(actual, df1)

        actual = df1 + df1
        assert_frame_equal(actual, df2)
        actual = df2 - df1
        assert_frame_equal(actual, df1)

        actual = df1 + scalar1
        assert_frame_equal(actual, df2)
        actual = df2 - scalar1
        assert_frame_equal(actual, df1)

        actual = df1 + timedelta_NaT
        assert_frame_equal(actual, dfn)
        actual = df1 - timedelta_NaT
        assert_frame_equal(actual, dfn)

        actual = df1 + NA
        assert_frame_equal(actual, dfn)
        actual = df1 - NA
        assert_frame_equal(actual, dfn)

        actual = df1 + pd.NaT  # NaT is datetime, not timedelta
        assert_frame_equal(actual, dfn)
        actual = df1 - pd.NaT
        assert_frame_equal(actual, dfn)

    def test_compare_timedelta_series(self):
        # regresssion test for GH5963
        s = pd.Series([timedelta(days=1), timedelta(days=2)])
        actual = s > timedelta(days=1)
        expected = pd.Series([False, True])
        tm.assert_series_equal(actual, expected)

    def test_compare_timedelta_ndarray(self):
        # GH11835
        periods = [Timedelta('0 days 01:00:00'), Timedelta('0 days 01:00:00')]
        arr = np.array(periods)
        result = arr[0] > arr
        expected = np.array([False, False])
        self.assert_numpy_array_equal(result, expected)


class TestSlicing(tm.TestCase):

    def test_tdi_ops_attributes(self):
        rng = timedelta_range('2 days', periods=5, freq='2D', name='x')

        result = rng + 1
        exp = timedelta_range('4 days', periods=5, freq='2D', name='x')
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.freq, '2D')

        result = rng - 2
        exp = timedelta_range('-2 days', periods=5, freq='2D', name='x')
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.freq, '2D')

        result = rng * 2
        exp = timedelta_range('4 days', periods=5, freq='4D', name='x')
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.freq, '4D')

        result = rng / 2
        exp = timedelta_range('1 days', periods=5, freq='D', name='x')
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.freq, 'D')

        result = -rng
        exp = timedelta_range('-2 days', periods=5, freq='-2D', name='x')
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.freq, '-2D')

        rng = pd.timedelta_range('-2 days', periods=5, freq='D', name='x')

        result = abs(rng)
        exp = TimedeltaIndex(['2 days', '1 days', '0 days', '1 days',
                              '2 days'], name='x')
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.freq, None)

    def test_add_overflow(self):
        # see gh-14068
        msg = "too (big|large) to convert"
        with tm.assertRaisesRegexp(OverflowError, msg):
            to_timedelta(106580, 'D') + Timestamp('2000')
        with tm.assertRaisesRegexp(OverflowError, msg):
            Timestamp('2000') + to_timedelta(106580, 'D')

        _NaT = int(pd.NaT) + 1
        msg = "Overflow in int64 addition"
        with tm.assertRaisesRegexp(OverflowError, msg):
            to_timedelta([106580], 'D') + Timestamp('2000')
        with tm.assertRaisesRegexp(OverflowError, msg):
            Timestamp('2000') + to_timedelta([106580], 'D')
        with tm.assertRaisesRegexp(OverflowError, msg):
            to_timedelta([_NaT]) - Timedelta('1 days')
        with tm.assertRaisesRegexp(OverflowError, msg):
            to_timedelta(['5 days', _NaT]) - Timedelta('1 days')
        with tm.assertRaisesRegexp(OverflowError, msg):
            (to_timedelta([_NaT, '5 days', '1 hours']) -
             to_timedelta(['7 seconds', _NaT, '4 hours']))

        # These should not overflow!
        exp = TimedeltaIndex([pd.NaT])
        result = to_timedelta([pd.NaT]) - Timedelta('1 days')
        tm.assert_index_equal(result, exp)

        exp = TimedeltaIndex(['4 days', pd.NaT])
        result = to_timedelta(['5 days', pd.NaT]) - Timedelta('1 days')
        tm.assert_index_equal(result, exp)

        exp = TimedeltaIndex([pd.NaT, pd.NaT, '5 hours'])
        result = (to_timedelta([pd.NaT, '5 days', '1 hours']) +
                  to_timedelta(['7 seconds', pd.NaT, '4 hours']))
        tm.assert_index_equal(result, exp)
