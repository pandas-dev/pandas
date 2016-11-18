from __future__ import print_function
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
from pandas import (Series, Index, Int64Index, Timestamp, Period,
                    DatetimeIndex, PeriodIndex, TimedeltaIndex,
                    Timedelta, timedelta_range, date_range, Float64Index,
                    _np_version_under1p10)
import pandas.tslib as tslib
import pandas.tseries.period as period

import pandas.util.testing as tm

from pandas.tests.test_base import Ops


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


class TestPeriodIndexOps(Ops):
    def setUp(self):
        super(TestPeriodIndexOps, self).setUp()
        mask = lambda x: (isinstance(x, DatetimeIndex) or
                          isinstance(x, PeriodIndex))
        self.is_valid_objs = [o for o in self.objs if mask(o)]
        self.not_valid_objs = [o for o in self.objs if not mask(o)]

    def test_ops_properties(self):
        self.check_ops_properties(
            ['year', 'month', 'day', 'hour', 'minute', 'second', 'weekofyear',
             'week', 'dayofweek', 'dayofyear', 'quarter'])
        self.check_ops_properties(['qyear'],
                                  lambda x: isinstance(x, PeriodIndex))

    def test_asobject_tolist(self):
        idx = pd.period_range(start='2013-01-01', periods=4, freq='M',
                              name='idx')
        expected_list = [pd.Period('2013-01-31', freq='M'),
                         pd.Period('2013-02-28', freq='M'),
                         pd.Period('2013-03-31', freq='M'),
                         pd.Period('2013-04-30', freq='M')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.asobject
        self.assertTrue(isinstance(result, Index))
        self.assertEqual(result.dtype, object)
        self.assert_index_equal(result, expected)
        self.assertEqual(result.name, expected.name)
        self.assertEqual(idx.tolist(), expected_list)

        idx = PeriodIndex(['2013-01-01', '2013-01-02', 'NaT',
                           '2013-01-04'], freq='D', name='idx')
        expected_list = [pd.Period('2013-01-01', freq='D'),
                         pd.Period('2013-01-02', freq='D'),
                         pd.Period('NaT', freq='D'),
                         pd.Period('2013-01-04', freq='D')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.asobject
        self.assertTrue(isinstance(result, Index))
        self.assertEqual(result.dtype, object)
        tm.assert_index_equal(result, expected)
        for i in [0, 1, 3]:
            self.assertEqual(result[i], expected[i])
        self.assertIs(result[2], pd.NaT)
        self.assertEqual(result.name, expected.name)

        result_list = idx.tolist()
        for i in [0, 1, 3]:
            self.assertEqual(result_list[i], expected_list[i])
        self.assertIs(result_list[2], pd.NaT)

    def test_minmax(self):

        # monotonic
        idx1 = pd.PeriodIndex([pd.NaT, '2011-01-01', '2011-01-02',
                               '2011-01-03'], freq='D')
        self.assertTrue(idx1.is_monotonic)

        # non-monotonic
        idx2 = pd.PeriodIndex(['2011-01-01', pd.NaT, '2011-01-03',
                               '2011-01-02', pd.NaT], freq='D')
        self.assertFalse(idx2.is_monotonic)

        for idx in [idx1, idx2]:
            self.assertEqual(idx.min(), pd.Period('2011-01-01', freq='D'))
            self.assertEqual(idx.max(), pd.Period('2011-01-03', freq='D'))
        self.assertEqual(idx1.argmin(), 1)
        self.assertEqual(idx2.argmin(), 0)
        self.assertEqual(idx1.argmax(), 3)
        self.assertEqual(idx2.argmax(), 2)

        for op in ['min', 'max']:
            # Return NaT
            obj = PeriodIndex([], freq='M')
            result = getattr(obj, op)()
            self.assertIs(result, tslib.NaT)

            obj = PeriodIndex([pd.NaT], freq='M')
            result = getattr(obj, op)()
            self.assertIs(result, tslib.NaT)

            obj = PeriodIndex([pd.NaT, pd.NaT, pd.NaT], freq='M')
            result = getattr(obj, op)()
            self.assertIs(result, tslib.NaT)

    def test_numpy_minmax(self):
        pr = pd.period_range(start='2016-01-15', end='2016-01-20')

        self.assertEqual(np.min(pr), Period('2016-01-15', freq='D'))
        self.assertEqual(np.max(pr), Period('2016-01-20', freq='D'))

        errmsg = "the 'out' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, errmsg, np.min, pr, out=0)
        tm.assertRaisesRegexp(ValueError, errmsg, np.max, pr, out=0)

        self.assertEqual(np.argmin(pr), 0)
        self.assertEqual(np.argmax(pr), 5)

        if not _np_version_under1p10:
            errmsg = "the 'out' parameter is not supported"
            tm.assertRaisesRegexp(ValueError, errmsg, np.argmin, pr, out=0)
            tm.assertRaisesRegexp(ValueError, errmsg, np.argmax, pr, out=0)

    def test_representation(self):
        # GH 7601
        idx1 = PeriodIndex([], freq='D')
        idx2 = PeriodIndex(['2011-01-01'], freq='D')
        idx3 = PeriodIndex(['2011-01-01', '2011-01-02'], freq='D')
        idx4 = PeriodIndex(['2011-01-01', '2011-01-02', '2011-01-03'],
                           freq='D')
        idx5 = PeriodIndex(['2011', '2012', '2013'], freq='A')
        idx6 = PeriodIndex(['2011-01-01 09:00', '2012-02-01 10:00',
                            'NaT'], freq='H')
        idx7 = pd.period_range('2013Q1', periods=1, freq="Q")
        idx8 = pd.period_range('2013Q1', periods=2, freq="Q")
        idx9 = pd.period_range('2013Q1', periods=3, freq="Q")
        idx10 = PeriodIndex(['2011-01-01', '2011-02-01'], freq='3D')

        exp1 = """PeriodIndex([], dtype='period[D]', freq='D')"""

        exp2 = """PeriodIndex(['2011-01-01'], dtype='period[D]', freq='D')"""

        exp3 = ("PeriodIndex(['2011-01-01', '2011-01-02'], dtype='period[D]', "
                "freq='D')")

        exp4 = ("PeriodIndex(['2011-01-01', '2011-01-02', '2011-01-03'], "
                "dtype='period[D]', freq='D')")

        exp5 = ("PeriodIndex(['2011', '2012', '2013'], dtype='period[A-DEC]', "
                "freq='A-DEC')")

        exp6 = ("PeriodIndex(['2011-01-01 09:00', '2012-02-01 10:00', 'NaT'], "
                "dtype='period[H]', freq='H')")

        exp7 = ("PeriodIndex(['2013Q1'], dtype='period[Q-DEC]', "
                "freq='Q-DEC')")

        exp8 = ("PeriodIndex(['2013Q1', '2013Q2'], dtype='period[Q-DEC]', "
                "freq='Q-DEC')")

        exp9 = ("PeriodIndex(['2013Q1', '2013Q2', '2013Q3'], "
                "dtype='period[Q-DEC]', freq='Q-DEC')")

        exp10 = ("PeriodIndex(['2011-01-01', '2011-02-01'], "
                 "dtype='period[3D]', freq='3D')")

        for idx, expected in zip([idx1, idx2, idx3, idx4, idx5,
                                  idx6, idx7, idx8, idx9, idx10],
                                 [exp1, exp2, exp3, exp4, exp5,
                                  exp6, exp7, exp8, exp9, exp10]):
            for func in ['__repr__', '__unicode__', '__str__']:
                result = getattr(idx, func)()
                self.assertEqual(result, expected)

    def test_representation_to_series(self):
        # GH 10971
        idx1 = PeriodIndex([], freq='D')
        idx2 = PeriodIndex(['2011-01-01'], freq='D')
        idx3 = PeriodIndex(['2011-01-01', '2011-01-02'], freq='D')
        idx4 = PeriodIndex(['2011-01-01', '2011-01-02',
                            '2011-01-03'], freq='D')
        idx5 = PeriodIndex(['2011', '2012', '2013'], freq='A')
        idx6 = PeriodIndex(['2011-01-01 09:00', '2012-02-01 10:00',
                            'NaT'], freq='H')

        idx7 = pd.period_range('2013Q1', periods=1, freq="Q")
        idx8 = pd.period_range('2013Q1', periods=2, freq="Q")
        idx9 = pd.period_range('2013Q1', periods=3, freq="Q")

        exp1 = """Series([], dtype: object)"""

        exp2 = """0   2011-01-01
dtype: object"""

        exp3 = """0   2011-01-01
1   2011-01-02
dtype: object"""

        exp4 = """0   2011-01-01
1   2011-01-02
2   2011-01-03
dtype: object"""

        exp5 = """0   2011
1   2012
2   2013
dtype: object"""

        exp6 = """0   2011-01-01 09:00
1   2012-02-01 10:00
2                NaT
dtype: object"""

        exp7 = """0   2013Q1
dtype: object"""

        exp8 = """0   2013Q1
1   2013Q2
dtype: object"""

        exp9 = """0   2013Q1
1   2013Q2
2   2013Q3
dtype: object"""

        for idx, expected in zip([idx1, idx2, idx3, idx4, idx5,
                                  idx6, idx7, idx8, idx9],
                                 [exp1, exp2, exp3, exp4, exp5,
                                  exp6, exp7, exp8, exp9]):
            result = repr(pd.Series(idx))
            self.assertEqual(result, expected)

    def test_summary(self):
        # GH9116
        idx1 = PeriodIndex([], freq='D')
        idx2 = PeriodIndex(['2011-01-01'], freq='D')
        idx3 = PeriodIndex(['2011-01-01', '2011-01-02'], freq='D')
        idx4 = PeriodIndex(
            ['2011-01-01', '2011-01-02', '2011-01-03'], freq='D')
        idx5 = PeriodIndex(['2011', '2012', '2013'], freq='A')
        idx6 = PeriodIndex(
            ['2011-01-01 09:00', '2012-02-01 10:00', 'NaT'], freq='H')

        idx7 = pd.period_range('2013Q1', periods=1, freq="Q")
        idx8 = pd.period_range('2013Q1', periods=2, freq="Q")
        idx9 = pd.period_range('2013Q1', periods=3, freq="Q")

        exp1 = """PeriodIndex: 0 entries
Freq: D"""

        exp2 = """PeriodIndex: 1 entries, 2011-01-01 to 2011-01-01
Freq: D"""

        exp3 = """PeriodIndex: 2 entries, 2011-01-01 to 2011-01-02
Freq: D"""

        exp4 = """PeriodIndex: 3 entries, 2011-01-01 to 2011-01-03
Freq: D"""

        exp5 = """PeriodIndex: 3 entries, 2011 to 2013
Freq: A-DEC"""

        exp6 = """PeriodIndex: 3 entries, 2011-01-01 09:00 to NaT
Freq: H"""

        exp7 = """PeriodIndex: 1 entries, 2013Q1 to 2013Q1
Freq: Q-DEC"""

        exp8 = """PeriodIndex: 2 entries, 2013Q1 to 2013Q2
Freq: Q-DEC"""

        exp9 = """PeriodIndex: 3 entries, 2013Q1 to 2013Q3
Freq: Q-DEC"""

        for idx, expected in zip([idx1, idx2, idx3, idx4, idx5,
                                  idx6, idx7, idx8, idx9],
                                 [exp1, exp2, exp3, exp4, exp5,
                                  exp6, exp7, exp8, exp9]):
            result = idx.summary()
            self.assertEqual(result, expected)

    def test_resolution(self):
        for freq, expected in zip(['A', 'Q', 'M', 'D', 'H',
                                   'T', 'S', 'L', 'U'],
                                  ['day', 'day', 'day', 'day',
                                   'hour', 'minute', 'second',
                                   'millisecond', 'microsecond']):

            idx = pd.period_range(start='2013-04-01', periods=30, freq=freq)
            self.assertEqual(idx.resolution, expected)

    def test_union(self):
        # union
        rng1 = pd.period_range('1/1/2000', freq='D', periods=5)
        other1 = pd.period_range('1/6/2000', freq='D', periods=5)
        expected1 = pd.period_range('1/1/2000', freq='D', periods=10)

        rng2 = pd.period_range('1/1/2000', freq='D', periods=5)
        other2 = pd.period_range('1/4/2000', freq='D', periods=5)
        expected2 = pd.period_range('1/1/2000', freq='D', periods=8)

        rng3 = pd.period_range('1/1/2000', freq='D', periods=5)
        other3 = pd.PeriodIndex([], freq='D')
        expected3 = pd.period_range('1/1/2000', freq='D', periods=5)

        rng4 = pd.period_range('2000-01-01 09:00', freq='H', periods=5)
        other4 = pd.period_range('2000-01-02 09:00', freq='H', periods=5)
        expected4 = pd.PeriodIndex(['2000-01-01 09:00', '2000-01-01 10:00',
                                    '2000-01-01 11:00', '2000-01-01 12:00',
                                    '2000-01-01 13:00', '2000-01-02 09:00',
                                    '2000-01-02 10:00', '2000-01-02 11:00',
                                    '2000-01-02 12:00', '2000-01-02 13:00'],
                                   freq='H')

        rng5 = pd.PeriodIndex(['2000-01-01 09:01', '2000-01-01 09:03',
                               '2000-01-01 09:05'], freq='T')
        other5 = pd.PeriodIndex(['2000-01-01 09:01', '2000-01-01 09:05'
                                                     '2000-01-01 09:08'],
                                freq='T')
        expected5 = pd.PeriodIndex(['2000-01-01 09:01', '2000-01-01 09:03',
                                    '2000-01-01 09:05', '2000-01-01 09:08'],
                                   freq='T')

        rng6 = pd.period_range('2000-01-01', freq='M', periods=7)
        other6 = pd.period_range('2000-04-01', freq='M', periods=7)
        expected6 = pd.period_range('2000-01-01', freq='M', periods=10)

        rng7 = pd.period_range('2003-01-01', freq='A', periods=5)
        other7 = pd.period_range('1998-01-01', freq='A', periods=8)
        expected7 = pd.period_range('1998-01-01', freq='A', periods=10)

        for rng, other, expected in [(rng1, other1, expected1),
                                     (rng2, other2, expected2),
                                     (rng3, other3, expected3), (rng4, other4,
                                                                 expected4),
                                     (rng5, other5, expected5), (rng6, other6,
                                                                 expected6),
                                     (rng7, other7, expected7)]:

            result_union = rng.union(other)
            tm.assert_index_equal(result_union, expected)

    def test_add_iadd(self):
        rng = pd.period_range('1/1/2000', freq='D', periods=5)
        other = pd.period_range('1/6/2000', freq='D', periods=5)

        # previously performed setop union, now raises TypeError (GH14164)
        with tm.assertRaises(TypeError):
            rng + other

        with tm.assertRaises(TypeError):
            rng += other

        # offset
        # DateOffset
        rng = pd.period_range('2014', '2024', freq='A')
        result = rng + pd.offsets.YearEnd(5)
        expected = pd.period_range('2019', '2029', freq='A')
        tm.assert_index_equal(result, expected)
        rng += pd.offsets.YearEnd(5)
        tm.assert_index_equal(rng, expected)

        for o in [pd.offsets.YearBegin(2), pd.offsets.MonthBegin(1),
                  pd.offsets.Minute(), np.timedelta64(365, 'D'),
                  timedelta(365), Timedelta(days=365)]:
            msg = ('Input has different freq(=.+)? '
                   'from PeriodIndex\\(freq=A-DEC\\)')
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                rng + o

        rng = pd.period_range('2014-01', '2016-12', freq='M')
        result = rng + pd.offsets.MonthEnd(5)
        expected = pd.period_range('2014-06', '2017-05', freq='M')
        tm.assert_index_equal(result, expected)
        rng += pd.offsets.MonthEnd(5)
        tm.assert_index_equal(rng, expected)

        for o in [pd.offsets.YearBegin(2), pd.offsets.MonthBegin(1),
                  pd.offsets.Minute(), np.timedelta64(365, 'D'),
                  timedelta(365), Timedelta(days=365)]:
            rng = pd.period_range('2014-01', '2016-12', freq='M')
            msg = 'Input has different freq(=.+)? from PeriodIndex\\(freq=M\\)'
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                rng + o

        # Tick
        offsets = [pd.offsets.Day(3), timedelta(days=3),
                   np.timedelta64(3, 'D'), pd.offsets.Hour(72),
                   timedelta(minutes=60 * 24 * 3), np.timedelta64(72, 'h'),
                   Timedelta('72:00:00')]
        for delta in offsets:
            rng = pd.period_range('2014-05-01', '2014-05-15', freq='D')
            result = rng + delta
            expected = pd.period_range('2014-05-04', '2014-05-18', freq='D')
            tm.assert_index_equal(result, expected)
            rng += delta
            tm.assert_index_equal(rng, expected)

        for o in [pd.offsets.YearBegin(2), pd.offsets.MonthBegin(1),
                  pd.offsets.Minute(), np.timedelta64(4, 'h'),
                  timedelta(hours=23), Timedelta('23:00:00')]:
            rng = pd.period_range('2014-05-01', '2014-05-15', freq='D')
            msg = 'Input has different freq(=.+)? from PeriodIndex\\(freq=D\\)'
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                rng + o

        offsets = [pd.offsets.Hour(2), timedelta(hours=2),
                   np.timedelta64(2, 'h'), pd.offsets.Minute(120),
                   timedelta(minutes=120), np.timedelta64(120, 'm'),
                   Timedelta(minutes=120)]
        for delta in offsets:
            rng = pd.period_range('2014-01-01 10:00', '2014-01-05 10:00',
                                  freq='H')
            result = rng + delta
            expected = pd.period_range('2014-01-01 12:00', '2014-01-05 12:00',
                                       freq='H')
            tm.assert_index_equal(result, expected)
            rng += delta
            tm.assert_index_equal(rng, expected)

        for delta in [pd.offsets.YearBegin(2), timedelta(minutes=30),
                      np.timedelta64(30, 's'), Timedelta(seconds=30)]:
            rng = pd.period_range('2014-01-01 10:00', '2014-01-05 10:00',
                                  freq='H')
            msg = 'Input has different freq(=.+)? from PeriodIndex\\(freq=H\\)'
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                result = rng + delta
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                rng += delta

        # int
        rng = pd.period_range('2000-01-01 09:00', freq='H', periods=10)
        result = rng + 1
        expected = pd.period_range('2000-01-01 10:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)
        rng += 1
        tm.assert_index_equal(rng, expected)

    def test_difference(self):
        # diff
        rng1 = pd.period_range('1/1/2000', freq='D', periods=5)
        other1 = pd.period_range('1/6/2000', freq='D', periods=5)
        expected1 = pd.period_range('1/1/2000', freq='D', periods=5)

        rng2 = pd.period_range('1/1/2000', freq='D', periods=5)
        other2 = pd.period_range('1/4/2000', freq='D', periods=5)
        expected2 = pd.period_range('1/1/2000', freq='D', periods=3)

        rng3 = pd.period_range('1/1/2000', freq='D', periods=5)
        other3 = pd.PeriodIndex([], freq='D')
        expected3 = pd.period_range('1/1/2000', freq='D', periods=5)

        rng4 = pd.period_range('2000-01-01 09:00', freq='H', periods=5)
        other4 = pd.period_range('2000-01-02 09:00', freq='H', periods=5)
        expected4 = rng4

        rng5 = pd.PeriodIndex(['2000-01-01 09:01', '2000-01-01 09:03',
                               '2000-01-01 09:05'], freq='T')
        other5 = pd.PeriodIndex(
            ['2000-01-01 09:01', '2000-01-01 09:05'], freq='T')
        expected5 = pd.PeriodIndex(['2000-01-01 09:03'], freq='T')

        rng6 = pd.period_range('2000-01-01', freq='M', periods=7)
        other6 = pd.period_range('2000-04-01', freq='M', periods=7)
        expected6 = pd.period_range('2000-01-01', freq='M', periods=3)

        rng7 = pd.period_range('2003-01-01', freq='A', periods=5)
        other7 = pd.period_range('1998-01-01', freq='A', periods=8)
        expected7 = pd.period_range('2006-01-01', freq='A', periods=2)

        for rng, other, expected in [(rng1, other1, expected1),
                                     (rng2, other2, expected2),
                                     (rng3, other3, expected3),
                                     (rng4, other4, expected4),
                                     (rng5, other5, expected5),
                                     (rng6, other6, expected6),
                                     (rng7, other7, expected7), ]:
            result_union = rng.difference(other)
            tm.assert_index_equal(result_union, expected)

    def test_sub_isub(self):

        # previously performed setop, now raises TypeError (GH14164)
        # TODO needs to wait on #13077 for decision on result type
        rng = pd.period_range('1/1/2000', freq='D', periods=5)
        other = pd.period_range('1/6/2000', freq='D', periods=5)

        with tm.assertRaises(TypeError):
            rng - other

        with tm.assertRaises(TypeError):
            rng -= other

        # offset
        # DateOffset
        rng = pd.period_range('2014', '2024', freq='A')
        result = rng - pd.offsets.YearEnd(5)
        expected = pd.period_range('2009', '2019', freq='A')
        tm.assert_index_equal(result, expected)
        rng -= pd.offsets.YearEnd(5)
        tm.assert_index_equal(rng, expected)

        for o in [pd.offsets.YearBegin(2), pd.offsets.MonthBegin(1),
                  pd.offsets.Minute(), np.timedelta64(365, 'D'),
                  timedelta(365)]:
            rng = pd.period_range('2014', '2024', freq='A')
            msg = ('Input has different freq(=.+)? '
                   'from PeriodIndex\\(freq=A-DEC\\)')
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                rng - o

        rng = pd.period_range('2014-01', '2016-12', freq='M')
        result = rng - pd.offsets.MonthEnd(5)
        expected = pd.period_range('2013-08', '2016-07', freq='M')
        tm.assert_index_equal(result, expected)
        rng -= pd.offsets.MonthEnd(5)
        tm.assert_index_equal(rng, expected)

        for o in [pd.offsets.YearBegin(2), pd.offsets.MonthBegin(1),
                  pd.offsets.Minute(), np.timedelta64(365, 'D'),
                  timedelta(365)]:
            rng = pd.period_range('2014-01', '2016-12', freq='M')
            msg = 'Input has different freq(=.+)? from PeriodIndex\\(freq=M\\)'
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                rng - o

        # Tick
        offsets = [pd.offsets.Day(3), timedelta(days=3),
                   np.timedelta64(3, 'D'), pd.offsets.Hour(72),
                   timedelta(minutes=60 * 24 * 3), np.timedelta64(72, 'h')]
        for delta in offsets:
            rng = pd.period_range('2014-05-01', '2014-05-15', freq='D')
            result = rng - delta
            expected = pd.period_range('2014-04-28', '2014-05-12', freq='D')
            tm.assert_index_equal(result, expected)
            rng -= delta
            tm.assert_index_equal(rng, expected)

        for o in [pd.offsets.YearBegin(2), pd.offsets.MonthBegin(1),
                  pd.offsets.Minute(), np.timedelta64(4, 'h'),
                  timedelta(hours=23)]:
            rng = pd.period_range('2014-05-01', '2014-05-15', freq='D')
            msg = 'Input has different freq(=.+)? from PeriodIndex\\(freq=D\\)'
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                rng - o

        offsets = [pd.offsets.Hour(2), timedelta(hours=2),
                   np.timedelta64(2, 'h'), pd.offsets.Minute(120),
                   timedelta(minutes=120), np.timedelta64(120, 'm')]
        for delta in offsets:
            rng = pd.period_range('2014-01-01 10:00', '2014-01-05 10:00',
                                  freq='H')
            result = rng - delta
            expected = pd.period_range('2014-01-01 08:00', '2014-01-05 08:00',
                                       freq='H')
            tm.assert_index_equal(result, expected)
            rng -= delta
            tm.assert_index_equal(rng, expected)

        for delta in [pd.offsets.YearBegin(2), timedelta(minutes=30),
                      np.timedelta64(30, 's')]:
            rng = pd.period_range('2014-01-01 10:00', '2014-01-05 10:00',
                                  freq='H')
            msg = 'Input has different freq(=.+)? from PeriodIndex\\(freq=H\\)'
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                result = rng + delta
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                rng += delta

        # int
        rng = pd.period_range('2000-01-01 09:00', freq='H', periods=10)
        result = rng - 1
        expected = pd.period_range('2000-01-01 08:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)
        rng -= 1
        tm.assert_index_equal(rng, expected)

    def test_comp_nat(self):
        left = pd.PeriodIndex([pd.Period('2011-01-01'), pd.NaT,
                               pd.Period('2011-01-03')])
        right = pd.PeriodIndex([pd.NaT, pd.NaT, pd.Period('2011-01-03')])

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
        idx = pd.period_range('2011-01-01 09:00', freq='H', periods=10)
        # create repeated values, 'n'th element is repeated by n+1 times
        idx = PeriodIndex(np.repeat(idx.values, range(1, len(idx) + 1)),
                          freq='H')

        exp_idx = PeriodIndex(['2011-01-01 18:00', '2011-01-01 17:00',
                               '2011-01-01 16:00', '2011-01-01 15:00',
                               '2011-01-01 14:00', '2011-01-01 13:00',
                               '2011-01-01 12:00', '2011-01-01 11:00',
                               '2011-01-01 10:00',
                               '2011-01-01 09:00'], freq='H')
        expected = Series(range(10, 0, -1), index=exp_idx, dtype='int64')

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(), expected)

        expected = pd.period_range('2011-01-01 09:00', freq='H',
                                   periods=10)
        tm.assert_index_equal(idx.unique(), expected)

        idx = PeriodIndex(['2013-01-01 09:00', '2013-01-01 09:00',
                           '2013-01-01 09:00', '2013-01-01 08:00',
                           '2013-01-01 08:00', pd.NaT], freq='H')

        exp_idx = PeriodIndex(['2013-01-01 09:00', '2013-01-01 08:00'],
                              freq='H')
        expected = Series([3, 2], index=exp_idx)

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(), expected)

        exp_idx = PeriodIndex(['2013-01-01 09:00', '2013-01-01 08:00',
                               pd.NaT], freq='H')
        expected = Series([3, 2, 1], index=exp_idx)

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(dropna=False), expected)

        tm.assert_index_equal(idx.unique(), exp_idx)

    def test_drop_duplicates_metadata(self):
        # GH 10115
        idx = pd.period_range('2011-01-01', '2011-01-31', freq='D', name='idx')
        result = idx.drop_duplicates()
        self.assert_index_equal(idx, result)
        self.assertEqual(idx.freq, result.freq)

        idx_dup = idx.append(idx)  # freq will not be reset
        result = idx_dup.drop_duplicates()
        self.assert_index_equal(idx, result)
        self.assertEqual(idx.freq, result.freq)

    def test_drop_duplicates(self):
        # to check Index/Series compat
        base = pd.period_range('2011-01-01', '2011-01-31', freq='D',
                               name='idx')
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

    def test_order_compat(self):
        def _check_freq(index, expected_index):
            if isinstance(index, PeriodIndex):
                self.assertEqual(index.freq, expected_index.freq)

        pidx = PeriodIndex(['2011', '2012', '2013'], name='pidx', freq='A')
        # for compatibility check
        iidx = Index([2011, 2012, 2013], name='idx')
        for idx in [pidx, iidx]:
            ordered = idx.sort_values()
            self.assert_index_equal(ordered, idx)
            _check_freq(ordered, idx)

            ordered = idx.sort_values(ascending=False)
            self.assert_index_equal(ordered, idx[::-1])
            _check_freq(ordered, idx[::-1])

            ordered, indexer = idx.sort_values(return_indexer=True)
            self.assert_index_equal(ordered, idx)
            self.assert_numpy_array_equal(indexer,
                                          np.array([0, 1, 2]),
                                          check_dtype=False)
            _check_freq(ordered, idx)

            ordered, indexer = idx.sort_values(return_indexer=True,
                                               ascending=False)
            self.assert_index_equal(ordered, idx[::-1])
            self.assert_numpy_array_equal(indexer,
                                          np.array([2, 1, 0]),
                                          check_dtype=False)
            _check_freq(ordered, idx[::-1])

        pidx = PeriodIndex(['2011', '2013', '2015', '2012',
                            '2011'], name='pidx', freq='A')
        pexpected = PeriodIndex(
            ['2011', '2011', '2012', '2013', '2015'], name='pidx', freq='A')
        # for compatibility check
        iidx = Index([2011, 2013, 2015, 2012, 2011], name='idx')
        iexpected = Index([2011, 2011, 2012, 2013, 2015], name='idx')
        for idx, expected in [(pidx, pexpected), (iidx, iexpected)]:
            ordered = idx.sort_values()
            self.assert_index_equal(ordered, expected)
            _check_freq(ordered, idx)

            ordered = idx.sort_values(ascending=False)
            self.assert_index_equal(ordered, expected[::-1])
            _check_freq(ordered, idx)

            ordered, indexer = idx.sort_values(return_indexer=True)
            self.assert_index_equal(ordered, expected)

            exp = np.array([0, 4, 3, 1, 2])
            self.assert_numpy_array_equal(indexer, exp, check_dtype=False)
            _check_freq(ordered, idx)

            ordered, indexer = idx.sort_values(return_indexer=True,
                                               ascending=False)
            self.assert_index_equal(ordered, expected[::-1])

            exp = np.array([2, 1, 3, 4, 0])
            self.assert_numpy_array_equal(indexer, exp,
                                          check_dtype=False)
            _check_freq(ordered, idx)

        pidx = PeriodIndex(['2011', '2013', 'NaT', '2011'], name='pidx',
                           freq='D')

        result = pidx.sort_values()
        expected = PeriodIndex(['NaT', '2011', '2011', '2013'],
                               name='pidx', freq='D')
        self.assert_index_equal(result, expected)
        self.assertEqual(result.freq, 'D')

        result = pidx.sort_values(ascending=False)
        expected = PeriodIndex(
            ['2013', '2011', '2011', 'NaT'], name='pidx', freq='D')
        self.assert_index_equal(result, expected)
        self.assertEqual(result.freq, 'D')

    def test_order(self):
        for freq in ['D', '2D', '4D']:
            idx = PeriodIndex(['2011-01-01', '2011-01-02', '2011-01-03'],
                              freq=freq, name='idx')

            ordered = idx.sort_values()
            self.assert_index_equal(ordered, idx)
            self.assertEqual(ordered.freq, idx.freq)

            ordered = idx.sort_values(ascending=False)
            expected = idx[::-1]
            self.assert_index_equal(ordered, expected)
            self.assertEqual(ordered.freq, expected.freq)
            self.assertEqual(ordered.freq, freq)

            ordered, indexer = idx.sort_values(return_indexer=True)
            self.assert_index_equal(ordered, idx)
            self.assert_numpy_array_equal(indexer,
                                          np.array([0, 1, 2]),
                                          check_dtype=False)
            self.assertEqual(ordered.freq, idx.freq)
            self.assertEqual(ordered.freq, freq)

            ordered, indexer = idx.sort_values(return_indexer=True,
                                               ascending=False)
            expected = idx[::-1]
            self.assert_index_equal(ordered, expected)
            self.assert_numpy_array_equal(indexer,
                                          np.array([2, 1, 0]),
                                          check_dtype=False)
            self.assertEqual(ordered.freq, expected.freq)
            self.assertEqual(ordered.freq, freq)

        idx1 = PeriodIndex(['2011-01-01', '2011-01-03', '2011-01-05',
                            '2011-01-02', '2011-01-01'], freq='D', name='idx1')
        exp1 = PeriodIndex(['2011-01-01', '2011-01-01', '2011-01-02',
                            '2011-01-03', '2011-01-05'], freq='D', name='idx1')

        idx2 = PeriodIndex(['2011-01-01', '2011-01-03', '2011-01-05',
                            '2011-01-02', '2011-01-01'],
                           freq='D', name='idx2')
        exp2 = PeriodIndex(['2011-01-01', '2011-01-01', '2011-01-02',
                            '2011-01-03', '2011-01-05'],
                           freq='D', name='idx2')

        idx3 = PeriodIndex([pd.NaT, '2011-01-03', '2011-01-05',
                            '2011-01-02', pd.NaT], freq='D', name='idx3')
        exp3 = PeriodIndex([pd.NaT, pd.NaT, '2011-01-02', '2011-01-03',
                            '2011-01-05'], freq='D', name='idx3')

        for idx, expected in [(idx1, exp1), (idx2, exp2), (idx3, exp3)]:
            ordered = idx.sort_values()
            self.assert_index_equal(ordered, expected)
            self.assertEqual(ordered.freq, 'D')

            ordered = idx.sort_values(ascending=False)
            self.assert_index_equal(ordered, expected[::-1])
            self.assertEqual(ordered.freq, 'D')

            ordered, indexer = idx.sort_values(return_indexer=True)
            self.assert_index_equal(ordered, expected)

            exp = np.array([0, 4, 3, 1, 2])
            self.assert_numpy_array_equal(indexer, exp, check_dtype=False)
            self.assertEqual(ordered.freq, 'D')

            ordered, indexer = idx.sort_values(return_indexer=True,
                                               ascending=False)
            self.assert_index_equal(ordered, expected[::-1])

            exp = np.array([2, 1, 3, 4, 0])
            self.assert_numpy_array_equal(indexer, exp, check_dtype=False)
            self.assertEqual(ordered.freq, 'D')

    def test_getitem(self):
        idx1 = pd.period_range('2011-01-01', '2011-01-31', freq='D',
                               name='idx')

        for idx in [idx1]:
            result = idx[0]
            self.assertEqual(result, pd.Period('2011-01-01', freq='D'))

            result = idx[-1]
            self.assertEqual(result, pd.Period('2011-01-31', freq='D'))

            result = idx[0:5]
            expected = pd.period_range('2011-01-01', '2011-01-05', freq='D',
                                       name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)
            self.assertEqual(result.freq, 'D')

            result = idx[0:10:2]
            expected = pd.PeriodIndex(['2011-01-01', '2011-01-03',
                                       '2011-01-05',
                                       '2011-01-07', '2011-01-09'],
                                      freq='D', name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)
            self.assertEqual(result.freq, 'D')

            result = idx[-20:-5:3]
            expected = pd.PeriodIndex(['2011-01-12', '2011-01-15',
                                       '2011-01-18',
                                       '2011-01-21', '2011-01-24'],
                                      freq='D', name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)
            self.assertEqual(result.freq, 'D')

            result = idx[4::-1]
            expected = PeriodIndex(['2011-01-05', '2011-01-04', '2011-01-03',
                                    '2011-01-02', '2011-01-01'],
                                   freq='D', name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)
            self.assertEqual(result.freq, 'D')

    def test_take(self):
        # GH 10295
        idx1 = pd.period_range('2011-01-01', '2011-01-31', freq='D',
                               name='idx')

        for idx in [idx1]:
            result = idx.take([0])
            self.assertEqual(result, pd.Period('2011-01-01', freq='D'))

            result = idx.take([5])
            self.assertEqual(result, pd.Period('2011-01-06', freq='D'))

            result = idx.take([0, 1, 2])
            expected = pd.period_range('2011-01-01', '2011-01-03', freq='D',
                                       name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, 'D')
            self.assertEqual(result.freq, expected.freq)

            result = idx.take([0, 2, 4])
            expected = pd.PeriodIndex(['2011-01-01', '2011-01-03',
                                       '2011-01-05'], freq='D', name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)
            self.assertEqual(result.freq, 'D')

            result = idx.take([7, 4, 1])
            expected = pd.PeriodIndex(['2011-01-08', '2011-01-05',
                                       '2011-01-02'],
                                      freq='D', name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)
            self.assertEqual(result.freq, 'D')

            result = idx.take([3, 2, 5])
            expected = PeriodIndex(['2011-01-04', '2011-01-03', '2011-01-06'],
                                   freq='D', name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)
            self.assertEqual(result.freq, 'D')

            result = idx.take([-3, 2, 5])
            expected = PeriodIndex(['2011-01-29', '2011-01-03', '2011-01-06'],
                                   freq='D', name='idx')
            self.assert_index_equal(result, expected)
            self.assertEqual(result.freq, expected.freq)
            self.assertEqual(result.freq, 'D')

    def test_nat_new(self):

        idx = pd.period_range('2011-01', freq='M', periods=5, name='x')
        result = idx._nat_new()
        exp = pd.PeriodIndex([pd.NaT] * 5, freq='M', name='x')
        tm.assert_index_equal(result, exp)

        result = idx._nat_new(box=False)
        exp = np.array([tslib.iNaT] * 5, dtype=np.int64)
        tm.assert_numpy_array_equal(result, exp)

    def test_shift(self):
        # GH 9903
        idx = pd.PeriodIndex([], name='xxx', freq='H')

        with tm.assertRaises(TypeError):
            # period shift doesn't accept freq
            idx.shift(1, freq='H')

        tm.assert_index_equal(idx.shift(0), idx)
        tm.assert_index_equal(idx.shift(3), idx)

        idx = pd.PeriodIndex(['2011-01-01 10:00', '2011-01-01 11:00'
                              '2011-01-01 12:00'], name='xxx', freq='H')
        tm.assert_index_equal(idx.shift(0), idx)
        exp = pd.PeriodIndex(['2011-01-01 13:00', '2011-01-01 14:00'
                              '2011-01-01 15:00'], name='xxx', freq='H')
        tm.assert_index_equal(idx.shift(3), exp)
        exp = pd.PeriodIndex(['2011-01-01 07:00', '2011-01-01 08:00'
                              '2011-01-01 09:00'], name='xxx', freq='H')
        tm.assert_index_equal(idx.shift(-3), exp)

    def test_repeat(self):
        index = pd.period_range('2001-01-01', periods=2, freq='D')
        exp = pd.PeriodIndex(['2001-01-01', '2001-01-01',
                              '2001-01-02', '2001-01-02'], freq='D')
        for res in [index.repeat(2), np.repeat(index, 2)]:
            tm.assert_index_equal(res, exp)

        index = pd.period_range('2001-01-01', periods=2, freq='2D')
        exp = pd.PeriodIndex(['2001-01-01', '2001-01-01',
                              '2001-01-03', '2001-01-03'], freq='2D')
        for res in [index.repeat(2), np.repeat(index, 2)]:
            tm.assert_index_equal(res, exp)

        index = pd.PeriodIndex(['2001-01', 'NaT', '2003-01'], freq='M')
        exp = pd.PeriodIndex(['2001-01', '2001-01', '2001-01',
                              'NaT', 'NaT', 'NaT',
                              '2003-01', '2003-01', '2003-01'], freq='M')
        for res in [index.repeat(3), np.repeat(index, 3)]:
            tm.assert_index_equal(res, exp)

    def test_nat(self):
        self.assertIs(pd.PeriodIndex._na_value, pd.NaT)
        self.assertIs(pd.PeriodIndex([], freq='M')._na_value, pd.NaT)

        idx = pd.PeriodIndex(['2011-01-01', '2011-01-02'], freq='D')
        self.assertTrue(idx._can_hold_na)

        tm.assert_numpy_array_equal(idx._isnan, np.array([False, False]))
        self.assertFalse(idx.hasnans)
        tm.assert_numpy_array_equal(idx._nan_idxs,
                                    np.array([], dtype=np.intp))

        idx = pd.PeriodIndex(['2011-01-01', 'NaT'], freq='D')
        self.assertTrue(idx._can_hold_na)

        tm.assert_numpy_array_equal(idx._isnan, np.array([False, True]))
        self.assertTrue(idx.hasnans)
        tm.assert_numpy_array_equal(idx._nan_idxs,
                                    np.array([1], dtype=np.intp))

    def test_equals(self):
        # GH 13107
        for freq in ['D', 'M']:
            idx = pd.PeriodIndex(['2011-01-01', '2011-01-02', 'NaT'],
                                 freq=freq)
            self.assertTrue(idx.equals(idx))
            self.assertTrue(idx.equals(idx.copy()))
            self.assertTrue(idx.equals(idx.asobject))
            self.assertTrue(idx.asobject.equals(idx))
            self.assertTrue(idx.asobject.equals(idx.asobject))
            self.assertFalse(idx.equals(list(idx)))
            self.assertFalse(idx.equals(pd.Series(idx)))

            idx2 = pd.PeriodIndex(['2011-01-01', '2011-01-02', 'NaT'],
                                  freq='H')
            self.assertFalse(idx.equals(idx2))
            self.assertFalse(idx.equals(idx2.copy()))
            self.assertFalse(idx.equals(idx2.asobject))
            self.assertFalse(idx.asobject.equals(idx2))
            self.assertFalse(idx.equals(list(idx2)))
            self.assertFalse(idx.equals(pd.Series(idx2)))

            # same internal, different tz
            idx3 = pd.PeriodIndex._simple_new(idx.asi8, freq='H')
            tm.assert_numpy_array_equal(idx.asi8, idx3.asi8)
            self.assertFalse(idx.equals(idx3))
            self.assertFalse(idx.equals(idx3.copy()))
            self.assertFalse(idx.equals(idx3.asobject))
            self.assertFalse(idx.asobject.equals(idx3))
            self.assertFalse(idx.equals(list(idx3)))
            self.assertFalse(idx.equals(pd.Series(idx3)))


if __name__ == '__main__':
    import nose

    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
