from __future__ import print_function
from datetime import timedelta
import numpy as np
import pandas as pd
from pandas import (Series, Index, Period, DatetimeIndex, PeriodIndex,
                    Timedelta, _np_version_under1p10)
import pandas.tslib as tslib
import pandas.tseries.period as period

import pandas.util.testing as tm

from pandas.tests.test_base import Ops


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
