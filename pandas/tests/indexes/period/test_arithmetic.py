# -*- coding: utf-8 -*-
from datetime import timedelta
import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas import (Timedelta,
                    period_range, Period, PeriodIndex,
                    _np_version_under1p10)
import pandas.core.indexes.period as period


class TestPeriodIndexArithmetic(object):
    def test_pi_add_offset_array(self):
        # GH#18849
        pi = pd.PeriodIndex([pd.Period('2015Q1'), pd.Period('2016Q2')])
        offs = np.array([pd.offsets.QuarterEnd(n=1, startingMonth=12),
                         pd.offsets.QuarterEnd(n=-2, startingMonth=12)])
        res = pi + offs
        expected = pd.PeriodIndex([pd.Period('2015Q2'), pd.Period('2015Q4')])
        tm.assert_index_equal(res, expected)

        unanchored = np.array([pd.offsets.Hour(n=1),
                               pd.offsets.Minute(n=-2)])
        with pytest.raises(period.IncompatibleFrequency):
            pi + unanchored
        with pytest.raises(TypeError):
            unanchored + pi

    @pytest.mark.xfail(reason='GH#18824 radd doesnt implement this case')
    def test_pi_radd_offset_array(self):
        # GH#18849
        pi = pd.PeriodIndex([pd.Period('2015Q1'), pd.Period('2016Q2')])
        offs = np.array([pd.offsets.QuarterEnd(n=1, startingMonth=12),
                         pd.offsets.QuarterEnd(n=-2, startingMonth=12)])
        res = offs + pi
        expected = pd.PeriodIndex([pd.Period('2015Q2'), pd.Period('2015Q4')])
        tm.assert_index_equal(res, expected)

    def test_add_iadd(self):
        rng = pd.period_range('1/1/2000', freq='D', periods=5)
        other = pd.period_range('1/6/2000', freq='D', periods=5)

        # previously performed setop union, now raises TypeError (GH14164)
        with pytest.raises(TypeError):
            rng + other

        with pytest.raises(TypeError):
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
            with tm.assert_raises_regex(
                    period.IncompatibleFrequency, msg):
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
            with tm.assert_raises_regex(
                    period.IncompatibleFrequency, msg):
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
            with tm.assert_raises_regex(
                    period.IncompatibleFrequency, msg):
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
            with tm.assert_raises_regex(
                    period.IncompatibleFrequency, msg):
                rng + delta
            with tm.assert_raises_regex(
                    period.IncompatibleFrequency, msg):
                rng += delta

    def test_pi_add_int(self, one):
        # Variants of `one` for #19012
        rng = pd.period_range('2000-01-01 09:00', freq='H', periods=10)
        result = rng + one
        expected = pd.period_range('2000-01-01 10:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)
        rng += one
        tm.assert_index_equal(rng, expected)

    @pytest.mark.parametrize('five', [5, np.array(5, dtype=np.int64)])
    def test_sub(self, five):
        rng = period_range('2007-01', periods=50)

        result = rng - five
        exp = rng + (-five)
        tm.assert_index_equal(result, exp)

    def test_sub_isub(self):

        # previously performed setop, now raises TypeError (GH14164)
        # TODO needs to wait on #13077 for decision on result type
        rng = pd.period_range('1/1/2000', freq='D', periods=5)
        other = pd.period_range('1/6/2000', freq='D', periods=5)

        with pytest.raises(TypeError):
            rng - other

        with pytest.raises(TypeError):
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
            with tm.assert_raises_regex(
                    period.IncompatibleFrequency, msg):
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
            with tm.assert_raises_regex(
                    period.IncompatibleFrequency, msg):
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
            with tm.assert_raises_regex(
                    period.IncompatibleFrequency, msg):
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
            with tm.assert_raises_regex(
                    period.IncompatibleFrequency, msg):
                rng + delta
            with tm.assert_raises_regex(
                    period.IncompatibleFrequency, msg):
                rng += delta

        # int
        rng = pd.period_range('2000-01-01 09:00', freq='H', periods=10)
        result = rng - 1
        expected = pd.period_range('2000-01-01 08:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)
        rng -= 1
        tm.assert_index_equal(rng, expected)


class TestPeriodIndexSeriesMethods(object):
    """ Test PeriodIndex and Period Series Ops consistency """

    def _check(self, values, func, expected):
        idx = pd.PeriodIndex(values)
        result = func(idx)
        if isinstance(expected, pd.Index):
            tm.assert_index_equal(result, expected)
        else:
            # comp op results in bool
            tm.assert_numpy_array_equal(result, expected)

        s = pd.Series(values)
        result = func(s)

        exp = pd.Series(expected, name=values.name)
        tm.assert_series_equal(result, exp)

    def test_pi_ops(self):
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03',
                           '2011-04'], freq='M', name='idx')

        expected = PeriodIndex(['2011-03', '2011-04',
                                '2011-05', '2011-06'], freq='M', name='idx')
        self._check(idx, lambda x: x + 2, expected)
        self._check(idx, lambda x: 2 + x, expected)

        self._check(idx + 2, lambda x: x - 2, idx)
        result = idx - Period('2011-01', freq='M')
        exp = pd.Index([0, 1, 2, 3], name='idx')
        tm.assert_index_equal(result, exp)

        result = Period('2011-01', freq='M') - idx
        exp = pd.Index([0, -1, -2, -3], name='idx')
        tm.assert_index_equal(result, exp)

    def test_pi_ops_errors(self):
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03',
                           '2011-04'], freq='M', name='idx')
        s = pd.Series(idx)

        msg = r"unsupported operand type\(s\)"

        for obj in [idx, s]:
            for ng in ["str", 1.5]:
                with tm.assert_raises_regex(TypeError, msg):
                    obj + ng

                with pytest.raises(TypeError):
                    # error message differs between PY2 and 3
                    ng + obj

                with tm.assert_raises_regex(TypeError, msg):
                    obj - ng

                with pytest.raises(TypeError):
                    np.add(obj, ng)

                if _np_version_under1p10:
                    assert np.add(ng, obj) is NotImplemented
                else:
                    with pytest.raises(TypeError):
                        np.add(ng, obj)

                with pytest.raises(TypeError):
                    np.subtract(obj, ng)

                if _np_version_under1p10:
                    assert np.subtract(ng, obj) is NotImplemented
                else:
                    with pytest.raises(TypeError):
                        np.subtract(ng, obj)

    def test_pi_ops_nat(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT',
                           '2011-04'], freq='M', name='idx')
        expected = PeriodIndex(['2011-03', '2011-04',
                                'NaT', '2011-06'], freq='M', name='idx')
        self._check(idx, lambda x: x + 2, expected)
        self._check(idx, lambda x: 2 + x, expected)
        self._check(idx, lambda x: np.add(x, 2), expected)

        self._check(idx + 2, lambda x: x - 2, idx)
        self._check(idx + 2, lambda x: np.subtract(x, 2), idx)

        # freq with mult
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT',
                           '2011-04'], freq='2M', name='idx')
        expected = PeriodIndex(['2011-07', '2011-08',
                                'NaT', '2011-10'], freq='2M', name='idx')
        self._check(idx, lambda x: x + 3, expected)
        self._check(idx, lambda x: 3 + x, expected)
        self._check(idx, lambda x: np.add(x, 3), expected)

        self._check(idx + 3, lambda x: x - 3, idx)
        self._check(idx + 3, lambda x: np.subtract(x, 3), idx)

    def test_pi_ops_array_int(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT',
                           '2011-04'], freq='M', name='idx')
        f = lambda x: x + np.array([1, 2, 3, 4])
        exp = PeriodIndex(['2011-02', '2011-04', 'NaT',
                           '2011-08'], freq='M', name='idx')
        self._check(idx, f, exp)

        f = lambda x: np.add(x, np.array([4, -1, 1, 2]))
        exp = PeriodIndex(['2011-05', '2011-01', 'NaT',
                           '2011-06'], freq='M', name='idx')
        self._check(idx, f, exp)

        f = lambda x: x - np.array([1, 2, 3, 4])
        exp = PeriodIndex(['2010-12', '2010-12', 'NaT',
                           '2010-12'], freq='M', name='idx')
        self._check(idx, f, exp)

        f = lambda x: np.subtract(x, np.array([3, 2, 3, -2]))
        exp = PeriodIndex(['2010-10', '2010-12', 'NaT',
                           '2011-06'], freq='M', name='idx')
        self._check(idx, f, exp)

    def test_pi_ops_offset(self):
        idx = PeriodIndex(['2011-01-01', '2011-02-01', '2011-03-01',
                           '2011-04-01'], freq='D', name='idx')
        f = lambda x: x + pd.offsets.Day()
        exp = PeriodIndex(['2011-01-02', '2011-02-02', '2011-03-02',
                           '2011-04-02'], freq='D', name='idx')
        self._check(idx, f, exp)

        f = lambda x: x + pd.offsets.Day(2)
        exp = PeriodIndex(['2011-01-03', '2011-02-03', '2011-03-03',
                           '2011-04-03'], freq='D', name='idx')
        self._check(idx, f, exp)

        f = lambda x: x - pd.offsets.Day(2)
        exp = PeriodIndex(['2010-12-30', '2011-01-30', '2011-02-27',
                           '2011-03-30'], freq='D', name='idx')
        self._check(idx, f, exp)

    def test_pi_offset_errors(self):
        idx = PeriodIndex(['2011-01-01', '2011-02-01', '2011-03-01',
                           '2011-04-01'], freq='D', name='idx')
        s = pd.Series(idx)

        # Series op is applied per Period instance, thus error is raised
        # from Period
        msg_idx = r"Input has different freq from PeriodIndex\(freq=D\)"
        msg_s = r"Input cannot be converted to Period\(freq=D\)"
        for obj, msg in [(idx, msg_idx), (s, msg_s)]:
            with tm.assert_raises_regex(
                    period.IncompatibleFrequency, msg):
                obj + pd.offsets.Hour(2)

            with tm.assert_raises_regex(
                    period.IncompatibleFrequency, msg):
                pd.offsets.Hour(2) + obj

            with tm.assert_raises_regex(
                    period.IncompatibleFrequency, msg):
                obj - pd.offsets.Hour(2)

    def test_pi_sub_period(self):
        # GH 13071
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03',
                           '2011-04'], freq='M', name='idx')

        result = idx - pd.Period('2012-01', freq='M')
        exp = pd.Index([-12, -11, -10, -9], name='idx')
        tm.assert_index_equal(result, exp)

        result = np.subtract(idx, pd.Period('2012-01', freq='M'))
        tm.assert_index_equal(result, exp)

        result = pd.Period('2012-01', freq='M') - idx
        exp = pd.Index([12, 11, 10, 9], name='idx')
        tm.assert_index_equal(result, exp)

        result = np.subtract(pd.Period('2012-01', freq='M'), idx)
        if _np_version_under1p10:
            assert result is NotImplemented
        else:
            tm.assert_index_equal(result, exp)

        exp = pd.TimedeltaIndex([np.nan, np.nan, np.nan, np.nan], name='idx')
        tm.assert_index_equal(idx - pd.Period('NaT', freq='M'), exp)
        tm.assert_index_equal(pd.Period('NaT', freq='M') - idx, exp)

    def test_pi_sub_pdnat(self):
        # GH 13071
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT',
                           '2011-04'], freq='M', name='idx')
        exp = pd.TimedeltaIndex([pd.NaT] * 4, name='idx')
        tm.assert_index_equal(pd.NaT - idx, exp)
        tm.assert_index_equal(idx - pd.NaT, exp)

    def test_pi_sub_period_nat(self):
        # GH 13071
        idx = PeriodIndex(['2011-01', 'NaT', '2011-03',
                           '2011-04'], freq='M', name='idx')

        result = idx - pd.Period('2012-01', freq='M')
        exp = pd.Index([-12, np.nan, -10, -9], name='idx')
        tm.assert_index_equal(result, exp)

        result = pd.Period('2012-01', freq='M') - idx
        exp = pd.Index([12, np.nan, 10, 9], name='idx')
        tm.assert_index_equal(result, exp)

        exp = pd.TimedeltaIndex([np.nan, np.nan, np.nan, np.nan], name='idx')
        tm.assert_index_equal(idx - pd.Period('NaT', freq='M'), exp)
        tm.assert_index_equal(pd.Period('NaT', freq='M') - idx, exp)
