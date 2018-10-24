# -*- coding: utf-8 -*-
# Arithmetc tests for DataFrame/Series/Index/Array classes that should
# behave identically.
# Specifically for Period dtype
import operator

import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm

from pandas.errors import PerformanceWarning

import pandas.core.indexes.period as period
from pandas.core import ops
from pandas import Period, PeriodIndex, period_range, Series
from pandas.tseries.frequencies import to_offset


# ------------------------------------------------------------------
# Comparisons

class TestPeriodIndexComparisons(object):
    def test_pi_cmp_period(self):
        idx = period_range('2007-01', periods=20, freq='M')

        result = idx < idx[10]
        exp = idx.values < idx.values[10]
        tm.assert_numpy_array_equal(result, exp)

    @pytest.mark.parametrize('freq', ['M', '2M', '3M'])
    def test_parr_cmp_period_scalar(self, freq, box):
        # GH#13200
        xbox = np.ndarray if box is pd.Index else box

        base = PeriodIndex(['2011-01', '2011-02', '2011-03', '2011-04'],
                           freq=freq)
        base = tm.box_expected(base, box)
        per = Period('2011-02', freq=freq)

        exp = np.array([False, True, False, False])
        exp = tm.box_expected(exp, xbox)
        tm.assert_equal(base == per, exp)
        tm.assert_equal(per == base, exp)

        exp = np.array([True, False, True, True])
        exp = tm.box_expected(exp, xbox)
        tm.assert_equal(base != per, exp)
        tm.assert_equal(per != base, exp)

        exp = np.array([False, False, True, True])
        exp = tm.box_expected(exp, xbox)
        tm.assert_equal(base > per, exp)
        tm.assert_equal(per < base, exp)

        exp = np.array([True, False, False, False])
        exp = tm.box_expected(exp, xbox)
        tm.assert_equal(base < per, exp)
        tm.assert_equal(per > base, exp)

        exp = np.array([False, True, True, True])
        exp = tm.box_expected(exp, xbox)
        tm.assert_equal(base >= per, exp)
        tm.assert_equal(per <= base, exp)

        exp = np.array([True, True, False, False])
        exp = tm.box_expected(exp, xbox)
        tm.assert_equal(base <= per, exp)
        tm.assert_equal(per >= base, exp)

    @pytest.mark.parametrize('freq', ['M', '2M', '3M'])
    def test_parr_cmp_pi(self, freq, box_df_fail):
        # GH#13200
        box = box_df_fail
        xbox = np.ndarray if box is pd.Index else box

        base = PeriodIndex(['2011-01', '2011-02', '2011-03', '2011-04'],
                           freq=freq)
        base = tm.box_expected(base, box)

        # TODO: could also box idx?
        idx = PeriodIndex(['2011-02', '2011-01', '2011-03', '2011-05'],
                          freq=freq)

        exp = np.array([False, False, True, False])
        exp = tm.box_expected(exp, xbox)
        tm.assert_equal(base == idx, exp)

        exp = np.array([True, True, False, True])
        exp = tm.box_expected(exp, xbox)
        tm.assert_equal(base != idx, exp)

        exp = np.array([False, True, False, False])
        exp = tm.box_expected(exp, xbox)
        tm.assert_equal(base > idx, exp)

        exp = np.array([True, False, False, True])
        exp = tm.box_expected(exp, xbox)
        tm.assert_equal(base < idx, exp)

        exp = np.array([False, True, True, False])
        exp = tm.box_expected(exp, xbox)
        tm.assert_equal(base >= idx, exp)

        exp = np.array([True, False, True, True])
        exp = tm.box_expected(exp, xbox)
        tm.assert_equal(base <= idx, exp)

    @pytest.mark.parametrize('freq', ['M', '2M', '3M'])
    def test_parr_cmp_pi_mismatched_freq_raises(self, freq, box_df_fail):
        # GH#13200
        # different base freq
        box = box_df_fail

        base = PeriodIndex(['2011-01', '2011-02', '2011-03', '2011-04'],
                           freq=freq)
        base = tm.box_expected(base, box)

        msg = "Input has different freq=A-DEC from "
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            base <= Period('2011', freq='A')

        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            Period('2011', freq='A') >= base

        # TODO: Could parametrize over boxes for idx?
        idx = PeriodIndex(['2011', '2012', '2013', '2014'], freq='A')
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            base <= idx

        # Different frequency
        msg = "Input has different freq=4M from "
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            base <= Period('2011', freq='4M')

        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            Period('2011', freq='4M') >= base

        idx = PeriodIndex(['2011', '2012', '2013', '2014'], freq='4M')
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            base <= idx

    @pytest.mark.parametrize('freq', ['M', '2M', '3M'])
    def test_pi_cmp_nat(self, freq):
        idx1 = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-05'], freq=freq)

        result = idx1 > Period('2011-02', freq=freq)
        exp = np.array([False, False, False, True])
        tm.assert_numpy_array_equal(result, exp)
        result = Period('2011-02', freq=freq) < idx1
        tm.assert_numpy_array_equal(result, exp)

        result = idx1 == Period('NaT', freq=freq)
        exp = np.array([False, False, False, False])
        tm.assert_numpy_array_equal(result, exp)
        result = Period('NaT', freq=freq) == idx1
        tm.assert_numpy_array_equal(result, exp)

        result = idx1 != Period('NaT', freq=freq)
        exp = np.array([True, True, True, True])
        tm.assert_numpy_array_equal(result, exp)
        result = Period('NaT', freq=freq) != idx1
        tm.assert_numpy_array_equal(result, exp)

        idx2 = PeriodIndex(['2011-02', '2011-01', '2011-04', 'NaT'], freq=freq)
        result = idx1 < idx2
        exp = np.array([True, False, False, False])
        tm.assert_numpy_array_equal(result, exp)

        result = idx1 == idx2
        exp = np.array([False, False, False, False])
        tm.assert_numpy_array_equal(result, exp)

        result = idx1 != idx2
        exp = np.array([True, True, True, True])
        tm.assert_numpy_array_equal(result, exp)

        result = idx1 == idx1
        exp = np.array([True, True, False, True])
        tm.assert_numpy_array_equal(result, exp)

        result = idx1 != idx1
        exp = np.array([False, False, True, False])
        tm.assert_numpy_array_equal(result, exp)

    @pytest.mark.parametrize('freq', ['M', '2M', '3M'])
    def test_pi_cmp_nat_mismatched_freq_raises(self, freq):
        idx1 = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-05'], freq=freq)

        diff = PeriodIndex(['2011-02', '2011-01', '2011-04', 'NaT'], freq='4M')
        msg = "Input has different freq=4M from PeriodIndex"
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            idx1 > diff

        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            idx1 == diff

    # TODO: De-duplicate with test_pi_cmp_nat
    @pytest.mark.parametrize('dtype', [object, None])
    def test_comp_nat(self, dtype):
        left = pd.PeriodIndex([pd.Period('2011-01-01'), pd.NaT,
                               pd.Period('2011-01-03')])
        right = pd.PeriodIndex([pd.NaT, pd.NaT, pd.Period('2011-01-03')])

        if dtype is not None:
            left = left.astype(dtype)
            right = right.astype(dtype)

        result = left == right
        expected = np.array([False, False, True])
        tm.assert_numpy_array_equal(result, expected)

        result = left != right
        expected = np.array([True, True, False])
        tm.assert_numpy_array_equal(result, expected)

        expected = np.array([False, False, False])
        tm.assert_numpy_array_equal(left == pd.NaT, expected)
        tm.assert_numpy_array_equal(pd.NaT == right, expected)

        expected = np.array([True, True, True])
        tm.assert_numpy_array_equal(left != pd.NaT, expected)
        tm.assert_numpy_array_equal(pd.NaT != left, expected)

        expected = np.array([False, False, False])
        tm.assert_numpy_array_equal(left < pd.NaT, expected)
        tm.assert_numpy_array_equal(pd.NaT > left, expected)


class TestPeriodSeriesComparisons(object):
    def test_cmp_series_period_series_mixed_freq(self):
        # GH#13200
        base = Series([Period('2011', freq='A'),
                       Period('2011-02', freq='M'),
                       Period('2013', freq='A'),
                       Period('2011-04', freq='M')])

        ser = Series([Period('2012', freq='A'),
                      Period('2011-01', freq='M'),
                      Period('2013', freq='A'),
                      Period('2011-05', freq='M')])

        exp = Series([False, False, True, False])
        tm.assert_series_equal(base == ser, exp)

        exp = Series([True, True, False, True])
        tm.assert_series_equal(base != ser, exp)

        exp = Series([False, True, False, False])
        tm.assert_series_equal(base > ser, exp)

        exp = Series([True, False, False, True])
        tm.assert_series_equal(base < ser, exp)

        exp = Series([False, True, True, False])
        tm.assert_series_equal(base >= ser, exp)

        exp = Series([True, False, True, True])
        tm.assert_series_equal(base <= ser, exp)


# ------------------------------------------------------------------
# Arithmetic

class TestPeriodFrameArithmetic(object):

    def test_ops_frame_period(self):
        # GH#13043
        df = pd.DataFrame({'A': [pd.Period('2015-01', freq='M'),
                                 pd.Period('2015-02', freq='M')],
                           'B': [pd.Period('2014-01', freq='M'),
                                 pd.Period('2014-02', freq='M')]})
        assert df['A'].dtype == object
        assert df['B'].dtype == object

        p = pd.Period('2015-03', freq='M')
        off = p.freq
        # dtype will be object because of original dtype
        exp = pd.DataFrame({'A': np.array([2 * off, 1 * off], dtype=object),
                            'B': np.array([14 * off, 13 * off], dtype=object)})
        tm.assert_frame_equal(p - df, exp)
        tm.assert_frame_equal(df - p, -1 * exp)

        df2 = pd.DataFrame({'A': [pd.Period('2015-05', freq='M'),
                                  pd.Period('2015-06', freq='M')],
                            'B': [pd.Period('2015-05', freq='M'),
                                  pd.Period('2015-06', freq='M')]})
        assert df2['A'].dtype == object
        assert df2['B'].dtype == object

        exp = pd.DataFrame({'A': np.array([4 * off, 4 * off], dtype=object),
                            'B': np.array([16 * off, 16 * off], dtype=object)})
        tm.assert_frame_equal(df2 - df, exp)
        tm.assert_frame_equal(df - df2, -1 * exp)


class TestPeriodIndexArithmetic(object):
    # ---------------------------------------------------------------
    # __add__/__sub__ with PeriodIndex
    # PeriodIndex + other is defined for integers and timedelta-like others
    # PeriodIndex - other is defined for integers, timedelta-like others,
    #   and PeriodIndex (with matching freq)

    def test_parr_add_iadd_parr_raises(self, box_df_broadcast_failure):
        box = box_df_broadcast_failure

        rng = pd.period_range('1/1/2000', freq='D', periods=5)
        other = pd.period_range('1/6/2000', freq='D', periods=5)
        # TODO: parametrize over boxes for other?

        rng = tm.box_expected(rng, box)
        # An earlier implementation of PeriodIndex addition performed
        # a set operation (union).  This has since been changed to
        # raise a TypeError. See GH#14164 and GH#13077 for historical
        # reference.
        with pytest.raises(TypeError):
            rng + other

        with pytest.raises(TypeError):
            rng += other

    def test_pi_sub_isub_pi(self):
        # GH#20049
        # For historical reference see GH#14164, GH#13077.
        # PeriodIndex subtraction originally performed set difference,
        # then changed to raise TypeError before being implemented in GH#20049
        rng = pd.period_range('1/1/2000', freq='D', periods=5)
        other = pd.period_range('1/6/2000', freq='D', periods=5)

        off = rng.freq
        expected = pd.Index([-5 * off] * 5)
        result = rng - other
        tm.assert_index_equal(result, expected)

        rng -= other
        tm.assert_index_equal(rng, expected)

    def test_pi_sub_pi_with_nat(self):
        rng = pd.period_range('1/1/2000', freq='D', periods=5)
        other = rng[1:].insert(0, pd.NaT)
        assert other[1:].equals(rng[1:])

        result = rng - other
        off = rng.freq
        expected = pd.Index([pd.NaT, 0 * off, 0 * off, 0 * off, 0 * off])
        tm.assert_index_equal(result, expected)

    def test_parr_sub_pi_mismatched_freq(self, box_df_broadcast_failure):
        box = box_df_broadcast_failure

        rng = pd.period_range('1/1/2000', freq='D', periods=5)
        other = pd.period_range('1/6/2000', freq='H', periods=5)
        # TODO: parametrize over boxes for other?

        rng = tm.box_expected(rng, box)
        with pytest.raises(period.IncompatibleFrequency):
            rng - other

    # -------------------------------------------------------------
    # Invalid Operations

    @pytest.mark.parametrize('other', [3.14, np.array([2.0, 3.0])])
    @pytest.mark.parametrize('op', [operator.add, ops.radd,
                                    operator.sub, ops.rsub])
    def test_pi_add_sub_float(self, op, other, box):
        if box is pd.DataFrame and isinstance(other, np.ndarray):
            pytest.xfail(reason="Tries to broadcast incorrectly")

        dti = pd.DatetimeIndex(['2011-01-01', '2011-01-02'], freq='D')
        pi = dti.to_period('D')
        pi = tm.box_expected(pi, box)
        with pytest.raises(TypeError):
            op(pi, other)

    @pytest.mark.parametrize('other', [pd.Timestamp.now(),
                                       pd.Timestamp.now().to_pydatetime(),
                                       pd.Timestamp.now().to_datetime64()])
    def test_pi_add_sub_datetime(self, other):
        # GH#23215
        rng = pd.period_range('1/1/2000', freq='D', periods=3)

        with pytest.raises(TypeError):
            rng + other
        with pytest.raises(TypeError):
            other + rng
        with pytest.raises(TypeError):
            rng - other
        with pytest.raises(TypeError):
            other - rng

    # -----------------------------------------------------------------
    # __add__/__sub__ with ndarray[datetime64] and ndarray[timedelta64]

    def test_pi_add_sub_dt64_array_raises(self):
        rng = pd.period_range('1/1/2000', freq='D', periods=3)
        dti = pd.date_range('2016-01-01', periods=3)
        dtarr = dti.values

        with pytest.raises(TypeError):
            rng + dtarr
        with pytest.raises(TypeError):
            dtarr + rng

        with pytest.raises(TypeError):
            rng - dtarr
        with pytest.raises(TypeError):
            dtarr - rng

    def test_pi_add_sub_td64_array_non_tick_raises(self):
        rng = pd.period_range('1/1/2000', freq='Q', periods=3)
        tdi = pd.TimedeltaIndex(['-1 Day', '-1 Day', '-1 Day'])
        tdarr = tdi.values

        with pytest.raises(period.IncompatibleFrequency):
            rng + tdarr
        with pytest.raises(period.IncompatibleFrequency):
            tdarr + rng

        with pytest.raises(period.IncompatibleFrequency):
            rng - tdarr
        with pytest.raises(period.IncompatibleFrequency):
            tdarr - rng

    def test_pi_add_sub_td64_array_tick(self):
        # PeriodIndex + Timedelta-like is allowed only with
        #   tick-like frequencies
        rng = pd.period_range('1/1/2000', freq='90D', periods=3)
        tdi = pd.TimedeltaIndex(['-1 Day', '-1 Day', '-1 Day'])
        tdarr = tdi.values

        expected = pd.period_range('12/31/1999', freq='90D', periods=3)
        result = rng + tdi
        tm.assert_index_equal(result, expected)
        result = rng + tdarr
        tm.assert_index_equal(result, expected)
        result = tdi + rng
        tm.assert_index_equal(result, expected)
        result = tdarr + rng
        tm.assert_index_equal(result, expected)

        expected = pd.period_range('1/2/2000', freq='90D', periods=3)

        result = rng - tdi
        tm.assert_index_equal(result, expected)
        result = rng - tdarr
        tm.assert_index_equal(result, expected)

        with pytest.raises(TypeError):
            tdarr - rng

        with pytest.raises(TypeError):
            tdi - rng

    # -----------------------------------------------------------------
    # operations with array/Index of DateOffset objects

    @pytest.mark.parametrize('box', [np.array, pd.Index])
    def test_pi_add_offset_array(self, box):
        # GH#18849
        pi = pd.PeriodIndex([pd.Period('2015Q1'), pd.Period('2016Q2')])
        offs = box([pd.offsets.QuarterEnd(n=1, startingMonth=12),
                    pd.offsets.QuarterEnd(n=-2, startingMonth=12)])
        expected = pd.PeriodIndex([pd.Period('2015Q2'), pd.Period('2015Q4')])

        with tm.assert_produces_warning(PerformanceWarning):
            res = pi + offs
        tm.assert_index_equal(res, expected)

        with tm.assert_produces_warning(PerformanceWarning):
            res2 = offs + pi
        tm.assert_index_equal(res2, expected)

        unanchored = np.array([pd.offsets.Hour(n=1),
                               pd.offsets.Minute(n=-2)])
        # addition/subtraction ops with incompatible offsets should issue
        # a PerformanceWarning and _then_ raise a TypeError.
        with pytest.raises(period.IncompatibleFrequency):
            with tm.assert_produces_warning(PerformanceWarning):
                pi + unanchored
        with pytest.raises(period.IncompatibleFrequency):
            with tm.assert_produces_warning(PerformanceWarning):
                unanchored + pi

    @pytest.mark.parametrize('box', [np.array, pd.Index])
    def test_pi_sub_offset_array(self, box):
        # GH#18824
        pi = pd.PeriodIndex([pd.Period('2015Q1'), pd.Period('2016Q2')])
        other = box([pd.offsets.QuarterEnd(n=1, startingMonth=12),
                     pd.offsets.QuarterEnd(n=-2, startingMonth=12)])

        expected = PeriodIndex([pi[n] - other[n] for n in range(len(pi))])

        with tm.assert_produces_warning(PerformanceWarning):
            res = pi - other
        tm.assert_index_equal(res, expected)

        anchored = box([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)])

        # addition/subtraction ops with anchored offsets should issue
        # a PerformanceWarning and _then_ raise a TypeError.
        with pytest.raises(period.IncompatibleFrequency):
            with tm.assert_produces_warning(PerformanceWarning):
                pi - anchored
        with pytest.raises(period.IncompatibleFrequency):
            with tm.assert_produces_warning(PerformanceWarning):
                anchored - pi

    def test_pi_add_iadd_int(self, one):
        # Variants of `one` for #19012
        rng = pd.period_range('2000-01-01 09:00', freq='H', periods=10)
        result = rng + one
        expected = pd.period_range('2000-01-01 10:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)
        rng += one
        tm.assert_index_equal(rng, expected)

    def test_pi_sub_isub_int(self, one):
        """
        PeriodIndex.__sub__ and __isub__ with several representations of
        the integer 1, e.g. int, long, np.int64, np.uint8, ...
        """
        rng = pd.period_range('2000-01-01 09:00', freq='H', periods=10)
        result = rng - one
        expected = pd.period_range('2000-01-01 08:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)
        rng -= one
        tm.assert_index_equal(rng, expected)

    @pytest.mark.parametrize('five', [5, np.array(5, dtype=np.int64)])
    def test_pi_sub_intlike(self, five):
        rng = period_range('2007-01', periods=50)

        result = rng - five
        exp = rng + (-five)
        tm.assert_index_equal(result, exp)

    def test_pi_sub_isub_offset(self):
        # offset
        # DateOffset
        rng = pd.period_range('2014', '2024', freq='A')
        result = rng - pd.offsets.YearEnd(5)
        expected = pd.period_range('2009', '2019', freq='A')
        tm.assert_index_equal(result, expected)
        rng -= pd.offsets.YearEnd(5)
        tm.assert_index_equal(rng, expected)

        rng = pd.period_range('2014-01', '2016-12', freq='M')
        result = rng - pd.offsets.MonthEnd(5)
        expected = pd.period_range('2013-08', '2016-07', freq='M')
        tm.assert_index_equal(result, expected)

        rng -= pd.offsets.MonthEnd(5)
        tm.assert_index_equal(rng, expected)

    def test_pi_add_offset_n_gt1(self, box):
        # GH#23215
        # add offset to PeriodIndex with freq.n > 1
        per = pd.Period('2016-01', freq='2M')
        pi = pd.PeriodIndex([per])

        expected = pd.PeriodIndex(['2016-03'], freq='2M')
        pi = tm.box_expected(pi, box)
        expected = tm.box_expected(expected, box)

        result = pi + per.freq
        tm.assert_equal(result, expected)

        result = per.freq + pi
        tm.assert_equal(result, expected)

    def test_pi_add_offset_n_gt1_not_divisible(self, box):
        # GH#23215
        # PeriodIndex with freq.n > 1 add offset with offset.n % freq.n != 0

        pi = pd.PeriodIndex(['2016-01'], freq='2M')
        pi = tm.box_expected(pi, box)

        expected = pd.PeriodIndex(['2016-04'], freq='2M')
        expected = tm.box_expected(expected, box)

        result = pi + to_offset('3M')
        tm.assert_equal(result, expected)

        result = to_offset('3M') + pi
        tm.assert_equal(result, expected)

    # ---------------------------------------------------------------
    # __add__/__sub__ with integer arrays

    @pytest.mark.parametrize('box', [np.array, pd.Index])
    @pytest.mark.parametrize('op', [operator.add, ops.radd])
    def test_pi_add_intarray(self, box, op):
        # GH#19959
        pi = pd.PeriodIndex([pd.Period('2015Q1'), pd.Period('NaT')])
        other = box([4, -1])
        result = op(pi, other)
        expected = pd.PeriodIndex([pd.Period('2016Q1'), pd.Period('NaT')])
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('box', [np.array, pd.Index])
    def test_pi_sub_intarray(self, box):
        # GH#19959
        pi = pd.PeriodIndex([pd.Period('2015Q1'), pd.Period('NaT')])
        other = box([4, -1])
        result = pi - other
        expected = pd.PeriodIndex([pd.Period('2014Q1'), pd.Period('NaT')])
        tm.assert_index_equal(result, expected)

        with pytest.raises(TypeError):
            other - pi

    # ---------------------------------------------------------------
    # Timedelta-like (timedelta, timedelta64, Timedelta, Tick)
    # TODO: Some of these are misnomers because of non-Tick DateOffsets

    def test_pi_add_timedeltalike_minute_gt1(self, three_days):
        # GH#23031 adding a time-delta-like offset to a PeriodArray that has
        # minute frequency with n != 1.  A more general case is tested below
        # in test_pi_add_timedeltalike_tick_gt1, but here we write out the
        # expected result more explicitly.
        other = three_days
        rng = pd.period_range('2014-05-01', periods=3, freq='2D')

        expected = pd.PeriodIndex(['2014-05-04', '2014-05-06', '2014-05-08'],
                                  freq='2D')

        result = rng + other
        tm.assert_index_equal(result, expected)

        result = other + rng
        tm.assert_index_equal(result, expected)

        # subtraction
        expected = pd.PeriodIndex(['2014-04-28', '2014-04-30', '2014-05-02'],
                                  freq='2D')
        result = rng - other
        tm.assert_index_equal(result, expected)

        with pytest.raises(TypeError):
            other - rng

    @pytest.mark.parametrize('freqstr', ['5ns', '5us', '5ms',
                                         '5s', '5T', '5h', '5d'])
    def test_pi_add_timedeltalike_tick_gt1(self, three_days, freqstr):
        # GH#23031 adding a time-delta-like offset to a PeriodArray that has
        # tick-like frequency with n != 1
        other = three_days
        rng = pd.period_range('2014-05-01', periods=6, freq=freqstr)

        expected = pd.period_range(rng[0] + other, periods=6, freq=freqstr)

        result = rng + other
        tm.assert_index_equal(result, expected)

        result = other + rng
        tm.assert_index_equal(result, expected)

        # subtraction
        expected = pd.period_range(rng[0] - other, periods=6, freq=freqstr)
        result = rng - other
        tm.assert_index_equal(result, expected)

        with pytest.raises(TypeError):
            other - rng

    def test_pi_add_iadd_timedeltalike_daily(self, three_days):
        # Tick
        other = three_days
        rng = pd.period_range('2014-05-01', '2014-05-15', freq='D')
        expected = pd.period_range('2014-05-04', '2014-05-18', freq='D')

        result = rng + other
        tm.assert_index_equal(result, expected)

        rng += other
        tm.assert_index_equal(rng, expected)

    def test_pi_sub_isub_timedeltalike_daily(self, three_days):
        # Tick-like 3 Days
        other = three_days
        rng = pd.period_range('2014-05-01', '2014-05-15', freq='D')
        expected = pd.period_range('2014-04-28', '2014-05-12', freq='D')

        result = rng - other
        tm.assert_index_equal(result, expected)

        rng -= other
        tm.assert_index_equal(rng, expected)

    def test_pi_add_sub_timedeltalike_freq_mismatch_daily(self, not_daily):
        other = not_daily
        rng = pd.period_range('2014-05-01', '2014-05-15', freq='D')
        msg = 'Input has different freq(=.+)? from PeriodIndex\\(freq=D\\)'
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng + other
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng += other
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng - other
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng -= other

    def test_pi_add_iadd_timedeltalike_hourly(self, two_hours):
        other = two_hours
        rng = pd.period_range('2014-01-01 10:00', '2014-01-05 10:00', freq='H')
        expected = pd.period_range('2014-01-01 12:00', '2014-01-05 12:00',
                                   freq='H')

        result = rng + other
        tm.assert_index_equal(result, expected)

        rng += other
        tm.assert_index_equal(rng, expected)

    def test_pi_add_timedeltalike_mismatched_freq_hourly(self, not_hourly):
        other = not_hourly
        rng = pd.period_range('2014-01-01 10:00', '2014-01-05 10:00', freq='H')
        msg = 'Input has different freq(=.+)? from PeriodIndex\\(freq=H\\)'

        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng + other

        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng += other

    def test_pi_sub_isub_timedeltalike_hourly(self, two_hours):
        other = two_hours
        rng = pd.period_range('2014-01-01 10:00', '2014-01-05 10:00', freq='H')
        expected = pd.period_range('2014-01-01 08:00', '2014-01-05 08:00',
                                   freq='H')

        result = rng - other
        tm.assert_index_equal(result, expected)

        rng -= other
        tm.assert_index_equal(rng, expected)

    def test_add_iadd_timedeltalike_annual(self):
        # offset
        # DateOffset
        rng = pd.period_range('2014', '2024', freq='A')
        result = rng + pd.offsets.YearEnd(5)
        expected = pd.period_range('2019', '2029', freq='A')
        tm.assert_index_equal(result, expected)
        rng += pd.offsets.YearEnd(5)
        tm.assert_index_equal(rng, expected)

    def test_pi_add_sub_timedeltalike_freq_mismatch_annual(self,
                                                           mismatched_freq):
        other = mismatched_freq
        rng = pd.period_range('2014', '2024', freq='A')
        msg = ('Input has different freq(=.+)? '
               'from PeriodIndex\\(freq=A-DEC\\)')
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng + other
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng += other
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng - other
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng -= other

    def test_pi_add_iadd_timedeltalike_M(self):
        rng = pd.period_range('2014-01', '2016-12', freq='M')
        expected = pd.period_range('2014-06', '2017-05', freq='M')

        result = rng + pd.offsets.MonthEnd(5)
        tm.assert_index_equal(result, expected)

        rng += pd.offsets.MonthEnd(5)
        tm.assert_index_equal(rng, expected)

    def test_pi_add_sub_timedeltalike_freq_mismatch_monthly(self,
                                                            mismatched_freq):
        other = mismatched_freq
        rng = pd.period_range('2014-01', '2016-12', freq='M')
        msg = 'Input has different freq(=.+)? from PeriodIndex\\(freq=M\\)'
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng + other
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng += other
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng - other
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            rng -= other


class TestPeriodSeriesArithmetic(object):
    def test_ops_series_timedelta(self):
        # GH 13043
        ser = pd.Series([pd.Period('2015-01-01', freq='D'),
                         pd.Period('2015-01-02', freq='D')], name='xxx')
        assert ser.dtype == object

        expected = pd.Series([pd.Period('2015-01-02', freq='D'),
                              pd.Period('2015-01-03', freq='D')], name='xxx')

        result = ser + pd.Timedelta('1 days')
        tm.assert_series_equal(result, expected)

        result = pd.Timedelta('1 days') + ser
        tm.assert_series_equal(result, expected)

        result = ser + pd.tseries.offsets.Day()
        tm.assert_series_equal(result, expected)

        result = pd.tseries.offsets.Day() + ser
        tm.assert_series_equal(result, expected)

    def test_ops_series_period(self):
        # GH 13043
        ser = pd.Series([pd.Period('2015-01-01', freq='D'),
                         pd.Period('2015-01-02', freq='D')], name='xxx')
        assert ser.dtype == object

        per = pd.Period('2015-01-10', freq='D')
        off = per.freq
        # dtype will be object because of original dtype
        expected = pd.Series([9 * off, 8 * off], name='xxx', dtype=object)
        tm.assert_series_equal(per - ser, expected)
        tm.assert_series_equal(ser - per, -1 * expected)

        s2 = pd.Series([pd.Period('2015-01-05', freq='D'),
                        pd.Period('2015-01-04', freq='D')], name='xxx')
        assert s2.dtype == object

        expected = pd.Series([4 * off, 2 * off], name='xxx', dtype=object)
        tm.assert_series_equal(s2 - ser, expected)
        tm.assert_series_equal(ser - s2, -1 * expected)


class TestPeriodIndexSeriesMethods(object):
    """ Test PeriodIndex and Period Series Ops consistency """

    def _check(self, values, func, expected):
        idx = pd.PeriodIndex(values)
        result = func(idx)
        tm.assert_equal(result, expected)

        ser = pd.Series(values)
        result = func(ser)

        exp = pd.Series(expected, name=values.name)
        tm.assert_series_equal(result, exp)

    def test_pi_ops(self):
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03', '2011-04'],
                          freq='M', name='idx')

        expected = PeriodIndex(['2011-03', '2011-04', '2011-05', '2011-06'],
                               freq='M', name='idx')
        self._check(idx, lambda x: x + 2, expected)
        self._check(idx, lambda x: 2 + x, expected)

        self._check(idx + 2, lambda x: x - 2, idx)
        result = idx - Period('2011-01', freq='M')
        off = idx.freq
        exp = pd.Index([0 * off, 1 * off, 2 * off, 3 * off], name='idx')
        tm.assert_index_equal(result, exp)

        result = Period('2011-01', freq='M') - idx
        exp = pd.Index([0 * off, -1 * off, -2 * off, -3 * off], name='idx')
        tm.assert_index_equal(result, exp)

    @pytest.mark.parametrize('ng', ["str", 1.5])
    def test_pi_ops_errors(self, ng, box):
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03', '2011-04'],
                          freq='M', name='idx')
        obj = tm.box_expected(idx, box)

        msg = r"unsupported operand type\(s\)"

        with tm.assert_raises_regex(TypeError, msg):
            obj + ng

        with pytest.raises(TypeError):
            # error message differs between PY2 and 3
            ng + obj

        with tm.assert_raises_regex(TypeError, msg):
            obj - ng

        with pytest.raises(TypeError):
            np.add(obj, ng)

        with pytest.raises(TypeError):
            np.add(ng, obj)

        with pytest.raises(TypeError):
            np.subtract(obj, ng)

        with pytest.raises(TypeError):
            np.subtract(ng, obj)

    def test_pi_ops_nat(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='M', name='idx')
        expected = PeriodIndex(['2011-03', '2011-04', 'NaT', '2011-06'],
                               freq='M', name='idx')
        self._check(idx, lambda x: x + 2, expected)
        self._check(idx, lambda x: 2 + x, expected)
        self._check(idx, lambda x: np.add(x, 2), expected)

        self._check(idx + 2, lambda x: x - 2, idx)
        self._check(idx + 2, lambda x: np.subtract(x, 2), idx)

        # freq with mult
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='2M', name='idx')
        expected = PeriodIndex(['2011-07', '2011-08', 'NaT', '2011-10'],
                               freq='2M', name='idx')
        self._check(idx, lambda x: x + 3, expected)
        self._check(idx, lambda x: 3 + x, expected)
        self._check(idx, lambda x: np.add(x, 3), expected)

        self._check(idx + 3, lambda x: x - 3, idx)
        self._check(idx + 3, lambda x: np.subtract(x, 3), idx)

    def test_pi_ops_array_int(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='M', name='idx')
        f = lambda x: x + np.array([1, 2, 3, 4])
        exp = PeriodIndex(['2011-02', '2011-04', 'NaT', '2011-08'],
                          freq='M', name='idx')
        self._check(idx, f, exp)

        f = lambda x: np.add(x, np.array([4, -1, 1, 2]))
        exp = PeriodIndex(['2011-05', '2011-01', 'NaT', '2011-06'],
                          freq='M', name='idx')
        self._check(idx, f, exp)

        f = lambda x: x - np.array([1, 2, 3, 4])
        exp = PeriodIndex(['2010-12', '2010-12', 'NaT', '2010-12'],
                          freq='M', name='idx')
        self._check(idx, f, exp)

        f = lambda x: np.subtract(x, np.array([3, 2, 3, -2]))
        exp = PeriodIndex(['2010-10', '2010-12', 'NaT', '2011-06'],
                          freq='M', name='idx')
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
        ser = pd.Series(idx)

        # Series op is applied per Period instance, thus error is raised
        # from Period
        msg_idx = r"Input has different freq from PeriodIndex\(freq=D\)"
        msg_s = r"Input cannot be converted to Period\(freq=D\)"
        for obj, msg in [(idx, msg_idx), (ser, msg_s)]:
            with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
                obj + pd.offsets.Hour(2)

            with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
                pd.offsets.Hour(2) + obj

            with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
                obj - pd.offsets.Hour(2)

    def test_pi_sub_period(self):
        # GH 13071
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03', '2011-04'],
                          freq='M', name='idx')

        result = idx - pd.Period('2012-01', freq='M')
        off = idx.freq
        exp = pd.Index([-12 * off, -11 * off, -10 * off, -9 * off], name='idx')
        tm.assert_index_equal(result, exp)

        result = np.subtract(idx, pd.Period('2012-01', freq='M'))
        tm.assert_index_equal(result, exp)

        result = pd.Period('2012-01', freq='M') - idx
        exp = pd.Index([12 * off, 11 * off, 10 * off, 9 * off], name='idx')
        tm.assert_index_equal(result, exp)

        result = np.subtract(pd.Period('2012-01', freq='M'), idx)
        tm.assert_index_equal(result, exp)

        exp = pd.TimedeltaIndex([np.nan, np.nan, np.nan, np.nan], name='idx')
        tm.assert_index_equal(idx - pd.Period('NaT', freq='M'), exp)
        tm.assert_index_equal(pd.Period('NaT', freq='M') - idx, exp)

    def test_pi_sub_pdnat(self):
        # GH 13071
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'],
                          freq='M', name='idx')
        exp = pd.TimedeltaIndex([pd.NaT] * 4, name='idx')
        tm.assert_index_equal(pd.NaT - idx, exp)
        tm.assert_index_equal(idx - pd.NaT, exp)

    def test_pi_sub_period_nat(self):
        # GH 13071
        idx = PeriodIndex(['2011-01', 'NaT', '2011-03', '2011-04'],
                          freq='M', name='idx')

        result = idx - pd.Period('2012-01', freq='M')
        off = idx.freq
        exp = pd.Index([-12 * off, pd.NaT, -10 * off, -9 * off], name='idx')
        tm.assert_index_equal(result, exp)

        result = pd.Period('2012-01', freq='M') - idx
        exp = pd.Index([12 * off, pd.NaT, 10 * off, 9 * off], name='idx')
        tm.assert_index_equal(result, exp)

        exp = pd.TimedeltaIndex([np.nan, np.nan, np.nan, np.nan], name='idx')
        tm.assert_index_equal(idx - pd.Period('NaT', freq='M'), exp)
        tm.assert_index_equal(pd.Period('NaT', freq='M') - idx, exp)
