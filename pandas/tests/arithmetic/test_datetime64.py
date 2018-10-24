# -*- coding: utf-8 -*-
# Arithmetc tests for DataFrame/Series/Index/Array classes that should
# behave identically.
# Specifically for datetime64 and datetime64tz dtypes
import operator
from datetime import datetime, timedelta
import warnings
from itertools import product, starmap

import numpy as np
import pytest
import pytz

import pandas as pd
import pandas.util.testing as tm

from pandas.compat.numpy import np_datetime64_compat
from pandas.errors import PerformanceWarning, NullFrequencyError

from pandas._libs.tslibs.conversion import localize_pydatetime
from pandas._libs.tslibs.offsets import shift_months

from pandas.core import ops

from pandas import (
    Timestamp, Timedelta, Period, Series, date_range, NaT,
    DatetimeIndex, TimedeltaIndex)


# ------------------------------------------------------------------
# Comparisons

class TestDatetime64DataFrameComparison(object):
    @pytest.mark.parametrize('timestamps', [
        [pd.Timestamp('2012-01-01 13:00:00+00:00')] * 2,
        [pd.Timestamp('2012-01-01 13:00:00')] * 2])
    def test_tz_aware_scalar_comparison(self, timestamps):
        # GH#15966
        df = pd.DataFrame({'test': timestamps})
        expected = pd.DataFrame({'test': [False, False]})
        tm.assert_frame_equal(df == -1, expected)

    def test_dt64_nat_comparison(self):
        # GH#22242, GH#22163 DataFrame considered NaT == ts incorrectly
        ts = pd.Timestamp.now()
        df = pd.DataFrame([ts, pd.NaT])
        expected = pd.DataFrame([True, False])

        result = df == ts
        tm.assert_frame_equal(result, expected)


class TestDatetime64SeriesComparison(object):
    # TODO: moved from tests.series.test_operators; needs cleanup
    def test_comparison_invalid(self):
        # GH#4968
        # invalid date/int comparisons
        ser = Series(range(5))
        ser2 = Series(pd.date_range('20010101', periods=5))

        for (x, y) in [(ser, ser2), (ser2, ser)]:

            result = x == y
            expected = Series([False] * 5)
            tm.assert_series_equal(result, expected)

            result = x != y
            expected = Series([True] * 5)
            tm.assert_series_equal(result, expected)

            with pytest.raises(TypeError):
                x >= y
            with pytest.raises(TypeError):
                x > y
            with pytest.raises(TypeError):
                x < y
            with pytest.raises(TypeError):
                x <= y

    @pytest.mark.parametrize('data', [
        [Timestamp('2011-01-01'), NaT, Timestamp('2011-01-03')],
        [Timedelta('1 days'), NaT, Timedelta('3 days')],
        [Period('2011-01', freq='M'), NaT, Period('2011-03', freq='M')]
    ])
    @pytest.mark.parametrize('dtype', [None, object])
    def test_nat_comparisons_scalar(self, dtype, data):
        left = Series(data, dtype=dtype)

        expected = Series([False, False, False])
        tm.assert_series_equal(left == NaT, expected)
        tm.assert_series_equal(NaT == left, expected)

        expected = Series([True, True, True])
        tm.assert_series_equal(left != NaT, expected)
        tm.assert_series_equal(NaT != left, expected)

        expected = Series([False, False, False])
        tm.assert_series_equal(left < NaT, expected)
        tm.assert_series_equal(NaT > left, expected)
        tm.assert_series_equal(left <= NaT, expected)
        tm.assert_series_equal(NaT >= left, expected)

        tm.assert_series_equal(left > NaT, expected)
        tm.assert_series_equal(NaT < left, expected)
        tm.assert_series_equal(left >= NaT, expected)
        tm.assert_series_equal(NaT <= left, expected)

    def test_series_comparison_scalars(self):
        series = Series(date_range('1/1/2000', periods=10))

        val = datetime(2000, 1, 4)
        result = series > val
        expected = Series([x > val for x in series])
        tm.assert_series_equal(result, expected)

        val = series[5]
        result = series > val
        expected = Series([x > val for x in series])
        tm.assert_series_equal(result, expected)

    def test_dt64_ser_cmp_date_warning(self):
        # https://github.com/pandas-dev/pandas/issues/21359
        # Remove this test and enble invalid test below
        ser = pd.Series(pd.date_range('20010101', periods=10), name='dates')
        date = ser.iloc[0].to_pydatetime().date()

        with tm.assert_produces_warning(FutureWarning) as m:
            result = ser == date
        expected = pd.Series([True] + [False] * 9, name='dates')
        tm.assert_series_equal(result, expected)
        assert "Comparing Series of datetimes " in str(m[0].message)
        assert "will not compare equal" in str(m[0].message)

        with tm.assert_produces_warning(FutureWarning) as m:
            result = ser != date
        tm.assert_series_equal(result, ~expected)
        assert "will not compare equal" in str(m[0].message)

        with tm.assert_produces_warning(FutureWarning) as m:
            result = ser <= date
        tm.assert_series_equal(result, expected)
        assert "a TypeError will be raised" in str(m[0].message)

        with tm.assert_produces_warning(FutureWarning) as m:
            result = ser < date
        tm.assert_series_equal(result, pd.Series([False] * 10, name='dates'))
        assert "a TypeError will be raised" in str(m[0].message)

        with tm.assert_produces_warning(FutureWarning) as m:
            result = ser >= date
        tm.assert_series_equal(result, pd.Series([True] * 10, name='dates'))
        assert "a TypeError will be raised" in str(m[0].message)

        with tm.assert_produces_warning(FutureWarning) as m:
            result = ser > date
        tm.assert_series_equal(result, pd.Series([False] + [True] * 9,
                                                 name='dates'))
        assert "a TypeError will be raised" in str(m[0].message)

    @pytest.mark.skip(reason="GH#21359")
    def test_dt64ser_cmp_date_invalid(self):
        # GH#19800 datetime.date comparison raises to
        # match DatetimeIndex/Timestamp.  This also matches the behavior
        # of stdlib datetime.datetime
        ser = pd.Series(pd.date_range('20010101', periods=10), name='dates')
        date = ser.iloc[0].to_pydatetime().date()
        assert not (ser == date).any()
        assert (ser != date).all()
        with pytest.raises(TypeError):
            ser > date
        with pytest.raises(TypeError):
            ser < date
        with pytest.raises(TypeError):
            ser >= date
        with pytest.raises(TypeError):
            ser <= date

    def test_dt64ser_cmp_period_scalar(self):
        ser = Series(pd.period_range('2000-01-01', periods=10, freq='D'))
        val = Period('2000-01-04', freq='D')
        result = ser > val
        expected = Series([x > val for x in ser])
        tm.assert_series_equal(result, expected)

        val = ser[5]
        result = ser > val
        expected = Series([x > val for x in ser])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("left,right", [
        ("lt", "gt"),
        ("le", "ge"),
        ("eq", "eq"),
        ("ne", "ne"),
    ])
    def test_timestamp_compare_series(self, left, right):
        # see gh-4982
        # Make sure we can compare Timestamps on the right AND left hand side.
        ser = pd.Series(pd.date_range("20010101", periods=10), name="dates")
        s_nat = ser.copy(deep=True)

        ser[0] = pd.Timestamp("nat")
        ser[3] = pd.Timestamp("nat")

        left_f = getattr(operator, left)
        right_f = getattr(operator, right)

        # No NaT
        expected = left_f(ser, pd.Timestamp("20010109"))
        result = right_f(pd.Timestamp("20010109"), ser)
        tm.assert_series_equal(result, expected)

        # NaT
        expected = left_f(ser, pd.Timestamp("nat"))
        result = right_f(pd.Timestamp("nat"), ser)
        tm.assert_series_equal(result, expected)

        # Compare to Timestamp with series containing NaT
        expected = left_f(s_nat, pd.Timestamp("20010109"))
        result = right_f(pd.Timestamp("20010109"), s_nat)
        tm.assert_series_equal(result, expected)

        # Compare to NaT with series containing NaT
        expected = left_f(s_nat, pd.Timestamp("nat"))
        result = right_f(pd.Timestamp("nat"), s_nat)
        tm.assert_series_equal(result, expected)

    def test_timestamp_equality(self):
        # GH#11034
        ser = pd.Series([pd.Timestamp('2000-01-29 01:59:00'), 'NaT'])
        result = ser != ser
        tm.assert_series_equal(result, pd.Series([False, True]))
        result = ser != ser[0]
        tm.assert_series_equal(result, pd.Series([False, True]))
        result = ser != ser[1]
        tm.assert_series_equal(result, pd.Series([True, True]))

        result = ser == ser
        tm.assert_series_equal(result, pd.Series([True, False]))
        result = ser == ser[0]
        tm.assert_series_equal(result, pd.Series([True, False]))
        result = ser == ser[1]
        tm.assert_series_equal(result, pd.Series([False, False]))


class TestDatetimeIndexComparisons(object):
    @pytest.mark.parametrize('other', [datetime(2016, 1, 1),
                                       Timestamp('2016-01-01'),
                                       np.datetime64('2016-01-01')])
    def test_dti_cmp_datetimelike(self, other, tz_naive_fixture):
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=2, tz=tz)
        if tz is not None:
            if isinstance(other, np.datetime64):
                # no tzaware version available
                return
            other = localize_pydatetime(other, dti.tzinfo)

        result = dti == other
        expected = np.array([True, False])
        tm.assert_numpy_array_equal(result, expected)

        result = dti > other
        expected = np.array([False, True])
        tm.assert_numpy_array_equal(result, expected)

        result = dti >= other
        expected = np.array([True, True])
        tm.assert_numpy_array_equal(result, expected)

        result = dti < other
        expected = np.array([False, False])
        tm.assert_numpy_array_equal(result, expected)

        result = dti <= other
        expected = np.array([True, False])
        tm.assert_numpy_array_equal(result, expected)

    def dti_cmp_non_datetime(self, tz_naive_fixture):
        # GH#19301 by convention datetime.date is not considered comparable
        # to Timestamp or DatetimeIndex.  This may change in the future.
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=2, tz=tz)

        other = datetime(2016, 1, 1).date()
        assert not (dti == other).any()
        assert (dti != other).all()
        with pytest.raises(TypeError):
            dti < other
        with pytest.raises(TypeError):
            dti <= other
        with pytest.raises(TypeError):
            dti > other
        with pytest.raises(TypeError):
            dti >= other

    @pytest.mark.parametrize('other', [None, np.nan, pd.NaT])
    def test_dti_eq_null_scalar(self, other, tz_naive_fixture):
        # GH#19301
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=2, tz=tz)
        assert not (dti == other).any()

    @pytest.mark.parametrize('other', [None, np.nan, pd.NaT])
    def test_dti_ne_null_scalar(self, other, tz_naive_fixture):
        # GH#19301
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=2, tz=tz)
        assert (dti != other).all()

    @pytest.mark.parametrize('other', [None, np.nan])
    def test_dti_cmp_null_scalar_inequality(self, tz_naive_fixture, other):
        # GH#19301
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=2, tz=tz)

        with pytest.raises(TypeError):
            dti < other
        with pytest.raises(TypeError):
            dti <= other
        with pytest.raises(TypeError):
            dti > other
        with pytest.raises(TypeError):
            dti >= other

    def test_dti_cmp_nat(self):
        left = pd.DatetimeIndex([pd.Timestamp('2011-01-01'), pd.NaT,
                                 pd.Timestamp('2011-01-03')])
        right = pd.DatetimeIndex([pd.NaT, pd.NaT, pd.Timestamp('2011-01-03')])

        for lhs, rhs in [(left, right),
                         (left.astype(object), right.astype(object))]:
            result = rhs == lhs
            expected = np.array([False, False, True])
            tm.assert_numpy_array_equal(result, expected)

            result = lhs != rhs
            expected = np.array([True, True, False])
            tm.assert_numpy_array_equal(result, expected)

            expected = np.array([False, False, False])
            tm.assert_numpy_array_equal(lhs == pd.NaT, expected)
            tm.assert_numpy_array_equal(pd.NaT == rhs, expected)

            expected = np.array([True, True, True])
            tm.assert_numpy_array_equal(lhs != pd.NaT, expected)
            tm.assert_numpy_array_equal(pd.NaT != lhs, expected)

            expected = np.array([False, False, False])
            tm.assert_numpy_array_equal(lhs < pd.NaT, expected)
            tm.assert_numpy_array_equal(pd.NaT > lhs, expected)

    def test_dti_cmp_nat_behaves_like_float_cmp_nan(self):
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

        cases = [(fidx1, fidx2), (didx1, didx2), (didx1, darr)]

        # Check pd.NaT is handles as the same as np.nan
        with tm.assert_produces_warning(None):
            for idx1, idx2 in cases:

                result = idx1 < idx2
                expected = np.array([True, False, False, False, True, False])
                tm.assert_numpy_array_equal(result, expected)

                result = idx2 > idx1
                expected = np.array([True, False, False, False, True, False])
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 <= idx2
                expected = np.array([True, False, False, False, True, True])
                tm.assert_numpy_array_equal(result, expected)

                result = idx2 >= idx1
                expected = np.array([True, False, False, False, True, True])
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 == idx2
                expected = np.array([False, False, False, False, False, True])
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 != idx2
                expected = np.array([True, True, True, True, True, False])
                tm.assert_numpy_array_equal(result, expected)

        with tm.assert_produces_warning(None):
            for idx1, val in [(fidx1, np.nan), (didx1, pd.NaT)]:
                result = idx1 < val
                expected = np.array([False, False, False, False, False, False])
                tm.assert_numpy_array_equal(result, expected)
                result = idx1 > val
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 <= val
                tm.assert_numpy_array_equal(result, expected)
                result = idx1 >= val
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 == val
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 != val
                expected = np.array([True, True, True, True, True, True])
                tm.assert_numpy_array_equal(result, expected)

        # Check pd.NaT is handles as the same as np.nan
        with tm.assert_produces_warning(None):
            for idx1, val in [(fidx1, 3), (didx1, datetime(2014, 3, 1))]:
                result = idx1 < val
                expected = np.array([True, False, False, False, False, False])
                tm.assert_numpy_array_equal(result, expected)
                result = idx1 > val
                expected = np.array([False, False, False, False, True, True])
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 <= val
                expected = np.array([True, False, True, False, False, False])
                tm.assert_numpy_array_equal(result, expected)
                result = idx1 >= val
                expected = np.array([False, False, True, False, True, True])
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 == val
                expected = np.array([False, False, True, False, False, False])
                tm.assert_numpy_array_equal(result, expected)

                result = idx1 != val
                expected = np.array([True, True, False, True, True, True])
                tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize('op', [operator.eq, operator.ne,
                                    operator.gt, operator.ge,
                                    operator.lt, operator.le])
    def test_comparison_tzawareness_compat(self, op):
        # GH#18162
        dr = pd.date_range('2016-01-01', periods=6)
        dz = dr.tz_localize('US/Pacific')

        with pytest.raises(TypeError):
            op(dr, dz)
        with pytest.raises(TypeError):
            op(dr, list(dz))
        with pytest.raises(TypeError):
            op(dz, dr)
        with pytest.raises(TypeError):
            op(dz, list(dr))

        # Check that there isn't a problem aware-aware and naive-naive do not
        # raise
        assert (dr == dr).all()
        assert (dr == list(dr)).all()
        assert (dz == dz).all()
        assert (dz == list(dz)).all()

        # Check comparisons against scalar Timestamps
        ts = pd.Timestamp('2000-03-14 01:59')
        ts_tz = pd.Timestamp('2000-03-14 01:59', tz='Europe/Amsterdam')

        assert (dr > ts).all()
        with pytest.raises(TypeError):
            op(dr, ts_tz)

        assert (dz > ts_tz).all()
        with pytest.raises(TypeError):
            op(dz, ts)

        # GH#12601: Check comparison against Timestamps and DatetimeIndex
        with pytest.raises(TypeError):
            op(ts, dz)

    @pytest.mark.parametrize('op', [operator.eq, operator.ne,
                                    operator.gt, operator.ge,
                                    operator.lt, operator.le])
    @pytest.mark.parametrize('other', [datetime(2016, 1, 1),
                                       Timestamp('2016-01-01'),
                                       np.datetime64('2016-01-01')])
    def test_scalar_comparison_tzawareness(self, op, other, tz_aware_fixture):
        tz = tz_aware_fixture
        dti = pd.date_range('2016-01-01', periods=2, tz=tz)
        with pytest.raises(TypeError):
            op(dti, other)
        with pytest.raises(TypeError):
            op(other, dti)

    @pytest.mark.parametrize('op', [operator.eq, operator.ne,
                                    operator.gt, operator.ge,
                                    operator.lt, operator.le])
    def test_nat_comparison_tzawareness(self, op):
        # GH#19276
        # tzaware DatetimeIndex should not raise when compared to NaT
        dti = pd.DatetimeIndex(['2014-01-01', pd.NaT, '2014-03-01', pd.NaT,
                                '2014-05-01', '2014-07-01'])
        expected = np.array([op == operator.ne] * len(dti))
        result = op(dti, pd.NaT)
        tm.assert_numpy_array_equal(result, expected)

        result = op(dti.tz_localize('US/Pacific'), pd.NaT)
        tm.assert_numpy_array_equal(result, expected)

    def test_dti_cmp_str(self, tz_naive_fixture):
        # GH#22074
        # regardless of tz, we expect these comparisons are valid
        tz = tz_naive_fixture
        rng = date_range('1/1/2000', periods=10, tz=tz)
        other = '1/1/2000'

        result = rng == other
        expected = np.array([True] + [False] * 9)
        tm.assert_numpy_array_equal(result, expected)

        result = rng != other
        expected = np.array([False] + [True] * 9)
        tm.assert_numpy_array_equal(result, expected)

        result = rng < other
        expected = np.array([False] * 10)
        tm.assert_numpy_array_equal(result, expected)

        result = rng <= other
        expected = np.array([True] + [False] * 9)
        tm.assert_numpy_array_equal(result, expected)

        result = rng > other
        expected = np.array([False] + [True] * 9)
        tm.assert_numpy_array_equal(result, expected)

        result = rng >= other
        expected = np.array([True] * 10)
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize('other', ['foo', 99, 4.0,
                                       object(), timedelta(days=2)])
    def test_dti_cmp_scalar_invalid(self, other, tz_naive_fixture):
        # GH#22074
        tz = tz_naive_fixture
        rng = date_range('1/1/2000', periods=10, tz=tz)

        result = rng == other
        expected = np.array([False] * 10)
        tm.assert_numpy_array_equal(result, expected)

        result = rng != other
        expected = np.array([True] * 10)
        tm.assert_numpy_array_equal(result, expected)

        with pytest.raises(TypeError):
            rng < other
        with pytest.raises(TypeError):
            rng <= other
        with pytest.raises(TypeError):
            rng > other
        with pytest.raises(TypeError):
            rng >= other

    def test_dti_cmp_list(self):
        rng = date_range('1/1/2000', periods=10)

        result = rng == list(rng)
        expected = rng == rng
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize('other', [
        pd.timedelta_range('1D', periods=10),
        pd.timedelta_range('1D', periods=10).to_series(),
        pd.timedelta_range('1D', periods=10).asi8.view('m8[ns]')
    ], ids=lambda x: type(x).__name__)
    def test_dti_cmp_tdi_tzawareness(self, other):
        # GH#22074
        # reversion test that we _don't_ call _assert_tzawareness_compat
        # when comparing against TimedeltaIndex
        dti = date_range('2000-01-01', periods=10, tz='Asia/Tokyo')

        result = dti == other
        expected = np.array([False] * 10)
        tm.assert_numpy_array_equal(result, expected)

        result = dti != other
        expected = np.array([True] * 10)
        tm.assert_numpy_array_equal(result, expected)

        with pytest.raises(TypeError):
            dti < other
        with pytest.raises(TypeError):
            dti <= other
        with pytest.raises(TypeError):
            dti > other
        with pytest.raises(TypeError):
            dti >= other

    def test_dti_cmp_object_dtype(self):
        # GH#22074
        dti = date_range('2000-01-01', periods=10, tz='Asia/Tokyo')

        other = dti.astype('O')

        result = dti == other
        expected = np.array([True] * 10)
        tm.assert_numpy_array_equal(result, expected)

        other = dti.tz_localize(None)
        with pytest.raises(TypeError):
            # tzawareness failure
            dti != other

        other = np.array(list(dti[:5]) + [Timedelta(days=1)] * 5)
        result = dti == other
        expected = np.array([True] * 5 + [False] * 5)
        tm.assert_numpy_array_equal(result, expected)

        with pytest.raises(TypeError):
            dti >= other


# ------------------------------------------------------------------
# Arithmetic

class TestFrameArithmetic(object):
    def test_dt64arr_sub_dtscalar(self, box):
        # GH#8554, GH#22163 DataFrame op should _not_ return dt64 dtype
        idx = pd.date_range('2013-01-01', periods=3)
        idx = tm.box_expected(idx, box)

        ts = pd.Timestamp('2013-01-01')
        # TODO: parametrize over scalar types

        expected = pd.TimedeltaIndex(['0 Days', '1 Day', '2 Days'])
        expected = tm.box_expected(expected, box)

        result = idx - ts
        tm.assert_equal(result, expected)

    def test_df_sub_datetime64_not_ns(self):
        # GH#7996, GH#22163 ensure non-nano datetime64 is converted to nano
        df = pd.DataFrame(pd.date_range('20130101', periods=3))
        dt64 = np.datetime64('2013-01-01')
        assert dt64.dtype == 'datetime64[D]'
        res = df - dt64
        expected = pd.DataFrame([pd.Timedelta(days=0), pd.Timedelta(days=1),
                                 pd.Timedelta(days=2)])
        tm.assert_frame_equal(res, expected)


class TestTimestampSeriesArithmetic(object):

    def test_timestamp_sub_series(self):
        ser = pd.Series(pd.date_range('2014-03-17', periods=2, freq='D',
                                      tz='US/Eastern'))
        ts = ser[0]

        delta_series = pd.Series([np.timedelta64(0, 'D'),
                                  np.timedelta64(1, 'D')])
        tm.assert_series_equal(ser - ts, delta_series)
        tm.assert_series_equal(ts - ser, -delta_series)

    def test_dt64ser_sub_datetime_dtype(self):
        ts = Timestamp(datetime(1993, 1, 7, 13, 30, 00))
        dt = datetime(1993, 6, 22, 13, 30)
        ser = Series([ts])
        result = pd.to_timedelta(np.abs(ser - dt))
        assert result.dtype == 'timedelta64[ns]'

    # -------------------------------------------------------------
    # TODO: This next block of tests came from tests.series.test_operators,
    # needs to be de-duplicated and parametrized over `box` classes

    @pytest.mark.parametrize('klass', [Series, pd.Index])
    def test_sub_datetime64_not_ns(self, klass):
        # GH#7996
        dt64 = np.datetime64('2013-01-01')
        assert dt64.dtype == 'datetime64[D]'

        obj = klass(date_range('20130101', periods=3))
        res = obj - dt64
        expected = klass([Timedelta(days=0), Timedelta(days=1),
                          Timedelta(days=2)])
        tm.assert_equal(res, expected)

        res = dt64 - obj
        tm.assert_equal(res, -expected)

    def test_sub_single_tz(self):
        # GH12290
        s1 = Series([pd.Timestamp('2016-02-10', tz='America/Sao_Paulo')])
        s2 = Series([pd.Timestamp('2016-02-08', tz='America/Sao_Paulo')])
        result = s1 - s2
        expected = Series([Timedelta('2days')])
        tm.assert_series_equal(result, expected)
        result = s2 - s1
        expected = Series([Timedelta('-2days')])
        tm.assert_series_equal(result, expected)

    def test_dt64tz_series_sub_dtitz(self):
        # GH#19071 subtracting tzaware DatetimeIndex from tzaware Series
        # (with same tz) raises, fixed by #19024
        dti = pd.date_range('1999-09-30', periods=10, tz='US/Pacific')
        ser = pd.Series(dti)
        expected = pd.Series(pd.TimedeltaIndex(['0days'] * 10))

        res = dti - ser
        tm.assert_series_equal(res, expected)
        res = ser - dti
        tm.assert_series_equal(res, expected)

    def test_sub_datetime_compat(self):
        # see gh-14088
        s = Series([datetime(2016, 8, 23, 12, tzinfo=pytz.utc), pd.NaT])
        dt = datetime(2016, 8, 22, 12, tzinfo=pytz.utc)
        exp = Series([Timedelta('1 days'), pd.NaT])
        tm.assert_series_equal(s - dt, exp)
        tm.assert_series_equal(s - Timestamp(dt), exp)

    def test_dt64_series_addsub_timedelta(self):
        # scalar timedeltas/np.timedelta64 objects
        # operate with np.timedelta64 correctly
        s = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])

        result = s + np.timedelta64(1, 's')
        result2 = np.timedelta64(1, 's') + s
        expected = Series([Timestamp('20130101 9:01:01'),
                           Timestamp('20130101 9:02:01')])
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result2, expected)

        result = s + np.timedelta64(5, 'ms')
        result2 = np.timedelta64(5, 'ms') + s
        expected = Series([Timestamp('20130101 9:01:00.005'),
                           Timestamp('20130101 9:02:00.005')])
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result2, expected)

    def test_dt64_series_add_tick_DateOffset(self):
        # GH 4532
        # operate with pd.offsets
        ser = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])
        expected = Series([Timestamp('20130101 9:01:05'),
                           Timestamp('20130101 9:02:05')])

        result = ser + pd.offsets.Second(5)
        tm.assert_series_equal(result, expected)

        result2 = pd.offsets.Second(5) + ser
        tm.assert_series_equal(result2, expected)

    def test_dt64_series_sub_tick_DateOffset(self):
        # GH 4532
        # operate with pd.offsets
        ser = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])
        expected = Series([Timestamp('20130101 9:00:55'),
                           Timestamp('20130101 9:01:55')])

        result = ser - pd.offsets.Second(5)
        tm.assert_series_equal(result, expected)

        result2 = -pd.offsets.Second(5) + ser
        tm.assert_series_equal(result2, expected)

        with pytest.raises(TypeError):
            pd.offsets.Second(5) - ser

    @pytest.mark.parametrize('cls_name', ['Day', 'Hour', 'Minute', 'Second',
                                          'Milli', 'Micro', 'Nano'])
    def test_dt64_series_add_tick_DateOffset_smoke(self, cls_name):
        # GH 4532
        # smoke tests for valid DateOffsets
        ser = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])

        offset_cls = getattr(pd.offsets, cls_name)
        ser + offset_cls(5)
        offset_cls(5) + ser

    def test_dt64_series_add_mixed_tick_DateOffset(self):
        # GH 4532
        # operate with pd.offsets
        s = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])

        result = s + pd.offsets.Milli(5)
        result2 = pd.offsets.Milli(5) + s
        expected = Series([Timestamp('20130101 9:01:00.005'),
                           Timestamp('20130101 9:02:00.005')])
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result2, expected)

        result = s + pd.offsets.Minute(5) + pd.offsets.Milli(5)
        expected = Series([Timestamp('20130101 9:06:00.005'),
                           Timestamp('20130101 9:07:00.005')])
        tm.assert_series_equal(result, expected)

    def test_dt64_series_sub_NaT(self):
        # GH#18808
        dti = pd.DatetimeIndex([pd.NaT, pd.Timestamp('19900315')])
        ser = pd.Series(dti)
        res = ser - pd.NaT
        expected = pd.Series([pd.NaT, pd.NaT], dtype='timedelta64[ns]')
        tm.assert_series_equal(res, expected)

        dti_tz = dti.tz_localize('Asia/Tokyo')
        ser_tz = pd.Series(dti_tz)
        res = ser_tz - pd.NaT
        expected = pd.Series([pd.NaT, pd.NaT], dtype='timedelta64[ns]')
        tm.assert_series_equal(res, expected)

    def test_dt64_series_arith_overflow(self):
        # GH#12534, fixed by #19024
        dt = pd.Timestamp('1700-01-31')
        td = pd.Timedelta('20000 Days')
        dti = pd.date_range('1949-09-30', freq='100Y', periods=4)
        ser = pd.Series(dti)
        with pytest.raises(OverflowError):
            ser - dt
        with pytest.raises(OverflowError):
            dt - ser
        with pytest.raises(OverflowError):
            ser + td
        with pytest.raises(OverflowError):
            td + ser

        ser.iloc[-1] = pd.NaT
        expected = pd.Series(['2004-10-03', '2104-10-04', '2204-10-04', 'NaT'],
                             dtype='datetime64[ns]')
        res = ser + td
        tm.assert_series_equal(res, expected)
        res = td + ser
        tm.assert_series_equal(res, expected)

        ser.iloc[1:] = pd.NaT
        expected = pd.Series(['91279 Days', 'NaT', 'NaT', 'NaT'],
                             dtype='timedelta64[ns]')
        res = ser - dt
        tm.assert_series_equal(res, expected)
        res = dt - ser
        tm.assert_series_equal(res, -expected)

    def test_datetime64_ops_nat(self):
        # GH 11349
        datetime_series = Series([NaT, Timestamp('19900315')])
        nat_series_dtype_timestamp = Series([NaT, NaT], dtype='datetime64[ns]')
        single_nat_dtype_datetime = Series([NaT], dtype='datetime64[ns]')

        # subtraction
        tm.assert_series_equal(-NaT + datetime_series,
                               nat_series_dtype_timestamp)
        with pytest.raises(TypeError):
            -single_nat_dtype_datetime + datetime_series

        tm.assert_series_equal(-NaT + nat_series_dtype_timestamp,
                               nat_series_dtype_timestamp)
        with pytest.raises(TypeError):
            -single_nat_dtype_datetime + nat_series_dtype_timestamp

        # addition
        tm.assert_series_equal(nat_series_dtype_timestamp + NaT,
                               nat_series_dtype_timestamp)
        tm.assert_series_equal(NaT + nat_series_dtype_timestamp,
                               nat_series_dtype_timestamp)

        tm.assert_series_equal(nat_series_dtype_timestamp + NaT,
                               nat_series_dtype_timestamp)
        tm.assert_series_equal(NaT + nat_series_dtype_timestamp,
                               nat_series_dtype_timestamp)

    # -------------------------------------------------------------
    # Invalid Operations
    # TODO: this block also needs to be de-duplicated and parametrized

    @pytest.mark.parametrize('dt64_series', [
        Series([Timestamp('19900315'), Timestamp('19900315')]),
        Series([pd.NaT, Timestamp('19900315')]),
        Series([pd.NaT, pd.NaT], dtype='datetime64[ns]')])
    @pytest.mark.parametrize('one', [1, 1.0, np.array(1)])
    def test_dt64_mul_div_numeric_invalid(self, one, dt64_series):
        # multiplication
        with pytest.raises(TypeError):
            dt64_series * one
        with pytest.raises(TypeError):
            one * dt64_series

        # division
        with pytest.raises(TypeError):
            dt64_series / one
        with pytest.raises(TypeError):
            one / dt64_series

    @pytest.mark.parametrize('op', ['__add__', '__radd__',
                                    '__sub__', '__rsub__'])
    @pytest.mark.parametrize('tz', [None, 'Asia/Tokyo'])
    def test_dt64_series_add_intlike(self, tz, op):
        # GH#19123
        dti = pd.DatetimeIndex(['2016-01-02', '2016-02-03', 'NaT'], tz=tz)
        ser = Series(dti)

        other = Series([20, 30, 40], dtype='uint8')

        pytest.raises(TypeError, getattr(ser, op), 1)

        pytest.raises(TypeError, getattr(ser, op), other)

        pytest.raises(TypeError, getattr(ser, op), other.values)

        pytest.raises(TypeError, getattr(ser, op), pd.Index(other))

    # -------------------------------------------------------------
    # Timezone-Centric Tests

    def test_operators_datetimelike_with_timezones(self):
        tz = 'US/Eastern'
        dt1 = Series(date_range('2000-01-01 09:00:00', periods=5,
                                tz=tz), name='foo')
        dt2 = dt1.copy()
        dt2.iloc[2] = np.nan

        td1 = Series(pd.timedelta_range('1 days 1 min', periods=5, freq='H'))
        td2 = td1.copy()
        td2.iloc[1] = np.nan

        result = dt1 + td1[0]
        exp = (dt1.dt.tz_localize(None) + td1[0]).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        result = dt2 + td2[0]
        exp = (dt2.dt.tz_localize(None) + td2[0]).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        # odd numpy behavior with scalar timedeltas
        result = td1[0] + dt1
        exp = (dt1.dt.tz_localize(None) + td1[0]).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        result = td2[0] + dt2
        exp = (dt2.dt.tz_localize(None) + td2[0]).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        result = dt1 - td1[0]
        exp = (dt1.dt.tz_localize(None) - td1[0]).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)
        with pytest.raises(TypeError):
            td1[0] - dt1

        result = dt2 - td2[0]
        exp = (dt2.dt.tz_localize(None) - td2[0]).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)
        with pytest.raises(TypeError):
            td2[0] - dt2

        result = dt1 + td1
        exp = (dt1.dt.tz_localize(None) + td1).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        result = dt2 + td2
        exp = (dt2.dt.tz_localize(None) + td2).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        result = dt1 - td1
        exp = (dt1.dt.tz_localize(None) - td1).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        result = dt2 - td2
        exp = (dt2.dt.tz_localize(None) - td2).dt.tz_localize(tz)
        tm.assert_series_equal(result, exp)

        with pytest.raises(TypeError):
            td1 - dt1
        with pytest.raises(TypeError):
            td2 - dt2


class TestDatetimeIndexArithmetic(object):

    # -------------------------------------------------------------
    # Invalid Operations

    @pytest.mark.parametrize('other', [3.14, np.array([2.0, 3.0])])
    @pytest.mark.parametrize('op', [operator.add, ops.radd,
                                    operator.sub, ops.rsub])
    def test_dti_add_sub_float(self, op, other):
        dti = DatetimeIndex(['2011-01-01', '2011-01-02'], freq='D')
        with pytest.raises(TypeError):
            op(dti, other)

    def test_dti_add_timestamp_raises(self, box):
        # GH#22163 ensure DataFrame doesn't cast Timestamp to i8
        idx = DatetimeIndex(['2011-01-01', '2011-01-02'])
        idx = tm.box_expected(idx, box)
        msg = "cannot add"
        with tm.assert_raises_regex(TypeError, msg):
            idx + Timestamp('2011-01-01')

    def test_dti_radd_timestamp_raises(self):
        idx = DatetimeIndex(['2011-01-01', '2011-01-02'])
        msg = "cannot add DatetimeIndex and Timestamp"
        with tm.assert_raises_regex(TypeError, msg):
            Timestamp('2011-01-01') + idx

    # -------------------------------------------------------------
    # Binary operations DatetimeIndex and int

    def test_dti_add_int(self, tz_naive_fixture, one):
        # Variants of `one` for #19012
        tz = tz_naive_fixture
        rng = pd.date_range('2000-01-01 09:00', freq='H',
                            periods=10, tz=tz)
        result = rng + one
        expected = pd.date_range('2000-01-01 10:00', freq='H',
                                 periods=10, tz=tz)
        tm.assert_index_equal(result, expected)

    def test_dti_iadd_int(self, tz_naive_fixture, one):
        tz = tz_naive_fixture
        rng = pd.date_range('2000-01-01 09:00', freq='H',
                            periods=10, tz=tz)
        expected = pd.date_range('2000-01-01 10:00', freq='H',
                                 periods=10, tz=tz)
        rng += one
        tm.assert_index_equal(rng, expected)

    def test_dti_sub_int(self, tz_naive_fixture, one):
        tz = tz_naive_fixture
        rng = pd.date_range('2000-01-01 09:00', freq='H',
                            periods=10, tz=tz)
        result = rng - one
        expected = pd.date_range('2000-01-01 08:00', freq='H',
                                 periods=10, tz=tz)
        tm.assert_index_equal(result, expected)

    def test_dti_isub_int(self, tz_naive_fixture, one):
        tz = tz_naive_fixture
        rng = pd.date_range('2000-01-01 09:00', freq='H',
                            periods=10, tz=tz)
        expected = pd.date_range('2000-01-01 08:00', freq='H',
                                 periods=10, tz=tz)
        rng -= one
        tm.assert_index_equal(rng, expected)

    # -------------------------------------------------------------
    # __add__/__sub__ with integer arrays

    @pytest.mark.parametrize('freq', ['H', 'D'])
    @pytest.mark.parametrize('box', [np.array, pd.Index])
    def test_dti_add_intarray_tick(self, box, freq):
        # GH#19959
        dti = pd.date_range('2016-01-01', periods=2, freq=freq)
        other = box([4, -1])
        expected = DatetimeIndex([dti[n] + other[n] for n in range(len(dti))])
        result = dti + other
        tm.assert_index_equal(result, expected)
        result = other + dti
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('freq', ['W', 'M', 'MS', 'Q'])
    @pytest.mark.parametrize('box', [np.array, pd.Index])
    def test_dti_add_intarray_non_tick(self, box, freq):
        # GH#19959
        dti = pd.date_range('2016-01-01', periods=2, freq=freq)
        other = box([4, -1])
        expected = DatetimeIndex([dti[n] + other[n] for n in range(len(dti))])
        with tm.assert_produces_warning(PerformanceWarning):
            result = dti + other
        tm.assert_index_equal(result, expected)
        with tm.assert_produces_warning(PerformanceWarning):
            result = other + dti
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('box', [np.array, pd.Index])
    def test_dti_add_intarray_no_freq(self, box):
        # GH#19959
        dti = pd.DatetimeIndex(['2016-01-01', 'NaT', '2017-04-05 06:07:08'])
        other = box([9, 4, -1])
        with pytest.raises(NullFrequencyError):
            dti + other
        with pytest.raises(NullFrequencyError):
            other + dti
        with pytest.raises(NullFrequencyError):
            dti - other
        with pytest.raises(TypeError):
            other - dti

    # -------------------------------------------------------------
    # Binary operations DatetimeIndex and timedelta-like

    def test_dti_add_timedeltalike(self, tz_naive_fixture, two_hours, box):
        # GH#22005, GH#22163 check DataFrame doesn't raise TypeError
        tz = tz_naive_fixture
        rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
        rng = tm.box_expected(rng, box)

        result = rng + two_hours
        expected = pd.date_range('2000-01-01 02:00',
                                 '2000-02-01 02:00', tz=tz)
        expected = tm.box_expected(expected, box)
        tm.assert_equal(result, expected)

    def test_dti_iadd_timedeltalike(self, tz_naive_fixture, two_hours):
        tz = tz_naive_fixture
        rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
        expected = pd.date_range('2000-01-01 02:00',
                                 '2000-02-01 02:00', tz=tz)
        rng += two_hours
        tm.assert_index_equal(rng, expected)

    def test_dti_sub_timedeltalike(self, tz_naive_fixture, two_hours):
        tz = tz_naive_fixture
        rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
        expected = pd.date_range('1999-12-31 22:00',
                                 '2000-01-31 22:00', tz=tz)
        result = rng - two_hours
        tm.assert_index_equal(result, expected)

    def test_dti_isub_timedeltalike(self, tz_naive_fixture, two_hours):
        tz = tz_naive_fixture
        rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
        expected = pd.date_range('1999-12-31 22:00',
                                 '2000-01-31 22:00', tz=tz)
        rng -= two_hours
        tm.assert_index_equal(rng, expected)

    # -------------------------------------------------------------
    # Binary operations DatetimeIndex and TimedeltaIndex/array
    def test_dti_add_tdi(self, tz_naive_fixture):
        # GH#17558
        tz = tz_naive_fixture
        dti = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        tdi = pd.timedelta_range('0 days', periods=10)
        expected = pd.date_range('2017-01-01', periods=10, tz=tz)

        # add with TimdeltaIndex
        result = dti + tdi
        tm.assert_index_equal(result, expected)

        result = tdi + dti
        tm.assert_index_equal(result, expected)

        # add with timedelta64 array
        result = dti + tdi.values
        tm.assert_index_equal(result, expected)

        result = tdi.values + dti
        tm.assert_index_equal(result, expected)

    def test_dti_iadd_tdi(self, tz_naive_fixture):
        # GH#17558
        tz = tz_naive_fixture
        dti = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        tdi = pd.timedelta_range('0 days', periods=10)
        expected = pd.date_range('2017-01-01', periods=10, tz=tz)

        # iadd with TimdeltaIndex
        result = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        result += tdi
        tm.assert_index_equal(result, expected)

        result = pd.timedelta_range('0 days', periods=10)
        result += dti
        tm.assert_index_equal(result, expected)

        # iadd with timedelta64 array
        result = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        result += tdi.values
        tm.assert_index_equal(result, expected)

        result = pd.timedelta_range('0 days', periods=10)
        result += dti
        tm.assert_index_equal(result, expected)

    def test_dti_sub_tdi(self, tz_naive_fixture):
        # GH#17558
        tz = tz_naive_fixture
        dti = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        tdi = pd.timedelta_range('0 days', periods=10)
        expected = pd.date_range('2017-01-01', periods=10, tz=tz, freq='-1D')

        # sub with TimedeltaIndex
        result = dti - tdi
        tm.assert_index_equal(result, expected)

        msg = 'cannot subtract .*TimedeltaIndex'
        with tm.assert_raises_regex(TypeError, msg):
            tdi - dti

        # sub with timedelta64 array
        result = dti - tdi.values
        tm.assert_index_equal(result, expected)

        msg = 'cannot subtract DatetimeIndex from'
        with tm.assert_raises_regex(TypeError, msg):
            tdi.values - dti

    def test_dti_isub_tdi(self, tz_naive_fixture):
        # GH#17558
        tz = tz_naive_fixture
        dti = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        tdi = pd.timedelta_range('0 days', periods=10)
        expected = pd.date_range('2017-01-01', periods=10, tz=tz, freq='-1D')

        # isub with TimedeltaIndex
        result = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        result -= tdi
        tm.assert_index_equal(result, expected)

        msg = 'cannot subtract .*TimedeltaIndex'
        with tm.assert_raises_regex(TypeError, msg):
            tdi -= dti

        # isub with timedelta64 array
        result = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        result -= tdi.values
        tm.assert_index_equal(result, expected)

        msg = '|'.join(['cannot perform __neg__ with this index type:',
                        'ufunc subtract cannot use operands with types',
                        'cannot subtract DatetimeIndex from'])
        with tm.assert_raises_regex(TypeError, msg):
            tdi.values -= dti

    # -------------------------------------------------------------
    # Binary Operations DatetimeIndex and datetime-like
    # TODO: A couple other tests belong in this section.  Move them in
    # A PR where there isn't already a giant diff.

    @pytest.mark.parametrize('addend', [
        datetime(2011, 1, 1),
        DatetimeIndex(['2011-01-01', '2011-01-02']),
        DatetimeIndex(['2011-01-01', '2011-01-02']).tz_localize('US/Eastern'),
        np.datetime64('2011-01-01'),
        Timestamp('2011-01-01')
    ], ids=lambda x: type(x).__name__)
    @pytest.mark.parametrize('tz', [None, 'US/Eastern'])
    def test_add_datetimelike_and_dti(self, addend, tz):
        # GH#9631
        dti = DatetimeIndex(['2011-01-01', '2011-01-02']).tz_localize(tz)
        msg = 'cannot add DatetimeIndex and {0}'.format(type(addend).__name__)
        with tm.assert_raises_regex(TypeError, msg):
            dti + addend
        with tm.assert_raises_regex(TypeError, msg):
            addend + dti

    # -------------------------------------------------------------
    # __add__/__sub__ with ndarray[datetime64] and ndarray[timedelta64]

    def test_dti_add_dt64_array_raises(self, tz_naive_fixture):
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        dtarr = dti.values

        with pytest.raises(TypeError):
            dti + dtarr
        with pytest.raises(TypeError):
            dtarr + dti

    def test_dti_sub_dt64_array_naive(self):
        dti = pd.date_range('2016-01-01', periods=3, tz=None)
        dtarr = dti.values

        expected = dti - dti
        result = dti - dtarr
        tm.assert_index_equal(result, expected)
        result = dtarr - dti
        tm.assert_index_equal(result, expected)

    def test_dti_sub_dt64_array_aware_raises(self, tz_naive_fixture):
        tz = tz_naive_fixture
        if tz is None:
            return
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        dtarr = dti.values

        with pytest.raises(TypeError):
            dti - dtarr
        with pytest.raises(TypeError):
            dtarr - dti

    def test_dti_add_td64_array(self, tz_naive_fixture):
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        tdi = pd.TimedeltaIndex(['-1 Day', '-1 Day', '-1 Day'])
        tdarr = tdi.values

        expected = dti + tdi
        result = dti + tdarr
        tm.assert_index_equal(result, expected)
        result = tdarr + dti
        tm.assert_index_equal(result, expected)

    def test_dti_sub_td64_array(self, tz_naive_fixture):
        tz = tz_naive_fixture
        dti = pd.date_range('2016-01-01', periods=3, tz=tz)
        tdi = pd.TimedeltaIndex(['-1 Day', '-1 Day', '-1 Day'])
        tdarr = tdi.values

        expected = dti - tdi
        result = dti - tdarr
        tm.assert_index_equal(result, expected)

        with pytest.raises(TypeError):
            tdarr - dti

    # -------------------------------------------------------------

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

        with pytest.raises(TypeError):
            dti_tz - dti

        with pytest.raises(TypeError):
            dti - dti_tz

        with pytest.raises(TypeError):
            dti_tz - dti_tz2

        # isub
        dti -= dti
        tm.assert_index_equal(dti, expected)

        # different length raises ValueError
        dti1 = date_range('20130101', periods=3)
        dti2 = date_range('20130101', periods=4)
        with pytest.raises(ValueError):
            dti1 - dti2

        # NaN propagation
        dti1 = DatetimeIndex(['2012-01-01', np.nan, '2012-01-03'])
        dti2 = DatetimeIndex(['2012-01-02', '2012-01-03', np.nan])
        expected = TimedeltaIndex(['1 days', np.nan, np.nan])
        result = dti2 - dti1
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('freq', [None, 'D'])
    def test_sub_period(self, freq, box):
        # GH#13078
        # not supported, check TypeError
        p = pd.Period('2011-01-01', freq='D')

        idx = pd.DatetimeIndex(['2011-01-01', '2011-01-02'], freq=freq)
        idx = tm.box_expected(idx, box)

        with pytest.raises(TypeError):
            idx - p

        with pytest.raises(TypeError):
            p - idx

    @pytest.mark.parametrize('op', [operator.add, ops.radd,
                                    operator.sub, ops.rsub])
    @pytest.mark.parametrize('pi_freq', ['D', 'W', 'Q', 'H'])
    @pytest.mark.parametrize('dti_freq', [None, 'D'])
    def test_dti_sub_pi(self, dti_freq, pi_freq, op, box):
        # GH#20049 subtracting PeriodIndex should raise TypeError
        dti = pd.DatetimeIndex(['2011-01-01', '2011-01-02'], freq=dti_freq)
        pi = dti.to_period(pi_freq)

        dti = tm.box_expected(dti, box, transpose=True)
        # TODO: Also box pi?
        with pytest.raises(TypeError):
            op(dti, pi)

    # -------------------------------------------------------------------
    # TODO: Most of this block is moved from series or frame tests, needs
    # cleanup, box-parametrization, and de-duplication

    @pytest.mark.parametrize('op', [operator.add, operator.sub])
    def test_timedelta64_equal_timedelta_supported_ops(self, op):
        ser = Series([Timestamp('20130301'),
                      Timestamp('20130228 23:00:00'),
                      Timestamp('20130228 22:00:00'),
                      Timestamp('20130228 21:00:00')])

        intervals = ['D', 'h', 'm', 's', 'us']

        # TODO: unused
        # npy16_mappings = {'D': 24 * 60 * 60 * 1000000,
        #                   'h': 60 * 60 * 1000000,
        #                   'm': 60 * 1000000,
        #                   's': 1000000,
        #                   'us': 1}

        def timedelta64(*args):
            return sum(starmap(np.timedelta64, zip(args, intervals)))

        for d, h, m, s, us in product(*([range(2)] * 5)):
            nptd = timedelta64(d, h, m, s, us)
            pytd = timedelta(days=d, hours=h, minutes=m, seconds=s,
                             microseconds=us)
            lhs = op(ser, nptd)
            rhs = op(ser, pytd)

            tm.assert_series_equal(lhs, rhs)

    def test_ops_nat_mixed_datetime64_timedelta64(self):
        # GH#11349
        timedelta_series = Series([NaT, Timedelta('1s')])
        datetime_series = Series([NaT, Timestamp('19900315')])
        nat_series_dtype_timedelta = Series([NaT, NaT],
                                            dtype='timedelta64[ns]')
        nat_series_dtype_timestamp = Series([NaT, NaT], dtype='datetime64[ns]')
        single_nat_dtype_datetime = Series([NaT], dtype='datetime64[ns]')
        single_nat_dtype_timedelta = Series([NaT], dtype='timedelta64[ns]')

        # subtraction
        tm.assert_series_equal(datetime_series - single_nat_dtype_datetime,
                               nat_series_dtype_timedelta)

        tm.assert_series_equal(datetime_series - single_nat_dtype_timedelta,
                               nat_series_dtype_timestamp)
        tm.assert_series_equal(-single_nat_dtype_timedelta + datetime_series,
                               nat_series_dtype_timestamp)

        # without a Series wrapping the NaT, it is ambiguous
        # whether it is a datetime64 or timedelta64
        # defaults to interpreting it as timedelta64
        tm.assert_series_equal(nat_series_dtype_timestamp -
                               single_nat_dtype_datetime,
                               nat_series_dtype_timedelta)

        tm.assert_series_equal(nat_series_dtype_timestamp -
                               single_nat_dtype_timedelta,
                               nat_series_dtype_timestamp)
        tm.assert_series_equal(-single_nat_dtype_timedelta +
                               nat_series_dtype_timestamp,
                               nat_series_dtype_timestamp)

        with pytest.raises(TypeError):
            timedelta_series - single_nat_dtype_datetime

        # addition
        tm.assert_series_equal(nat_series_dtype_timestamp +
                               single_nat_dtype_timedelta,
                               nat_series_dtype_timestamp)
        tm.assert_series_equal(single_nat_dtype_timedelta +
                               nat_series_dtype_timestamp,
                               nat_series_dtype_timestamp)

        tm.assert_series_equal(nat_series_dtype_timestamp +
                               single_nat_dtype_timedelta,
                               nat_series_dtype_timestamp)
        tm.assert_series_equal(single_nat_dtype_timedelta +
                               nat_series_dtype_timestamp,
                               nat_series_dtype_timestamp)

        tm.assert_series_equal(nat_series_dtype_timedelta +
                               single_nat_dtype_datetime,
                               nat_series_dtype_timestamp)
        tm.assert_series_equal(single_nat_dtype_datetime +
                               nat_series_dtype_timedelta,
                               nat_series_dtype_timestamp)

    def test_ufunc_coercions(self):
        idx = date_range('2011-01-01', periods=3, freq='2D', name='x')

        delta = np.timedelta64(1, 'D')
        for result in [idx + delta, np.add(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = date_range('2011-01-02', periods=3, freq='2D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == '2D'

        for result in [idx - delta, np.subtract(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = date_range('2010-12-31', periods=3, freq='2D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == '2D'

        delta = np.array([np.timedelta64(1, 'D'), np.timedelta64(2, 'D'),
                          np.timedelta64(3, 'D')])
        for result in [idx + delta, np.add(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = DatetimeIndex(['2011-01-02', '2011-01-05', '2011-01-08'],
                                freq='3D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == '3D'

        for result in [idx - delta, np.subtract(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = DatetimeIndex(['2010-12-31', '2011-01-01', '2011-01-02'],
                                freq='D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == 'D'

    def test_datetimeindex_sub_timestamp_overflow(self):
        dtimax = pd.to_datetime(['now', pd.Timestamp.max])
        dtimin = pd.to_datetime(['now', pd.Timestamp.min])

        tsneg = Timestamp('1950-01-01')
        ts_neg_variants = [tsneg,
                           tsneg.to_pydatetime(),
                           tsneg.to_datetime64().astype('datetime64[ns]'),
                           tsneg.to_datetime64().astype('datetime64[D]')]

        tspos = Timestamp('1980-01-01')
        ts_pos_variants = [tspos,
                           tspos.to_pydatetime(),
                           tspos.to_datetime64().astype('datetime64[ns]'),
                           tspos.to_datetime64().astype('datetime64[D]')]

        for variant in ts_neg_variants:
            with pytest.raises(OverflowError):
                dtimax - variant

        expected = pd.Timestamp.max.value - tspos.value
        for variant in ts_pos_variants:
            res = dtimax - variant
            assert res[1].value == expected

        expected = pd.Timestamp.min.value - tsneg.value
        for variant in ts_neg_variants:
            res = dtimin - variant
            assert res[1].value == expected

        for variant in ts_pos_variants:
            with pytest.raises(OverflowError):
                dtimin - variant

    def test_datetimeindex_sub_datetimeindex_overflow(self):
        # GH#22492, GH#22508
        dtimax = pd.to_datetime(['now', pd.Timestamp.max])
        dtimin = pd.to_datetime(['now', pd.Timestamp.min])

        ts_neg = pd.to_datetime(['1950-01-01', '1950-01-01'])
        ts_pos = pd.to_datetime(['1980-01-01', '1980-01-01'])

        # General tests
        expected = pd.Timestamp.max.value - ts_pos[1].value
        result = dtimax - ts_pos
        assert result[1].value == expected

        expected = pd.Timestamp.min.value - ts_neg[1].value
        result = dtimin - ts_neg
        assert result[1].value == expected

        with pytest.raises(OverflowError):
            dtimax - ts_neg

        with pytest.raises(OverflowError):
            dtimin - ts_pos

        # Edge cases
        tmin = pd.to_datetime([pd.Timestamp.min])
        t1 = tmin + pd.Timedelta.max + pd.Timedelta('1us')
        with pytest.raises(OverflowError):
            t1 - tmin

        tmax = pd.to_datetime([pd.Timestamp.max])
        t2 = tmax + pd.Timedelta.min - pd.Timedelta('1us')
        with pytest.raises(OverflowError):
            tmax - t2

    @pytest.mark.parametrize('names', [('foo', None, None),
                                       ('baz', 'bar', None),
                                       ('bar', 'bar', 'bar')])
    @pytest.mark.parametrize('tz', [None, 'America/Chicago'])
    def test_dti_add_series(self, tz, names):
        # GH#13905
        index = DatetimeIndex(['2016-06-28 05:30', '2016-06-28 05:31'],
                              tz=tz, name=names[0])
        ser = Series([Timedelta(seconds=5)] * 2,
                     index=index, name=names[1])
        expected = Series(index + Timedelta(seconds=5),
                          index=index, name=names[2])

        # passing name arg isn't enough when names[2] is None
        expected.name = names[2]
        assert expected.dtype == index.dtype
        result = ser + index
        tm.assert_series_equal(result, expected)
        result2 = index + ser
        tm.assert_series_equal(result2, expected)

        expected = index + Timedelta(seconds=5)
        result3 = ser.values + index
        tm.assert_index_equal(result3, expected)
        result4 = index + ser.values
        tm.assert_index_equal(result4, expected)

    def test_dti_add_offset_array(self, tz_naive_fixture):
        # GH#18849
        tz = tz_naive_fixture
        dti = pd.date_range('2017-01-01', periods=2, tz=tz)
        other = np.array([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti + other
        expected = DatetimeIndex([dti[n] + other[n] for n in range(len(dti))],
                                 name=dti.name, freq='infer')
        tm.assert_index_equal(res, expected)

        with tm.assert_produces_warning(PerformanceWarning):
            res2 = other + dti
        tm.assert_index_equal(res2, expected)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('foo', 'foo', 'foo')])
    def test_dti_add_offset_index(self, tz_naive_fixture, names):
        # GH#18849, GH#19744
        tz = tz_naive_fixture
        dti = pd.date_range('2017-01-01', periods=2, tz=tz, name=names[0])
        other = pd.Index([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)],
                         name=names[1])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti + other
        expected = DatetimeIndex([dti[n] + other[n] for n in range(len(dti))],
                                 name=names[2], freq='infer')
        tm.assert_index_equal(res, expected)

        with tm.assert_produces_warning(PerformanceWarning):
            res2 = other + dti
        tm.assert_index_equal(res2, expected)

    def test_dti_sub_offset_array(self, tz_naive_fixture):
        # GH#18824
        tz = tz_naive_fixture
        dti = pd.date_range('2017-01-01', periods=2, tz=tz)
        other = np.array([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti - other
        expected = DatetimeIndex([dti[n] - other[n] for n in range(len(dti))],
                                 name=dti.name, freq='infer')
        tm.assert_index_equal(res, expected)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('foo', 'foo', 'foo')])
    def test_dti_sub_offset_index(self, tz_naive_fixture, names):
        # GH#18824, GH#19744
        tz = tz_naive_fixture
        dti = pd.date_range('2017-01-01', periods=2, tz=tz, name=names[0])
        other = pd.Index([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)],
                         name=names[1])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti - other
        expected = DatetimeIndex([dti[n] - other[n] for n in range(len(dti))],
                                 name=names[2], freq='infer')
        tm.assert_index_equal(res, expected)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('foo', 'foo', 'foo')])
    def test_dti_with_offset_series(self, tz_naive_fixture, names):
        # GH#18849
        tz = tz_naive_fixture
        dti = pd.date_range('2017-01-01', periods=2, tz=tz, name=names[0])
        other = Series([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)],
                       name=names[1])

        expected_add = Series([dti[n] + other[n] for n in range(len(dti))],
                              name=names[2])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti + other
        tm.assert_series_equal(res, expected_add)

        with tm.assert_produces_warning(PerformanceWarning):
            res2 = other + dti
        tm.assert_series_equal(res2, expected_add)

        expected_sub = Series([dti[n] - other[n] for n in range(len(dti))],
                              name=names[2])

        with tm.assert_produces_warning(PerformanceWarning):
            res3 = dti - other
        tm.assert_series_equal(res3, expected_sub)

    def test_dti_add_offset_tzaware(self, tz_aware_fixture, box):
        # GH#21610, GH#22163 ensure DataFrame doesn't return object-dtype
        timezone = tz_aware_fixture
        if timezone == 'US/Pacific':
            dates = date_range('2012-11-01', periods=3, tz=timezone)
            offset = dates + pd.offsets.Hour(5)
            assert dates[0] + pd.offsets.Hour(5) == offset[0]

        dates = date_range('2010-11-01 00:00',
                           periods=3, tz=timezone, freq='H')
        expected = DatetimeIndex(['2010-11-01 05:00', '2010-11-01 06:00',
                                  '2010-11-01 07:00'], freq='H', tz=timezone)

        dates = tm.box_expected(dates, box)
        expected = tm.box_expected(expected, box)

        # TODO: parametrize over the scalar being added?  radd?  sub?
        offset = dates + pd.offsets.Hour(5)
        tm.assert_equal(offset, expected)
        offset = dates + np.timedelta64(5, 'h')
        tm.assert_equal(offset, expected)
        offset = dates + timedelta(hours=5)
        tm.assert_equal(offset, expected)


@pytest.mark.parametrize('klass', [Series, DatetimeIndex])
def test_dt64_with_offset_array(klass):
    # GH#10699
    # array of offsets
    box = Series if klass is Series else pd.Index
    dti = DatetimeIndex([Timestamp('2000-1-1'), Timestamp('2000-2-1')])

    s = klass(dti)

    with tm.assert_produces_warning(PerformanceWarning):
        result = s + box([pd.offsets.DateOffset(years=1),
                          pd.offsets.MonthEnd()])
        exp = klass([Timestamp('2001-1-1'), Timestamp('2000-2-29')])
        tm.assert_equal(result, exp)

        # same offset
        result = s + box([pd.offsets.DateOffset(years=1),
                          pd.offsets.DateOffset(years=1)])
        exp = klass([Timestamp('2001-1-1'), Timestamp('2001-2-1')])
        tm.assert_equal(result, exp)


@pytest.mark.parametrize('klass', [Series, DatetimeIndex])
def test_dt64_with_DateOffsets_relativedelta(klass):
    # GH#10699
    vec = klass([Timestamp('2000-01-05 00:15:00'),
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
        tm.assert_equal(klass([x + op for x in vec]), vec + op)
        tm.assert_equal(klass([x - op for x in vec]), vec - op)
        op = pd.DateOffset(**dict(relative_kwargs[:i + 1]))
        tm.assert_equal(klass([x + op for x in vec]), vec + op)
        tm.assert_equal(klass([x - op for x in vec]), vec - op)


@pytest.mark.parametrize('cls_and_kwargs', [
    'YearBegin', ('YearBegin', {'month': 5}),
    'YearEnd', ('YearEnd', {'month': 5}),
    'MonthBegin', 'MonthEnd',
    'SemiMonthEnd', 'SemiMonthBegin',
    'Week', ('Week', {'weekday': 3}),
    'Week', ('Week', {'weekday': 6}),
    'BusinessDay', 'BDay', 'QuarterEnd', 'QuarterBegin',
    'CustomBusinessDay', 'CDay', 'CBMonthEnd',
    'CBMonthBegin', 'BMonthBegin', 'BMonthEnd',
    'BusinessHour', 'BYearBegin', 'BYearEnd',
    'BQuarterBegin', ('LastWeekOfMonth', {'weekday': 2}),
    ('FY5253Quarter', {'qtr_with_extra_week': 1,
                       'startingMonth': 1,
                       'weekday': 2,
                       'variation': 'nearest'}),
    ('FY5253', {'weekday': 0, 'startingMonth': 2, 'variation': 'nearest'}),
    ('WeekOfMonth', {'weekday': 2, 'week': 2}),
    'Easter', ('DateOffset', {'day': 4}),
    ('DateOffset', {'month': 5})])
@pytest.mark.parametrize('normalize', [True, False])
@pytest.mark.parametrize('klass', [Series, DatetimeIndex])
def test_dt64_with_DateOffsets(klass, normalize, cls_and_kwargs):
    # GH#10699
    # assert these are equal on a piecewise basis
    vec = klass([Timestamp('2000-01-05 00:15:00'),
                 Timestamp('2000-01-31 00:23:00'),
                 Timestamp('2000-01-01'),
                 Timestamp('2000-03-31'),
                 Timestamp('2000-02-29'),
                 Timestamp('2000-12-31'),
                 Timestamp('2000-05-15'),
                 Timestamp('2001-06-15')])

    if isinstance(cls_and_kwargs, tuple):
        # If cls_name param is a tuple, then 2nd entry is kwargs for
        # the offset constructor
        cls_name, kwargs = cls_and_kwargs
    else:
        cls_name = cls_and_kwargs
        kwargs = {}

    offset_cls = getattr(pd.offsets, cls_name)

    with warnings.catch_warnings(record=True):
        # pandas.errors.PerformanceWarning: Non-vectorized DateOffset being
        # applied to Series or DatetimeIndex
        # we aren't testing that here, so ignore.
        warnings.simplefilter("ignore", PerformanceWarning)
        for n in [0, 5]:
            if (cls_name in ['WeekOfMonth', 'LastWeekOfMonth',
                             'FY5253Quarter', 'FY5253'] and n == 0):
                # passing n = 0 is invalid for these offset classes
                continue

            offset = offset_cls(n, normalize=normalize, **kwargs)
            tm.assert_equal(klass([x + offset for x in vec]), vec + offset)
            tm.assert_equal(klass([x - offset for x in vec]), vec - offset)
            tm.assert_equal(klass([offset + x for x in vec]), offset + vec)


@pytest.mark.parametrize('klass', [Series, DatetimeIndex])
def test_datetime64_with_DateOffset(klass):
    # GH#10699
    s = klass(date_range('2000-01-01', '2000-01-31'), name='a')
    result = s + pd.DateOffset(years=1)
    result2 = pd.DateOffset(years=1) + s
    exp = klass(date_range('2001-01-01', '2001-01-31'), name='a')
    tm.assert_equal(result, exp)
    tm.assert_equal(result2, exp)

    result = s - pd.DateOffset(years=1)
    exp = klass(date_range('1999-01-01', '1999-01-31'), name='a')
    tm.assert_equal(result, exp)

    s = klass([Timestamp('2000-01-15 00:15:00', tz='US/Central'),
               pd.Timestamp('2000-02-15', tz='US/Central')], name='a')
    result = s + pd.offsets.Day()
    result2 = pd.offsets.Day() + s
    exp = klass([Timestamp('2000-01-16 00:15:00', tz='US/Central'),
                 Timestamp('2000-02-16', tz='US/Central')], name='a')
    tm.assert_equal(result, exp)
    tm.assert_equal(result2, exp)

    s = klass([Timestamp('2000-01-15 00:15:00', tz='US/Central'),
               pd.Timestamp('2000-02-15', tz='US/Central')], name='a')
    result = s + pd.offsets.MonthEnd()
    result2 = pd.offsets.MonthEnd() + s
    exp = klass([Timestamp('2000-01-31 00:15:00', tz='US/Central'),
                 Timestamp('2000-02-29', tz='US/Central')], name='a')
    tm.assert_equal(result, exp)
    tm.assert_equal(result2, exp)


@pytest.mark.parametrize('years', [-1, 0, 1])
@pytest.mark.parametrize('months', [-2, 0, 2])
def test_shift_months(years, months):
    dti = DatetimeIndex([Timestamp('2000-01-05 00:15:00'),
                         Timestamp('2000-01-31 00:23:00'),
                         Timestamp('2000-01-01'),
                         Timestamp('2000-02-29'),
                         Timestamp('2000-12-31')])
    actual = DatetimeIndex(shift_months(dti.asi8, years * 12 + months))

    raw = [x + pd.offsets.DateOffset(years=years, months=months)
           for x in dti]
    expected = DatetimeIndex(raw)
    tm.assert_index_equal(actual, expected)
