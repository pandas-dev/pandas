# -*- coding: utf-8 -*-
# Arithmetc tests for DataFrame/Series/Index/Array classes that should
# behave identically.
# Specifically for numeric dtypes
from datetime import timedelta
from decimal import Decimal

import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas import Timedelta, Series, TimedeltaIndex


# ------------------------------------------------------------------
# Comparisons


# ------------------------------------------------------------------
# Numeric dtypes Arithmetic with Timedelta Scalar

class TestNumericArraylikeArithmeticWithTimedeltaScalar(object):

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="block.eval incorrect",
                                             strict=True))
    ])
    @pytest.mark.parametrize('index', [
        pd.Int64Index(range(1, 11)),
        pd.UInt64Index(range(1, 11)),
        pd.Float64Index(range(1, 11)),
        pd.RangeIndex(1, 11)],
        ids=lambda x: type(x).__name__)
    @pytest.mark.parametrize('scalar_td', [
        Timedelta(days=1),
        Timedelta(days=1).to_timedelta64(),
        Timedelta(days=1).to_pytimedelta()],
        ids=lambda x: type(x).__name__)
    def test_numeric_arr_mul_tdscalar(self, scalar_td, index, box):
        # GH#19333

        if (box is Series and
                type(scalar_td) is timedelta and index.dtype == 'f8'):
            raise pytest.xfail(reason="Cannot multiply timedelta by float")

        expected = pd.timedelta_range('1 days', '10 days')

        index = tm.box_expected(index, box)
        expected = tm.box_expected(expected, box)

        result = index * scalar_td
        tm.assert_equal(result, expected)

        commute = scalar_td * index
        tm.assert_equal(commute, expected)

    @pytest.mark.parametrize('index', [
        pd.Int64Index(range(1, 3)),
        pd.UInt64Index(range(1, 3)),
        pd.Float64Index(range(1, 3)),
        pd.RangeIndex(1, 3)],
        ids=lambda x: type(x).__name__)
    @pytest.mark.parametrize('scalar_td', [
        Timedelta(days=1),
        Timedelta(days=1).to_timedelta64(),
        Timedelta(days=1).to_pytimedelta()],
        ids=lambda x: type(x).__name__)
    def test_numeric_arr_rdiv_tdscalar(self, scalar_td, index, box):

        if box is Series and type(scalar_td) is timedelta:
            raise pytest.xfail(reason="TODO: Figure out why this case fails")
        if box is pd.DataFrame and isinstance(scalar_td, timedelta):
            raise pytest.xfail(reason="TODO: Figure out why this case fails")

        expected = TimedeltaIndex(['1 Day', '12 Hours'])

        index = tm.box_expected(index, box)
        expected = tm.box_expected(expected, box)

        result = scalar_td / index
        tm.assert_equal(result, expected)

        with pytest.raises(TypeError):
            index / scalar_td


# ------------------------------------------------------------------
# Arithmetic

class TestDivisionByZero(object):

    def test_div_zero(self, zero, idx):
        expected = pd.Index([np.nan, np.inf, np.inf, np.inf, np.inf],
                            dtype=np.float64)
        result = idx / zero
        tm.assert_index_equal(result, expected)
        ser_compat = Series(idx).astype('i8') / np.array(zero).astype('i8')
        tm.assert_series_equal(ser_compat, Series(result))

    def test_floordiv_zero(self, zero, idx):
        expected = pd.Index([np.nan, np.inf, np.inf, np.inf, np.inf],
                            dtype=np.float64)

        result = idx // zero
        tm.assert_index_equal(result, expected)
        ser_compat = Series(idx).astype('i8') // np.array(zero).astype('i8')
        tm.assert_series_equal(ser_compat, Series(result))

    def test_mod_zero(self, zero, idx):
        expected = pd.Index([np.nan, np.nan, np.nan, np.nan, np.nan],
                            dtype=np.float64)
        result = idx % zero
        tm.assert_index_equal(result, expected)
        ser_compat = Series(idx).astype('i8') % np.array(zero).astype('i8')
        tm.assert_series_equal(ser_compat, Series(result))

    def test_divmod_zero(self, zero, idx):

        exleft = pd.Index([np.nan, np.inf, np.inf, np.inf, np.inf],
                          dtype=np.float64)
        exright = pd.Index([np.nan, np.nan, np.nan, np.nan, np.nan],
                           dtype=np.float64)

        result = divmod(idx, zero)
        tm.assert_index_equal(result[0], exleft)
        tm.assert_index_equal(result[1], exright)

    # ------------------------------------------------------------------

    @pytest.mark.parametrize('dtype2', [
        np.int64, np.int32, np.int16, np.int8,
        np.float64, np.float32, np.float16,
        np.uint64, np.uint32, np.uint16, np.uint8])
    @pytest.mark.parametrize('dtype1', [np.int64, np.float64, np.uint64])
    def test_ser_div_ser(self, dtype1, dtype2):
        # no longer do integer div for any ops, but deal with the 0's
        first = Series([3, 4, 5, 8], name='first').astype(dtype1)
        second = Series([0, 0, 0, 3], name='second').astype(dtype2)

        with np.errstate(all='ignore'):
            expected = Series(first.values.astype(np.float64) / second.values,
                              dtype='float64', name=None)
        expected.iloc[0:3] = np.inf

        result = first / second
        tm.assert_series_equal(result, expected)
        assert not result.equals(second / first)

    def test_rdiv_zero_compat(self):
        # GH#8674
        zero_array = np.array([0] * 5)
        data = np.random.randn(5)
        expected = Series([0.] * 5)

        result = zero_array / Series(data)
        tm.assert_series_equal(result, expected)

        result = Series(zero_array) / data
        tm.assert_series_equal(result, expected)

        result = Series(zero_array) / Series(data)
        tm.assert_series_equal(result, expected)

    def test_div_zero_inf_signs(self):
        # GH#9144, inf signing
        ser = Series([-1, 0, 1], name='first')
        expected = Series([-np.inf, np.nan, np.inf], name='first')

        result = ser / 0
        tm.assert_series_equal(result, expected)

    def test_rdiv_zero(self):
        # GH#9144
        ser = Series([-1, 0, 1], name='first')
        expected = Series([0.0, np.nan, 0.0], name='first')

        result = 0 / ser
        tm.assert_series_equal(result, expected)

    def test_floordiv_div(self):
        # GH#9144
        ser = Series([-1, 0, 1], name='first')

        result = ser // 0
        expected = Series([-np.inf, np.nan, np.inf], name='first')
        tm.assert_series_equal(result, expected)

    def test_df_div_zero_df(self):
        # integer div, but deal with the 0's (GH#9144)
        df = pd.DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})
        result = df / df

        first = pd.Series([1.0, 1.0, 1.0, 1.0])
        second = pd.Series([np.nan, np.nan, np.nan, 1])
        expected = pd.DataFrame({'first': first, 'second': second})
        tm.assert_frame_equal(result, expected)

    def test_df_div_zero_array(self):
        # integer div, but deal with the 0's (GH#9144)
        df = pd.DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})

        first = pd.Series([1.0, 1.0, 1.0, 1.0])
        second = pd.Series([np.nan, np.nan, np.nan, 1])
        expected = pd.DataFrame({'first': first, 'second': second})

        with np.errstate(all='ignore'):
            arr = df.values.astype('float') / df.values
        result = pd.DataFrame(arr, index=df.index,
                              columns=df.columns)
        tm.assert_frame_equal(result, expected)

    def test_df_div_zero_int(self):
        # integer div, but deal with the 0's (GH#9144)
        df = pd.DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})

        result = df / 0
        expected = pd.DataFrame(np.inf, index=df.index, columns=df.columns)
        expected.iloc[0:3, 1] = np.nan
        tm.assert_frame_equal(result, expected)

        # numpy has a slightly different (wrong) treatment
        with np.errstate(all='ignore'):
            arr = df.values.astype('float64') / 0
        result2 = pd.DataFrame(arr, index=df.index,
                               columns=df.columns)
        tm.assert_frame_equal(result2, expected)

    def test_df_div_zero_series_does_not_commute(self):
        # integer div, but deal with the 0's (GH#9144)
        df = pd.DataFrame(np.random.randn(10, 5))
        ser = df[0]
        res = ser / df
        res2 = df / ser
        assert not res.fillna(0).equals(res2.fillna(0))

    # ------------------------------------------------------------------
    # Mod By Zero

    def test_df_mod_zero_df(self):
        # GH#3590, modulo as ints
        df = pd.DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})

        # this is technically wrong, as the integer portion is coerced to float
        # ###
        first = pd.Series([0, 0, 0, 0], dtype='float64')
        second = pd.Series([np.nan, np.nan, np.nan, 0])
        expected = pd.DataFrame({'first': first, 'second': second})
        result = df % df
        tm.assert_frame_equal(result, expected)

    def test_df_mod_zero_array(self):
        # GH#3590, modulo as ints
        df = pd.DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})

        # this is technically wrong, as the integer portion is coerced to float
        # ###
        first = pd.Series([0, 0, 0, 0], dtype='float64')
        second = pd.Series([np.nan, np.nan, np.nan, 0])
        expected = pd.DataFrame({'first': first, 'second': second})

        # numpy has a slightly different (wrong) treatment
        with np.errstate(all='ignore'):
            arr = df.values % df.values
        result2 = pd.DataFrame(arr, index=df.index,
                               columns=df.columns, dtype='float64')
        result2.iloc[0:3, 1] = np.nan
        tm.assert_frame_equal(result2, expected)

    def test_df_mod_zero_int(self):
        # GH#3590, modulo as ints
        df = pd.DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})

        result = df % 0
        expected = pd.DataFrame(np.nan, index=df.index, columns=df.columns)
        tm.assert_frame_equal(result, expected)

        # numpy has a slightly different (wrong) treatment
        with np.errstate(all='ignore'):
            arr = df.values.astype('float64') % 0
        result2 = pd.DataFrame(arr, index=df.index, columns=df.columns)
        tm.assert_frame_equal(result2, expected)

    def test_df_mod_zero_series_does_not_commute(self):
        # GH#3590, modulo as ints
        # not commutative with series
        df = pd.DataFrame(np.random.randn(10, 5))
        ser = df[0]
        res = ser % df
        res2 = df % ser
        assert not res.fillna(0).equals(res2.fillna(0))


class TestDivision(object):
    # __div__, __rdiv__, __floordiv__, __rfloordiv__
    # for non-timestamp/timedelta/period dtypes

    @pytest.mark.parametrize('box', [
        pytest.param(pd.Index,
                     marks=pytest.mark.xfail(reason="Index.__div__ always "
                                                    "raises",
                                             raises=TypeError, strict=True)),
        pd.Series,
        pd.DataFrame
    ], ids=lambda x: x.__name__)
    def test_divide_decimal(self, box):
        # resolves issue GH#9787
        ser = Series([Decimal(10)])
        expected = Series([Decimal(5)])

        ser = tm.box_expected(ser, box)
        expected = tm.box_expected(expected, box)

        result = ser / Decimal(2)

        tm.assert_equal(result, expected)

        result = ser // Decimal(2)
        tm.assert_equal(result, expected)

    def test_div_equiv_binop(self):
        # Test Series.div as well as Series.__div__
        # float/integer issue
        # GH#7785
        first = Series([1, 0], name='first')
        second = Series([-0.01, -0.02], name='second')
        expected = Series([-0.01, -np.inf])

        result = second.div(first)
        tm.assert_series_equal(result, expected, check_names=False)

        result = second / first
        tm.assert_series_equal(result, expected)


class TestAdditionSubtraction(object):
    # __add__, __sub__, __radd__, __rsub__, __iadd__, __isub__
    # for non-timestamp/timedelta/period dtypes
    pass


class TestObjectDtypeEquivalence(object):
    # Tests that arithmetic operations match operations executed elementwise

    @pytest.mark.parametrize('dtype', [None, object])
    def test_series_with_dtype_radd_nan(self, dtype):
        ser = pd.Series([1, 2, 3], dtype=dtype)
        expected = pd.Series([np.nan, np.nan, np.nan], dtype=dtype)

        result = np.nan + ser
        tm.assert_series_equal(result, expected)

        result = ser + np.nan
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('dtype', [None, object])
    def test_series_with_dtype_radd_int(self, dtype):
        ser = pd.Series([1, 2, 3], dtype=dtype)
        expected = pd.Series([2, 3, 4], dtype=dtype)

        result = 1 + ser
        tm.assert_series_equal(result, expected)

        result = ser + 1
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('dtype', [None, object])
    def test_df_with_dtype_radd_nan(self, dtype):
        df = pd.DataFrame([1, 2, 3], dtype=dtype)
        expected = pd.DataFrame([np.nan, np.nan, np.nan], dtype=dtype)

        result = np.nan + df
        tm.assert_frame_equal(result, expected)

        result = df + np.nan
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize('dtype', [None, object])
    def test_df_with_dtype_radd_int(self, dtype):
        df = pd.DataFrame([1, 2, 3], dtype=dtype)
        expected = pd.DataFrame([2, 3, 4], dtype=dtype)

        result = 1 + df
        tm.assert_frame_equal(result, expected)

        result = df + 1
        tm.assert_frame_equal(result, expected)
