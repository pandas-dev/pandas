# -*- coding: utf-8 -*-
# Arithmetc tests for DataFrame/Series/Index/Array classes that should
# behave identically.
# Specifically for numeric dtypes
from decimal import Decimal
import operator

import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm

from pandas.compat import PY3, Iterable
from pandas.core import ops
from pandas import Timedelta, Series, Index, TimedeltaIndex


# ------------------------------------------------------------------
# Comparisons

class TestNumericComparisons(object):
    def test_operator_series_comparison_zerorank(self):
        # GH#13006
        result = np.float64(0) > pd.Series([1, 2, 3])
        expected = 0.0 > pd.Series([1, 2, 3])
        tm.assert_series_equal(result, expected)
        result = pd.Series([1, 2, 3]) < np.float64(0)
        expected = pd.Series([1, 2, 3]) < 0.0
        tm.assert_series_equal(result, expected)
        result = np.array([0, 1, 2])[0] > pd.Series([0, 1, 2])
        expected = 0.0 > pd.Series([1, 2, 3])
        tm.assert_series_equal(result, expected)

    def test_df_numeric_cmp_dt64_raises(self):
        # GH#8932, GH#22163
        ts = pd.Timestamp.now()
        df = pd.DataFrame({'x': range(5)})
        with pytest.raises(TypeError):
            df > ts
        with pytest.raises(TypeError):
            df < ts
        with pytest.raises(TypeError):
            ts < df
        with pytest.raises(TypeError):
            ts > df

        assert not (df == ts).any().any()
        assert (df != ts).all().all()

    def test_compare_invalid(self):
        # GH#8058
        # ops testing
        a = pd.Series(np.random.randn(5), name=0)
        b = pd.Series(np.random.randn(5))
        b.name = pd.Timestamp('2000-01-01')
        tm.assert_series_equal(a / b, 1 / (b / a))


# ------------------------------------------------------------------
# Numeric dtypes Arithmetic with Timedelta Scalar

class TestNumericArraylikeArithmeticWithTimedeltaLike(object):

    # TODO: also check name retentention
    @pytest.mark.parametrize('box_cls', [np.array, pd.Index, pd.Series])
    @pytest.mark.parametrize('left', [
        pd.RangeIndex(10, 40, 10)] + [cls([10, 20, 30], dtype=dtype)
                                      for dtype in ['i1', 'i2', 'i4', 'i8',
                                                    'u1', 'u2', 'u4', 'u8',
                                                    'f2', 'f4', 'f8']
                                      for cls in [pd.Series, pd.Index]],
        ids=lambda x: type(x).__name__ + str(x.dtype))
    def test_mul_td64arr(self, left, box_cls):
        # GH#22390
        right = np.array([1, 2, 3], dtype='m8[s]')
        right = box_cls(right)

        expected = pd.TimedeltaIndex(['10s', '40s', '90s'])
        if isinstance(left, pd.Series) or box_cls is pd.Series:
            expected = pd.Series(expected)

        result = left * right
        tm.assert_equal(result, expected)

        result = right * left
        tm.assert_equal(result, expected)

    # TODO: also check name retentention
    @pytest.mark.parametrize('box_cls', [np.array, pd.Index, pd.Series])
    @pytest.mark.parametrize('left', [
        pd.RangeIndex(10, 40, 10)] + [cls([10, 20, 30], dtype=dtype)
                                      for dtype in ['i1', 'i2', 'i4', 'i8',
                                                    'u1', 'u2', 'u4', 'u8',
                                                    'f2', 'f4', 'f8']
                                      for cls in [pd.Series, pd.Index]],
        ids=lambda x: type(x).__name__ + str(x.dtype))
    def test_div_td64arr(self, left, box_cls):
        # GH#22390
        right = np.array([10, 40, 90], dtype='m8[s]')
        right = box_cls(right)

        expected = pd.TimedeltaIndex(['1s', '2s', '3s'])
        if isinstance(left, pd.Series) or box_cls is pd.Series:
            expected = pd.Series(expected)

        result = right / left
        tm.assert_equal(result, expected)

        result = right // left
        tm.assert_equal(result, expected)

        with pytest.raises(TypeError):
            left / right

        with pytest.raises(TypeError):
            left // right

    # TODO: de-duplicate with test_numeric_arr_mul_tdscalar
    def test_ops_series(self):
        # regression test for G#H8813
        td = Timedelta('1 day')
        other = pd.Series([1, 2])
        expected = pd.Series(pd.to_timedelta(['1 day', '2 days']))
        tm.assert_series_equal(expected, td * other)
        tm.assert_series_equal(expected, other * td)

    # TODO: also test non-nanosecond timedelta64 and Tick objects;
    #  see test_numeric_arr_rdiv_tdscalar for note on these failing
    @pytest.mark.parametrize('scalar_td', [
        Timedelta(days=1),
        Timedelta(days=1).to_timedelta64(),
        Timedelta(days=1).to_pytimedelta()],
        ids=lambda x: type(x).__name__)
    def test_numeric_arr_mul_tdscalar(self, scalar_td, numeric_idx, box):
        # GH#19333
        index = numeric_idx

        expected = pd.timedelta_range('0 days', '4 days')

        index = tm.box_expected(index, box)
        expected = tm.box_expected(expected, box)

        result = index * scalar_td
        tm.assert_equal(result, expected)

        commute = scalar_td * index
        tm.assert_equal(commute, expected)

    def test_numeric_arr_rdiv_tdscalar(self, three_days, numeric_idx, box):
        index = numeric_idx[1:3]

        broken = (isinstance(three_days, np.timedelta64) and
                  three_days.dtype != 'm8[ns]')
        broken = broken or isinstance(three_days, pd.offsets.Tick)
        if box is not pd.Index and broken:
            # np.timedelta64(3, 'D') / 2 == np.timedelta64(1, 'D')
            raise pytest.xfail("timedelta64 not converted to nanos; "
                               "Tick division not implemented")

        expected = TimedeltaIndex(['3 Days', '36 Hours'])

        index = tm.box_expected(index, box)
        expected = tm.box_expected(expected, box)

        result = three_days / index
        tm.assert_equal(result, expected)

        with pytest.raises(TypeError):
            index / three_days


# ------------------------------------------------------------------
# Arithmetic

class TestDivisionByZero(object):

    def test_div_zero(self, zero, numeric_idx):
        idx = numeric_idx

        expected = pd.Index([np.nan, np.inf, np.inf, np.inf, np.inf],
                            dtype=np.float64)
        result = idx / zero
        tm.assert_index_equal(result, expected)
        ser_compat = Series(idx).astype('i8') / np.array(zero).astype('i8')
        tm.assert_series_equal(ser_compat, Series(result))

    def test_floordiv_zero(self, zero, numeric_idx):
        idx = numeric_idx

        expected = pd.Index([np.nan, np.inf, np.inf, np.inf, np.inf],
                            dtype=np.float64)

        result = idx // zero
        tm.assert_index_equal(result, expected)
        ser_compat = Series(idx).astype('i8') // np.array(zero).astype('i8')
        tm.assert_series_equal(ser_compat, Series(result))

    def test_mod_zero(self, zero, numeric_idx):
        idx = numeric_idx

        expected = pd.Index([np.nan, np.nan, np.nan, np.nan, np.nan],
                            dtype=np.float64)
        result = idx % zero
        tm.assert_index_equal(result, expected)
        ser_compat = Series(idx).astype('i8') % np.array(zero).astype('i8')
        tm.assert_series_equal(ser_compat, Series(result))

    def test_divmod_zero(self, zero, numeric_idx):
        idx = numeric_idx

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


class TestMultiplicationDivision(object):
    # __mul__, __rmul__, __div__, __rdiv__, __floordiv__, __rfloordiv__
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

    def test_div_int(self, numeric_idx):
        # truediv under PY3
        idx = numeric_idx
        result = idx / 1
        expected = idx
        if PY3:
            expected = expected.astype('float64')
        tm.assert_index_equal(result, expected)

        result = idx / 2
        if PY3:
            expected = expected.astype('float64')
        expected = Index(idx.values / 2)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('op', [operator.mul, ops.rmul, operator.floordiv])
    def test_mul_int_identity(self, op, numeric_idx, box):
        idx = numeric_idx
        idx = tm.box_expected(idx, box)

        result = op(idx, 1)
        tm.assert_equal(result, idx)

    def test_mul_int_array(self, numeric_idx):
        idx = numeric_idx
        didx = idx * idx

        result = idx * np.array(5, dtype='int64')
        tm.assert_index_equal(result, idx * 5)

        arr_dtype = 'uint64' if isinstance(idx, pd.UInt64Index) else 'int64'
        result = idx * np.arange(5, dtype=arr_dtype)
        tm.assert_index_equal(result, didx)

    def test_mul_int_series(self, numeric_idx):
        idx = numeric_idx
        didx = idx * idx

        arr_dtype = 'uint64' if isinstance(idx, pd.UInt64Index) else 'int64'
        result = idx * Series(np.arange(5, dtype=arr_dtype))
        tm.assert_series_equal(result, Series(didx))

    def test_mul_float_series(self, numeric_idx):
        idx = numeric_idx
        rng5 = np.arange(5, dtype='float64')

        result = idx * Series(rng5 + 0.1)
        expected = Series(rng5 * (rng5 + 0.1))
        tm.assert_series_equal(result, expected)

    def test_mul_index(self, numeric_idx):
        # in general not true for RangeIndex
        idx = numeric_idx
        if not isinstance(idx, pd.RangeIndex):
            result = idx * idx
            tm.assert_index_equal(result, idx ** 2)

    def test_mul_datelike_raises(self, numeric_idx):
        idx = numeric_idx
        with pytest.raises(TypeError):
            idx * pd.date_range('20130101', periods=5)

    def test_mul_size_mismatch_raises(self, numeric_idx):
        idx = numeric_idx
        with pytest.raises(ValueError):
            idx * idx[0:3]
        with pytest.raises(ValueError):
            idx * np.array([1, 2])

    @pytest.mark.parametrize('op', [operator.pow, ops.rpow])
    def test_pow_float(self, op, numeric_idx, box):
        # test power calculations both ways, GH#14973
        idx = numeric_idx
        expected = pd.Float64Index(op(idx.values, 2.0))

        idx = tm.box_expected(idx, box)
        expected = tm.box_expected(expected, box)

        result = op(idx, 2.0)
        tm.assert_equal(result, expected)

    def test_modulo(self, numeric_idx, box):
        # GH#9244
        idx = numeric_idx
        expected = Index(idx.values % 2)

        idx = tm.box_expected(idx, box)
        expected = tm.box_expected(expected, box)

        result = idx % 2
        tm.assert_equal(result, expected)

    def test_divmod_scalar(self, numeric_idx):
        idx = numeric_idx

        result = divmod(idx, 2)
        with np.errstate(all='ignore'):
            div, mod = divmod(idx.values, 2)

        expected = Index(div), Index(mod)
        for r, e in zip(result, expected):
            tm.assert_index_equal(r, e)

    def test_divmod_ndarray(self, numeric_idx):
        idx = numeric_idx
        other = np.ones(idx.values.shape, dtype=idx.values.dtype) * 2

        result = divmod(idx, other)
        with np.errstate(all='ignore'):
            div, mod = divmod(idx.values, other)

        expected = Index(div), Index(mod)
        for r, e in zip(result, expected):
            tm.assert_index_equal(r, e)

    def test_divmod_series(self, numeric_idx):
        idx = numeric_idx
        other = np.ones(idx.values.shape, dtype=idx.values.dtype) * 2

        result = divmod(idx, Series(other))
        with np.errstate(all='ignore'):
            div, mod = divmod(idx.values, other)

        expected = Series(div), Series(mod)
        for r, e in zip(result, expected):
            tm.assert_series_equal(r, e)

    @pytest.mark.parametrize('other', [np.nan, 7, -23, 2.718, -3.14, np.inf])
    def test_ops_np_scalar(self, other):
        vals = np.random.randn(5, 3)
        f = lambda x: pd.DataFrame(x, index=list('ABCDE'),
                                   columns=['jim', 'joe', 'jolie'])

        df = f(vals)

        tm.assert_frame_equal(df / np.array(other), f(vals / other))
        tm.assert_frame_equal(np.array(other) * df, f(vals * other))
        tm.assert_frame_equal(df + np.array(other), f(vals + other))
        tm.assert_frame_equal(np.array(other) - df, f(other - vals))

    # TODO: This came from series.test.test_operators, needs cleanup
    def test_operators_frame(self):
        # rpow does not work with DataFrame
        ts = tm.makeTimeSeries()
        ts.name = 'ts'

        df = pd.DataFrame({'A': ts})

        tm.assert_series_equal(ts + ts, ts + df['A'],
                               check_names=False)
        tm.assert_series_equal(ts ** ts, ts ** df['A'],
                               check_names=False)
        tm.assert_series_equal(ts < ts, ts < df['A'],
                               check_names=False)
        tm.assert_series_equal(ts / ts, ts / df['A'],
                               check_names=False)


class TestAdditionSubtraction(object):
    # __add__, __sub__, __radd__, __rsub__, __iadd__, __isub__
    # for non-timestamp/timedelta/period dtypes

    # TODO: This came from series.test.test_operators, needs cleanup
    def test_arith_ops_df_compat(self):
        # GH#1134
        s1 = pd.Series([1, 2, 3], index=list('ABC'), name='x')
        s2 = pd.Series([2, 2, 2], index=list('ABD'), name='x')

        exp = pd.Series([3.0, 4.0, np.nan, np.nan],
                        index=list('ABCD'), name='x')
        tm.assert_series_equal(s1 + s2, exp)
        tm.assert_series_equal(s2 + s1, exp)

        exp = pd.DataFrame({'x': [3.0, 4.0, np.nan, np.nan]},
                           index=list('ABCD'))
        tm.assert_frame_equal(s1.to_frame() + s2.to_frame(), exp)
        tm.assert_frame_equal(s2.to_frame() + s1.to_frame(), exp)

        # different length
        s3 = pd.Series([1, 2, 3], index=list('ABC'), name='x')
        s4 = pd.Series([2, 2, 2, 2], index=list('ABCD'), name='x')

        exp = pd.Series([3, 4, 5, np.nan],
                        index=list('ABCD'), name='x')
        tm.assert_series_equal(s3 + s4, exp)
        tm.assert_series_equal(s4 + s3, exp)

        exp = pd.DataFrame({'x': [3, 4, 5, np.nan]},
                           index=list('ABCD'))
        tm.assert_frame_equal(s3.to_frame() + s4.to_frame(), exp)
        tm.assert_frame_equal(s4.to_frame() + s3.to_frame(), exp)

    # TODO: This came from series.test.test_operators, needs cleanup
    def test_series_frame_radd_bug(self):
        # GH#353
        vals = pd.Series(tm.rands_array(5, 10))
        result = 'foo_' + vals
        expected = vals.map(lambda x: 'foo_' + x)
        tm.assert_series_equal(result, expected)

        frame = pd.DataFrame({'vals': vals})
        result = 'foo_' + frame
        expected = pd.DataFrame({'vals': vals.map(lambda x: 'foo_' + x)})
        tm.assert_frame_equal(result, expected)

        ts = tm.makeTimeSeries()
        ts.name = 'ts'

        # really raise this time
        now = pd.Timestamp.now().to_pydatetime()
        with pytest.raises(TypeError):
            now + ts

        with pytest.raises(TypeError):
            ts + now

    # TODO: This came from series.test.test_operators, needs cleanup
    def test_datetime64_with_index(self):
        # arithmetic integer ops with an index
        ser = pd.Series(np.random.randn(5))
        expected = ser - ser.index.to_series()
        result = ser - ser.index
        tm.assert_series_equal(result, expected)

        # GH#4629
        # arithmetic datetime64 ops with an index
        ser = pd.Series(pd.date_range('20130101', periods=5),
                        index=pd.date_range('20130101', periods=5))
        expected = ser - ser.index.to_series()
        result = ser - ser.index
        tm.assert_series_equal(result, expected)

        with pytest.raises(TypeError):
            # GH#18850
            result = ser - ser.index.to_period()

        df = pd.DataFrame(np.random.randn(5, 2),
                          index=pd.date_range('20130101', periods=5))
        df['date'] = pd.Timestamp('20130102')
        df['expected'] = df['date'] - df.index.to_series()
        df['result'] = df['date'] - df.index
        tm.assert_series_equal(df['result'], df['expected'], check_names=False)

    # TODO: taken from tests.frame.test_operators, needs cleanup
    def test_frame_operators(self):
        seriesd = tm.getSeriesData()
        frame = pd.DataFrame(seriesd)
        frame2 = pd.DataFrame(seriesd, columns=['D', 'C', 'B', 'A'])

        garbage = np.random.random(4)
        colSeries = pd.Series(garbage, index=np.array(frame.columns))

        idSum = frame + frame
        seriesSum = frame + colSeries

        for col, series in idSum.items():
            for idx, val in series.items():
                origVal = frame[col][idx] * 2
                if not np.isnan(val):
                    assert val == origVal
                else:
                    assert np.isnan(origVal)

        for col, series in seriesSum.items():
            for idx, val in series.items():
                origVal = frame[col][idx] + colSeries[col]
                if not np.isnan(val):
                    assert val == origVal
                else:
                    assert np.isnan(origVal)

        added = frame2 + frame2
        expected = frame2 * 2
        tm.assert_frame_equal(added, expected)

        df = pd.DataFrame({'a': ['a', None, 'b']})
        tm.assert_frame_equal(df + df,
                              pd.DataFrame({'a': ['aa', np.nan, 'bb']}))

        # Test for issue #10181
        for dtype in ('float', 'int64'):
            frames = [
                pd.DataFrame(dtype=dtype),
                pd.DataFrame(columns=['A'], dtype=dtype),
                pd.DataFrame(index=[0], dtype=dtype),
            ]
            for df in frames:
                assert (df + df).equals(df)
                tm.assert_frame_equal(df + df, df)

    # TODO: taken from tests.series.test_operators; needs cleanup
    def test_series_operators(self):
        def _check_op(series, other, op, pos_only=False, check_dtype=True):
            left = np.abs(series) if pos_only else series
            right = np.abs(other) if pos_only else other

            cython_or_numpy = op(left, right)
            python = left.combine(right, op)
            tm.assert_series_equal(cython_or_numpy, python,
                                   check_dtype=check_dtype)

        def check(series, other):
            simple_ops = ['add', 'sub', 'mul', 'truediv', 'floordiv', 'mod']

            for opname in simple_ops:
                _check_op(series, other, getattr(operator, opname))

            _check_op(series, other, operator.pow, pos_only=True)

            _check_op(series, other, lambda x, y: operator.add(y, x))
            _check_op(series, other, lambda x, y: operator.sub(y, x))
            _check_op(series, other, lambda x, y: operator.truediv(y, x))
            _check_op(series, other, lambda x, y: operator.floordiv(y, x))
            _check_op(series, other, lambda x, y: operator.mul(y, x))
            _check_op(series, other, lambda x, y: operator.pow(y, x),
                      pos_only=True)
            _check_op(series, other, lambda x, y: operator.mod(y, x))

        tser = tm.makeTimeSeries().rename('ts')
        check(tser, tser * 2)
        check(tser, tser * 0)
        check(tser, tser[::2])
        check(tser, 5)

        def check_comparators(series, other, check_dtype=True):
            _check_op(series, other, operator.gt, check_dtype=check_dtype)
            _check_op(series, other, operator.ge, check_dtype=check_dtype)
            _check_op(series, other, operator.eq, check_dtype=check_dtype)
            _check_op(series, other, operator.lt, check_dtype=check_dtype)
            _check_op(series, other, operator.le, check_dtype=check_dtype)

        check_comparators(tser, 5)
        check_comparators(tser, tser + 1, check_dtype=False)

    # TODO: taken from tests.series.test_operators; needs cleanup
    def test_divmod(self):
        def check(series, other):
            results = divmod(series, other)
            if isinstance(other, Iterable) and len(series) != len(other):
                # if the lengths don't match, this is the test where we use
                # `tser[::2]`. Pad every other value in `other_np` with nan.
                other_np = []
                for n in other:
                    other_np.append(n)
                    other_np.append(np.nan)
            else:
                other_np = other
            other_np = np.asarray(other_np)
            with np.errstate(all='ignore'):
                expecteds = divmod(series.values, np.asarray(other_np))

            for result, expected in zip(results, expecteds):
                # check the values, name, and index separately
                tm.assert_almost_equal(np.asarray(result), expected)

                assert result.name == series.name
                tm.assert_index_equal(result.index, series.index)

        tser = tm.makeTimeSeries().rename('ts')
        check(tser, tser * 2)
        check(tser, tser * 0)
        check(tser, tser[::2])
        check(tser, 5)


class TestUFuncCompat(object):
    @pytest.mark.parametrize('holder', [pd.Int64Index, pd.UInt64Index,
                                        pd.Float64Index, pd.Series])
    def test_ufunc_coercions(self, holder):
        idx = holder([1, 2, 3, 4, 5], name='x')
        box = pd.Series if holder is pd.Series else pd.Index

        result = np.sqrt(idx)
        assert result.dtype == 'f8' and isinstance(result, box)
        exp = pd.Float64Index(np.sqrt(np.array([1, 2, 3, 4, 5])), name='x')
        exp = tm.box_expected(exp, box)
        tm.assert_equal(result, exp)

        result = np.divide(idx, 2.)
        assert result.dtype == 'f8' and isinstance(result, box)
        exp = pd.Float64Index([0.5, 1., 1.5, 2., 2.5], name='x')
        exp = tm.box_expected(exp, box)
        tm.assert_equal(result, exp)

        # _evaluate_numeric_binop
        result = idx + 2.
        assert result.dtype == 'f8' and isinstance(result, box)
        exp = pd.Float64Index([3., 4., 5., 6., 7.], name='x')
        exp = tm.box_expected(exp, box)
        tm.assert_equal(result, exp)

        result = idx - 2.
        assert result.dtype == 'f8' and isinstance(result, box)
        exp = pd.Float64Index([-1., 0., 1., 2., 3.], name='x')
        exp = tm.box_expected(exp, box)
        tm.assert_equal(result, exp)

        result = idx * 1.
        assert result.dtype == 'f8' and isinstance(result, box)
        exp = pd.Float64Index([1., 2., 3., 4., 5.], name='x')
        exp = tm.box_expected(exp, box)
        tm.assert_equal(result, exp)

        result = idx / 2.
        assert result.dtype == 'f8' and isinstance(result, box)
        exp = pd.Float64Index([0.5, 1., 1.5, 2., 2.5], name='x')
        exp = tm.box_expected(exp, box)
        tm.assert_equal(result, exp)

    @pytest.mark.parametrize('holder', [pd.Int64Index, pd.UInt64Index,
                                        pd.Float64Index, pd.Series])
    def test_ufunc_multiple_return_values(self, holder):
        obj = holder([1, 2, 3], name='x')
        box = pd.Series if holder is pd.Series else pd.Index

        result = np.modf(obj)
        assert isinstance(result, tuple)
        exp1 = pd.Float64Index([0., 0., 0.], name='x')
        exp2 = pd.Float64Index([1., 2., 3.], name='x')
        tm.assert_equal(result[0], tm.box_expected(exp1, box))
        tm.assert_equal(result[1], tm.box_expected(exp2, box))

    def test_ufunc_at(self):
        s = pd.Series([0, 1, 2], index=[1, 2, 3], name='x')
        np.add.at(s, [0, 2], 10)
        expected = pd.Series([10, 1, 12], index=[1, 2, 3], name='x')
        tm.assert_series_equal(s, expected)


class TestObjectDtypeEquivalence(object):
    # Tests that arithmetic operations match operations executed elementwise

    @pytest.mark.parametrize('dtype', [None, object])
    def test_numarr_with_dtype_add_nan(self, dtype, box):
        ser = pd.Series([1, 2, 3], dtype=dtype)
        expected = pd.Series([np.nan, np.nan, np.nan], dtype=dtype)

        ser = tm.box_expected(ser, box)
        expected = tm.box_expected(expected, box)

        result = np.nan + ser
        tm.assert_equal(result, expected)

        result = ser + np.nan
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize('dtype', [None, object])
    def test_numarr_with_dtype_add_int(self, dtype, box):
        ser = pd.Series([1, 2, 3], dtype=dtype)
        expected = pd.Series([2, 3, 4], dtype=dtype)

        ser = tm.box_expected(ser, box)
        expected = tm.box_expected(expected, box)

        result = 1 + ser
        tm.assert_equal(result, expected)

        result = ser + 1
        tm.assert_equal(result, expected)

    # TODO: moved from tests.series.test_operators; needs cleanup
    @pytest.mark.parametrize('op', [operator.add, operator.sub, operator.mul,
                                    operator.truediv, operator.floordiv])
    def test_operators_reverse_object(self, op):
        # GH#56
        arr = pd.Series(np.random.randn(10), index=np.arange(10), dtype=object)

        result = op(1., arr)
        expected = op(1., arr.astype(float))
        tm.assert_series_equal(result.astype(float), expected)
