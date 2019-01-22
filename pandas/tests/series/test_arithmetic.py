# -*- coding: utf-8 -*-
import operator

import numpy as np
import pytest

import pandas as pd
from pandas import Series, compat
from pandas.core.indexes.period import IncompatibleFrequency
import pandas.util.testing as tm


def _permute(obj):
    return obj.take(np.random.permutation(len(obj)))


class TestSeriesFlexArithmetic(object):
    @pytest.mark.parametrize(
        'ts',
        [
            (lambda x: x, lambda x: x * 2, False),
            (lambda x: x, lambda x: x[::2], False),
            (lambda x: x, lambda x: 5, True),
            (lambda x: tm.makeFloatSeries(),
             lambda x: tm.makeFloatSeries(),
             True)
        ])
    @pytest.mark.parametrize('opname', ['add', 'sub', 'mul', 'floordiv',
                                        'truediv', 'div', 'pow'])
    def test_flex_method_equivalence(self, opname, ts):
        # check that Series.{opname} behaves like Series.__{opname}__,
        tser = tm.makeTimeSeries().rename('ts')

        series = ts[0](tser)
        other = ts[1](tser)
        check_reverse = ts[2]

        if opname == 'div' and compat.PY3:
            pytest.skip('div test only for Py3')

        op = getattr(Series, opname)

        if op == 'div':
            alt = operator.truediv
        else:
            alt = getattr(operator, opname)

        result = op(series, other)
        expected = alt(series, other)
        tm.assert_almost_equal(result, expected)
        if check_reverse:
            rop = getattr(Series, "r" + opname)
            result = rop(series, other)
            expected = alt(other, series)
            tm.assert_almost_equal(result, expected)


class TestSeriesArithmetic(object):
    # Some of these may end up in tests/arithmetic, but are not yet sorted

    def test_add_series_with_period_index(self):
        rng = pd.period_range('1/1/2000', '1/1/2010', freq='A')
        ts = Series(np.random.randn(len(rng)), index=rng)

        result = ts + ts[::2]
        expected = ts + ts
        expected[1::2] = np.nan
        tm.assert_series_equal(result, expected)

        result = ts + _permute(ts[::2])
        tm.assert_series_equal(result, expected)

        msg = "Input has different freq=D from PeriodIndex\\(freq=A-DEC\\)"
        with pytest.raises(IncompatibleFrequency, match=msg):
            ts + ts.asfreq('D', how="end")


# ------------------------------------------------------------------
# Comparisons

class TestSeriesFlexComparison(object):
    def test_comparison_flex_basic(self):
        left = pd.Series(np.random.randn(10))
        right = pd.Series(np.random.randn(10))

        tm.assert_series_equal(left.eq(right), left == right)
        tm.assert_series_equal(left.ne(right), left != right)
        tm.assert_series_equal(left.le(right), left < right)
        tm.assert_series_equal(left.lt(right), left <= right)
        tm.assert_series_equal(left.gt(right), left > right)
        tm.assert_series_equal(left.ge(right), left >= right)

        # axis
        for axis in [0, None, 'index']:
            tm.assert_series_equal(left.eq(right, axis=axis), left == right)
            tm.assert_series_equal(left.ne(right, axis=axis), left != right)
            tm.assert_series_equal(left.le(right, axis=axis), left < right)
            tm.assert_series_equal(left.lt(right, axis=axis), left <= right)
            tm.assert_series_equal(left.gt(right, axis=axis), left > right)
            tm.assert_series_equal(left.ge(right, axis=axis), left >= right)

        #
        msg = 'No axis named 1 for object type'
        for op in ['eq', 'ne', 'le', 'le', 'gt', 'ge']:
            with pytest.raises(ValueError, match=msg):
                getattr(left, op)(right, axis=1)


class TestSeriesComparison(object):
    def test_comparison_different_length(self):
        a = Series(['a', 'b', 'c'])
        b = Series(['b', 'a'])
        with pytest.raises(ValueError):
            a < b

        a = Series([1, 2])
        b = Series([2, 3, 4])
        with pytest.raises(ValueError):
            a == b

    @pytest.mark.parametrize('opname', ['eq', 'ne', 'gt', 'lt', 'ge', 'le'])
    def test_ser_flex_cmp_return_dtypes(self, opname):
        # GH#15115
        ser = Series([1, 3, 2], index=range(3))
        const = 2

        result = getattr(ser, opname)(const).get_dtype_counts()
        tm.assert_series_equal(result, Series([1], ['bool']))

    @pytest.mark.parametrize('opname', ['eq', 'ne', 'gt', 'lt', 'ge', 'le'])
    def test_ser_flex_cmp_return_dtypes_empty(self, opname):
        # GH#15115 empty Series case
        ser = Series([1, 3, 2], index=range(3))
        empty = ser.iloc[:0]
        const = 2

        result = getattr(empty, opname)(const).get_dtype_counts()
        tm.assert_series_equal(result, Series([1], ['bool']))

    @pytest.mark.parametrize('op', [operator.eq, operator.ne,
                                    operator.le, operator.lt,
                                    operator.ge, operator.gt])
    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('baz', 'baz', 'baz')])
    def test_ser_cmp_result_names(self, names, op):
        # datetime64 dtype
        dti = pd.date_range('1949-06-07 03:00:00',
                            freq='H', periods=5, name=names[0])
        ser = Series(dti).rename(names[1])
        result = op(ser, dti)
        assert result.name == names[2]

        # datetime64tz dtype
        dti = dti.tz_localize('US/Central')
        ser = Series(dti).rename(names[1])
        result = op(ser, dti)
        assert result.name == names[2]

        # timedelta64 dtype
        tdi = dti - dti.shift(1)
        ser = Series(tdi).rename(names[1])
        result = op(ser, tdi)
        assert result.name == names[2]

        # categorical
        if op in [operator.eq, operator.ne]:
            # categorical dtype comparisons raise for inequalities
            cidx = tdi.astype('category')
            ser = Series(cidx).rename(names[1])
            result = op(ser, cidx)
            assert result.name == names[2]

# ------------------------------------------------------------------
# Unary


class TestSeriesUnary(object):

    def test_neg(self):
        ser = tm.makeStringSeries()
        ser.name = 'series'
        tm.assert_series_equal(-ser, -1 * ser)

    def test_pos(self):
        ser = tm.makeStringSeries()
        ser.name = 'series'
        tm.assert_series_equal(+ser, ser)

    def test_invert(self):
        ser = tm.makeStringSeries()
        ser.name = 'series'
        tm.assert_series_equal(-(ser < 0), ~(ser < 0))

    @pytest.mark.parametrize('op', [operator.pos, operator.neg,
                                    operator.inv, operator.abs])
    def test_integer(self, any_numpy_ea_int_dtype, op):
        # GH#23087
        typ = any_numpy_ea_int_dtype

        ser = Series([1, 3, 2], index=range(3), dtype=typ)
        result = op(ser)
        exp = Series(op(ser.values), index=ser.index)
        tm.assert_series_equal(result, exp)

    @pytest.mark.parametrize('op', [operator.pos, operator.neg,
                                    operator.abs])
    def test_float(self, float_dtype, op):
        # GH#23087
        ser = Series([1.1, 3.1, 2.1], index=range(3), dtype=float_dtype)
        result = op(ser)
        exp = Series(op(ser.values), index=ser.index)
        tm.assert_series_equal(result, exp)

    @pytest.mark.parametrize('op', [operator.inv])
    def test_float_inv_raise(self, float_dtype, op):
        # GH#23087
        ser = Series([1.1, 3.1, 2.1], index=range(3), dtype=float_dtype)
        with pytest.raises(TypeError):
            op(ser)
        # consistent with numpy
        with pytest.raises(TypeError):
            op(ser.values)

    @pytest.mark.parametrize('typ', ['bool'])
    @pytest.mark.parametrize('op', [operator.pos, operator.neg,
                                    operator.inv, operator.abs])
    def test_bool(self, typ, op):
        # GH#23087
        ser = Series([True, False, True], index=range(3), dtype=typ)

        if op is operator.neg:
            result = op(ser)
            exp = Series(operator.inv(ser.values), index=ser.index)
            tm.assert_series_equal(result, exp)

            if pd.compat.numpy._np_version_under1p13:
                # older NumPy support this
                tm.assert_numpy_array_equal(~ser.values, op(ser.values))
            else:
                # inconsistent with NumPy
                with pytest.raises(TypeError):
                    op(ser.values)
        else:
            result = op(ser)
            exp = Series(op(ser.values), index=ser.index)
            tm.assert_series_equal(result, exp)

    @pytest.mark.parametrize('typ', ['datetime64[ns]', 'datetime64[ns, GMT]',
                                     'datetime64[ns, US/Eastern]'])
    @pytest.mark.parametrize('op', [operator.pos])
    def test_datetime_pos_raises(self, typ, op):
        # GH#23087
        ser = Series([pd.Timestamp('2011-01-01'), pd.Timestamp('2011-01-02'),
                      pd.Timestamp('2011-01-03')], index=range(3), dtype=typ)
        # pandas all raises
        with pytest.raises(TypeError):
            op(ser)

        # inconsistent with NumPy
        tm.assert_numpy_array_equal(ser.values, op(ser.values))

    @pytest.mark.parametrize('typ', ['datetime64[ns]', 'datetime64[ns, GMT]',
                                     'datetime64[ns, US/Eastern]'])
    @pytest.mark.parametrize('op', [operator.neg,
                                    operator.inv, operator.abs])
    def test_datetime_raises(self, typ, op):
        # GH#23087
        ser = Series([pd.Timestamp('2011-01-01'), pd.Timestamp('2011-01-02'),
                      pd.Timestamp('2011-01-03')], index=range(3), dtype=typ)
        # pandas all raises
        with pytest.raises(TypeError):
            op(ser)
        with pytest.raises(TypeError):
            op(ser.values)

    @pytest.mark.parametrize('typ', ['timedelta64[ns]'])
    @pytest.mark.parametrize('op', [operator.pos, operator.neg,
                                    operator.abs])
    def test_timedelta(self, typ, op):
        # GH#23087
        ser = Series([pd.Timedelta('1 days'), pd.Timedelta('2 days'),
                      pd.Timedelta('3 days')], index=range(3), dtype=typ)
        result = op(ser)
        exp = Series(op(ser.values), index=ser.index)
        tm.assert_series_equal(result, exp)

    @pytest.mark.parametrize('typ', ['timedelta64[ns]'])
    @pytest.mark.parametrize('op', [operator.inv])
    def test_timedelta_inv_raise(self, typ, op):
        # GH#23087
        ser = Series([pd.Timedelta('1 days'), pd.Timedelta('2 days'),
                      pd.Timedelta('3 days')], index=range(3), dtype=typ)
        with pytest.raises(TypeError):
            op(ser)
        # consistent with numpy
        with pytest.raises(TypeError):
            op(ser.values)

    @pytest.mark.skipif(not pd.compat.numpy._np_version_under1p16,
                        reason='NumPy 1.6 or later shows warning when op '
                        'performed against non-numeric dtype')
    @pytest.mark.parametrize('typ', ['str', 'object'])
    @pytest.mark.parametrize('op', [operator.pos])
    def test_object(self, typ, op):
        # GH#23087
        ser = Series(['a', 'b', 'c'], index=range(3), dtype=typ)
        result = op(ser)
        exp = Series(op(ser.values), index=ser.index)
        tm.assert_series_equal(result, exp)

    @pytest.mark.parametrize('typ', ['str', 'object'])
    @pytest.mark.parametrize('op', [operator.neg,
                                    operator.inv, operator.abs])
    def test_object_raises(self, typ, op):
        # GH#23087
        ser = Series(['a', 'b', 'c'], index=range(3), dtype=typ)
        with pytest.raises(TypeError):
            op(ser)
        # consistent with numpy
        with pytest.raises(TypeError):
            op(ser.values)
