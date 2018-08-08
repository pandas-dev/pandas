# -*- coding: utf-8 -*-
from datetime import timedelta
import operator
from decimal import Decimal

import numpy as np
import pytest

from pandas import Series, Timestamp

import pandas as pd
import pandas.util.testing as tm


# ------------------------------------------------------------------
# Comparisons

class TestSeriesComparison(object):
    def test_compare_invalid(self):
        # GH#8058
        # ops testing
        a = pd.Series(np.random.randn(5), name=0)
        b = pd.Series(np.random.randn(5))
        b.name = pd.Timestamp('2000-01-01')
        tm.assert_series_equal(a / b, 1 / (b / a))

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


class TestTimedeltaSeriesComparisons(object):
    def test_compare_timedelta_series(self):
        # regresssion test for GH5963
        s = pd.Series([timedelta(days=1), timedelta(days=2)])
        actual = s > timedelta(days=1)
        expected = pd.Series([False, True])
        tm.assert_series_equal(actual, expected)


# ------------------------------------------------------------------
# Arithmetic

class TestSeriesDivision(object):
    # __div__, __rdiv__, __floordiv__, __rfloordiv__
    # for non-timestamp/timedelta/period dtypes

    def test_divide_decimal(self):
        # resolves issue GH#9787
        expected = Series([Decimal(5)])

        ser = Series([Decimal(10)])
        result = ser / Decimal(2)

        tm.assert_series_equal(result, expected)

        ser = Series([Decimal(10)])
        result = ser // Decimal(2)

        tm.assert_series_equal(result, expected)

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


class TestSeriesArithmetic(object):
    # Standard, numeric, or otherwise not-Timestamp/Timedelta/Period dtypes
    @pytest.mark.parametrize('data', [
        [1, 2, 3],
        [1.1, 2.2, 3.3],
        [Timestamp('2011-01-01'), Timestamp('2011-01-02'), pd.NaT],
        ['x', 'y', 1]])
    @pytest.mark.parametrize('dtype', [None, object])
    def test_series_radd_str_invalid(self, dtype, data):
        ser = Series(data, dtype=dtype)
        with pytest.raises(TypeError):
            'foo_' + ser

    # TODO: parametrize, better name
    def test_object_ser_add_invalid(self):
        # invalid ops
        obj_ser = tm.makeObjectSeries()
        obj_ser.name = 'objects'
        with pytest.raises(Exception):
            obj_ser + 1
        with pytest.raises(Exception):
            obj_ser + np.array(1, dtype=np.int64)
        with pytest.raises(Exception):
            obj_ser - 1
        with pytest.raises(Exception):
            obj_ser - np.array(1, dtype=np.int64)

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

    def test_series_radd_str(self):
        ser = pd.Series(['x', np.nan, 'x'])
        tm.assert_series_equal('a' + ser, pd.Series(['ax', np.nan, 'ax']))
        tm.assert_series_equal(ser + 'a', pd.Series(['xa', np.nan, 'xa']))

    @pytest.mark.parametrize('dtype', [None, object])
    def test_series_with_dtype_radd_timedelta(self, dtype):
        # note this test is _not_ aimed at timedelta64-dtyped Series
        ser = pd.Series([pd.Timedelta('1 days'), pd.Timedelta('2 days'),
                         pd.Timedelta('3 days')], dtype=dtype)
        expected = pd.Series([pd.Timedelta('4 days'), pd.Timedelta('5 days'),
                              pd.Timedelta('6 days')])

        result = pd.Timedelta('3 days') + ser
        tm.assert_series_equal(result, expected)

        result = ser + pd.Timedelta('3 days')
        tm.assert_series_equal(result, expected)
