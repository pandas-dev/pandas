# -*- coding: utf-8 -*-
import operator
from datetime import timedelta

import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm
from pandas import Series, compat
from pandas.core.indexes.period import IncompatibleFrequency


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

    def test_empty_series_add_sub(self):
        # GH#13844
        a = Series(dtype='M8[ns]')
        b = Series(dtype='m8[ns]')
        tm.assert_series_equal(a, a + b)
        tm.assert_series_equal(a, a - b)
        tm.assert_series_equal(a, b + a)
        with pytest.raises(TypeError):
            b - a

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
        with tm.assert_raises_regex(IncompatibleFrequency, msg):
            ts + ts.asfreq('D', how="end")

    def test_operators_datetimelike(self):

        # ## timedelta64 ###
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        # ## datetime64 ###
        dt1 = Series([pd.Timestamp('20111230'), pd.Timestamp('20120101'),
                      pd.Timestamp('20120103')])
        dt1.iloc[2] = np.nan
        dt2 = Series([pd.Timestamp('20111231'), pd.Timestamp('20120102'),
                      pd.Timestamp('20120104')])
        dt1 - dt2
        dt2 - dt1

        # ## datetime64 with timetimedelta ###
        dt1 + td1
        td1 + dt1
        dt1 - td1
        # TODO: Decide if this ought to work.
        # td1 - dt1

        # ## timetimedelta with datetime64 ###
        td1 + dt1
        dt1 + td1


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
            with tm.assert_raises_regex(ValueError, msg):
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


def test_pow_ops_object():
    # 22922
    # pow is weird with masking & 1, so testing here
    a = Series([1, np.nan, 1, np.nan], dtype=object)
    b = Series([1, np.nan, np.nan, 1], dtype=object)
    result = a ** b
    expected = Series(a.values ** b.values, dtype=object)
    tm.assert_series_equal(result, expected)

    result = b ** a
    expected = Series(b.values ** a.values, dtype=object)

    tm.assert_series_equal(result, expected)
