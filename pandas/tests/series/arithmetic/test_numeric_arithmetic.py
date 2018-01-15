# coding=utf-8

import pytest

import numpy as np
import pandas as pd

from pandas import Series
import pandas.util.testing as tm


class TestSeriesArithmeticDtypeCompat(object):

    def test_object_series_add_int_invalid(self):
        # invalid ops
        obj_series = tm.makeObjectSeries()
        obj_series.name = 'objects'

        with pytest.raises(Exception):
            obj_series + 1
        with pytest.raises(Exception):
            obj_series + np.array(1, dtype=np.int64)

    def test_object_series_sub_int_invalid(self):
        # invalid ops
        obj_series = tm.makeObjectSeries()
        obj_series.name = 'objects'

        with pytest.raises(Exception):
            obj_series - 1
        with pytest.raises(Exception):
            obj_series - np.array(1, dtype=np.int64)

    def test_series_radd_str(self):
        ser = pd.Series(['x', np.nan, 'x'])
        tm.assert_series_equal('a' + ser, pd.Series(['ax', np.nan, 'ax']))
        tm.assert_series_equal(ser + 'a', pd.Series(['xa', np.nan, 'xa']))

    @pytest.mark.parametrize('data', [
        [1, 2, 3],
        [1.1, 2.2, 3.3],
        [pd.Timestamp('2011-01-01'), pd.Timestamp('2011-01-02'), pd.NaT],
        ['x', 'y', 1]])
    @pytest.mark.parametrize('dtype', [None, object])
    def test_series_radd_str_invalid(self, dtype, data):
        ser = Series(data, dtype=dtype)
        with pytest.raises(TypeError):
            'foo_' + ser

    @pytest.mark.parametrize('dtype', [None, object])
    def test_series_with_dtype_radd_int(self, dtype):
        ser = pd.Series([1, 2, 3], dtype=dtype)
        expected = pd.Series([2, 3, 4], dtype=dtype)

        result = 1 + ser
        tm.assert_series_equal(result, expected)

        result = ser + 1
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('dtype', [None, object])
    def test_series_with_dtype_radd_nan(self, dtype):
        ser = pd.Series([1, 2, 3], dtype=dtype)
        expected = pd.Series([np.nan, np.nan, np.nan], dtype=dtype)

        result = np.nan + ser
        tm.assert_series_equal(result, expected)

        result = ser + np.nan
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('dtype', [None, object])
    def test_series_with_dtype_radd_timedelta(self, dtype):
        ser = pd.Series([pd.Timedelta('1 days'), pd.Timedelta('2 days'),
                         pd.Timedelta('3 days')], dtype=dtype)
        expected = pd.Series([pd.Timedelta('4 days'), pd.Timedelta('5 days'),
                              pd.Timedelta('6 days')])

        result = pd.Timedelta('3 days') + ser
        tm.assert_series_equal(result, expected)

        result = ser + pd.Timedelta('3 days')
        tm.assert_series_equal(result, expected)


class TestSeriesArithmetic(object):
    def test_divide_decimal(self):
        """ resolves issue #9787 """
        from decimal import Decimal

        expected = Series([Decimal(5)])

        s = Series([Decimal(10)])
        s = s / Decimal(2)

        tm.assert_series_equal(expected, s)

        s = Series([Decimal(10)])
        s = s // Decimal(2)

        tm.assert_series_equal(expected, s)

    def test_div(self):
        with np.errstate(all='ignore'):
            # no longer do integer div for any ops, but deal with the 0's
            p = pd.DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})
            result = p['first'] / p['second']
            expected = Series(
                p['first'].values.astype(float) / p['second'].values,
                dtype='float64')
            expected.iloc[0:3] = np.inf
            tm.assert_series_equal(result, expected)

            result = p['first'] / 0
            expected = Series(np.inf, index=p.index, name='first')
            tm.assert_series_equal(result, expected)

            p = p.astype('float64')
            result = p['first'] / p['second']
            expected = Series(p['first'].values / p['second'].values)
            tm.assert_series_equal(result, expected)

            p = pd.DataFrame({'first': [3, 4, 5, 8], 'second': [1, 1, 1, 1]})
            result = p['first'] / p['second']
            tm.assert_series_equal(result, p['first'].astype('float64'),
                                   check_names=False)
            assert result.name is None
            assert not result.equals(p['second'] / p['first'])

            # inf signing
            s = Series([np.nan, 1., -1.])
            result = s / 0
            expected = Series([np.nan, np.inf, -np.inf])
            tm.assert_series_equal(result, expected)

            # float/integer issue
            # GH 7785
            p = pd.DataFrame({'first': (1, 0), 'second': (-0.01, -0.02)})
            expected = Series([-0.01, -np.inf])

            result = p['second'].div(p['first'])
            tm.assert_series_equal(result, expected, check_names=False)

            result = p['second'] / p['first']
            tm.assert_series_equal(result, expected)

            # GH 9144
            s = Series([-1, 0, 1])

            result = 0 / s
            expected = Series([0.0, np.nan, 0.0])
            tm.assert_series_equal(result, expected)

            result = s / 0
            expected = Series([-np.inf, np.nan, np.inf])
            tm.assert_series_equal(result, expected)

            result = s // 0
            expected = Series([-np.inf, np.nan, np.inf])
            tm.assert_series_equal(result, expected)

            # GH 8674
            zero_array = np.array([0] * 5)
            data = np.random.randn(5)
            expected = pd.Series([0.] * 5)
            result = zero_array / pd.Series(data)
            tm.assert_series_equal(result, expected)

            result = pd.Series(zero_array) / data
            tm.assert_series_equal(result, expected)

            result = pd.Series(zero_array) / pd.Series(data)
            tm.assert_series_equal(result, expected)
