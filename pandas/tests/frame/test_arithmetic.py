# -*- coding: utf-8 -*-

import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm


# -------------------------------------------------------------------
# Comparisons

class TestFrameComparisons(object):
    def test_df_boolean_comparison_error(self):
        # GH#4576
        # boolean comparisons with a tuple/list give unexpected results
        df = pd.DataFrame(np.arange(6).reshape((3, 2)))

        # not shape compatible
        with pytest.raises(ValueError):
            df == (2, 2)
        with pytest.raises(ValueError):
            df == [2, 2]

    def test_df_float_none_comparison(self):
        df = pd.DataFrame(np.random.randn(8, 3), index=range(8),
                          columns=['A', 'B', 'C'])

        with pytest.raises(TypeError):
            df.__eq__(None)

    def test_df_string_comparison(self):
        df = pd.DataFrame([{"a": 1, "b": "foo"}, {"a": 2, "b": "bar"}])
        mask_a = df.a > 1
        tm.assert_frame_equal(df[mask_a], df.loc[1:1, :])
        tm.assert_frame_equal(df[-mask_a], df.loc[0:0, :])

        mask_b = df.b == "foo"
        tm.assert_frame_equal(df[mask_b], df.loc[0:0, :])
        tm.assert_frame_equal(df[-mask_b], df.loc[1:1, :])

    @pytest.mark.parametrize('opname', ['eq', 'ne', 'gt', 'lt', 'ge', 'le'])
    def test_df_flex_cmp_constant_return_types(self, opname):
        # GH#15077, non-empty DataFrame
        df = pd.DataFrame({'x': [1, 2, 3], 'y': [1., 2., 3.]})
        const = 2

        result = getattr(df, opname)(const).get_dtype_counts()
        tm.assert_series_equal(result, pd.Series([2], ['bool']))

    @pytest.mark.parametrize('opname', ['eq', 'ne', 'gt', 'lt', 'ge', 'le'])
    def test_df_flex_cmp_constant_return_types_empty(self, opname):
        # GH#15077 empty DataFrame
        df = pd.DataFrame({'x': [1, 2, 3], 'y': [1., 2., 3.]})
        const = 2

        empty = df.iloc[:0]
        result = getattr(empty, opname)(const).get_dtype_counts()
        tm.assert_series_equal(result, pd.Series([2], ['bool']))


# -------------------------------------------------------------------
# Arithmetic

class TestPeriodFrameArithmetic(object):

    def test_ops_frame_period(self):
        # GH 13043
        df = pd.DataFrame({'A': [pd.Period('2015-01', freq='M'),
                                 pd.Period('2015-02', freq='M')],
                           'B': [pd.Period('2014-01', freq='M'),
                                 pd.Period('2014-02', freq='M')]})
        assert df['A'].dtype == object
        assert df['B'].dtype == object

        p = pd.Period('2015-03', freq='M')
        # dtype will be object because of original dtype
        exp = pd.DataFrame({'A': np.array([2, 1], dtype=object),
                            'B': np.array([14, 13], dtype=object)})
        tm.assert_frame_equal(p - df, exp)
        tm.assert_frame_equal(df - p, -exp)

        df2 = pd.DataFrame({'A': [pd.Period('2015-05', freq='M'),
                                  pd.Period('2015-06', freq='M')],
                            'B': [pd.Period('2015-05', freq='M'),
                                  pd.Period('2015-06', freq='M')]})
        assert df2['A'].dtype == object
        assert df2['B'].dtype == object

        exp = pd.DataFrame({'A': np.array([4, 4], dtype=object),
                            'B': np.array([16, 16], dtype=object)})
        tm.assert_frame_equal(df2 - df, exp)
        tm.assert_frame_equal(df - df2, -exp)
