# -*- coding: utf-8 -*-
import pytest
import numpy as np

from pandas.compat import range

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

class TestFrameFlexArithmetic(object):
    def test_df_add_flex_filled_mixed_dtypes(self):
        # GH#19611
        dti = pd.date_range('2016-01-01', periods=3)
        ser = pd.Series(['1 Day', 'NaT', '2 Days'], dtype='timedelta64[ns]')
        df = pd.DataFrame({'A': dti, 'B': ser})
        other = pd.DataFrame({'A': ser, 'B': ser})
        fill = pd.Timedelta(days=1).to_timedelta64()
        result = df.add(other, fill_value=fill)

        expected = pd.DataFrame(
            {'A': pd.Series(['2016-01-02', '2016-01-03', '2016-01-05'],
                            dtype='datetime64[ns]'),
             'B': ser * 2})
        tm.assert_frame_equal(result, expected)
