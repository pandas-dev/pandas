# -*- coding: utf-8 -*-
# Arithmetc tests for DataFrame/Series/Index/Array classes that should
# behave identically.
# Specifically for object dtype

import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm

from pandas import Series, Timestamp


# ------------------------------------------------------------------
# Comparisons


# ------------------------------------------------------------------
# Arithmetic

class TestArithmetic(object):
    def test_df_radd_str(self):
        df = pd.DataFrame(['x', np.nan, 'x'])

        expected = pd.DataFrame(['ax', np.nan, 'ax'])
        result = 'a' + df
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame(['xa', np.nan, 'xa'])
        result = df + 'a'
        tm.assert_frame_equal(result, expected)

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
    def test_df_radd_str_invalid(self, dtype, data):
        df = pd.DataFrame(data, dtype=dtype)
        with pytest.raises(TypeError):
            'foo_' + df

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
