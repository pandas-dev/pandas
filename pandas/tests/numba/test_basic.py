"""
Test our numba extension types work correctly
"""

import numpy as np
import pytest

from pandas.compat._optional import import_optional_dependency
import pandas.util._test_decorators as td

from pandas import (
    DataFrame,
    Index,
    Series,
)
import pandas._testing as tm

pytestmark = pytest.mark.single_cpu


@pytest.fixture
def test_df():
    return tm.makeDataFrame()


@pytest.fixture()
def test_series():
    return tm.makeFloatSeries()


@td.skip_if_no("numba")
class TestBasic:
    numba = import_optional_dependency("numba", errors="ignore")

    def test_identity(self, test_df):
        @self.numba.njit
        def f(df):
            return df

        tm.assert_frame_equal(f(test_df), test_df)

    def test_transformation(self, test_df):
        def f(df):
            new_values = (df.values - df.values.mean()) / df.values.std()
            return DataFrame(new_values, index=df.index, columns=df.columns)

        njit_f = self.numba.njit(f)
        tm.assert_frame_equal(njit_f(test_df), f(test_df))

    def test_reduce_scalar(self, test_df):
        def f(df):
            return df.values.sum()

        njit_f = self.numba.njit(f)
        tm.assert_almost_equal(njit_f(test_df), f(test_df))

    def test_reduction_axis_0(self, test_df):
        def f(df):
            vals = df.values.sum(axis=0)
            # TODO: Don't need to convert the columns to a ndarray manually after
            # https://github.com/numba/numba/issues/6803 is resolved
            cols_array = np.empty(len(df.columns), dtype=np.int64)
            for i, v in enumerate(df.columns):
                cols_array[i] = v
            return Series(vals, index=Index(cols_array))

        # Use integer column names
        # since numpy/numba doesn't support string arrays
        test_df.columns = range(len(test_df.columns))
        njit_f = self.numba.njit(f)
        tm.assert_series_equal(njit_f(test_df), f(test_df))

    def test_reduction_axis_1(self, test_df):
        def f(df):
            vals = df.values.sum(axis=1)
            return Series(vals, index=df.index)

        njit_f = self.numba.njit(f)
        tm.assert_series_equal(njit_f(test_df), f(test_df))

    def test_expand_cols(self, test_series):
        # Can't use self.numba since numba doesn't recognize that inside jitted function
        import numba

        def f(series):
            demeaned = series.values - series.values.mean()
            new_vals = np.concatenate(
                (series.values.astype(np.float64), demeaned), axis=1
            )
            # TODO: Can remove use of typed.List once that
            # becomes the default list container in numba
            new_cols = numba.typed.List(["original", "demeaned"])
            return DataFrame(new_vals, index=series.index, columns=new_cols)

        njit_f = numba.njit(f)
        test_df = DataFrame(test_series)
        tm.assert_frame_equal(njit_f(test_df), f(test_df))
