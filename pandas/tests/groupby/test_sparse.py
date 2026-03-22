import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    Series,
    SparseDtype,
)
import pandas._testing as tm


@pytest.fixture(params=[0, np.nan])
def fill_value(request):
    return request.param


@pytest.fixture
def sparse_df(fill_value):
    """DataFrame with sparse int columns and a dense grouping column."""
    dense = DataFrame(
        {
            "key": ["a", "a", "b", "b", "a", "b"],
            "val1": [1, 0, 0, 3, 0, 5],
            "val2": [0, 2, 0, 0, 4, 0],
        }
    )
    sparse = dense.copy()
    sparse["val1"] = sparse["val1"].astype(SparseDtype(np.int64, fill_value))
    sparse["val2"] = sparse["val2"].astype(SparseDtype(np.int64, fill_value))
    return dense, sparse


class TestSparseGroupby:
    """Tests for groupby operations on SparseArray columns (GH#36123)."""

    @pytest.mark.parametrize("op", ["sum", "mean", "min", "max", "std", "var"])
    def test_sparse_groupby_agg(self, sparse_df, op):
        dense, sparse = sparse_df
        result = getattr(sparse.groupby("key"), op)()
        expected = getattr(dense.groupby("key"), op)()
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("op", ["first", "last"])
    def test_sparse_groupby_first_last(self, sparse_df, op):
        # first/last preserve the SparseArray dtype
        dense, sparse = sparse_df
        result = getattr(sparse.groupby("key"), op)()
        expected = getattr(dense.groupby("key"), op)()
        tm.assert_frame_equal(result, expected, check_dtype=False)

    @pytest.mark.parametrize("op", ["any", "all"])
    def test_sparse_groupby_any_all(self, sparse_df, op):
        dense, sparse = sparse_df
        result = getattr(sparse.groupby("key"), op)()
        expected = getattr(dense.groupby("key"), op)()
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("op", ["sem", "prod", "median"])
    def test_sparse_groupby_other_aggs(self, sparse_df, op):
        dense, sparse = sparse_df
        result = getattr(sparse.groupby("key"), op)()
        expected = getattr(dense.groupby("key"), op)()
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("op", ["cumsum", "cummin", "cummax", "cumprod"])
    def test_sparse_groupby_transform(self, sparse_df, op):
        dense, sparse = sparse_df
        result = getattr(sparse.groupby("key"), op)()
        expected = getattr(dense.groupby("key"), op)()
        tm.assert_frame_equal(result, expected)

    def test_sparse_groupby_rank(self, sparse_df):
        dense, sparse = sparse_df
        result = sparse.groupby("key").rank()
        expected = dense.groupby("key").rank()
        tm.assert_frame_equal(result, expected)

    def test_sparse_groupby_idxmin_idxmax(self, sparse_df):
        dense, sparse = sparse_df
        for op in ["idxmin", "idxmax"]:
            result = getattr(sparse.groupby("key"), op)()
            expected = getattr(dense.groupby("key"), op)()
            tm.assert_frame_equal(result, expected)

    def test_sparse_groupby_nan_fill_value(self):
        # When fill_value=NaN, gap positions are NaN and should be
        # excluded from aggregations.
        dense = DataFrame(
            {
                "key": ["a", "a", "b", "b"],
                "val": [1.0, np.nan, np.nan, 4.0],
            }
        )
        sparse = dense.copy()
        sparse["val"] = pd.array(dense["val"], dtype=SparseDtype(float, np.nan))
        result = sparse.groupby("key").mean()
        expected = dense.groupby("key").mean()
        tm.assert_frame_equal(result, expected)

    def test_sparse_groupby_series(self, fill_value):
        # Groupby on a single sparse Series.
        vals = [1, 0, 0, 3, 0, 5]
        keys = ["a", "a", "b", "b", "a", "b"]
        dense_ser = Series(vals, name="val")
        sparse_ser = dense_ser.astype(SparseDtype(np.int64, fill_value))
        result = sparse_ser.groupby(keys).mean()
        expected = dense_ser.groupby(keys).mean()
        tm.assert_series_equal(result, expected)
