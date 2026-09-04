import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.fixture(params=[0, np.nan])
def fill_value(request):
    return request.param


@pytest.fixture
def sparse_df(fill_value):
    """DataFrame with sparse int columns and a dense grouping column."""
    dense = pd.DataFrame(
        {
            "key": ["a", "a", "b", "b", "a", "b"],
            "val1": [1, 0, 0, 3, 0, 5],
            "val2": [0, 2, 0, 0, 4, 0],
        }
    )
    sparse = dense.copy()
    sparse["val1"] = sparse["val1"].astype(pd.SparseDtype(np.int64, fill_value))
    sparse["val2"] = sparse["val2"].astype(pd.SparseDtype(np.int64, fill_value))
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
        dense = pd.DataFrame(
            {
                "key": ["a", "a", "b", "b"],
                "val": [1.0, np.nan, np.nan, 4.0],
            }
        )
        sparse = dense.copy()
        sparse["val"] = pd.array(dense["val"], dtype=pd.SparseDtype(float, np.nan))
        result = sparse.groupby("key").mean()
        expected = dense.groupby("key").mean()
        tm.assert_frame_equal(result, expected)

    def test_sparse_groupby_series(self, fill_value):
        # Groupby on a single sparse Series.
        vals = [1, 0, 0, 3, 0, 5]
        keys = ["a", "a", "b", "b", "a", "b"]
        dense_ser = pd.Series(vals, name="val")
        sparse_ser = dense_ser.astype(pd.SparseDtype(np.int64, fill_value))
        result = sparse_ser.groupby(keys).mean()
        expected = dense_ser.groupby(keys).mean()
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("op", ["min", "max", "sum", "mean", "cumsum", "any"])
    def test_sparse_groupby_nan_fill_value_gaps(self, op):
        # GH#64758 the NaN fill_value must be materialized rather than
        #  unsafely cast to the int64 subtype
        keys = ["a", "a", "b", "b"]
        sparse_ser = pd.Series([1, 0, 1, 0]).astype("Sparse[int64]")
        sparse_ser = sparse_ser.astype(pd.SparseDtype(np.int64, np.nan))
        assert sparse_ser.array.sp_index.ngaps == 2

        dense_ser = pd.Series([1.0, np.nan, 1.0, np.nan])
        result = getattr(sparse_ser.groupby(keys), op)()
        expected = getattr(dense_ser.groupby(keys), op)()
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("op", ["sum", "any", "all"])
    def test_sparse_groupby_bool_nan_fill_value(self, op):
        # GH#64758 a NaN fill_value on a bool subtype was cast to True; the gaps
        #  are missing values, so the dense equivalent is object dtype
        keys = ["a", "a", "b", "b"]
        sparse_ser = pd.Series([True, False, True, False]).astype(
            pd.SparseDtype(bool, False)
        )
        sparse_ser = sparse_ser.astype(pd.SparseDtype(bool, np.nan))
        assert sparse_ser.array.sp_index.ngaps == 2

        dense_ser = pd.Series([True, np.nan, True, np.nan])
        result = getattr(sparse_ser.groupby(keys), op)()
        expected = getattr(dense_ser.groupby(keys), op)()
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "subtype", ["int8", "int32", "uint16", "uint64", "float32"]
    )
    @pytest.mark.parametrize("op", ["min", "max", "sum", "cumsum"])
    def test_sparse_groupby_preserves_subtype(self, subtype, op):
        # GH#64758 densifying must give the subtype's own dense dtype; widening
        #  uint64 to float64 would round values above 2**53
        keys = ["a", "a", "b", "b"]
        big = 2**63 + 12345 if subtype == "uint64" else 3
        dense_ser = pd.Series(np.array([1, 0, big, 0], dtype=subtype))
        sparse_ser = dense_ser.astype(pd.SparseDtype(subtype, 0))
        assert sparse_ser.array.sp_index.ngaps == 2

        result = getattr(sparse_ser.groupby(keys), op)()
        expected = getattr(dense_ser.groupby(keys), op)()
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("op", ["min", "max", "sum", "cumsum"])
    def test_sparse_groupby_wider_fill_value(self, op):
        # GH#64758 SparseDtype allows a fill_value whose type is wider than the
        #  subtype as long as it round-trips; the result still follows the subtype
        keys = ["a", "a", "b", "b"]
        dense_ser = pd.Series(np.array([1, 0, 3, 0], dtype="int64"))
        sparse_ser = dense_ser.astype(pd.SparseDtype("int64", 0.0))
        assert sparse_ser.array.sp_index.ngaps == 2

        result = getattr(sparse_ser.groupby(keys), op)()
        expected = getattr(dense_ser.groupby(keys), op)()
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("nat_fill", [None, pd.NaT])
    @pytest.mark.parametrize(
        "unit",
        ["datetime64[ns]", "datetime64[s]", "timedelta64[ns]", "timedelta64[us]"],
    )
    @pytest.mark.parametrize("op", ["min", "max", "mean", "median", "cummax", "rank"])
    def test_sparse_groupby_datetimelike(self, unit, op, nat_fill):
        # GH#64758 datetimelike subtypes must not be treated as raw integers,
        #  and the result must keep the subtype's unit
        keys = ["a", "a", "b", "b"]
        dense_ser = pd.Series(np.array([1, 2, 3, "NaT"], dtype=unit))
        sparse_ser = dense_ser.astype(pd.SparseDtype(unit, nat_fill))
        assert sparse_ser.array.sp_index.ngaps == 1

        result = getattr(sparse_ser.groupby(keys), op)()
        expected = getattr(dense_ser.groupby(keys), op)()
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("op", ["min", "max", "mean", "cummax"])
    def test_sparse_groupby_datetimelike_subsecond_fill_value(self, op):
        # GH#64758 a sub-microsecond Timestamp fill_value must not be truncated
        #  when the gaps are materialized
        keys = ["a", "a", "b", "b"]
        dense_ser = pd.Series(np.array([1, 2, 3, 4], dtype="M8[ns]"))
        sparse_ser = dense_ser.astype(pd.SparseDtype("M8[ns]", pd.Timestamp(3)))
        assert sparse_ser.array.sp_index.ngaps == 1

        result = getattr(sparse_ser.groupby(keys), op)()
        expected = getattr(dense_ser.groupby(keys), op)()
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("op", ["min", "sum", "cumsum"])
    def test_sparse_groupby_timedelta_subsecond_fill_value(self, op):
        # GH#64758 same for a sub-microsecond Timedelta fill_value
        keys = ["a", "a", "b", "b"]
        dense_ser = pd.Series(np.array([1, 2, 3, 4], dtype="m8[ns]"))
        sparse_ser = dense_ser.astype(pd.SparseDtype("m8[ns]", pd.Timedelta(3)))
        assert sparse_ser.array.sp_index.ngaps == 1

        result = getattr(sparse_ser.groupby(keys), op)()
        expected = getattr(dense_ser.groupby(keys), op)()
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("nat_fill", [None, pd.NaT])
    @pytest.mark.parametrize("op", ["sum", "cumsum"])
    def test_sparse_groupby_timedelta_sum(self, op, nat_fill):
        # GH#64758 adding timedeltas is valid, unlike adding datetimes
        keys = ["a", "a", "b", "b"]
        dense_ser = pd.Series(np.array([1, 2, 3, "NaT"], dtype="m8[ns]"))
        sparse_ser = dense_ser.astype(pd.SparseDtype("m8[ns]", nat_fill))

        result = getattr(sparse_ser.groupby(keys), op)()
        expected = getattr(dense_ser.groupby(keys), op)()
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("nat_fill", [None, pd.NaT])
    @pytest.mark.parametrize(
        "unit, op",
        [
            ("datetime64[ns]", "sum"),
            ("datetime64[ns]", "prod"),
            ("datetime64[ns]", "var"),
            ("datetime64[ns]", "skew"),
            ("datetime64[ns]", "kurt"),
            ("datetime64[ns]", "cumsum"),
            ("datetime64[ns]", "cumprod"),
            # any/all raise with a different message than the arithmetic ops
            ("datetime64[ns]", "any"),
            ("datetime64[ns]", "all"),
            ("timedelta64[ns]", "prod"),
            ("timedelta64[ns]", "var"),
            ("timedelta64[ns]", "skew"),
            ("timedelta64[ns]", "kurt"),
            ("timedelta64[ns]", "cumprod"),
        ],
    )
    def test_sparse_groupby_datetimelike_invalid(self, unit, op, nat_fill):
        # GH#64758 operations that are invalid for the subtype must raise
        #  instead of silently returning garbage
        keys = ["a", "a", "b", "b"]
        dense_ser = pd.Series(np.array([1, 2, 3, "NaT"], dtype=unit))
        sparse_ser = dense_ser.astype(pd.SparseDtype(unit, nat_fill))

        with pytest.raises(TypeError) as sparse_err:
            getattr(sparse_ser.groupby(keys), op)()
        with pytest.raises(TypeError) as dense_err:
            getattr(dense_ser.groupby(keys), op)()
        assert str(sparse_err.value) == str(dense_err.value)
