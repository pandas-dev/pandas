import numpy as np
import pytest

from pandas import Categorical, DataFrame, Index, Series, Timestamp, date_range
import pandas._testing as tm
from pandas.core.arrays import SparseArray


class TestDataFrameSetItem:
    def test_setitem_error_msmgs(self):

        # GH 7432
        df = DataFrame(
            {"bar": [1, 2, 3], "baz": ["d", "e", "f"]},
            index=Index(["a", "b", "c"], name="foo"),
        )
        ser = Series(
            ["g", "h", "i", "j"],
            index=Index(["a", "b", "c", "a"], name="foo"),
            name="fiz",
        )
        msg = "cannot reindex from a duplicate axis"
        with pytest.raises(ValueError, match=msg):
            df["newcol"] = ser

        # GH 4107, more descriptive error message
        df = DataFrame(np.random.randint(0, 2, (4, 4)), columns=["a", "b", "c", "d"])

        msg = "incompatible index of inserted column with frame index"
        with pytest.raises(TypeError, match=msg):
            df["gr"] = df.groupby(["b", "c"]).count()

    def test_setitem_benchmark(self):
        # from the vb_suite/frame_methods/frame_insert_columns
        N = 10
        K = 5
        df = DataFrame(index=range(N))
        new_col = np.random.randn(N)
        for i in range(K):
            df[i] = new_col
        expected = DataFrame(np.repeat(new_col, K).reshape(N, K), index=range(N))
        tm.assert_frame_equal(df, expected)

    def test_setitem_different_dtype(self):
        df = DataFrame(
            np.random.randn(5, 3), index=np.arange(5), columns=["c", "b", "a"]
        )
        df.insert(0, "foo", df["a"])
        df.insert(2, "bar", df["c"])

        # diff dtype

        # new item
        df["x"] = df["a"].astype("float32")
        result = df.dtypes
        expected = Series(
            [np.dtype("float64")] * 5 + [np.dtype("float32")],
            index=["foo", "c", "bar", "b", "a", "x"],
        )
        tm.assert_series_equal(result, expected)

        # replacing current (in different block)
        df["a"] = df["a"].astype("float32")
        result = df.dtypes
        expected = Series(
            [np.dtype("float64")] * 4 + [np.dtype("float32")] * 2,
            index=["foo", "c", "bar", "b", "a", "x"],
        )
        tm.assert_series_equal(result, expected)

        df["y"] = df["a"].astype("int32")
        result = df.dtypes
        expected = Series(
            [np.dtype("float64")] * 4 + [np.dtype("float32")] * 2 + [np.dtype("int32")],
            index=["foo", "c", "bar", "b", "a", "x", "y"],
        )
        tm.assert_series_equal(result, expected)

    def test_setitem_empty_columns(self):
        # GH 13522
        df = DataFrame(index=["A", "B", "C"])
        df["X"] = df.index
        df["X"] = ["x", "y", "z"]
        exp = DataFrame(data={"X": ["x", "y", "z"]}, index=["A", "B", "C"])
        tm.assert_frame_equal(df, exp)

    def test_setitem_dt64_index_empty_columns(self):
        rng = date_range("1/1/2000 00:00:00", "1/1/2000 1:59:50", freq="10s")
        df = DataFrame(index=np.arange(len(rng)))

        df["A"] = rng
        assert df["A"].dtype == np.dtype("M8[ns]")

    def test_setitem_timestamp_empty_columns(self):
        # GH#19843
        df = DataFrame(index=range(3))
        df["now"] = Timestamp("20130101", tz="UTC")

        expected = DataFrame(
            [[Timestamp("20130101", tz="UTC")]] * 3, index=[0, 1, 2], columns=["now"],
        )
        tm.assert_frame_equal(df, expected)

    def test_setitem_wrong_length_categorical_dtype_raises(self):
        # GH#29523
        cat = Categorical.from_codes([0, 1, 1, 0, 1, 2], ["a", "b", "c"])
        df = DataFrame(range(10), columns=["bar"])

        msg = "Length of values does not match length of index"
        with pytest.raises(ValueError, match=msg):
            df["foo"] = cat

    def test_setitem_with_sparse_value(self):
        # GH#8131
        df = DataFrame({"c_1": ["a", "b", "c"], "n_1": [1.0, 2.0, 3.0]})
        sp_array = SparseArray([0, 0, 1])
        df["new_column"] = sp_array

        expected = Series(sp_array, name="new_column")
        tm.assert_series_equal(df["new_column"], expected)

    def test_setitem_with_unaligned_sparse_value(self):
        df = DataFrame({"c_1": ["a", "b", "c"], "n_1": [1.0, 2.0, 3.0]})
        sp_series = Series(SparseArray([0, 0, 1]), index=[2, 1, 0])

        df["new_column"] = sp_series
        expected = Series(SparseArray([1, 0, 0]), name="new_column")
        tm.assert_series_equal(df["new_column"], expected)

    def test_setitem_dict_preserves_dtypes(self):
        # https://github.com/pandas-dev/pandas/issues/34573
        expected = DataFrame(
            {
                "a": Series([0, 1, 2], dtype="int64"),
                "b": Series([1, 2, 3], dtype=float),
                "c": Series([1, 2, 3], dtype=float),
            }
        )
        df = DataFrame(
            {
                "a": Series([], dtype="int64"),
                "b": Series([], dtype=float),
                "c": Series([], dtype=float),
            }
        )
        for idx, b in enumerate([1, 2, 3]):
            df.loc[df.shape[0]] = {
                "a": int(idx),
                "b": float(b),
                "c": float(b),
            }
        tm.assert_frame_equal(df, expected)
