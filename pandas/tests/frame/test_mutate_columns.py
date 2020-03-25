import re

import numpy as np
import pytest

from pandas import DataFrame, MultiIndex, Series
import pandas._testing as tm

# Column add, remove, delete.


class TestDataFrameMutateColumns:
    def test_setitem_error_msmgs(self):

        # GH 7432
        df = DataFrame(
            {"foo": ["a", "b", "c"], "bar": [1, 2, 3], "baz": ["d", "e", "f"]}
        ).set_index("foo")
        s = DataFrame(
            {"foo": ["a", "b", "c", "a"], "fiz": ["g", "h", "i", "j"]}
        ).set_index("foo")
        msg = "cannot reindex from a duplicate axis"
        with pytest.raises(ValueError, match=msg):
            df["newcol"] = s

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

    def test_delitem(self, float_frame):
        del float_frame["A"]
        assert "A" not in float_frame

    def test_delitem_multiindex(self):
        midx = MultiIndex.from_product([["A", "B"], [1, 2]])
        df = DataFrame(np.random.randn(4, 4), columns=midx)
        assert len(df.columns) == 4
        assert ("A",) in df.columns
        assert "A" in df.columns

        result = df["A"]
        assert isinstance(result, DataFrame)
        del df["A"]

        assert len(df.columns) == 2

        # A still in the levels, BUT get a KeyError if trying
        # to delete
        assert ("A",) not in df.columns
        with pytest.raises(KeyError, match=re.escape("('A',)")):
            del df[("A",)]

        # behavior of dropped/deleted MultiIndex levels changed from
        # GH 2770 to GH 19027: MultiIndex no longer '.__contains__'
        # levels which are dropped/deleted
        assert "A" not in df.columns
        with pytest.raises(KeyError, match=re.escape("('A',)")):
            del df["A"]

    def test_pop(self, float_frame):
        float_frame.columns.name = "baz"

        float_frame.pop("A")
        assert "A" not in float_frame

        float_frame["foo"] = "bar"
        float_frame.pop("foo")
        assert "foo" not in float_frame
        assert float_frame.columns.name == "baz"

        # gh-10912: inplace ops cause caching issue
        a = DataFrame([[1, 2, 3], [4, 5, 6]], columns=["A", "B", "C"], index=["X", "Y"])
        b = a.pop("B")
        b += 1

        # original frame
        expected = DataFrame([[1, 3], [4, 6]], columns=["A", "C"], index=["X", "Y"])
        tm.assert_frame_equal(a, expected)

        # result
        expected = Series([2, 5], index=["X", "Y"], name="B") + 1
        tm.assert_series_equal(b, expected)

    def test_pop_non_unique_cols(self):
        df = DataFrame({0: [0, 1], 1: [0, 1], 2: [4, 5]})
        df.columns = ["a", "b", "a"]

        res = df.pop("a")
        assert type(res) == DataFrame
        assert len(res) == 2
        assert len(df.columns) == 1
        assert "b" in df.columns
        assert "a" not in df.columns
        assert len(df.index) == 2
