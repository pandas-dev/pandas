"""
test_insert is specifically for the DataFrame.insert method; not to be
confused with tests with "insert" in their names that are really testing
__setitem__.
"""
import numpy as np
import pytest

from pandas import DataFrame, Index
import pandas._testing as tm


class TestDataFrameInsert:
    def test_insert(self):
        df = DataFrame(
            np.random.randn(5, 3), index=np.arange(5), columns=["c", "b", "a"]
        )

        df.insert(0, "foo", df["a"])
        tm.assert_index_equal(df.columns, Index(["foo", "c", "b", "a"]))
        tm.assert_series_equal(df["a"], df["foo"], check_names=False)

        df.insert(2, "bar", df["c"])
        tm.assert_index_equal(df.columns, Index(["foo", "c", "bar", "b", "a"]))
        tm.assert_almost_equal(df["c"], df["bar"], check_names=False)

        with pytest.raises(ValueError, match="already exists"):
            df.insert(1, "a", df["b"])

        msg = "cannot insert c, already exists"
        with pytest.raises(ValueError, match=msg):
            df.insert(1, "c", df["b"])

        df.columns.name = "some_name"
        # preserve columns name field
        df.insert(0, "baz", df["c"])
        assert df.columns.name == "some_name"

    def test_insert_column_bug_4032(self):

        # GH4032, inserting a column and renaming causing errors
        df = DataFrame({"b": [1.1, 2.2]})

        df = df.rename(columns={})
        df.insert(0, "a", [1, 2])
        result = df.rename(columns={})

        str(result)
        expected = DataFrame([[1, 1.1], [2, 2.2]], columns=["a", "b"])
        tm.assert_frame_equal(result, expected)

        df.insert(0, "c", [1.3, 2.3])
        result = df.rename(columns={})

        str(result)
        expected = DataFrame([[1.3, 1, 1.1], [2.3, 2, 2.2]], columns=["c", "a", "b"])
        tm.assert_frame_equal(result, expected)
