import pytest

import pandas as pd
from pandas import DataFrame
import pandas._testing as tm


@pytest.fixture
def regular_df():
    return DataFrame({"a": [1, 2], "b": [3, 4], "c": [5, 6], "d": [7, 8]})


@pytest.fixture
def multiindex_df():
    return DataFrame(
        [(0, 2, 4), (1, 3, 5)],
        columns=pd.MultiIndex.from_tuples([("A", "c"), ("A", "d"), ("B", "e")]),
    )


class TestSelect:
    def test_select_subset_cols(self, regular_df):
        expected = DataFrame({"a": [1, 2], "c": [5, 6]})
        result = regular_df.select("a", "c")
        tm.assert_frame_equal(result, expected)

    def test_single_value(self, regular_df):
        expected = DataFrame({"a": [1, 2]})
        result = regular_df.select("a")
        assert isinstance(result, DataFrame)
        tm.assert_frame_equal(result, expected)

    def test_select_change_order(self, regular_df):
        expected = DataFrame({"b": [3, 4], "d": [7, 8], "a": [1, 2], "c": [5, 6]})
        result = regular_df.select("b", "d", "a", "c")
        tm.assert_frame_equal(result, expected)

    def test_select_none(self, regular_df):
        result = regular_df.select()
        assert result.empty

    def test_select_duplicated(self, regular_df):
        expected = ["a", "d", "a"]
        result = regular_df.select("a", "d", "a")
        assert result.columns.tolist() == expected

    def test_select_list(self, regular_df):
        with pytest.raises(ValueError, match="does not support a list"):
            regular_df.select(["a", "b"])

    def test_select_missing(self, regular_df):
        with pytest.raises(KeyError, match=r"None of .* are in the \[columns\]"):
            regular_df.select("z")

    def test_select_not_hashable(self, regular_df):
        with pytest.raises(TypeError, match="unhashable type"):
            regular_df.select(set())

    def test_select_multiindex_one_level(self, multiindex_df):
        expected = DataFrame(
            [(0, 2), (1, 3)],
            columns=pd.MultiIndex.from_tuples([("A", "c"), ("A", "d")]),
        )
        result = multiindex_df.select("A")
        tm.assert_frame_equal(result, expected)

    def test_select_multiindex_single_column(self, multiindex_df):
        expected = DataFrame(
            [(2,), (3,)], columns=pd.MultiIndex.from_tuples([("A", "d")])
        )
        result = multiindex_df.select(("A", "d"))
        assert isinstance(result, DataFrame)
        tm.assert_frame_equal(result, expected)

    def test_select_multiindex_multiple_columns(self, multiindex_df):
        expected = DataFrame(
            [(0, 4), (1, 5)],
            columns=pd.MultiIndex.from_tuples([("A", "c"), ("B", "e")]),
        )
        result = multiindex_df.select(("A", "c"), ("B", "e"))
        tm.assert_frame_equal(result, expected)

    def test_select_multiindex_missing(self, multiindex_df):
        with pytest.raises(KeyError, match="not in index"):
            multiindex_df.select("Z")
