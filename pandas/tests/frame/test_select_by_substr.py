import pytest
from pandas import DataFrame
import pandas._testing as tm

class TestSelectBySubstr:
    @pytest.fixture
    def df(self):
        return DataFrame({
            "first_name": ["Alice", "Bob"],
            "last_name": ["Smith", "Jones"],
            "age": [25, 30],
            "City": ["NY", "LA"],
            "FOO_column": [1.1, 2.2],
            "bar_col": [True, False]
        })

    def test_basic_substring(self, df):
        result = df.select_by_substr("name")
        expected = df[["first_name", "last_name"]]
        tm.assert_frame_equal(result, expected)

    def test_multiple_substrings(self, df):
        result = df.select_by_substr(["name", "city"])
        expected = df[["first_name", "last_name", "City"]]
        tm.assert_frame_equal(result, expected)

    def test_case_sensitivity(self, df):
        # Case-sensitive search (no match)
        result = df.select_by_substr("city", ignore_case=False)
        expected = df[[]]
        tm.assert_frame_equal(result, expected)

        # Case-insensitive search (matches "City")
        result = df.select_by_substr("city")
        expected = df[["City"]]
        tm.assert_frame_equal(result, expected)

    def test_no_matches(self, df):
        result = df.select_by_substr("xyz")
        expected = df[[]]
        tm.assert_frame_equal(result, expected)

    def test_empty_dataframe(self):
        df = DataFrame()
        result = df.select_by_substr("test")
        expected = df[[]]
        tm.assert_frame_equal(result, expected)

    def test_duplicate_columns(self):
        df = DataFrame([[1, 2]], columns=["foo", "foo"])
        result = df.select_by_substr("foo")
        expected = df[["foo", "foo"]]
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("substr,expected_cols", [
        ("", []),  # Empty substring is ignored
        (["", "name"], ["first_name", "last_name"]),  # "" is ignored
        (["  ", "age"], ["age"]),  # "  " is ignored
    ])
    def test_edge_cases(self, df, substr, expected_cols):
        result = df.select_by_substr(substr)
        expected = df[expected_cols] if expected_cols else df[[]]
        tm.assert_frame_equal(result, expected)

    def test_mixed_case_columns(self, df):
        result = df.select_by_substr("foo")
        expected = df[["FOO_column"]]
        tm.assert_frame_equal(result, expected)

    def test_return_type(self, df):
        # With matches
        result = df.select_by_substr("name")
        assert isinstance(result, DataFrame)
        # Without matches
        result = df.select_by_substr("invalid")
        assert isinstance(result, DataFrame)
        assert result.empty

    def test_original_frame_unmodified(self, df):
        original_cols = df.columns.copy()
        df.select_by_substr("name")
        tm.assert_index_equal(df.columns, original_cols)
