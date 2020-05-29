import pytest

import pandas as pd
import pandas._testing as tm


class TestFormat:
    @pytest.mark.parametrize("format_str", ["{}-{}", "{A}-{B}", "{}-{B}"])
    @pytest.mark.parametrize("name", [None, "X"])
    @pytest.mark.parametrize("how_na", ["all", "any"])
    def test_basic(self, format_str, name, how_na):
        df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        expected = pd.Series(["1-4", "2-5", "3-6"], dtype="string", name=name)

        result = df.format(format_str, name=name, how_na=how_na)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("format_str", ["{Index}-{}-{}", "{Index}-{A}-{B}"])
    def test_with_index(self, format_str):
        df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        expected = pd.Series(["0-1-4", "1-2-5", "2-3-6"], dtype="string")

        result = df.format(format_str)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("format_str", ["{}-{}"])
    @pytest.mark.parametrize("positional_only", [True, False])
    def test_positional_only(self, format_str, positional_only):
        df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        expected = pd.Series(["1-4", "2-5", "3-6"], dtype="string")

        result = df.format(format_str, positional_only=positional_only)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("format_str", ["{A}-{B}", "{A}-{}", "{Index}-{}"])
    def test_positional_only_raises(self, format_str):
        df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        with pytest.raises(KeyError):
            df.format(format_str, positional_only=True)

    @pytest.mark.parametrize(
        "how_na, expected",
        [("any", ["1-4", pd.NA, pd.NA]), ("all", ["1-4", "nan-5", pd.NA])],
    )
    def test_na_how(self, how_na, expected):
        df = pd.DataFrame({"A": [1, None, None], "B": [4, 5, None]})
        expected = pd.Series(expected, dtype="string")

        result = df.format("{:.0f}-{:.0f}", how_na=how_na)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("format_string", ["{}-{}-{}", "{0}-{1}-{2}"])
    def test_too_many_positional_args(self, format_string):
        df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        with pytest.raises(IndexError):
            df.format(format_string)

    @pytest.mark.parametrize("format_string", ["{A}-{B}-{C}", "{C}"])
    def test_too_many_named_args(self, format_string):
        df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        with pytest.raises(KeyError):
            df.format(format_string)
