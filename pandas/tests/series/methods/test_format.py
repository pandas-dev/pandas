import pytest

import pandas as pd
import pandas._testing as tm


class TestFormat:
    @pytest.mark.parametrize("format_str", ["Value: {}", "Value: {A}"])
    @pytest.mark.parametrize("how_na", ["all", "any"])
    def test_basic(self, format_str, how_na):
        ser = pd.Series([1, 2, 3], name="A")
        expected = pd.Series(
            ["Value: 1", "Value: 2", "Value: 3"], dtype="string", name="X"
        )

        result = ser.format(format_str, how_na=how_na, name="X")
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("format_str", ["{Index}-{}", "{Index}-{A}"])
    def test_with_index(self, format_str):
        ser = pd.Series([1, 2, 3], name="A")
        expected = pd.Series(["0-1", "1-2", "2-3"], dtype="string", name="X")

        result = ser.format(format_str, name="X")
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("format_str", ["Value: {}"])
    @pytest.mark.parametrize("positional_only", [True, False])
    def test_positional_only(self, format_str, positional_only):
        ser = pd.Series([1, 2, 3], name="A")
        expected = pd.Series(["Value: 1", "Value: 2", "Value: 3"], dtype="string")

        result = ser.format(format_str, positional_only=positional_only)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("format_str", ["{A}-{}", "{Index}-{}"])
    def test_positional_only_raises(self, format_str):
        ser = pd.Series([1, 2, 3], name="A")
        with pytest.raises(KeyError):
            ser.format(format_str, positional_only=True)

    @pytest.mark.parametrize(
        "how_na, expected",
        [("any", ["Value: 1", pd.NA, pd.NA]), ("all", ["Value: 1", pd.NA, pd.NA])],
    )
    @pytest.mark.parametrize("format_str", ["Value: {}", "Value: {A}"])
    def test_na_how(self, how_na, expected, format_str):
        ser = pd.Series([1, pd.NA, pd.NA], name="A")
        expected = pd.Series(expected, dtype="string")

        result = ser.format("Value: {}", how_na=how_na)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("format_string", ["{}-{}", "{0}-{1}"])
    def test_too_many_positional_args(self, format_string):
        ser = pd.Series([1, 2, 3], name="A")
        with pytest.raises(IndexError):
            ser.format(format_string)

    @pytest.mark.parametrize("format_string", ["{A}-{B}", "{B}"])
    def test_unknown_named_args(self, format_string):
        ser = pd.Series([1, 2, 3], name="A")
        with pytest.raises(KeyError):
            ser.format(format_string)
