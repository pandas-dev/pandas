"""
Tests for Series.str.splitlines and Series.str.expandtabs.

These methods mirror Python's str.splitlines() and str.expandtabs().
"""
import numpy as np
import pytest

import pandas as pd
from pandas import Series, _testing as tm
from pandas.tests.strings import _convert_na_value


# --- splitlines tests ---


@pytest.mark.parametrize("keepends", [False, True])
def test_str_splitlines_basic(any_string_dtype, keepends):
    values = Series(
        ["line1\nline2", "line1\r\nline2", "no newline", np.nan],
        dtype=any_string_dtype,
    )
    result = values.str.splitlines(keepends=keepends)
    if keepends:
        exp = Series(
            [
                ["line1\n", "line2"],
                ["line1\r\n", "line2"],
                ["no newline"],
                np.nan,
            ]
        )
    else:
        exp = Series(
            [
                ["line1", "line2"],
                ["line1", "line2"],
                ["no newline"],
                np.nan,
            ]
        )
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


def test_str_splitlines_empty_string(any_string_dtype):
    values = Series(["", "a\nb", np.nan], dtype=any_string_dtype)
    result = values.str.splitlines()
    exp = Series([[], ["a", "b"], np.nan])
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


def test_str_splitlines_keepends_true(any_string_dtype):
    values = Series(["a\nb\nc"], dtype=any_string_dtype)
    result = values.str.splitlines(keepends=True)
    exp = Series([["a\n", "b\n", "c"]])
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


def test_str_splitlines_universal_newlines(any_string_dtype):
    # Test \n, \r, \r\n, \v, \f
    values = Series(
        ["a\nb", "a\rb", "a\r\nb", "a\vb", "a\fb"],
        dtype=any_string_dtype,
    )
    result = values.str.splitlines()
    exp = Series([["a", "b"], ["a", "b"], ["a", "b"], ["a", "b"], ["a", "b"]])
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


def test_str_splitlines_index(any_string_dtype):
    idx = pd.Index(["line1\nline2", "single"], dtype=any_string_dtype)
    result = idx.str.splitlines()
    exp = pd.Index([["line1", "line2"], ["single"]], dtype=object)
    tm.assert_index_equal(result, exp)


def test_str_splitlines_all_na(any_string_dtype):
    values = Series([np.nan, np.nan], dtype=any_string_dtype)
    result = values.str.splitlines()
    exp = Series([np.nan, np.nan], dtype=object)
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


def test_str_splitlines_empty_series(any_string_dtype):
    values = Series([], dtype=any_string_dtype)
    result = values.str.splitlines()
    exp = Series([], dtype=object)
    tm.assert_series_equal(result, exp)


# --- expandtabs tests ---


def test_str_expandtabs_tabsize_four(any_string_dtype):
    values = Series(["1\t2", "a\t\tb", np.nan], dtype=any_string_dtype)
    result = values.str.expandtabs(tabsize=4)
    exp = Series(["1   2", "a       b", np.nan])
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


def test_str_expandtabs_tabsize_eight(any_string_dtype):
    values = Series(["1\t2", "a\tb", np.nan], dtype=any_string_dtype)
    result = values.str.expandtabs(tabsize=8)
    exp = Series(["1       2", "a       b", np.nan])
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


def test_str_expandtabs_default(any_string_dtype):
    values = Series(["1\t2", "a\tb", np.nan], dtype=any_string_dtype)
    result = values.str.expandtabs()
    exp = Series(["1       2", "a       b", np.nan])
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


def test_str_expandtabs_no_tabs(any_string_dtype):
    values = Series(["no tabs", "here", np.nan], dtype=any_string_dtype)
    result = values.str.expandtabs(8)
    exp = Series(["no tabs", "here", np.nan])
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


def test_str_expandtabs_empty_string(any_string_dtype):
    values = Series(["", "\t", "a\tb"], dtype=any_string_dtype)
    result = values.str.expandtabs(4)
    exp = Series(["", "    ", "a   b"])
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


def test_str_expandtabs_index(any_string_dtype):
    idx = pd.Index(["1\t2", "a\tb"], dtype=any_string_dtype)
    result = idx.str.expandtabs(4)
    exp = pd.Index(["1   2", "a   b"])
    tm.assert_index_equal(result, exp)


def test_str_expandtabs_tabsize_zero(any_string_dtype):
    values = Series(["1\t2"], dtype=any_string_dtype)
    result = values.str.expandtabs(0)
    exp = Series(["12"])
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


def test_str_expandtabs_invalid_tabsize():
    values = Series(["a\tb"])
    with pytest.raises(TypeError, match="tabsize must be of integer type"):
        values.str.expandtabs(tabsize="8")
    with pytest.raises(ValueError, match="tabsize must be >= 0"):
        values.str.expandtabs(tabsize=-1)


def test_str_expandtabs_all_na(any_string_dtype):
    values = Series([np.nan, np.nan], dtype=any_string_dtype)
    result = values.str.expandtabs(8)
    exp = Series([np.nan, np.nan], dtype=any_string_dtype)
    exp = _convert_na_value(values, exp)
    tm.assert_series_equal(result, exp)


def test_str_expandtabs_empty_series(any_string_dtype):
    values = Series([], dtype=any_string_dtype)
    result = values.str.expandtabs(8)
    exp = Series([], dtype=any_string_dtype)
    tm.assert_series_equal(result, exp)


# --- equivalence with Python str ---


def test_str_splitlines_matches_str_splitlines():
    test_cases = ["a\nb\nc", "a\r\nb", "", "no newline", "one\ntwo\nthree"]
    for s in test_cases:
        series = Series([s])
        result = series.str.splitlines(keepends=False).iloc[0]
        expected = s.splitlines(keepends=False)
        assert result == expected
    for s in test_cases:
        series = Series([s])
        result = series.str.splitlines(keepends=True).iloc[0]
        expected = s.splitlines(keepends=True)
        assert result == expected


def test_str_expandtabs_matches_str_expandtabs():
    test_cases = ["1\t2", "a\t\tb", "\t", "", " \t "]
    for s in test_cases:
        for tabsize in [0, 4, 8]:
            series = Series([s])
            result = series.str.expandtabs(tabsize=tabsize).iloc[0]
            expected = s.expandtabs(tabsize=tabsize)
            assert result == expected
