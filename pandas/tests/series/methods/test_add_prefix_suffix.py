import pytest

import pandas as pd
import pandas._testing as tm


def test_add_prefix_suffix(string_series):
    with_prefix = string_series.add_prefix("foo#")
    expected = pd.Index([f"foo#{c}" for c in string_series.index])
    tm.assert_index_equal(with_prefix.index, expected)

    with_suffix = string_series.add_suffix("#foo")
    expected = pd.Index([f"{c}#foo" for c in string_series.index])
    tm.assert_index_equal(with_suffix.index, expected)

    with_pct_prefix = string_series.add_prefix("%")
    expected = pd.Index([f"%{c}" for c in string_series.index])
    tm.assert_index_equal(with_pct_prefix.index, expected)

    with_pct_suffix = string_series.add_suffix("%")
    expected = pd.Index([f"{c}%" for c in string_series.index])
    tm.assert_index_equal(with_pct_suffix.index, expected)


def test_add_prefix_suffix_axis(string_series):
    # GH 47819
    with_prefix = string_series.add_prefix("foo#", axis=0)
    expected = pd.Index([f"foo#{c}" for c in string_series.index])
    tm.assert_index_equal(with_prefix.index, expected)

    with_pct_suffix = string_series.add_suffix("#foo", axis=0)
    expected = pd.Index([f"{c}#foo" for c in string_series.index])
    tm.assert_index_equal(with_pct_suffix.index, expected)


def test_add_prefix_suffix_invalid_axis(string_series):
    with pytest.raises(ValueError, match="No axis named 1 for object type Series"):
        string_series.add_prefix("foo#", axis=1)

    with pytest.raises(ValueError, match="No axis named 1 for object type Series"):
        string_series.add_suffix("foo#", axis=1)
