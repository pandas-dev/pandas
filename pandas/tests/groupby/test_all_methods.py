"""
Tests that apply to all groupby operation methods.

The only tests that should appear here are those that use the `groupby_func` fixture.
Even if it does use that fixture, prefer a more specific test file if it available
such as:

 - test_categorical
 - test_groupby_dropna
 - test_groupby_subclass
 - test_raises
"""

import pytest

from pandas import DataFrame
import pandas._testing as tm
from pandas.tests.groupby import get_groupby_method_args


def test_multiindex_group_all_columns_when_empty(groupby_func):
    # GH 32464
    df = DataFrame({"a": [], "b": [], "c": []}).set_index(["a", "b", "c"])
    gb = df.groupby(["a", "b", "c"], group_keys=False)
    method = getattr(gb, groupby_func)
    args = get_groupby_method_args(groupby_func, df)

    warn = FutureWarning if groupby_func == "fillna" else None
    warn_msg = "DataFrameGroupBy.fillna is deprecated"
    with tm.assert_produces_warning(warn, match=warn_msg):
        result = method(*args).index
    expected = df.index
    tm.assert_index_equal(result, expected)


def test_duplicate_columns(request, groupby_func, as_index):
    # GH#50806
    if groupby_func == "corrwith":
        msg = "GH#50845 - corrwith fails when there are duplicate columns"
        request.applymarker(pytest.mark.xfail(reason=msg))
    df = DataFrame([[1, 3, 6], [1, 4, 7], [2, 5, 8]], columns=list("abb"))
    args = get_groupby_method_args(groupby_func, df)
    gb = df.groupby("a", as_index=as_index)
    warn = FutureWarning if groupby_func == "fillna" else None
    warn_msg = "DataFrameGroupBy.fillna is deprecated"
    with tm.assert_produces_warning(warn, match=warn_msg):
        result = getattr(gb, groupby_func)(*args)

    expected_df = df.set_axis(["a", "b", "c"], axis=1)
    expected_args = get_groupby_method_args(groupby_func, expected_df)
    expected_gb = expected_df.groupby("a", as_index=as_index)
    warn = FutureWarning if groupby_func == "fillna" else None
    warn_msg = "DataFrameGroupBy.fillna is deprecated"
    with tm.assert_produces_warning(warn, match=warn_msg):
        expected = getattr(expected_gb, groupby_func)(*expected_args)
    if groupby_func not in ("size", "ngroup", "cumcount"):
        expected = expected.rename(columns={"c": "b"})
    tm.assert_equal(result, expected)
