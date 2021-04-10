import builtins

import numpy as np
import pytest

from pandas import (
    DataFrame,
    Index,
    isna,
)
import pandas._testing as tm


@pytest.mark.parametrize("agg_func", ["any", "all"])
@pytest.mark.parametrize("skipna", [True, False])
@pytest.mark.parametrize(
    "vals",
    [
        ["foo", "bar", "baz"],
        ["foo", "", ""],
        ["", "", ""],
        [1, 2, 3],
        [1, 0, 0],
        [0, 0, 0],
        [1.0, 2.0, 3.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [True, True, True],
        [True, False, False],
        [False, False, False],
        [np.nan, np.nan, np.nan],
    ],
)
def test_groupby_bool_aggs(agg_func, skipna, vals):
    df = DataFrame({"key": ["a"] * 3 + ["b"] * 3, "val": vals * 2})

    # Figure out expectation using Python builtin
    exp = getattr(builtins, agg_func)(vals)

    # edge case for missing data with skipna and 'any'
    if skipna and all(isna(vals)) and agg_func == "any":
        exp = False

    exp_df = DataFrame([exp] * 2, columns=["val"], index=Index(["a", "b"], name="key"))
    result = getattr(df.groupby("key"), agg_func)(skipna=skipna)
    tm.assert_frame_equal(result, exp_df)


def test_any():
    df = DataFrame(
        [[1, 2, "foo"], [1, np.nan, "bar"], [3, np.nan, "baz"]],
        columns=["A", "B", "C"],
    )
    expected = DataFrame(
        [[True, True], [False, True]], columns=["B", "C"], index=[1, 3]
    )
    expected.index.name = "A"
    result = df.groupby("A").any()
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("bool_agg_func", ["any", "all"])
def test_bool_aggs_dup_column_labels(bool_agg_func):
    # 21668
    df = DataFrame([[True, True]], columns=["a", "a"])
    grp_by = df.groupby([0])
    result = getattr(grp_by, bool_agg_func)()

    expected = df
    tm.assert_frame_equal(result, expected)
