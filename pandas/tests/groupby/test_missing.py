import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    Index,
)
import pandas._testing as tm


@pytest.mark.parametrize("func", ["ffill", "bfill"])
def test_groupby_column_index_name_lost_fill_funcs(func):
    # GH: 29764 groupby loses index sometimes
    df = DataFrame(
        [[1, 1.0, -1.0], [1, np.nan, np.nan], [1, 2.0, -2.0]],
        columns=Index(["type", "a", "b"], name="idx"),
    )
    df_grouped = df.groupby(["type"])[["a", "b"]]
    result = getattr(df_grouped, func)().columns
    expected = Index(["a", "b"], name="idx")
    tm.assert_index_equal(result, expected)


@pytest.mark.parametrize("func", ["ffill", "bfill"])
def test_groupby_fill_duplicate_column_names(func):
    # GH: 25610 ValueError with duplicate column names
    df1 = DataFrame({"field1": [1, 3, 4], "field2": [1, 3, 4]})
    df2 = DataFrame({"field1": [1, np.nan, 4]})
    df_grouped = pd.concat([df1, df2], axis=1).groupby(by=["field2"])
    expected = DataFrame(
        [[1, 1.0], [3, np.nan], [4, 4.0]], columns=["field1", "field1"]
    )
    result = getattr(df_grouped, func)()
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("method", ["ffill", "bfill"])
@pytest.mark.parametrize("has_nan_group", [True, False])
def test_ffill_handles_nan_groups(dropna, method, has_nan_group):
    # GH 34725

    df_without_nan_rows = DataFrame([(1, 0.1), (2, 0.2)])

    ridx = [-1, 0, -1, -1, 1, -1]
    df = df_without_nan_rows.reindex(ridx).reset_index(drop=True)

    group_b = np.nan if has_nan_group else "b"
    df["group_col"] = pd.Series(["a"] * 3 + [group_b] * 3)

    grouped = df.groupby(by="group_col", dropna=dropna)
    result = getattr(grouped, method)(limit=None)

    expected_rows = {
        ("ffill", True, True): [-1, 0, 0, -1, -1, -1],
        ("ffill", True, False): [-1, 0, 0, -1, 1, 1],
        ("ffill", False, True): [-1, 0, 0, -1, 1, 1],
        ("ffill", False, False): [-1, 0, 0, -1, 1, 1],
        ("bfill", True, True): [0, 0, -1, -1, -1, -1],
        ("bfill", True, False): [0, 0, -1, 1, 1, -1],
        ("bfill", False, True): [0, 0, -1, 1, 1, -1],
        ("bfill", False, False): [0, 0, -1, 1, 1, -1],
    }

    ridx = expected_rows.get((method, dropna, has_nan_group))
    expected = df_without_nan_rows.reindex(ridx).reset_index(drop=True)
    # columns are a 'take' on df.columns, which are object dtype
    expected.columns = expected.columns.astype(object)

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("min_count, value", [(2, np.nan), (-1, 1.0)])
@pytest.mark.parametrize("func", ["first", "last", "max", "min"])
def test_min_count(func, min_count, value):
    # GH#37821
    df = DataFrame({"a": [1] * 3, "b": [1, np.nan, np.nan], "c": [np.nan] * 3})
    result = getattr(df.groupby("a"), func)(min_count=min_count)
    expected = DataFrame({"b": [value], "c": [np.nan]}, index=Index([1], name="a"))
    tm.assert_frame_equal(result, expected)


def test_indices_with_missing():
    # GH 9304
    df = DataFrame({"a": [1, 1, np.nan], "b": [2, 3, 4], "c": [5, 6, 7]})
    g = df.groupby(["a", "b"])
    result = g.indices
    expected = {(1.0, 2): np.array([0]), (1.0, 3): np.array([1])}
    assert result == expected
