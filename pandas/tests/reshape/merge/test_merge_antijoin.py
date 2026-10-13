import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm
from pandas.core.reshape.merge import merge


def test_merge_antijoin():
    # GH#42916
    left = pd.DataFrame({"A": [1, 2, 3]}, index=["a", "b", "c"])
    right = pd.DataFrame({"B": [1, 2, 4]}, index=["a", "b", "d"])

    result = merge(left, right, how="left_anti", left_index=True, right_index=True)
    expected = pd.DataFrame({"A": [3], "B": [np.nan]}, index=["c"])
    tm.assert_frame_equal(result, expected)

    result = merge(left, right, how="right_anti", left_index=True, right_index=True)
    expected = pd.DataFrame({"A": [np.nan], "B": [4]}, index=["d"])
    tm.assert_frame_equal(result, expected)


def test_merge_antijoin_on_different_columns():
    left = pd.DataFrame({"A": [1.0, 2.0, 3.0], "B": ["a", "b", "c"]}).astype(
        {"B": object}
    )
    right = pd.DataFrame({"C": [1.0, 2.0, 4.0], "D": ["a", "d", "b"]}).astype(
        {"D": object}
    )

    result = merge(left, right, how="left_anti", left_on="B", right_on="D")
    expected = pd.DataFrame(
        {
            "A": [3.0],
            "B": ["c"],
            "C": [np.nan],
            "D": [np.nan],
        },
        index=[2],
    ).astype({"B": object, "D": object})
    tm.assert_frame_equal(result, expected)

    result = merge(left, right, how="right_anti", left_on="B", right_on="D")
    expected = pd.DataFrame(
        {
            "A": [np.nan],
            "B": [np.nan],
            "C": [2.0],
            "D": ["d"],
        },
        index=[1],
    ).astype({"B": object, "D": object})
    tm.assert_frame_equal(result, expected)


def test_merge_antijoin_nonunique_keys():
    left = pd.DataFrame({"A": [1.0, 2.0, 3.0], "B": ["a", "b", "b"]}).astype(
        {"B": object}
    )
    right = pd.DataFrame({"C": [1.0, 2.0, 4.0], "D": ["b", "d", "d"]}).astype(
        {"D": object}
    )

    result = merge(left, right, how="left_anti", left_on="B", right_on="D")
    expected = pd.DataFrame(
        {
            "A": [1.0],
            "B": ["a"],
            "C": [np.nan],
            "D": [np.nan],
        },
        index=[0],
    ).astype({"B": object, "D": object})
    tm.assert_frame_equal(result, expected)

    result = merge(left, right, how="right_anti", left_on="B", right_on="D")
    expected = pd.DataFrame(
        {
            "A": [np.nan, np.nan],
            "B": [np.nan, np.nan],
            "C": [2.0, 4.0],
            "D": ["d", "d"],
        },
        index=[2, 3],
    ).astype({"B": object, "D": object})
    tm.assert_frame_equal(result, expected)


def test_merge_antijoin_same_df():
    left = pd.DataFrame({"A": [1, 2, 3]}, index=["a", "b", "c"], dtype=np.int64)
    result = merge(left, left, how="left_anti", left_index=True, right_index=True)
    expected = pd.DataFrame([], columns=["A_x", "A_y"], dtype=np.int64)
    tm.assert_frame_equal(result, expected, check_index_type=False)


def test_merge_antijoin_nans():
    left = pd.DataFrame({"A": [1.0, 2.0, np.nan], "C": ["a", "b", "c"]}).astype(
        {"C": object}
    )
    right = pd.DataFrame({"A": [3.0, 2.0, np.nan], "D": ["d", "e", "f"]}).astype(
        {"D": object}
    )
    result = merge(left, right, how="left_anti", on="A")
    expected = pd.DataFrame({"A": [1.0], "C": ["a"], "D": [np.nan]}).astype(
        {"C": object, "D": object}
    )
    tm.assert_frame_equal(result, expected)


def test_merge_antijoin_on_datetime64tz():
    # GH11405
    left = pd.DataFrame(
        {
            "key": pd.date_range("20151010", periods=2, tz="US/Eastern"),
            "value": [1.0, 2.0],
        }
    )
    right = pd.DataFrame(
        {
            "key": pd.date_range("20151011", periods=3, tz="US/Eastern"),
            "value": [1.0, 2.0, 3.0],
        }
    )

    expected = pd.DataFrame(
        {
            "key": pd.date_range("20151010", periods=1, tz="US/Eastern"),
            "value_x": [1.0],
            "value_y": [np.nan],
        },
        index=[0],
    )
    result = merge(left, right, on="key", how="left_anti")
    tm.assert_frame_equal(result, expected)

    expected = pd.DataFrame(
        {
            "key": pd.date_range("20151012", periods=2, tz="US/Eastern"),
            "value_x": [np.nan, np.nan],
            "value_y": [2.0, 3.0],
        },
        index=[1, 2],
    )
    result = merge(left, right, on="key", how="right_anti")
    tm.assert_frame_equal(result, expected)


def test_merge_antijoin_multiindex():
    left = pd.DataFrame(
        {
            "A": [1, 2, 3],
            "B": [4, 5, 6],
        },
        index=pd.MultiIndex.from_tuples(
            [("a", "x"), ("b", "y"), ("c", "z")], names=["first", "second"]
        ),
    )
    right = pd.DataFrame(
        {
            "C": [7, 8, 9],
            "D": [10, 11, 12],
        },
        index=pd.MultiIndex.from_tuples(
            [("a", "x"), ("b", "y"), ("c", "w")], names=["first", "second"]
        ),
    )

    result = merge(left, right, how="left_anti", left_index=True, right_index=True)
    expected = pd.DataFrame(
        {
            "A": [3],
            "B": [6],
            "C": [np.nan],
            "D": [np.nan],
        },
        index=pd.MultiIndex.from_tuples([("c", "z")], names=["first", "second"]),
    )
    tm.assert_frame_equal(result, expected)

    result = merge(left, right, how="right_anti", left_index=True, right_index=True)
    expected = pd.DataFrame(
        {
            "A": [np.nan],
            "B": [np.nan],
            "C": [9],
            "D": [12],
        },
        index=pd.MultiIndex.from_tuples([("c", "w")], names=["first", "second"]),
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "dtype",
    [
        "Int64",
        pytest.param("int64[pyarrow]", marks=td.skip_if_no("pyarrow")),
        pytest.param("timestamp[s][pyarrow]", marks=td.skip_if_no("pyarrow")),
        pytest.param("string[pyarrow]", marks=td.skip_if_no("pyarrow")),
    ],
)
def test_merge_antijoin_extension_dtype(dtype):
    left = pd.DataFrame(
        {
            "join_col": [1, 3, 5],
            "left_val": [1, 2, 3],
        }
    )
    right = pd.DataFrame(
        {
            "join_col": [2, 3, 4],
            "right_val": [1, 2, 3],
        }
    )
    left = left.astype({"join_col": dtype})
    right = right.astype({"join_col": dtype})
    result = merge(left, right, how="left_anti", on="join_col")
    expected = pd.DataFrame(
        {
            "join_col": [1, 5],
            "left_val": [1, 3],
            "right_val": [np.nan, np.nan],
        },
        index=[0, 2],
    )
    expected = expected.astype({"join_col": dtype})
    tm.assert_frame_equal(result, expected)


def test_merge_antijoin_empty_dataframe():
    left = pd.DataFrame({"A": [], "B": []})
    right = pd.DataFrame({"C": [], "D": []})

    result = merge(left, right, how="left_anti", left_on="A", right_on="C")
    expected = pd.DataFrame({"A": [], "B": [], "C": [], "D": []})
    tm.assert_frame_equal(result, expected)

    result = merge(left, right, how="right_anti", left_on="A", right_on="C")
    tm.assert_frame_equal(result, expected)


def test_merge_antijoin_no_common_elements():
    left = pd.DataFrame({"A": [1, 2, 3]})
    right = pd.DataFrame({"B": [4, 5, 6]})

    result = merge(left, right, how="left_anti", left_on="A", right_on="B")
    expected = pd.DataFrame({"A": [1, 2, 3], "B": [np.nan, np.nan, np.nan]})
    tm.assert_frame_equal(result, expected)

    result = merge(left, right, how="right_anti", left_on="A", right_on="B")
    expected = pd.DataFrame({"A": [np.nan, np.nan, np.nan], "B": [4, 5, 6]})
    tm.assert_frame_equal(result, expected)


def test_merge_antijoin_with_null_values():
    left = pd.DataFrame({"A": [1.0, 2.0, None, 4.0]})
    right = pd.DataFrame({"B": [2.0, None, 5.0]})

    result = merge(left, right, how="left_anti", left_on="A", right_on="B")
    expected = pd.DataFrame({"A": [1.0, 4.0], "B": [np.nan, np.nan]}, index=[0, 3])
    tm.assert_frame_equal(result, expected)

    result = merge(left, right, how="right_anti", left_on="A", right_on="B")
    expected = pd.DataFrame({"A": [np.nan], "B": [5.0]}, index=[2])
    tm.assert_frame_equal(result, expected)


def test_merge_antijoin_with_mixed_dtypes():
    left = pd.DataFrame({"A": [1, "2", 3.0]})
    right = pd.DataFrame({"B": ["2", 3.0, 4]})

    result = merge(left, right, how="left_anti", left_on="A", right_on="B")
    expected = pd.DataFrame({"A": [1], "B": [np.nan]}, dtype=object)
    tm.assert_frame_equal(result, expected)

    result = merge(left, right, how="right_anti", left_on="A", right_on="B")
    expected = pd.DataFrame({"A": [np.nan], "B": [4]}, dtype=object, index=[2])
    tm.assert_frame_equal(result, expected)
