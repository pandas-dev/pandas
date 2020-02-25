import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize("axis", [0, 1, "index", "columns"])
def test_differences_axis(axis):
    df = pd.DataFrame(
        {"col1": ["a", "b", "c"], "col2": [1.0, 2.0, np.nan], "col3": [1.0, 2.0, 3.0]},
        columns=["col1", "col2", "col3"],
    )
    df2 = df.copy()
    df2.loc[0, "col1"] = "c"
    df2.loc[2, "col3"] = 4.0

    result = df.differences(df2, axis=axis)

    if axis in (1, "columns"):
        indices = pd.Index([0, 2])
        columns = pd.MultiIndex.from_product([["col1", "col3"], ["self", "other"]])
        expected = pd.DataFrame(
            [["a", "c", np.nan, np.nan], [np.nan, np.nan, 3.0, 4.0]],
            index=indices,
            columns=columns,
        )
    else:
        indices = pd.MultiIndex.from_product([[0, 2], ["self", "other"]])
        columns = pd.Index(["col1", "col3"])
        expected = pd.DataFrame(
            [["a", np.nan], ["c", np.nan], [np.nan, 3.0], [np.nan, 4.0]],
            index=indices,
            columns=columns,
        )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "keep_shape, keep_equal",
    [
        (True, False),
        (False, True),
        (True, True),
        # False, False case is already covered in test_differences_axis
    ],
)
def test_differences_various_formats(keep_shape, keep_equal):
    df = pd.DataFrame(
        {"col1": ["a", "b", "c"], "col2": [1.0, 2.0, np.nan], "col3": [1.0, 2.0, 3.0]},
        columns=["col1", "col2", "col3"],
    )
    df2 = df.copy()
    df2.loc[0, "col1"] = "c"
    df2.loc[2, "col3"] = 4.0

    result = df.differences(df2, keep_shape=keep_shape, keep_equal=keep_equal)

    if keep_shape:
        indices = pd.Index([0, 1, 2])
        columns = pd.MultiIndex.from_product(
            [["col1", "col2", "col3"], ["self", "other"]]
        )
        if keep_equal:
            expected = pd.DataFrame(
                [
                    ["a", "c", 1.0, 1.0, 1.0, 1.0],
                    ["b", "b", 2.0, 2.0, 2.0, 2.0],
                    ["c", "c", np.nan, np.nan, 3.0, 4.0],
                ],
                index=indices,
                columns=columns,
            )
        else:
            expected = pd.DataFrame(
                [
                    ["a", "c", np.nan, np.nan, np.nan, np.nan],
                    [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                    [np.nan, np.nan, np.nan, np.nan, 3.0, 4.0],
                ],
                index=indices,
                columns=columns,
            )
    else:
        indices = pd.Index([0, 2])
        columns = pd.MultiIndex.from_product([["col1", "col3"], ["self", "other"]])
        expected = pd.DataFrame(
            [["a", "c", 1.0, 1.0], ["c", "c", 3.0, 4.0]], index=indices, columns=columns
        )
    tm.assert_frame_equal(result, expected)


def test_differences_with_equal_nulls():
    # We want to make sure two NaNs are considered the same
    # and dropped where applicable
    df = pd.DataFrame(
        {"col1": ["a", "b", "c"], "col2": [1.0, 2.0, np.nan], "col3": [1.0, 2.0, 3.0]},
        columns=["col1", "col2", "col3"],
    )
    df2 = df.copy()
    df2.loc[0, "col1"] = "c"

    result = df.differences(df2)
    indices = pd.Index([0])
    columns = pd.MultiIndex.from_product([["col1"], ["self", "other"]])
    expected = pd.DataFrame([["a", "c"]], index=indices, columns=columns)
    tm.assert_frame_equal(result, expected)


def test_differences_with_non_equal_nulls():
    # We want to make sure the relevant NaNs do not get dropped
    # even if the entire row or column are NaNs
    df = pd.DataFrame(
        {"col1": ["a", "b", "c"], "col2": [1.0, 2.0, np.nan], "col3": [1.0, 2.0, 3.0]},
        columns=["col1", "col2", "col3"],
    )
    df2 = df.copy()
    df2.loc[0, "col1"] = "c"
    df2.loc[2, "col3"] = np.nan

    result = df.differences(df2)

    indices = pd.Index([0, 2])
    columns = pd.MultiIndex.from_product([["col1", "col3"], ["self", "other"]])
    expected = pd.DataFrame(
        [["a", "c", np.nan, np.nan], [np.nan, np.nan, 3.0, np.nan]],
        index=indices,
        columns=columns,
    )
    tm.assert_frame_equal(result, expected)
