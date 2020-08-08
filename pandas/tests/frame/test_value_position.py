import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize("keep", ["first", "last", "all", "invalid_parameter"])
def test_idxmax(keep):
    # GH#35257
    df = pd.DataFrame(
        {
            "a": [1, 5, 2, 4, 3, 5, 4, 2],
            "b": [9, 4, 3, 4, 2, 2, 4, 9],
            "c": [1, 8, 2, 4, 1, 5, 8, 3],
        },
        index=["A", "B", "C", "D", "E", "F", "G", "H"],
    )

    if keep == "first":
        expected = pd.Series(["B", "A", "B"], index=["a", "b", "c"])
    elif keep == "last":
        expected = pd.Series(["F", "H", "G"], index=["a", "b", "c"])
    elif keep == "all":
        expected = pd.Series(
            [["B", "F"], ["A", "H"], ["B", "G"]], index=["a", "b", "c"]
        )
    else:
        try:
            result = df.idxmax(keep=keep)
        except ValueError:
            return True
    result = df.idxmax(keep=keep)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("keep", ["first", "last", "all", "invalid_parameter"])
def test_idxmin(keep):
    # GH#35257
    df = pd.DataFrame(
        {
            "a": [1, 5, 2, 4, 3, 5, 4, 2],
            "b": [9, 4, 3, 4, 2, 2, 4, 9],
            "c": [1, 8, 2, 4, 1, 5, 8, 3],
        },
        index=["A", "B", "C", "D", "E", "F", "G", "H"],
    )
    if keep == "first":
        expected = pd.Series(["A", "E", "A"], index=["a", "b", "c"])
    elif keep == "last":
        expected = pd.Series(["A", "F", "E"], index=["a", "b", "c"])
    elif keep == "all":
        expected = pd.Series([["A"], ["E", "F"], ["A", "E"]], index=["a", "b", "c"])
    else:
        try:
            result = df.idxmin(keep=keep)
        except ValueError:
            return True
    result = df.idxmin(keep=keep)
    tm.assert_series_equal(result, expected)
