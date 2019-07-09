import numpy as np
import pytest

import pandas as pd
from pandas.util import testing as tm


def test_error():
    df = pd.DataFrame(
        {"A": pd.Series([[0, 1, 2], np.nan, [], (3, 4)], index=list("abcd")), "B": 1}
    )
    df.columns = list("AA")
    with pytest.raises(ValueError):
        df.explode(subset=list("AA"))


def test_basic():
    df = pd.DataFrame(
        {"A": pd.Series([[0, 1, 2], np.nan, [], (3, 4)], index=list("abcd")), "B": 1}
    )
    result = df.explode(subset=["A"])
    expected = pd.DataFrame(
        {
            "A": pd.Series(
                [0, 1, 2, np.nan, np.nan, 3, 4], index=list("aaabcdd"), dtype=object
            ),
            "B": 1,
        }
    )
    tm.assert_frame_equal(result, expected)


def test_all_columns():
    df = pd.DataFrame(
        {"A": pd.Series([[0, 1, 2], np.nan, [], (3, 4)], index=list("abcd")), "B": 1}
    )
    result = df.explode(subset=["A", "B"])
    expected = pd.DataFrame(
        {
            "A": pd.Series(
                [0, 1, 2, np.nan, np.nan, 3, 4], index=list("aaabcdd"), dtype=object
            ),
            "B": 1,
        }
    )
    tm.assert_frame_equal(result, expected)


def test_multiple_columns():
    df = pd.DataFrame(
        {
            "A": pd.Series([[0, 1, 2], np.nan, [], (3, 4)], index=list("abcd")),
            "B": pd.Series([[0, 1, 2], np.nan, np.nan, 3], index=list("abcd")),
        }
    )
    result = df.explode(subset=["A", "B"])
    expected = pd.DataFrame(
        {
            "A": [0, 0, 0, 1, 1, 1, 2, 2, 2, np.nan, np.nan, 3, 4],
            "B": [0, 1, 2, 0, 1, 2, 0, 1, 2, np.nan, np.nan, 3, 3],
        },
        dtype=object,
        index=list("aaaaaaaaabcdd"),
    )
    tm.assert_frame_equal(result, expected)


def test_usecase():
    # explode a single column
    # gh-10511
    df = pd.DataFrame(
        [[11, range(5), 10], [22, range(3), 20]], columns=["A", "B", "C"]
    ).set_index("C")
    result = df.explode(["B"])

    expected = pd.DataFrame(
        {
            "A": [11, 11, 11, 11, 11, 22, 22, 22],
            "B": np.array([0, 1, 2, 3, 4, 0, 1, 2], dtype=object),
            "C": [10, 10, 10, 10, 10, 20, 20, 20],
        },
        columns=list("ABC"),
    ).set_index("C")

    tm.assert_frame_equal(result, expected)

    # gh-8517
    df = pd.DataFrame(
        [["2014-01-01", "Alice", "A B"], ["2014-01-02", "Bob", "C D"]],
        columns=["dt", "name", "text"],
    )
    result = df.assign(text=df.text.str.split(" ")).explode(["text"])
    expected = pd.DataFrame(
        [
            ["2014-01-01", "Alice", "A"],
            ["2014-01-01", "Alice", "B"],
            ["2014-01-02", "Bob", "C"],
            ["2014-01-02", "Bob", "D"],
        ],
        columns=["dt", "name", "text"],
        index=[0, 0, 1, 1],
    )
    tm.assert_frame_equal(result, expected)
