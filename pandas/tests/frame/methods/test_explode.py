import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


def test_error():
    df = pd.DataFrame(
        {"A": pd.Series([[0, 1, 2], np.nan, [], (3, 4)], index=list("abcd")), "B": 1}
    )
    with pytest.raises(ValueError):
        df.explode(list("AA"))

    df.columns = list("AA")
    with pytest.raises(ValueError):
        df.explode("A")


def test_basic():
    df = pd.DataFrame(
        {"A": pd.Series([[0, 1, 2], np.nan, [], (3, 4)], index=list("abcd")), "B": 1}
    )
    result = df.explode("A")
    expected = pd.DataFrame(
        {
            "A": pd.Series(
                [0, 1, 2, np.nan, np.nan, 3, 4], index=list("aaabcdd"), dtype=object
            ),
            "B": 1,
        }
    )
    tm.assert_frame_equal(result, expected)


def test_multi_index_rows():
    df = pd.DataFrame(
        {"A": np.array([[0, 1, 2], np.nan, [], (3, 4)], dtype=object), "B": 1},
        index=pd.MultiIndex.from_tuples([("a", 1), ("a", 2), ("b", 1), ("b", 2)]),
    )

    result = df.explode("A")
    expected = pd.DataFrame(
        {
            "A": pd.Series(
                [0, 1, 2, np.nan, np.nan, 3, 4],
                index=pd.MultiIndex.from_tuples(
                    [
                        ("a", 1),
                        ("a", 1),
                        ("a", 1),
                        ("a", 2),
                        ("b", 1),
                        ("b", 2),
                        ("b", 2),
                    ]
                ),
                dtype=object,
            ),
            "B": 1,
        }
    )
    tm.assert_frame_equal(result, expected)


def test_multi_index_columns():
    df = pd.DataFrame(
        {("A", 1): np.array([[0, 1, 2], np.nan, [], (3, 4)], dtype=object), ("A", 2): 1}
    )

    result = df.explode(("A", 1))
    expected = pd.DataFrame(
        {
            ("A", 1): pd.Series(
                [0, 1, 2, np.nan, np.nan, 3, 4],
                index=pd.Index([0, 0, 0, 1, 2, 3, 3]),
                dtype=object,
            ),
            ("A", 2): 1,
        }
    )
    tm.assert_frame_equal(result, expected)


def test_usecase():
    # explode a single column
    # gh-10511
    df = pd.DataFrame(
        [[11, range(5), 10], [22, range(3), 20]], columns=list("ABC")
    ).set_index("C")
    result = df.explode("B")

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
    result = df.assign(text=df.text.str.split(" ")).explode("text")
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


@pytest.mark.parametrize(
    "input_dict, input_index, expected_dict, expected_index",
    [
        (
            {"col1": [[1, 2], [3, 4]], "col2": ["foo", "bar"]},
            [0, 0],
            {"col1": [1, 2, 3, 4], "col2": ["foo", "foo", "bar", "bar"]},
            [0, 0, 0, 0],
        ),
        (
            {"col1": [[1, 2], [3, 4]], "col2": ["foo", "bar"]},
            pd.Index([0, 0], name="my_index"),
            {"col1": [1, 2, 3, 4], "col2": ["foo", "foo", "bar", "bar"]},
            pd.Index([0, 0, 0, 0], name="my_index"),
        ),
        (
            {"col1": [[1, 2], [3, 4]], "col2": ["foo", "bar"]},
            pd.MultiIndex.from_arrays(
                [[0, 0], [1, 1]], names=["my_first_index", "my_second_index"]
            ),
            {"col1": [1, 2, 3, 4], "col2": ["foo", "foo", "bar", "bar"]},
            pd.MultiIndex.from_arrays(
                [[0, 0, 0, 0], [1, 1, 1, 1]],
                names=["my_first_index", "my_second_index"],
            ),
        ),
        (
            {"col1": [[1, 2], [3, 4]], "col2": ["foo", "bar"]},
            pd.MultiIndex.from_arrays([[0, 0], [1, 1]], names=["my_index", None]),
            {"col1": [1, 2, 3, 4], "col2": ["foo", "foo", "bar", "bar"]},
            pd.MultiIndex.from_arrays(
                [[0, 0, 0, 0], [1, 1, 1, 1]], names=["my_index", None]
            ),
        ),
    ],
)
def test_duplicate_index(input_dict, input_index, expected_dict, expected_index):
    # GH 28005
    df = pd.DataFrame(input_dict, index=input_index)
    result = df.explode("col1")
    expected = pd.DataFrame(expected_dict, index=expected_index, dtype=object)
    tm.assert_frame_equal(result, expected)
