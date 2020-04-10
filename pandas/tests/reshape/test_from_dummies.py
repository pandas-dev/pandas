import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "dtype, expected_dict",
    [
        ("str", {"col1": ["a", "a", "b"]}),
        (str, {"col1": ["a", "a", "b"]},),
        ("category", {"col1": ["a", "a", "b"]}),
    ],
)
def test_dtype(dtype, expected_dict):
    df = pd.DataFrame({"col1_a": [1, 1, 0], "col1_b": [0, 0, 1]})
    result = pd.from_dummies(df, dtype=dtype)
    expected = pd.DataFrame(expected_dict, dtype=dtype)
    tm.assert_frame_equal(result, expected)


def test_malformed():
    df = pd.DataFrame({"col1_a": [1, 1, 0], "col1_b": [1, 0, 1]})
    msg = (
        "Data cannot be decoded! Each row must contain only 0s and 1s"
        ", and each row may have at most one 1"
    )
    with pytest.raises(ValueError, match=msg):
        pd.from_dummies(df)


@pytest.mark.parametrize(
    "prefix_sep, input_dict",
    [
        ("_", {"col1_a": [1, 1, 0], "col1_b": [0, 0, 1]}),
        ("*", {"col1*a": [1, 1, 0], "col1*b": [0, 0, 1]}),
        (".", {"col1.a": [1, 1, 0], "col1.b": [0, 0, 1]}),
    ],
)
def test_prefix_sep(prefix_sep, input_dict):
    df = pd.DataFrame(input_dict)
    result = pd.from_dummies(df, prefix_sep=prefix_sep)
    expected = pd.DataFrame({"col1": ["a", "a", "b"]}, dtype="category")
    tm.assert_frame_equal(result, expected)


def test_no_prefix():
    df = pd.DataFrame({"a": [1, 1, 0], "b": [0, 0, 1]})
    result = pd.from_dummies(df, prefix="letter")
    expected = pd.DataFrame({"letter": ["a", "a", "b"]}, dtype="category")
    tm.assert_frame_equal(result, expected)


def test_multiple_columns():
    df = pd.DataFrame(
        {"col1_a": [1, 0], "col1_b": [0, 1], "col2_a": [0, 0], "col2_c": [1, 1]}
    )
    result = pd.from_dummies(df)
    expected = pd.DataFrame({"col1": ["a", "b"], "col2": ["c", "c"]}, dtype="category")
    tm.assert_frame_equal(result, expected)
