import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "to_replace,value,expected,flip_categories",
    [
        # one-to-one
        (1, 2, [2, 2, 3], False),
        (1, 4, [4, 2, 3], False),
        (4, 1, [1, 2, 3], False),
        (5, 6, [1, 2, 3], False),
        # many-to-one
        ([1], 2, [2, 2, 3], False),
        ([1, 2], 3, [3, 3, 3], False),
        ([1, 2], 4, [4, 4, 3], False),
        ((1, 2, 4), 5, [5, 5, 3], False),
        ((5, 6), 2, [1, 2, 3], False),
        # many-to-many, handled outside of Categorical and results in separate dtype
        ([1], [2], [2, 2, 3], True),
        ([1, 4], [5, 2], [5, 2, 3], True),
        # check_categorical sorts categories, which crashes on mixed dtypes
        (3, "4", [1, 2, "4"], False),
        ([1, 2, "3"], "5", ["5", "5", 3], True),
    ],
)
def test_replace(to_replace, value, expected, flip_categories):
    # GH 31720
    stays_categorical = not isinstance(value, list)

    s = pd.Series([1, 2, 3], dtype="category")
    result = s.replace(to_replace, value)
    expected = pd.Series(expected, dtype="category")
    s.replace(to_replace, value, inplace=True)

    if flip_categories:
        expected = expected.cat.set_categories(expected.cat.categories[::-1])

    if not stays_categorical:
        # the replace call loses categorical dtype
        expected = pd.Series(np.asarray(expected))

    tm.assert_series_equal(
        expected, result, check_category_order=False,
    )
    tm.assert_series_equal(
        expected, s, check_category_order=False,
    )


@pytest.mark.parametrize(
    "to_replace,value,input_data,expected_data,inplace",
    [
        (r"^\s*$", pd.NA, ["d", "ee", "f", ""], ["d", "ee", "f", pd.NA], False),
        (r"e{2}", "replace", ["d", "ee", "f", ""], ["d", "replace", "f", ""], False),
        (r"f", "replace", ["d", "ee", "f", ""], ["d", "ee", "replace", ""], False),
        (r"^\s*$", pd.NA, ["d", "ee", "f", ""], ["d", "ee", "f", pd.NA], True),
        (r"e{2}", "replace", ["d", "ee", "f", ""], ["d", "replace", "f", ""], True),
        (r"f", "replace", ["d", "ee", "f", ""], ["d", "ee", "replace", ""], True),
    ],
)
def test_replace_regex(to_replace, value, input_data, expected_data, inplace):
    # GH35977
    input_df = pd.DataFrame({"col1": input_data}, dtype="string")
    expected_df = pd.DataFrame({"col1": expected_data}, dtype="string")
    converted = input_df.replace(to_replace, value, inplace=inplace, regex=True)
    if inplace:
        tm.assert_frame_equal(expected_df, input_df)
    else:
        tm.assert_frame_equal(expected_df, converted)


@pytest.mark.parametrize(
    "to_replace,value,input_data,expected_data",
    [
        ("", pd.NA, ["d", "ee", "f", ""], ["d", "ee", "f", pd.NA]),
        ("ee", "replace", ["d", "ee", "f", ""], ["d", "replace", "f", ""]),
        ("f", "replace", ["d", "ee", "f", ""], ["d", "ee", "replace", ""]),
    ],
)
def test_replace_string(to_replace, value, input_data, expected_data):
    # GH35977
    input_df = pd.DataFrame({"col1": input_data}, dtype="string")
    expected_df = pd.DataFrame({"col1": expected_data}, dtype="string")
    input_df = input_df.replace(to_replace, value, inplace=False, regex=False)
    tm.assert_frame_equal(expected_df, input_df)
