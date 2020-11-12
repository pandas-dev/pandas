import numpy as np
import pytest

import pandas as pd
from pandas import Categorical
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
        #  except for cases with only 1 unique entry in `value`
        ([1], [2], [2, 2, 3], True),
        ([1, 4], [5, 2], [5, 2, 3], True),
        # check_categorical sorts categories, which crashes on mixed dtypes
        (3, "4", [1, 2, "4"], False),
        ([1, 2, "3"], "5", ["5", "5", 3], True),
    ],
)
def test_replace(to_replace, value, expected, flip_categories):
    # GH 31720
    stays_categorical = not isinstance(value, list) or len(pd.unique(value)) == 1

    s = pd.Series([1, 2, 3], dtype="category")
    result = s.replace(to_replace, value)
    expected = pd.Series(expected, dtype="category")
    s.replace(to_replace, value, inplace=True)

    if flip_categories:
        expected = expected.cat.set_categories(expected.cat.categories[::-1])

    if not stays_categorical:
        # the replace call loses categorical dtype
        expected = pd.Series(np.asarray(expected))

    tm.assert_series_equal(expected, result, check_category_order=False)
    tm.assert_series_equal(expected, s, check_category_order=False)


@pytest.mark.parametrize(
    "to_replace, value, result, expected_error_msg",
    [
        ("b", "c", ["a", "c"], "Categorical.categories are different"),
        ("c", "d", ["a", "b"], None),
        # https://github.com/pandas-dev/pandas/issues/33288
        ("a", "a", ["a", "b"], None),
        ("b", None, ["a", None], "Categorical.categories length are different"),
    ],
)
def test_replace2(to_replace, value, result, expected_error_msg):
    # TODO: better name
    # GH#26988
    cat = Categorical(["a", "b"])
    expected = Categorical(result)
    result = cat.replace(to_replace, value)
    tm.assert_categorical_equal(result, expected)
    if to_replace == "b":  # the "c" test is supposed to be unchanged
        with pytest.raises(AssertionError, match=expected_error_msg):
            # ensure non-inplace call does not affect original
            tm.assert_categorical_equal(cat, expected)
    cat.replace(to_replace, value, inplace=True)
    tm.assert_categorical_equal(cat, expected)


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
    df = pd.DataFrame({"col1": input_data}, dtype="string")
    expected = pd.DataFrame({"col1": expected_data}, dtype="string")
    df_replaced = df.replace(to_replace, value, inplace=inplace, regex=True)
    result = df if inplace else df_replaced
    tm.assert_frame_equal(result, expected)


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
    df = pd.DataFrame({"col1": input_data}, dtype="string")
    expected = pd.DataFrame({"col1": expected_data}, dtype="string")
    result = df.replace(to_replace, value, inplace=False, regex=False)
    tm.assert_frame_equal(result, expected)
