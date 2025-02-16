import pytest

import pandas as pd
from pandas import Categorical
import pandas._testing as tm


@pytest.mark.parametrize(
    "to_replace,value,expected",
    [
        # one-to-one
        (4, 1, [1, 2, 3]),
        (3, 1, [1, 2, 1]),
        # many-to-one
        ((5, 6), 2, [1, 2, 3]),
        ((3, 2), 1, [1, 1, 1]),
    ],
)
def test_replace_categorical_series(to_replace, value, expected):
    # GH 31720
    ser = pd.Series([1, 2, 3], dtype="category")
    result = ser.replace(to_replace, value)
    expected = pd.Series(Categorical(expected, categories=[1, 2, 3]))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "to_replace,value",
    [
        # one-to-one
        (3, 5),
        # many-to-one
        ((3, 2), 5),
    ],
)
def test_replace_categorical_series_new_category_raises(to_replace, value):
    # GH 31720
    ser = pd.Series([1, 2, 3], dtype="category")
    with pytest.raises(
        TypeError, match="Cannot setitem on a Categorical with a new category"
    ):
        ser.replace(to_replace, value)


def test_replace_maintain_ordering():
    # GH51016
    dtype = pd.CategoricalDtype([0, 1, 2], ordered=True)
    ser = pd.Series([0, 1, 2], dtype=dtype)
    result = ser.replace(0, 2)
    expected = pd.Series([2, 1, 2], dtype=dtype)
    tm.assert_series_equal(expected, result, check_category_order=True)


def test_replace_categorical_ea_dtype():
    # GH49404
    cat = Categorical(pd.array(["a", "b", "c"], dtype="string"))
    result = pd.Series(cat).replace(["a", "b"], ["c", "c"])._values
    expected = Categorical(
        pd.array(["c"] * 3, dtype="string"),
        categories=pd.array(["a", "b", "c"], dtype="string"),
    )
    tm.assert_categorical_equal(result, expected)


def test_replace_categorical_ea_dtype_different_cats_raises():
    # GH49404
    cat = Categorical(pd.array(["a", "b"], dtype="string"))
    with pytest.raises(
        TypeError, match="Cannot setitem on a Categorical with a new category"
    ):
        pd.Series(cat).replace(["a", "b"], ["c", pd.NA])
