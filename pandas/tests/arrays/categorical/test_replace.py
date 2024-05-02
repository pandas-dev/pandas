import pytest

import pandas as pd
from pandas import Categorical
import pandas._testing as tm


@pytest.mark.parametrize(
    "to_replace,value,expected,flip_categories",
    [
        # one-to-one
        (4, 1, [1, 2, 3], False),
        # many-to-one
        ((5, 6), 2, [1, 2, 3], False),
    ],
)
def test_replace_categorical_series(to_replace, value, expected, flip_categories):
    # GH 31720

    ser = pd.Series([1, 2, 3], dtype="category")
    msg = "with CategoricalDtype is not supported"
    with pytest.raises(TypeError, match=msg):
        ser.replace(to_replace, value)


@pytest.mark.parametrize(
    "to_replace, value, result, expected_error_msg",
    [
        ("b", "c", ["a", "c"], "Categorical.categories are different"),
        # https://github.com/pandas-dev/pandas/issues/33288
        ("b", None, ["a", None], "Categorical.categories length are different"),
    ],
)
def test_replace_categorical(to_replace, value, result, expected_error_msg):
    # GH#26988
    cat = Categorical(["a", "b"])
    expected = Categorical(result)

    if expected_error_msg is None:
        result = pd.Series(cat, copy=False).replace(to_replace, value)._values
        tm.assert_categorical_equal(result, expected)
    elif value is not None:
        result = (
            pd.Series(cat, copy=False)
            .cat.rename_categories({to_replace: value})
            ._values
        )
        tm.assert_categorical_equal(result, expected)

    if to_replace == "b":  # the "c" test is supposed to be unchanged
        with pytest.raises(AssertionError, match=expected_error_msg):
            # ensure non-inplace call does not affect original
            tm.assert_categorical_equal(cat, expected)

    ser = pd.Series(cat, copy=False)
    if expected_error_msg is None:
        ser.replace(to_replace, value, inplace=True)
        tm.assert_categorical_equal(cat, expected)
    else:
        msg2 = "with CategoricalDtype is not supported"
        with pytest.raises(TypeError, match=msg2):
            ser.replace(to_replace, value, inplace=True)


def test_replace_categorical_ea_dtype_raises():
    # GH49404
    cat = Categorical(pd.array(["a", "b"], dtype="string"))
    msg2 = "with CategoricalDtype is not supported"
    with pytest.raises(TypeError, match=msg2):
        pd.Series(cat).replace(["a", "b"], ["c", pd.NA])._values
