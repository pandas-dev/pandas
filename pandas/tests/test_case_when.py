import numpy as np
import pytest

from pandas import (
    DataFrame,
    Series,
    case_when,
)
import pandas._testing as tm


@pytest.fixture
def df():
    """
    base dataframe for testing
    """
    return DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})


# use fixture and parametrize
def test_case_when_no_args():
    """
    Raise ValueError if no args is provided.
    """
    msg = "Kindly provide at least one boolean condition, "
    msg += "with a corresponding replacement."
    with pytest.raises(ValueError, match=msg):
        case_when()


def test_case_when_odd_args(df):
    """
    Raise ValueError if no of args is odd.
    """
    msg = "The number of boolean conditions should be equal "
    msg += "to the number of replacements. "
    msg += "However, the total number of conditions and replacements "
    msg += "is 3, which is an odd number."
    with pytest.raises(ValueError, match=msg):
        case_when(df["a"].eq(1), 1, df.a.gt(1))


def test_case_when_boolean_ndim(df):
    """
    Raise ValueError if boolean array ndim is greater than 1.
    """
    with pytest.raises(ValueError, match="condition0 is not a one dimensional array."):
        case_when(df, 2)


def test_case_when_not_boolean(df):
    """
    Raise TypeError if condition is not a boolean array.
    """
    with pytest.raises(TypeError, match="condition1 is not a boolean array."):
        case_when(df["a"].eq(1), 1, df["a"], 2)


def test_case_when_multiple_bool_lengths(df):
    """
    Raise ValueError if the boolean conditions do not have the same length.
    """
    with pytest.raises(
        ValueError, match="All boolean conditions should have the same length."
    ):
        case_when(df["a"].eq(1), 1, [True, False], 2)


def test_case_when_default_ndim(df):
    """
    Raise ValueError if the default is not a 1D array.
    """
    msg = "The provided default argument should "
    msg += "either be a scalar or a 1-D array."
    with pytest.raises(ValueError, match=msg):
        case_when(df["a"].eq(1), 1, default=df)


def test_case_when_default_length(df):
    """
    Raise ValueError if the default is not
    the same length as any of the boolean conditions.
    """
    msg = "The length of the default argument should "
    msg += "be the same as the length of any "
    msg += "of the boolean conditions."
    with pytest.raises(ValueError, match=msg):
        case_when(df["a"].eq(1), 1, default=[2])


def test_case_when_replacement_ndim(df):
    """
    Raise ValueError if the ndim of the replacement value is greater than 1.
    """
    with pytest.raises(ValueError, match="replacement0 should be a 1-D array."):
        case_when(df["a"].eq(1), df)


def test_case_when_replacement_length(df):
    """
    Raise ValueError if the replacement size is not the same as the boolean array.
    """
    msg = "The size of replacement0 array"
    msg += "does not match the size of condition0 array."
    with pytest.raises(ValueError, match=msg):
        case_when(df["a"].eq(1), np.array([2]))


def test_case_when_multiple_conditions(df):
    """
    Test output when booleans are derived from a computation
    """
    result = case_when(df.a.eq(1), 1, Series([False, True, False]), 2)
    expected = Series([1, 2, np.nan])
    tm.assert_series_equal(result, expected)


def test_case_when_multiple_conditions_replacement_list(df):
    """
    Test output when replacement is a list
    """
    result = case_when(
        [True, False, False], 1, df["a"].gt(1) & df["b"].eq(5), [1, 2, 3]
    )
    expected = Series([1, 2, np.nan], dtype="Int64")
    tm.assert_series_equal(result, expected)


def test_case_when_multiple_conditions_replacement_series(df):
    """
    Test output when replacement is a Series
    """
    result = case_when(
        np.array([True, False, False]),
        1,
        df["a"].gt(1) & df["b"].eq(5),
        Series([1, 2, 3]),
    )
    expected = Series([1, 2, np.nan])
    tm.assert_series_equal(result, expected)


def test_case_when_multiple_conditions_default_is_not_none(df):
    """
    Test output when default is not None
    """
    result = case_when(
        [True, False, False],
        1,
        df["a"].gt(1) & df["b"].eq(5),
        Series([1, 2, 3]),
        default=-1,
    )
    expected = Series([1, 2, -1])
    tm.assert_series_equal(result, expected)
