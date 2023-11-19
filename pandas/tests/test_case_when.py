import numpy as np
import pytest

from pandas import (
    DataFrame,
    Series,
    array as pd_array,
    case_when,
    date_range,
)
import pandas._testing as tm


@pytest.fixture
def df():
    """
    base dataframe for testing
    """
    return DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})


def test_case_when_no_args():
    """
    Raise ValueError if no args is provided.
    """
    msg = "provide at least one boolean condition, "
    msg += "with a corresponding replacement."
    with pytest.raises(ValueError, match=msg):
        case_when()


def test_case_when_odd_args(df):
    """
    Raise ValueError if no of args is odd.
    """
    msg = "Argument 0 must have length 2; "
    msg += "a condition and replacement; got length 3."

    with pytest.raises(ValueError, match=msg):
        case_when((df["a"].eq(1), 1, df.a.gt(1)))


def test_case_when_raise_error_from_mask(df):
    """
    Raise Error from within Series.mask
    """
    msg = "condition0 and replacement0 failed to evaluate."
    with pytest.raises(ValueError, match=msg):
        case_when((df["a"].eq(1), df))


def test_case_when_error_multiple_replacements_series(df):
    """
    Test output when the replacements indices are different.
    """
    with pytest.raises(
        AssertionError, match="All replacement objects must have the same index."
    ):
        case_when(
            ([True, False, False], Series(1)),
            (df["a"].gt(1) & df["b"].eq(5), Series([1, 2, 3])),
        )


def test_case_when_error_multiple_conditions_series(df):
    """
    Test output when the conditions indices are different.
    """
    with pytest.raises(
        AssertionError, match="All condition objects must have the same index."
    ):
        case_when(
            (Series([True, False, False], index=[2, 3, 4]), 1),
            (df["a"].gt(1) & df["b"].eq(5), Series([1, 2, 3])),
        )


def test_case_when_raise_error_different_index_condition_and_replacements(df):
    """
    Raise if the replacement index and condition index are different.
    """
    msg = "All replacement objects and condition objects "
    msg += "should have the same index."
    with pytest.raises(AssertionError, match=msg):
        case_when(
            (df.a.eq(1), 1), (Series([False, True, False], index=["a", "b", "c"]), 2)
        )


def test_case_when_single_condition(df):
    """
    Test output on a single condition.
    """
    result = case_when((df.a.eq(1), 1))
    expected = Series([1, np.nan, np.nan])
    tm.assert_series_equal(result, expected)


def test_case_when_multiple_conditions(df):
    """
    Test output when booleans are derived from a computation
    """
    result = case_when((df.a.eq(1), 1), (Series([False, True, False]), 2))
    expected = Series([1, 2, np.nan])
    tm.assert_series_equal(result, expected)


def test_case_when_multiple_conditions_replacement_list(df):
    """
    Test output when replacement is a list
    """
    result = case_when(
        ([True, False, False], 1), (df["a"].gt(1) & df["b"].eq(5), [1, 2, 3])
    )
    expected = Series([1, 2, np.nan])
    tm.assert_series_equal(result, expected)


def test_case_when_multiple_conditions_replacement_extension_dtype(df):
    """
    Test output when replacement has an extension dtype
    """
    result = case_when(
        ([True, False, False], 1),
        (df["a"].gt(1) & df["b"].eq(5), pd_array([1, 2, 3], dtype="Int64")),
    )
    expected = Series([1, 2, np.nan], dtype="Int64")
    tm.assert_series_equal(result, expected)


def test_case_when_multiple_conditions_replacement_series(df):
    """
    Test output when replacement is a Series
    """
    result = case_when(
        (np.array([True, False, False]), 1),
        (df["a"].gt(1) & df["b"].eq(5), Series([1, 2, 3])),
    )
    expected = Series([1, 2, np.nan])
    tm.assert_series_equal(result, expected)


def test_case_when_multiple_conditions_default_is_not_none(df):
    """
    Test output when default is not None
    """
    result = case_when(
        ([True, False, False], 1),
        (df["a"].gt(1) & df["b"].eq(5), Series([1, 2, 3])),
        default=-1,
    )
    expected = Series([1, 2, -1])
    tm.assert_series_equal(result, expected)


def test_case_when_non_range_index():
    """
    Test output if index is not RangeIndex
    """
    rng = np.random.default_rng(seed=123)
    dates = date_range("1/1/2000", periods=8)
    df = DataFrame(
        rng.standard_normal(size=(8, 4)), index=dates, columns=["A", "B", "C", "D"]
    )
    result = case_when((df.A.gt(0), df.B), default=5)
    result = Series(result, name="A")
    expected = df.A.mask(df.A.gt(0), df.B).where(df.A.gt(0), 5)
    tm.assert_series_equal(result, expected)
