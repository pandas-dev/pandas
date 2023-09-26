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
def test_case_when_multiple_conditions(df):
    """
    Test output when booleans are derived from a computation
    """
    result = case_when(df.a.eq(1), 1, df.a.gt(1) & df.b.eq(5), 2)
    expected = Series([1, 2, np.nan])
    tm.assert_series_equal(result, expected)


def test_case_when_multiple_conditions_array_series(df):
    """
    Test output when boolean arrays are passed
    """
    result = case_when([True, False, False], 1, Series([False, True, False]), 2)
    expected = Series([1, 2, np.nan])
    tm.assert_series_equal(result, expected)
