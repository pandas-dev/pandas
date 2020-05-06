import numpy as np
import pytest

from pandas import Series
import pandas._testing as tm


@pytest.mark.parametrize(
    "keep, expected",
    [
        ("first", Series([False, False, True, False, True], name="name")),
        ("last", Series([True, True, False, False, False], name="name")),
        (False, Series([True, True, True, False, True], name="name")),
    ],
)
def test_duplicated_keep(keep, expected):
    ser = Series(["a", "b", "b", "c", "a"], name="name")

    result = ser.duplicated(keep=keep)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "keep, expected",
    [
        ("first", Series([False, False, True, False, True])),
        ("last", Series([True, True, False, False, False])),
        (False, Series([True, True, True, False, True])),
    ],
)
def test_duplicated_nan_none(keep, expected):
    ser = Series([np.nan, 3, 3, None, np.nan], dtype=object)

    result = ser.duplicated(keep=keep)
    tm.assert_series_equal(result, expected)
