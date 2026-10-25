import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "keep, expected",
    [
        ("first", [False, False, True, False, True]),
        ("last", [True, True, False, False, False]),
        (False, [True, True, True, False, True]),
    ],
)
def test_duplicated_keep(keep, expected):
    ser = pd.Series(["a", "b", "b", "c", "a"], name="name")

    result = ser.duplicated(keep=keep)
    expected = pd.Series(expected, name="name")
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "keep, expected",
    [
        ("first", [False, False, True, False, True]),
        ("last", [True, True, False, False, False]),
        (False, [True, True, True, False, True]),
    ],
)
def test_duplicated_nan_none(keep, expected):
    ser = pd.Series([np.nan, 3, 3, None, np.nan], dtype=object)

    result = ser.duplicated(keep=keep)
    expected = pd.Series(expected)
    tm.assert_series_equal(result, expected)


def test_duplicated_categorical_bool_na(nulls_fixture):
    # GH#44351
    ser = pd.Series(
        pd.Categorical(
            [True, False, True, False, nulls_fixture],
            categories=[True, False],
            ordered=True,
        )
    )
    result = ser.duplicated()
    expected = pd.Series([False, False, True, True, False])
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "keep, vals",
    [
        ("last", [True, True, False]),
        ("first", [False, True, True]),
        (False, [True, True, True]),
    ],
)
def test_duplicated_mask(keep, vals):
    # GH#48150
    ser = pd.Series([1, 2, pd.NA, pd.NA, pd.NA], dtype="Int64")
    result = ser.duplicated(keep=keep)
    expected = pd.Series([False, False, *vals])
    tm.assert_series_equal(result, expected)


def test_duplicated_mask_no_duplicated_na(keep):
    # GH#48150
    ser = pd.Series([1, 2, pd.NA], dtype="Int64")
    result = ser.duplicated(keep=keep)
    expected = pd.Series([False, False, False])
    tm.assert_series_equal(result, expected)
