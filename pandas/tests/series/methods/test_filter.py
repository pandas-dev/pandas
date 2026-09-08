import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.fixture
def ser():
    return pd.Series([1, 2, 3], index=["x", "y", "z"])


@pytest.mark.parametrize(
    "mask",
    [
        [False, True, True],
        np.array([False, True, True]),
        pd.array([False, True, True], dtype="boolean"),
        lambda ser: ser > 1,
    ],
)
def test_filter_mask(ser, mask):
    # GH#61317
    result = ser.filter(mask)
    expected = ser.iloc[1:]
    tm.assert_series_equal(result, expected)


def test_filter_mask_series_aligns(ser):
    # GH#61317
    mask = pd.Series([True, True, False], index=["z", "y", "x"])
    result = ser.filter(mask)
    expected = ser.iloc[1:]
    tm.assert_series_equal(result, expected)


def test_filter_expression_raises(ser):
    # GH#61317
    msg = "Expressions such as pd.col\\(...\\) are only supported by DataFrame.filter"
    with pytest.raises(TypeError, match=msg):
        ser.filter(pd.col("a") > 1)


def test_filter_callable_must_return_mask(ser):
    # GH#61317
    msg = "The callable passed to Series.filter must evaluate to a boolean mask"
    with pytest.raises(TypeError, match=msg):
        ser.filter(lambda ser: ["x"])


def test_filter_mask_2d_raises(ser):
    # GH#61317
    msg = "The mask passed to Series.filter must be one-dimensional"
    with pytest.raises(ValueError, match=msg):
        ser.filter(ser.to_frame() > 1)


@pytest.mark.parametrize(
    "mask",
    [
        [True, None, False],
        pd.array([True, None, False], dtype="boolean"),
    ],
)
def test_filter_mask_na(ser, mask):
    # GH#61317
    result = ser.filter(mask)
    expected = ser.iloc[[0]]
    tm.assert_series_equal(result, expected)

    result = ser.filter(mask, na=False)
    tm.assert_series_equal(result, expected)

    msg = "The mask contains missing values"
    with pytest.raises(ValueError, match=msg):
        ser.filter(mask, na="raise")

    result = ser.filter(mask, na=True)
    expected = ser.iloc[[0, 1]]
    tm.assert_series_equal(result, expected)


def test_filter_bool_labels_select_labels():
    # GH#61317
    # Boolean values on an axis with boolean labels keep selecting labels
    ser = pd.Series([1, 2, 3], index=[True, False, "c"])
    result = ser.filter([True, False])
    expected = ser.iloc[:2]
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("kwargs", [{"items": ["x"]}, {"like": "x"}, {"regex": "x"}])
def test_filter_labels(ser, kwargs):
    # GH#61317
    result = ser.filter(**kwargs)
    expected = ser.iloc[[0]]
    tm.assert_series_equal(result, expected)
