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
        pd.Series([False, True, True], index=["x", "y", "z"]),
        lambda ser: ser > 1,
    ],
)
def test_filter_cond(ser, mask):
    # GH#61317
    result = ser.filter(cond=mask)
    expected = ser.iloc[1:]
    tm.assert_series_equal(result, expected)


def test_filter_positional_callable_is_mask(ser):
    # GH#61317
    result = ser.filter(lambda ser: ser > 1)
    expected = ser.iloc[1:]
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "mask",
    [
        [True, False, True],
        np.array([True, False, True]),
        pd.Series([True, False, True], index=["x", "y", "z"]),
    ],
)
def test_filter_positional_bools_select_labels(ser, mask):
    # GH#61317
    msg = "A list-like of booleans passed positionally to Series.filter"
    with tm.assert_produces_warning(UserWarning, match=msg):
        result = ser.filter(mask)
    expected = ser.iloc[[]]
    tm.assert_series_equal(result, expected)

    with tm.assert_produces_warning(None):
        result = ser.filter(items=mask)
    tm.assert_series_equal(result, expected)


def test_filter_bool_labels():
    # GH#61317
    ser = pd.Series([1, 2, 3], index=[True, False, "c"])
    result = ser.filter(items=[True, False])
    expected = ser.iloc[:2]
    tm.assert_series_equal(result, expected)

    msg = "A list-like of booleans passed positionally to Series.filter"
    with tm.assert_produces_warning(UserWarning, match=msg):
        result = ser.filter([True, False])
    tm.assert_series_equal(result, expected)

    result = ser.filter(cond=[False, True, True])
    expected = ser.iloc[1:]
    tm.assert_series_equal(result, expected)


def test_filter_cond_series_aligns(ser):
    # GH#61317
    mask = pd.Series([True, True, False], index=["z", "y", "x"])
    result = ser.filter(cond=mask)
    expected = ser.iloc[1:]
    tm.assert_series_equal(result, expected)


def test_filter_expression_raises(ser):
    # GH#61317
    msg = "Expressions such as pd.col\\(...\\) are only supported by DataFrame.filter"
    with pytest.raises(TypeError, match=msg):
        ser.filter(pd.col("a") > 1)
    with pytest.raises(TypeError, match=msg):
        ser.filter(cond=pd.col("a") > 1)


def test_filter_callable_must_return_mask(ser):
    # GH#61317
    msg = "The callable passed to Series.filter must evaluate to a boolean mask"
    with pytest.raises(TypeError, match=msg):
        ser.filter(lambda ser: ["x"])


def test_filter_cond_not_mask_raises(ser):
    # GH#61317
    msg = "cond passed to Series.filter must be a boolean mask"
    with pytest.raises(TypeError, match=msg):
        ser.filter(cond=["x"])


def test_filter_cond_2d_raises(ser):
    # GH#61317
    mask = pd.DataFrame({"a": [False, True, True]}, index=["x", "y", "z"])
    msg = "The mask passed to Series.filter must be one-dimensional"
    with pytest.raises(ValueError, match=msg):
        ser.filter(cond=mask)


@pytest.mark.parametrize(
    "mask",
    [
        [True, None, False],
        pd.array([True, None, False], dtype="boolean"),
    ],
)
def test_filter_cond_na(ser, mask):
    # GH#61317
    result = ser.filter(cond=mask)
    expected = ser.iloc[[0]]
    tm.assert_series_equal(result, expected)

    result = ser.filter(cond=mask, na=False)
    tm.assert_series_equal(result, expected)

    msg = "The mask contains missing values"
    with pytest.raises(ValueError, match=msg):
        ser.filter(cond=mask, na="raise")

    result = ser.filter(cond=mask, na=True)
    expected = ser.iloc[[0, 1]]
    tm.assert_series_equal(result, expected)


def test_filter_na_label():
    # GH#61317
    ser = pd.Series([1, 2], index=[np.nan, "x"])
    with tm.assert_produces_warning(None):
        result = ser.filter([np.nan])
    expected = ser.iloc[[0]]
    tm.assert_series_equal(result, expected, check_index_type=False)


def test_filter_tuple_labels_multiindex():
    # GH#61317
    mi = pd.MultiIndex.from_tuples([(True, False), (False, True)])
    ser = pd.Series([1, 2], index=mi)
    with tm.assert_produces_warning(None):
        result = ser.filter([(True, False)])
    expected = ser.iloc[[0]]
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("kwargs", [{"items": ["x"]}, {"like": "x"}, {"regex": "x"}])
def test_filter_labels(ser, kwargs):
    # GH#61317
    result = ser.filter(**kwargs)
    expected = ser.iloc[[0]]
    tm.assert_series_equal(result, expected)
