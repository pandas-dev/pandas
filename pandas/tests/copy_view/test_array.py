import numpy as np
import pytest

from pandas.compat.numpy import np_version_gt2

from pandas import (
    DataFrame,
    Series,
    date_range,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array

# -----------------------------------------------------------------------------
# Copy/view behaviour for accessing underlying array of Series/DataFrame


@pytest.mark.parametrize(
    "method",
    [
        lambda ser: ser.values,
        lambda ser: np.asarray(ser.array),
        lambda ser: np.asarray(ser),
        lambda ser: np.array(ser, copy=False),
    ],
    ids=["values", "array", "np.asarray", "np.array"],
)
def test_series_values(request, method):
    ser = Series([1, 2, 3], name="name")
    ser_orig = ser.copy()

    arr = method(ser)

    if request.node.callspec.id == "array":
        # https://github.com/pandas-dev/pandas/issues/63099
        # .array for now does not return a read-only view
        assert arr.flags.writeable is True
        # updating the array updates the series
        arr[0] = 0
        assert ser.iloc[0] == 0
        return

    # .values still gives a view but is read-only
    assert np.shares_memory(arr, get_array(ser, "name"))
    assert arr.flags.writeable is False

    # mutating series through arr therefore doesn't work
    with pytest.raises(ValueError, match="read-only"):
        arr[0] = 0
    tm.assert_series_equal(ser, ser_orig)

    # mutating the series itself still works
    ser.iloc[0] = 0
    assert ser.values[0] == 0


@pytest.mark.parametrize(
    "method",
    [
        lambda df: df.values,
        lambda df: np.asarray(df),
        lambda ser: np.array(ser, copy=False),
    ],
    ids=["values", "asarray", "array"],
)
def test_dataframe_values(method):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    df_orig = df.copy()

    arr = method(df)

    # .values still gives a view but is read-only
    assert np.shares_memory(arr, get_array(df, "a"))
    assert arr.flags.writeable is False

    # mutating series through arr therefore doesn't work
    with pytest.raises(ValueError, match="read-only"):
        arr[0, 0] = 0
    tm.assert_frame_equal(df, df_orig)

    # mutating the series itself still works
    df.iloc[0, 0] = 0
    assert df.values[0, 0] == 0


def test_series_to_numpy():
    ser = Series([1, 2, 3], name="name")
    ser_orig = ser.copy()

    # default: copy=False, no dtype or NAs
    arr = ser.to_numpy()
    # to_numpy still gives a view but is read-only
    assert np.shares_memory(arr, get_array(ser, "name"))
    assert arr.flags.writeable is False

    # mutating series through arr therefore doesn't work
    with pytest.raises(ValueError, match="read-only"):
        arr[0] = 0
    tm.assert_series_equal(ser, ser_orig)

    # mutating the series itself still works
    ser.iloc[0] = 0
    assert ser.values[0] == 0

    # specify copy=True gives a writeable array
    ser = Series([1, 2, 3], name="name")
    arr = ser.to_numpy(copy=True)
    assert not np.shares_memory(arr, get_array(ser, "name"))
    assert arr.flags.writeable is True

    # specifying a dtype that already causes a copy also gives a writeable array
    ser = Series([1, 2, 3], name="name")
    arr = ser.to_numpy(dtype="float64")
    assert not np.shares_memory(arr, get_array(ser, "name"))
    assert arr.flags.writeable is True


@pytest.mark.parametrize(
    "method",
    [
        lambda ser: np.asarray(ser.values),
        lambda ser: np.asarray(ser.array),
        lambda ser: np.asarray(ser),
        lambda ser: np.asarray(ser, dtype="int64"),
        lambda ser: np.array(ser, copy=False),
    ],
    ids=["values", "array", "np.asarray", "np.asarray-dtype", "np.array"],
)
def test_series_values_ea_dtypes(request, method):
    ser = Series([1, 2, 3], dtype="Int64")
    ser_orig = ser.copy()

    arr = method(ser)

    if request.node.callspec.id in ("values", "array"):
        # https://github.com/pandas-dev/pandas/issues/63099
        # .array/values for now does not return a read-only view
        assert arr.flags.writeable is True
        # updating the array updates the series
        arr[0] = 0
        assert ser.iloc[0] == 0
        return

    # conversion to ndarray gives a view but is read-only
    assert np.shares_memory(arr, get_array(ser))
    assert arr.flags.writeable is False

    # mutating series through arr therefore doesn't work
    with pytest.raises(ValueError, match="read-only"):
        arr[0] = 0
    tm.assert_series_equal(ser, ser_orig)

    # mutating the series itself still works
    ser.iloc[0] = 0
    assert ser.values[0] == 0


@pytest.mark.parametrize(
    "method",
    [
        lambda df: df.values,
        lambda df: np.asarray(df),
        lambda df: np.asarray(df, dtype="int64"),
        lambda df: np.array(df, copy=False),
    ],
    ids=["values", "np.asarray", "np.asarray-dtype", "np.array"],
)
def test_dataframe_array_ea_dtypes(method):
    df = DataFrame({"a": [1, 2, 3]}, dtype="Int64")
    arr = method(df)

    assert np.shares_memory(arr, get_array(df, "a"))
    assert arr.flags.writeable is False


def test_dataframe_array_string_dtype():
    df = DataFrame({"a": ["a", "b"]}, dtype="string[python]")
    arr = np.asarray(df)
    assert np.shares_memory(arr, get_array(df, "a"))
    assert arr.flags.writeable is False


def test_series_array_string_dtype(any_string_dtype):
    ser = Series(["a", "b"], dtype=any_string_dtype)
    arr = np.asarray(ser)
    if any_string_dtype == "string" and any_string_dtype.storage == "pyarrow":
        # for pyarrow strings, the numpy arrays is not a view, so also does
        # not need to be read-only (https://github.com/pandas-dev/pandas/pull/64035)
        assert not np.shares_memory(arr, get_array(ser))
        assert arr.flags.writeable is True
    else:
        assert np.shares_memory(arr, get_array(ser))
        assert arr.flags.writeable is False


def test_dataframe_multiple_numpy_dtypes():
    df = DataFrame({"a": [1, 2, 3], "b": 1.5})
    arr = np.asarray(df)
    assert not np.shares_memory(arr, get_array(df, "a"))
    assert arr.flags.writeable is True

    if np_version_gt2:
        # copy=False semantics are only supported in NumPy>=2.

        with pytest.raises(ValueError, match="Unable to avoid copy while creating"):
            arr = np.array(df, copy=False)

    arr = np.array(df, copy=True)
    assert arr.flags.writeable is True


def test_dataframe_single_block_copy_true():
    # the copy=False/None cases are tested above in test_dataframe_values
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    arr = np.array(df, copy=True)
    assert not np.shares_memory(arr, get_array(df, "a"))
    assert arr.flags.writeable is True


def test_values_is_ea():
    df = DataFrame({"a": date_range("2012-01-01", periods=3)})
    arr = np.asarray(df)
    assert arr.flags.writeable is False


def test_empty_dataframe():
    df = DataFrame()
    arr = np.asarray(df)
    assert arr.flags.writeable is True
