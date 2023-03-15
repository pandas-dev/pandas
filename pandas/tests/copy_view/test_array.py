import numpy as np
import pytest

from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array

# -----------------------------------------------------------------------------
# Copy/view behaviour for accessing underlying array of Series/DataFrame


@pytest.mark.parametrize(
    "method",
    [lambda ser: ser.values, lambda ser: np.asarray(ser)],
    ids=["values", "asarray"],
)
def test_series_values(using_copy_on_write, method):
    ser = Series([1, 2, 3], name="name")
    ser_orig = ser.copy()

    arr = method(ser)

    if using_copy_on_write:
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
    else:
        assert arr.flags.writeable is True
        arr[0] = 0
        assert ser.iloc[0] == 0


@pytest.mark.parametrize(
    "method",
    [lambda df: df.values, lambda df: np.asarray(df)],
    ids=["values", "asarray"],
)
def test_dataframe_values(using_copy_on_write, using_array_manager, method):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    df_orig = df.copy()

    arr = method(df)

    if using_copy_on_write:
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
    else:
        assert arr.flags.writeable is True
        arr[0, 0] = 0
        if not using_array_manager:
            assert df.iloc[0, 0] == 0
        else:
            tm.assert_frame_equal(df, df_orig)


def test_series_to_numpy(using_copy_on_write):
    ser = Series([1, 2, 3], name="name")
    ser_orig = ser.copy()

    # default: copy=False, no dtype or NAs
    arr = ser.to_numpy()
    if using_copy_on_write:
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
    else:
        assert arr.flags.writeable is True
        arr[0] = 0
        assert ser.iloc[0] == 0

    # specify copy=False gives a writeable array
    ser = Series([1, 2, 3], name="name")
    arr = ser.to_numpy(copy=True)
    assert not np.shares_memory(arr, get_array(ser, "name"))
    assert arr.flags.writeable is True

    # specifying a dtype that already causes a copy also gives a writeable array
    ser = Series([1, 2, 3], name="name")
    arr = ser.to_numpy(dtype="float64")
    assert not np.shares_memory(arr, get_array(ser, "name"))
    assert arr.flags.writeable is True
