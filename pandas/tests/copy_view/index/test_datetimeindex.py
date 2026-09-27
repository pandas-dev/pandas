import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


@pytest.mark.parametrize("box", [lambda x: x, pd.DatetimeIndex])
def test_datetimeindex(box):
    dt = pd.date_range("2019-12-31", periods=3, freq="D")
    ser = pd.Series(dt)
    idx = box(pd.DatetimeIndex(ser))
    expected = idx.copy(deep=True)
    ser.iloc[0] = pd.Timestamp("2020-12-31")
    tm.assert_index_equal(idx, expected)


def test_datetimeindex_tz_convert():
    dt = pd.date_range("2019-12-31", periods=3, freq="D", tz="Europe/Berlin")
    ser = pd.Series(dt)
    idx = pd.DatetimeIndex(ser).tz_convert("US/Eastern")
    expected = idx.copy(deep=True)
    ser.iloc[0] = pd.Timestamp("2020-12-31", tz="Europe/Berlin")
    tm.assert_index_equal(idx, expected)


def test_datetimeindex_tz_localize():
    dt = pd.date_range("2019-12-31", periods=3, freq="D")
    ser = pd.Series(dt)
    idx = pd.DatetimeIndex(ser).tz_localize("Europe/Berlin")
    expected = idx.copy(deep=True)
    ser.iloc[0] = pd.Timestamp("2020-12-31")
    tm.assert_index_equal(idx, expected)


def test_datetimeindex_isocalendar():
    dt = pd.date_range("2019-12-31", periods=3, freq="D")
    ser = pd.Series(dt)
    df = pd.DatetimeIndex(ser).isocalendar()
    expected = df.index.copy(deep=True)
    ser.iloc[0] = pd.Timestamp("2020-12-31")
    tm.assert_index_equal(df.index, expected)


def test_index_values():
    idx = pd.date_range("2019-12-31", periods=3, freq="D")
    result = idx.values
    assert result.flags.writeable is False


def test_constructor_copy_input_datetime_ndarray_default():
    # GH 63388
    arr = np.array(["2020-01-01", "2020-01-02"], dtype="datetime64[ns]")
    idx = pd.DatetimeIndex(arr)
    assert not np.shares_memory(arr, get_array(idx))


def test_constructor_copy_input_datetime_ea_default():
    # GH 63388
    arr = pd.array(["2020-01-01", "2020-01-02"], dtype="datetime64[ns]")
    idx = pd.DatetimeIndex(arr)
    assert not tm.shares_memory(arr, idx.array)


def test_series_from_temporary_datetimeindex_readonly_data():
    # GH 63388
    arr = np.array(["2020-01-01", "2020-01-02"], dtype="datetime64[ns]")
    arr.flags.writeable = False
    ser = pd.Series(pd.DatetimeIndex(arr))
    assert not np.shares_memory(arr, get_array(ser))
    ser.iloc[0] = pd.Timestamp("2020-01-01")
    expected = pd.Series(
        [pd.Timestamp("2020-01-01"), pd.Timestamp("2020-01-02")], dtype="datetime64[ns]"
    )
    tm.assert_series_equal(ser, expected)
