import numpy as np
import pytest

from pandas import (
    DatetimeIndex,
    Series,
    Timestamp,
    array,
    date_range,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array

pytestmark = pytest.mark.filterwarnings(
    "ignore:Setting a value on a view:FutureWarning"
)


@pytest.mark.parametrize("box", [lambda x: x, DatetimeIndex])
def test_datetimeindex(box):
    dt = date_range("2019-12-31", periods=3, freq="D")
    ser = Series(dt)
    idx = box(DatetimeIndex(ser))
    expected = idx.copy(deep=True)
    ser.iloc[0] = Timestamp("2020-12-31")
    tm.assert_index_equal(idx, expected)


def test_datetimeindex_tz_convert():
    dt = date_range("2019-12-31", periods=3, freq="D", tz="Europe/Berlin")
    ser = Series(dt)
    idx = DatetimeIndex(ser).tz_convert("US/Eastern")
    expected = idx.copy(deep=True)
    ser.iloc[0] = Timestamp("2020-12-31", tz="Europe/Berlin")
    tm.assert_index_equal(idx, expected)


def test_datetimeindex_tz_localize():
    dt = date_range("2019-12-31", periods=3, freq="D")
    ser = Series(dt)
    idx = DatetimeIndex(ser).tz_localize("Europe/Berlin")
    expected = idx.copy(deep=True)
    ser.iloc[0] = Timestamp("2020-12-31")
    tm.assert_index_equal(idx, expected)


def test_datetimeindex_isocalendar():
    dt = date_range("2019-12-31", periods=3, freq="D")
    ser = Series(dt)
    df = DatetimeIndex(ser).isocalendar()
    expected = df.index.copy(deep=True)
    ser.iloc[0] = Timestamp("2020-12-31")
    tm.assert_index_equal(df.index, expected)


def test_index_values():
    idx = date_range("2019-12-31", periods=3, freq="D")
    result = idx.values
    assert result.flags.writeable is False


def test_constructor_copy_input_datetime_ndarray_default():
    # GH 63388
    arr = np.array(["2020-01-01", "2020-01-02"], dtype="datetime64[ns]")
    idx = DatetimeIndex(arr)
    assert not np.shares_memory(arr, get_array(idx))


def test_constructor_copy_input_datetime_ea_default():
    # GH 63388
    arr = array(["2020-01-01", "2020-01-02"], dtype="datetime64[ns]")
    idx = DatetimeIndex(arr)
    assert not tm.shares_memory(arr, idx.array)


def test_series_from_temporary_datetimeindex_readonly_data():
    # GH 63388
    arr = np.array(["2020-01-01", "2020-01-02"], dtype="datetime64[ns]")
    arr.flags.writeable = False
    ser = Series(DatetimeIndex(arr))
    assert not np.shares_memory(arr, get_array(ser))
    ser.iloc[0] = Timestamp("2020-01-01")
    expected = Series(
        [Timestamp("2020-01-01"), Timestamp("2020-01-02")], dtype="datetime64[ns]"
    )
    tm.assert_series_equal(ser, expected)
