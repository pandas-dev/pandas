import pytest

from pandas import (
    DatetimeIndex,
    Series,
    Timestamp,
    date_range,
)
import pandas._testing as tm

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
