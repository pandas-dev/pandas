from pandas import (
    DatetimeIndex,
    Series,
    Timestamp,
    date_range,
)
import pandas._testing as tm


def test_datetimeindex_series(using_copy_on_write):
    dt = date_range("2019-12-31", periods=3, freq="D")
    ser = Series(dt)
    idx = DatetimeIndex(ser)
    expected = idx.copy(deep=True)
    ser.iloc[0] = Timestamp("2020-12-31")
    if using_copy_on_write:
        tm.assert_index_equal(idx, expected)


def test_datetimeindex_index(using_copy_on_write):
    dt = date_range("2019-12-31", periods=3, freq="D")
    ser = Series(dt)
    idx = DatetimeIndex(DatetimeIndex(ser))
    expected = idx.copy(deep=True)
    ser.iloc[0] = Timestamp("2020-12-31")
    if using_copy_on_write:
        tm.assert_index_equal(idx, expected)


def test_datetimeindex_tz_convert(using_copy_on_write):
    dt = date_range("2019-12-31", periods=3, freq="D", tz="Europe/Berlin")
    ser = Series(dt)
    idx = DatetimeIndex(ser).tz_convert("US/Eastern")
    expected = idx.copy(deep=True)
    ser.iloc[0] = Timestamp("2020-12-31", tz="Europe/Berlin")
    if using_copy_on_write:
        tm.assert_index_equal(idx, expected)


def test_datetimeindex_tz_localize(using_copy_on_write):
    dt = date_range("2019-12-31", periods=3, freq="D")
    ser = Series(dt)
    idx = DatetimeIndex(ser).tz_localize("Europe/Berlin")
    expected = idx.copy(deep=True)
    ser.iloc[0] = Timestamp("2020-12-31")
    if using_copy_on_write:
        tm.assert_index_equal(idx, expected)


def test_datetimeindex_isocalendar(using_copy_on_write):
    dt = date_range("2019-12-31", periods=3, freq="D")
    ser = Series(dt)
    df = DatetimeIndex(ser).isocalendar()
    expected = df.index.copy(deep=True)
    ser.iloc[0] = Timestamp("2020-12-31")
    if using_copy_on_write:
        tm.assert_index_equal(df.index, expected)
