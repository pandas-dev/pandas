import numpy as np

from pandas import (
    DataFrame,
    Series,
    concat,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


def test_concat_frames(using_copy_on_write):
    df = DataFrame({"b": ["a"] * 3})
    df2 = DataFrame({"a": ["a"] * 3})
    df_orig = df.copy()
    result = concat([df, df2], axis=1)

    if using_copy_on_write:
        assert np.shares_memory(get_array(result, "b"), get_array(df, "b"))
        assert np.shares_memory(get_array(result, "a"), get_array(df2, "a"))
    else:
        assert not np.shares_memory(get_array(result, "b"), get_array(df, "b"))
        assert not np.shares_memory(get_array(result, "a"), get_array(df2, "a"))

    result.iloc[0, 0] = "d"
    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "b"), get_array(df, "b"))
        assert np.shares_memory(get_array(result, "a"), get_array(df2, "a"))

    result.iloc[0, 1] = "d"
    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "a"), get_array(df2, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_concat_frames_updating_input(using_copy_on_write):
    df = DataFrame({"b": ["a"] * 3})
    df2 = DataFrame({"a": ["a"] * 3})
    result = concat([df, df2], axis=1)

    if using_copy_on_write:
        assert np.shares_memory(get_array(result, "b"), get_array(df, "b"))
        assert np.shares_memory(get_array(result, "a"), get_array(df2, "a"))
    else:
        assert not np.shares_memory(get_array(result, "b"), get_array(df, "b"))
        assert not np.shares_memory(get_array(result, "a"), get_array(df2, "a"))

    expected = result.copy()
    df.iloc[0, 0] = "d"
    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "b"), get_array(df, "b"))
        assert np.shares_memory(get_array(result, "a"), get_array(df2, "a"))

    df2.iloc[0, 0] = "d"
    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "a"), get_array(df2, "a"))
    tm.assert_frame_equal(result, expected)


def test_concat_series(using_copy_on_write):
    ser = Series([1, 2], name="a")
    ser2 = Series([3, 4], name="b")
    ser_orig = ser.copy()
    ser2_orig = ser2.copy()
    result = concat([ser, ser2], axis=1)

    if using_copy_on_write:
        assert np.shares_memory(get_array(result, "a"), ser.values)
        assert np.shares_memory(get_array(result, "b"), ser2.values)
    else:
        assert not np.shares_memory(get_array(result, "a"), ser.values)
        assert not np.shares_memory(get_array(result, "b"), ser2.values)

    result.iloc[0, 0] = 100
    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "a"), ser.values)
        assert np.shares_memory(get_array(result, "b"), ser2.values)

    result.iloc[0, 1] = 1000
    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "b"), ser2.values)
    tm.assert_series_equal(ser, ser_orig)
    tm.assert_series_equal(ser2, ser2_orig)


def test_concat_frames_chained(using_copy_on_write):
    df1 = DataFrame({"a": [1, 2, 3], "b": [0.1, 0.2, 0.3]})
    df2 = DataFrame({"c": [4, 5, 6]})
    df3 = DataFrame({"d": [4, 5, 6]})
    result = concat([concat([df1, df2], axis=1), df3], axis=1)
    expected = result.copy()

    if using_copy_on_write:
        assert np.shares_memory(get_array(result, "a"), get_array(df1, "a"))
        assert np.shares_memory(get_array(result, "c"), get_array(df2, "c"))
        assert np.shares_memory(get_array(result, "d"), get_array(df3, "d"))
    else:
        assert not np.shares_memory(get_array(result, "a"), get_array(df1, "a"))
        assert not np.shares_memory(get_array(result, "c"), get_array(df2, "c"))
        assert not np.shares_memory(get_array(result, "d"), get_array(df3, "d"))

    df1.iloc[0, 0] = 100
    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "a"), get_array(df1, "a"))

    tm.assert_frame_equal(result, expected)


def test_concat_series_chained(using_copy_on_write):
    ser1 = Series([1, 2, 3], name="a")
    ser2 = Series([4, 5, 6], name="c")
    ser3 = Series([4, 5, 6], name="d")
    result = concat([concat([ser1, ser2], axis=1), ser3], axis=1)
    expected = result.copy()

    if using_copy_on_write:
        assert np.shares_memory(get_array(result, "a"), get_array(ser1, "a"))
        assert np.shares_memory(get_array(result, "c"), get_array(ser2, "c"))
        assert np.shares_memory(get_array(result, "d"), get_array(ser3, "d"))
    else:
        assert not np.shares_memory(get_array(result, "a"), get_array(ser1, "a"))
        assert not np.shares_memory(get_array(result, "c"), get_array(ser2, "c"))
        assert not np.shares_memory(get_array(result, "d"), get_array(ser3, "d"))

    ser1.iloc[0] = 100
    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "a"), get_array(ser1, "a"))

    tm.assert_frame_equal(result, expected)


def test_concat_series_updating_input(using_copy_on_write):
    ser = Series([1, 2], name="a")
    ser2 = Series([3, 4], name="b")
    expected = DataFrame({"a": [1, 2], "b": [3, 4]})
    result = concat([ser, ser2], axis=1)

    if using_copy_on_write:
        assert np.shares_memory(get_array(result, "a"), get_array(ser, "a"))
        assert np.shares_memory(get_array(result, "b"), get_array(ser2, "b"))
    else:
        assert not np.shares_memory(get_array(result, "a"), get_array(ser, "a"))
        assert not np.shares_memory(get_array(result, "b"), get_array(ser2, "b"))

    ser.iloc[0] = 100
    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "a"), get_array(ser, "a"))
        assert np.shares_memory(get_array(result, "b"), get_array(ser2, "b"))
    tm.assert_frame_equal(result, expected)

    ser2.iloc[0] = 1000
    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "b"), get_array(ser2, "b"))
    tm.assert_frame_equal(result, expected)


def test_concat_mixed_series_frame(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "c": 1})
    ser = Series([4, 5, 6], name="d")
    result = concat([df, ser], axis=1)
    expected = result.copy()

    if using_copy_on_write:
        assert np.shares_memory(get_array(result, "a"), get_array(df, "a"))
        assert np.shares_memory(get_array(result, "c"), get_array(df, "c"))
        assert np.shares_memory(get_array(result, "d"), get_array(ser, "d"))
    else:
        assert not np.shares_memory(get_array(result, "a"), get_array(df, "a"))
        assert not np.shares_memory(get_array(result, "c"), get_array(df, "c"))
        assert not np.shares_memory(get_array(result, "d"), get_array(ser, "d"))

    ser.iloc[0] = 100
    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "d"), get_array(ser, "d"))

    df.iloc[0, 0] = 100
    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "a"), get_array(df, "a"))
    tm.assert_frame_equal(result, expected)
