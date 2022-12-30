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


def test_concat_series(using_copy_on_write):
    ser = Series([1, 2], name="a")
    ser2 = Series([3, 4], name="b")
    ser_orig = ser.copy()
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
