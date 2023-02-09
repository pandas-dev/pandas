import numpy as np
import pytest

from pandas import (
    DataFrame,
    Interval,
    Series,
    interval_range,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


def test_fillna(using_copy_on_write):
    df = DataFrame({"a": [1.5, np.nan], "b": 1})
    df_orig = df.copy()

    df2 = df.fillna(5.5)
    if using_copy_on_write:
        assert np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
    else:
        assert not np.shares_memory(get_array(df, "b"), get_array(df2, "b"))

    df2.iloc[0, 1] = 100
    tm.assert_frame_equal(df_orig, df)


@pytest.mark.parametrize("downcast", [None, False])
def test_fillna_inplace(using_copy_on_write, downcast):
    df = DataFrame({"a": [1.5, np.nan], "b": 1})
    arr_a = get_array(df, "a")
    arr_b = get_array(df, "b")

    df.fillna(5.5, inplace=True, downcast=downcast)
    assert np.shares_memory(get_array(df, "a"), arr_a)
    assert np.shares_memory(get_array(df, "b"), arr_b)
    if using_copy_on_write:
        assert df._mgr._has_no_reference(0)
        assert df._mgr._has_no_reference(1)


def test_fillna_inplace_reference(using_copy_on_write):
    df = DataFrame({"a": [1.5, np.nan], "b": 1})
    df_orig = df.copy()
    arr_a = get_array(df, "a")
    arr_b = get_array(df, "b")
    view = df[:]

    df.fillna(5.5, inplace=True)
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df, "a"), arr_a)
        assert np.shares_memory(get_array(df, "b"), arr_b)
        assert view._mgr._has_no_reference(0)
        assert df._mgr._has_no_reference(0)
        tm.assert_frame_equal(view, df_orig)
    else:
        assert np.shares_memory(get_array(df, "a"), arr_a)
        assert np.shares_memory(get_array(df, "b"), arr_b)
    expected = DataFrame({"a": [1.5, 5.5], "b": 1})
    tm.assert_frame_equal(df, expected)


def test_fillna_interval_inplace_reference(using_copy_on_write):
    ser = Series(interval_range(start=0, end=5), name="a")
    ser.iloc[1] = np.nan

    ser_orig = ser.copy()
    view = ser[:]
    ser.fillna(value=Interval(left=0, right=5), inplace=True)

    if using_copy_on_write:
        assert not np.shares_memory(
            get_array(ser, "a").left.values, get_array(view, "a").left.values
        )
        tm.assert_series_equal(view, ser_orig)
    else:
        assert np.shares_memory(
            get_array(ser, "a").left.values, get_array(view, "a").left.values
        )
