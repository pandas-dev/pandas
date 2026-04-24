import numpy as np
import pytest

from pandas import (
    NA,
    DataFrame,
    Interval,
    NaT,
    Series,
    Timestamp,
    interval_range,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


@pytest.mark.parametrize("method", ["pad", "nearest", "linear"])
def test_interpolate_no_op(method):
    df = DataFrame({"a": [1, 2]})
    df_orig = df.copy()

    if method == "pad":
        msg = f"Can not interpolate with method={method}"
        with pytest.raises(ValueError, match=msg):
            df.interpolate(method=method)
    else:
        result = df.interpolate(method=method)
        assert np.shares_memory(get_array(result, "a"), get_array(df, "a"))
        assert result.index is not df.index
        assert result.columns is not df.columns

        result.iloc[0, 0] = 100

        assert not np.shares_memory(get_array(result, "a"), get_array(df, "a"))
        tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("func", ["ffill", "bfill"])
def test_interp_fill_functions(func):
    # Check that these takes the same code paths as interpolate
    df = DataFrame({"a": [1, 2]})
    df_orig = df.copy()

    result = getattr(df, func)()

    assert np.shares_memory(get_array(result, "a"), get_array(df, "a"))
    assert result.index is not df.index
    assert result.columns is not df.columns

    result.iloc[0, 0] = 100
    assert not np.shares_memory(get_array(result, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("func", ["ffill", "bfill"])
@pytest.mark.parametrize(
    "vals", [[1, np.nan, 2], [Timestamp("2019-12-31"), NaT, Timestamp("2020-12-31")]]
)
def test_interpolate_triggers_copy(vals, func):
    df = DataFrame({"a": vals})
    result = getattr(df, func)()

    assert not np.shares_memory(get_array(result, "a"), get_array(df, "a"))
    # Check that we don't have references when triggering a copy
    assert result._mgr._has_no_reference(0)


@pytest.mark.parametrize(
    "vals", [[1, np.nan, 2], [Timestamp("2019-12-31"), NaT, Timestamp("2020-12-31")]]
)
def test_interpolate_inplace_no_reference_no_copy(vals):
    df = DataFrame({"a": vals})
    arr = get_array(df, "a")
    df.interpolate(method="linear", inplace=True)

    assert np.shares_memory(arr, get_array(df, "a"))
    # Check that we don't have references when triggering a copy
    assert df._mgr._has_no_reference(0)


@pytest.mark.parametrize(
    "vals", [[1, np.nan, 2], [Timestamp("2019-12-31"), NaT, Timestamp("2020-12-31")]]
)
def test_interpolate_inplace_with_refs(vals):
    df = DataFrame({"a": [1, np.nan, 2]})
    df_orig = df.copy()
    arr = get_array(df, "a")
    view = df[:]
    df.interpolate(method="linear", inplace=True)
    # Check that copy was triggered in interpolate and that we don't
    # have any references left
    assert not np.shares_memory(arr, get_array(df, "a"))
    tm.assert_frame_equal(df_orig, view)
    assert df._mgr._has_no_reference(0)
    assert view._mgr._has_no_reference(0)


@pytest.mark.parametrize("func", ["ffill", "bfill"])
@pytest.mark.parametrize("dtype", ["float64", "Float64"])
def test_interp_fill_functions_inplace(func, dtype):
    # Check that these takes the same code paths as interpolate
    df = DataFrame({"a": [1, np.nan, 2]}, dtype=dtype)
    df_orig = df.copy()
    arr = get_array(df, "a")
    view = df[:]

    getattr(df, func)(inplace=True)

    # Check that copy was triggered in interpolate and that we don't
    # have any references left
    assert not np.shares_memory(arr, get_array(df, "a"))
    tm.assert_frame_equal(df_orig, view)
    assert df._mgr._has_no_reference(0)
    assert view._mgr._has_no_reference(0)


def test_interpolate_cannot_with_object_dtype():
    df = DataFrame({"a": ["a", np.nan, "c"], "b": 1})
    df["a"] = df["a"].astype(object)

    msg = "DataFrame cannot interpolate with object dtype"
    with pytest.raises(TypeError, match=msg):
        df.interpolate()


def test_interpolate_object_convert_no_op():
    df = DataFrame({"a": ["a", "b", "c"], "b": 1})
    df["a"] = df["a"].astype(object)
    arr_a = get_array(df, "a")

    # Now CoW makes a copy, it should not!
    assert df._mgr._has_no_reference(0)
    assert np.shares_memory(arr_a, get_array(df, "a"))


def test_interpolate_object_convert_copies():
    df = DataFrame({"a": [1, np.nan, 2.5], "b": 1})
    arr_a = get_array(df, "a")
    msg = "Can not interpolate with method=pad"
    with pytest.raises(ValueError, match=msg):
        df.interpolate(method="pad", inplace=True)

    assert df._mgr._has_no_reference(0)
    assert np.shares_memory(arr_a, get_array(df, "a"))


def test_interpolate_downcast_reference_triggers_copy():
    df = DataFrame({"a": [1, np.nan, 2.5], "b": 1})
    df_orig = df.copy()
    arr_a = get_array(df, "a")
    view = df[:]

    msg = "Can not interpolate with method=pad"
    with pytest.raises(ValueError, match=msg):
        df.interpolate(method="pad", inplace=True)
        assert df._mgr._has_no_reference(0)
        assert not np.shares_memory(arr_a, get_array(df, "a"))

    tm.assert_frame_equal(df_orig, view)


def test_fillna():
    df = DataFrame({"a": [1.5, np.nan], "b": 1})
    df_orig = df.copy()

    df2 = df.fillna(5.5)
    assert np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
    assert df2.index is not df.index
    assert df2.columns is not df.columns

    df2.iloc[0, 1] = 100
    tm.assert_frame_equal(df_orig, df)


def test_fillna_dict():
    df = DataFrame({"a": [1.5, np.nan], "b": 1})
    df_orig = df.copy()

    df2 = df.fillna({"a": 100.5})
    assert np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
    assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    df2.iloc[0, 1] = 100
    tm.assert_frame_equal(df_orig, df)


def test_fillna_inplace():
    df = DataFrame({"a": [1.5, np.nan], "b": 1})
    arr_a = get_array(df, "a")
    arr_b = get_array(df, "b")

    df.fillna(5.5, inplace=True)
    assert np.shares_memory(get_array(df, "a"), arr_a)
    assert np.shares_memory(get_array(df, "b"), arr_b)
    assert df._mgr._has_no_reference(0)
    assert df._mgr._has_no_reference(1)


def test_fillna_inplace_reference():
    df = DataFrame({"a": [1.5, np.nan], "b": 1})
    df_orig = df.copy()
    arr_a = get_array(df, "a")
    arr_b = get_array(df, "b")
    view = df[:]

    df.fillna(5.5, inplace=True)
    assert not np.shares_memory(get_array(df, "a"), arr_a)
    assert np.shares_memory(get_array(df, "b"), arr_b)
    assert view._mgr._has_no_reference(0)
    assert df._mgr._has_no_reference(0)
    tm.assert_frame_equal(view, df_orig)
    expected = DataFrame({"a": [1.5, 5.5], "b": 1})
    tm.assert_frame_equal(df, expected)


def test_fillna_interval_inplace_reference():
    # Set dtype explicitly to avoid implicit cast when setting nan
    ser = Series(
        interval_range(start=0, end=5), name="a", dtype="interval[float64, right]"
    )
    ser.iloc[1] = np.nan

    ser_orig = ser.copy()
    view = ser[:]
    ser.fillna(value=Interval(left=0, right=5), inplace=True)

    assert not np.shares_memory(
        get_array(ser, "a").left.values, get_array(view, "a").left.values
    )
    tm.assert_series_equal(view, ser_orig)


def test_fillna_series_empty_arg():
    ser = Series([1, np.nan, 2])
    ser_orig = ser.copy()
    result = ser.fillna({})
    assert np.shares_memory(get_array(ser), get_array(result))

    ser.iloc[0] = 100.5
    tm.assert_series_equal(ser_orig, result)


def test_fillna_series_empty_arg_inplace():
    ser = Series([1, np.nan, 2])
    arr = get_array(ser)
    ser.fillna({}, inplace=True)

    assert np.shares_memory(get_array(ser), arr)
    assert ser._mgr._has_no_reference(0)


def test_fillna_ea_noop_shares_memory(any_numeric_ea_and_arrow_dtype):
    df = DataFrame({"a": [1, NA, 3], "b": 1}, dtype=any_numeric_ea_and_arrow_dtype)
    df_orig = df.copy()
    df2 = df.fillna(100)

    assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    assert np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
    assert not df2._mgr._has_no_reference(1)
    tm.assert_frame_equal(df_orig, df)

    df2.iloc[0, 1] = 100
    assert not np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
    assert df2._mgr._has_no_reference(1)
    assert df._mgr._has_no_reference(1)
    tm.assert_frame_equal(df_orig, df)


def test_fillna_inplace_ea_noop_shares_memory(any_numeric_ea_and_arrow_dtype):
    df = DataFrame({"a": [1, NA, 3], "b": 1}, dtype=any_numeric_ea_and_arrow_dtype)
    df_orig = df.copy()
    view = df[:]
    df.fillna(100, inplace=True)
    assert not np.shares_memory(get_array(df, "a"), get_array(view, "a"))

    assert np.shares_memory(get_array(df, "b"), get_array(view, "b"))
    assert not df._mgr._has_no_reference(1)
    assert not view._mgr._has_no_reference(1)

    df.iloc[0, 1] = 100
    tm.assert_frame_equal(df_orig, view)


def test_fillna_chained_assignment():
    df = DataFrame({"a": [1, np.nan, 2], "b": 1})
    df_orig = df.copy()
    with tm.raises_chained_assignment_error():
        df["a"].fillna(100, inplace=True)
    tm.assert_frame_equal(df, df_orig)

    with tm.raises_chained_assignment_error():
        df[["a"]].fillna(100, inplace=True)
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("func", ["interpolate", "ffill", "bfill"])
def test_interpolate_chained_assignment(func):
    df = DataFrame({"a": [1, np.nan, 2], "b": 1})
    df_orig = df.copy()
    with tm.raises_chained_assignment_error():
        getattr(df["a"], func)(inplace=True)
    tm.assert_frame_equal(df, df_orig)

    with tm.raises_chained_assignment_error():
        getattr(df[["a"]], func)(inplace=True)
    tm.assert_frame_equal(df, df_orig)
