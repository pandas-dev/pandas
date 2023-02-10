import numpy as np
import pytest

from pandas import (
    DataFrame,
    NaT,
    Series,
    Timestamp,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


@pytest.mark.parametrize("method", ["pad", "nearest", "linear"])
def test_interpolate_no_op(using_copy_on_write, method):
    df = DataFrame({"a": [1, 2]})
    df_orig = df.copy()

    result = df.interpolate(method=method)

    if using_copy_on_write:
        assert np.shares_memory(get_array(result, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(result, "a"), get_array(df, "a"))

    result.iloc[0, 0] = 100

    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("func", ["ffill", "bfill"])
def test_interp_fill_functions(using_copy_on_write, func):
    # Check that these takes the same code paths as interpolate
    df = DataFrame({"a": [1, 2]})
    df_orig = df.copy()

    result = getattr(df, func)()

    if using_copy_on_write:
        assert np.shares_memory(get_array(result, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(result, "a"), get_array(df, "a"))

    result.iloc[0, 0] = 100

    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("func", ["ffill", "bfill"])
@pytest.mark.parametrize(
    "vals", [[1, np.nan, 2], [Timestamp("2019-12-31"), NaT, Timestamp("2020-12-31")]]
)
def test_interpolate_triggers_copy(using_copy_on_write, vals, func):
    df = DataFrame({"a": vals})
    result = getattr(df, func)()

    assert not np.shares_memory(get_array(result, "a"), get_array(df, "a"))
    if using_copy_on_write:
        # Check that we don't have references when triggering a copy
        assert result._mgr._has_no_reference(0)


@pytest.mark.parametrize(
    "vals", [[1, np.nan, 2], [Timestamp("2019-12-31"), NaT, Timestamp("2020-12-31")]]
)
def test_interpolate_inplace_no_reference_no_copy(using_copy_on_write, vals):
    df = DataFrame({"a": vals})
    arr = get_array(df, "a")
    df.interpolate(method="linear", inplace=True)

    assert np.shares_memory(arr, get_array(df, "a"))
    if using_copy_on_write:
        # Check that we don't have references when triggering a copy
        assert df._mgr._has_no_reference(0)


@pytest.mark.parametrize(
    "vals", [[1, np.nan, 2], [Timestamp("2019-12-31"), NaT, Timestamp("2020-12-31")]]
)
def test_interpolate_inplace_with_refs(using_copy_on_write, vals):
    df = DataFrame({"a": [1, np.nan, 2]})
    df_orig = df.copy()
    arr = get_array(df, "a")
    view = df[:]
    df.interpolate(method="linear", inplace=True)

    if using_copy_on_write:
        # Check that copy was triggered in interpolate and that we don't
        # have any references left
        assert not np.shares_memory(arr, get_array(df, "a"))
        tm.assert_frame_equal(df_orig, view)
        assert df._mgr._has_no_reference(0)
        assert view._mgr._has_no_reference(0)
    else:
        assert np.shares_memory(arr, get_array(df, "a"))


def test_interpolate_cleaned_fill_method(using_copy_on_write):
    # Check that "method is set to None" case works correctly
    df = DataFrame({"a": ["a", np.nan, "c"], "b": 1})
    df_orig = df.copy()

    result = df.interpolate(method="asfreq")

    if using_copy_on_write:
        assert np.shares_memory(get_array(result, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(result, "a"), get_array(df, "a"))

    result.iloc[0, 0] = Timestamp("2021-12-31")

    if using_copy_on_write:
        assert not np.shares_memory(get_array(result, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_interpolate_object_convert_no_op(using_copy_on_write):
    df = DataFrame({"a": ["a", "b", "c"], "b": 1})
    arr_a = get_array(df, "a")
    df.interpolate(method="pad", inplace=True)

    # Now CoW makes a copy, it should not!
    if using_copy_on_write:
        assert df._mgr._has_no_reference(0)
        assert np.shares_memory(arr_a, get_array(df, "a"))


def test_interpolate_object_convert_copies(using_copy_on_write):
    df = DataFrame({"a": Series([1, 2], dtype=object), "b": 1})
    arr_a = get_array(df, "a")
    df.interpolate(method="pad", inplace=True)

    if using_copy_on_write:
        assert df._mgr._has_no_reference(0)
        assert not np.shares_memory(arr_a, get_array(df, "a"))


def test_interpolate_downcast(using_copy_on_write):
    df = DataFrame({"a": [1, np.nan, 2.5], "b": 1})
    arr_a = get_array(df, "a")
    df.interpolate(method="pad", inplace=True, downcast="infer")

    if using_copy_on_write:
        assert df._mgr._has_no_reference(0)
    assert np.shares_memory(arr_a, get_array(df, "a"))


def test_interpolate_downcast_reference_triggers_copy(using_copy_on_write):
    df = DataFrame({"a": [1, np.nan, 2.5], "b": 1})
    df_orig = df.copy()
    arr_a = get_array(df, "a")
    view = df[:]
    df.interpolate(method="pad", inplace=True, downcast="infer")

    if using_copy_on_write:
        assert df._mgr._has_no_reference(0)
        assert not np.shares_memory(arr_a, get_array(df, "a"))
        tm.assert_frame_equal(df_orig, view)
    else:
        tm.assert_frame_equal(df, view)
