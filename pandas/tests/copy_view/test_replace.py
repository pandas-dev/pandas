import numpy as np
import pytest

from pandas import (
    Categorical,
    DataFrame,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


@pytest.mark.parametrize(
    "replace_kwargs",
    [
        {"to_replace": {"a": 1, "b": 4}, "value": -1},
        # Test CoW splits blocks to avoid copying unchanged columns
        {"to_replace": {"a": 1}, "value": -1},
        {"to_replace": {"b": 4}, "value": -1},
        {"to_replace": {"b": {4: 1}}},
        # TODO: Add these in a further optimization
        # We would need to see which columns got replaced in the mask
        # which could be expensive
        # {"to_replace": {"b": 1}},
        # 1
    ],
)
def test_replace(using_copy_on_write, replace_kwargs):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": ["foo", "bar", "baz"]})
    df_orig = df.copy()

    df_replaced = df.replace(**replace_kwargs)

    if using_copy_on_write:
        if (df_replaced["b"] == df["b"]).all():
            assert np.shares_memory(get_array(df_replaced, "b"), get_array(df, "b"))
        assert np.shares_memory(get_array(df_replaced, "c"), get_array(df, "c"))

    # mutating squeezed df triggers a copy-on-write for that column/block
    df_replaced.loc[0, "c"] = -1
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df_replaced, "c"), get_array(df, "c"))

    if "a" in replace_kwargs["to_replace"]:
        arr = get_array(df_replaced, "a")
        df_replaced.loc[0, "a"] = 100
        assert np.shares_memory(get_array(df_replaced, "a"), arr)
    tm.assert_frame_equal(df, df_orig)


def test_replace_mask_all_false_second_block(using_copy_on_write):
    df = DataFrame({"a": [1.5, 2, 3], "b": 100.5, "c": 1, "d": 2})
    df_orig = df.copy()

    df2 = df.replace(to_replace=1.5, value=55.5)

    if using_copy_on_write:
        # TODO: Block splitting would allow us to avoid copying b
        assert np.shares_memory(get_array(df, "c"), get_array(df2, "c"))
        assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    else:
        assert not np.shares_memory(get_array(df, "c"), get_array(df2, "c"))
        assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    df2.loc[0, "c"] = 1
    tm.assert_frame_equal(df, df_orig)  # Original is unchanged

    if using_copy_on_write:
        assert not np.shares_memory(get_array(df, "c"), get_array(df2, "c"))
        # TODO: This should split and not copy the whole block
        # assert np.shares_memory(get_array(df, "d"), get_array(df2, "d"))


def test_replace_coerce_single_column(using_copy_on_write, using_array_manager):
    df = DataFrame({"a": [1.5, 2, 3], "b": 100.5})
    df_orig = df.copy()

    df2 = df.replace(to_replace=1.5, value="a")

    if using_copy_on_write:
        assert np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
        assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    elif not using_array_manager:
        assert np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
        assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    if using_copy_on_write:
        df2.loc[0, "b"] = 0.5
        tm.assert_frame_equal(df, df_orig)  # Original is unchanged
        assert not np.shares_memory(get_array(df, "b"), get_array(df2, "b"))


def test_replace_to_replace_wrong_dtype(using_copy_on_write):
    df = DataFrame({"a": [1.5, 2, 3], "b": 100.5})
    df_orig = df.copy()

    df2 = df.replace(to_replace="xxx", value=1.5)

    if using_copy_on_write:
        assert np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
        assert np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    else:
        assert not np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
        assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    df2.loc[0, "b"] = 0.5
    tm.assert_frame_equal(df, df_orig)  # Original is unchanged

    if using_copy_on_write:
        assert not np.shares_memory(get_array(df, "b"), get_array(df2, "b"))


def test_replace_inplace(using_copy_on_write):
    df = DataFrame({"a": [1.5, 2, 3]})
    arr_a = get_array(df, "a")
    df.replace(to_replace=1.5, value=15.5, inplace=True)

    assert np.shares_memory(get_array(df, "a"), arr_a)
    if using_copy_on_write:
        assert df._mgr._has_no_reference(0)


@pytest.mark.parametrize("to_replace", [1.5, [1.5]])
def test_replace_inplace_reference(using_copy_on_write, to_replace):
    df = DataFrame({"a": [1.5, 2, 3]})
    arr_a = get_array(df, "a")
    view = df[:]
    df.replace(to_replace=to_replace, value=15.5, inplace=True)

    if using_copy_on_write:
        assert not np.shares_memory(get_array(df, "a"), arr_a)
        assert df._mgr._has_no_reference(0)
        assert view._mgr._has_no_reference(0)
    else:
        assert np.shares_memory(get_array(df, "a"), arr_a)


@pytest.mark.parametrize("to_replace", ["a", 100.5])
def test_replace_inplace_reference_no_op(using_copy_on_write, to_replace):
    df = DataFrame({"a": [1.5, 2, 3]})
    arr_a = get_array(df, "a")
    view = df[:]
    df.replace(to_replace=to_replace, value=15.5, inplace=True)

    assert np.shares_memory(get_array(df, "a"), arr_a)
    if using_copy_on_write:
        assert not df._mgr._has_no_reference(0)
        assert not view._mgr._has_no_reference(0)


@pytest.mark.parametrize("to_replace", [1, [1]])
@pytest.mark.parametrize("val", [1, 1.5])
def test_replace_categorical_inplace_reference(using_copy_on_write, val, to_replace):
    df = DataFrame({"a": Categorical([1, 2, 3])})
    df_orig = df.copy()
    arr_a = get_array(df, "a")
    view = df[:]
    df.replace(to_replace=to_replace, value=val, inplace=True)

    if using_copy_on_write:
        assert not np.shares_memory(get_array(df, "a").codes, arr_a.codes)
        assert df._mgr._has_no_reference(0)
        assert view._mgr._has_no_reference(0)
        tm.assert_frame_equal(view, df_orig)
    else:
        assert np.shares_memory(get_array(df, "a").codes, arr_a.codes)


@pytest.mark.parametrize("val", [1, 1.5])
def test_replace_categorical_inplace(using_copy_on_write, val):
    df = DataFrame({"a": Categorical([1, 2, 3])})
    arr_a = get_array(df, "a")
    df.replace(to_replace=1, value=val, inplace=True)

    assert np.shares_memory(get_array(df, "a").codes, arr_a.codes)
    if using_copy_on_write:
        assert df._mgr._has_no_reference(0)

    expected = DataFrame({"a": Categorical([val, 2, 3])})
    tm.assert_frame_equal(df, expected)


@pytest.mark.parametrize("val", [1, 1.5])
def test_replace_categorical(using_copy_on_write, val):
    df = DataFrame({"a": Categorical([1, 2, 3])})
    df_orig = df.copy()
    df2 = df.replace(to_replace=1, value=val)

    if using_copy_on_write:
        assert df._mgr._has_no_reference(0)
        assert df2._mgr._has_no_reference(0)
    assert not np.shares_memory(get_array(df, "a").codes, get_array(df2, "a").codes)
    tm.assert_frame_equal(df, df_orig)

    arr_a = get_array(df2, "a").codes
    df2.iloc[0, 0] = 2.0
    assert np.shares_memory(get_array(df2, "a").codes, arr_a)


@pytest.mark.parametrize("method", ["where", "mask"])
def test_masking_inplace(using_copy_on_write, method):
    df = DataFrame({"a": [1.5, 2, 3]})
    df_orig = df.copy()
    arr_a = get_array(df, "a")
    view = df[:]

    method = getattr(df, method)
    method(df["a"] > 1.6, -1, inplace=True)

    if using_copy_on_write:
        assert not np.shares_memory(get_array(df, "a"), arr_a)
        assert df._mgr._has_no_reference(0)
        assert view._mgr._has_no_reference(0)
        tm.assert_frame_equal(view, df_orig)
    else:
        assert np.shares_memory(get_array(df, "a"), arr_a)
