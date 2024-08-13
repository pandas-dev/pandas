import numpy as np
import pytest

from pandas._config import using_string_dtype

from pandas.compat import HAS_PYARROW

from pandas import (
    Categorical,
    DataFrame,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


@pytest.mark.xfail(using_string_dtype(), reason="TODO(infer_string)")
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
def test_replace(replace_kwargs):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": ["foo", "bar", "baz"]})
    df_orig = df.copy()

    df_replaced = df.replace(**replace_kwargs)

    if (df_replaced["b"] == df["b"]).all():
        assert np.shares_memory(get_array(df_replaced, "b"), get_array(df, "b"))
    assert np.shares_memory(get_array(df_replaced, "c"), get_array(df, "c"))

    # mutating squeezed df triggers a copy-on-write for that column/block
    df_replaced.loc[0, "c"] = -1
    assert not np.shares_memory(get_array(df_replaced, "c"), get_array(df, "c"))

    if "a" in replace_kwargs["to_replace"]:
        arr = get_array(df_replaced, "a")
        df_replaced.loc[0, "a"] = 100
        assert np.shares_memory(get_array(df_replaced, "a"), arr)
    tm.assert_frame_equal(df, df_orig)


def test_replace_regex_inplace_refs():
    df = DataFrame({"a": ["aaa", "bbb"]})
    df_orig = df.copy()
    view = df[:]
    arr = get_array(df, "a")
    df.replace(to_replace=r"^a.*$", value="new", inplace=True, regex=True)
    assert not np.shares_memory(arr, get_array(df, "a"))
    assert df._mgr._has_no_reference(0)
    tm.assert_frame_equal(view, df_orig)


@pytest.mark.xfail(using_string_dtype() and HAS_PYARROW, reason="TODO(infer_string)")
def test_replace_regex_inplace():
    df = DataFrame({"a": ["aaa", "bbb"]})
    arr = get_array(df, "a")
    df.replace(to_replace=r"^a.*$", value="new", inplace=True, regex=True)
    assert df._mgr._has_no_reference(0)
    assert np.shares_memory(arr, get_array(df, "a"))

    df_orig = df.copy()
    df2 = df.replace(to_replace=r"^b.*$", value="new", regex=True)
    tm.assert_frame_equal(df_orig, df)
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))


def test_replace_regex_inplace_no_op():
    df = DataFrame({"a": [1, 2]})
    arr = get_array(df, "a")
    df.replace(to_replace=r"^a.$", value="new", inplace=True, regex=True)
    assert df._mgr._has_no_reference(0)
    assert np.shares_memory(arr, get_array(df, "a"))

    df_orig = df.copy()
    df2 = df.replace(to_replace=r"^x.$", value="new", regex=True)
    tm.assert_frame_equal(df_orig, df)
    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))


def test_replace_mask_all_false_second_block():
    df = DataFrame({"a": [1.5, 2, 3], "b": 100.5, "c": 1, "d": 2})
    df_orig = df.copy()

    df2 = df.replace(to_replace=1.5, value=55.5)

    # TODO: Block splitting would allow us to avoid copying b
    assert np.shares_memory(get_array(df, "c"), get_array(df2, "c"))
    assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    df2.loc[0, "c"] = 1
    tm.assert_frame_equal(df, df_orig)  # Original is unchanged

    assert not np.shares_memory(get_array(df, "c"), get_array(df2, "c"))
    assert np.shares_memory(get_array(df, "d"), get_array(df2, "d"))


def test_replace_coerce_single_column():
    df = DataFrame({"a": [1.5, 2, 3], "b": 100.5})
    df_orig = df.copy()

    df2 = df.replace(to_replace=1.5, value="a")
    assert np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
    assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    df2.loc[0, "b"] = 0.5
    tm.assert_frame_equal(df, df_orig)  # Original is unchanged
    assert not np.shares_memory(get_array(df, "b"), get_array(df2, "b"))


def test_replace_to_replace_wrong_dtype():
    df = DataFrame({"a": [1.5, 2, 3], "b": 100.5})
    df_orig = df.copy()

    df2 = df.replace(to_replace="xxx", value=1.5)

    assert np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
    assert np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    df2.loc[0, "b"] = 0.5
    tm.assert_frame_equal(df, df_orig)  # Original is unchanged
    assert not np.shares_memory(get_array(df, "b"), get_array(df2, "b"))


def test_replace_list_categorical():
    df = DataFrame({"a": ["a", "b", "c"]}, dtype="category")
    arr = get_array(df, "a")

    df.replace(["c"], value="a", inplace=True)
    assert np.shares_memory(arr.codes, get_array(df, "a").codes)
    assert df._mgr._has_no_reference(0)

    df_orig = df.copy()
    df.replace(["b"], value="a")
    df2 = df.apply(lambda x: x.cat.rename_categories({"b": "d"}))
    assert not np.shares_memory(arr.codes, get_array(df2, "a").codes)

    tm.assert_frame_equal(df, df_orig)


def test_replace_list_inplace_refs_categorical():
    df = DataFrame({"a": ["a", "b", "c"]}, dtype="category")
    view = df[:]
    df_orig = df.copy()
    df.replace(["c"], value="a", inplace=True)
    tm.assert_frame_equal(df_orig, view)


@pytest.mark.parametrize("to_replace", [1.5, [1.5], []])
def test_replace_inplace(to_replace):
    df = DataFrame({"a": [1.5, 2, 3]})
    arr_a = get_array(df, "a")
    df.replace(to_replace=1.5, value=15.5, inplace=True)

    assert np.shares_memory(get_array(df, "a"), arr_a)
    assert df._mgr._has_no_reference(0)


@pytest.mark.parametrize("to_replace", [1.5, [1.5]])
def test_replace_inplace_reference(to_replace):
    df = DataFrame({"a": [1.5, 2, 3]})
    arr_a = get_array(df, "a")
    view = df[:]
    df.replace(to_replace=to_replace, value=15.5, inplace=True)

    assert not np.shares_memory(get_array(df, "a"), arr_a)
    assert df._mgr._has_no_reference(0)
    assert view._mgr._has_no_reference(0)


@pytest.mark.parametrize("to_replace", ["a", 100.5])
def test_replace_inplace_reference_no_op(to_replace):
    df = DataFrame({"a": [1.5, 2, 3]})
    arr_a = get_array(df, "a")
    view = df[:]
    df.replace(to_replace=to_replace, value=15.5, inplace=True)

    assert np.shares_memory(get_array(df, "a"), arr_a)
    assert not df._mgr._has_no_reference(0)
    assert not view._mgr._has_no_reference(0)


@pytest.mark.parametrize("to_replace", [1, [1]])
def test_replace_categorical_inplace_reference(to_replace):
    df = DataFrame({"a": Categorical([1, 2, 3])})
    df_orig = df.copy()
    arr_a = get_array(df, "a")
    view = df[:]
    df.replace(to_replace=to_replace, value=1, inplace=True)
    assert not np.shares_memory(get_array(df, "a").codes, arr_a.codes)
    assert df._mgr._has_no_reference(0)
    assert view._mgr._has_no_reference(0)
    tm.assert_frame_equal(view, df_orig)


def test_replace_categorical_inplace():
    df = DataFrame({"a": Categorical([1, 2, 3])})
    arr_a = get_array(df, "a")
    df.replace(to_replace=1, value=1, inplace=True)

    assert np.shares_memory(get_array(df, "a").codes, arr_a.codes)
    assert df._mgr._has_no_reference(0)

    expected = DataFrame({"a": Categorical([1, 2, 3])})
    tm.assert_frame_equal(df, expected)


def test_replace_categorical():
    df = DataFrame({"a": Categorical([1, 2, 3])})
    df_orig = df.copy()
    df2 = df.replace(to_replace=1, value=1)

    assert df._mgr._has_no_reference(0)
    assert df2._mgr._has_no_reference(0)
    assert not np.shares_memory(get_array(df, "a").codes, get_array(df2, "a").codes)
    tm.assert_frame_equal(df, df_orig)

    arr_a = get_array(df2, "a").codes
    df2.iloc[0, 0] = 2.0
    assert np.shares_memory(get_array(df2, "a").codes, arr_a)


@pytest.mark.parametrize("method", ["where", "mask"])
def test_masking_inplace(method):
    df = DataFrame({"a": [1.5, 2, 3]})
    df_orig = df.copy()
    arr_a = get_array(df, "a")
    view = df[:]

    method = getattr(df, method)
    method(df["a"] > 1.6, -1, inplace=True)

    assert not np.shares_memory(get_array(df, "a"), arr_a)
    assert df._mgr._has_no_reference(0)
    assert view._mgr._has_no_reference(0)
    tm.assert_frame_equal(view, df_orig)


def test_replace_empty_list():
    df = DataFrame({"a": [1, 2]})

    df2 = df.replace([], [])
    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    assert not df._mgr._has_no_reference(0)
    arr_a = get_array(df, "a")
    df.replace([], [])
    assert np.shares_memory(get_array(df, "a"), arr_a)
    assert not df._mgr._has_no_reference(0)
    assert not df2._mgr._has_no_reference(0)


@pytest.mark.parametrize("value", ["d", None])
def test_replace_object_list_inplace(value):
    df = DataFrame({"a": ["a", "b", "c"]}, dtype=object)
    arr = get_array(df, "a")
    df.replace(["c"], value, inplace=True)
    assert np.shares_memory(arr, get_array(df, "a"))
    assert df._mgr._has_no_reference(0)


def test_replace_list_multiple_elements_inplace():
    df = DataFrame({"a": [1, 2, 3]})
    arr = get_array(df, "a")
    df.replace([1, 2], 4, inplace=True)
    assert np.shares_memory(arr, get_array(df, "a"))
    assert df._mgr._has_no_reference(0)


def test_replace_list_none():
    df = DataFrame({"a": ["a", "b", "c"]})

    df_orig = df.copy()
    df2 = df.replace(["b"], value=None)
    tm.assert_frame_equal(df, df_orig)

    assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))


def test_replace_list_none_inplace_refs():
    df = DataFrame({"a": ["a", "b", "c"]})
    arr = get_array(df, "a")
    df_orig = df.copy()
    view = df[:]
    df.replace(["a"], value=None, inplace=True)
    assert df._mgr._has_no_reference(0)
    assert not np.shares_memory(arr, get_array(df, "a"))
    tm.assert_frame_equal(df_orig, view)


def test_replace_columnwise_no_op_inplace():
    df = DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]})
    view = df[:]
    df_orig = df.copy()
    df.replace({"a": 10}, 100, inplace=True)
    assert np.shares_memory(get_array(view, "a"), get_array(df, "a"))
    df.iloc[0, 0] = 100
    tm.assert_frame_equal(view, df_orig)


def test_replace_columnwise_no_op():
    df = DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]})
    df_orig = df.copy()
    df2 = df.replace({"a": 10}, 100)
    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    df2.iloc[0, 0] = 100
    tm.assert_frame_equal(df, df_orig)


def test_replace_chained_assignment():
    df = DataFrame({"a": [1, np.nan, 2], "b": 1})
    df_orig = df.copy()
    with tm.raises_chained_assignment_error():
        df["a"].replace(1, 100, inplace=True)
    tm.assert_frame_equal(df, df_orig)

    with tm.raises_chained_assignment_error():
        df[["a"]].replace(1, 100, inplace=True)
    tm.assert_frame_equal(df, df_orig)


def test_replace_listlike():
    df = DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]})
    df_orig = df.copy()

    result = df.replace([200, 201], [11, 11])
    assert np.shares_memory(get_array(result, "a"), get_array(df, "a"))

    result.iloc[0, 0] = 100
    tm.assert_frame_equal(df, df)

    result = df.replace([200, 2], [10, 10])
    assert not np.shares_memory(get_array(df, "a"), get_array(result, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_replace_listlike_inplace():
    df = DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]})
    arr = get_array(df, "a")
    df.replace([200, 2], [10, 11], inplace=True)
    assert np.shares_memory(get_array(df, "a"), arr)

    view = df[:]
    df_orig = df.copy()
    df.replace([200, 3], [10, 11], inplace=True)
    assert not np.shares_memory(get_array(df, "a"), arr)
    tm.assert_frame_equal(view, df_orig)
