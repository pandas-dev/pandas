import numpy as np
import pytest

import pandas as pd
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
def test_replace(replace_kwargs):
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()

    df_replaced = df.replace(**replace_kwargs)

    if (df_replaced["b"] == df["b"]).all():
        assert np.shares_memory(get_array(df_replaced, "b"), get_array(df, "b"))
    assert tm.shares_memory(get_array(df_replaced, "c"), get_array(df, "c"))

    # mutating squeezed df triggers a copy-on-write for that column/block
    df_replaced.loc[0, "c"] = -1
    assert not np.shares_memory(get_array(df_replaced, "c"), get_array(df, "c"))

    if "a" in replace_kwargs["to_replace"]:
        arr = get_array(df_replaced, "a")
        df_replaced.loc[0, "a"] = 100
        assert np.shares_memory(get_array(df_replaced, "a"), arr)
    tm.assert_frame_equal(df, df_orig)


def test_replace_regex_inplace_refs():
    df = pd.DataFrame({"a": ["aaa", "bbb"]})
    df_orig = df.copy()
    view = df[:]
    arr = get_array(df, "a")
    df.replace(to_replace=r"^a.*$", value="new", inplace=True, regex=True)
    assert not np.shares_memory(arr, get_array(df, "a"))
    assert df._mgr._has_no_reference(0)
    tm.assert_frame_equal(view, df_orig)


def test_replace_regex_inplace():
    df = pd.DataFrame({"a": ["aaa", "bbb"]})
    arr = get_array(df, "a")
    df.replace(to_replace=r"^a.*$", value="new", inplace=True, regex=True)
    assert df._mgr._has_no_reference(0)
    assert tm.shares_memory(arr, get_array(df, "a"))

    df_orig = df.copy()
    df2 = df.replace(to_replace=r"^b.*$", value="new", regex=True)
    tm.assert_frame_equal(df_orig, df)
    assert not tm.shares_memory(get_array(df2, "a"), get_array(df, "a"))


@pytest.mark.parametrize("dtype", ["object", "string[python]"])
def test_replace_regex_non_inplace_does_not_mutate_original(dtype):
    # GH#57733 - replace_regex with non-string value triggers astype to object,
    # which for StringArray (backed by object ndarray) returned a view, causing
    # in-place mutation of the original
    df = pd.DataFrame({"b": list("ab..")}, dtype=dtype)
    df_orig = df.copy()
    result = df.replace([r"\s*\.\s*", "b"], 0, regex=True)

    expected = pd.DataFrame({"b": ["a", 0, 0, 0]}, dtype=object)
    tm.assert_frame_equal(result, expected)
    # Original must be unchanged
    tm.assert_frame_equal(df, df_orig)


def test_replace_regex_inplace_no_op():
    df = pd.DataFrame({"a": [1, 2]})
    arr = get_array(df, "a")
    df.replace(to_replace=r"^a.$", value="new", inplace=True, regex=True)
    assert df._mgr._has_no_reference(0)
    assert np.shares_memory(arr, get_array(df, "a"))

    df_orig = df.copy()
    df2 = df.replace(to_replace=r"^x.$", value="new", regex=True)
    tm.assert_frame_equal(df_orig, df)
    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))


def test_replace_mask_all_false_second_block():
    df = pd.DataFrame({"a": [1.5, 2, 3], "b": 100.5, "c": 1, "d": 2})
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
    df = pd.DataFrame({"a": [1.5, 2, 3], "b": 100.5})
    df_orig = df.copy()

    df2 = df.replace(to_replace=1.5, value="a")
    assert np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
    assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    df2.loc[0, "b"] = 0.5
    tm.assert_frame_equal(df, df_orig)  # Original is unchanged
    assert not np.shares_memory(get_array(df, "b"), get_array(df2, "b"))


def test_replace_to_replace_wrong_dtype():
    df = pd.DataFrame({"a": [1.5, 2, 3], "b": 100.5})
    df_orig = df.copy()

    df2 = df.replace(to_replace="xxx", value=1.5)

    assert np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
    assert np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    df2.loc[0, "b"] = 0.5
    tm.assert_frame_equal(df, df_orig)  # Original is unchanged
    assert not np.shares_memory(get_array(df, "b"), get_array(df2, "b"))


def test_replace_list_categorical():
    df = pd.DataFrame({"a": ["a", "b", "c"]}, dtype="category")
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
    df = pd.DataFrame({"a": ["a", "b", "c"]}, dtype="category")
    view = df[:]
    df_orig = df.copy()
    df.replace(["c"], value="a", inplace=True)
    tm.assert_frame_equal(df_orig, view)


@pytest.mark.parametrize("to_replace", [1.5, [1.5], []])
def test_replace_inplace(to_replace):
    df = pd.DataFrame({"a": [1.5, 2, 3]})
    arr_a = get_array(df, "a")
    df.replace(to_replace=1.5, value=15.5, inplace=True)

    assert np.shares_memory(get_array(df, "a"), arr_a)
    assert df._mgr._has_no_reference(0)


@pytest.mark.parametrize("to_replace", [1.5, [1.5]])
def test_replace_inplace_reference(to_replace):
    df = pd.DataFrame({"a": [1.5, 2, 3]})
    arr_a = get_array(df, "a")
    view = df[:]
    df.replace(to_replace=to_replace, value=15.5, inplace=True)

    assert not np.shares_memory(get_array(df, "a"), arr_a)
    assert df._mgr._has_no_reference(0)
    assert view._mgr._has_no_reference(0)


@pytest.mark.parametrize("to_replace", ["a", 100.5])
def test_replace_inplace_reference_no_op(to_replace):
    df = pd.DataFrame({"a": [1.5, 2, 3]})
    arr_a = get_array(df, "a")
    view = df[:]
    df.replace(to_replace=to_replace, value=15.5, inplace=True)

    assert np.shares_memory(get_array(df, "a"), arr_a)
    assert not df._mgr._has_no_reference(0)
    assert not view._mgr._has_no_reference(0)


@pytest.mark.parametrize("to_replace", [1, [1]])
def test_replace_categorical_inplace_reference(to_replace):
    df = pd.DataFrame({"a": pd.Categorical([1, 2, 3])})
    df_orig = df.copy()
    arr_a = get_array(df, "a")
    view = df[:]
    df.replace(to_replace=to_replace, value=1, inplace=True)
    assert not np.shares_memory(get_array(df, "a").codes, arr_a.codes)
    assert df._mgr._has_no_reference(0)
    assert view._mgr._has_no_reference(0)
    tm.assert_frame_equal(view, df_orig)


def test_replace_categorical_inplace():
    df = pd.DataFrame({"a": pd.Categorical([1, 2, 3])})
    arr_a = get_array(df, "a")
    df.replace(to_replace=1, value=1, inplace=True)

    assert np.shares_memory(get_array(df, "a").codes, arr_a.codes)
    assert df._mgr._has_no_reference(0)

    expected = pd.DataFrame({"a": pd.Categorical([1, 2, 3])})
    tm.assert_frame_equal(df, expected)


def test_replace_categorical():
    df = pd.DataFrame({"a": pd.Categorical([1, 2, 3])})
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
    df = pd.DataFrame({"a": [1.5, 2, 3]})
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
    df = pd.DataFrame({"a": [1, 2]})

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
    df = pd.DataFrame({"a": ["a", "b", "c"]}, dtype=object)
    arr = get_array(df, "a")
    df.replace(["c"], value, inplace=True)
    assert np.shares_memory(arr, get_array(df, "a"))
    assert df._mgr._has_no_reference(0)


def test_replace_list_multiple_elements_inplace():
    df = pd.DataFrame({"a": [1, 2, 3]})
    arr = get_array(df, "a")
    df.replace([1, 2], 4, inplace=True)
    assert np.shares_memory(arr, get_array(df, "a"))
    assert df._mgr._has_no_reference(0)


def test_replace_list_none():
    df = pd.DataFrame({"a": ["a", "b", "c"]})

    df_orig = df.copy()
    df2 = df.replace(["b"], value=None)
    tm.assert_frame_equal(df, df_orig)

    assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    # replace multiple values that don't actually replace anything with None
    # https://github.com/pandas-dev/pandas/issues/59770
    df3 = df.replace(["d", "e", "f"], value=None)
    tm.assert_frame_equal(df3, df_orig)
    assert tm.shares_memory(get_array(df, "a"), get_array(df3, "a"))


def test_replace_list_none_inplace_refs():
    df = pd.DataFrame({"a": ["a", "b", "c"]})
    arr = get_array(df, "a")
    df_orig = df.copy()
    view = df[:]
    df.replace(["a"], value=None, inplace=True)
    assert df._mgr._has_no_reference(0)
    assert not np.shares_memory(arr, get_array(df, "a"))
    tm.assert_frame_equal(df_orig, view)


def test_replace_columnwise_no_op_inplace():
    df = pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]})
    view = df[:]
    df_orig = df.copy()
    df.replace({"a": 10}, 100, inplace=True)
    assert np.shares_memory(get_array(view, "a"), get_array(df, "a"))
    df.iloc[0, 0] = 100
    tm.assert_frame_equal(view, df_orig)


def test_replace_columnwise_no_op():
    df = pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]})
    df_orig = df.copy()
    df2 = df.replace({"a": 10}, 100)
    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    df2.iloc[0, 0] = 100
    tm.assert_frame_equal(df, df_orig)


def test_replace_chained_assignment():
    df = pd.DataFrame({"a": [1, np.nan, 2], "b": 1})
    df_orig = df.copy()
    with tm.raises_chained_assignment_error():
        df["a"].replace(1, 100, inplace=True)
    tm.assert_frame_equal(df, df_orig)

    with tm.raises_chained_assignment_error():
        df[["a"]].replace(1, 100, inplace=True)
    tm.assert_frame_equal(df, df_orig)


def test_replace_listlike():
    df = pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]})
    df_orig = df.copy()

    result = df.replace([200, 201], [11, 11])
    assert np.shares_memory(get_array(result, "a"), get_array(df, "a"))

    result.iloc[0, 0] = 100
    tm.assert_frame_equal(df, df)

    result = df.replace([200, 2], [10, 10])
    assert not np.shares_memory(get_array(df, "a"), get_array(result, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_replace_listlike_inplace():
    df = pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]})
    arr = get_array(df, "a")
    df.replace([200, 2], [10, 11], inplace=True)
    assert np.shares_memory(get_array(df, "a"), arr)

    view = df[:]
    df_orig = df.copy()
    df.replace([200, 3], [10, 11], inplace=True)
    assert not np.shares_memory(get_array(df, "a"), arr)
    tm.assert_frame_equal(view, df_orig)


@pytest.mark.parametrize(
    "to_replace, value",
    [
        # matching pair before a no-op pair
        ([2, 200], [20, 10]),
        # every pair matches
        ([2, 3], [20, 30]),
        # more than two pairs
        ([2, 200, 201], [20, 10, 10]),
    ],
)
def test_replace_listlike_inplace_keeps_tracking_copies(to_replace, value):
    # GH#58966 replace_list dropped the block's self-reference once it had
    #  copied, so views taken afterwards were no longer copy-on-write protected
    df = pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]})
    view = df[:]
    df.replace(to_replace, value, inplace=True)
    # the copy triggered by ``view`` must not have leaked the original values
    tm.assert_frame_equal(view, pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]}))

    df_orig = df.copy()
    view2 = df[:]
    df.iloc[0, 0] = 999
    assert not np.shares_memory(get_array(df, "a"), get_array(view2, "a"))
    tm.assert_frame_equal(view2, df_orig)


def test_replace_dictlike_inplace_keeps_tracking_copies():
    # GH#58966
    df = pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]})
    view = df[:]
    df.replace({2: 20, 3: 30}, inplace=True)
    tm.assert_frame_equal(view, pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]}))

    df_orig = df.copy()
    view2 = df[:]
    df.iloc[0, 0] = 999
    assert not np.shares_memory(get_array(df, "a"), get_array(view2, "a"))
    tm.assert_frame_equal(view2, df_orig)


def test_replace_object_inplace_keeps_tracking_copies():
    # GH#58966
    ser = pd.Series(["a", "b", "c"], dtype=object)
    view = ser[:]
    ser.replace(["b", "z"], ["B", "Z"], inplace=True)
    tm.assert_series_equal(view, pd.Series(["a", "b", "c"], dtype=object))

    ser_orig = ser.copy()
    view2 = ser[:]
    ser.iloc[0] = "x"
    assert not np.shares_memory(get_array(ser), get_array(view2))
    tm.assert_series_equal(view2, ser_orig)


def test_replace_listlike_inplace_upcast_keeps_tracking_copies():
    # GH#58966 no pre-existing reference is needed: the upcast to object mints
    #  fresh refs, which were then emptied
    ser = pd.Series([1, 2, 3])
    ser.replace([2, 200], ["x", "y"], inplace=True)

    ser_orig = ser.copy()
    view = ser[:]
    ser.iloc[0] = "Q"
    assert not np.shares_memory(get_array(ser), get_array(view))
    tm.assert_series_equal(view, ser_orig)


def test_replace_regex_inplace_keeps_tracking_copies():
    # GH#58966
    ser = pd.Series(["ab", "bc", "cd"], dtype=object)
    view = ser[:]
    ser.replace(["^a", "^z"], ["A", "Z"], regex=True, inplace=True)
    tm.assert_series_equal(view, pd.Series(["ab", "bc", "cd"], dtype=object))

    ser_orig = ser.copy()
    view2 = ser[:]
    ser.iloc[0] = "Q"
    assert not np.shares_memory(get_array(ser), get_array(view2))
    tm.assert_series_equal(view2, ser_orig)


@pytest.mark.parametrize("dtype", ["Int64", "datetime64[ns]"])
def test_replace_listlike_inplace_ea_keeps_tracking_copies(dtype):
    # GH#58966 EABackedBlock, not just NumpyBlock
    to_replace = [2, 200]
    value = [20, 10]
    new_val = 999
    if dtype == "datetime64[ns]":
        to_replace = [pd.Timestamp(val) for val in to_replace]
        value = [pd.Timestamp(val) for val in value]
        new_val = pd.Timestamp(new_val)

    ser = pd.Series([1, 2, 3]).astype(dtype)
    orig = ser.copy()
    view = ser[:]
    ser.replace(to_replace, value, inplace=True)
    tm.assert_series_equal(view, orig)

    ser_orig = ser.copy()
    view2 = ser[:]
    ser.iloc[0] = new_val
    assert not np.shares_memory(get_array(ser), get_array(view2))
    tm.assert_series_equal(view2, ser_orig)


@pytest.mark.parametrize(
    "to_replace, value",
    [
        ([1, 2], ["s0", "s1"]),
        ([1, 2, 3], ["s0", "s1", "s2"]),
    ],
)
def test_replace_listlike_inplace_split_keeps_tracking_copies(to_replace, value):
    # GH#58966 upcasting one column of a multi-column block splits it; the
    #  untouched sibling survives the pair loop and must stay registered
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    df.replace(to_replace, value, inplace=True)

    df_orig = df.copy()
    view = df[:]
    df.iloc[0, 1] = 999
    assert not np.shares_memory(get_array(df, "b"), get_array(view, "b"))
    tm.assert_frame_equal(view, df_orig)


def test_replace_dictlike_inplace_split_keeps_tracking_copies():
    # GH#58966
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    df.replace({1: "s0", 2: "s1"}, inplace=True)

    df_orig = df.copy()
    view = df[:]
    df.iloc[0, 1] = 999
    assert not np.shares_memory(get_array(df, "b"), get_array(view, "b"))
    tm.assert_frame_equal(view, df_orig)
