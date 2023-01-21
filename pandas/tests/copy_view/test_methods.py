import numpy as np
import pytest

from pandas import (
    DataFrame,
    Index,
    MultiIndex,
    Period,
    Series,
    Timestamp,
    date_range,
    period_range,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


def test_copy(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_copy = df.copy()

    # the deep copy doesn't share memory
    assert not np.shares_memory(get_array(df_copy, "a"), get_array(df, "a"))
    if using_copy_on_write:
        assert df_copy._mgr.refs is None

    # mutating copy doesn't mutate original
    df_copy.iloc[0, 0] = 0
    assert df.iloc[0, 0] == 1


def test_copy_shallow(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_copy = df.copy(deep=False)

    # the shallow copy still shares memory
    assert np.shares_memory(get_array(df_copy, "a"), get_array(df, "a"))
    if using_copy_on_write:
        assert df_copy._mgr.refs is not None

    if using_copy_on_write:
        # mutating shallow copy doesn't mutate original
        df_copy.iloc[0, 0] = 0
        assert df.iloc[0, 0] == 1
        # mutating triggered a copy-on-write -> no longer shares memory
        assert not np.shares_memory(get_array(df_copy, "a"), get_array(df, "a"))
        # but still shares memory for the other columns/blocks
        assert np.shares_memory(get_array(df_copy, "c"), get_array(df, "c"))
    else:
        # mutating shallow copy does mutate original
        df_copy.iloc[0, 0] = 0
        assert df.iloc[0, 0] == 0
        # and still shares memory
        assert np.shares_memory(get_array(df_copy, "a"), get_array(df, "a"))


@pytest.mark.parametrize("copy", [True, None, False])
@pytest.mark.parametrize(
    "method",
    [
        lambda df, copy: df.rename(columns=str.lower, copy=copy),
        lambda df, copy: df.reindex(columns=["a", "c"], copy=copy),
        lambda df, copy: df.reindex_like(df, copy=copy),
        lambda df, copy: df.set_axis(["a", "b", "c"], axis="index", copy=copy),
        lambda df, copy: df.rename_axis(index="test", copy=copy),
        lambda df, copy: df.rename_axis(columns="test", copy=copy),
        # lambda df, copy: df.astype({'b': 'int64'}, copy=copy),
        # lambda df, copy: df.swaplevel(0, 0, copy=copy),
        lambda df, copy: df.swapaxes(0, 0, copy=copy),
        lambda df, copy: df.truncate(0, 5, copy=copy),
        # lambda df, copy: df.infer_objects(copy=copy)
        lambda df, copy: df.to_timestamp(copy=copy),
        lambda df, copy: df.to_period(freq="D", copy=copy),
        lambda df, copy: df.tz_localize("US/Central", copy=copy),
        lambda df, copy: df.tz_convert("US/Central", copy=copy),
        lambda df, copy: df.set_flags(allows_duplicate_labels=False, copy=copy),
    ],
    ids=[
        "rename",
        "reindex",
        "reindex_like",
        "set_axis",
        "rename_axis0",
        "rename_axis1",
        # "astype",  # CoW not yet implemented
        # "swaplevel",  # only series
        "swapaxes",
        "truncate",
        # "infer_objects",  # CoW not yet implemented
        "to_timestamp",
        "to_period",
        "tz_localize",
        "tz_convert",
        "set_flags",
    ],
)
def test_methods_copy_keyword(
    request, method, copy, using_copy_on_write, using_array_manager
):
    index = None
    if "to_timestamp" in request.node.callspec.id:
        index = period_range("2012-01-01", freq="D", periods=3)
    elif "to_period" in request.node.callspec.id:
        index = date_range("2012-01-01", freq="D", periods=3)
    elif "tz_localize" in request.node.callspec.id:
        index = date_range("2012-01-01", freq="D", periods=3)
    elif "tz_convert" in request.node.callspec.id:
        index = date_range("2012-01-01", freq="D", periods=3, tz="Europe/Brussels")

    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]}, index=index)
    df2 = method(df, copy=copy)

    share_memory = (using_copy_on_write and copy is not True) or copy is False

    if request.node.callspec.id.startswith("reindex-"):
        # TODO copy=False without CoW still returns a copy in this case
        if not using_copy_on_write and not using_array_manager and copy is False:
            share_memory = False
        # TODO copy=True with CoW still returns a view
        if using_copy_on_write:
            share_memory = True

    if share_memory:
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))


# -----------------------------------------------------------------------------
# DataFrame methods returning new DataFrame using shallow copy


def test_reset_index(using_copy_on_write):
    # Case: resetting the index (i.e. adding a new column) + mutating the
    # resulting dataframe
    df = DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]}, index=[10, 11, 12]
    )
    df_orig = df.copy()
    df2 = df.reset_index()
    df2._mgr._verify_integrity()

    if using_copy_on_write:
        # still shares memory (df2 is a shallow copy)
        assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
        assert np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    # mutating df2 triggers a copy-on-write for that column / block
    df2.iloc[0, 2] = 0
    assert not np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    tm.assert_frame_equal(df, df_orig)


def test_rename_columns(using_copy_on_write):
    # Case: renaming columns returns a new dataframe
    # + afterwards modifying the result
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.rename(columns=str.upper)

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "A"), get_array(df, "a"))
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "A"), get_array(df, "a"))
    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "C"), get_array(df, "c"))
    expected = DataFrame({"A": [0, 2, 3], "B": [4, 5, 6], "C": [0.1, 0.2, 0.3]})
    tm.assert_frame_equal(df2, expected)
    tm.assert_frame_equal(df, df_orig)


def test_rename_columns_modify_parent(using_copy_on_write):
    # Case: renaming columns returns a new dataframe
    # + afterwards modifying the original (parent) dataframe
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df2 = df.rename(columns=str.upper)
    df2_orig = df2.copy()

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "A"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "A"), get_array(df, "a"))
    df.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "A"), get_array(df, "a"))
    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "C"), get_array(df, "c"))
    expected = DataFrame({"a": [0, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    tm.assert_frame_equal(df, expected)
    tm.assert_frame_equal(df2, df2_orig)


def test_pipe(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "b": 1.5})
    df_orig = df.copy()

    def testfunc(df):
        return df

    df2 = df.pipe(testfunc)

    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    # mutating df2 triggers a copy-on-write for that column
    df2.iloc[0, 0] = 0
    if using_copy_on_write:
        tm.assert_frame_equal(df, df_orig)
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        expected = DataFrame({"a": [0, 2, 3], "b": 1.5})
        tm.assert_frame_equal(df, expected)

        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))


def test_pipe_modify_df(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "b": 1.5})
    df_orig = df.copy()

    def testfunc(df):
        df.iloc[0, 0] = 100
        return df

    df2 = df.pipe(testfunc)

    assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))

    if using_copy_on_write:
        tm.assert_frame_equal(df, df_orig)
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        expected = DataFrame({"a": [100, 2, 3], "b": 1.5})
        tm.assert_frame_equal(df, expected)

        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))


def test_reindex_columns(using_copy_on_write):
    # Case: reindexing the column returns a new dataframe
    # + afterwards modifying the result
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.reindex(columns=["a", "c"])

    if using_copy_on_write:
        # still shares memory (df2 is a shallow copy)
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    # mutating df2 triggers a copy-on-write for that column
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    tm.assert_frame_equal(df, df_orig)


def test_drop_on_column(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.drop(columns="a")
    df2._mgr._verify_integrity()

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
        assert np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    else:
        assert not np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
        assert not np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    tm.assert_frame_equal(df, df_orig)


def test_select_dtypes(using_copy_on_write):
    # Case: selecting columns using `select_dtypes()` returns a new dataframe
    # + afterwards modifying the result
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.select_dtypes("int64")
    df2._mgr._verify_integrity()

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    # mutating df2 triggers a copy-on-write for that column/block
    df2.iloc[0, 0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "filter_kwargs", [{"items": ["a"]}, {"like": "a"}, {"regex": "a"}]
)
def test_filter(using_copy_on_write, filter_kwargs):
    # Case: selecting columns using `filter()` returns a new dataframe
    # + afterwards modifying the result
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.filter(**filter_kwargs)
    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    # mutating df2 triggers a copy-on-write for that column/block
    if using_copy_on_write:
        df2.iloc[0, 0] = 0
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_shift_no_op(using_copy_on_write):
    df = DataFrame(
        [[1, 2], [3, 4], [5, 6]],
        index=date_range("2020-01-01", "2020-01-03"),
        columns=["a", "b"],
    )
    df_orig = df.copy()
    df2 = df.shift(periods=0)

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    df.iloc[0, 0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
    tm.assert_frame_equal(df2, df_orig)


def test_shift_index(using_copy_on_write):
    df = DataFrame(
        [[1, 2], [3, 4], [5, 6]],
        index=date_range("2020-01-01", "2020-01-03"),
        columns=["a", "b"],
    )
    df2 = df.shift(periods=1, axis=0)

    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))


def test_shift_rows_freq(using_copy_on_write):
    df = DataFrame(
        [[1, 2], [3, 4], [5, 6]],
        index=date_range("2020-01-01", "2020-01-03"),
        columns=["a", "b"],
    )
    df_orig = df.copy()
    df_orig.index = date_range("2020-01-02", "2020-01-04")
    df2 = df.shift(periods=1, freq="1D")

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    df.iloc[0, 0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))
    tm.assert_frame_equal(df2, df_orig)


def test_shift_columns(using_copy_on_write):
    df = DataFrame(
        [[1, 2], [3, 4], [5, 6]], columns=date_range("2020-01-01", "2020-01-02")
    )
    df2 = df.shift(periods=1, axis=1)

    assert np.shares_memory(get_array(df2, "2020-01-02"), get_array(df, "2020-01-01"))
    df.iloc[0, 1] = 0
    if using_copy_on_write:
        assert not np.shares_memory(
            get_array(df2, "2020-01-02"), get_array(df, "2020-01-01")
        )
    expected = DataFrame(
        [[np.nan, 1], [np.nan, 3], [np.nan, 5]],
        columns=date_range("2020-01-01", "2020-01-02"),
    )
    tm.assert_frame_equal(df2, expected)


def test_pop(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    view_original = df[:]
    result = df.pop("a")

    assert np.shares_memory(result.values, get_array(view_original, "a"))
    assert np.shares_memory(get_array(df, "b"), get_array(view_original, "b"))

    if using_copy_on_write:
        result.iloc[0] = 0
        assert not np.shares_memory(result.values, get_array(view_original, "a"))
    df.iloc[0, 0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df, "b"), get_array(view_original, "b"))
        tm.assert_frame_equal(view_original, df_orig)
    else:
        expected = DataFrame({"a": [1, 2, 3], "b": [0, 5, 6], "c": [0.1, 0.2, 0.3]})
        tm.assert_frame_equal(view_original, expected)


@pytest.mark.parametrize(
    "func",
    [
        lambda x, y: x.align(y),
        lambda x, y: x.align(y.a, axis=0),
        lambda x, y: x.align(y.a.iloc[slice(0, 1)], axis=1),
    ],
)
def test_align_frame(using_copy_on_write, func):
    df = DataFrame({"a": [1, 2, 3], "b": "a"})
    df_orig = df.copy()
    df_changed = df[["b", "a"]].copy()
    df2, _ = func(df, df_changed)

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    df2.iloc[0, 0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_align_series(using_copy_on_write):
    ser = Series([1, 2])
    ser_orig = ser.copy()
    ser_other = ser.copy()
    ser2, ser_other_result = ser.align(ser_other)

    if using_copy_on_write:
        assert np.shares_memory(ser2.values, ser.values)
        assert np.shares_memory(ser_other_result.values, ser_other.values)
    else:
        assert not np.shares_memory(ser2.values, ser.values)
        assert not np.shares_memory(ser_other_result.values, ser_other.values)

    ser2.iloc[0] = 0
    ser_other_result.iloc[0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(ser2.values, ser.values)
        assert not np.shares_memory(ser_other_result.values, ser_other.values)
    tm.assert_series_equal(ser, ser_orig)
    tm.assert_series_equal(ser_other, ser_orig)


def test_to_frame(using_copy_on_write):
    # Case: converting a Series to a DataFrame with to_frame
    ser = Series([1, 2, 3])
    ser_orig = ser.copy()

    df = ser[:].to_frame()

    # currently this always returns a "view"
    assert np.shares_memory(ser.values, get_array(df, 0))

    df.iloc[0, 0] = 0

    if using_copy_on_write:
        # mutating df triggers a copy-on-write for that column
        assert not np.shares_memory(ser.values, get_array(df, 0))
        tm.assert_series_equal(ser, ser_orig)
    else:
        # but currently select_dtypes() actually returns a view -> mutates parent
        expected = ser_orig.copy()
        expected.iloc[0] = 0
        tm.assert_series_equal(ser, expected)

    # modify original series -> don't modify dataframe
    df = ser[:].to_frame()
    ser.iloc[0] = 0

    if using_copy_on_write:
        tm.assert_frame_equal(df, ser_orig.to_frame())
    else:
        expected = ser_orig.copy().to_frame()
        expected.iloc[0, 0] = 0
        tm.assert_frame_equal(df, expected)


@pytest.mark.parametrize("ax", ["index", "columns"])
def test_swapaxes_noop(using_copy_on_write, ax):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    df_orig = df.copy()
    df2 = df.swapaxes(ax, ax)

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    # mutating df2 triggers a copy-on-write for that column/block
    df2.iloc[0, 0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_swapaxes_single_block(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]}, index=["x", "y", "z"])
    df_orig = df.copy()
    df2 = df.swapaxes("index", "columns")

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "x"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "x"), get_array(df, "a"))

    # mutating df2 triggers a copy-on-write for that column/block
    df2.iloc[0, 0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df2, "x"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "method, idx",
    [
        (lambda df: df.copy(deep=False).copy(deep=False), 0),
        (lambda df: df.reset_index().reset_index(), 2),
        (lambda df: df.rename(columns=str.upper).rename(columns=str.lower), 0),
        (lambda df: df.copy(deep=False).select_dtypes(include="number"), 0),
    ],
    ids=["shallow-copy", "reset_index", "rename", "select_dtypes"],
)
def test_chained_methods(request, method, idx, using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()

    # when not using CoW, only the copy() variant actually gives a view
    df2_is_view = not using_copy_on_write and request.node.callspec.id == "shallow-copy"

    # modify df2 -> don't modify df
    df2 = method(df)
    df2.iloc[0, idx] = 0
    if not df2_is_view:
        tm.assert_frame_equal(df, df_orig)

    # modify df -> don't modify df2
    df2 = method(df)
    df.iloc[0, 0] = 0
    if not df2_is_view:
        tm.assert_frame_equal(df2.iloc[:, idx:], df_orig)


@pytest.mark.parametrize("obj", [Series([1, 2], name="a"), DataFrame({"a": [1, 2]})])
def test_to_timestamp(using_copy_on_write, obj):
    obj.index = Index([Period("2012-1-1", freq="D"), Period("2012-1-2", freq="D")])

    obj_orig = obj.copy()
    obj2 = obj.to_timestamp()

    if using_copy_on_write:
        assert np.shares_memory(get_array(obj2, "a"), get_array(obj, "a"))
    else:
        assert not np.shares_memory(get_array(obj2, "a"), get_array(obj, "a"))

    # mutating obj2 triggers a copy-on-write for that column / block
    obj2.iloc[0] = 0
    assert not np.shares_memory(get_array(obj2, "a"), get_array(obj, "a"))
    tm.assert_equal(obj, obj_orig)


@pytest.mark.parametrize("obj", [Series([1, 2], name="a"), DataFrame({"a": [1, 2]})])
def test_to_period(using_copy_on_write, obj):
    obj.index = Index([Timestamp("2019-12-31"), Timestamp("2020-12-31")])

    obj_orig = obj.copy()
    obj2 = obj.to_period(freq="Y")

    if using_copy_on_write:
        assert np.shares_memory(get_array(obj2, "a"), get_array(obj, "a"))
    else:
        assert not np.shares_memory(get_array(obj2, "a"), get_array(obj, "a"))

    # mutating obj2 triggers a copy-on-write for that column / block
    obj2.iloc[0] = 0
    assert not np.shares_memory(get_array(obj2, "a"), get_array(obj, "a"))
    tm.assert_equal(obj, obj_orig)


def test_set_index(using_copy_on_write):
    # GH 49473
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.set_index("a")

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
    else:
        assert not np.shares_memory(get_array(df2, "b"), get_array(df, "b"))

    # mutating df2 triggers a copy-on-write for that column / block
    df2.iloc[0, 1] = 0
    assert not np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    tm.assert_frame_equal(df, df_orig)


def test_add_prefix(using_copy_on_write):
    # GH 49473
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.add_prefix("CoW_")

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "CoW_a"), get_array(df, "a"))
    df2.iloc[0, 0] = 0

    assert not np.shares_memory(get_array(df2, "CoW_a"), get_array(df, "a"))

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "CoW_c"), get_array(df, "c"))
    expected = DataFrame(
        {"CoW_a": [0, 2, 3], "CoW_b": [4, 5, 6], "CoW_c": [0.1, 0.2, 0.3]}
    )
    tm.assert_frame_equal(df2, expected)
    tm.assert_frame_equal(df, df_orig)


def test_add_suffix(using_copy_on_write):
    # GH 49473
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.add_suffix("_CoW")
    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "a_CoW"), get_array(df, "a"))
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "a_CoW"), get_array(df, "a"))
    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "c_CoW"), get_array(df, "c"))
    expected = DataFrame(
        {"a_CoW": [0, 2, 3], "b_CoW": [4, 5, 6], "c_CoW": [0.1, 0.2, 0.3]}
    )
    tm.assert_frame_equal(df2, expected)
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("axis, val", [(0, 5.5), (1, np.nan)])
def test_dropna(using_copy_on_write, axis, val):
    df = DataFrame({"a": [1, 2, 3], "b": [4, val, 6], "c": "d"})
    df_orig = df.copy()
    df2 = df.dropna(axis=axis)

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    df2.iloc[0, 0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("val", [5, 5.5])
def test_dropna_series(using_copy_on_write, val):
    ser = Series([1, val, 4])
    ser_orig = ser.copy()
    ser2 = ser.dropna()

    if using_copy_on_write:
        assert np.shares_memory(ser2.values, ser.values)
    else:
        assert not np.shares_memory(ser2.values, ser.values)

    ser2.iloc[0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(ser2.values, ser.values)
    tm.assert_series_equal(ser, ser_orig)


@pytest.mark.parametrize(
    "method",
    [
        lambda df: df.head(),
        lambda df: df.head(2),
        lambda df: df.tail(),
        lambda df: df.tail(3),
    ],
)
def test_head_tail(method, using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "b": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = method(df)
    df2._mgr._verify_integrity()

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
        assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))

    # modify df2 to trigger CoW for that block
    df2.iloc[0, 0] = 0
    assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        # without CoW enabled, head and tail return views. Mutating df2 also mutates df.
        df2.iloc[0, 0] = 1
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"before": "a", "after": "b", "axis": 1},
        {"before": 0, "after": 1, "axis": 0},
    ],
)
def test_truncate(using_copy_on_write, kwargs):
    df = DataFrame({"a": [1, 2, 3], "b": 1, "c": 2})
    df_orig = df.copy()
    df2 = df.truncate(**kwargs)
    df2._mgr._verify_integrity()

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    df2.iloc[0, 0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("method", ["assign", "drop_duplicates"])
def test_assign_drop_duplicates(using_copy_on_write, method):
    df = DataFrame({"a": [1, 2, 3]})
    df_orig = df.copy()
    df2 = getattr(df, method)()
    df2._mgr._verify_integrity()

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    df2.iloc[0, 0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("obj", [Series([1, 2]), DataFrame({"a": [1, 2]})])
def test_take(using_copy_on_write, obj):
    # Check that no copy is made when we take all rows in original order
    obj_orig = obj.copy()
    obj2 = obj.take([0, 1])

    if using_copy_on_write:
        assert np.shares_memory(obj2.values, obj.values)
    else:
        assert not np.shares_memory(obj2.values, obj.values)

    obj2.iloc[0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(obj2.values, obj.values)
    tm.assert_equal(obj, obj_orig)


@pytest.mark.parametrize("obj", [Series([1, 2]), DataFrame({"a": [1, 2]})])
def test_between_time(using_copy_on_write, obj):
    obj.index = date_range("2018-04-09", periods=2, freq="1D20min")
    obj_orig = obj.copy()
    obj2 = obj.between_time("0:00", "1:00")

    if using_copy_on_write:
        assert np.shares_memory(obj2.values, obj.values)
    else:
        assert not np.shares_memory(obj2.values, obj.values)

    obj2.iloc[0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(obj2.values, obj.values)
    tm.assert_equal(obj, obj_orig)


def test_reindex_like(using_copy_on_write):
    df = DataFrame({"a": [1, 2], "b": "a"})
    other = DataFrame({"b": "a", "a": [1, 2]})

    df_orig = df.copy()
    df2 = df.reindex_like(other)

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    df2.iloc[0, 1] = 0
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_sort_index(using_copy_on_write):
    # GH 49473
    ser = Series([1, 2, 3])
    ser_orig = ser.copy()
    ser2 = ser.sort_index()

    if using_copy_on_write:
        assert np.shares_memory(ser.values, ser2.values)
    else:
        assert not np.shares_memory(ser.values, ser2.values)

    # mutating ser triggers a copy-on-write for the column / block
    ser2.iloc[0] = 0
    assert not np.shares_memory(ser2.values, ser.values)
    tm.assert_series_equal(ser, ser_orig)


@pytest.mark.parametrize(
    "obj, kwargs",
    [(Series([1, 2, 3], name="a"), {}), (DataFrame({"a": [1, 2, 3]}), {"by": "a"})],
)
def test_sort_values(using_copy_on_write, obj, kwargs):
    obj_orig = obj.copy()
    obj2 = obj.sort_values(**kwargs)

    if using_copy_on_write:
        assert np.shares_memory(get_array(obj2, "a"), get_array(obj, "a"))
    else:
        assert not np.shares_memory(get_array(obj2, "a"), get_array(obj, "a"))

    # mutating df triggers a copy-on-write for the column / block
    obj2.iloc[0] = 0
    assert not np.shares_memory(get_array(obj2, "a"), get_array(obj, "a"))
    tm.assert_equal(obj, obj_orig)


@pytest.mark.parametrize(
    "obj, kwargs",
    [(Series([1, 2, 3], name="a"), {}), (DataFrame({"a": [1, 2, 3]}), {"by": "a"})],
)
def test_sort_values_inplace(using_copy_on_write, obj, kwargs, using_array_manager):
    obj_orig = obj.copy()
    view = obj[:]
    obj.sort_values(inplace=True, **kwargs)

    assert np.shares_memory(get_array(obj, "a"), get_array(view, "a"))

    # mutating obj triggers a copy-on-write for the column / block
    obj.iloc[0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(get_array(obj, "a"), get_array(view, "a"))
        tm.assert_equal(view, obj_orig)
    else:
        assert np.shares_memory(get_array(obj, "a"), get_array(view, "a"))


def test_round(using_copy_on_write):
    df = DataFrame({"a": [1, 2], "b": "c"})
    df2 = df.round()
    df_orig = df.copy()

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "b"), get_array(df, "b"))

    df2.iloc[0, 1] = "d"
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
    tm.assert_frame_equal(df, df_orig)


def test_reorder_levels(using_copy_on_write):
    index = MultiIndex.from_tuples(
        [(1, 1), (1, 2), (2, 1), (2, 2)], names=["one", "two"]
    )
    df = DataFrame({"a": [1, 2, 3, 4]}, index=index)
    df_orig = df.copy()
    df2 = df.reorder_levels(order=["two", "one"])

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    df2.iloc[0, 0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_series_reorder_levels(using_copy_on_write):
    index = MultiIndex.from_tuples(
        [(1, 1), (1, 2), (2, 1), (2, 2)], names=["one", "two"]
    )
    ser = Series([1, 2, 3, 4], index=index)
    ser_orig = ser.copy()
    ser2 = ser.reorder_levels(order=["two", "one"])

    if using_copy_on_write:
        assert np.shares_memory(ser2.values, ser.values)
    else:
        assert not np.shares_memory(ser2.values, ser.values)

    ser2.iloc[0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(ser2.values, ser.values)
    tm.assert_series_equal(ser, ser_orig)


@pytest.mark.parametrize("obj", [Series([1, 2, 3]), DataFrame({"a": [1, 2, 3]})])
def test_swaplevel(using_copy_on_write, obj):
    index = MultiIndex.from_tuples([(1, 1), (1, 2), (2, 1)], names=["one", "two"])
    obj.index = index
    obj_orig = obj.copy()
    obj2 = obj.swaplevel()

    if using_copy_on_write:
        assert np.shares_memory(obj2.values, obj.values)
    else:
        assert not np.shares_memory(obj2.values, obj.values)

    obj2.iloc[0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(obj2.values, obj.values)
    tm.assert_equal(obj, obj_orig)


def test_frame_set_axis(using_copy_on_write):
    # GH 49473
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.set_axis(["a", "b", "c"], axis="index")

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    # mutating df2 triggers a copy-on-write for that column / block
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_series_set_axis(using_copy_on_write):
    # GH 49473
    ser = Series([1, 2, 3])
    ser_orig = ser.copy()
    ser2 = ser.set_axis(["a", "b", "c"], axis="index")

    if using_copy_on_write:
        assert np.shares_memory(ser, ser2)
    else:
        assert not np.shares_memory(ser, ser2)

    # mutating ser triggers a copy-on-write for the column / block
    ser2.iloc[0] = 0
    assert not np.shares_memory(ser2, ser)
    tm.assert_series_equal(ser, ser_orig)


def test_set_flags(using_copy_on_write):
    ser = Series([1, 2, 3])
    ser_orig = ser.copy()
    ser2 = ser.set_flags(allows_duplicate_labels=False)

    assert np.shares_memory(ser, ser2)

    # mutating ser triggers a copy-on-write for the column / block
    ser2.iloc[0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(ser2, ser)
        tm.assert_series_equal(ser, ser_orig)
    else:
        assert np.shares_memory(ser2, ser)
        expected = Series([0, 2, 3])
        tm.assert_series_equal(ser, expected)


@pytest.mark.parametrize("copy_kwargs", [{"copy": True}, {}])
@pytest.mark.parametrize("kwargs", [{"mapper": "test"}, {"index": "test"}])
def test_rename_axis(using_copy_on_write, kwargs, copy_kwargs):
    df = DataFrame({"a": [1, 2, 3, 4]}, index=Index([1, 2, 3, 4], name="a"))
    df_orig = df.copy()
    df2 = df.rename_axis(**kwargs, **copy_kwargs)

    if using_copy_on_write and not copy_kwargs:
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    df2.iloc[0, 0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "func, tz", [("tz_convert", "Europe/Berlin"), ("tz_localize", None)]
)
def test_tz_convert_localize(using_copy_on_write, func, tz):
    # GH 49473
    ser = Series(
        [1, 2], index=date_range(start="2014-08-01 09:00", freq="H", periods=2, tz=tz)
    )
    ser_orig = ser.copy()
    ser2 = getattr(ser, func)("US/Central")

    if using_copy_on_write:
        assert np.shares_memory(ser.values, ser2.values)
    else:
        assert not np.shares_memory(ser.values, ser2.values)

    # mutating ser triggers a copy-on-write for the column / block
    ser2.iloc[0] = 0
    assert not np.shares_memory(ser2.values, ser.values)
    tm.assert_series_equal(ser, ser_orig)


def test_droplevel(using_copy_on_write):
    # GH 49473
    index = MultiIndex.from_tuples([(1, 1), (1, 2), (2, 1)], names=["one", "two"])
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}, index=index)
    df_orig = df.copy()
    df2 = df.droplevel(0)

    if using_copy_on_write:
        assert np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    else:
        assert not np.shares_memory(get_array(df2, "c"), get_array(df, "c"))

    # mutating df2 triggers a copy-on-write for that column / block
    df2.iloc[0, 0] = 0

    assert not np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    tm.assert_frame_equal(df, df_orig)


def test_squeeze(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3]})
    df_orig = df.copy()
    series = df.squeeze()

    # Should share memory regardless of CoW since squeeze is just an iloc
    assert np.shares_memory(series.values, get_array(df, "a"))

    # mutating squeezed df triggers a copy-on-write for that column/block
    series.iloc[0] = 0
    if using_copy_on_write:
        assert not np.shares_memory(series.values, get_array(df, "a"))
        tm.assert_frame_equal(df, df_orig)
    else:
        # Without CoW the original will be modified
        assert np.shares_memory(series.values, get_array(df, "a"))
        assert df.loc[0, "a"] == 0


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


def test_putmask(using_copy_on_write):
    df = DataFrame({"a": [1, 2], "b": 1, "c": 2})
    view = df[:]
    df_orig = df.copy()
    df[df == df] = 5

    if using_copy_on_write:
        assert not np.shares_memory(get_array(view, "a"), get_array(df, "a"))
        tm.assert_frame_equal(view, df_orig)
    else:
        # Without CoW the original will be modified
        assert np.shares_memory(get_array(view, "a"), get_array(df, "a"))
        assert view.iloc[0, 0] == 5


def test_isetitem(using_copy_on_write):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})
    df_orig = df.copy()
    df2 = df.copy(deep=None)  # Trigger a CoW
    df2.isetitem(1, np.array([-1, -2, -3]))  # This is inplace

    if using_copy_on_write:
        assert np.shares_memory(get_array(df, "c"), get_array(df2, "c"))
        assert np.shares_memory(get_array(df, "a"), get_array(df2, "a"))
    else:
        assert not np.shares_memory(get_array(df, "c"), get_array(df2, "c"))
        assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    df2.loc[0, "a"] = 0
    tm.assert_frame_equal(df, df_orig)  # Original is unchanged

    if using_copy_on_write:
        assert np.shares_memory(get_array(df, "c"), get_array(df2, "c"))
    else:
        assert not np.shares_memory(get_array(df, "c"), get_array(df2, "c"))
