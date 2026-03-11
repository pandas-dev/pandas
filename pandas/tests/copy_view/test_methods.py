import numpy as np
import pytest

from pandas.compat import HAS_PYARROW

import pandas as pd
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
from pandas.util.version import Version


def test_copy():
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_copy = df.copy()

    # the deep copy by defaults takes a shallow copy of the Index
    assert df_copy.index is not df.index
    assert df_copy.columns is not df.columns
    assert df_copy.index.is_(df.index)
    assert df_copy.columns.is_(df.columns)

    # the deep copy doesn't share memory
    assert not np.shares_memory(get_array(df_copy, "a"), get_array(df, "a"))
    assert not df_copy._mgr.blocks[0].refs.has_reference()
    assert not df_copy._mgr.blocks[1].refs.has_reference()

    assert df_copy.index is not df.index
    assert df_copy.columns is not df.columns

    # mutating copy doesn't mutate original
    df_copy.iloc[0, 0] = 0
    assert df.iloc[0, 0] == 1


def test_copy_shallow():
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_copy = df.copy(deep=False)

    # the shallow copy also makes a shallow copy of the index
    assert df_copy.index is not df.index
    assert df_copy.columns is not df.columns
    assert df_copy.index.is_(df.index)
    assert df_copy.columns.is_(df.columns)

    # the shallow copy still shares memory
    assert np.shares_memory(get_array(df_copy, "a"), get_array(df, "a"))
    assert df_copy._mgr.blocks[0].refs.has_reference()
    assert df_copy._mgr.blocks[1].refs.has_reference()

    # mutating shallow copy doesn't mutate original
    df_copy.iloc[0, 0] = 0
    assert df.iloc[0, 0] == 1
    # mutating triggered a copy-on-write -> no longer shares memory
    assert not np.shares_memory(get_array(df_copy, "a"), get_array(df, "a"))
    # but still shares memory for the other columns/blocks
    assert np.shares_memory(get_array(df_copy, "c"), get_array(df, "c"))


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
@pytest.mark.parametrize("copy", [True, None, False])
@pytest.mark.parametrize(
    "method",
    [
        lambda df, copy: df.rename(columns=str.lower, copy=copy),
        lambda df, copy: df.reindex(columns=["a", "c"], copy=copy),
        lambda df, copy: df.reindex_like(df, copy=copy),
        lambda df, copy: df.align(df, copy=copy)[0],
        lambda df, copy: df.set_axis(["a", "b", "c"], axis="index", copy=copy),
        lambda df, copy: df.rename_axis(index="test", copy=copy),
        lambda df, copy: df.rename_axis(columns="test", copy=copy),
        lambda df, copy: df.astype({"b": "int64"}, copy=copy),
        # lambda df, copy: df.swaplevel(0, 0, copy=copy),
        lambda df, copy: df.truncate(0, 5, copy=copy),
        lambda df, copy: df.infer_objects(copy=copy),
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
        "align",
        "set_axis",
        "rename_axis0",
        "rename_axis1",
        "astype",
        # "swaplevel",  # only series
        "truncate",
        "infer_objects",
        "to_timestamp",
        "to_period",
        "tz_localize",
        "tz_convert",
        "set_flags",
    ],
)
def test_methods_copy_keyword(request, method, copy):
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
    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))


@pytest.mark.filterwarnings("ignore::DeprecationWarning")
@pytest.mark.parametrize("copy", [True, None, False])
@pytest.mark.parametrize(
    "method",
    [
        lambda ser, copy: ser.rename(index={0: 100}, copy=copy),
        lambda ser, copy: ser.rename(None, copy=copy),
        lambda ser, copy: ser.reindex(index=ser.index, copy=copy),
        lambda ser, copy: ser.reindex_like(ser, copy=copy),
        lambda ser, copy: ser.align(ser, copy=copy)[0],
        lambda ser, copy: ser.set_axis(["a", "b", "c"], axis="index", copy=copy),
        lambda ser, copy: ser.rename_axis(index="test", copy=copy),
        lambda ser, copy: ser.astype("int64", copy=copy),
        lambda ser, copy: ser.swaplevel(0, 1, copy=copy),
        lambda ser, copy: ser.truncate(0, 5, copy=copy),
        lambda ser, copy: ser.infer_objects(copy=copy),
        lambda ser, copy: ser.to_timestamp(copy=copy),
        lambda ser, copy: ser.to_period(freq="D", copy=copy),
        lambda ser, copy: ser.tz_localize("US/Central", copy=copy),
        lambda ser, copy: ser.tz_convert("US/Central", copy=copy),
        lambda ser, copy: ser.set_flags(allows_duplicate_labels=False, copy=copy),
    ],
    ids=[
        "rename (dict)",
        "rename",
        "reindex",
        "reindex_like",
        "align",
        "set_axis",
        "rename_axis0",
        "astype",
        "swaplevel",
        "truncate",
        "infer_objects",
        "to_timestamp",
        "to_period",
        "tz_localize",
        "tz_convert",
        "set_flags",
    ],
)
def test_methods_series_copy_keyword(request, method, copy):
    index = None
    if "to_timestamp" in request.node.callspec.id:
        index = period_range("2012-01-01", freq="D", periods=3)
    elif "to_period" in request.node.callspec.id:
        index = date_range("2012-01-01", freq="D", periods=3)
    elif "tz_localize" in request.node.callspec.id:
        index = date_range("2012-01-01", freq="D", periods=3)
    elif "tz_convert" in request.node.callspec.id:
        index = date_range("2012-01-01", freq="D", periods=3, tz="Europe/Brussels")
    elif "swaplevel" in request.node.callspec.id:
        index = MultiIndex.from_arrays([[1, 2, 3], [4, 5, 6]])

    ser = Series([1, 2, 3], index=index)
    ser2 = method(ser, copy=copy)
    assert np.shares_memory(get_array(ser2), get_array(ser))


# -----------------------------------------------------------------------------
# DataFrame methods returning new DataFrame using shallow copy


def test_reset_index():
    # Case: resetting the index (i.e. adding a new column) + mutating the
    # resulting dataframe
    df = DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]}, index=[10, 11, 12]
    )
    df_orig = df.copy()
    df2 = df.reset_index()
    df2._mgr._verify_integrity()

    # still shares memory (df2 is a shallow copy)
    assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
    assert np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    # mutating df2 triggers a copy-on-write for that column / block
    df2.iloc[0, 2] = 0
    assert not np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
    assert np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("index", [pd.RangeIndex(0, 2), Index([1, 2])])
def test_reset_index_series_drop(index):
    ser = Series([1, 2], index=index)
    ser_orig = ser.copy()
    ser2 = ser.reset_index(drop=True)
    assert np.shares_memory(get_array(ser), get_array(ser2))
    assert not ser._mgr._has_no_reference(0)

    ser2.iloc[0] = 100
    tm.assert_series_equal(ser, ser_orig)


def test_groupby_column_index_in_references():
    df = DataFrame(
        {"A": ["a", "b", "c", "d"], "B": [1, 2, 3, 4], "C": ["a", "a", "b", "b"]}
    )
    df = df.set_index("A")
    key = df["C"]
    result = df.groupby(key, observed=True).sum()
    expected = df.groupby("C", observed=True).sum()
    tm.assert_frame_equal(result, expected)


def test_groupby_modify_series():
    # https://github.com/pandas-dev/pandas/issues/63219
    # Modifying a Series after using it to groupby should not impact
    # the groupby operation.
    ser = Series([1, 2, 1])
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    gb = df.groupby(ser)
    ser.iloc[0] = 100
    result = gb.sum()
    expected = DataFrame({"a": [4, 2], "b": [10, 5]}, index=[1, 2])
    tm.assert_frame_equal(result, expected)


def test_rename_columns():
    # Case: renaming columns returns a new dataframe
    # + afterwards modifying the result
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.rename(columns=str.upper)

    assert np.shares_memory(get_array(df2, "A"), get_array(df, "a"))
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "A"), get_array(df, "a"))
    assert np.shares_memory(get_array(df2, "C"), get_array(df, "c"))
    expected = DataFrame({"A": [0, 2, 3], "B": [4, 5, 6], "C": [0.1, 0.2, 0.3]})
    tm.assert_frame_equal(df2, expected)
    tm.assert_frame_equal(df, df_orig)


def test_rename_columns_modify_parent():
    # Case: renaming columns returns a new dataframe
    # + afterwards modifying the original (parent) dataframe
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df2 = df.rename(columns=str.upper)
    df2_orig = df2.copy()

    assert np.shares_memory(get_array(df2, "A"), get_array(df, "a"))
    df.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "A"), get_array(df, "a"))
    assert np.shares_memory(get_array(df2, "C"), get_array(df, "c"))
    expected = DataFrame({"a": [0, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    tm.assert_frame_equal(df, expected)
    tm.assert_frame_equal(df2, df2_orig)


def test_pipe():
    df = DataFrame({"a": [1, 2, 3], "b": 1.5})
    df_orig = df.copy()

    def testfunc(df):
        return df

    df2 = df.pipe(testfunc)

    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    # mutating df2 triggers a copy-on-write for that column
    df2.iloc[0, 0] = 0
    tm.assert_frame_equal(df, df_orig)
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))


def test_pipe_modify_df():
    df = DataFrame({"a": [1, 2, 3], "b": 1.5})
    df_orig = df.copy()

    def testfunc(df):
        df.iloc[0, 0] = 100
        return df

    df2 = df.pipe(testfunc)

    assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))

    tm.assert_frame_equal(df, df_orig)
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))


def test_reindex_columns():
    # Case: reindexing the column returns a new dataframe
    # + afterwards modifying the result
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.reindex(columns=["a", "c"])

    # still shares memory (df2 is a shallow copy)
    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    # mutating df2 triggers a copy-on-write for that column
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    assert np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "index",
    [
        lambda idx: idx,
        lambda idx: idx.view(),
        lambda idx: idx.copy(),
        lambda idx: list(idx),
    ],
    ids=["identical", "view", "copy", "values"],
)
def test_reindex_rows(index):
    # Case: reindexing the rows with an index that matches the current index
    # can use a shallow copy
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.reindex(index=index(df.index))

    # still shares memory (df2 is a shallow copy)
    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    # mutating df2 triggers a copy-on-write for that column
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    assert np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    tm.assert_frame_equal(df, df_orig)


def test_drop_on_column():
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.drop(columns="a")
    df2._mgr._verify_integrity()

    assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
    assert np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
    assert np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    tm.assert_frame_equal(df, df_orig)


def test_select_dtypes():
    # Case: selecting columns using `select_dtypes()` returns a new dataframe
    # + afterwards modifying the result
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.select_dtypes("int64")
    df2._mgr._verify_integrity()

    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    # mutating df2 triggers a copy-on-write for that column/block
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "filter_kwargs", [{"items": ["a"]}, {"like": "a"}, {"regex": "a"}]
)
def test_filter(filter_kwargs):
    # Case: selecting columns using `filter()` returns a new dataframe
    # + afterwards modifying the result
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.filter(**filter_kwargs)
    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    # mutating df2 triggers a copy-on-write for that column/block
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_shift_no_op():
    df = DataFrame(
        [[1, 2], [3, 4], [5, 6]],
        index=date_range("2020-01-01", "2020-01-03"),
        columns=["a", "b"],
    )
    df_orig = df.copy()
    df2 = df.shift(periods=0)
    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    assert df2.index is not df.index
    assert df2.columns is not df.columns

    df.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))
    assert np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
    tm.assert_frame_equal(df2, df_orig)


def test_shift_index():
    df = DataFrame(
        [[1, 2], [3, 4], [5, 6]],
        index=date_range("2020-01-01", "2020-01-03"),
        columns=["a", "b"],
    )
    df2 = df.shift(periods=1, axis=0)

    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    assert df2.index is not df.index
    assert df2.columns is not df.columns


def test_shift_rows_freq():
    df = DataFrame(
        [[1, 2], [3, 4], [5, 6]],
        index=date_range("2020-01-01", "2020-01-03"),
        columns=["a", "b"],
    )
    df_orig = df.copy()
    df_orig.index = date_range("2020-01-02", "2020-01-04")
    df2 = df.shift(periods=1, freq="1D")

    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    df.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df, "a"), get_array(df2, "a"))
    tm.assert_frame_equal(df2, df_orig)


def test_shift_columns():
    df = DataFrame(
        [[1, 2], [3, 4], [5, 6]], columns=date_range("2020-01-01", "2020-01-02")
    )
    df2 = df.shift(periods=1, axis=1)

    assert np.shares_memory(get_array(df2, "2020-01-02"), get_array(df, "2020-01-01"))
    df.iloc[0, 0] = 0
    assert not np.shares_memory(
        get_array(df2, "2020-01-02"), get_array(df, "2020-01-01")
    )
    expected = DataFrame(
        [[np.nan, 1], [np.nan, 3], [np.nan, 5]],
        columns=date_range("2020-01-01", "2020-01-02"),
    )
    tm.assert_frame_equal(df2, expected)


def test_pop():
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    view_original = df[:]
    result = df.pop("a")

    assert np.shares_memory(result.values, get_array(view_original, "a"))
    assert np.shares_memory(get_array(df, "b"), get_array(view_original, "b"))

    result.iloc[0] = 0
    assert not np.shares_memory(result.values, get_array(view_original, "a"))
    df.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df, "b"), get_array(view_original, "b"))
    tm.assert_frame_equal(view_original, df_orig)


@pytest.mark.parametrize(
    "func",
    [
        lambda x, y: x.align(y),
        lambda x, y: x.align(y.a, axis=0),
        lambda x, y: x.align(y.a.iloc[slice(0, 1)], axis=1),
    ],
)
def test_align_frame(func):
    df = DataFrame({"a": [1, 2, 3], "b": "a"})
    df_orig = df.copy()
    df_changed = df[["b", "a"]].copy()
    df2, _ = func(df, df_changed)

    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_align_series():
    ser = Series([1, 2])
    ser_orig = ser.copy()
    ser_other = ser.copy()
    ser2, ser_other_result = ser.align(ser_other)

    assert np.shares_memory(ser2.values, ser.values)
    assert np.shares_memory(ser_other_result.values, ser_other.values)
    ser2.iloc[0] = 0
    ser_other_result.iloc[0] = 0
    assert not np.shares_memory(ser2.values, ser.values)
    assert not np.shares_memory(ser_other_result.values, ser_other.values)
    tm.assert_series_equal(ser, ser_orig)
    tm.assert_series_equal(ser_other, ser_orig)


def test_align_copy_false():
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    df_orig = df.copy()
    df2, df3 = df.align(df)

    assert np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
    assert np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    df2.loc[0, "a"] = 0
    tm.assert_frame_equal(df, df_orig)  # Original is unchanged

    df3.loc[0, "a"] = 0
    tm.assert_frame_equal(df, df_orig)  # Original is unchanged


def test_align_with_series_copy_false():
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    ser = Series([1, 2, 3], name="x")
    ser_orig = ser.copy()
    df_orig = df.copy()
    df2, ser2 = df.align(ser, axis=0)

    assert np.shares_memory(get_array(df, "b"), get_array(df2, "b"))
    assert np.shares_memory(get_array(df, "a"), get_array(df2, "a"))
    assert np.shares_memory(get_array(ser, "x"), get_array(ser2, "x"))

    df2.loc[0, "a"] = 0
    tm.assert_frame_equal(df, df_orig)  # Original is unchanged

    ser2.loc[0] = 0
    tm.assert_series_equal(ser, ser_orig)  # Original is unchanged


def test_to_frame():
    # Case: converting a Series to a DataFrame with to_frame
    ser = Series([1, 2, 3])
    ser_orig = ser.copy()

    df = ser[:].to_frame()

    # currently this always returns a "view"
    assert np.shares_memory(ser.values, get_array(df, 0))

    df.iloc[0, 0] = 0

    # mutating df triggers a copy-on-write for that column
    assert not np.shares_memory(ser.values, get_array(df, 0))
    tm.assert_series_equal(ser, ser_orig)

    # modify original series -> don't modify dataframe
    df = ser[:].to_frame()
    ser.iloc[0] = 0

    tm.assert_frame_equal(df, ser_orig.to_frame())

    df = ser.to_frame()
    assert df.index is not ser.index


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
def test_chained_methods(method, idx):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()

    # modify df2 -> don't modify df
    df2 = method(df)
    df2.iloc[0, idx] = 0
    tm.assert_frame_equal(df, df_orig)

    # modify df -> don't modify df2
    df2 = method(df)
    df.iloc[0, 0] = 0
    tm.assert_frame_equal(df2.iloc[:, idx:], df_orig)


@pytest.mark.parametrize("obj", [Series([1, 2], name="a"), DataFrame({"a": [1, 2]})])
def test_to_timestamp(obj):
    obj.index = Index([Period("2012-1-1", freq="D"), Period("2012-1-2", freq="D")])

    obj_orig = obj.copy()
    obj2 = obj.to_timestamp()

    assert np.shares_memory(get_array(obj2, "a"), get_array(obj, "a"))

    # mutating obj2 triggers a copy-on-write for that column / block
    obj2.iloc[0] = 0
    assert not np.shares_memory(get_array(obj2, "a"), get_array(obj, "a"))
    tm.assert_equal(obj, obj_orig)


@pytest.mark.parametrize("obj", [Series([1, 2], name="a"), DataFrame({"a": [1, 2]})])
def test_to_period(obj):
    obj.index = Index([Timestamp("2019-12-31"), Timestamp("2020-12-31")])

    obj_orig = obj.copy()
    obj2 = obj.to_period(freq="Y")

    assert np.shares_memory(get_array(obj2, "a"), get_array(obj, "a"))

    # mutating obj2 triggers a copy-on-write for that column / block
    obj2.iloc[0] = 0
    assert not np.shares_memory(get_array(obj2, "a"), get_array(obj, "a"))
    tm.assert_equal(obj, obj_orig)


def test_set_index():
    # GH 49473
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.set_index("a")

    assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))

    # mutating df2 triggers a copy-on-write for that column / block
    df2.iloc[0, 1] = 0
    assert not np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    tm.assert_frame_equal(df, df_orig)


def test_set_index_mutating_parent_does_not_mutate_index():
    df = DataFrame({"a": [1, 2, 3], "b": 1})
    result = df.set_index("a")
    expected = result.copy()

    df.iloc[0, 0] = 100
    tm.assert_frame_equal(result, expected)


def test_add_prefix():
    # GH 49473
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.add_prefix("CoW_")

    assert np.shares_memory(get_array(df2, "CoW_a"), get_array(df, "a"))
    df2.iloc[0, 0] = 0

    assert not np.shares_memory(get_array(df2, "CoW_a"), get_array(df, "a"))

    assert np.shares_memory(get_array(df2, "CoW_c"), get_array(df, "c"))
    expected = DataFrame(
        {"CoW_a": [0, 2, 3], "CoW_b": [4, 5, 6], "CoW_c": [0.1, 0.2, 0.3]}
    )
    tm.assert_frame_equal(df2, expected)
    tm.assert_frame_equal(df, df_orig)


def test_add_suffix():
    # GH 49473
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.add_suffix("_CoW")
    assert np.shares_memory(get_array(df2, "a_CoW"), get_array(df, "a"))
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "a_CoW"), get_array(df, "a"))
    assert np.shares_memory(get_array(df2, "c_CoW"), get_array(df, "c"))
    expected = DataFrame(
        {"a_CoW": [0, 2, 3], "b_CoW": [4, 5, 6], "c_CoW": [0.1, 0.2, 0.3]}
    )
    tm.assert_frame_equal(df2, expected)
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("axis, val", [(0, 5.5), (1, np.nan)])
def test_dropna(axis, val):
    df = DataFrame({"a": [1, 2, 3], "b": [4, val, 6], "c": "d"})
    df_orig = df.copy()
    df2 = df.dropna(axis=axis)

    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("val", [5, 5.5])
def test_dropna_series(val):
    ser = Series([1, val, 4])
    ser_orig = ser.copy()
    ser2 = ser.dropna()
    assert np.shares_memory(ser2.values, ser.values)

    ser2.iloc[0] = 0
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
def test_head_tail(method):
    df = DataFrame({"a": [1, 2, 3], "b": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = method(df)
    df2._mgr._verify_integrity()

    # We are explicitly deviating for CoW here to make an eager copy (avoids
    # tracking references for very cheap ops)
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    assert not np.shares_memory(get_array(df2, "b"), get_array(df, "b"))

    # modify df2 to trigger CoW for that block
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_infer_objects(using_infer_string):
    df = DataFrame(
        {"a": [1, 2], "b": Series(["x", "y"], dtype=object), "c": 1, "d": "x"}
    )
    df_orig = df.copy()
    df2 = df.infer_objects()

    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    if using_infer_string and HAS_PYARROW:
        assert not tm.shares_memory(get_array(df2, "b"), get_array(df, "b"))
    else:
        assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))

    df2.iloc[0, 0] = 0
    df2.iloc[0, 1] = "d"
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    assert not np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
    tm.assert_frame_equal(df, df_orig)


def test_infer_objects_no_reference(using_infer_string):
    df = DataFrame(
        {
            "a": [1, 2],
            "b": Series(["x", "y"], dtype=object),
            "c": 1,
            "d": Series(
                [Timestamp("2019-12-31"), Timestamp("2020-12-31")], dtype="object"
            ),
            "e": Series(["z", "w"], dtype=object),
        }
    )
    df = df.infer_objects()

    arr_a = get_array(df, "a")
    arr_b = get_array(df, "b")
    arr_d = get_array(df, "d")

    df.iloc[0, 0] = 0
    df.iloc[0, 1] = "d"
    df.iloc[0, 3] = Timestamp("2018-12-31")
    assert np.shares_memory(arr_a, get_array(df, "a"))
    if using_infer_string and HAS_PYARROW:
        # note that the underlying memory of arr_b has been copied anyway
        # because of the assignment, but the EA is updated inplace so still
        # appears the share memory
        assert tm.shares_memory(arr_b, get_array(df, "b"))
    else:
        # TODO(CoW): Block splitting causes references here
        assert not np.shares_memory(arr_b, get_array(df, "b"))
    assert np.shares_memory(arr_d, get_array(df, "d"))


def test_infer_objects_reference():
    df = DataFrame(
        {
            "a": [1, 2],
            "b": Series(["x", "y"], dtype=object),
            "c": 1,
            "d": Series(
                [Timestamp("2019-12-31"), Timestamp("2020-12-31")], dtype="object"
            ),
        }
    )
    view = df[:]  # noqa: F841
    df = df.infer_objects()

    arr_a = get_array(df, "a")
    arr_b = get_array(df, "b")
    arr_d = get_array(df, "d")

    df.iloc[0, 0] = 0
    df.iloc[0, 1] = "d"
    df.iloc[0, 3] = Timestamp("2018-12-31")
    assert not np.shares_memory(arr_a, get_array(df, "a"))
    assert not np.shares_memory(arr_b, get_array(df, "b"))
    assert np.shares_memory(arr_d, get_array(df, "d"))


@pytest.mark.parametrize(
    "kwargs",
    [
        {"before": "a", "after": "b", "axis": 1},
        {"before": 0, "after": 1, "axis": 0},
    ],
)
def test_truncate(kwargs):
    df = DataFrame({"a": [1, 2, 3], "b": 1, "c": 2})
    df_orig = df.copy()
    df2 = df.truncate(**kwargs)
    df2._mgr._verify_integrity()

    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("method", ["assign", "drop_duplicates"])
def test_assign_drop_duplicates(method):
    df = DataFrame({"a": [1, 2, 3]})
    df_orig = df.copy()
    df2 = getattr(df, method)()
    df2._mgr._verify_integrity()

    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("obj", [Series([1, 2]), DataFrame({"a": [1, 2]})])
def test_take(obj):
    # Check that no copy is made when we take all rows in original order
    obj_orig = obj.copy()
    obj2 = obj.take([0, 1])
    assert np.shares_memory(obj2.values, obj.values)

    obj2.iloc[0] = 0
    assert not np.shares_memory(obj2.values, obj.values)
    tm.assert_equal(obj, obj_orig)


@pytest.mark.parametrize("obj", [Series([1, 2]), DataFrame({"a": [1, 2]})])
def test_between_time(obj):
    obj.index = date_range("2018-04-09", periods=2, freq="1D20min")
    obj_orig = obj.copy()
    obj2 = obj.between_time("0:00", "1:00")
    assert np.shares_memory(obj2.values, obj.values)

    obj2.iloc[0] = 0
    assert not np.shares_memory(obj2.values, obj.values)
    tm.assert_equal(obj, obj_orig)


def test_reindex_like():
    df = DataFrame({"a": [1, 2], "b": "a"})
    other = DataFrame({"b": "a", "a": [1, 2]})

    df_orig = df.copy()
    df2 = df.reindex_like(other)
    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    df2.iloc[0, 1] = 0
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_sort_index():
    # GH 49473
    ser = Series([1, 2, 3])
    ser_orig = ser.copy()
    ser2 = ser.sort_index()
    assert np.shares_memory(ser.values, ser2.values)

    # mutating ser triggers a copy-on-write for the column / block
    ser2.iloc[0] = 0
    assert not np.shares_memory(ser2.values, ser.values)
    tm.assert_series_equal(ser, ser_orig)


@pytest.mark.parametrize(
    "obj, kwargs",
    [(Series([1, 2, 3], name="a"), {}), (DataFrame({"a": [1, 2, 3]}), {"by": "a"})],
)
def test_sort_values(obj, kwargs):
    obj_orig = obj.copy()
    obj2 = obj.sort_values(**kwargs)
    assert np.shares_memory(get_array(obj2, "a"), get_array(obj, "a"))

    # mutating df triggers a copy-on-write for the column / block
    obj2.iloc[0] = 0
    assert not np.shares_memory(get_array(obj2, "a"), get_array(obj, "a"))
    tm.assert_equal(obj, obj_orig)


@pytest.mark.parametrize(
    "obj, kwargs",
    [(Series([1, 2, 3], name="a"), {}), (DataFrame({"a": [1, 2, 3]}), {"by": "a"})],
)
def test_sort_values_inplace(obj, kwargs):
    obj_orig = obj.copy()
    view = obj[:]
    obj.sort_values(inplace=True, **kwargs)

    assert np.shares_memory(get_array(obj, "a"), get_array(view, "a"))

    # mutating obj triggers a copy-on-write for the column / block
    obj.iloc[0] = 0
    assert not np.shares_memory(get_array(obj, "a"), get_array(view, "a"))
    tm.assert_equal(view, obj_orig)


@pytest.mark.parametrize("decimals", [-1, 0, 1])
def test_round(decimals):
    df = DataFrame({"a": [1, 2], "b": "c"})
    df_orig = df.copy()
    df2 = df.round(decimals=decimals)

    assert tm.shares_memory(get_array(df2, "b"), get_array(df, "b"))
    # TODO: Make inplace by using out parameter of ndarray.round?
    if decimals >= 0 and Version(np.__version__) < Version("2.4.0.dev0"):
        # Ensure lazy copy if no-op
        # TODO: Cannot rely on Numpy returning view after version 2.3
        assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    else:
        assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    assert df2.index is not df.index
    assert df2.columns is not df.columns

    df2.iloc[0, 1] = "d"
    df2.iloc[0, 0] = 4
    assert not np.shares_memory(get_array(df2, "b"), get_array(df, "b"))
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_reorder_levels():
    index = MultiIndex.from_tuples(
        [(1, 1), (1, 2), (2, 1), (2, 2)], names=["one", "two"]
    )
    df = DataFrame({"a": [1, 2, 3, 4]}, index=index)
    df_orig = df.copy()
    df2 = df.reorder_levels(order=["two", "one"])
    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_series_reorder_levels():
    index = MultiIndex.from_tuples(
        [(1, 1), (1, 2), (2, 1), (2, 2)], names=["one", "two"]
    )
    ser = Series([1, 2, 3, 4], index=index)
    ser_orig = ser.copy()
    ser2 = ser.reorder_levels(order=["two", "one"])
    assert np.shares_memory(ser2.values, ser.values)

    ser2.iloc[0] = 0
    assert not np.shares_memory(ser2.values, ser.values)
    tm.assert_series_equal(ser, ser_orig)


@pytest.mark.parametrize("obj", [Series([1, 2, 3]), DataFrame({"a": [1, 2, 3]})])
def test_swaplevel(obj):
    index = MultiIndex.from_tuples([(1, 1), (1, 2), (2, 1)], names=["one", "two"])
    obj.index = index
    obj_orig = obj.copy()
    obj2 = obj.swaplevel()
    assert np.shares_memory(obj2.values, obj.values)

    obj2.iloc[0] = 0
    assert not np.shares_memory(obj2.values, obj.values)
    tm.assert_equal(obj, obj_orig)


def test_frame_set_axis():
    # GH 49473
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [0.1, 0.2, 0.3]})
    df_orig = df.copy()
    df2 = df.set_axis(["a", "b", "c"], axis="index")

    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    # mutating df2 triggers a copy-on-write for that column / block
    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_series_set_axis():
    # GH 49473
    ser = Series([1, 2, 3])
    ser_orig = ser.copy()
    ser2 = ser.set_axis(["a", "b", "c"], axis="index")
    assert np.shares_memory(ser, ser2)

    # mutating ser triggers a copy-on-write for the column / block
    ser2.iloc[0] = 0
    assert not np.shares_memory(ser2, ser)
    tm.assert_series_equal(ser, ser_orig)


def test_set_flags():
    ser = Series([1, 2, 3])
    ser_orig = ser.copy()
    ser2 = ser.set_flags(allows_duplicate_labels=False)

    assert np.shares_memory(ser, ser2)

    # mutating ser triggers a copy-on-write for the column / block
    ser2.iloc[0] = 0
    assert not np.shares_memory(ser2, ser)
    tm.assert_series_equal(ser, ser_orig)


@pytest.mark.parametrize("kwargs", [{"mapper": "test"}, {"index": "test"}])
def test_rename_axis(kwargs):
    df = DataFrame({"a": [1, 2, 3, 4]}, index=Index([1, 2, 3, 4], name="a"))
    df_orig = df.copy()
    df2 = df.rename_axis(**kwargs)
    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    df2.iloc[0, 0] = 0
    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize(
    "func, tz", [("tz_convert", "Europe/Berlin"), ("tz_localize", None)]
)
def test_tz_convert_localize(func, tz):
    # GH 49473
    ser = Series(
        [1, 2], index=date_range(start="2014-08-01 09:00", freq="h", periods=2, tz=tz)
    )
    ser_orig = ser.copy()
    ser2 = getattr(ser, func)("US/Central")
    assert np.shares_memory(ser.values, ser2.values)

    # mutating ser triggers a copy-on-write for the column / block
    ser2.iloc[0] = 0
    assert not np.shares_memory(ser2.values, ser.values)
    tm.assert_series_equal(ser, ser_orig)


def test_droplevel():
    # GH 49473
    index = MultiIndex.from_tuples([(1, 1), (1, 2), (2, 1)], names=["one", "two"])
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]}, index=index)
    df_orig = df.copy()
    df2 = df.droplevel(0)

    assert np.shares_memory(get_array(df2, "c"), get_array(df, "c"))
    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    # mutating df2 triggers a copy-on-write for that column / block
    df2.iloc[0, 0] = 0

    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    assert np.shares_memory(get_array(df2, "b"), get_array(df, "b"))

    tm.assert_frame_equal(df, df_orig)


def test_squeeze():
    df = DataFrame({"a": [1, 2, 3]})
    df_orig = df.copy()
    series = df.squeeze()

    # Should share memory regardless of CoW since squeeze is just an iloc
    assert np.shares_memory(series.values, get_array(df, "a"))

    # mutating squeezed df triggers a copy-on-write for that column/block
    series.iloc[0] = 0
    assert not np.shares_memory(series.values, get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_items():
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})
    df_orig = df.copy()

    # Test this twice, since the second time, the item cache will be
    # triggered, and we want to make sure it still works then.
    for i in range(2):
        for name, ser in df.items():
            assert np.shares_memory(get_array(ser, name), get_array(df, name))

            # mutating df triggers a copy-on-write for that column / block
            ser.iloc[0] = 0

            assert not np.shares_memory(get_array(ser, name), get_array(df, name))
            tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("dtype", ["int64", "Int64"])
def test_putmask(dtype):
    df = DataFrame({"a": [1, 2], "b": 1, "c": 2}, dtype=dtype)
    view = df[:]
    df_orig = df.copy()
    df[df == df] = 5

    assert not np.shares_memory(get_array(view, "a"), get_array(df, "a"))
    tm.assert_frame_equal(view, df_orig)


@pytest.mark.parametrize("dtype", ["int64", "Int64"])
def test_putmask_no_reference(dtype):
    df = DataFrame({"a": [1, 2], "b": 1, "c": 2}, dtype=dtype)
    arr_a = get_array(df, "a")
    df[df == df] = 5
    assert np.shares_memory(arr_a, get_array(df, "a"))


@pytest.mark.parametrize("dtype", ["float64", "Float64"])
def test_putmask_aligns_rhs_no_reference(dtype):
    df = DataFrame({"a": [1.5, 2], "b": 1.5}, dtype=dtype)
    arr_a = get_array(df, "a")
    df[df == df] = DataFrame({"a": [5.5, 5]})
    assert np.shares_memory(arr_a, get_array(df, "a"))


@pytest.mark.parametrize("val, exp, raises", [(5.5, True, True), (5, False, False)])
def test_putmask_dont_copy_some_blocks(val, exp, raises: bool):
    df = DataFrame({"a": [1, 2], "b": 1, "c": 1.5})
    view = df[:]
    df_orig = df.copy()
    indexer = DataFrame(
        [[True, False, False], [True, False, False]], columns=list("abc")
    )
    if raises:
        with pytest.raises(TypeError, match="Invalid value"):
            df[indexer] = val
    else:
        df[indexer] = val
        assert not np.shares_memory(get_array(view, "a"), get_array(df, "a"))
        # TODO(CoW): Could split blocks to avoid copying the whole block
        assert np.shares_memory(get_array(view, "b"), get_array(df, "b")) is exp
        assert np.shares_memory(get_array(view, "c"), get_array(df, "c"))
        assert df._mgr._has_no_reference(1) is not exp
        assert not df._mgr._has_no_reference(2)
        tm.assert_frame_equal(view, df_orig)


@pytest.mark.parametrize("dtype", ["int64", "Int64"])
@pytest.mark.parametrize(
    "func",
    [
        lambda ser: ser.where(ser > 0, 10),
        lambda ser: ser.mask(ser <= 0, 10),
    ],
)
def test_where_mask_noop(dtype, func):
    ser = Series([1, 2, 3], dtype=dtype)
    ser_orig = ser.copy()

    result = func(ser)
    assert np.shares_memory(get_array(ser), get_array(result))
    assert result.index is not ser.index

    result.iloc[0] = 10
    assert not np.shares_memory(get_array(ser), get_array(result))
    tm.assert_series_equal(ser, ser_orig)


@pytest.mark.parametrize("dtype", ["int64", "Int64"])
@pytest.mark.parametrize(
    "func",
    [
        lambda ser: ser.where(ser < 0, 10),
        lambda ser: ser.mask(ser >= 0, 10),
    ],
)
def test_where_mask(dtype, func):
    ser = Series([1, 2, 3], dtype=dtype)
    ser_orig = ser.copy()

    result = func(ser)

    assert not np.shares_memory(get_array(ser), get_array(result))
    assert result.index is not ser.index
    tm.assert_series_equal(ser, ser_orig)


@pytest.mark.parametrize("dtype, val", [("int64", 10.5), ("Int64", 10)])
@pytest.mark.parametrize(
    "func",
    [
        lambda df, val: df.where(df < 0, val),
        lambda df, val: df.mask(df >= 0, val),
    ],
)
def test_where_mask_noop_on_single_column(dtype, val, func):
    df = DataFrame({"a": [1, 2, 3], "b": [-4, -5, -6]}, dtype=dtype)
    df_orig = df.copy()

    result = func(df, val)
    assert np.shares_memory(get_array(df, "b"), get_array(result, "b"))
    assert not np.shares_memory(get_array(df, "a"), get_array(result, "a"))

    result.iloc[0, 1] = 10
    assert not np.shares_memory(get_array(df, "b"), get_array(result, "b"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("func", ["mask", "where"])
def test_chained_where_mask(func):
    df = DataFrame({"a": [1, 4, 2], "b": 1})
    df_orig = df.copy()
    with tm.raises_chained_assignment_error():
        getattr(df["a"], func)(df["a"] > 2, 5, inplace=True)
    tm.assert_frame_equal(df, df_orig)

    with tm.raises_chained_assignment_error():
        getattr(df[["a"]], func)(df["a"] > 2, 5, inplace=True)
    tm.assert_frame_equal(df, df_orig)


def test_asfreq_noop():
    df = DataFrame(
        {"a": [0.0, None, 2.0, 3.0]},
        index=date_range("1/1/2000", periods=4, freq="min"),
    )
    df_orig = df.copy()
    df2 = df.asfreq(freq="min")
    assert np.shares_memory(get_array(df2, "a"), get_array(df, "a"))

    # mutating df2 triggers a copy-on-write for that column / block
    df2.iloc[0, 0] = 0

    assert not np.shares_memory(get_array(df2, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


def test_iterrows():
    df = DataFrame({"a": 0, "b": 1}, index=[1, 2, 3])
    df_orig = df.copy()

    for _, sub in df.iterrows():
        sub.iloc[0] = 100
    tm.assert_frame_equal(df, df_orig)


def test_interpolate_creates_copy():
    # GH#51126
    df = DataFrame({"a": [1.5, np.nan, 3]})
    view = df[:]
    expected = df.copy()

    df.ffill(inplace=True)
    df.iloc[0, 0] = 100.5
    tm.assert_frame_equal(view, expected)


def test_isetitem():
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})
    df_orig = df.copy()
    df2 = df.copy(deep=False)  # Trigger a CoW
    df2.isetitem(1, np.array([-1, -2, -3]))  # This is inplace
    assert np.shares_memory(get_array(df, "c"), get_array(df2, "c"))
    assert np.shares_memory(get_array(df, "a"), get_array(df2, "a"))

    df2.loc[0, "a"] = 0
    tm.assert_frame_equal(df, df_orig)  # Original is unchanged
    assert np.shares_memory(get_array(df, "c"), get_array(df2, "c"))


@pytest.mark.parametrize(
    "dtype", ["int64", "float64"], ids=["single-block", "mixed-block"]
)
def test_isetitem_series(dtype):
    df = DataFrame({"a": [1, 2, 3], "b": np.array([4, 5, 6], dtype=dtype)})
    ser = Series([7, 8, 9])
    ser_orig = ser.copy()
    df.isetitem(0, ser)

    assert np.shares_memory(get_array(df, "a"), get_array(ser))
    assert not df._mgr._has_no_reference(0)

    # mutating dataframe doesn't update series
    df.loc[0, "a"] = 0
    tm.assert_series_equal(ser, ser_orig)

    # mutating series doesn't update dataframe
    df = DataFrame({"a": [1, 2, 3], "b": np.array([4, 5, 6], dtype=dtype)})
    ser = Series([7, 8, 9])
    df.isetitem(0, ser)

    ser.loc[0] = 0
    expected = DataFrame({"a": [7, 8, 9], "b": np.array([4, 5, 6], dtype=dtype)})
    tm.assert_frame_equal(df, expected)


def test_isetitem_frame():
    df = DataFrame({"a": [1, 2, 3], "b": 1, "c": 2})
    rhs = DataFrame({"a": [4, 5, 6], "b": 2})
    df.isetitem([0, 1], rhs)
    assert np.shares_memory(get_array(df, "a"), get_array(rhs, "a"))
    assert np.shares_memory(get_array(df, "b"), get_array(rhs, "b"))
    assert not df._mgr._has_no_reference(0)
    expected = df.copy()
    rhs.iloc[0, 0] = 100
    rhs.iloc[0, 1] = 100
    tm.assert_frame_equal(df, expected)


@pytest.mark.parametrize("key", ["a", ["a"]])
def test_get(key):
    df = DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    df_orig = df.copy()

    result = df.get(key)

    assert np.shares_memory(get_array(result, "a"), get_array(df, "a"))
    result.iloc[0] = 0
    assert not np.shares_memory(get_array(result, "a"), get_array(df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("axis, key", [(0, 0), (1, "a")])
@pytest.mark.parametrize(
    "dtype", ["int64", "float64"], ids=["single-block", "mixed-block"]
)
def test_xs(axis, key, dtype):
    single_block = dtype == "int64"
    df = DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6], "c": np.array([7, 8, 9], dtype=dtype)}
    )
    df_orig = df.copy()

    result = df.xs(key, axis=axis)

    if axis == 1 or single_block:
        assert np.shares_memory(get_array(df, "a"), get_array(result))
    else:
        assert result._mgr._has_no_reference(0)
    if axis == 0:
        assert result.index is not df.columns
    else:
        assert result.index is not df.index

    result.iloc[0] = 0
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("axis", [0, 1])
@pytest.mark.parametrize("key, level", [("l1", 0), (2, 1)])
def test_xs_multiindex(key, level, axis):
    arr = np.arange(18).reshape(6, 3)
    index = MultiIndex.from_product([["l1", "l2"], [1, 2, 3]], names=["lev1", "lev2"])
    df = DataFrame(arr, index=index, columns=list("abc"))
    if axis == 1:
        df = df.transpose().copy()
    df_orig = df.copy()

    result = df.xs(key, level=level, axis=axis)

    if level == 0:
        assert np.shares_memory(
            get_array(df, df.columns[0]), get_array(result, result.columns[0])
        )
    assert result.index is not df.index
    assert result.columns is not df.columns

    result.iloc[0, 0] = 0
    tm.assert_frame_equal(df, df_orig)


def test_update_frame():
    df1 = DataFrame({"a": [1.0, 2.0, 3.0], "b": [4.0, 5.0, 6.0]})
    df2 = DataFrame({"b": [100.0]}, index=[1])
    df1_orig = df1.copy()
    view = df1[:]
    df1.update(df2)

    expected = DataFrame({"a": [1.0, 2.0, 3.0], "b": [4.0, 100.0, 6.0]})
    tm.assert_frame_equal(df1, expected)
    # df1 is updated, but its view not
    tm.assert_frame_equal(view, df1_orig)
    assert np.shares_memory(get_array(df1, "a"), get_array(view, "a"))
    assert not np.shares_memory(get_array(df1, "b"), get_array(view, "b"))


def test_update_series():
    ser1 = Series([1.0, 2.0, 3.0])
    ser2 = Series([100.0], index=[1])
    ser1_orig = ser1.copy()
    view = ser1[:]

    ser1.update(ser2)

    expected = Series([1.0, 100.0, 3.0])
    tm.assert_series_equal(ser1, expected)
    # ser1 is updated, but its view not
    tm.assert_series_equal(view, ser1_orig)


def test_update_chained_assignment():
    df = DataFrame({"a": [1, 2, 3]})
    ser2 = Series([100.0], index=[1])
    df_orig = df.copy()
    with tm.raises_chained_assignment_error():
        df["a"].update(ser2)
    tm.assert_frame_equal(df, df_orig)

    with tm.raises_chained_assignment_error():
        df[["a"]].update(ser2.to_frame())
    tm.assert_frame_equal(df, df_orig)


def test_inplace_arithmetic_series():
    ser = Series([1, 2, 3])
    ser_orig = ser.copy()
    data = get_array(ser)
    ser *= 2
    # https://github.com/pandas-dev/pandas/pull/55745
    # changed to NOT update inplace because there is no benefit (actual
    # operation already done non-inplace). This was only for the optics
    # of updating the backing array inplace, but we no longer want to make
    # that guarantee
    assert not np.shares_memory(get_array(ser), data)
    tm.assert_numpy_array_equal(data, get_array(ser_orig))


def test_inplace_arithmetic_series_with_reference():
    ser = Series([1, 2, 3])
    ser_orig = ser.copy()
    view = ser[:]
    ser *= 2
    assert not np.shares_memory(get_array(ser), get_array(view))
    tm.assert_series_equal(ser_orig, view)


def test_transpose():
    df = DataFrame({"a": [1, 2, 3], "b": 1})
    df_orig = df.copy()
    result = df.transpose()
    assert np.shares_memory(get_array(df, "a"), get_array(result, 0))

    result.iloc[0, 0] = 100
    tm.assert_frame_equal(df, df_orig)


def test_transpose_different_dtypes():
    df = DataFrame({"a": [1, 2, 3], "b": 1.5})
    df_orig = df.copy()
    result = df.T

    assert not np.shares_memory(get_array(df, "a"), get_array(result, 0))
    result.iloc[0, 0] = 100
    tm.assert_frame_equal(df, df_orig)


def test_transpose_ea_single_column():
    df = DataFrame({"a": [1, 2, 3]}, dtype="Int64")
    result = df.T

    assert not np.shares_memory(get_array(df, "a"), get_array(result, 0))


def test_transform_frame():
    df = DataFrame({"a": [1, 2, 3], "b": 1})
    df_orig = df.copy()

    def func(ser):
        ser.iloc[0] = 100
        return ser

    df.transform(func)
    tm.assert_frame_equal(df, df_orig)


def test_transform_series():
    ser = Series([1, 2, 3])
    ser_orig = ser.copy()

    def func(ser):
        ser.iloc[0] = 100
        return ser

    ser.transform(func)
    tm.assert_series_equal(ser, ser_orig)


def test_count_read_only_array():
    df = DataFrame({"a": [1, 2], "b": 3})
    result = df.count()
    result.iloc[0] = 100
    expected = Series([100, 2], index=["a", "b"])
    tm.assert_series_equal(result, expected)


def test_insert_series():
    df = DataFrame({"a": [1, 2, 3]})
    ser = Series([1, 2, 3])
    ser_orig = ser.copy()
    df.insert(loc=1, value=ser, column="b")
    assert np.shares_memory(get_array(ser), get_array(df, "b"))
    assert not df._mgr._has_no_reference(1)

    df.iloc[0, 1] = 100
    tm.assert_series_equal(ser, ser_orig)


def test_eval():
    df = DataFrame({"a": [1, 2, 3], "b": 1})
    df_orig = df.copy()

    result = df.eval("c = a+b")
    assert np.shares_memory(get_array(df, "a"), get_array(result, "a"))

    result.iloc[0, 0] = 100
    tm.assert_frame_equal(df, df_orig)


def test_eval_inplace():
    df = DataFrame({"a": [1, 2, 3], "b": 1})
    df_orig = df.copy()
    df_view = df[:]

    df.eval("c = a+b", inplace=True)
    assert np.shares_memory(get_array(df, "a"), get_array(df_view, "a"))

    df.iloc[0, 0] = 100
    tm.assert_frame_equal(df_view, df_orig)


def test_apply_modify_row():
    # Case: applying a function on each row as a Series object, where the
    # function mutates the row object (which needs to trigger CoW if row is a view)
    df = DataFrame({"A": [1, 2], "B": [3, 4]})
    df_orig = df.copy()

    def transform(row):
        row["B"] = 100
        return row

    df.apply(transform, axis=1)

    tm.assert_frame_equal(df, df_orig)

    # row Series is a copy
    df = DataFrame({"A": [1, 2], "B": ["b", "c"]})
    df_orig = df.copy()

    with tm.assert_produces_warning(None):
        df.apply(transform, axis=1)

    tm.assert_frame_equal(df, df_orig)


def test_reduce():
    df = DataFrame({"a": [1, 2, 3], "b": 1.5})

    result = df.sum()
    assert result.index is not df.columns

    result = df.groupby([0, 0, 1]).sum()
    assert result.columns is not df.columns

    result = df.quantile(0.5)
    assert result.index is not df.columns
    result = df.quantile([0.25, 0.5, 0.75])
    assert result.columns is not df.columns


def test_diff():
    df = DataFrame({"a": [1, 2, 3], "b": 1.5})

    result = df.diff()
    assert result.index is not df.index
    assert result.columns is not df.columns

    ser = Series([1, 2, 3])
    result = ser.diff()
    assert result.index is not ser.index
