import numpy as np
import pytest

from pandas._config import using_string_dtype

import pandas as pd
from pandas import (
    DataFrame,
    DatetimeIndex,
    Index,
    Period,
    PeriodIndex,
    Series,
    Timedelta,
    TimedeltaIndex,
    Timestamp,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array

# -----------------------------------------------------------------------------
# Copy/view behaviour for Series / DataFrame constructors


@pytest.mark.parametrize("dtype", [None, "int64"])
def test_series_from_series(dtype):
    # Case: constructing a Series from another Series object follows CoW rules:
    # a new object is returned and thus mutations are not propagated
    ser = Series([1, 2, 3], name="name")

    # default is copy=False -> new Series is a shallow copy / view of original
    result = Series(ser, dtype=dtype)

    # the shallow copy still shares memory
    assert np.shares_memory(get_array(ser), get_array(result))

    assert result._mgr.blocks[0].refs.has_reference()

    # mutating new series copy doesn't mutate original
    result.iloc[0] = 0
    assert ser.iloc[0] == 1
    # mutating triggered a copy-on-write -> no longer shares memory
    assert not np.shares_memory(get_array(ser), get_array(result))

    # the same when modifying the parent
    result = Series(ser, dtype=dtype)

    # mutating original doesn't mutate new series
    ser.iloc[0] = 0
    assert result.iloc[0] == 1

    # forcing copy=False still gives a CoW shallow copy
    result = Series(ser, dtype=dtype, copy=False)
    assert np.shares_memory(get_array(ser), get_array(result))
    assert result._mgr.blocks[0].refs.has_reference()

    # forcing copy=True still results in an actual hard copy up front
    result = Series(ser, dtype=dtype, copy=True)
    assert not np.shares_memory(get_array(ser), get_array(result))
    assert ser._mgr._has_no_reference(0)


def test_series_from_series_with_reindex():
    # Case: constructing a Series from another Series with specifying an index
    # that potentially requires a reindex of the values
    ser = Series([1, 2, 3], name="name")

    # passing an index that doesn't actually require a reindex of the values
    # -> still getting a CoW shallow copy
    for index in [
        ser.index,
        ser.index.copy(),
        list(ser.index),
        ser.index.rename("idx"),
    ]:
        result = Series(ser, index=index)
        assert np.shares_memory(ser.values, result.values)
        result.iloc[0] = 0
        assert ser.iloc[0] == 1

        # forcing copy=True still results in an actual hard copy up front
        result = Series(ser, index=index, copy=True)
        assert not np.shares_memory(ser.values, result.values)
        assert not result._mgr.blocks[0].refs.has_reference()

    # ensure that if an actual reindex is needed, we don't have any refs
    # (mutating the result wouldn't trigger CoW)
    result = Series(ser, index=[0, 1, 2, 3])
    assert not np.shares_memory(ser.values, result.values)
    assert not result._mgr.blocks[0].refs.has_reference()


@pytest.mark.parametrize("dtype", [None, "int64"])
@pytest.mark.parametrize("idx", [None, pd.RangeIndex(start=0, stop=3, step=1)])
@pytest.mark.parametrize(
    "arr", [np.array([1, 2, 3], dtype="int64"), pd.array([1, 2, 3], dtype="Int64")]
)
def test_series_from_array(idx, dtype, arr):
    ser = Series(arr, dtype=dtype, index=idx)
    ser_orig = ser.copy()
    data = getattr(arr, "_data", arr)
    assert not np.shares_memory(get_array(ser), data)

    arr[0] = 100
    tm.assert_series_equal(ser, ser_orig)

    # if the user explicitly passes copy=False, we get an actual view
    # not protected by CoW
    ser = Series(arr, dtype=dtype, index=idx, copy=False)
    assert np.shares_memory(get_array(ser), data)
    arr[0] = 50
    assert ser.iloc[0] == 50


@pytest.mark.parametrize("copy", [True, False, None])
def test_series_from_array_different_dtype(copy):
    arr = np.array([1, 2, 3], dtype="int64")
    ser = Series(arr, dtype="int32", copy=copy)
    assert not np.shares_memory(get_array(ser), arr)


@pytest.mark.parametrize(
    "idx",
    [
        Index([1, 2]),
        DatetimeIndex([Timestamp("2019-12-31"), Timestamp("2020-12-31")]),
        PeriodIndex([Period("2019-12-31"), Period("2020-12-31")]),
        TimedeltaIndex([Timedelta("1 days"), Timedelta("2 days")]),
    ],
)
def test_series_from_index(idx):
    ser = Series(idx)
    expected = idx.copy(deep=True)
    assert np.shares_memory(get_array(ser), get_array(idx))
    assert not ser._mgr._has_no_reference(0)
    ser.iloc[0] = ser.iloc[1]
    tm.assert_index_equal(idx, expected)

    # forcing copy=False still gives a CoW shallow copy
    ser = Series(idx, copy=False)
    assert np.shares_memory(get_array(ser), get_array(idx))
    assert not ser._mgr._has_no_reference(0)
    ser.iloc[0] = ser.iloc[1]
    tm.assert_index_equal(idx, expected)

    # forcing copy=True still results in a copy
    ser = Series(idx, copy=True)
    assert not np.shares_memory(get_array(ser), get_array(idx))
    assert ser._mgr._has_no_reference(0)


@pytest.mark.parametrize("copy", [True, False, None])
def test_series_from_index_different_dtypes(copy):
    idx = Index([1, 2, 3], dtype="int64", copy=copy)
    ser = Series(idx, dtype="int32")
    assert not np.shares_memory(get_array(ser), get_array(idx))
    assert ser._mgr._has_no_reference(0)


def test_series_from_block_manager_different_dtype():
    ser = Series([1, 2, 3], dtype="int64")
    msg = "Passing a SingleBlockManager to Series"
    with tm.assert_produces_warning(DeprecationWarning, match=msg):
        ser2 = Series(ser._mgr, dtype="int32")
    assert not np.shares_memory(get_array(ser), get_array(ser2))
    assert ser2._mgr._has_no_reference(0)


@pytest.mark.parametrize("use_mgr", [True, False])
@pytest.mark.parametrize("columns", [None, ["a"]])
def test_dataframe_constructor_mgr_or_df(columns, use_mgr):
    df = DataFrame({"a": [1, 2, 3]})
    df_orig = df.copy()

    if use_mgr:
        data = df._mgr
        warn = DeprecationWarning
    else:
        data = df
        warn = None
    msg = "Passing a BlockManager to DataFrame"
    with tm.assert_produces_warning(warn, match=msg, check_stacklevel=False):
        new_df = DataFrame(data)

    assert np.shares_memory(get_array(df, "a"), get_array(new_df, "a"))
    new_df.iloc[0] = 100

    assert not np.shares_memory(get_array(df, "a"), get_array(new_df, "a"))
    tm.assert_frame_equal(df, df_orig)


@pytest.mark.parametrize("dtype", [None, "int64", "Int64"])
@pytest.mark.parametrize("index", [None, [0, 1, 2]])
@pytest.mark.parametrize("columns", [None, ["a", "b"], ["a", "b", "c"]])
def test_dataframe_from_dict_of_series(columns, index, dtype):
    # Case: constructing a DataFrame from Series objects with copy=False
    # has to do a lazy following CoW rules
    # (the default for DataFrame(dict) is still to copy to ensure consolidation)
    s1 = Series([1, 2, 3])
    s2 = Series([4, 5, 6])
    s1_orig = s1.copy()
    expected = DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6]}, index=index, columns=columns, dtype=dtype
    )

    result = DataFrame(
        {"a": s1, "b": s2}, index=index, columns=columns, dtype=dtype, copy=False
    )

    # the shallow copy still shares memory
    assert np.shares_memory(get_array(result, "a"), get_array(s1))

    # mutating the new dataframe doesn't mutate original
    result.iloc[0, 0] = 10
    assert not np.shares_memory(get_array(result, "a"), get_array(s1))
    tm.assert_series_equal(s1, s1_orig)

    # the same when modifying the parent series
    s1 = Series([1, 2, 3])
    s2 = Series([4, 5, 6])
    result = DataFrame(
        {"a": s1, "b": s2}, index=index, columns=columns, dtype=dtype, copy=False
    )
    s1.iloc[0] = 10
    assert not np.shares_memory(get_array(result, "a"), get_array(s1))
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("dtype", [None, "int64"])
def test_dataframe_from_dict_of_series_with_reindex(dtype):
    # Case: constructing a DataFrame from Series objects with copy=False
    # and passing an index that requires an actual (no-view) reindex -> need
    # to ensure the result doesn't have refs set up to unnecessarily trigger
    # a copy on write
    s1 = Series([1, 2, 3])
    s2 = Series([4, 5, 6])
    df = DataFrame({"a": s1, "b": s2}, index=[1, 2, 3], dtype=dtype, copy=False)

    # df should own its memory, so mutating shouldn't trigger a copy
    arr_before = get_array(df, "a")
    assert not np.shares_memory(arr_before, get_array(s1))
    df.iloc[0, 0] = 100
    arr_after = get_array(df, "a")
    assert np.shares_memory(arr_before, arr_after)


@pytest.mark.parametrize(
    "data, dtype",
    [
        ([1, 2], "int64"),
        # 1D-only EA
        ([1, 2], "Int64"),
        pytest.param(
            ["a", "b"],
            "str",
            marks=pytest.mark.xfail(
                reason="TODO bug with infer_string=False and specifying dtype='str'"
            )
            if not using_string_dtype()
            else [],
        ),
        (["a", "b"], object),
        # 2D EA
        (
            [Timestamp("2020", tz="UTC"), Timestamp("2021", tz="UTC")],
            "datetime64[ns, UTC]",
        ),
    ],
    ids=["int", "int-ea", "str", "object", "datetime64tz"],
)
def test_dataframe_from_series_or_index(data, dtype, index_or_series):
    obj = index_or_series(data, dtype=dtype)
    obj_orig = obj.copy(deep=True)  # deep=True needed for Index

    # default is copy=False -> DataFrame holds a shallow copy of original Index/Series
    df = DataFrame(obj)
    assert tm.shares_memory(get_array(obj), get_array(df, 0))
    assert not df._mgr._has_no_reference(0)

    df.iloc[0, 0] = data[-1]
    tm.assert_equal(obj, obj_orig)

    # with passing the (identical) dtype -> same
    df = DataFrame(obj, dtype=dtype)
    assert tm.shares_memory(get_array(obj), get_array(df, 0))
    assert not df._mgr._has_no_reference(0)

    df.iloc[0, 0] = data[-1]
    tm.assert_equal(obj, obj_orig)

    # forcing copy=True still results in an actual hard copy up front
    df = DataFrame(obj, copy=True)
    if not (obj.dtype == "str" and obj.dtype.storage == "pyarrow"):
        # ArrowExtensionArray deep copy still points to the same underlying data
        assert not tm.shares_memory(get_array(obj), get_array(df, 0))
        assert df._mgr._has_no_reference(0)

    df.iloc[0, 0] = data[-1]
    tm.assert_equal(obj, obj_orig)


def test_dataframe_from_series_or_index_different_dtype(index_or_series):
    obj = index_or_series([1, 2], dtype="int64")
    df = DataFrame(obj, dtype="int32")
    assert not np.shares_memory(get_array(obj), get_array(df, 0))
    assert df._mgr._has_no_reference(0)


def test_dataframe_from_series_dont_infer_datetime():
    ser = Series([Timestamp("2019-12-31"), Timestamp("2020-12-31")], dtype=object)
    df = DataFrame(ser)
    assert df.dtypes.iloc[0] == np.dtype(object)
    assert np.shares_memory(get_array(ser), get_array(df, 0))
    assert not df._mgr._has_no_reference(0)


@pytest.mark.parametrize("index", [None, [0, 1, 2]])
def test_dataframe_from_dict_of_series_with_dtype(index):
    # Variant of above, but now passing a dtype that causes a copy
    # -> need to ensure the result doesn't have refs set up to unnecessarily
    # trigger a copy on write
    s1 = Series([1.0, 2.0, 3.0])
    s2 = Series([4, 5, 6])
    df = DataFrame({"a": s1, "b": s2}, index=index, dtype="int64", copy=False)

    # df should own its memory, so mutating shouldn't trigger a copy
    arr_before = get_array(df, "a")
    assert not np.shares_memory(arr_before, get_array(s1))
    df.iloc[0, 0] = 100
    arr_after = get_array(df, "a")
    assert np.shares_memory(arr_before, arr_after)


@pytest.mark.parametrize("copy", [False, None, True])
def test_dataframe_from_numpy_array(copy):
    arr = np.array([[1, 2], [3, 4]])
    df = DataFrame(arr, copy=copy)

    if copy is not False or copy is True:
        assert not np.shares_memory(get_array(df, 0), arr)
    else:
        assert np.shares_memory(get_array(df, 0), arr)


@pytest.mark.parametrize(
    "data, dtype",
    [
        # 1D-only EA
        ([1, 2], "Int64"),
        # 2D EA
        (
            [Timestamp("2020", tz="UTC"), Timestamp("2021", tz="UTC")],
            "datetime64[ns, UTC]",
        ),
    ],
    ids=["int-ea", "datetime64tz"],
)
@pytest.mark.parametrize("copy", [False, None, True])
def test_dataframe_from_extension_array(copy, data, dtype):
    arr = pd.array(data, dtype=dtype)
    df = DataFrame(arr, copy=copy)

    if arr.dtype == "Int64":
        # to ensure tm.shares_memory works correctly
        # TODO fix in tm.shares_memory or get_array?
        arr = arr._data

    if copy is None or copy is True:
        assert not tm.shares_memory(get_array(df, 0), arr)
    else:
        assert tm.shares_memory(get_array(df, 0), arr)


def test_frame_from_dict_of_index():
    idx = Index([1, 2, 3])
    expected = idx.copy(deep=True)
    df = DataFrame({"a": idx}, copy=False)
    assert np.shares_memory(get_array(df, "a"), idx._values)
    assert not df._mgr._has_no_reference(0)

    df.iloc[0, 0] = 100
    tm.assert_index_equal(idx, expected)
