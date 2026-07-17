import numpy as np
import pytest

from pandas import (
    DataFrame,
    Index,
    Interval,
    Series,
    array,
    interval_range,
)
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


def index_view(index_data):
    df = DataFrame({"a": index_data, "b": 1.5})
    view = df[:]
    df = df.set_index("a", drop=True)
    idx = df.index
    # df = None
    return idx, view


def test_set_index_update_column():
    df = DataFrame({"a": [1, 2], "b": 1})
    df = df.set_index("a", drop=False)
    expected = df.index.copy(deep=True)
    df.iloc[0, 0] = 100
    tm.assert_index_equal(df.index, expected)


def test_set_index_drop_update_column():
    df = DataFrame({"a": [1, 2], "b": 1.5})
    view = df[:]
    df = df.set_index("a", drop=True)
    expected = df.index.copy(deep=True)
    view.iloc[0, 0] = 100
    tm.assert_index_equal(df.index, expected)


def test_set_index_series():
    df = DataFrame({"a": [1, 2], "b": 1.5})
    ser = Series([10, 11])
    df = df.set_index(ser)
    expected = df.index.copy(deep=True)
    ser.iloc[0] = 100
    tm.assert_index_equal(df.index, expected)


def test_assign_index_as_series():
    df = DataFrame({"a": [1, 2], "b": 1.5})
    ser = Series([10, 11])
    df.index = ser
    expected = df.index.copy(deep=True)
    ser.iloc[0] = 100
    tm.assert_index_equal(df.index, expected)


def test_assign_index_as_index():
    df = DataFrame({"a": [1, 2], "b": 1.5})
    ser = Series([10, 11])
    rhs_index = Index(ser)
    df.index = rhs_index
    rhs_index = None  # overwrite to clear reference
    expected = df.index.copy(deep=True)
    ser.iloc[0] = 100
    tm.assert_index_equal(df.index, expected)


def test_index_from_series():
    ser = Series([1, 2])
    idx = Index(ser)
    expected = idx.copy(deep=True)
    ser.iloc[0] = 100
    tm.assert_index_equal(idx, expected)


def test_index_from_series_copy():
    ser = Series([1, 2])
    idx = Index(ser, copy=True)  # noqa: F841
    arr = get_array(ser)
    ser.iloc[0] = 100
    assert np.shares_memory(get_array(ser), arr)


def test_index_from_index():
    ser = Series([1, 2])
    idx = Index(ser)
    idx = Index(idx)
    expected = idx.copy(deep=True)
    ser.iloc[0] = 100
    tm.assert_index_equal(idx, expected)


@pytest.mark.parametrize(
    "func",
    [
        lambda x: x._shallow_copy(x._values),
        lambda x: x.view(),
        lambda x: x.take([0, 1]),
        lambda x: x.repeat([1, 1]),
        lambda x: x[slice(0, 2)],
        lambda x: x[[0, 1]],
        lambda x: x._getitem_slice(slice(0, 2)),
        lambda x: x.delete([]),
        lambda x: x.rename("b"),
        lambda x: x.astype("Int64", copy=False),
    ],
    ids=[
        "_shallow_copy",
        "view",
        "take",
        "repeat",
        "getitem_slice",
        "getitem_list",
        "_getitem_slice",
        "delete",
        "rename",
        "astype",
    ],
)
def test_index_ops(func, request):
    idx, view_ = index_view([1, 2])
    expected = idx.copy(deep=True)
    if "astype" in request.node.callspec.id:
        expected = expected.astype("Int64")
    idx = func(idx)
    view_.iloc[0, 0] = 100
    tm.assert_index_equal(idx, expected, exact="equiv", check_names=False)


def test_infer_objects():
    idx, view_ = index_view(["a", "b"])
    expected = idx.copy(deep=True)
    idx = idx.infer_objects(copy=False)
    view_.iloc[0, 0] = "aaaa"
    tm.assert_index_equal(idx, expected, check_names=False)


def test_index_to_frame():
    idx = Index([1, 2, 3], name="a")
    expected = idx.copy(deep=True)
    df = idx.to_frame()
    assert np.shares_memory(get_array(df, "a"), idx._values)
    assert not df._mgr._has_no_reference(0)

    df.iloc[0, 0] = 100
    tm.assert_index_equal(idx, expected)


def test_index_where_noop():
    # GH#65265 CoW references should be tracked through Index.where
    idx = Index([1, 2, 3])
    result = idx.where([True] * 3, 100)
    assert np.shares_memory(get_array(idx), get_array(result))
    assert result._references.has_reference()

    expected = idx.copy(deep=True)
    result = Series(result)
    result.iloc[0] = 100
    tm.assert_index_equal(idx, expected)


@pytest.mark.parametrize(
    "data, dtype",
    [
        ([1, 2, 3], None),
        (["a", "b", "c"], None),
        ([1, 2, 3], "Int64"),
        (["a", "b", "c"], "category"),
    ],
)
def test_index_values(data, dtype):
    # GH#38547 mutating the array returned by Index.values must not be able
    #  to silently corrupt the (immutable) Index
    idx = Index(data, dtype=dtype)
    result = idx.values
    if isinstance(result, np.ndarray):
        assert result.flags.writeable is False
    else:
        assert result._readonly is True
    with pytest.raises(ValueError, match="read-only"):
        result[0] = data[1]


@pytest.mark.parametrize(
    "data, dtype",
    [
        ([1, 2, 3], None),
        (["a", "b", "c"], None),
        ([1, 2, 3], "Int64"),
        (["a", "b", "c"], "category"),
    ],
)
def test_index_array_readonly(data, dtype):
    # GH#38547 Index.array is read-only so the immutable Index cannot be
    #  mutated through it
    idx = Index(data, dtype=dtype)
    result = idx.array
    assert result._readonly is True
    with pytest.raises(ValueError, match="read-only"):
        result[0] = data[1]


@pytest.mark.parametrize(
    "data, dtype",
    [
        (["a", "b", "a"], "category"),
        (["2020-01-01", "2020-01-02"], "M8[ns]"),
        (["2020-01-01", "2020-01-02"], "M8[ns, UTC]"),
        (["1 day", "2 days"], "m8[ns]"),
        (["2020-01-01", "2020-01-02"], "period[D]"),
        ([1, 2, 3], "Int64"),
    ],
)
def test_index_array_inplace_op_raises(data, dtype):
    # GH#38547 in-place Series ops must not be able to write into the
    #  zero-copy readonly array escaping through Index.array
    idx = Index(data, dtype=dtype)
    expected = idx.copy(deep=True)
    ser = Series(idx.array, copy=False)

    msg = "Cannot modify read-only array"
    with pytest.raises(ValueError, match=msg):
        ser.iloc[0] = idx[1]
    with pytest.raises(ValueError, match=msg):
        ser.mask(ser == idx[0], idx[1], inplace=True)

    tm.assert_index_equal(idx, expected)


def test_index_array_inplace_op_raises_interval():
    # GH#38547 IntervalArray._putmask has its own path mutating _left/_right
    idx = interval_range(0, 2)
    expected = idx.copy(deep=True)
    ser = Series(idx.array, copy=False)

    msg = "Cannot modify read-only array"
    with pytest.raises(ValueError, match=msg):
        ser.iloc[0] = Interval(8, 9)
    with pytest.raises(ValueError, match=msg):
        ser.mask(ser == idx[0], Interval(8, 9), inplace=True)

    tm.assert_index_equal(idx, expected)


def test_columns_values_setitem_raises():
    # GH#38547 directly mutating df.columns.values used to silently corrupt
    #  the columns Index (it appeared renamed but lookups failed)
    df = DataFrame(np.eye(2), columns=["a", "b"])
    msg = "read-only"
    with pytest.raises(ValueError, match=msg):
        df.columns.values[0] = "x"
    with pytest.raises(ValueError, match=msg):
        df.columns.array[0] = "x"
    # the Index is left intact and still usable under its original labels
    tm.assert_index_equal(df.columns, Index(["a", "b"]))
    tm.assert_series_equal(df["a"], Series([1.0, 0.0], name="a"))


def test_constructor_copy_input_ndarray_default():
    arr = np.array([0, 1])
    idx = Index(arr)
    assert not np.shares_memory(arr, get_array(idx))


def test_constructor_copy_input_ea_default():
    arr = array([0, 1], dtype="Int64")
    idx = Index(arr)
    assert not tm.shares_memory(arr, idx.array)


def test_series_from_temporary_index_readonly_data():
    # GH 63370
    arr = np.array([0, 1], dtype=np.dtype(np.int8))
    arr.flags.writeable = False
    ser = Series(Index(arr))
    assert not np.shares_memory(arr, get_array(ser))
    assert ser._mgr._has_no_reference(0)
    ser[[False, True]] = np.array([0, 2], dtype=np.dtype(np.int8))
    expected = Series([0, 2], dtype=np.dtype(np.int8))
    tm.assert_series_equal(ser, expected)
