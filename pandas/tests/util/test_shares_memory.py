import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm


def test_shares_memory_numpy():
    arr = np.arange(10)
    view = arr[:5]
    assert tm.shares_memory(arr, view)
    arr2 = np.arange(10)
    assert not tm.shares_memory(arr, arr2)


def test_shares_memory_ndarray_and_ea():
    # GH#55372 the ndarray may be passed as either argument
    obj = pd.Categorical(["a", "b", "a"])
    assert tm.shares_memory(obj._ndarray, obj)
    assert tm.shares_memory(obj, obj._ndarray)

    assert not tm.shares_memory(np.arange(3, dtype=np.int8), obj)


def test_shares_memory_rangeindex():
    idx = pd.RangeIndex(10)
    arr = np.arange(10)
    assert not tm.shares_memory(idx, arr)
    assert not tm.shares_memory(idx, idx)


def test_shares_memory_multiindex():
    # GH#55372
    mi = pd.MultiIndex.from_product([[1, 2], ["a", "b"]])

    assert tm.shares_memory(mi, mi)
    for codes in mi.codes:
        assert tm.shares_memory(mi, codes)

    assert not tm.shares_memory(mi, mi.copy(deep=True))
    # the levels are not what backs a MultiIndex's values
    assert not tm.shares_memory(mi, mi.levels[0])


def test_shares_memory_index_and_series():
    ser = pd.Series(np.arange(10))
    idx = pd.Index(ser._values, copy=False)

    assert tm.shares_memory(ser, idx)
    assert tm.shares_memory(idx, ser)
    assert tm.shares_memory(ser, ser._values)

    assert not tm.shares_memory(ser, ser.copy(deep=True))


@pytest.mark.parametrize(
    "obj",
    [
        pd.Categorical(["a", "b", "a"]),
        pd.array(["a", "b"], dtype=pd.StringDtype("python")),
        pd.date_range("2016-01-01", periods=3, tz="US/Pacific")._data,
        pd.period_range("2016-01-01", periods=3, freq="D")._data,
        pd.timedelta_range("1 Day", periods=3)._data,
    ],
)
def test_shares_memory_ndarray_backed(obj):
    # GH#55372
    assert tm.shares_memory(obj, obj)
    assert tm.shares_memory(obj, obj[:2])
    assert tm.shares_memory(obj, obj._ndarray)

    assert not tm.shares_memory(obj, obj.copy())


def test_shares_memory_sparse():
    # GH#55372
    obj = pd.arrays.SparseArray([0, 1, 0, 2])

    assert tm.shares_memory(obj, obj)
    assert tm.shares_memory(obj, obj.sp_values)

    assert not tm.shares_memory(obj, obj.copy())
    # the sparse index is not memory we track
    assert not tm.shares_memory(obj, obj.sp_index.indices)


def test_shares_memory_interval():
    obj = pd.interval_range(1, 5)

    assert tm.shares_memory(obj, obj)
    assert tm.shares_memory(obj, obj._data)
    assert tm.shares_memory(obj, obj[::-1])
    assert tm.shares_memory(obj, obj[:2])

    assert not tm.shares_memory(obj, obj._data.copy())


def test_shares_memory_interval_one_side():
    # GH#55372 sharing either side counts as sharing
    obj = pd.arrays.IntervalArray.from_breaks(np.arange(5))
    other = pd.arrays.IntervalArray.from_arrays(obj.left.copy(), obj.right)

    assert tm.shares_memory(obj, other)
    assert tm.shares_memory(obj, obj.right)


def test_shares_memory_masked():
    # GH#55372
    obj = pd.array([1, 2, None], dtype="Int64")

    assert tm.shares_memory(obj, obj)
    assert not tm.shares_memory(obj, obj.copy())

    # a masked array can be compared against something other than a masked array
    assert tm.shares_memory(obj, obj._data)
    assert tm.shares_memory(obj, obj._mask)
    assert tm.shares_memory(obj._data, obj)
    assert not tm.shares_memory(obj, np.arange(3))

    assert tm.shares_memory(pd.Series(obj, copy=False), obj._data)


def test_shares_memory_masked_shared_mask_only():
    # GH#55372 sharing either the values or the mask counts as sharing
    obj = pd.array([1, 2, None], dtype="Int64")
    other = pd.arrays.IntegerArray(obj._data.copy(), obj._mask)

    assert not tm.shares_memory(obj._data, other._data)
    assert tm.shares_memory(obj, other)


@td.skip_if_no("pyarrow")
def test_shares_memory_string():
    # GH#55823
    import pyarrow as pa

    obj = pd.array(["a", "b"], dtype=pd.StringDtype("pyarrow", na_value=pd.NA))
    assert tm.shares_memory(obj, obj)

    obj = pd.array(["a", "b"], dtype=pd.StringDtype("pyarrow", na_value=np.nan))
    assert tm.shares_memory(obj, obj)

    obj = pd.array(["a", "b"], dtype=pd.ArrowDtype(pa.string()))
    assert tm.shares_memory(obj, obj)


@td.skip_if_no("pyarrow")
def test_shares_memory_arrow():
    obj = pd.array([1, 2, 3], dtype="int64[pyarrow]")

    # pyarrow buffers are immutable, so a copy is not a defensive one
    assert tm.shares_memory(obj, obj)
    assert tm.shares_memory(obj, obj.copy())

    assert not tm.shares_memory(obj, pd.array([1, 2, 3], dtype="int64[pyarrow]"))


@td.skip_if_no("pyarrow")
def test_shares_memory_arrow_chunked():
    # GH#55372 sharing in any chunk counts as sharing, not just the first
    import pyarrow as pa

    shared = pa.array([1, 2, 3])
    obj = pd.arrays.ArrowExtensionArray(pa.chunked_array([pa.array([9, 9]), shared]))
    other = pd.arrays.ArrowExtensionArray(pa.chunked_array([pa.array([8, 8]), shared]))

    assert tm.shares_memory(obj, other)


@td.skip_if_no("pyarrow")
@pytest.mark.parametrize(
    "arrow_array",
    [
        lambda pa: pa.chunked_array([], type=pa.int64()),  # no chunks at all
        lambda pa: pa.chunked_array([pa.nulls(3)]),  # null type has no data buffer
        lambda pa: pa.chunked_array([pa.array([], type=pa.int64())]),  # empty buffers
    ],
)
def test_shares_memory_arrow_no_buffers(arrow_array):
    # GH#55372 these used to raise instead of returning a bool
    import pyarrow as pa

    obj = pd.arrays.ArrowExtensionArray(arrow_array(pa))
    other = pd.arrays.ArrowExtensionArray(arrow_array(pa))

    # nothing is actually backed by memory here, so nothing can be shared
    assert not tm.shares_memory(obj, obj)
    assert not tm.shares_memory(obj, other)


@td.skip_if_no("pyarrow")
def test_shares_memory_arrow_nested():
    # GH#55372 nested types have None entries in .buffers()
    import pyarrow as pa

    obj = pd.arrays.ArrowExtensionArray(pa.array([{"a": 1}, {"a": 2}]))

    assert tm.shares_memory(obj, obj)
    assert not tm.shares_memory(
        obj, pd.arrays.ArrowExtensionArray(pa.array([{"a": 3}, {"a": 4}]))
    )


@td.skip_if_no("pyarrow")
def test_shares_memory_arrow_and_numpy():
    # pa.array is zero-copy for a numpy int64 array
    import pyarrow as pa

    arr = np.arange(10, dtype=np.int64)
    obj = pd.arrays.ArrowExtensionArray(pa.array(arr))

    assert tm.shares_memory(obj, arr)
    assert not tm.shares_memory(obj, arr.copy())


def test_shares_memory_dataframe():
    # GH#55372
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    arr = df._mgr.blocks[0].values

    assert tm.shares_memory(df, df)
    assert tm.shares_memory(df, arr)
    assert tm.shares_memory(arr, df)

    assert not tm.shares_memory(df, df.copy(deep=True))


def test_shares_memory_not_implemented():
    # GH#55372 unrecognized objects raise rather than silently returning False
    with pytest.raises(NotImplementedError, match="int"):
        tm.shares_memory(1, 2)

    # more than one block is not supported
    df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
    with pytest.raises(NotImplementedError, match="DataFrame"):
        tm.shares_memory(df, df)
