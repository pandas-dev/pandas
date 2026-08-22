import numpy as np

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm


def test_shares_memory_interval():
    obj = pd.interval_range(1, 5)

    assert tm.shares_memory(obj, obj)
    assert tm.shares_memory(obj, obj._data)
    assert tm.shares_memory(obj, obj[::-1])
    assert tm.shares_memory(obj, obj[:2])

    assert not tm.shares_memory(obj, obj._data.copy())


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


def test_shares_memory_numpy():
    arr = np.arange(10)
    view = arr[:5]
    assert tm.shares_memory(arr, view)
    arr2 = np.arange(10)
    assert not tm.shares_memory(arr, arr2)


def test_shares_memory_series_index_and_frame():
    arr = np.arange(4)
    ser = pd.Series(arr, copy=False)
    idx = pd.Index(np.arange(4))
    df = pd.DataFrame(arr.reshape(2, 2), copy=False)

    assert tm.shares_memory(ser, arr)
    assert tm.shares_memory(idx, idx[:2])
    assert tm.shares_memory(df, arr)

    arr2 = np.arange(4)
    assert not tm.shares_memory(ser, arr2)
    assert not tm.shares_memory(idx, arr2)
    assert not tm.shares_memory(df, arr2)


def test_shares_memory_masked_array():
    obj = pd.array([1, None], dtype="Int64")

    assert tm.shares_memory(obj, obj[:])
    assert not tm.shares_memory(obj, obj.copy())


def test_shares_memory_sparse_array():
    obj = pd.arrays.SparseArray([0, 1, 0])

    assert tm.shares_memory(obj, obj)
    assert not tm.shares_memory(obj, obj.copy())


def test_shares_memory_rangeindex():
    idx = pd.RangeIndex(10)
    arr = np.arange(10)
    assert not tm.shares_memory(idx, arr)
