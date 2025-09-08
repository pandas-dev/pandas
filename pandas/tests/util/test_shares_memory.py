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


# Unit tests for tm.shares_memory (#55372)
def test_shares_memory_numpy():
    arr = np.arange(10)
    view = arr[:5]
    assert tm.shares_memory(arr, view)
    arr2 = np.arange(10)
    assert not tm.shares_memory(arr, arr2)


def test_shares_memory_series():
    arr = np.arange(10)
    s = pd.Series(arr)
    assert tm.shares_memory(arr, s)
    s2 = pd.Series(np.arange(10))
    assert not tm.shares_memory(s, s2)


def test_shares_memory_dataframe_single_block():
    arr = np.arange(10)
    df = pd.DataFrame({"a": arr})
    assert tm.shares_memory(arr, df)
    df2 = pd.DataFrame({"a": np.arange(10)})
    assert not tm.shares_memory(df, df2)


def test_shares_memory_rangeindex():
    idx = pd.RangeIndex(10)
    arr = np.arange(10)
    assert not tm.shares_memory(idx, arr)


def test_shares_memory_multiindex():
    idx = pd.MultiIndex.from_arrays([np.arange(10), np.arange(10, 20)])
    arr = idx.codes[0]
    assert tm.shares_memory(idx, arr)
    arr2 = np.arange(10)
    assert not tm.shares_memory(idx, arr2)
