"""
Tests for sharing one sorted copy of the data between group-by-group operations.
"""

import gc
import weakref

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


def _frame():
    return pd.DataFrame(
        {"k": [1, 0, 1, 0], "x": [1.0, 2.0, 3.0, 4.0], "y": [5.0, 6.0, 7.0, 8.0]}
    )


@pytest.mark.parametrize(
    "first, second",
    [
        (lambda gb, f: gb["x"].agg(f), lambda gb, f: gb["x"].agg(f)),
        (lambda gb, f: gb["x"].agg([f, f]), lambda gb, f: gb["x"].apply(f)),
        (lambda gb, f: gb["x"].agg(f), lambda gb, f: gb["x"].transform(f)),
        (lambda gb, f: gb.agg(f), lambda gb, f: gb.agg(f)),
    ],
)
def test_repeated_calls_share_sorted_data(first, second):
    # GH#51077
    gb = _frame().groupby("k")
    seen = []

    def record(group):
        seen.append(group.to_numpy())
        return 0.0

    first(gb, record)
    n_first = len(seen)
    second(gb, record)
    assert np.shares_memory(seen[0], seen[n_first])


def test_cache_shared_with_selection():
    # GH#51077 the selection is a new Series on each access
    gb = _frame().groupby("k")
    gb["x"].agg(lambda x: 0.0)
    gb["x"].agg(lambda x: 0.0)
    gb["y"].agg(lambda x: 0.0)
    assert len(gb._grouper._splitter_cache) == 2


@pytest.mark.parametrize(
    "values",
    [
        np.array([1.0, 2.0, 3.0, 4.0]),
        pd.date_range("2016-01-01", periods=4),
    ],
)
def test_duplicate_names_in_one_block_not_confused(values):
    # GH#51077
    df = pd.DataFrame({"a": values, "b": values[::-1]})
    df.columns = ["a", "a"]
    assert len(df._mgr.blocks) == 1
    gb = df.groupby(np.array([1, 0, 1, 0]))
    expected = gb.first()

    for _ in range(2):
        result = gb.agg(lambda x: x.iloc[0])
        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "write",
    [
        lambda ser: ser.iloc.__setitem__(0, 100),
        lambda ser: ser.loc.__setitem__(0, 100),
        lambda ser: ser.__iadd__(100),
        lambda ser: ser.mask(ser == ser.iloc[0], 100, inplace=True),
    ],
)
@pytest.mark.parametrize("dtype", ["float64", "Int64", "Float64", "int64[pyarrow]"])
def test_series_source_write_seen(write, dtype):
    # GH#51077
    if dtype == "int64[pyarrow]":
        pytest.importorskip("pyarrow")
    ser = pd.Series([1, 2, 3, 4], dtype=dtype)
    keys = np.array([1, 0, 1, 0])
    gb = ser.groupby(keys)
    func = lambda x: x.iloc[0]
    gb.agg(func)

    write(ser)

    tm.assert_series_equal(gb.agg(func), ser.groupby(keys).agg(func))


@pytest.mark.parametrize(
    "write",
    [
        lambda df: df.loc.__setitem__((0, "x"), 100.0),
        lambda df: df.iloc.__setitem__((0, 0), 100.0),
        lambda df: df.__setitem__("x", 100.0),
        lambda df: df.__setitem__("z", 100.0),
    ],
)
def test_frame_source_write_seen(write):
    # GH#51077 keys come from outside, so nothing is excluded from the frame
    df = _frame().drop(columns="k")
    keys = np.array([1, 0, 1, 0])
    gb = df.groupby(keys)
    func = lambda frame: frame.iloc[0, 0]
    gb.apply(func)

    write(df)

    tm.assert_series_equal(gb.apply(func), df.groupby(keys).apply(func))


def test_source_relabel_seen():
    # GH#51077
    ser = pd.Series([1.0, 2.0, 3.0, 4.0], name="a")
    keys = np.array([1, 0, 1, 0])
    gb = ser.groupby(keys)
    seen = []

    def record(group):
        seen.append(group)
        return 0.0

    gb.agg(record)
    ser.name = "b"
    gb.agg(record)
    assert seen[-1].name == "b"

    ser.index = list("wxyz")
    gb.agg(record)
    assert list(seen[-1].index) == ["w", "y"]

    ser.index.name = "idx"
    gb.agg(record)
    assert seen[-1].index.name == "idx"


def test_source_attrs_seen():
    # GH#51077
    ser = pd.Series([1.0, 2.0, 3.0, 4.0])
    gb = ser.groupby([1, 0, 1, 0])
    seen = []

    def record(group):
        seen.append(group)
        return 0.0

    gb.agg(record)
    ser.attrs["a"] = 1
    gb.agg(record)
    assert seen[-1].attrs == {"a": 1}


def test_sorted_copy_freed_with_groupby():
    # GH#51077
    gb = _frame().groupby("k")
    gb["x"].agg(lambda x: 0.0)
    sorted_data = next(iter(gb._grouper._splitter_cache.values()))[1]._sorted_data
    ref = weakref.ref(sorted_data)
    del sorted_data, gb
    gc.collect()
    assert ref() is None


def test_stale_entry_dropped():
    # GH#51077
    ser = pd.Series([1.0, 2.0, 3.0, 4.0])
    gb = ser.groupby([1, 0, 1, 0])
    for value in [100.0, 200.0, 300.0]:
        gb.agg(lambda x: 0.0)
        ser.iloc[0] = value
    assert len(gb._grouper._splitter_cache) == 1


def test_sorted_input_not_cached():
    # GH#51077 no copy to retain, and no reason to make writes to the source copy
    gb = _frame().sort_values("k", kind="stable").groupby("k")
    gb.agg(lambda x: 0.0)
    assert gb._grouper._splitter_cache == {}


def test_group_setitem_not_shared():
    # GH#51077 Copy-on-Write keeps setitem on a group out of the shared sorted copy
    def mutate(group):
        group.iloc[0] = 100.0
        return group.iloc[0]

    ser = pd.Series([3.0, 2.0, 1.0, 0.0])
    result = ser.groupby([1, 0, 1, 0]).agg([mutate, lambda x: x.iloc[0]])
    expected = pd.DataFrame(
        {"mutate": [100.0, 100.0], "<lambda_0>": [2.0, 3.0]}, index=[0, 1]
    )
    tm.assert_frame_equal(result, expected)


def test_other_column_in_udf():
    # GH#51077 a nested call on the same grouper must not get another column's data
    df = pd.DataFrame(
        {"k": [1, 0, 1, 0], "a": [1.0, 2.0, 3.0, 4.0], "b": [10.0, 20.0, 30.0, 40.0]}
    )
    gb = df.groupby("k")

    def other(group):
        return gb["b"].agg(lambda x: x.iloc[0]).sum()

    result = gb["a"].agg([other, lambda x: x.iloc[0]])
    expected = pd.DataFrame(
        {"other": [30.0, 30.0], "<lambda_0>": [2.0, 1.0]},
        index=pd.Index([0, 1], name="k"),
    )
    tm.assert_frame_equal(result, expected)
