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


@pytest.fixture(
    params=["float64", "Int64", "int64[pyarrow]", pd.StringDtype(na_value=np.nan)],
    ids=str,
)
def ser(request):
    dtype = request.param
    if dtype == "int64[pyarrow]":
        pytest.importorskip("pyarrow")
    values = list("1234") if isinstance(dtype, pd.StringDtype) else [1, 2, 3, 4]
    return pd.Series(values, dtype=dtype, name="a")


@pytest.mark.parametrize(
    "dtype",
    ["float64", "int64[pyarrow]", pd.StringDtype("pyarrow", na_value=np.nan)],
    ids=str,
)
def test_repeated_calls_hit_cache_by_dtype(dtype):
    # GH#51077
    pytest.importorskip("pyarrow")
    values = list("1234") if isinstance(dtype, pd.StringDtype) else [1, 2, 3, 4]
    df = pd.DataFrame({"a": pd.Series(values, dtype=dtype)})
    gb = df.groupby(np.array([1, 0, 1, 0]))
    gb["a"].agg(lambda x: x.iloc[0])
    entries = list(gb._grouper._splitter_cache.values())
    gb["a"].agg(lambda x: x.iloc[0])
    new_entries = list(gb._grouper._splitter_cache.values())
    assert len(new_entries) == 1
    assert new_entries[0][1] is entries[0][1]


@pytest.mark.parametrize(
    "values, dtype",
    [
        ([5, 4, 3, 2, 1, 9], "Int64"),
        ([5.0, 4.0, 3.0, 2.0, 1.0, 9.0], "Float64"),
        ([True, True, False, False, False, False], "boolean"),
        (list("543219"), pd.StringDtype("python", na_value=np.nan)),
        (list("543219"), pd.StringDtype("python", na_value=pd.NA)),
    ],
    ids=str,
)
def test_writable_to_numpy_not_cached(values, dtype):
    # GH#51077 sorting to_numpy in a UDF must not reach later calls, see GH#63099
    ser = pd.Series(values, dtype=dtype)
    keys = [1, 0, 1, 0, 1, 0]
    gb = ser.groupby(keys)
    first = lambda x: x.iloc[0]

    def sort_in_place(group):
        arr = group.to_numpy()
        arr.sort()
        return arr[0]

    gb.agg(sort_in_place)

    assert len(gb._grouper._splitter_cache) == 0
    tm.assert_series_equal(gb.agg(first), ser.groupby(keys).agg(first))


def test_writable_categorical_values_not_cached():
    # GH#51077 sorting Series.values in a UDF must not reach later calls
    ser = pd.Series(pd.Categorical([5, 4, 3, 2, 1, 9]))
    keys = [1, 0, 1, 0, 1, 0]
    gb = ser.groupby(keys)
    first = lambda x: x.iloc[0]

    def sort_in_place(group):
        arr = group.values
        arr.sort()
        return arr[0]

    gb.agg(sort_in_place)

    assert len(gb._grouper._splitter_cache) == 0
    tm.assert_series_equal(gb.agg(first), ser.groupby(keys).agg(first))


def test_masked_source_to_numpy_write_seen():
    # GH#51077
    ser = pd.Series([1, 2, 3, 4], dtype="Int64")
    keys = np.array([1, 0, 1, 0])
    gb = ser.groupby(keys)
    first = lambda x: x.iloc[0]
    gb.agg(first)

    ser.to_numpy()[0] = 100

    tm.assert_series_equal(gb.agg(first), ser.groupby(keys).agg(first))


def test_frame_with_masked_block_not_cached():
    # GH#51077
    df = pd.DataFrame({"a": [1.0, 2.0, 3.0, 4.0]})
    df["b"] = pd.array([1, 2, 3, 4], dtype="Int64")
    gb = df.groupby(np.array([1, 0, 1, 0]))

    gb.apply(lambda group: 0.0)

    assert len(gb._grouper._splitter_cache) == 0


def test_duplicate_names_sharing_refs_not_confused(ser):
    # GH#51077 the second column is a slice of the first
    df = pd.concat([ser, ser.iloc[::-1].reset_index(drop=True)], axis=1)
    gb = df.groupby(np.array([1, 0, 1, 0]))
    result = gb.agg(lambda x: x.iloc[0])
    tm.assert_frame_equal(result, gb.first())


def test_reassigned_reversed_column_seen(ser):
    # GH#51077 the reversed column is a slice, so it shares refs with the old one
    df = ser.to_frame()
    keys = np.array([1, 0, 1, 0])
    gb = df.groupby(keys)
    func = lambda x: x.iloc[0]
    gb["a"].agg(func)

    df["a"] = df["a"].iloc[::-1].reset_index(drop=True)

    tm.assert_series_equal(gb["a"].agg(func), df.groupby(keys)["a"].agg(func))


@pytest.mark.parametrize(
    "values, new_dtype",
    [
        (
            pd.date_range("2016-01-01", periods=4, tz="UTC"),
            "datetime64[us, Asia/Tokyo]",
        ),
        (
            pd.array(list("wxyz"), dtype=pd.StringDtype(na_value=np.nan)),
            pd.StringDtype(na_value=pd.NA),
        ),
        (pd.Categorical(list("wxyz")), pd.CategoricalDtype(ordered=True)),
    ],
)
def test_dtype_change_sharing_data_seen(values, new_dtype):
    # GH#51077
    df = pd.DataFrame({"a": values})
    gb = df.groupby(np.array([1, 0, 1, 0]))
    seen = []

    def record(group):
        seen.append(group.dtype)
        return 0.0

    gb["a"].agg(record)
    df["a"] = df["a"].astype(new_dtype)
    gb["a"].agg(record)
    assert seen[-1] == df["a"].dtype


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


def test_nan_label_not_cached():
    # GH#51077 otherwise each call adds an entry
    df = pd.DataFrame(np.arange(8.0).reshape(4, 2), columns=[1.0, np.nan])
    gb = df.groupby(np.array([1, 0, 1, 0]))
    for _ in range(3):
        result = gb.agg(lambda x: x.iloc[0])
    tm.assert_frame_equal(result, gb.first())
    assert len(gb._grouper._splitter_cache) == 1


def test_nan_in_multiindex_label_not_cached():
    # GH#51077
    columns = pd.MultiIndex.from_tuples([("a", 1.0), ("a", np.nan)])
    df = pd.DataFrame(np.arange(8.0).reshape(4, 2), columns=columns)
    gb = df.groupby(np.array([1, 0, 1, 0]))
    for _ in range(3):
        gb.apply(lambda x: 0.0)
    assert len(gb._grouper._splitter_cache) == 0


def test_unhashable_label_not_cached():
    # GH#51077
    df = pd.DataFrame(np.arange(8.0).reshape(4, 2))
    df.columns = pd.Index([[1], [2]], tupleize_cols=False)
    result = df.groupby(np.array([1, 0, 1, 0])).apply(lambda x: x.iloc[0, 0])
    expected = pd.Series([2.0, 0.0], index=np.array([0, 1]))
    tm.assert_series_equal(result, expected)


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
        {"mutate": [100.0, 100.0], "<lambda_0>": [2.0, 3.0]}, index=np.array([0, 1])
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
