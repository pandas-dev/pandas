"""
Tests for the deferred gather of string and binary Arrow arrays: ``take``,
slices and boolean masks of such an ``ArrowExtensionArray`` record the
positions or the mask and gather on first access to ``_pa_array``.
"""

from __future__ import annotations

import copy
import pickle
import threading

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")


def _is_pending(arr) -> bool:
    return arr._pending_take is not None


@pytest.fixture(
    params=[pa.string(), pa.large_string(), pa.binary(), pa.large_binary(), "str"],
    ids=["string", "large_string", "binary", "large_binary", "str"],
)
def data(request):
    if request.param == "str":
        values = [f"v{i}" if i % 7 else None for i in range(100)]
        return pd.array(values, dtype=pd.StringDtype("pyarrow", na_value=np.nan))
    return _arrow_array(request.param)


def _arrow_array(pa_type, n=100):
    if pa.types.is_binary(pa_type) or pa.types.is_large_binary(pa_type):
        values = [f"v{i}".encode() if i % 7 else None for i in range(n)]
    else:
        values = [f"v{i}" if i % 7 else None for i in range(n)]
    chunked = pa.chunked_array([values[:40], values[40:]], type=pa_type)
    return pd.arrays.ArrowExtensionArray(chunked)


def _expected(data, positions):
    # the same take, run eagerly on the Arrow array
    return data._from_pyarrow_array(data._pa_array.take(pa.array(positions)))


def test_take_is_deferred(data):
    positions = np.arange(99, 9, -1)
    result = data.take(positions)
    assert _is_pending(result)
    assert len(result) == len(positions)
    assert result.dtype == data.dtype
    # len and dtype do not run the gather
    assert _is_pending(result)
    tm.assert_extension_array_equal(result, _expected(data, positions))
    assert not _is_pending(result)


def test_take_does_not_alias_indexer(data):
    positions = np.arange(50)
    result = data.take(positions)
    positions[:] = 0
    tm.assert_extension_array_equal(result, _expected(data, np.arange(50)))


def test_take_keeps_immutable_indexer(data):
    # an array over ``bytes`` can never be made writeable: kept as it is
    positions = np.frombuffer(np.arange(50, dtype=np.intp).tobytes(), dtype=np.intp)
    result = data.take(positions)
    assert result._pending_take.selection is positions
    result = data.take(positions[10:])
    assert result._pending_take.selection.base is positions


@pytest.mark.parametrize("how", ["owner", "view"])
def test_take_copies_readonly_indexer(data, how):
    # a read-only indexer can still change: its owner can make it writeable
    # again, and a view of it changes with its base
    base = np.arange(50)
    positions = base if how == "owner" else base[:]
    positions.flags.writeable = False
    result = data.take(positions)
    assert result._pending_take.selection is not positions
    base.flags.writeable = True
    base[:] = 0
    tm.assert_extension_array_equal(result, _expected(data, np.arange(50)))


def test_frame_take_shares_indexer_across_columns():
    df = pd.DataFrame(
        {
            "a": pd.array([f"a{i}" for i in range(100)], dtype="str"),
            "b": pd.array([f"b{i}" for i in range(100)], dtype="str"),
            "c": np.arange(100),
        }
    )
    indexer = np.arange(99, -1, -1)
    result = df.take(indexer)
    a = result["a"].array._pending_take
    b = result["b"].array._pending_take
    # one copy of the indexer for the whole frame, not one per column
    assert a.selection is b.selection
    assert a.selection is not indexer
    # the caller's indexer is left writeable, and changing it does not
    # change the result
    assert indexer.flags.writeable
    indexer[:] = 0
    tm.assert_frame_equal(result, df.iloc[::-1])


@pytest.mark.parametrize(
    "select",
    [
        lambda df: df.take([5, 3, 1, 7]),
        lambda df: df[df["c"] % 3 != 0],
        lambda df: df.loc[df["c"] % 3 != 0],
        lambda df: df.iloc[[5, 3, 1, 7]],
        lambda df: df.reindex([5, 3, 1, 7]),
        lambda df: df.sort_values("c"),
        lambda df: df.sort_index(ascending=False),
        lambda df: df.drop_duplicates("d"),
        lambda df: df.query("c > 10", engine="python"),
        lambda df: df.sample(10, random_state=2),
        lambda df: df.nlargest(10, "c"),
        lambda df: df.dropna(subset=["e"]),
        lambda df: df.join(pd.DataFrame({"k": [3, 1]}, index=[5, 1]), how="inner"),
    ],
    ids=[
        "take",
        "mask",
        "loc_mask",
        "iloc",
        "reindex",
        "sort_values",
        "sort_index",
        "dedup",
        "query",
        "sample",
        "nlargest",
        "dropna",
        "join",
    ],
)
def test_frame_row_selection_shares_indexer(select):
    # every row selection of a frame reaches the columns through
    # BlockManager.reindex_indexer, which hands them one read-only indexer
    n = 100
    df = pd.DataFrame(
        {
            "a": pd.array([f"a{i}" for i in range(n)], dtype="str"),
            "b": pd.array([f"b{i}" for i in range(n)], dtype="str"),
            "c": np.arange(n)[::-1].copy(),
            "d": pd.array(np.arange(n) // 2, dtype="Int64"),
            "e": np.where(np.arange(n) % 3 == 0, np.nan, 1.0),
        }
    )
    result = select(df)
    a = result["a"].array._pending_take
    b = result["b"].array._pending_take
    assert a.selection is b.selection
    assert not a.selection.flags.writeable
    expected = select(df.astype({"a": object, "b": object}))
    tm.assert_frame_equal(result, expected.astype({"a": "str", "b": "str"}))


@pytest.mark.parametrize(
    "df, indexer",
    [
        # no column can keep the indexer
        (
            pd.DataFrame(
                {
                    "a": pd.array(range(10), dtype="Int64"),
                    "b": pd.Categorical(list("abcdefghij")),
                    "c": np.arange(10),
                }
            ),
            [3, 1, 2],
        ),
        # a -1 makes every take eager
        (
            pd.DataFrame(
                {"a": pd.array(list("abcdefghij"), dtype="str"), "c": np.arange(10)}
            ),
            [3, -1, 2],
        ),
    ],
    ids=["no_string_column", "missing_rows"],
)
def test_frame_take_copies_indexer_only_if_kept(df, indexer, monkeypatch):
    copies = []
    frombuffer = np.frombuffer

    def spy(*args, **kwargs):
        copies.append(args)
        return frombuffer(*args, **kwargs)

    monkeypatch.setattr(np, "frombuffer", spy)
    df._mgr.reindex_indexer(df.index[indexer], np.array(indexer, dtype=np.intp), axis=1)
    assert not copies


def test_merge_shares_indexer_per_side():
    left = pd.DataFrame(
        {
            "key": np.arange(100),
            "a": pd.array([f"a{i}" for i in range(100)], dtype="str"),
            "b": pd.array([f"b{i}" for i in range(100)], dtype="str"),
        }
    )
    right = pd.DataFrame(
        {
            "key": np.arange(0, 100, 2)[::-1].copy(),
            "c": pd.array([f"c{i}" for i in range(50)], dtype="str"),
            "d": pd.array([f"d{i}" for i in range(50)], dtype="str"),
        }
    )
    result = left.merge(right, on="key")
    pending = {col: result[col].array._pending_take for col in "abcd"}
    # one indexer per side, shared by that side's columns
    assert pending["a"].selection is pending["b"].selection
    assert pending["c"].selection is pending["d"].selection
    expected = pd.merge(
        left.astype({"a": object, "b": object}),
        right.astype({"c": object, "d": object}),
        on="key",
    ).astype(dict.fromkeys("abcd", "str"))
    tm.assert_frame_equal(result, expected)


def _deferred_frame(n=1000):
    # every string column deferred, by a sort that reorders all rows
    df = pd.DataFrame(
        {
            "a": pd.array([f"a{i}" for i in range(n)], dtype="str"),
            "b": pd.array([f"b{i}" for i in range(n)], dtype="str"),
            "c": np.arange(n)[::-1].copy(),
        }
    )
    return df.sort_values("c")


def _gathered(df, col="a") -> bool:
    pending = df[col].array._pending_take
    return pending is None or pending.result is not None


@pytest.mark.parametrize(
    "op",
    [
        lambda df: (len(df), df.shape, df.dtypes),
        lambda df: df[["a", "c"]],
        lambda df: df.assign(z=1),
        lambda df: df.rename(columns={"c": "cc"}),
        lambda df: df.drop(columns="b"),
        lambda df: df["a"].iloc[3],
        lambda df: df.at[df.index[3], "a"],
        lambda df: df[df["c"] > 5],
        lambda df: repr(df),
        lambda df: df.copy(),
    ],
    ids=[
        "shape",
        "select_columns",
        "assign",
        "rename",
        "drop_column",
        "iloc_scalar",
        "at",
        "select_rows",
        "repr",
        "copy",
    ],
)
def test_frame_op_does_not_gather(op):
    df = _deferred_frame()
    op(df)
    assert not _gathered(df)


def test_frame_repr_gathers_only_shown_rows():
    df = _deferred_frame()
    # copy gathers into new memory and leaves ``df`` deferred
    expected = repr(df.copy())
    assert repr(df) == expected
    assert not _gathered(df)


def test_frame_copy_gathers_into_new_memory():
    df = _deferred_frame()
    copied = df.copy()
    assert not _is_pending(copied["a"].array)
    assert not tm.shares_memory(copied["a"].array, df["a"].array)
    tm.assert_frame_equal(copied, df)


@pytest.mark.parametrize(
    "select",
    [
        lambda df: df.merge(
            pd.DataFrame({"c": [3, 1, -1], "k": [0, 1, 2]}), on="c", how="right"
        ),
        lambda df: df.reindex([5, 3, -1]),
    ],
    ids=["merge_with_unmatched", "reindex_new_label"],
)
def test_frame_selection_with_missing_rows_is_eager(select):
    df = pd.DataFrame(
        {
            "a": pd.array([f"a{i}" for i in range(10)], dtype="str"),
            "c": np.arange(10),
        }
    )
    result = select(df)
    assert not _is_pending(result["a"].array)
    assert result["a"].isna().iloc[-1]


def test_compose_take_slice_mask(data):
    first = data.take(np.arange(99, -1, -1))
    second = first.take(np.arange(0, 100, 2))
    third = second[5:45]
    mask = np.zeros(len(third), dtype=bool)
    mask[::3] = True
    mask[1] = True
    fourth = third[mask]
    for arr in [first, second, third, fourth]:
        assert _is_pending(arr)
    # the chain is one deferred take from the original data
    parent, _ = fourth._take_source()
    assert parent is data._pa_array

    expected = data[::-1][::2][5:45][mask]
    tm.assert_extension_array_equal(fourth, expected)
    # gathering the end of the chain leaves the rest deferred
    assert _is_pending(first)
    tm.assert_extension_array_equal(first, data[::-1])


def test_mask_is_deferred(data):
    mask = np.arange(len(data)) % 3 != 0
    result = data[mask]
    assert _is_pending(result)
    tm.assert_extension_array_equal(result, _expected(data, np.flatnonzero(mask)))


def test_mask_is_kept_as_mask(data):
    # gathered with filter, like an eager mask, which keeps the chunking
    mask = np.arange(len(data)) % 3 != 0
    result = data[mask]
    selection = result._pending_take.selection
    assert selection.dtype == np.bool_
    assert selection is not mask
    assert len(result) == mask.sum()
    gathered = result._pa_array
    assert gathered.num_chunks == data._pa_array.num_chunks
    tm.assert_extension_array_equal(result, _expected(data, np.flatnonzero(mask)))


def test_mask_does_not_alias_indexer(data):
    mask = np.arange(len(data)) % 2 == 0
    result = data[mask]
    mask[:] = True
    tm.assert_extension_array_equal(result, data[::2])


def test_mask_of_mask_stays_mask(data):
    first = data[np.arange(len(data)) % 2 == 0]
    second_mask = np.arange(len(first)) % 5 != 0
    second = first[second_mask]
    assert second._pending_take.selection.dtype == np.bool_
    parent, _ = second._take_source()
    assert parent is data._pa_array
    tm.assert_extension_array_equal(second, data[::2][second_mask])
    assert _is_pending(first)


@pytest.mark.parametrize(
    "select",
    [
        lambda arr: arr.take(np.arange(len(arr) - 1, -1, -1)),
        lambda arr: arr[3:30],
        lambda arr: arr[::-2],
    ],
    ids=["take", "slice", "step-slice"],
)
def test_take_or_slice_of_mask_composes_positions(data, select):
    mask = np.arange(len(data)) % 3 != 1
    pending = data[mask]
    result = select(pending)
    assert _is_pending(result)
    assert result._pending_take.selection.dtype == np.intp
    tm.assert_extension_array_equal(result, select(data[mask].copy()))
    # the mask itself is still gathered with filter
    assert pending._pending_take.selection.dtype == np.bool_


@pytest.mark.parametrize("i", [0, 7, -1, -66])
def test_scalar_getitem_of_mask_does_not_gather(data, i):
    mask = np.arange(len(data)) % 3 != 1
    pending = data[mask]
    expected = data[mask].copy()[i]
    result = pending[i]
    if pd.isna(expected):
        assert result is data.dtype.na_value
    else:
        assert result == expected
    assert _is_pending(pending)


def test_mask_of_take_composes_positions(data):
    pending = data.take(np.arange(99, -1, -1))
    mask = np.arange(len(pending)) % 4 == 0
    result = pending[mask]
    assert result._pending_take.selection.dtype == np.intp
    tm.assert_extension_array_equal(result, data[::-1][mask])


def test_copy_of_mask_gathers(data):
    mask = np.arange(len(data)) % 2 == 1
    pending = data[mask]
    result = pending.copy()
    assert not _is_pending(result)
    assert _is_pending(pending)
    tm.assert_extension_array_equal(result, data[1::2])


def test_small_result_is_deferred(data):
    result = data.take(np.array([3, 1, 4]))
    assert _is_pending(result)
    tm.assert_extension_array_equal(result, _expected(data, [3, 1, 4]))

    pending = data.take(np.arange(99, -1, -1))
    result = pending.take(np.array([0, 1]))
    assert _is_pending(result)
    tm.assert_extension_array_equal(result, _expected(data, [99, 98]))
    assert _is_pending(pending)

    mask = np.zeros(len(data), dtype=bool)
    mask[[3, 50]] = True
    result = data[mask]
    assert _is_pending(result)
    tm.assert_extension_array_equal(result, _expected(data, [3, 50]))


def test_empty_result_is_eager(data):
    result = data[np.zeros(len(data), dtype=bool)]
    assert not _is_pending(result)
    assert len(result) == 0
    assert result.dtype == data.dtype

    pending = data.take(np.arange(99, -1, -1))
    result = pending[np.zeros(len(pending), dtype=bool)]
    assert not _is_pending(result)
    assert len(result) == 0


def test_slice_of_materialized_array_is_not_deferred(data):
    result = data[10:90]
    assert not _is_pending(result)


def test_full_slice_shares_pending(data):
    pending = data.take(np.arange(99, -1, -1))
    result = pending[:]
    assert result._pending_take is pending._pending_take


def _deferred(data, how):
    # a deferred selection of all 100 rows of ``data``, in reverse order
    if how == "take":
        return data.take(np.arange(99, -1, -1)), data[::-1]
    mask = np.ones(len(data), dtype=bool)
    mask[::7] = False
    # copy gathers, so the expected values do not go through the deferral
    return data[mask], data[mask].copy()


@pytest.mark.parametrize("how", ["take", "mask"])
def test_small_slice_gathers_only_its_rows(data, how):
    pending, expected = _deferred(data, how)
    head = pending[:10]
    tm.assert_extension_array_equal(head, expected[:10])
    # the slice gathered its own rows and left the selection deferred
    assert _is_pending(pending)
    assert pending._pending_take.result is None


@pytest.mark.parametrize("how", ["take", "mask"])
def test_overlapping_slices_share_one_gather(data, how):
    # once the slices would gather more rows than the selection holds, the
    # selection is gathered and the remaining slices are views of it
    pending, expected = _deferred(data, how)
    slices = [pending[:-i] for i in range(1, 6)]
    for i, sliced in enumerate(slices, start=1):
        tm.assert_extension_array_equal(sliced, expected[:-i])
    assert pending._pending_take.result is not None
    assert not tm.shares_memory(slices[0], pending)
    assert all(tm.shares_memory(sliced, pending) for sliced in slices[1:])
    tm.assert_extension_array_equal(pending, expected)


def test_slice_after_gather_is_a_view(data):
    pending = data.take(np.arange(99, -1, -1))
    sliced = pending[5:50]
    assert _is_pending(sliced)
    pending._pa_array
    tm.assert_extension_array_equal(sliced, data[::-1][5:50])
    assert tm.shares_memory(sliced, pending)


def test_slices_of_a_slice(data):
    pending = data.take(np.arange(99, -1, -1))
    outer = pending[:90]
    first = outer[:80]
    second = outer[5:85]
    # slices of ``outer`` are cut from ``pending``, not from ``outer``
    assert first._pending_take.origin is pending._pending_take
    assert second._pending_take.origin_start == 5
    tm.assert_extension_array_equal(first, data[::-1][:80])
    tm.assert_extension_array_equal(second, data[::-1][5:85])
    # 80 + 80 rows exceed the 100 of ``pending``, which then gathered whole;
    # ``second`` and ``outer`` are views of it
    assert pending._pending_take.result is not None
    assert tm.shares_memory(second, pending)
    assert not tm.shares_memory(first, pending)
    tm.assert_extension_array_equal(outer, data[::-1][:90])
    assert tm.shares_memory(outer, pending)


def test_long_chain_of_slices():
    # a sliding window: each step slices the previous slice, and reads a
    # window that overlaps the next one
    data = _arrow_array(pa.string(), n=5000)
    rest = data.take(np.arange(len(data))[::-1])
    root = rest._pending_take
    for _ in range(2000):
        window = rest[:3]
        window._pa_array
        rest = rest[1:]
        # the links stay one deep, so reading ``rest`` does not recurse
        assert rest._pending_take.origin in (root, None)
    tm.assert_extension_array_equal(rest, data[::-1][2000:])
    # the windows count against the rows of ``root``, so once they would
    # gather more than it holds, ``root`` is gathered whole
    assert root.result is not None
    assert tm.shares_memory(rest, window)


def test_step_slice_gathers_on_its_own(data):
    pending = data.take(np.arange(99, -1, -1))
    sliced = pending[::2]
    assert sliced._pending_take.origin is None
    tm.assert_extension_array_equal(sliced, data[::-1][::2])


def test_take_with_fill_is_eager(data):
    positions = np.array([0, -1, 5] * 20)
    result = data.take(positions, allow_fill=True)
    assert not _is_pending(result)
    expected = pd.array(
        [data[0], None, data[5]] * 20,
        dtype=data.dtype,
    )
    tm.assert_extension_array_equal(result, expected)


def test_take_of_pending_with_fill(data):
    pending = data.take(np.arange(99, -1, -1))
    result = pending.take(np.array([0, -1] * 30), allow_fill=True)
    expected = pd.array([data[99], None] * 30, dtype=data.dtype)
    tm.assert_extension_array_equal(result, expected)


def test_take_negative_without_fill(data):
    pending = data.take(np.arange(99, -1, -1))
    result = pending.take(np.array([-1] * 20))
    tm.assert_extension_array_equal(result, _expected(data, [0] * 20))


def test_take_out_of_bounds_raises(data):
    pending = data.take(np.arange(50))
    with pytest.raises(IndexError, match="out of bounds"):
        pending.take(np.array([50] * 10))


def test_fixed_width_not_deferred():
    # fixed-width values are cheap to gather
    arr = pd.array(np.arange(100), dtype="int64[pyarrow]")
    positions = np.arange(99, -1, -1)
    result = arr.take(positions)
    assert not _is_pending(result)
    tm.assert_extension_array_equal(result, _expected(arr, positions))
    mask = np.arange(100) % 2 == 0
    result = arr[mask]
    assert not _is_pending(result)
    tm.assert_extension_array_equal(result, _expected(arr, np.flatnonzero(mask)))


@pytest.mark.parametrize("sort", [True, False], ids=["sorted", "unsorted"])
def test_gather_across_chunks(data, sort):
    rng = np.random.default_rng(2)
    positions = rng.choice(len(data), size=60, replace=False)
    if sort:
        positions = np.sort(positions)
    result = data.take(positions)
    assert _is_pending(result)
    gathered = result._pa_array
    tm.assert_extension_array_equal(result, _expected(data, positions))
    if sort:
        # sorted positions are gathered chunk by chunk, like filter
        assert gathered.num_chunks == data._pa_array.num_chunks


@pytest.mark.slow
@pytest.mark.parametrize("pa_type", [pa.string(), pa.binary()], ids=str)
def test_gather_over_2gib(pa_type):
    # ChunkedArray.take concatenates the chunks, which overflows 32-bit
    # offsets past 2 GiB. The data buffer is never written, so this needs
    # little real memory.
    width, n = 1000, 1_150_000
    data = pa.py_buffer(np.zeros(width * n, dtype=np.uint8))
    offsets = pa.py_buffer((np.arange(n + 1, dtype=np.int64) * width).astype(np.int32))
    chunk = pa.Array.from_buffers(pa_type, n, [None, offsets, data])
    arr = pd.arrays.ArrowExtensionArray(pa.chunked_array([chunk, chunk]))

    # every other row: deferred, and never gathered in full
    pending = arr.take(np.arange(0, 2 * n, 2))
    assert _is_pending(pending)
    half = n // 2
    expected = pa.array([b"\0" * width] * 3, type=pa_type)

    # small results composed onto the parent, gathered from it: sorted
    # positions ...
    result = pending.take(np.array([0, half - 1, half]))
    assert pa.chunked_array(result._pa_array).combine_chunks().equals(expected)
    mask = np.zeros(len(pending), dtype=bool)
    mask[[1, half + 1, half + 2]] = True
    result = pending[mask]
    assert result._pa_array.combine_chunks().equals(expected)
    # ... and unsorted ones
    result = pending.take(np.array([half + 3, 2, half]))
    assert result._pa_array.combine_chunks().equals(expected)
    # a mask of the parent itself is gathered with filter, per chunk
    mask = np.zeros(2 * n, dtype=bool)
    mask[[0, n - 1, n]] = True
    result = arr[mask]
    assert result._pending_take.selection.dtype == np.bool_
    assert result._pa_array.combine_chunks().equals(expected)


@pytest.mark.parametrize("i", [0, 7, 42, -1, -100])
def test_scalar_getitem_does_not_gather(data, i):
    pending = data.take(np.arange(99, -1, -1))
    result = pending[i]
    expected = data[::-1][i]
    if pd.isna(expected):
        assert result is data.dtype.na_value
    else:
        assert result == expected
    assert _is_pending(pending)


def test_scalar_getitem_out_of_bounds(data):
    pending = data.take(np.arange(99, -1, -1))
    with pytest.raises(IndexError, match="out of bounds"):
        pending[100]


def test_materialize_releases_parent(data):
    result = data.take(np.arange(99, -1, -1))
    pending = result._pending_take
    result._pa_array
    assert pending.parent is None
    assert pending.selection is None


def test_view_shares_buffers(data):
    pending = data.take(np.arange(99, -1, -1))
    view = pending.view()
    assert _is_pending(view)
    assert tm.shares_memory(pending, view)
    assert view._pa_array is pending._pa_array


def test_copy_gathers(data):
    pending = data.take(np.arange(99, -1, -1))
    copied = pending.copy()
    # the copy owns its buffers, so it does not keep the parent alive
    assert not _is_pending(copied)
    assert _is_pending(pending)
    assert not tm.shares_memory(pending, copied)
    tm.assert_extension_array_equal(copied, pending)


@pytest.mark.parametrize("how", ["take", "mask"])
@pytest.mark.parametrize("select", ["copy", "take", "slice"])
def test_select_after_a_view_gathered(data, how, select):
    pending, expected = _deferred(data, how)
    view = pending[:]
    view._pa_array
    # ``pending`` still holds the take its view gathered
    assert pending._pending_take.result is not None
    if select == "copy":
        result = pending.copy()
        assert not tm.shares_memory(result, view)
        tm.assert_extension_array_equal(result, expected)
    elif select == "take":
        # deferred again, from the gathered values
        result = pending.take([3, 1, 2])
        assert _is_pending(result)
        tm.assert_extension_array_equal(result, expected.take([3, 1, 2]))
    else:
        # a view of the gathered values, as on an eager array
        result = pending[10:20]
        assert not _is_pending(result)
        assert tm.shares_memory(result, view)
        tm.assert_extension_array_equal(result, expected[10:20])


def test_getattr_after_another_thread_gathered(data):
    # the state a reader finds if another thread finishes the gather between
    # the failed lookup of ``_pa_array`` and the call to ``__getattr__``
    pending = data.take(np.arange(10))
    gathered = pending._pending_take.materialize()
    pending.__dict__["_pa_array"] = gathered
    del pending.__dict__["_pending_take"]
    assert type(pending).__getattr__(pending, "_pa_array") is gathered

    # neither the values nor a deferred take: a plain AttributeError
    empty = type(data).__new__(type(data))
    with pytest.raises(AttributeError, match="no attribute '_pa_array'"):
        type(data).__getattr__(empty, "_pa_array")
    with pytest.raises(AttributeError, match="no attribute 'missing'"):
        data.missing


def test_gather_reraises_with_64bit_offsets():
    # only 32-bit offsets can overflow, so with 64-bit offsets an error of
    # the take is raised as it is, instead of retrying chunk by chunk
    from pandas.core.arrays.arrow.array import _gather

    class Parent:
        num_chunks = 2
        type = pa.large_string()

        def take(self, indices):
            raise pa.ArrowInvalid("not an overflow")

    with pytest.raises(pa.ArrowInvalid, match="not an overflow"):
        _gather(Parent(), np.array([1, 0], dtype=np.intp))


def test_setitem_on_pending(data):
    pending = data.take(np.arange(99, -1, -1))
    view = pending.view()
    expected = data[::-1].copy()
    pending[0] = data[1]
    expected[0] = data[1]
    tm.assert_extension_array_equal(pending, expected)
    tm.assert_extension_array_equal(view, data[::-1])


def test_pickle_pending(data):
    pending = data.take(np.arange(99, -1, -1))
    result = pickle.loads(pickle.dumps(pending))
    tm.assert_extension_array_equal(result, data[::-1])
    assert not _is_pending(result)


@pytest.mark.parametrize("dump", [pickle.dumps, copy.deepcopy])
def test_pickle_during_first_read(data, dump):
    # the state another thread's pickle finds if it runs while the first
    # read of ``_pa_array`` has stored the values but not yet dropped the
    # deferred take
    pending = data.take(np.arange(99, -1, -1))
    pending.__dict__["_pa_array"] = pending._pending_take.materialize()
    assert "_pending_take" in pending.__dict__
    assert "_pending_take" not in pending.__getstate__()
    result = dump(pending)
    if dump is pickle.dumps:
        result = pickle.loads(result)
    tm.assert_extension_array_equal(result, data[::-1])


@pytest.mark.parametrize("share", [True, False], ids=["same-array", "views"])
def test_concurrent_materialize(share):
    # large enough that the gather, which releases the GIL, takes a while
    n = 200_000
    data = pd.array([f"value{i}" for i in range(n)], dtype="str")
    pending = data.take(np.arange(n - 1, -1, -1))
    views = [pending if share else pending.view() for _ in range(8)]
    barrier = threading.Barrier(len(views))
    results: list = [None] * len(views)
    errors: list[BaseException] = []

    def run(i):
        barrier.wait()
        try:
            assert len(views[i]) == n
            results[i] = views[i]._pa_array
            assert len(views[i]) == n
        except BaseException as err:
            errors.append(err)

    threads = [threading.Thread(target=run, args=(i,)) for i in range(len(views))]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert not errors
    assert all(result is results[0] for result in results)
    tm.assert_extension_array_equal(pending, data[::-1])


def test_concurrent_overlapping_slices():
    n = 200_000
    data = pd.array([f"value{i}" for i in range(n)], dtype="str")
    pending = data.take(np.arange(n - 1, -1, -1))
    slices = [pending[i : n - i] for i in range(8)]
    barrier = threading.Barrier(len(slices))
    errors: list[BaseException] = []

    def run(i):
        barrier.wait()
        try:
            slices[i]._pa_array
        except BaseException as err:
            errors.append(err)

    threads = [threading.Thread(target=run, args=(i,)) for i in range(len(slices))]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert not errors
    expected = data[::-1]
    for i, sliced in enumerate(slices):
        tm.assert_extension_array_equal(sliced, expected[i : n - i])
    # at most one slice gathered on its own before the rest became views
    assert pending._pending_take.result is not None
    assert sum(not tm.shares_memory(sliced, pending) for sliced in slices) <= 1


def test_series_mask_and_merge():
    n = 1_000
    left = pd.DataFrame(
        {
            "key": np.arange(n) % 50,
            "name": [f"name{i}" for i in range(n)],
            "comment": [f"comment{i}" for i in range(n)],
        }
    )
    right = pd.DataFrame({"key": np.arange(50), "label": [f"l{i}" for i in range(50)]})
    filtered = left[left["key"] % 2 == 0]
    result = filtered.merge(right, on="key")
    assert len(result) == n // 2
    assert (result["label"] == "l" + result["key"].astype(str)).all()
    assert (result["name"].str.slice(4).astype(int) % 50 == result["key"]).all()
