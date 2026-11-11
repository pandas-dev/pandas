"""
Tests for the deferred gather of string and binary Arrow arrays: ``take``,
slices and boolean masks of such an ``ArrowExtensionArray`` record the
positions or the mask and gather on first access to ``_pa_array``.
"""

from __future__ import annotations

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


def test_take_keeps_owned_readonly_indexer(data):
    positions = np.arange(50)
    positions.flags.writeable = False
    result = data.take(positions)
    assert result._pending_take.selection is positions

    # a read-only view could still change through its base: copied
    base = np.arange(50)
    view = base[:]
    view.flags.writeable = False
    result = data.take(view)
    assert result._pending_take.selection is not view
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
