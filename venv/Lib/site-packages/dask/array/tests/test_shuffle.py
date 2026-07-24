from __future__ import annotations

import numpy as np
import pytest

import dask
import dask.array as da
from dask.array import assert_eq, shuffle
from dask.array._shuffle import _rechunk_other_dimensions
from dask.core import flatten


@pytest.fixture()
def arr():
    return np.arange(0, 24).reshape(8, 3).T.copy()


@pytest.fixture()
def darr(arr):
    return da.from_array(arr, chunks=((2, 1), (4, 4)))


@pytest.mark.parametrize(
    "indexer, chunks, other_chunks",
    [
        ([[1, 5, 6], [0, 2, 3, 4, 7]], (3, 5), (2, 1)),
        ([[1, 5, 6], [0, 3], [4, 2, 7]], (5, 3), (2, 1)),
        ([[1], [0, 6, 5, 3, 2, 4], [7]], (1, 6, 1), (1, 1, 1)),
        ([[1, 5, 1, 5, 1, 5], [1, 6, 4, 2, 7]], (6, 5), (1, 1, 1)),
    ],
)
def test_shuffle(arr, darr, indexer, chunks, other_chunks):
    result = darr.shuffle(indexer, axis=1)
    expected = arr[:, list(flatten(indexer))]
    assert_eq(result, expected)
    assert result.chunks[0] == other_chunks
    assert result.chunks[1] == chunks


@pytest.mark.parametrize("tol, chunks", ((1, (3, 2, 3)), (1.4, (5, 3))))
def test_shuffle_config_tolerance(arr, darr, tol, chunks):
    indexer = [[1, 5, 6], [0, 3], [4, 2, 7]]
    with dask.config.set({"array.chunk-size-tolerance": tol}):
        result = darr.shuffle(indexer, axis=1)
    expected = arr[:, [1, 5, 6, 0, 3, 4, 2, 7]]
    assert_eq(result, expected)
    assert result.chunks[0] == darr.chunks[0]
    assert result.chunks[1] == chunks


def test_shuffle_larger_array():
    arr = da.random.random((15, 15, 15), chunks=(5, 5, 5))
    indexer = np.arange(0, 15)
    np.random.shuffle(indexer)
    indexer = [indexer[0:6], indexer[6:8], indexer[8:9], indexer[9:]]
    indexer = list(map(list, indexer))
    take_indexer = list(flatten(indexer))
    assert_eq(shuffle(arr, indexer, axis=1), arr[..., take_indexer, :])


def test_incompatible_indexer(darr):
    with pytest.raises(ValueError, match="indexer must be a list of lists"):
        darr.shuffle("s", axis=1)

    with pytest.raises(ValueError, match="indexer must be a list of lists"):
        darr.shuffle([1], axis=1)


def test_unknown_chunk_sizes(darr):
    darr._chunks = ((np.nan, 1), (4, 4))
    with pytest.raises(
        ValueError, match="Shuffling only allowed with known chunk sizes"
    ):
        darr.shuffle([[1]], axis=1)


def test_oob_axis(darr):
    with pytest.raises(ValueError, match="is out of bounds"):
        darr.shuffle([[1]], axis=5)


def test_oob_indexer(darr):
    with pytest.raises(IndexError, match="Indexer contains out of bounds index"):
        darr.shuffle([[16]], axis=1)


def test_shuffle_no_op_with_correct_indexer():
    arr = da.ones((250, 100), chunks=((50, 100, 33, 67), 100))
    indexer = [
        list(range(0, 50)),
        list(range(50, 150)),
        list(range(150, 183)),
        list(range(183, 250)),
    ]
    result = arr.shuffle(indexer, axis=0)
    assert result.dask == arr.dask
    assert_eq(arr, result)


def test_resize_other_dimensions():
    arr = da.random.random((250, 50), chunks=((45, 100, 38, 67), 10))
    result = _rechunk_other_dimensions(arr, 20, 1, "auto")
    assert result.chunks == ((45, 50, 50, 38, 34, 33), (10,) * 5)
    assert_eq(arr, result)

    arr = da.random.random((250, 50, 20), chunks=((45, 100, 38, 67), 10, 10))
    result = _rechunk_other_dimensions(arr, 20, 1, "auto")
    assert result.chunks == ((45, 50, 50, 38, 67), (10,) * 5, (5,) * 4)
    assert_eq(arr, result)

    result = _rechunk_other_dimensions(arr, 40, 1, "auto")
    assert result.chunks == ((45, 50, 50, 38, 34, 33), (10,) * 5, (5,) * 4)
    assert_eq(arr, result)

    arr = da.random.random((250, 50, 5), chunks=((45, 100, 38, 67), 10, 1))
    result = _rechunk_other_dimensions(arr, 40, 1, "auto")
    assert result.chunks == (
        (23, 22, 25, 25, 25, 25, 19, 19, 23, 22, 22),
        (10,) * 5,
        (1,) * 5,
    )
    assert_eq(arr, result)

    arr = da.random.random((5, 50, 5), chunks=(5, 10, (2, 3)))
    result = _rechunk_other_dimensions(arr, 100, 1, "auto")
    assert result.chunks == ((2, 1, 1, 1), (10,) * 5, (1,) * 5)
    assert_eq(arr, result)

    arr = da.random.random((5, 50, 5), chunks=(5, 10, (2, 3)))
    result = _rechunk_other_dimensions(arr, 1000, 1, "auto")
    assert result.chunks == ((1,) * 5, (10,) * 5, (1,) * 5)
    assert_eq(arr, result)

    arr = da.random.random((5, 50, 5), chunks=(5, 10, (2, 3)))
    result = _rechunk_other_dimensions(arr, 250, 1, "auto")
    assert result.chunks == ((1,) * 5, (10,) * 5, (1,) * 5)
    assert_eq(arr, result)

    arr = da.random.random((2, 1, 2), chunks=(2, 1, 2))
    result = _rechunk_other_dimensions(arr, 4, 1, "auto")
    assert result.chunks == ((1, 1), (1,), (1, 1))
    assert_eq(arr, result)


def test_dtype_taker(arr, darr):
    indexer = [[1, 5, 6], [0, 3], [4, 2, 7]]
    result = darr.shuffle(indexer, axis=1)
    expected = arr[:, [1, 5, 6, 0, 3, 4, 2, 7]]
    assert_eq(result, expected)
    assert all(
        v.value[1][1].dtype == np.uint8
        for k, v in dict(result.dask).items()
        if "shuffle-taker" in k
    )
