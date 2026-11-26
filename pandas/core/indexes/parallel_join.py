"""
Run the sorted-join kernels in ``pandas._libs.join`` over value-range chunks
on a thread pool.  See GH#51364.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import (
    TYPE_CHECKING,
    TypeVar,
)

import numpy as np

from pandas.compat._cpu import get_max_threads

if TYPE_CHECKING:
    from collections.abc import Callable

    from pandas._typing import npt

    JoinResult = tuple[np.ndarray, npt.NDArray[np.intp], npt.NDArray[np.intp]]

T = TypeVar("T")

# Aim for at least this many elements (left + right) per chunk, so thread
# startup stays a small share of the kernel time even on a busy machine.
MIN_CHUNK_SIZE = 250_000


def _chunk_bounds(
    left: np.ndarray, right: np.ndarray
) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Return chunk start/stop positions in ``left`` and ``right``, or None to
    run serially.

    Both sides are cut at the leftmost occurrence of each boundary value, so a
    run of equal values never straddles two chunks.  The cuts are spaced evenly
    over the longer side only, so chunk sizes are uneven when the two sides'
    values are distributed differently.
    """
    max_chunks = (len(left) + len(right)) // MIN_CHUNK_SIZE
    if max_chunks < 2 or left.dtype == object:
        # object kernels hold the GIL
        return None
    n_chunks = min(get_max_threads(), max_chunks)
    if n_chunks < 2:
        return None
    larger = left if len(left) >= len(right) else right
    positions = np.arange(1, n_chunks) * len(larger) // n_chunks
    values = np.unique(larger[positions])
    lbounds = np.concatenate([[0], left.searchsorted(values, "left"), [len(left)]])
    rbounds = np.concatenate([[0], right.searchsorted(values, "left"), [len(right)]])
    # the first chunk is empty when a cut lands on the smallest value
    _, keep = np.unique(lbounds + rbounds, return_index=True)
    if len(keep) < 3:
        return None
    return lbounds[keep], rbounds[keep]


def _shift(indexer: np.ndarray, offset: int) -> np.ndarray:
    """Add ``offset`` in place, leaving the -1 "no match" sentinel as is."""
    if offset == 0:
        return indexer
    # sign is -1 at the sentinel and 0 elsewhere; much faster than ``where=``
    sign = indexer >> (indexer.dtype.itemsize * 8 - 1)
    indexer += offset
    indexer |= sign
    return indexer


def _run_chunks(run: Callable[[int], T], n_chunks: int) -> list[T] | None:
    """
    Return ``[run(k) for k in range(n_chunks)]`` computed on a thread pool, or
    None if no threads can be started, so the caller runs serially.
    """
    with ThreadPoolExecutor(n_chunks) as pool:
        try:
            futures = [pool.submit(run, k) for k in range(n_chunks)]
        except RuntimeError:
            # interpreter shutdown (e.g. inside an atexit handler) or the OS
            # refused a new thread; errors raised by ``run`` are not caught here
            return None
        return [future.result() for future in futures]


def join_indexer(
    kernel: Callable[[np.ndarray, np.ndarray], JoinResult],
    left: np.ndarray,
    right: np.ndarray,
) -> JoinResult:
    """
    Chunked equivalent of ``kernel(left, right)`` for ``left_join_indexer``,
    ``inner_join_indexer`` and ``outer_join_indexer``.
    """
    bounds = _chunk_bounds(left, right)
    if bounds is None:
        return kernel(left, right)
    lbounds, rbounds = bounds

    def run(k: int) -> JoinResult:
        lstart, rstart = lbounds[k], rbounds[k]
        joined, lidx, ridx = kernel(
            left[lstart : lbounds[k + 1]], right[rstart : rbounds[k + 1]]
        )
        return joined, _shift(lidx, lstart), _shift(ridx, rstart)

    pieces = _run_chunks(run, len(lbounds) - 1)
    if pieces is None:
        return kernel(left, right)
    columns = [list(column) for column in zip(*pieces, strict=True)]
    del pieces
    result = []
    for column in columns:
        result.append(np.concatenate(column))
        # free this output's pieces before building the next, to limit peak memory
        column.clear()
    return result[0], result[1], result[2]


def left_join_indexer_unique(
    kernel: Callable[[np.ndarray, np.ndarray], npt.NDArray[np.intp]],
    left: np.ndarray,
    right: np.ndarray,
) -> npt.NDArray[np.intp]:
    """
    Chunked equivalent of ``libjoin.left_join_indexer_unique(left, right)``.
    """
    bounds = _chunk_bounds(left, right)
    if bounds is None:
        return kernel(left, right)
    lbounds, rbounds = bounds

    def run(k: int) -> npt.NDArray[np.intp]:
        rstart = rbounds[k]
        ridx = kernel(left[lbounds[k] : lbounds[k + 1]], right[rstart : rbounds[k + 1]])
        return _shift(ridx, rstart)

    pieces = _run_chunks(run, len(lbounds) - 1)
    if pieces is None:
        return kernel(left, right)
    return np.concatenate(pieces)
