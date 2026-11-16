"""
Tests for running the sorted-join kernels on a thread pool (GH#51364).

Each test forces chunking on small inputs and compares against the serial
kernel, which is the oracle.
"""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

from pandas._libs import join as libjoin
from pandas.compat import WASM

import pandas as pd
import pandas._testing as tm
from pandas.core.indexes import parallel_join

pytestmark = pytest.mark.skipif(WASM, reason="cannot start threads")


@pytest.fixture
def chunked(monkeypatch):
    """Chunk even tiny inputs, and count how often the threaded path runs."""
    monkeypatch.setattr(parallel_join, "MIN_CHUNK_SIZE", 2)
    calls = []
    run_chunks = parallel_join._run_chunks

    def counting_run_chunks(run, n_chunks):
        calls.append(n_chunks)
        return run_chunks(run, n_chunks)

    monkeypatch.setattr(parallel_join, "_run_chunks", counting_run_chunks)
    with pd.option_context("mode.max_threads", 4):
        yield calls


# long runs of equal values that straddle the chunk boundaries, values present
# on only one side (-1 in the indexers), and empty or disjoint sides
CASES = {
    "runs": ([0, 0, 0, 1, 1, 2, 2, 2, 2, 4, 5, 5], [0, 1, 1, 2, 3, 3, 5, 5, 5, 6]),
    "left-only": ([1, 2, 3, 4, 5, 6, 7, 8], []),
    "right-only": ([], [1, 2, 3, 4, 5, 6, 7, 8]),
    "disjoint": ([0, 1, 2, 3, 4, 5], [6, 7, 8, 9, 10, 11]),
    "skewed": ([4], list(range(20))),
}


@pytest.mark.parametrize(
    "kernel",
    [
        libjoin.left_join_indexer,
        libjoin.inner_join_indexer,
        libjoin.outer_join_indexer,
    ],
)
@pytest.mark.parametrize("case", CASES)
@pytest.mark.parametrize("dtype", ["int64", "uint8", "float64"])
def test_join_indexer_matches_serial(chunked, kernel, case, dtype):
    left, right = (np.array(values, dtype=dtype) for values in CASES[case])
    expected = kernel(left, right)
    result = parallel_join.join_indexer(kernel, left, right)
    assert chunked
    for res, exp in zip(result, expected, strict=True):
        tm.assert_numpy_array_equal(res, exp)


@pytest.mark.parametrize("case", CASES)
def test_left_join_indexer_unique_matches_serial(chunked, case):
    left, right = (np.array(values, dtype="int64") for values in CASES[case])
    right = np.unique(right)
    kernel = libjoin.left_join_indexer_unique
    expected = kernel(left, right)
    result = parallel_join.left_join_indexer_unique(kernel, left, right)
    assert chunked
    tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize(
    "values",
    [
        np.array([-np.inf, -1.5, -0.0, 0.0, 0.0, 2.5, 2.5, np.inf]),
        pd.date_range("2016-01-01", periods=8, freq="D", tz="US/Eastern"),
        pd.array(range(8), dtype="Int64"),
        pd.Categorical(range(8)),
    ],
    ids=["float-signed-zero", "datetime-tz", "Int64", "categorical"],
)
@pytest.mark.parametrize("how", ["left", "right", "inner", "outer"])
def test_index_join_matches_serial(chunked, values, how):
    left = pd.Index(values[[0, 1, 2, 2, 3, 4, 4, 4, 6]])
    # unique, so the libjoin path applies; -0.0 on the left matches 0.0 here
    right = pd.Index(values[[1, 3, 5, 7]])
    with pd.option_context("mode.max_threads", 1):
        expected = left.join(right, how=how, return_indexers=True)
    calls_before = len(chunked)
    result = left.join(right, how=how, return_indexers=True)
    assert len(chunked) > calls_before
    tm.assert_index_equal(result[0], expected[0])
    for res, exp in zip(result[1:], expected[1:], strict=True):
        if exp is None:
            assert res is None
        else:
            tm.assert_numpy_array_equal(res, exp)


def test_merge_and_align_match_serial(chunked):
    keys = np.array([0, 0, 1, 2, 2, 2, 3, 5, 5, 8])
    left = pd.DataFrame({"key": keys, "a": range(len(keys))})
    right = pd.DataFrame({"key": [0, 2, 3, 4, 5, 9], "b": range(6)})
    ser1 = pd.Series(range(6), index=[0, 1, 2, 4, 6, 7])
    ser2 = pd.Series(range(6), index=[1, 2, 3, 5, 7, 9])
    with pd.option_context("mode.max_threads", 1):
        expected_merge = pd.merge(left, right, on="key", how="outer")
        expected_sum = ser1 + ser2
    tm.assert_frame_equal(pd.merge(left, right, on="key", how="outer"), expected_merge)
    assert len(chunked) == 1
    tm.assert_series_equal(ser1 + ser2, expected_sum)
    assert len(chunked) == 2


def test_object_dtype_stays_serial(chunked):
    # object kernels hold the GIL, so threads would gain nothing
    left = np.array(["a", "b", "c", "d"], dtype=object)
    right = np.array(["b", "c", "e", "f"], dtype=object)
    expected = libjoin.inner_join_indexer(left, right)
    result = parallel_join.join_indexer(libjoin.inner_join_indexer, left, right)
    assert not chunked
    for res, exp in zip(result, expected, strict=True):
        tm.assert_numpy_array_equal(res, exp)


def test_single_run_of_equal_values_stays_serial(chunked):
    # every cut lands at the start of the run, leaving one non-empty chunk
    left = np.array([3] * 10)
    right = np.array([3, 3])
    expected = libjoin.inner_join_indexer(left, right)
    result = parallel_join.join_indexer(libjoin.inner_join_indexer, left, right)
    assert not chunked
    for res, exp in zip(result, expected, strict=True):
        tm.assert_numpy_array_equal(res, exp)


def test_max_threads_caps_chunk_count(chunked):
    left = np.arange(40)
    with pd.option_context("mode.max_threads", 3):
        parallel_join.join_indexer(libjoin.outer_join_indexer, left, left + 5)
    assert chunked == [3]


def test_max_threads_one_stays_serial(chunked):
    left = np.arange(10)
    with pd.option_context("mode.max_threads", 1):
        parallel_join.join_indexer(libjoin.outer_join_indexer, left, left + 5)
    assert not chunked


def test_thread_start_failure_falls_back_to_serial(chunked, monkeypatch):
    # e.g. a join inside an atexit handler, after the interpreter has stopped
    # accepting new threads
    def refuse(self, *args, **kwargs):
        raise RuntimeError("cannot schedule new futures after interpreter shutdown")

    monkeypatch.setattr(ThreadPoolExecutor, "submit", refuse)
    left, right = np.arange(10), np.arange(5, 15)
    expected = libjoin.outer_join_indexer(left, right)
    result = parallel_join.join_indexer(libjoin.outer_join_indexer, left, right)
    assert chunked
    for res, exp in zip(result, expected, strict=True):
        tm.assert_numpy_array_equal(res, exp)


def test_error_in_chunk_propagates(chunked):
    # a RuntimeError from the kernel must not trigger the serial fallback,
    # which here would succeed and hide it
    left = np.arange(10)

    def kernel(lvals, rvals):
        if len(lvals) < len(left):
            raise RuntimeError("kernel failed")
        return libjoin.inner_join_indexer(lvals, rvals)

    with pytest.raises(RuntimeError, match="kernel failed"):
        parallel_join.join_indexer(kernel, left, left)
