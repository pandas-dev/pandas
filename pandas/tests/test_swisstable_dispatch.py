import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame
from pandas._libs import hashtable as htable
from pandas.core import algorithms


def _unexpected_khash_call(*args, **kwargs):
    raise AssertionError("khash fallback was called")


def test_isin_numeric_uses_swisstable(monkeypatch):
    monkeypatch.setattr(htable, "ismember", _unexpected_khash_call)

    with pd.option_context("compute.use_swisstable", True):
        result = algorithms.isin(
            np.array([1, 2, 3, 2], dtype=np.int64),
            np.array([2, 4], dtype=np.int64),
        )

    expected = np.array([False, True, False, True])
    np.testing.assert_array_equal(result, expected)


def test_duplicated_numeric_uses_swisstable(monkeypatch):
    monkeypatch.setattr(htable, "duplicated", _unexpected_khash_call)

    with pd.option_context("compute.use_swisstable", True):
        result = algorithms.duplicated(np.array([1, 2, 1, 3], dtype=np.int64))

    expected = np.array([False, False, True, False])
    np.testing.assert_array_equal(result, expected)


def test_value_counts_numeric_uses_swisstable(monkeypatch):
    monkeypatch.setattr(htable, "value_count", _unexpected_khash_call)

    with pd.option_context("compute.use_swisstable", True):
        keys, counts, na_count = algorithms.value_counts_arraylike(
            np.array([2, 1, 2, 3], dtype=np.int64),
            dropna=True,
        )

    np.testing.assert_array_equal(keys, np.array([2, 1, 3], dtype=np.int64))
    np.testing.assert_array_equal(counts, np.array([2, 1, 1], dtype=np.int64))
    assert na_count == 0


def test_series_numeric_methods_use_swisstable(monkeypatch):
    monkeypatch.setattr(htable, "ismember", _unexpected_khash_call)
    monkeypatch.setattr(htable, "duplicated", _unexpected_khash_call)
    monkeypatch.setattr(htable, "value_count", _unexpected_khash_call)
    series = pd.Series([3, 1, 3, 2], dtype=np.int64)

    with pd.option_context("compute.use_swisstable", True):
        isin_result = series.isin([2, 3])
        duplicated_result = series.duplicated()
        counts_result = series.value_counts(sort=False)

    expected_isin = pd.Series([True, False, True, True])
    expected_duplicated = pd.Series([False, False, True, False])
    expected_counts = pd.Series(
        [2, 1, 1],
        index=pd.Index([3, 1, 2], dtype=np.int64),
        name="count",
    )
    pd.testing.assert_series_equal(isin_result, expected_isin)
    pd.testing.assert_series_equal(duplicated_result, expected_duplicated)
    pd.testing.assert_series_equal(counts_result, expected_counts)


@pytest.mark.parametrize(
    "dtype",
    [
        np.int8,
        np.int16,
        np.int32,
        np.int64,
        np.uint8,
        np.uint16,
        np.uint32,
        np.uint64,
        np.float32,
        np.float64,
        np.complex64,
        np.complex128,
    ],
)
def test_numeric_algorithms_match_khash(dtype):
    values = np.array([3, 1, 2, 1, 3], dtype=dtype)
    targets = np.array([1, 4], dtype=dtype)

    with pd.option_context("compute.use_swisstable", False):
        expected_isin = algorithms.isin(values, targets)
        expected_duplicated = algorithms.duplicated(values, keep=False)
        expected_counts = algorithms.value_counts_arraylike(values, dropna=True)

    with pd.option_context("compute.use_swisstable", True):
        result_isin = algorithms.isin(values, targets)
        result_duplicated = algorithms.duplicated(values, keep=False)
        result_counts = algorithms.value_counts_arraylike(values, dropna=True)

    np.testing.assert_array_equal(result_isin, expected_isin)
    np.testing.assert_array_equal(result_duplicated, expected_duplicated)
    for result, expected in zip(result_counts, expected_counts, strict=True):
        np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.complex64, np.complex128])
def test_numeric_algorithms_nan_semantics_match_khash(dtype):
    values = np.array([1, np.nan, 2, np.nan], dtype=dtype)
    targets = np.array([np.nan, 3], dtype=dtype)

    with pd.option_context("compute.use_swisstable", False):
        expected_isin = algorithms.isin(values, targets)
        expected_duplicated = algorithms.duplicated(values)
        expected_counts = algorithms.value_counts_arraylike(values, dropna=False)

    with pd.option_context("compute.use_swisstable", True):
        result_isin = algorithms.isin(values, targets)
        result_duplicated = algorithms.duplicated(values)
        result_counts = algorithms.value_counts_arraylike(values, dropna=False)

    np.testing.assert_array_equal(result_isin, expected_isin)
    np.testing.assert_array_equal(result_duplicated, expected_duplicated)
    np.testing.assert_array_equal(result_counts[0], expected_counts[0])
    np.testing.assert_array_equal(result_counts[1], expected_counts[1])
    assert result_counts[2] == expected_counts[2]


@pytest.mark.parametrize("dropna", [False, True])
def test_masked_numeric_algorithms_match_khash(dropna):
    values = np.array([1, 9, 1, 9], dtype=np.int64)
    mask = np.array([False, True, False, True])

    with pd.option_context("compute.use_swisstable", False):
        expected_duplicated = algorithms.duplicated(values, mask=mask)
        expected_counts = algorithms.value_counts_arraylike(
            values, dropna=dropna, mask=mask
        )

    with pd.option_context("compute.use_swisstable", True):
        result_duplicated = algorithms.duplicated(values, mask=mask)
        result_counts = algorithms.value_counts_arraylike(
            values, dropna=dropna, mask=mask
        )

    np.testing.assert_array_equal(result_duplicated, expected_duplicated)
    np.testing.assert_array_equal(result_counts[0], expected_counts[0])
    np.testing.assert_array_equal(result_counts[1], expected_counts[1])
    assert result_counts[2] == expected_counts[2]


def test_numeric_inner_merge_uses_swisstable_hash_join():
    left = DataFrame({"key": [3, 1, 2], "left": ["c", "a", "b"]})
    right = DataFrame({"key": [2, 3], "right": ["B", "C"]})

    with pd.option_context("compute.use_swisstable", True):
        result = left.merge(right, on="key", how="inner", sort=False)

    expected = DataFrame(
        {"key": [3, 2], "left": ["c", "b"], "right": ["C", "B"]}
    )
    pd.testing.assert_frame_equal(result, expected)


def test_nullable_numeric_inner_merge_uses_swisstable_hash_join():
    left = DataFrame(
        {
            "key": pd.array([3, None, 2], dtype="Int64"),
            "left": ["c", "missing", "b"],
        }
    )
    right = DataFrame(
        {"key": pd.array([2, 3], dtype="Int64"), "right": ["B", "C"]}
    )

    with pd.option_context("compute.use_swisstable", True):
        result = left.merge(right, on="key", how="inner", sort=False)

    expected = DataFrame(
        {
            "key": pd.array([3, 2], dtype="Int64"),
            "left": ["c", "b"],
            "right": ["C", "B"],
        }
    )
    pd.testing.assert_frame_equal(result, expected)


def test_numeric_merge_sort_uses_swisstable_uniques():
    left = DataFrame({"key": [3, 1, 2], "left": ["c", "a", "b"]})
    right = DataFrame({"key": [2, 3], "right": ["B", "C"]})

    with pd.option_context("compute.use_swisstable", True):
        result = left.merge(right, on="key", how="left", sort=True)

    expected = DataFrame(
        {
            "key": [1, 2, 3],
            "left": ["a", "b", "c"],
            "right": [np.nan, "B", "C"],
        }
    )
    pd.testing.assert_frame_equal(result, expected)
