import numpy as np
import pytest

from pandas._libs import index as libindex

import pandas as pd
from pandas import (
    Index,
    NaT,
    RangeIndex,
)
import pandas._testing as tm


class TestGetSliceBounds:
    @pytest.mark.parametrize("side, expected", [("left", 4), ("right", 5)])
    def test_get_slice_bounds_within(self, side, expected):
        index = Index(list("abcdef"))
        result = index.get_slice_bound("e", side=side)
        assert result == expected

    @pytest.mark.parametrize("side", ["left", "right"])
    @pytest.mark.parametrize(
        "data, bound, expected", [(list("abcdef"), "x", 6), (list("bcdefg"), "a", 0)]
    )
    def test_get_slice_bounds_outside(self, side, expected, data, bound):
        index = Index(data)
        result = index.get_slice_bound(bound, side=side)
        assert result == expected

    def test_get_slice_bounds_invalid_side(self):
        with pytest.raises(ValueError, match="Invalid value for side kwarg"):
            Index([]).get_slice_bound("a", side="middle")


class TestGetIndexerNonUnique:
    def test_get_indexer_non_unique_dtype_mismatch(self):
        # GH#25459
        indexes, missing = Index(["A", "B"]).get_indexer_non_unique(Index([0]))
        tm.assert_numpy_array_equal(np.array([-1], dtype=np.intp), indexes)
        tm.assert_numpy_array_equal(np.array([0], dtype=np.intp), missing)

    @pytest.mark.parametrize(
        "idx_values,idx_non_unique",
        [
            ([np.nan, 100, 200, 100], [np.nan, 100]),
            ([np.nan, 100.0, 200.0, 100.0], [np.nan, 100.0]),
        ],
    )
    def test_get_indexer_non_unique_int_index(self, idx_values, idx_non_unique):
        indexes, missing = Index(idx_values).get_indexer_non_unique(Index([np.nan]))
        tm.assert_numpy_array_equal(np.array([0], dtype=np.intp), indexes)
        tm.assert_numpy_array_equal(np.array([], dtype=np.intp), missing)

        indexes, missing = Index(idx_values).get_indexer_non_unique(
            Index(idx_non_unique)
        )
        tm.assert_numpy_array_equal(np.array([0, 1, 3], dtype=np.intp), indexes)
        tm.assert_numpy_array_equal(np.array([], dtype=np.intp), missing)


class TestGetLoc:
    @pytest.mark.slow  # to_flat_index takes a while
    def test_get_loc_tuple_monotonic_above_size_cutoff(self, monkeypatch):
        # Go through the libindex path for which using
        # _bin_search vs ndarray.searchsorted makes a difference

        with monkeypatch.context():
            monkeypatch.setattr(libindex, "_SIZE_CUTOFF", 100)
            lev = list("ABCD")
            dti = pd.date_range("2016-01-01", periods=10)

            mi = pd.MultiIndex.from_product([lev, range(5), dti])
            oidx = mi.to_flat_index()

            loc = len(oidx) // 2
            tup = oidx[loc]

            res = oidx.get_loc(tup)
        assert res == loc

    def test_get_loc_tuple_monotonic_with_duplicates(self):
        # GH#37800 _get_loc_duplicates should not use ndarray.searchsorted
        # with tuple keys, as searchsorted interprets tuples as array-like
        values = np.empty(3, dtype=object)
        values[0] = (1, 1)
        values[1] = (1, 1)
        values[2] = (2, 2)
        idx = Index(values)
        result = idx.get_loc((1, 1))
        assert result == slice(0, 2)

        result = idx.get_loc((2, 2))
        assert result == 2

        with pytest.raises(KeyError, match=r"\(3, 3\)"):
            idx.get_loc((3, 3))

    def test_get_loc_nan_object_dtype_nonmonotonic_nonunique(self):
        # case that goes through _maybe_get_bool_indexer
        idx = Index(["foo", np.nan, None, "foo", 1.0, None], dtype=object)

        # we dont raise KeyError on nan
        res = idx.get_loc(np.nan)
        assert res == 1

        # we only match on None, not on np.nan
        res = idx.get_loc(None)
        expected = np.array([False, False, True, False, False, True])
        tm.assert_numpy_array_equal(res, expected)

        # we don't match at all on mismatched NA
        with pytest.raises(KeyError, match="NaT"):
            idx.get_loc(NaT)


def test_get_indexer_monotonic_above_size_cutoff(monkeypatch):
    # GH#14273 get_indexer should avoid building a hash table for large
    # monotonic indices.
    monkeypatch.setattr(libindex, "_SIZE_CUTOFF", 100)
    idx = Index(np.arange(200, dtype=np.int64))

    target = Index([10, 50, 150, 999])
    result = idx.get_indexer(target)
    expected = np.array([10, 50, 150, -1], dtype=np.intp)
    tm.assert_numpy_array_equal(result, expected)

    # Verify hash table was NOT built
    assert not idx._engine.is_mapping_populated


def test_get_indexer_monotonic_above_size_cutoff_datetime(monkeypatch):
    # GH#14273
    monkeypatch.setattr(libindex, "_SIZE_CUTOFF", 100)
    dti = pd.date_range("2000-01-01", periods=200)

    target = dti[[5, 100, 199]]
    result = dti.get_indexer(target)
    expected = np.array([5, 100, 199], dtype=np.intp)
    tm.assert_numpy_array_equal(result, expected)

    assert not dti._engine.is_mapping_populated


def test_get_indexer_monotonic_above_size_cutoff_not_found(monkeypatch):
    # GH#14273
    monkeypatch.setattr(libindex, "_SIZE_CUTOFF", 100)
    idx = Index(np.arange(0, 400, 2, dtype=np.int64))  # even numbers only

    target = Index([1, 3, 5])  # odd numbers, not in index
    result = idx.get_indexer(target)
    expected = np.array([-1, -1, -1], dtype=np.intp)
    tm.assert_numpy_array_equal(result, expected)

    assert not idx._engine.is_mapping_populated


def test_get_indexer_monotonic_above_size_cutoff_mixed(monkeypatch):
    # GH#14273
    monkeypatch.setattr(libindex, "_SIZE_CUTOFF", 100)
    idx = Index(np.arange(0, 400, 2, dtype=np.int64))

    # Mix of found and not found, including values beyond the range
    target = Index([0, 1, 398, 399, -1, 500])
    result = idx.get_indexer(target)
    expected = np.array([0, -1, 199, -1, -1, -1], dtype=np.intp)
    tm.assert_numpy_array_equal(result, expected)

    assert not idx._engine.is_mapping_populated


def test_getitem_boolean_ea_indexer():
    # GH#45806
    ser = pd.Series([True, False, pd.NA], dtype="boolean")
    result = ser.index[ser]
    expected = RangeIndex(start=0, stop=1, step=1)
    tm.assert_index_equal(result, expected, exact=True)
