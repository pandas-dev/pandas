"""
Tests for 2D support in BaseMaskedArray methods.

These tests exercise code paths that were previously gated with
``ndim > 1`` checks or assumed 1D shapes. Each section corresponds to
a method that was made 2D-compatible.
"""

import numpy as np
import pytest

import pandas as pd
from pandas.core.arrays import (
    BooleanArray,
    FloatingArray,
    IntegerArray,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def int_arr2d():
    # shape (2, 3):
    #   [[1,  NA, 5],
    #    [3,  4,  6]]
    data = np.array([[1, 2, 5], [3, 4, 6]], dtype=np.int64)
    mask = np.array([[False, True, False], [False, False, False]], dtype=bool)
    return IntegerArray._simple_new(data, mask)


# ---------------------------------------------------------------------------
# _rank
# ---------------------------------------------------------------------------


class TestRank2D:
    def test_rank_2d_axis0(self, int_arr2d):
        result = int_arr2d._rank(axis=0, method="average")
        assert result.shape == (2, 3)
        # col 0: values [1, 3] -> ranks [1, 2]
        assert result[0, 0] == 1.0
        assert result[1, 0] == 2.0
        # col 1: values [NA, 4] -> ranks [NA, 1]
        assert pd.isna(result[0, 1])
        assert result[1, 1] == 1.0
        # col 2: values [5, 6] -> ranks [1, 2]
        assert result[0, 2] == 1.0
        assert result[1, 2] == 2.0

    def test_rank_2d_preserves_dtype(self, int_arr2d):
        result = int_arr2d._rank(axis=0, method="average")
        assert isinstance(result, FloatingArray)

    def test_rank_2d_min_method(self, int_arr2d):
        result = int_arr2d._rank(axis=0, method="min")
        assert isinstance(result, IntegerArray)
        assert result[0, 0] == 1
        assert result[1, 0] == 2

    def test_rank_2d_na_option_top(self, int_arr2d):
        result = int_arr2d._rank(axis=0, method="average", na_option="top")
        # col 1: NA gets rank 1, 4 gets rank 2
        assert result[0, 1] == 1.0
        assert result[1, 1] == 2.0
        # No NAs in result mask
        assert not result._mask.any()

    def test_rank_2d_descending(self, int_arr2d):
        result = int_arr2d._rank(axis=0, method="average", ascending=False)
        # col 0: [1, 3] desc -> ranks [2, 1]
        assert result[0, 0] == 2.0
        assert result[1, 0] == 1.0

    def test_rank_2d_all_na_column(self):
        data = np.array([[1, 2], [3, 4]], dtype=np.int64)
        mask = np.array([[False, True], [False, True]], dtype=bool)
        arr = IntegerArray._simple_new(data, mask)
        result = arr._rank(axis=0, method="average")
        # col 0: [1, 3] -> [1, 2]
        assert result[0, 0] == 1.0
        assert result[1, 0] == 2.0
        # col 1: all NA
        assert pd.isna(result[0, 1])
        assert pd.isna(result[1, 1])

    def test_rank_2d_boolean(self):
        data = np.array([[True, False], [False, True]], dtype=bool)
        mask = np.zeros((2, 2), dtype=bool)
        arr = BooleanArray._simple_new(data, mask)
        result = arr._rank(axis=0, method="average")
        assert result.shape == (2, 2)

    def test_rank_2d_axis1_raises(self, int_arr2d):
        with pytest.raises(NotImplementedError):
            int_arr2d._rank(axis=1)


# ---------------------------------------------------------------------------
# _quantile
# ---------------------------------------------------------------------------


class TestQuantile2D:
    def test_quantile_2d_no_na(self):
        data = np.array([[1, 3], [2, 4]], dtype=np.int64)
        mask = np.zeros((2, 2), dtype=bool)
        arr = IntegerArray._simple_new(data, mask)
        result = arr._quantile(np.array([0.5]), "linear")
        assert result.shape == (2, 1)
        # mask should be all-False (no NAs)
        assert not result._mask.any()

    def test_quantile_2d_with_na(self):
        # row 0: [1, NA] -> quantile computed only from 1
        # row 1: [3, 4] -> quantile is 3.5
        data = np.array([[1, 2], [3, 4]], dtype=np.int64)
        mask = np.array([[False, True], [False, False]], dtype=bool)
        arr = IntegerArray._simple_new(data, mask)
        result = arr._quantile(np.array([0.5]), "linear")
        assert result.shape == (2, 1)
        # Row 0 has at least one non-NA value -> not masked
        assert not result._mask[0, 0]
        # Row 1 has no NAs -> not masked
        assert not result._mask[1, 0]

    def test_quantile_2d_all_na_row(self):
        # row 0: all NA, row 1: [3, 4]
        # Use float dtype to avoid RuntimeWarning on NaN -> int cast
        data = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
        mask = np.array([[True, True], [False, False]], dtype=bool)
        arr = FloatingArray._simple_new(data, mask)
        result = arr._quantile(np.array([0.5]), "linear")
        assert result.shape == (2, 1)
        # Row 0 is all-NA -> masked
        assert result._mask[0, 0]
        # Row 1 is not all-NA -> not masked
        assert not result._mask[1, 0]

    def test_quantile_2d_multiple_qs(self):
        data = np.array([[10, 20], [30, 40]], dtype=np.int64)
        mask = np.zeros((2, 2), dtype=bool)
        arr = IntegerArray._simple_new(data, mask)
        qs = np.array([0.25, 0.5, 0.75])
        result = arr._quantile(qs, "linear")
        # shape: (nrows=2, nqs=3)
        assert result.shape == (2, 3)
        assert not result._mask.any()

    def test_quantile_2d_float(self):
        data = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
        mask = np.array([[False, True], [False, False]], dtype=bool)
        arr = FloatingArray._simple_new(data, mask)
        result = arr._quantile(np.array([0.5]), "linear")
        assert result.shape == (2, 1)
        assert isinstance(result, FloatingArray)
