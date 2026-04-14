"""
Tests for 2D support in BaseMaskedArray methods.

These tests exercise code paths that were previously gated with
``ndim > 1`` checks or assumed 1D shapes. Each section corresponds to
a method that was made 2D-compatible.
"""

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
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


@pytest.fixture
def float_arr2d():
    # shape (2, 3):
    #   [[1.0,  NA, 5.0],
    #    [3.0, 4.0, 6.0]]
    data = np.array([[1.0, 2.0, 5.0], [3.0, 4.0, 6.0]], dtype=np.float64)
    mask = np.array([[False, True, False], [False, False, False]], dtype=bool)
    return FloatingArray._simple_new(data, mask)


@pytest.fixture
def bool_arr2d():
    # shape (2, 3):
    #   [[True,  NA,  False],
    #    [False, True, True ]]
    data = np.array([[True, False, False], [False, True, True]], dtype=bool)
    mask = np.array([[False, True, False], [False, False, False]], dtype=bool)
    return BooleanArray._simple_new(data, mask)


# ---------------------------------------------------------------------------
# __array_ufunc__ (PR 1-A)
# ---------------------------------------------------------------------------


class TestArrayUfunc2D:
    def test_add_2d(self, int_arr2d):
        result = np.add(int_arr2d, int_arr2d)
        assert result.shape == (2, 3)
        assert isinstance(result, IntegerArray)
        # Non-NA positions: doubled
        assert result[1, 0] == 6
        assert result[1, 1] == 8
        # NA propagation
        assert pd.isna(result[0, 1])
        assert not pd.isna(result[0, 0])

    def test_ufunc_2d_with_ndarray(self, int_arr2d):
        other = np.array([[10, 20, 30], [40, 50, 60]])
        result = np.add(int_arr2d, other)
        assert result.shape == (2, 3)
        assert result[0, 0] == 11
        assert result[1, 2] == 66
        assert pd.isna(result[0, 1])

    def test_ufunc_2d_with_scalar(self, int_arr2d):
        result = np.multiply(int_arr2d, 2)
        assert result.shape == (2, 3)
        assert result[0, 0] == 2
        assert result[1, 1] == 8
        assert pd.isna(result[0, 1])

    def test_ufunc_2d_mask_propagation(self):
        data1 = np.array([[1, 2], [3, 4]], dtype=np.int64)
        mask1 = np.array([[True, False], [False, False]], dtype=bool)
        arr1 = IntegerArray._simple_new(data1, mask1)

        data2 = np.array([[5, 6], [7, 8]], dtype=np.int64)
        mask2 = np.array([[False, True], [False, False]], dtype=bool)
        arr2 = IntegerArray._simple_new(data2, mask2)

        result = np.add(arr1, arr2)
        assert result.shape == (2, 2)
        # mask union: (0,0) and (0,1) are masked
        assert pd.isna(result[0, 0])
        assert pd.isna(result[0, 1])
        assert not pd.isna(result[1, 0])
        assert result[1, 0] == 10
        assert result[1, 1] == 12


# ---------------------------------------------------------------------------
# _arith_method (PR 1-B)
# ---------------------------------------------------------------------------


class TestArithMethod2D:
    def test_add_scalar(self, int_arr2d):
        result = int_arr2d + 10
        assert result.shape == (2, 3)
        assert result[0, 0] == 11
        assert result[1, 1] == 14
        assert pd.isna(result[0, 1])

    def test_add_same_shape(self, int_arr2d):
        result = int_arr2d + int_arr2d
        assert result.shape == (2, 3)
        assert result[0, 0] == 2
        assert result[1, 2] == 12
        assert pd.isna(result[0, 1])

    def test_sub_2d_ndarray(self, int_arr2d):
        other = np.ones((2, 3), dtype=np.int64)
        result = int_arr2d - other
        assert result.shape == (2, 3)
        assert result[0, 0] == 0
        assert result[1, 0] == 2
        assert pd.isna(result[0, 1])

    def test_mul_float_2d(self, float_arr2d):
        result = float_arr2d * 2.0
        assert result.shape == (2, 3)
        assert result[0, 0] == 2.0
        assert result[1, 1] == 8.0
        assert pd.isna(result[0, 1])

    def test_1d_rejects_2d_operand(self):
        arr = pd.array([1, 2, 3, 4], dtype="Int64")
        msg = "can only perform ops with 1-d structures"
        with pytest.raises(NotImplementedError, match=msg):
            arr + np.ones((1, 4), dtype=np.int64)

    def test_2d_rejects_3d_operand(self, int_arr2d):
        msg = "can only perform ops with 1-d structures"
        with pytest.raises(NotImplementedError, match=msg):
            int_arr2d + np.ones((2, 3, 1), dtype=np.int64)

    def test_pow_mask_handling_2d(self):
        # 1 ** x should be unmasked; x ** 0 should be unmasked
        data = np.array([[1, 2], [3, 4]], dtype=np.int64)
        mask = np.array([[False, True], [True, False]], dtype=bool)
        arr = IntegerArray._simple_new(data, mask)

        result = arr**0
        # anything ** 0 == 1, regardless of mask
        assert result[0, 0] == 1
        assert result[0, 1] == 1
        assert result[1, 0] == 1
        assert result[1, 1] == 1


# ---------------------------------------------------------------------------
# _cmp_method (PR 1-C)
# ---------------------------------------------------------------------------


class TestCmpMethod2D:
    def test_eq_scalar(self, int_arr2d):
        result = int_arr2d == 3
        assert isinstance(result, BooleanArray)
        assert result.shape == (2, 3)
        assert result[1, 0] is np.bool_(True)
        assert result[0, 0] is np.bool_(False)
        assert pd.isna(result[0, 1])

    def test_eq_same_shape(self, int_arr2d):
        result = int_arr2d == int_arr2d
        assert result.shape == (2, 3)
        # non-NA positions are equal to themselves
        assert result[0, 0] is np.bool_(True)
        assert result[1, 2] is np.bool_(True)
        # NA == NA -> NA
        assert pd.isna(result[0, 1])

    def test_lt_2d_ndarray(self, int_arr2d):
        other = np.array([[2, 3, 4], [2, 3, 4]], dtype=np.int64)
        result = int_arr2d < other
        assert result.shape == (2, 3)
        assert result[0, 0] is np.bool_(True)  # 1 < 2
        assert result[1, 0] is np.bool_(False)  # 3 < 2 = False
        assert pd.isna(result[0, 1])

    def test_1d_rejects_2d_operand(self):
        arr = pd.array([1, 2, 3, 4], dtype="Int64")
        msg = "can only perform ops with 1-d structures"
        with pytest.raises(NotImplementedError, match=msg):
            arr == np.ones((1, 4), dtype=np.int64)

    def test_shape_mismatch_raises(self, int_arr2d):
        other = np.ones((2, 2), dtype=np.int64)
        with pytest.raises(ValueError, match="Lengths must match"):
            int_arr2d == other

    def test_cmp_na(self, int_arr2d):
        result = int_arr2d == pd.NA
        assert result.shape == (2, 3)
        assert pd.isna(result[0, 0])
        assert pd.isna(result[1, 1])


# ---------------------------------------------------------------------------
# _rank (PR 1-D)
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
# _quantile (PR 1-E)
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


# ---------------------------------------------------------------------------
# fillna fast path (PR 1-F)
# ---------------------------------------------------------------------------


class TestFillna2D:
    def test_fillna_scalar_2d(self, int_arr2d):
        result = int_arr2d.fillna(99)
        assert result.shape == (2, 3)
        # The NA at (0,1) should be filled
        assert result[0, 1] == 99
        # Non-NA values unchanged
        assert result[0, 0] == 1
        assert result[1, 1] == 4
        # No remaining NAs
        assert not result._mask.any()

    def test_fillna_scalar_2d_no_na(self):
        data = np.array([[1, 2], [3, 4]], dtype=np.int64)
        mask = np.zeros((2, 2), dtype=bool)
        arr = IntegerArray._simple_new(data, mask)
        result = arr.fillna(99)
        tm.assert_extension_array_equal(result, arr)

    def test_fillna_scalar_2d_float(self, float_arr2d):
        result = float_arr2d.fillna(0.0)
        assert result.shape == (2, 3)
        assert result[0, 1] == 0.0
        assert not result._mask.any()

    def test_fillna_limit_2d(self):
        # 3 rows, 2 cols, multiple NAs per column
        data = np.array([[1, 2], [3, 4], [5, 6]], dtype=np.int64)
        mask = np.array([[True, False], [True, True], [False, True]], dtype=bool)
        arr = IntegerArray._simple_new(data, mask)
        result = arr.fillna(0, limit=1)
        # limit=1: only the first NA per column is filled
        assert result[0, 0] == 0  # first NA in col 0 -> filled
        assert pd.isna(result[1, 0])  # second NA in col 0 -> not filled
        assert not pd.isna(result[2, 0])  # was not NA
        assert result[1, 1] == 0  # first NA in col 1 -> filled
        assert pd.isna(result[2, 1])  # second NA in col 1 -> not filled

    def test_fillna_copy_false_2d(self, int_arr2d):
        result = int_arr2d.fillna(99, copy=False)
        assert result.shape == (2, 3)
        assert result[0, 1] == 99


# ---------------------------------------------------------------------------
# interpolate (PR 1-G)
# ---------------------------------------------------------------------------


class TestInterpolate2D:
    def test_interpolate_2d_linear(self):
        # shape (2, 4): two "columns" (block rows), 4 "rows" each
        # row 0: [1.0, NA, NA, 4.0]  -> after linear interp: [1, 2, 3, 4]
        # row 1: [10.0, NA, 30.0, NA] -> forward fill at end: [10, 20, 30, 30]
        data = np.array(
            [[1.0, 99.0, 99.0, 4.0], [10.0, 99.0, 30.0, 99.0]],
            dtype=np.float64,
        )
        mask = np.array(
            [[False, True, True, False], [False, True, False, True]], dtype=bool
        )
        arr = FloatingArray._simple_new(data, mask)
        result = arr.interpolate(
            method="linear",
            axis=0,
            index=pd.RangeIndex(4),
            limit=None,
            limit_direction="forward",
            limit_area=None,
            copy=True,
        )
        assert result.shape == (2, 4)
        assert isinstance(result, FloatingArray)
        # row 0: linear interpolation between 1.0 and 4.0
        assert result[0, 0] == 1.0
        assert result[0, 1] == 2.0
        assert result[0, 2] == 3.0
        assert result[0, 3] == 4.0
        # row 1: linear interpolation between 10.0 and 30.0
        assert result[1, 0] == 10.0
        assert result[1, 1] == 20.0
        assert result[1, 2] == 30.0

    def test_interpolate_2d_integer_returns_float(self):
        data = np.array([[1, 99, 3], [10, 99, 30]], dtype=np.int64)
        mask = np.array([[False, True, False], [False, True, False]], dtype=bool)
        arr = IntegerArray._simple_new(data, mask)
        result = arr.interpolate(
            method="linear",
            axis=0,
            index=pd.RangeIndex(3),
            limit=None,
            limit_direction="forward",
            limit_area=None,
            copy=True,
        )
        assert isinstance(result, FloatingArray)
        assert result[0, 1] == 2.0
        assert result[1, 1] == 20.0

    def test_interpolate_2d_no_na(self):
        data = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
        mask = np.zeros((2, 2), dtype=bool)
        arr = FloatingArray._simple_new(data, mask)
        result = arr.interpolate(
            method="linear",
            axis=0,
            index=pd.RangeIndex(2),
            limit=None,
            limit_direction="forward",
            limit_area=None,
            copy=True,
        )
        # No NAs, so values unchanged
        assert result[0, 0] == 1.0
        assert result[1, 1] == 4.0
