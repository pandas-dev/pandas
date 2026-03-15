import operator

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import SparseArray


@pytest.mark.filterwarnings("ignore:invalid value encountered in cast:RuntimeWarning")
@pytest.mark.parametrize("fill_value", [0, np.nan])
@pytest.mark.parametrize("op", [operator.pos, operator.neg])
def test_unary_op(op, fill_value):
    arr = np.array([0, 1, np.nan, 2])
    sparray = SparseArray(arr, fill_value=fill_value)
    result = op(sparray)
    expected = SparseArray(op(arr), fill_value=op(fill_value))
    tm.assert_sp_array_equal(result, expected)


@pytest.mark.parametrize("fill_value", [True, False])
def test_invert(fill_value):
    arr = np.array([True, False, False, True])
    sparray = SparseArray(arr, fill_value=fill_value)
    result = ~sparray
    expected = SparseArray(~arr, fill_value=not fill_value)
    tm.assert_sp_array_equal(result, expected)

    result = ~pd.Series(sparray)
    expected = pd.Series(expected)
    tm.assert_series_equal(result, expected)

    result = ~pd.DataFrame({"A": sparray})
    expected = pd.DataFrame({"A": expected})
    tm.assert_frame_equal(result, expected)


class TestUnaryMethodFillValueChange:
    """Tests for _unary_method when the fill_value changes after the operation.

    GH#60179 - the old code densified the entire array; the fix operates only
    on sp_values / sp_index so these tests verify the sparse structure stays
    canonical (no stored values equal to the new fill_value).
    """

    def test_invert_bool_fill_true(self):
        # Primary path hit by notna(SparseArray): isna returns
        # SparseDtype[bool, True], then ~ flips fill to False.
        arr = SparseArray([True, True, False, True, True], fill_value=True)
        result = ~arr
        expected = SparseArray([False, False, True, False, False], fill_value=False)
        tm.assert_sp_array_equal(result, expected)
        # Only the non-fill positions should be stored
        assert result.sp_index.npoints == 1
        assert list(result.sp_values) == [True]

    def test_invert_bool_fill_false(self):
        arr = SparseArray([False, True, True, False], fill_value=False)
        result = ~arr
        expected = SparseArray([True, False, False, True], fill_value=True)
        tm.assert_sp_array_equal(result, expected)
        assert result.sp_index.npoints == 2
        assert list(result.sp_values) == [False, False]

    def test_neg_fill_value_changes(self):
        # neg changes fill_value from -1 to 1; no stored values match new fill
        arr = SparseArray([1, -1, 2, -2], fill_value=-1)
        result = -arr
        expected = SparseArray([-1, 1, -2, 2], fill_value=1)
        tm.assert_sp_array_equal(result, expected)
        # stored values are [-1, -2, 2] (none equal new fill_value=1)
        assert result.sp_index.npoints == 3

    def test_neg_stored_collapses_to_fill(self):
        # neg changes fill_value from 1 to -1; stored value -1 becomes 1
        # which matches new fill_value (-1) — wait, neg(-1)=1, new fill=-1
        # Actually: stored value 1 at some position becomes -1 == new fill.
        arr = SparseArray([1, 1, 2, 3], fill_value=1)
        result = -arr
        expected = SparseArray([-1, -1, -2, -3], fill_value=-1)
        tm.assert_sp_array_equal(result, expected)
        # sp_values were [2, 3], neg -> [-2, -3], neither equals -1
        assert result.sp_index.npoints == 2

    def test_abs_fill_value_matches_stored(self):
        # abs changes fill_value from -1 to 1; stored value 1 matches new fill
        arr = SparseArray([1, -1, 2, -2], fill_value=-1)
        result = abs(arr)
        expected = SparseArray([1, 1, 2, 2], fill_value=1)
        tm.assert_sp_array_equal(result, expected)
        # positions 0 and 1 have value 1 == fill_value, only 2 and 3 stored
        assert result.sp_index.npoints == 2

    def test_neg_no_stored_matches_fill(self):
        # neg changes fill_value from 0 to 0 — fill unchanged, fast path
        arr = SparseArray([0, 1, 0, -1], fill_value=0)
        result = -arr
        expected = SparseArray([0, -1, 0, 1], fill_value=0)
        tm.assert_sp_array_equal(result, expected)
        assert result.sp_index.npoints == 2

    def test_invert_int_fill_value_changes(self):
        # bitwise invert on int: ~0 = -1
        # sp_values=[-1, 2] at [1, 3]; ~(-1)=0, ~2=-3; neither equals new fill -1
        arr = SparseArray([0, -1, 0, 2], fill_value=0)
        result = ~arr
        expected = SparseArray([-1, 0, -1, -3], fill_value=-1)
        tm.assert_sp_array_equal(result, expected)
        assert result.sp_index.npoints == 2


class TestUnaryMethods:
    @pytest.mark.filterwarnings(
        "ignore:invalid value encountered in cast:RuntimeWarning"
    )
    def test_neg_operator(self):
        arr = SparseArray([-1, -2, np.nan, 3], fill_value=np.nan, dtype=np.int8)
        res = -arr
        exp = SparseArray([1, 2, np.nan, -3], fill_value=np.nan, dtype=np.int8)
        tm.assert_sp_array_equal(exp, res)

        arr = SparseArray([-1, -2, 1, 3], fill_value=-1, dtype=np.int8)
        res = -arr
        exp = SparseArray([1, 2, -1, -3], fill_value=1, dtype=np.int8)
        tm.assert_sp_array_equal(exp, res)

    @pytest.mark.filterwarnings(
        "ignore:invalid value encountered in cast:RuntimeWarning"
    )
    def test_abs_operator(self):
        arr = SparseArray([-1, -2, np.nan, 3], fill_value=np.nan, dtype=np.int8)
        res = abs(arr)
        exp = SparseArray([1, 2, np.nan, 3], fill_value=np.nan, dtype=np.int8)
        tm.assert_sp_array_equal(exp, res)

        arr = SparseArray([-1, -2, 1, 3], fill_value=-1, dtype=np.int8)
        res = abs(arr)
        exp = SparseArray([1, 2, 1, 3], fill_value=1, dtype=np.int8)
        tm.assert_sp_array_equal(exp, res)

    def test_invert_operator(self):
        arr = SparseArray([False, True, False, True], fill_value=False, dtype=np.bool_)
        exp = SparseArray(
            np.invert([False, True, False, True]), fill_value=True, dtype=np.bool_
        )
        res = ~arr
        tm.assert_sp_array_equal(exp, res)

        arr = SparseArray([0, 1, 0, 2, 3, 0], fill_value=0, dtype=np.int32)
        res = ~arr
        exp = SparseArray([-1, -2, -1, -3, -4, -1], fill_value=-1, dtype=np.int32)
        tm.assert_sp_array_equal(exp, res)
