"""
Unit tests for private methods in pandas.core.nanops_numba module.

This module tests only the private methods (prefixed with underscore).
"""

import numpy as np
import pytest
from numba.typed import List as NumbaList

from pandas.core.nanops_numba import (
    MIN_INT,
    NumbaReductionOps,
    _cast_to_timelike,
    _chunk_arr_into_arr_list,
    _get_initial_value,
    _nb_reduce_arr_list_in_parallel,
    _nb_reduce_single_arr,
    _nanvar_std_sem,
    _nullify_below_mincount,
    _reduce_chunked_results,
    _reduce_empty_array,
    nb_reduce,
)

import pandas._testing as tm


class TestGetInitialValue:
    """Test the _get_initial_value private function."""

    @pytest.fixture
    def float_array_with_nan(self):
        return np.array([np.nan, 2.0, 3.0, np.nan])

    @pytest.fixture
    def valid_array(self):
        return np.array([1.0, 2.0, 3.0])

    def test_skipna_true_with_leading_nulls(self, float_array_with_nan):
        index, value = _get_initial_value(float_array_with_nan, skipna=True)
        assert index == 1
        assert value == 2.0

    def test_skipna_false_with_leading_null(self, float_array_with_nan):
        index, value = _get_initial_value(float_array_with_nan, skipna=False)
        assert index == 0
        assert np.isnan(value)

    def test_skipna_false_with_valid_first(self, valid_array):
        index, value = _get_initial_value(valid_array, skipna=False)
        assert index == 0
        assert value == 1.0

    def test_all_null_values(self):
        arr = np.array([np.nan, np.nan, np.nan])
        index, value = _get_initial_value(arr, skipna=True)
        assert index == -1
        assert np.isnan(value)

    def test_with_mask_skipna_true(self):
        arr = np.array([1.0, 2.0, 3.0])
        mask = np.array([True, False, False])  # First element masked
        index, value = _get_initial_value(arr, skipna=True, mask=mask)
        assert index == 1
        assert value == 2.0

    def test_with_mask_skipna_false(self):
        arr = np.array([1.0, 2.0, 3.0])
        mask = np.array([True, False, False])  # First element masked
        index, value = _get_initial_value(arr, skipna=False, mask=mask)
        assert index == 0
        assert np.isnan(value)


class TestNbReduceSingleArr:
    """Test _nb_reduce_single_arr private function."""

    @pytest.fixture
    def sample_array(self):
        return np.array([1.0, 2.0, 3.0, 4.0, 5.0])

    @pytest.fixture
    def array_with_nan(self):
        return np.array([1.0, np.nan, 3.0, 4.0, np.nan])

    def test_sum_no_nulls(self, sample_array):
        result, count = _nb_reduce_single_arr(
            NumbaReductionOps.sum, sample_array, skipna=True
        )
        assert result == 15.0
        assert count == 5

    def test_sum_with_nans_skipna_true(self, array_with_nan):
        result, count = _nb_reduce_single_arr(
            NumbaReductionOps.sum, array_with_nan, skipna=True
        )
        assert result == 8.0  # 1 + 3 + 4
        assert count == 3

    def test_sum_with_nans_skipna_false(self, array_with_nan):
        result, count = _nb_reduce_single_arr(
            NumbaReductionOps.sum, array_with_nan, skipna=False
        )
        assert np.isnan(result)
        assert count == 1

    def test_min_operation(self, sample_array):
        result, count = _nb_reduce_single_arr(
            NumbaReductionOps.min, sample_array, skipna=True
        )
        assert result == 1.0
        assert count == 5

    def test_max_operation(self, sample_array):
        result, count = _nb_reduce_single_arr(
            NumbaReductionOps.max, sample_array, skipna=True
        )
        assert result == 5.0
        assert count == 5

    def test_with_mask(self):
        arr = np.array([1.0, 2.0, 3.0, 4.0])
        mask = np.array([False, True, False, True])  # mask 2.0 and 4.0
        result, count = _nb_reduce_single_arr(
            NumbaReductionOps.sum, arr, skipna=True, mask=mask
        )
        assert result == 4.0  # 1 + 3
        assert count == 2

    def test_find_initial_value_false(self, sample_array):
        result, count = _nb_reduce_single_arr(
            NumbaReductionOps.sum, sample_array, find_initial_value=False
        )
        # Should start with 0 and add all values
        assert result == 15.0
        assert count == 5


class TestNullifyBelowMincount:
    """Test _nullify_below_mincount private function."""

    def test_float_array(self):
        result = np.array([1.0, 2.0, 3.0])
        count = np.array([2, 1, 3])
        min_count = 2

        modified = _nullify_below_mincount(result, count, min_count)

        expected = np.array([1.0, np.nan, 3.0])
        tm.assert_numpy_array_equal(modified, expected)

    def test_int_array(self):
        result = np.array([1, 2, 3], dtype=np.int64)
        count = np.array([2, 1, 3])
        min_count = 2

        modified = _nullify_below_mincount(result, count, min_count)

        expected = np.array([1, MIN_INT, 3])
        tm.assert_numpy_array_equal(modified, expected)


class TestReduceEmptyArray:
    """Test _reduce_empty_array private function."""

    def test_1d_empty_array(self):
        arr = np.array([], dtype=np.float64)
        result, count = _reduce_empty_array("sum", arr, axis=None, min_count=0)
        assert result == 0.0
        assert count == 0

    def test_1d_empty_array_with_min_count(self):
        arr = np.array([], dtype=np.float64)
        result, count = _reduce_empty_array("sum", arr, axis=None, min_count=1)
        assert np.isnan(result)
        assert count == 0

    def test_2d_empty_array_axis_0(self):
        arr = np.array([[], []], dtype=np.float64)
        result, count = _reduce_empty_array("sum", arr, axis=0, min_count=0)
        # Empty array along axis 0 means no elements to reduce
        tm.assert_numpy_array_equal(result, np.array([]))
        tm.assert_numpy_array_equal(count, np.array([]))


class TestChunkArrIntoArrList:
    """Test _chunk_arr_into_arr_list private function."""

    def test_1d_array_single_thread(self):
        arr = np.array([1, 2, 3, 4, 5])
        arr_list, _, final_length = _chunk_arr_into_arr_list(
            arr, multi_threading=False, axis=None
        )
        assert len(arr_list) == 1
        tm.assert_numpy_array_equal(arr_list[0], arr)
        assert final_length == 0

    def test_2d_array_axis_0(self):
        arr = np.array([[1, 2], [3, 4], [5, 6]])
        arr_list, _, final_length = _chunk_arr_into_arr_list(
            arr, multi_threading=False, axis=0
        )
        # For axis=0, returns transposed array, so 2 columns (arrays)
        assert len(arr_list) == 2
        assert final_length == 2
        # arr_list is the transpose of the original array
        expected_transposed = arr.T
        for i in range(len(arr_list)):
            tm.assert_numpy_array_equal(arr_list[i], expected_transposed[i])

    def test_2d_array_axis_1(self):
        arr = np.array([[1, 2], [3, 4], [5, 6]])
        arr_list, _, final_length = _chunk_arr_into_arr_list(
            arr, multi_threading=False, axis=1
        )
        # For axis=1, returns original array, so 3 rows
        assert len(arr_list) == 3
        assert final_length == 3
        # arr_list is the original array
        for i in range(len(arr_list)):
            tm.assert_numpy_array_equal(arr_list[i], arr[i])


class TestNbReduceArrListInParallel:
    """Test _nb_reduce_arr_list_in_parallel private function."""

    @pytest.fixture
    def array_list(self):
        # Create a NumbaList of arrays for parallel processing
        arr_list = NumbaList()
        arr_list.append(np.array([1.0, 2.0, 3.0]))
        arr_list.append(np.array([4.0, 5.0, 6.0]))
        arr_list.append(np.array([7.0, 8.0, 9.0]))
        return arr_list

    def test_parallel_sum_without_mask(self, array_list):
        target = np.zeros(len(array_list), dtype=np.float64)
        result, counts = _nb_reduce_arr_list_in_parallel(
            NumbaReductionOps.sum, array_list, target, mask_list=None, skipna=True
        )

        expected_results = np.array([6.0, 15.0, 24.0])  # [1+2+3, 4+5+6, 7+8+9]
        expected_counts = np.array([3, 3, 3])

        tm.assert_numpy_array_equal(result, expected_results)
        tm.assert_numpy_array_equal(counts, expected_counts)

    def test_parallel_with_mask(self):
        # Create array list with some elements that should be masked
        arr_list = NumbaList()
        arr_list.append(np.array([1.0, 2.0, 3.0]))
        arr_list.append(np.array([4.0, 5.0, 6.0]))

        # Create corresponding mask list
        mask_list = NumbaList()
        mask_list.append(np.array([False, True, False]))  # Mask middle element
        mask_list.append(np.array([True, False, False]))   # Mask first element

        target = np.zeros(len(arr_list), dtype=np.float64)
        result, counts = _nb_reduce_arr_list_in_parallel(
            NumbaReductionOps.sum, arr_list, target, mask_list=mask_list, skipna=True
        )

        expected_results = np.array([4.0, 11.0])  # [1+3, 5+6]
        expected_counts = np.array([2, 2])

        tm.assert_numpy_array_equal(result, expected_results)
        tm.assert_numpy_array_equal(counts, expected_counts)


class TestReduceChunkedResults:
    """Test _reduce_chunked_results private function."""

    def test_single_chunk_reduction(self):
        # Test when final_length == 0 (reduce both axes)
        chunk_results = np.array([1.0, 2.0, 3.0])
        counts = np.array([2, 2, 2])
        final_length = 0
        return_dtype = np.dtype("float64")

        result, count = _reduce_chunked_results(
            "sum", chunk_results, counts, final_length, return_dtype,
            skipna=True, find_initial_value=True
        )

        # Should reduce the chunk_results array itself
        expected_result = np.array([6.0])  # 1 + 2 + 3
        expected_count = np.array([6])     # 2 + 2 + 2

        tm.assert_numpy_array_equal(result, expected_result)
        tm.assert_numpy_array_equal(count, expected_count)

    def test_no_chunking_needed(self):
        # Test when chunk_results and counts are already in final form
        chunk_results = np.array([10.0, 20.0])
        counts = np.array([3, 4])
        final_length = 2
        return_dtype = np.dtype("float64")

        result, count = _reduce_chunked_results(
            "sum", chunk_results, counts, final_length, return_dtype,
            skipna=True, find_initial_value=True
        )

        # Should return results as-is (no further reduction needed)
        tm.assert_numpy_array_equal(result, chunk_results)
        tm.assert_numpy_array_equal(count, counts)


class TestCastToTimelike:
    """Test _cast_to_timelike private function."""

    def test_cast_to_datetime(self):
        arr = np.array([1000000000000000000, 2000000000000000000], dtype=np.int64)
        result = _cast_to_timelike(arr, np.dtype("datetime64[ns]"))
        assert result.dtype == np.dtype("datetime64[ns]")

    def test_cast_to_timedelta(self):
        arr = np.array([1000000000, 2000000000], dtype=np.int64)
        result = _cast_to_timelike(arr, np.dtype("timedelta64[ns]"))
        assert result.dtype == np.dtype("timedelta64[ns]")

    def test_cast_with_nan_values(self):
        arr = np.array([1000000000000000000, np.nan], dtype=np.float64)
        result = _cast_to_timelike(arr, np.dtype("datetime64[ns]"))
        assert result.dtype == np.dtype("datetime64[ns]")
        assert np.isnat(result[1])


class TestNanvarStdSem:
    """Test the _nanvar_std_sem private function."""

    @pytest.fixture
    def sample_data(self):
        return np.array([1.0, 2.0, 3.0, 4.0, 5.0])

    def test_variance_calculation(self, sample_data):
        result = _nanvar_std_sem(sample_data)
        # Sample variance with ddof=1
        # mean = 3.0, deviations = [-2,-1,0,1,2], sum of squares = 10
        # variance = 10/4 = 2.5
        expected = 2.5
        tm.assert_almost_equal(result, expected)

    def test_standard_deviation(self, sample_data):
        result = _nanvar_std_sem(sample_data, std=True)
        expected = np.sqrt(2.5)
        tm.assert_almost_equal(result, expected)

    def test_standard_error_of_mean(self, sample_data):
        result = _nanvar_std_sem(sample_data, sem=True)
        expected = np.sqrt(2.5) / np.sqrt(5)  # std / sqrt(n)
        tm.assert_almost_equal(result, expected)

    def test_with_ddof_0(self, sample_data):
        result = _nanvar_std_sem(sample_data, ddof=0)
        # Population variance: sum((x - mean)^2) / n = 10/5 = 2.0
        expected = 2.0
        tm.assert_almost_equal(result, expected)

    def test_insufficient_data_for_ddof(self):
        arr = np.array([1.0])  # Only one value
        result = _nanvar_std_sem(arr, ddof=1)  # ddof=1 requires at least 2 values
        assert np.isnan(result)

    def test_with_nan_values(self):
        arr = np.array([1.0, np.nan, 3.0, 4.0, np.nan])
        result = _nanvar_std_sem(arr, skipna=True)
        # Should not raise error and return finite result
        assert np.isfinite(result)

    def test_complex_array(self):
        arr = np.array([1+2j, 3+4j])
        result = _nanvar_std_sem(arr)
        # Should handle complex numbers by processing real and imag parts
        assert np.isfinite(result)

    def test_datetime_array(self):
        arr = np.array(
            ["2020-01-01", "2020-01-02", "2020-01-03"], dtype="datetime64[ns]"
        )
        result = _nanvar_std_sem(arr, std=True)
        # Should return timedelta type for datetime input with std=True
        assert result.dtype.kind == "m"  # timedelta type

    def test_2d_array_calculation(self):
        arr = np.array([[1.0, 2.0], [3.0, 4.0]])
        result = _nanvar_std_sem(arr, axis=0)
        # Variance along axis 0: var([1,3]) and var([2,4])
        expected = np.array([1.0, 1.0])  # Each column has variance 1.0 with ddof=1
        tm.assert_numpy_array_equal(result, expected)


class TestNbReduceMultithreading:
    """Test nb_reduce with large arrays to trigger multi-threading."""

    @pytest.fixture
    def large_2d_array(self):
        """Create a large 2D array to trigger multi-threading."""
        # Create array large enough to trigger multi-threading (> 1e6 elements)
        rows, cols = 2000, 600  # 1.2M elements
        rng = np.random.default_rng(42)  # For reproducible tests
        arr = rng.standard_normal((rows, cols)).astype(np.float64)
        return arr

    @pytest.fixture
    def large_2d_array_with_nans(self):
        """Create a large 2D array with some NaN values."""
        rows, cols = 2000, 600
        rng = np.random.default_rng(42)
        arr = rng.standard_normal((rows, cols)).astype(np.float64)
        # Add some NaN values at random positions
        nan_mask = rng.random((rows, cols)) < 0.05  # 5% NaN values
        arr[nan_mask] = np.nan
        return arr

    def test_nb_reduce_sum_axis_none_multithreaded(self, large_2d_array):
        """Test sum reduction with axis=None on large array (multi-threaded)."""
        result, count = nb_reduce(
            "sum", large_2d_array, axis=None, multi_threading=True
        )

        # Compare with numpy result
        expected = np.sum(large_2d_array)
        expected_count = large_2d_array.size

        tm.assert_almost_equal(result, expected, rtol=1e-10)
        assert count == expected_count

    def test_nb_reduce_sum_axis_0_multithreaded(self, large_2d_array):
        """Test sum reduction along axis 0 on large array (multi-threaded)."""
        result, count = nb_reduce("sum", large_2d_array, axis=0, multi_threading=True)

        # Compare with numpy result
        expected = np.sum(large_2d_array, axis=0)
        expected_count = np.full(large_2d_array.shape[1], large_2d_array.shape[0])

        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(count, expected_count)

    def test_nb_reduce_sum_axis_1_multithreaded(self, large_2d_array):
        """Test sum reduction along axis 1 on large array (multi-threaded)."""
        result, count = nb_reduce("sum", large_2d_array, axis=1, multi_threading=True)

        # Compare with numpy result
        expected = np.sum(large_2d_array, axis=1)
        expected_count = np.full(large_2d_array.shape[0], large_2d_array.shape[1])

        np.testing.assert_array_almost_equal(result, expected)
        tm.assert_numpy_array_equal(count, expected_count)

    def test_nb_reduce_min_axis_none_multithreaded(self, large_2d_array):
        """Test min reduction with axis=None on large array (multi-threaded)."""
        result, count = nb_reduce(
            "min", large_2d_array, axis=None, multi_threading=True
        )

        # Compare with numpy result
        expected = np.min(large_2d_array)
        expected_count = large_2d_array.size

        tm.assert_almost_equal(result, expected)
        assert count == expected_count

    def test_nb_reduce_min_axis_0_multithreaded(self, large_2d_array):
        """Test min reduction along axis 0 on large array (multi-threaded)."""
        result, count = nb_reduce("min", large_2d_array, axis=0, multi_threading=True)

        # Compare with numpy result
        expected = np.min(large_2d_array, axis=0)
        expected_count = np.full(large_2d_array.shape[1], large_2d_array.shape[0])

        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(count, expected_count)

    def test_nb_reduce_max_axis_1_multithreaded(self, large_2d_array):
        """Test max reduction along axis 1 on large array (multi-threaded)."""
        result, count = nb_reduce("max", large_2d_array, axis=1, multi_threading=True)

        # Compare with numpy result
        expected = np.max(large_2d_array, axis=1)
        expected_count = np.full(large_2d_array.shape[0], large_2d_array.shape[1])

        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(count, expected_count)

    def test_nb_reduce_with_nans_skipna_true_multithreaded(
        self, large_2d_array_with_nans
    ):
        """Test sum with NaN values and skipna=True on large array (multi-threaded)."""
        result, count = nb_reduce(
            "sum", large_2d_array_with_nans, axis=None,
            skipna=True, multi_threading=True
        )

        # Compare with numpy nansum
        expected = np.nansum(large_2d_array_with_nans)
        expected_count = np.sum(~np.isnan(large_2d_array_with_nans))

        tm.assert_almost_equal(result, expected, rtol=1e-10)
        assert count == expected_count

    def test_nb_reduce_with_nans_axis_0_multithreaded(self, large_2d_array_with_nans):
        """Test sum with NaN values along axis 0 (multi-threaded)."""
        result, count = nb_reduce("sum", large_2d_array_with_nans, axis=0,
                                  skipna=True, multi_threading=True)

        # Compare with numpy nansum
        expected = np.nansum(large_2d_array_with_nans, axis=0)
        expected_count = np.sum(~np.isnan(large_2d_array_with_nans), axis=0)

        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(count, expected_count)

    def test_nb_reduce_with_nans_axis_1_multithreaded(self, large_2d_array_with_nans):
        """Test sum with NaN values along axis 1 (multi-threaded)."""
        result, count = nb_reduce("sum", large_2d_array_with_nans, axis=1,
                                  skipna=True, multi_threading=True)

        # Compare with numpy nansum
        expected = np.nansum(large_2d_array_with_nans, axis=1)
        expected_count = np.sum(~np.isnan(large_2d_array_with_nans), axis=1)

        np.testing.assert_array_almost_equal(result, expected)
        tm.assert_numpy_array_equal(count, expected_count)

    def test_nb_reduce_single_thread_vs_multithread_consistency(self, large_2d_array):
        """Test that single-threaded and multi-threaded results are identical."""
        # Single-threaded result
        result_st, count_st = nb_reduce("sum", large_2d_array, axis=0,
                                        multi_threading=False)

        # Multi-threaded result
        result_mt, count_mt = nb_reduce("sum", large_2d_array, axis=0,
                                        multi_threading=True)

        # Results should be identical
        tm.assert_numpy_array_equal(result_st, result_mt)
        tm.assert_numpy_array_equal(count_st, count_mt)

    @pytest.mark.parametrize("op", ["sum", "min", "max"])
    @pytest.mark.parametrize("axis", [None, 0, 1])
    def test_nb_reduce_operations_multithreaded(self, large_2d_array, op, axis):
        """Test various operations with different axes on large array."""
        result, count = nb_reduce(op, large_2d_array, axis=axis, multi_threading=True)

        # Verify result shape is correct
        if axis is None:
            assert np.isscalar(result)
            assert np.isscalar(count)
        elif axis == 0:
            assert result.shape == (large_2d_array.shape[1],)
            assert count.shape == (large_2d_array.shape[1],)
        elif axis == 1:
            assert result.shape == (large_2d_array.shape[0],)
            assert count.shape == (large_2d_array.shape[0],)

        # Verify count values are reasonable
        if axis is None:
            assert count == large_2d_array.size
        elif axis == 0:
            assert np.all(count == large_2d_array.shape[0])
        elif axis == 1:
            assert np.all(count == large_2d_array.shape[1])

    def test_nb_reduce_min_count_multithreaded(self, large_2d_array_with_nans):
        """Test min_count parameter with large array (multi-threaded)."""
        min_count = 100  # Require at least 100 non-NaN values per column

        result, count = nb_reduce("sum", large_2d_array_with_nans, axis=0,
                                  skipna=True, min_count=min_count,
                                  multi_threading=True)

        # Check that columns with insufficient data are NaN
        valid_columns = count >= min_count
        assert np.all(np.isfinite(result[valid_columns]))

        # Some columns should be nullified due to min_count
        if not np.all(valid_columns):
            assert np.any(np.isnan(result[~valid_columns]))

    def test_nb_reduce_mean_axis_none_multithreaded(self, large_2d_array):
        """Test mean reduction with axis=None on large array (multi-threaded)."""
        result, count = nb_reduce("mean", large_2d_array, axis=None,
                                  multi_threading=True)

        # Compare with numpy result
        expected = np.mean(large_2d_array)
        expected_count = large_2d_array.size

        tm.assert_almost_equal(result, expected, rtol=1e-10)
        assert count == expected_count

    def test_nb_reduce_mean_axis_0_multithreaded(self, large_2d_array):
        """Test mean reduction along axis 0 on large array (multi-threaded)."""
        result, count = nb_reduce("mean", large_2d_array, axis=0,
                                  multi_threading=True)

        # Compare with numpy result
        expected = np.mean(large_2d_array, axis=0)
        expected_count = np.full(large_2d_array.shape[1], large_2d_array.shape[0])

        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(count, expected_count)

    def test_nb_reduce_mean_axis_1_multithreaded(self, large_2d_array):
        """Test mean reduction along axis 1 on large array (multi-threaded)."""
        result, count = nb_reduce("mean", large_2d_array, axis=1,
                                  multi_threading=True)

        # Compare with numpy result
        expected = np.mean(large_2d_array, axis=1)
        expected_count = np.full(large_2d_array.shape[0], large_2d_array.shape[1])

        np.testing.assert_array_almost_equal(result, expected)
        tm.assert_numpy_array_equal(count, expected_count)

    def test_nb_reduce_sum_square_axis_none_multithreaded(self, large_2d_array):
        """Test sum_square reduction with axis=None on large array."""
        result, count = nb_reduce("sum_square", large_2d_array, axis=None,
                                  multi_threading=True)

        # Compare with numpy result (sum of squares)
        expected = np.sum(large_2d_array ** 2)
        expected_count = large_2d_array.size

        tm.assert_almost_equal(result, expected, rtol=1e-10)
        assert count == expected_count

    def test_nb_reduce_sum_square_axis_0_multithreaded(self, large_2d_array):
        """Test sum_square reduction along axis 0 on large array."""
        result, count = nb_reduce("sum_square", large_2d_array, axis=0,
                                  multi_threading=True)

        # Compare with numpy result (sum of squares along axis 0)
        expected = np.sum(large_2d_array ** 2, axis=0)
        expected_count = np.full(large_2d_array.shape[1], large_2d_array.shape[0])

        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(count, expected_count)

    def test_nb_reduce_sum_square_axis_1_multithreaded(self, large_2d_array):
        """Test sum_square reduction along axis 1 on large array."""
        result, count = nb_reduce("sum_square", large_2d_array, axis=1,
                                  multi_threading=True)

        # Compare with numpy result (sum of squares along axis 1)
        expected = np.sum(large_2d_array ** 2, axis=1)
        expected_count = np.full(large_2d_array.shape[0], large_2d_array.shape[1])

        np.testing.assert_array_almost_equal(result, expected)
        tm.assert_numpy_array_equal(count, expected_count)


class TestNbReduceTimedelta64:
    """Test nb_reduce with timedelta64 dtype arrays."""

    @pytest.fixture
    def timedelta64_2d_array(self):
        """Create a 2D array of timedelta64 values."""
        # Create timedelta values in different units
        rows, cols = 1000, 300  # Large enough for multi-threading
        rng = np.random.default_rng(42)

        # Generate random days between 0 and 365
        days = rng.integers(0, 365, size=(rows, cols))
        td_array = np.array(days, dtype="timedelta64[D]")

        return td_array

    @pytest.fixture
    def timedelta64_2d_array_with_nat(self):
        """Create a 2D array of timedelta64 values with some NaT values."""
        rows, cols = 1000, 300
        rng = np.random.default_rng(42)

        # Generate random days
        days = rng.integers(0, 365, size=(rows, cols))
        td_array = np.array(days, dtype="timedelta64[D]")

        # Add some NaT (Not a Time) values
        nat_mask = rng.random((rows, cols)) < 0.05  # 5% NaT values
        td_array[nat_mask] = np.timedelta64("NaT")

        return td_array

    def test_nb_reduce_timedelta64_sum_axis_none(self, timedelta64_2d_array):
        """Test sum reduction on timedelta64 array with axis=None."""
        result, count = nb_reduce("sum", timedelta64_2d_array, axis=None,
                                  multi_threading=True)

        # Compare with numpy result
        expected = np.sum(timedelta64_2d_array)
        expected_count = timedelta64_2d_array.size

        assert result == expected
        assert count == expected_count

    def test_nb_reduce_timedelta64_sum_axis_0(self, timedelta64_2d_array):
        """Test sum reduction on timedelta64 array along axis 0."""
        result, count = nb_reduce("sum", timedelta64_2d_array, axis=0,
                                  multi_threading=True)

        # Compare with numpy result
        expected = np.sum(timedelta64_2d_array, axis=0)
        expected_count = np.full(timedelta64_2d_array.shape[1],
                               timedelta64_2d_array.shape[0])

        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(count, expected_count)

    def test_nb_reduce_timedelta64_sum_axis_1(self, timedelta64_2d_array):
        """Test sum reduction on timedelta64 array along axis 1."""
        result, count = nb_reduce("sum", timedelta64_2d_array, axis=1,
                                  multi_threading=True)

        # Compare with numpy result
        expected = np.sum(timedelta64_2d_array, axis=1)
        expected_count = np.full(timedelta64_2d_array.shape[0],
                               timedelta64_2d_array.shape[1])

        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(count, expected_count)

    def test_nb_reduce_timedelta64_min_max(self, timedelta64_2d_array):
        """Test min/max reduction on timedelta64 array."""
        # Test min
        result_min, count_min = nb_reduce("min", timedelta64_2d_array, axis=None,
                                         multi_threading=True)
        expected_min = np.min(timedelta64_2d_array)
        assert result_min == expected_min
        assert count_min == timedelta64_2d_array.size

        # Test max
        result_max, count_max = nb_reduce("max", timedelta64_2d_array, axis=None,
                                         multi_threading=True)
        expected_max = np.max(timedelta64_2d_array)
        assert result_max == expected_max
        assert count_max == timedelta64_2d_array.size

    def test_nb_reduce_timedelta64_with_nat_skipna_true(
        self, timedelta64_2d_array_with_nat
    ):
        """Test reduction on timedelta64 array with NaT values, skipna=True."""
        result, count = nb_reduce("sum", timedelta64_2d_array_with_nat, axis=None,
                                  skipna=True, multi_threading=True)

        # Compare with numpy result
        # For timedelta64 with NaT, we need to use nansum equivalent
        valid_mask = ~np.isnat(timedelta64_2d_array_with_nat)
        expected = np.sum(timedelta64_2d_array_with_nat[valid_mask])
        expected_count = np.sum(valid_mask)

        assert result == expected
        assert count == expected_count

    def test_nb_reduce_timedelta64_with_nat_skipna_false(
        self, timedelta64_2d_array_with_nat
    ):
        """Test reduction on timedelta64 array with NaT values, skipna=False."""
        result, count = nb_reduce("sum", timedelta64_2d_array_with_nat, axis=None,
                                  skipna=False, multi_threading=True)

        # When skipna=False and there are NaT values, result should be NaT
        assert np.isnat(result)

    def test_nb_reduce_timedelta64_mean_axis_0(self, timedelta64_2d_array):
        """Test mean reduction on timedelta64 array along axis 0."""
        result, count = nb_reduce("mean", timedelta64_2d_array, axis=0,
                                  multi_threading=True)

        # Compare with numpy result
        expected = np.mean(timedelta64_2d_array, axis=0)
        expected_count = np.full(timedelta64_2d_array.shape[1],
                               timedelta64_2d_array.shape[0])

        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(count, expected_count)
