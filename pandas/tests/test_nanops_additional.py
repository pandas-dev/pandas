"""
Additional unit test cases for pandas nanops module - edge cases and boundary conditions.
These tests are separate from the baseline test suite to avoid interference.
"""
import numpy as np
import pytest

from pandas.core import nanops


def test_nansum_empty_array_edge_cases():
    """Test nansum behavior with empty arrays and different dtypes - edge case coverage."""
    # Empty float array should return 0.0
    empty_float = np.array([], dtype=np.float64)
    result = nanops.nansum(empty_float)
    assert result == 0.0
    
    # Empty integer array should return 0
    empty_int = np.array([], dtype=np.int64)
    result = nanops.nansum(empty_int)
    assert result == 0
    
    # Empty array with min_count requirement should return NaN
    result = nanops.nansum(empty_float, min_count=1)
    assert np.isnan(result)


def test_nanmean_mask_edge_cases():
    """Test nanmean with different mask scenarios - uncovered mask logic."""
    values = np.array([1.0, 2.0, 3.0, 4.0])
    
    # All values masked should return NaN
    all_masked = np.array([True, True, True, True])
    result = nanops.nanmean(values, mask=all_masked)
    assert np.isnan(result)
    
    # Partial mask should return mean of unmasked values
    partial_mask = np.array([True, False, False, True])
    result = nanops.nanmean(values, mask=partial_mask)
    assert result == 2.5  # mean of [2.0, 3.0]
    
    # No mask should return regular mean
    result = nanops.nanmean(values, mask=None)
    assert result == 2.5  # mean of [1.0, 2.0, 3.0, 4.0]


def test_nanvar_ddof_boundary_conditions():
    """Test nanvar with boundary ddof values - statistical edge cases."""
    values = np.array([1.0, 2.0, 3.0])
    
    # ddof equal to sample size should return NaN
    result = nanops.nanvar(values, ddof=3)
    assert np.isnan(result)
    
    # ddof greater than sample size should return NaN
    result = nanops.nanvar(values, ddof=4)
    assert np.isnan(result)
    
    # ddof = 0 should work normally
    result = nanops.nanvar(values, ddof=0)
    assert not np.isnan(result) and not np.isinf(result)


def test_nanargmax_nanargmin_error_conditions():
    """Test error handling in nanargmax/nanargmin - error path coverage."""
    # All NaN array should raise ValueError
    all_nan = np.array([np.nan, np.nan, np.nan])
    
    with pytest.raises(ValueError):
        nanops.nanargmax(all_nan)
    
    with pytest.raises(ValueError):
        nanops.nanargmin(all_nan)
    
    # Empty array should raise ValueError
    empty_array = np.array([])
    
    with pytest.raises(ValueError):
        nanops.nanargmax(empty_array)
    
    with pytest.raises(ValueError):
        nanops.nanargmin(empty_array)


def test_nanskew_nankurt_insufficient_samples():
    """Test skewness/kurtosis with insufficient sample sizes - statistical boundary cases."""
    # Single value should return NaN for skewness
    single_value = np.array([1.0])
    result = nanops.nanskew(single_value)
    assert np.isnan(result)
    
    # Two values should return NaN for kurtosis (need at least 4)
    two_values = np.array([1.0, 2.0])
    result = nanops.nankurt(two_values)
    assert np.isnan(result)
    
    # All NaN should return NaN
    all_nan = np.array([np.nan, np.nan, np.nan])
    result = nanops.nanskew(all_nan)
    assert np.isnan(result)
    
    result = nanops.nankurt(all_nan)
    assert np.isnan(result)