"""
Additional unit test cases for pandas Series constructors - edge cases and boundary conditions.
These tests are separate from the baseline test suite to avoid interference.
"""
import numpy as np
import pytest

from pandas import Series
import pandas as pd


def test_series_constructor_invalid_key_types():
    """Test Series construction with invalid dictionary key types - error handling coverage."""
    # Test with unhashable keys should raise
    with pytest.raises(TypeError):
        Series({[1, 2]: 'value'})  # List as key should fail
    
    # Test with complex nested unhashable keys
    with pytest.raises(TypeError):
        Series({frozenset([1, {2: 3}]): 'value'})  # Nested unhashable
    
    # Test with mixed hashable and unhashable keys
    with pytest.raises(TypeError):
        Series({1: 'valid', [2, 3]: 'invalid'})


def test_series_constructor_empty_edge_cases():
    """Test Series construction with various empty inputs - boundary condition coverage."""
    # Empty list should create empty Series
    s1 = Series([])
    assert len(s1) == 0
    
    # None should create empty Series
    s2 = Series(None)
    assert len(s2) == 0
    
    # Empty dict should create empty Series
    s3 = Series({})
    assert len(s3) == 0
    
    # Empty string should create Series with single empty string
    s4 = Series('')
    assert len(s4) == 1 and s4.iloc[0] == ''


def test_series_constructor_mixed_dtype_edge_cases():
    """Test Series construction with mixed data types - dtype inference coverage."""
    # Mixed numeric and string should result in object dtype
    mixed_data = [1, 'two', 3.0, 'four']
    s = Series(mixed_data)
    assert s.dtype == object
    
    # Mixed with None values
    mixed_with_none = [1, None, 'three', 4.0]
    s2 = Series(mixed_with_none)
    assert s2.dtype == object
    
    # Boolean mixed with numeric should promote to object
    bool_numeric = [True, 1, False, 2.5]
    s3 = Series(bool_numeric)
    assert s3.dtype == object


def test_series_constructor_memory_intensive():
    """Test Series construction with large datasets - memory edge case coverage."""
    # Large array should not cause memory issues
    large_size = 100000
    large_array = np.arange(large_size)
    s = Series(large_array)
    assert len(s) == large_size
    assert s.iloc[0] == 0
    assert s.iloc[-1] == large_size - 1


def test_series_constructor_invalid_index_length():
    """Test Series construction with mismatched index lengths - validation coverage."""
    data = [1, 2, 3, 4]
    
    # Index longer than data should raise
    with pytest.raises(ValueError):
        Series(data, index=['a', 'b', 'c', 'd', 'e'])
    
    # Index shorter than data should raise
    with pytest.raises(ValueError):
        Series(data, index=['a', 'b', 'c'])
    
    # Matching lengths should work
    s = Series(data, index=['a', 'b', 'c', 'd'])
    assert len(s) == 4