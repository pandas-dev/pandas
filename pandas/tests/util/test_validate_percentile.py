"""
Tests for pandas.util._validators.validate_percentile

Student: Nithikesh Bobbili
Component: validate_percentile (lines 339-368 in pandas/util/_validators.py)
"""

import numpy as np
import pytest
from pandas.util._validators import validate_percentile


class TestValidatePercentile:
    """Test suite for validate_percentile function - Initial tests."""

    def test_valid_single_percentile(self):
        """Test that a single valid percentile is accepted."""
        result = validate_percentile(0.5)
        np.testing.assert_array_equal(result, np.array(0.5))

    def test_valid_zero(self):
        """Test that 0 is accepted as a valid percentile."""
        result = validate_percentile(0.0)
        np.testing.assert_array_equal(result, np.array(0.0))

    def test_valid_one(self):
        """Test that 1 is accepted as a valid percentile."""
        result = validate_percentile(1.0)
        np.testing.assert_array_equal(result, np.array(1.0))

    def test_valid_list(self):
        """Test that a list of valid percentiles is accepted."""
        result = validate_percentile([0.25, 0.5, 0.75])
        np.testing.assert_array_equal(result, np.array([0.25, 0.5, 0.75]))

    def test_valid_tuple(self):
        """Test that a tuple of valid percentiles is accepted."""
        result = validate_percentile((0.1, 0.9))
        np.testing.assert_array_equal(result, np.array([0.1, 0.9]))

    def test_valid_array(self):
        """Test that a numpy array of valid percentiles is accepted."""
        result = validate_percentile(np.array([0.1, 0.5, 0.9]))
        np.testing.assert_array_equal(result, np.array([0.1, 0.5, 0.9]))

    def test_valid_boundary_values(self):
        """Test boundary values [0, 1] are accepted."""
        result = validate_percentile([0, 1])
        np.testing.assert_array_equal(result, np.array([0, 1]))

    def test_invalid_below_zero(self):
        """Test that percentiles below 0 raise ValueError."""
        with pytest.raises(ValueError, match="percentiles should all be in the interval"):
            validate_percentile(-0.1)

    def test_invalid_above_one(self):
        """Test that percentiles above 1 raise ValueError."""
        with pytest.raises(ValueError, match="percentiles should all be in the interval"):
            validate_percentile(1.5)

    def test_invalid_in_list(self):
        """Test that invalid percentiles in a list raise ValueError."""
        with pytest.raises(ValueError, match="percentiles should all be in the interval"):
            validate_percentile([0.5, 1.5])

    def test_mixed_valid_invalid(self):
        """Test that mixed valid/invalid percentiles raise ValueError."""
        with pytest.raises(ValueError, match="percentiles should all be in the interval"):
            validate_percentile([0.25, -0.1, 0.75])

    # Improvement tests added after analyzing initial mutation results
    def test_returns_ndarray_type(self):
        """Test that the function always returns a numpy ndarray."""
        result = validate_percentile(0.5)
        assert isinstance(result, np.ndarray)
        
        result = validate_percentile([0.25, 0.75])
        assert isinstance(result, np.ndarray)

    def test_edge_case_just_above_one(self):
        """Test that 1.0000001 raises ValueError."""
        with pytest.raises(ValueError, match="percentiles should all be in the interval"):
            validate_percentile(1.0000001)

    def test_edge_case_just_below_zero(self):
        """Test that -0.0000001 raises ValueError."""
        with pytest.raises(ValueError, match="percentiles should all be in the interval"):
            validate_percentile(-0.0000001)
