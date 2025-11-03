"""
Tests for pandas.util._validators.validate_endpoints

Student: Sandeep Kumar
Component: validate_endpoints (lines 391-420 in pandas/util/_validators.py)
"""

import pytest
from pandas.util._validators import validate_endpoints


class TestValidateEndpoints:
    """Test suite for validate_endpoints function - Initial tests."""

    def test_closed_none(self):
        """Test with closed=None returns both True."""
        left, right = validate_endpoints(None)
        assert left is True
        assert right is True

    def test_closed_left(self):
        """Test with closed='left' returns left=True, right=False."""
        left, right = validate_endpoints("left")
        assert left is True
        assert right is False

    def test_closed_right(self):
        """Test with closed='right' returns left=False, right=True."""
        left, right = validate_endpoints("right")
        assert left is False
        assert right is True

    def test_invalid_string_raises_error(self):
        """Test that invalid strings raise ValueError."""
        with pytest.raises(ValueError, match="Closed has to be either"):
            validate_endpoints("invalid")

    def test_empty_string_raises_error(self):
        """Test that empty string raises ValueError."""
        with pytest.raises(ValueError, match="Closed has to be either"):
            validate_endpoints("")

    def test_uppercase_raises_error(self):
        """Test that uppercase 'LEFT' raises ValueError (case sensitive)."""
        with pytest.raises(ValueError, match="Closed has to be either"):
            validate_endpoints("LEFT")

    def test_integer_raises_error(self):
        """Test that integers raise ValueError."""
        with pytest.raises(ValueError, match="Closed has to be either"):
            validate_endpoints(1)  # type: ignore[arg-type]

    # Improvement tests added after analyzing initial mutation results
    def test_returns_tuple_type(self):
        """Test that the function returns a tuple of exactly 2 booleans."""
        result = validate_endpoints(None)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], bool)
        assert isinstance(result[1], bool)

    def test_left_and_right_mutually_exclusive(self):
        """Test that when left is True, right is False and vice versa."""
        left_closed, right_closed = validate_endpoints("left")
        assert left_closed is True
        assert right_closed is False
        
        left_closed, right_closed = validate_endpoints("right")
        assert left_closed is False
        assert right_closed is True
