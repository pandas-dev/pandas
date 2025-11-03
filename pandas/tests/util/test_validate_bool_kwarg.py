"""
Tests for pandas.util._validators.validate_bool_kwarg

Student: Mallikarjuna
Component: validate_bool_kwarg (lines 228-270 in pandas/util/_validators.py)
"""

import pytest
from pandas.util._validators import validate_bool_kwarg


class TestValidateBoolKwarg:
    """Test suite for validate_bool_kwarg function - Initial tests."""

    def test_valid_true(self):
        """Test that True is accepted as a valid boolean."""
        result = validate_bool_kwarg(True, "test_arg")
        assert result is True

    def test_valid_false(self):
        """Test that False is accepted as a valid boolean."""
        result = validate_bool_kwarg(False, "test_arg")
        assert result is False

    def test_none_allowed_default(self):
        """Test that None is allowed by default."""
        result = validate_bool_kwarg(None, "test_arg")
        assert result is None

    def test_none_disallowed(self):
        """Test that None raises ValueError when none_allowed=False."""
        with pytest.raises(ValueError, match='For argument "test_arg" expected type bool'):
            validate_bool_kwarg(None, "test_arg", none_allowed=False)

    def test_int_disallowed_default(self):
        """Test that integers raise ValueError by default."""
        with pytest.raises(ValueError, match='For argument "test_arg" expected type bool'):
            validate_bool_kwarg(1, "test_arg")

    def test_int_allowed(self):
        """Test that integers are accepted when int_allowed=True."""
        result = validate_bool_kwarg(1, "test_arg", int_allowed=True)
        assert result == 1
        
        result = validate_bool_kwarg(0, "test_arg", int_allowed=True)
        assert result == 0

    def test_string_raises_error(self):
        """Test that strings raise ValueError."""
        with pytest.raises(ValueError, match='For argument "my_param" expected type bool, received type str'):
            validate_bool_kwarg("true", "my_param")

    def test_list_raises_error(self):
        """Test that lists raise ValueError."""
        with pytest.raises(ValueError, match='For argument "test_arg" expected type bool'):
            validate_bool_kwarg([True], "test_arg")

    def test_float_raises_error(self):
        """Test that floats raise ValueError."""
        with pytest.raises(ValueError, match='For argument "test_arg" expected type bool'):
            validate_bool_kwarg(1.0, "test_arg")

    # Improvement tests added after analyzing initial mutation results
    def test_none_and_int_both_allowed(self):
        """Test that None and integers can both be allowed together."""
        result = validate_bool_kwarg(None, "test_arg", none_allowed=True, int_allowed=True)
        assert result is None
        
        result = validate_bool_kwarg(1, "test_arg", none_allowed=True, int_allowed=True)
        assert result == 1

    def test_none_and_int_both_disallowed(self):
        """Test that None and integers are both disallowed when both flags are False."""
        with pytest.raises(ValueError, match='For argument "test_arg" expected type bool'):
            validate_bool_kwarg(None, "test_arg", none_allowed=False, int_allowed=False)
        
        with pytest.raises(ValueError, match='For argument "test_arg" expected type bool'):
            validate_bool_kwarg(1, "test_arg", none_allowed=False, int_allowed=False)

    def test_zero_as_integer(self):
        """Test that zero is treated as an integer when int_allowed=True."""
        result = validate_bool_kwarg(0, "test_arg", int_allowed=True)
        assert result == 0
        assert isinstance(result, int)
