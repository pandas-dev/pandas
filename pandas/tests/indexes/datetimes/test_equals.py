import pandas as pd
import pytest


class TestDatetimeIndexEqualsCheckDtype:
    """Tests for DatetimeIndex.equals() with check_dtype parameter."""

    def test_equals_datetime_different_units_check_dtype_true(self):
        """Default behavior: different datetime units should not be equal."""
        idx1 = pd.date_range("2020-01-01", periods=10).as_unit("ns")
        idx2 = pd.date_range("2020-01-01", periods=10).as_unit("s")

        # Default check_dtype=True should return False
        assert idx1.equals(idx2, check_dtype=True) is False
        assert idx1.equals(idx2) is False  # default behavior

    def test_equals_datetime_different_units_check_dtype_false(self):
        """With check_dtype=False, same values but different units should be equal."""
        idx1 = pd.date_range("2020-01-01", periods=10).as_unit("ns")
        idx2 = pd.date_range("2020-01-01", periods=10).as_unit("s")

        # check_dtype=False should return True (same values)
        assert idx1.equals(idx2, check_dtype=False) is True

    def test_equals_datetime_ns_vs_ms(self):
        """Test nanosecond vs millisecond precision."""
        idx1 = pd.date_range("2020-01-01", periods=5).as_unit("ns")
        idx2 = pd.date_range("2020-01-01", periods=5).as_unit("ms")

        assert idx1.equals(idx2, check_dtype=True) is False
        assert idx1.equals(idx2, check_dtype=False) is True

    def test_equals_datetime_ns_vs_us(self):
        """Test nanosecond vs microsecond precision."""
        idx1 = pd.date_range("2020-01-01", periods=5).as_unit("ns")
        idx2 = pd.date_range("2020-01-01", periods=5).as_unit("us")

        assert idx1.equals(idx2, check_dtype=True) is False
        assert idx1.equals(idx2, check_dtype=False) is True

    def test_equals_datetime_different_values_check_dtype_false(self):
        """Even with check_dtype=False, different values should not be equal."""
        idx1 = pd.date_range("2020-01-01", periods=5)
        idx2 = pd.date_range("2020-01-02", periods=5)

        assert idx1.equals(idx2, check_dtype=False) is False

    def test_equals_datetime_different_tz_check_dtype_false(self):
        """check_dtype=False should still respect timezone differences."""
        idx1 = pd.date_range("2020-01-01", periods=3, tz="UTC")
        idx2 = pd.date_range("2020-01-01", periods=3, tz="US/Eastern")

        # Different timezones have different underlying values
        assert idx1.equals(idx2, check_dtype=False) is False

    def test_equals_datetime_same_tz_different_units(self):
        """Test timezone-aware indexes with different units."""
        idx1 = pd.date_range("2020-01-01", periods=5, tz="UTC").as_unit("ns")
        idx2 = pd.date_range("2020-01-01", periods=5, tz="UTC").as_unit("s")

        assert idx1.equals(idx2, check_dtype=True) is False
        assert idx1.equals(idx2, check_dtype=False) is True

    def test_equals_datetime_with_freq_different_units(self):
        """Test datetime indexes with frequency and different units."""
        idx1 = pd.date_range("2020-01-01", periods=100, freq="D").as_unit("ns")
        idx2 = pd.date_range("2020-01-01", periods=100, freq="D").as_unit("s")

        assert idx1.equals(idx2, check_dtype=True) is False
        assert idx1.equals(idx2, check_dtype=False) is True

    def test_equals_datetime_nat_values(self):
        """Test equals with NaT (Not-a-Time) values."""
        idx1 = pd.DatetimeIndex(["2020-01-01", "NaT", "2020-01-03"]).as_unit("ns")
        idx2 = pd.DatetimeIndex(["2020-01-01", "NaT", "2020-01-03"]).as_unit("s")

        assert idx1.equals(idx2, check_dtype=True) is False
        assert idx1.equals(idx2, check_dtype=False) is True