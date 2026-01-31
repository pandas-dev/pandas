import numpy as np
import pytest

import pandas as pd
from pandas import (
    Series,
    Timestamp,
    date_range,
)
import pandas._testing as tm
from pandas.api.types import is_scalar


class TestSeriesSearchSorted:
    def test_searchsorted(self):
        ser = Series([1, 2, 3])

        result = ser.searchsorted(1, side="left")
        assert is_scalar(result)
        assert result == 0

        result = ser.searchsorted(1, side="right")
        assert is_scalar(result)
        assert result == 1

    def test_searchsorted_numeric_dtypes_scalar(self):
        ser = Series([1, 2, 90, 1000, 3e9])
        res = ser.searchsorted(30)
        assert is_scalar(res)
        assert res == 2

        res = ser.searchsorted([30])
        exp = np.array([2], dtype=np.intp)
        tm.assert_numpy_array_equal(res, exp)

    def test_searchsorted_numeric_dtypes_vector(self):
        ser = Series([1, 2, 90, 1000, 3e9])
        res = ser.searchsorted([91, 2e6])
        exp = np.array([3, 4], dtype=np.intp)
        tm.assert_numpy_array_equal(res, exp)

    def test_searchsorted_datetime64_scalar(self):
        ser = Series(date_range("20120101", periods=10, freq="2D"))
        val = Timestamp("20120102")
        res = ser.searchsorted(val)
        assert is_scalar(res)
        assert res == 1

    def test_searchsorted_datetime64_scalar_mixed_timezones(self):
        # GH 30086
        ser = Series(date_range("20120101", periods=10, freq="2D", tz="UTC"))
        val = Timestamp("20120102", tz="America/New_York")
        res = ser.searchsorted(val)
        assert is_scalar(res)
        assert res == 1

    def test_searchsorted_datetime64_list(self):
        ser = Series(date_range("20120101", periods=10, freq="2D"))
        vals = [Timestamp("20120102"), Timestamp("20120104")]
        res = ser.searchsorted(vals)
        exp = np.array([1, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(res, exp)

    def test_searchsorted_sorter(self):
        # GH8490
        ser = Series([3, 1, 2])
        res = ser.searchsorted([0, 3], sorter=np.argsort(ser))
        exp = np.array([0, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(res, exp)

    def test_searchsorted_dataframe_fail(self):
        # GH#49620
        ser = Series([1, 2, 3, 4, 5])
        vals = pd.DataFrame([[1, 2], [3, 4]])
        msg = "Value must be 1-D array-like or scalar, DataFrame is not supported"
        with pytest.raises(ValueError, match=msg):
            ser.searchsorted(vals)


class TestSearchsortedWithNullableIntegers:
    """Comprehensive tests for searchsorted with nullable integer arrays."""
    def test_series_searchsorted_nullable_int_common_dtype(self):
        base = 2**53
        s = pd.Series(np.array([base, base+1, base+2, base+3], dtype=np.int64))
        value = pd.array([base+1, pd.NA], dtype="Int64")

        result = s.searchsorted(value)
        assert result.tolist() == [1, 4]

    def test_searchsorted_int64_with_single_na(self):
        """Test searchsorted with Int64 array containing one NA."""
        base = 2**53
        s = pd.Series(np.array([base, base+1, base+2, base+3], dtype=np.int64))
        value = pd.array([base+1, pd.NA], dtype="Int64")
        
        result = s.searchsorted(value, side="left")
        expected = np.array([1, 4], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_int64_with_single_na_right(self):
        """Test searchsorted with side='right'."""
        base = 2**53
        s = pd.Series(np.array([base, base+1, base+2, base+3], dtype=np.int64))
        value = pd.array([base+1, pd.NA], dtype="Int64")
        
        result = s.searchsorted(value, side="right")
        expected = np.array([2, 4], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_all_na_values(self):
        """Test searchsorted when all values are NA."""
        s = pd.Series([1, 2, 3, 4, 5], dtype="int64")
        value = pd.array([pd.NA, pd.NA, pd.NA], dtype="Int64")
        
        result = s.searchsorted(value)
        expected = np.array([5, 5, 5], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_no_na_values(self):
        """Test searchsorted with Int64 array but no NA values."""
        s = pd.Series([1, 2, 3, 4, 5], dtype="int64")
        value = pd.array([2, 4], dtype="Int64")
        
        result = s.searchsorted(value)
        expected = np.array([1, 3], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_mixed_valid_and_na(self):
        """Test searchsorted with mixture of valid values and NAs."""
        s = pd.Series([10, 20, 30, 40, 50], dtype="int64")
        value = pd.array([15, pd.NA, 25, pd.NA, 35], dtype="Int64")
        
        result = s.searchsorted(value)
        expected = np.array([1, 5, 2, 5, 3], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_na_at_start(self):
        """Test searchsorted with NA at the beginning."""
        s = pd.Series([1, 2, 3, 4, 5], dtype="int64")
        value = pd.array([pd.NA, 2, 3], dtype="Int64")
        
        result = s.searchsorted(value)
        expected = np.array([5, 1, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_na_at_end(self):
        """Test searchsorted with NA at the end."""
        s = pd.Series([1, 2, 3, 4, 5], dtype="int64")
        value = pd.array([2, 3, pd.NA], dtype="Int64")
        
        result = s.searchsorted(value)
        expected = np.array([1, 2, 5], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_empty_array_with_na(self):
        """Test searchsorted on empty array with NA values."""
        s = pd.Series([], dtype="int64")
        value = pd.array([pd.NA, 1, pd.NA], dtype="Int64")
        
        result = s.searchsorted(value)
        expected = np.array([0, 0, 0], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_single_element_with_na(self):
        """Test searchsorted on single element array."""
        s = pd.Series([5], dtype="int64")
        value = pd.array([pd.NA, 3, 5, 7], dtype="Int64")
        
        result = s.searchsorted(value)
        expected = np.array([1, 0, 0, 1], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_int32_with_na(self):
        """Test searchsorted with Int32 dtype."""
        s = pd.Series([1, 2, 3, 4, 5], dtype="int32")
        value = pd.array([2, pd.NA, 4], dtype="Int32")
        
        result = s.searchsorted(value)
        expected = np.array([1, 5, 3], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_int8_with_na(self):
        """Test searchsorted with Int8 dtype."""
        s = pd.Series([1, 2, 3, 4, 5], dtype="int8")
        value = pd.array([pd.NA, 3], dtype="Int8")
        
        result = s.searchsorted(value)
        expected = np.array([5, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_duplicates_with_na(self):
        """Test searchsorted with duplicate values in array."""
        s = pd.Series([1, 2, 2, 3, 3, 3, 4], dtype="int64")
        value = pd.array([2, pd.NA, 3], dtype="Int64")
        
        # side='left'
        result_left = s.searchsorted(value, side="left")
        expected_left = np.array([1, 7, 3], dtype=np.intp)
        tm.assert_numpy_array_equal(result_left, expected_left)
        
        # side='right'
        result_right = s.searchsorted(value, side="right")
        expected_right = np.array([3, 7, 6], dtype=np.intp)
        tm.assert_numpy_array_equal(result_right, expected_right)

    def test_searchsorted_uint64_with_na(self):
        """Test searchsorted with UInt64 dtype."""
        s = pd.Series([1, 2, 3, 4, 5], dtype="uint64")
        value = pd.array([2, pd.NA], dtype="UInt64")
        
        result = s.searchsorted(value)
        expected = np.array([1, 5], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_large_array_with_na(self):
        """Test searchsorted with larger array."""
        s = pd.Series(range(1000), dtype="int64")
        value = pd.array([100, pd.NA, 500, pd.NA, 999], dtype="Int64")
        
        result = s.searchsorted(value)
        expected = np.array([100, 1000, 500, 1000, 999], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_boundary_values_with_na(self):
        """Test searchsorted with values at boundaries."""
        s = pd.Series([10, 20, 30, 40, 50], dtype="int64")
        value = pd.array([0, pd.NA, 10, 50, 60, pd.NA], dtype="Int64")
        
        result = s.searchsorted(value)
        expected = np.array([0, 5, 0, 4, 5, 5], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_negative_values_with_na(self):
        """Test searchsorted with negative values."""
        s = pd.Series([-50, -30, -10, 0, 10, 30, 50], dtype="int64")
        value = pd.array([-40, pd.NA, 0, pd.NA, 20], dtype="Int64")
        
        result = s.searchsorted(value)
        expected = np.array([1, 7, 3, 7, 5], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_preserves_regular_behavior(self):
        """Ensure regular searchsorted still works without NA."""
        s = pd.Series([1, 2, 3, 4, 5], dtype="int64")
        
        # Regular list
        result = s.searchsorted([2, 4])
        expected = np.array([1, 3], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)
        
        # Regular numpy array
        result = s.searchsorted(np.array([2, 4]))
        expected = np.array([1, 3], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)
        
        # Scalar
        result = s.searchsorted(3)
        assert result == 2

    def test_searchsorted_float_array_with_na(self):
        """Test that Float64 arrays with NA also work."""
        s = pd.Series([1.5, 2.5, 3.5, 4.5, 5.5], dtype="float64")
        value = pd.array([2.5, pd.NA, 4.0], dtype="Float64")
        
        result = s.searchsorted(value)
        expected = np.array([1, 5, 3], dtype=np.intp)  # Changed from [1, 5, 2]
        tm.assert_numpy_array_equal(result, expected)
    def test_searchsorted_series_with_int64_array(self):
        """Test on Series (not just arrays)."""
        s = pd.Series([10, 20, 30, 40, 50])
        value = pd.array([25, pd.NA, 35], dtype="Int64")
        
        result = s.searchsorted(value)
        expected = np.array([2, 5, 3], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)

    def test_searchsorted_index_with_na(self):
        """Test searchsorted on Index."""
        idx = pd.Index([1, 2, 3, 4, 5], dtype="int64")
        value = pd.array([2, pd.NA, 4], dtype="Int64")
        
        result = idx.searchsorted(value)
        expected = np.array([1, 5, 3], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)
