"""
Comprehensive tests for Swiss Table hash maps

Tests all numeric types: int64, uint64, int32, uint32, float64, float32
Tests NaN handling for float types
Tests delete, clear, and iteration operations
"""

import math

import numpy as np
import pytest

from pandas._libs import swisstable

import pandas._testing as tm


class TestSwissInt64Map:
    """Test SwissInt64Map basic functionality"""

    def test_init_empty(self):
        """Test empty table initialization"""
        table = swisstable.SwissInt64Map()
        assert len(table) == 0
        assert table.size == 0
        assert table.capacity == 0

    def test_init_with_size_hint(self):
        """Test initialization with size hint"""
        table = swisstable.SwissInt64Map(size_hint=100)
        assert len(table) == 0
        assert table.capacity >= 100
        # Should be normalized to power of 2
        assert table.capacity & (table.capacity - 1) == 0

    def test_insert_single(self):
        """Test inserting a single key-value pair"""
        table = swisstable.SwissInt64Map()
        ret = table.insert(42, 100)
        assert ret is True  # Newly inserted
        assert len(table) == 1
        assert table.get(42) == 100

    def test_insert_update(self):
        """Test updating existing key"""
        table = swisstable.SwissInt64Map()
        table.insert(42, 100)
        ret = table.insert(42, 200)
        assert ret is False  # Updated existing
        assert len(table) == 1
        assert table.get(42) == 200

    def test_insert_multiple(self):
        """Test inserting multiple key-value pairs"""
        table = swisstable.SwissInt64Map()
        keys = [10, 20, 30, 40, 50]
        vals = [100, 200, 300, 400, 500]

        for k, v in zip(keys, vals, strict=True):
            table.insert(k, v)

        assert len(table) == 5
        for k, v in zip(keys, vals, strict=True):
            assert table.get(k) == v

    def test_get_missing_key(self):
        """Test getting non-existent key returns None"""
        table = swisstable.SwissInt64Map()
        table.insert(42, 100)
        assert table.get(999) is None

    def test_get_with_default(self):
        """Test get with custom default value"""
        table = swisstable.SwissInt64Map()
        assert table.get(999, default=-1) == -1

    def test_getitem(self):
        """Test __getitem__ operator"""
        table = swisstable.SwissInt64Map()
        table.insert(42, 100)
        assert table[42] == 100

    def test_getitem_missing_raises(self):
        """Test __getitem__ raises KeyError for missing key"""
        table = swisstable.SwissInt64Map()
        with pytest.raises(KeyError, match="999"):
            _ = table[999]

    def test_setitem(self):
        """Test __setitem__ operator"""
        table = swisstable.SwissInt64Map()
        table[42] = 100
        assert table[42] == 100
        assert len(table) == 1

    def test_setitem_update(self):
        """Test __setitem__ updates existing"""
        table = swisstable.SwissInt64Map()
        table[42] = 100
        table[42] = 200
        assert table[42] == 200
        assert len(table) == 1

    def test_contains(self):
        """Test __contains__ operator (in)"""
        table = swisstable.SwissInt64Map()
        table[42] = 100
        assert 42 in table
        assert 999 not in table

    def test_delete(self):
        """Test delete operation"""
        table = swisstable.SwissInt64Map()
        table.insert(42, 100)
        assert table.delete(42) is True
        assert 42 not in table
        assert table.delete(42) is False
        assert len(table) == 0

    def test_delitem(self):
        """Test __delitem__ operator"""
        table = swisstable.SwissInt64Map()
        table[42] = 100
        del table[42]
        assert 42 not in table

    def test_delitem_keyerror(self):
        """Test __delitem__ raises KeyError for missing key"""
        table = swisstable.SwissInt64Map()
        with pytest.raises(KeyError):
            del table[42]

    def test_clear(self):
        """Test clear operation"""
        table = swisstable.SwissInt64Map()
        for i in range(100):
            table.insert(i, i * 10)
        assert len(table) == 100
        table.clear()
        assert len(table) == 0
        assert table.capacity > 0  # Capacity preserved

    def test_iteration(self):
        """Test iteration over keys"""
        table = swisstable.SwissInt64Map()
        keys = set(range(100))
        for k in keys:
            table.insert(k, k * 10)
        found_keys = set(table)
        assert found_keys == keys

    def test_items(self):
        """Test items() iteration"""
        table = swisstable.SwissInt64Map()
        expected = {i: i * 10 for i in range(100)}
        for k, v in expected.items():
            table.insert(k, v)
        found = dict(table.items())
        assert found == expected

    def test_resize(self):
        """Test automatic resizing on growth"""
        table = swisstable.SwissInt64Map(size_hint=16)
        initial_capacity = table.capacity

        # Insert many elements to trigger resize
        for i in range(100):
            table[i] = i * 10

        assert len(table) == 100
        assert table.capacity > initial_capacity

        # Verify all values are still correct
        for i in range(100):
            assert table[i] == i * 10

    def test_load_factor(self):
        """Test load factor calculation"""
        table = swisstable.SwissInt64Map(size_hint=16)
        assert table.load_factor == 0.0

        # Insert some elements
        for i in range(10):
            table[i] = i

        expected_lf = 10.0 / table.capacity
        assert abs(table.load_factor - expected_lf) < 0.001

    def test_negative_keys(self):
        """Test negative int64 keys"""
        table = swisstable.SwissInt64Map()
        table[-42] = 100
        table[-999] = 200
        assert table[-42] == 100
        assert table[-999] == 200
        assert len(table) == 2

    def test_large_keys(self):
        """Test large int64 keys"""
        table = swisstable.SwissInt64Map()
        large_key = 2**60
        table[large_key] = 999
        assert table[large_key] == 999

    def test_collision_handling(self):
        """Test handling of hash collisions"""
        table = swisstable.SwissInt64Map()

        # Insert many keys to ensure some collisions
        n = 1000
        for i in range(n):
            table[i] = i * 2

        # Verify all keys are stored correctly
        for i in range(n):
            assert table[i] == i * 2
            assert i in table

        assert len(table) == n

    def test_repr(self):
        """Test __repr__ output"""
        table = swisstable.SwissInt64Map()
        table[42] = 100

        repr_str = repr(table)
        assert "SwissInt64Map" in repr_str

    def test_zero_key(self):
        """Test zero as a key"""
        table = swisstable.SwissInt64Map()
        table[0] = 999
        assert table[0] == 999
        assert 0 in table

    def test_min_max_int64(self):
        """Test minimum and maximum int64 values"""
        table = swisstable.SwissInt64Map()

        min_val = -(2**63)
        max_val = 2**63 - 1

        table[min_val] = 1
        table[max_val] = 2

        assert table[min_val] == 1
        assert table[max_val] == 2

    def test_large_scale(self):
        """Test with large number of elements"""
        table = swisstable.SwissInt64Map(size_hint=10000)
        n = 10000

        # Insert
        for i in range(n):
            table[i] = i * 3

        # Verify
        assert len(table) == n
        for i in range(0, n, 100):  # Spot check every 100th
            assert table[i] == i * 3

    def test_delete_and_iterate(self):
        """Test iteration after deletes"""
        table = swisstable.SwissInt64Map()
        for i in range(10):
            table.insert(i, i * 10)

        # Delete every other element
        for i in range(0, 10, 2):
            table.delete(i)

        remaining = set(table)
        assert remaining == {1, 3, 5, 7, 9}


class TestSwissUInt64Map:
    """Tests for SwissUInt64Map: uint64 -> size_t"""

    def test_basic(self):
        table = swisstable.SwissUInt64Map()
        table.insert(42, 100)
        assert table[42] == 100

    def test_max_uint64(self):
        table = swisstable.SwissUInt64Map()
        max_val = np.iinfo(np.uint64).max
        table.insert(max_val, 999)
        assert table[max_val] == 999

    def test_iteration(self):
        table = swisstable.SwissUInt64Map()
        keys = set(range(100))
        for k in keys:
            table.insert(k, k * 10)
        found_keys = set(table)
        assert found_keys == keys


class TestSwissInt32Map:
    """Tests for SwissInt32Map: int32 -> size_t"""

    def test_basic(self):
        table = swisstable.SwissInt32Map()
        table.insert(42, 100)
        assert table[42] == 100

    def test_range(self):
        table = swisstable.SwissInt32Map()
        min_val = np.iinfo(np.int32).min
        max_val = np.iinfo(np.int32).max
        table.insert(min_val, 1)
        table.insert(max_val, 2)
        assert table[min_val] == 1
        assert table[max_val] == 2

    def test_delete_clear(self):
        table = swisstable.SwissInt32Map()
        table.insert(1, 10)
        table.insert(2, 20)
        assert table.delete(1) is True
        assert 1 not in table
        table.clear()
        assert len(table) == 0


class TestSwissUInt32Map:
    """Tests for SwissUInt32Map: uint32 -> size_t"""

    def test_basic(self):
        table = swisstable.SwissUInt32Map()
        table.insert(42, 100)
        assert table[42] == 100


class TestSwissFloat64Map:
    """Tests for SwissFloat64Map: float64 -> size_t with NaN handling"""

    def test_basic(self):
        table = swisstable.SwissFloat64Map()
        table.insert(3.14, 100)
        assert table[3.14] == 100

    def test_nan_equality(self):
        """All NaN values should be treated as equal"""
        table = swisstable.SwissFloat64Map()
        nan1 = float("nan")
        nan2 = float("nan")
        assert nan1 is not nan2  # Different NaN objects

        table.insert(nan1, 100)
        assert nan1 in table
        assert nan2 in table  # Different NaN should find the same entry
        assert table[nan2] == 100

    def test_nan_update(self):
        """Inserting different NaN values should update the same entry"""
        table = swisstable.SwissFloat64Map()
        nan1 = float("nan")
        nan2 = float("nan")

        assert table.insert(nan1, 100) is True  # New
        assert table.insert(nan2, 200) is False  # Update
        assert len(table) == 1
        assert table[nan1] == 200

    def test_nan_from_numpy(self):
        """Test NaN from numpy"""
        table = swisstable.SwissFloat64Map()
        table.insert(np.nan, 100)
        assert np.nan in table
        assert float("nan") in table

    def test_zero_equality(self):
        """+0.0 and -0.0 should be treated as equal"""
        table = swisstable.SwissFloat64Map()
        table.insert(0.0, 100)
        assert 0.0 in table
        assert -0.0 in table
        assert table[-0.0] == 100
        assert len(table) == 1

    def test_inf(self):
        """Test infinity values"""
        table = swisstable.SwissFloat64Map()
        table.insert(float("inf"), 100)
        table.insert(float("-inf"), 200)
        assert table[float("inf")] == 100
        assert table[float("-inf")] == 200
        assert len(table) == 2

    def test_iteration_with_nan(self):
        """Iteration should include NaN values"""
        table = swisstable.SwissFloat64Map()
        table.insert(1.0, 10)
        table.insert(float("nan"), 20)
        table.insert(2.0, 30)

        keys = list(table)
        assert len(keys) == 3
        # Check that exactly one key is NaN
        nan_count = sum(1 for k in keys if math.isnan(k))
        assert nan_count == 1

    def test_delete_nan(self):
        """Test deleting NaN value"""
        table = swisstable.SwissFloat64Map()
        table.insert(float("nan"), 100)
        assert table.delete(float("nan")) is True
        assert float("nan") not in table


class TestSwissFloat32Map:
    """Tests for SwissFloat32Map: float32 -> size_t with NaN handling"""

    def test_basic(self):
        table = swisstable.SwissFloat32Map()
        table.insert(3.14, 100)
        assert table.get(3.14) is not None

    def test_nan_equality(self):
        table = swisstable.SwissFloat32Map()
        nan = float("nan")
        table.insert(nan, 100)
        assert nan in table

    def test_zero_equality(self):
        table = swisstable.SwissFloat32Map()
        table.insert(0.0, 100)
        assert -0.0 in table


class TestCollisions:
    """Test behavior under hash collisions"""

    def test_many_elements_same_hash_bucket(self):
        """Insert many elements that may collide"""
        table = swisstable.SwissInt64Map()
        # Insert 1000 sequential integers
        for i in range(1000):
            table.insert(i, i * 10)

        # Verify all are retrievable
        for i in range(1000):
            assert table[i] == i * 10

    def test_delete_reinsert(self):
        """Test delete followed by reinsert"""
        table = swisstable.SwissInt64Map()
        table.insert(42, 100)
        table.delete(42)
        table.insert(42, 200)
        assert table[42] == 200


class TestEdgeCases:
    """Edge case tests"""

    def test_empty_table(self):
        table = swisstable.SwissInt64Map()
        assert len(table) == 0
        assert 42 not in table
        assert table.get(42) is None
        assert list(table) == []

    def test_single_element(self):
        table = swisstable.SwissInt64Map()
        table.insert(42, 100)
        assert len(table) == 1
        assert list(table) == [42]
        assert list(table.items()) == [(42, 100)]

    def test_all_duplicates(self):
        """Inserting the same key repeatedly"""
        table = swisstable.SwissInt64Map()
        for i in range(100):
            table.insert(42, i)
        assert len(table) == 1
        assert table[42] == 99


class TestProperties:
    """Test table properties"""

    def test_size_property(self):
        table = swisstable.SwissInt64Map()
        assert table.size == 0
        table.insert(1, 10)
        assert table.size == 1
        table.insert(2, 20)
        assert table.size == 2
        table.delete(1)
        assert table.size == 1

    def test_capacity_property(self):
        table = swisstable.SwissInt64Map()
        assert table.capacity == 0
        table.insert(1, 10)
        assert table.capacity >= 16  # Minimum capacity is 16

    def test_load_factor_bounds(self):
        """Test load factor stays within reasonable bounds"""
        table = swisstable.SwissInt64Map(size_hint=16)

        for i in range(100):
            table[i] = i
            # Load factor should stay reasonable
            assert 0.0 <= table.load_factor <= 0.90


# =============================================================================
# Integration Tests - map_locations, lookup, ismember
# =============================================================================


class TestMapLocations:
    """Test map_locations() method for array-based indexing"""

    def test_int64_basic(self):
        """Test basic map_locations functionality"""
        table = swisstable.SwissInt64Map()
        values = np.array([10, 20, 30, 40, 50], dtype=np.int64)
        table.map_locations(values)

        assert len(table) == 5
        assert table[10] == 0
        assert table[20] == 1
        assert table[30] == 2
        assert table[40] == 3
        assert table[50] == 4

    def test_int64_duplicates(self):
        """Test map_locations with duplicate values (last occurrence wins)"""
        table = swisstable.SwissInt64Map()
        values = np.array([10, 20, 10, 30, 20], dtype=np.int64)
        table.map_locations(values)

        assert len(table) == 3
        assert table[10] == 2  # Last occurrence at index 2
        assert table[20] == 4  # Last occurrence at index 4
        assert table[30] == 3

    def test_int64_empty(self):
        """Test map_locations with empty array"""
        table = swisstable.SwissInt64Map()
        values = np.array([], dtype=np.int64)
        table.map_locations(values)
        assert len(table) == 0

    def test_float64_with_nan(self):
        """Test map_locations with NaN values"""
        table = swisstable.SwissFloat64Map()
        values = np.array([1.0, np.nan, 2.0, np.nan, 3.0], dtype=np.float64)
        table.map_locations(values)

        assert len(table) == 4  # Two NaNs count as one
        assert table[1.0] == 0
        assert table[np.nan] == 3  # Last NaN occurrence
        assert table[2.0] == 2
        assert table[3.0] == 4

    def test_all_types(self):
        """Test map_locations for all integer types"""
        # Int64
        t64 = swisstable.SwissInt64Map()
        t64.map_locations(np.array([1, 2, 3], dtype=np.int64))
        assert t64[2] == 1

        # UInt64
        tu64 = swisstable.SwissUInt64Map()
        tu64.map_locations(np.array([1, 2, 3], dtype=np.uint64))
        assert tu64[2] == 1

        # Int32
        t32 = swisstable.SwissInt32Map()
        t32.map_locations(np.array([1, 2, 3], dtype=np.int32))
        assert t32[2] == 1

        # UInt32
        tu32 = swisstable.SwissUInt32Map()
        tu32.map_locations(np.array([1, 2, 3], dtype=np.uint32))
        assert tu32[2] == 1


class TestLookup:
    """Test lookup() method for array-based lookups"""

    def test_int64_basic(self):
        """Test basic lookup functionality"""
        table = swisstable.SwissInt64Map()
        table.map_locations(np.array([10, 20, 30], dtype=np.int64))

        lookup_values = np.array([20, 10, 30], dtype=np.int64)
        result = table.lookup(lookup_values)

        tm.assert_numpy_array_equal(result, np.array([1, 0, 2], dtype=np.intp))

    def test_int64_misses(self):
        """Test lookup with missing values returns -1"""
        table = swisstable.SwissInt64Map()
        table.map_locations(np.array([10, 20, 30], dtype=np.int64))

        lookup_values = np.array([10, 99, 20, 88], dtype=np.int64)
        result = table.lookup(lookup_values)

        tm.assert_numpy_array_equal(result, np.array([0, -1, 1, -1], dtype=np.intp))

    def test_int64_all_misses(self):
        """Test lookup where all values are missing"""
        table = swisstable.SwissInt64Map()
        table.map_locations(np.array([10, 20, 30], dtype=np.int64))

        lookup_values = np.array([100, 200, 300], dtype=np.int64)
        result = table.lookup(lookup_values)

        tm.assert_numpy_array_equal(result, np.array([-1, -1, -1], dtype=np.intp))

    def test_int64_empty_table(self):
        """Test lookup on empty table"""
        table = swisstable.SwissInt64Map()
        lookup_values = np.array([1, 2, 3], dtype=np.int64)
        result = table.lookup(lookup_values)

        tm.assert_numpy_array_equal(result, np.array([-1, -1, -1], dtype=np.intp))

    def test_float64_with_nan(self):
        """Test lookup with NaN values"""
        table = swisstable.SwissFloat64Map()
        table.map_locations(np.array([1.0, np.nan, 2.0], dtype=np.float64))

        lookup_values = np.array([np.nan, 1.0, 3.0], dtype=np.float64)
        result = table.lookup(lookup_values)

        tm.assert_numpy_array_equal(result, np.array([1, 0, -1], dtype=np.intp))

    def test_all_types(self):
        """Test lookup for all integer types"""
        for MapClass, dtype in [
            (swisstable.SwissInt64Map, np.int64),
            (swisstable.SwissUInt64Map, np.uint64),
            (swisstable.SwissInt32Map, np.int32),
            (swisstable.SwissUInt32Map, np.uint32),
            (swisstable.SwissInt16Map, np.int16),
            (swisstable.SwissUInt16Map, np.uint16),
            (swisstable.SwissInt8Map, np.int8),
            (swisstable.SwissUInt8Map, np.uint8),
        ]:
            table = MapClass()
            table.map_locations(np.array([1, 2, 3], dtype=dtype))
            result = table.lookup(np.array([2, 99], dtype=dtype))
            tm.assert_numpy_array_equal(result, np.array([1, -1], dtype=np.intp))


class TestIsmember:
    """Test ismember_*() standalone functions for membership testing"""

    def test_int64_basic(self):
        """Test basic ismember functionality"""
        arr = np.array([1, 2, 3, 4, 5], dtype=np.int64)
        values = np.array([2, 4], dtype=np.int64)

        result = swisstable.ismember_int64(arr, values)

        expected = np.array([False, True, False, True, False])
        tm.assert_numpy_array_equal(result, expected)

    def test_int64_no_matches(self):
        """Test ismember with no matches"""
        arr = np.array([1, 2, 3], dtype=np.int64)
        values = np.array([10, 20, 30], dtype=np.int64)

        result = swisstable.ismember_int64(arr, values)

        tm.assert_numpy_array_equal(result, np.array([False, False, False]))

    def test_int64_all_matches(self):
        """Test ismember with all matches"""
        arr = np.array([1, 2, 3], dtype=np.int64)
        values = np.array([1, 2, 3, 4, 5], dtype=np.int64)

        result = swisstable.ismember_int64(arr, values)

        tm.assert_numpy_array_equal(result, np.array([True, True, True]))

    def test_int64_empty_values(self):
        """Test ismember with empty values array"""
        arr = np.array([1, 2, 3], dtype=np.int64)
        values = np.array([], dtype=np.int64)

        result = swisstable.ismember_int64(arr, values)

        tm.assert_numpy_array_equal(result, np.array([False, False, False]))

    def test_int64_empty_arr(self):
        """Test ismember with empty arr array"""
        arr = np.array([], dtype=np.int64)
        values = np.array([1, 2, 3], dtype=np.int64)

        result = swisstable.ismember_int64(arr, values)

        assert len(result) == 0

    def test_int64_duplicates_in_values(self):
        """Test ismember with duplicates in values"""
        arr = np.array([1, 2, 3, 4], dtype=np.int64)
        values = np.array([2, 2, 2, 4, 4], dtype=np.int64)

        result = swisstable.ismember_int64(arr, values)

        tm.assert_numpy_array_equal(result, np.array([False, True, False, True]))

    def test_float64_basic(self):
        """Test ismember for float64"""
        arr = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        values = np.array([2.0, 4.0], dtype=np.float64)

        result = swisstable.ismember_float64(arr, values)

        tm.assert_numpy_array_equal(result, np.array([False, True, False]))

    def test_float64_with_nan(self):
        """Test ismember with NaN values"""
        arr = np.array([1.0, np.nan, 2.0, np.nan], dtype=np.float64)
        values = np.array([np.nan, 2.0], dtype=np.float64)

        result = swisstable.ismember_float64(arr, values)

        tm.assert_numpy_array_equal(result, np.array([False, True, True, True]))

    def test_all_types(self):
        """Test ismember for all numeric types"""
        test_cases = [
            (swisstable.ismember_int64, np.int64),
            (swisstable.ismember_uint64, np.uint64),
            (swisstable.ismember_int32, np.int32),
            (swisstable.ismember_uint32, np.uint32),
            (swisstable.ismember_int16, np.int16),
            (swisstable.ismember_uint16, np.uint16),
            (swisstable.ismember_int8, np.int8),
            (swisstable.ismember_uint8, np.uint8),
            (swisstable.ismember_float64, np.float64),
            (swisstable.ismember_float32, np.float32),
        ]

        for ismember_fn, dtype in test_cases:
            arr = np.array([1, 2, 3, 4, 5], dtype=dtype)
            values = np.array([2, 4], dtype=dtype)
            result = ismember_fn(arr, values)
            expected = np.array([False, True, False, True, False])
            tm.assert_numpy_array_equal(result, expected, err_msg=f"Failed for {dtype}")

    def test_large_scale(self):
        """Test ismember with large arrays (performance-oriented)"""
        rng = np.random.default_rng(42)
        arr = rng.integers(0, 100000, size=100000, dtype=np.int64)
        values = rng.integers(0, 100000, size=1000, dtype=np.int64)

        result = swisstable.ismember_int64(arr, values)

        # Verify result type and shape
        assert result.dtype == np.bool_
        assert result.shape == arr.shape

        # Verify correctness with Python set (slow but correct)
        value_set = set(values)
        expected = np.array([x in value_set for x in arr])
        tm.assert_numpy_array_equal(result, expected)


class TestUnique:
    """Test unique() method for getting unique values"""

    def test_int64_basic(self):
        """Test basic unique functionality"""
        table = swisstable.SwissInt64Map()
        values = np.array([3, 1, 2, 1, 3, 4, 2], dtype=np.int64)
        uniques = table.unique(values)

        # Should preserve order of first occurrence
        tm.assert_numpy_array_equal(uniques, np.array([3, 1, 2, 4], dtype=np.intp))

    def test_int64_with_inverse(self):
        """Test unique with return_inverse"""
        table = swisstable.SwissInt64Map()
        values = np.array([3, 1, 2, 1, 3, 4, 2], dtype=np.int64)
        uniques, labels = table.unique(values, return_inverse=True)

        tm.assert_numpy_array_equal(uniques, np.array([3, 1, 2, 4], dtype=np.intp))
        expected_labels = np.array([0, 1, 2, 1, 0, 3, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(labels, expected_labels)

    def test_int64_all_unique(self):
        """Test unique with no duplicates"""
        table = swisstable.SwissInt64Map()
        values = np.array([1, 2, 3, 4, 5], dtype=np.int64)
        uniques = table.unique(values)

        tm.assert_numpy_array_equal(uniques, values)

    def test_int64_all_same(self):
        """Test unique with all same values"""
        table = swisstable.SwissInt64Map()
        values = np.array([42, 42, 42, 42], dtype=np.int64)
        uniques = table.unique(values)

        tm.assert_numpy_array_equal(uniques, np.array([42], dtype=np.intp))

    def test_int64_empty(self):
        """Test unique with empty array"""
        table = swisstable.SwissInt64Map()
        values = np.array([], dtype=np.int64)
        uniques = table.unique(values)

        assert len(uniques) == 0

    def test_float64_with_nan(self):
        """Test unique with NaN values"""
        table = swisstable.SwissFloat64Map()
        values = np.array([1.0, np.nan, 2.0, np.nan, 1.0], dtype=np.float64)
        uniques = table.unique(values)

        # NaN appears once, in first position
        assert len(uniques) == 3
        assert uniques[0] == 1.0
        assert np.isnan(uniques[1])
        assert uniques[2] == 2.0


class TestFactorize:
    """Test factorize() method for categorical encoding"""

    def test_int64_basic(self):
        """Test basic factorize functionality"""
        table = swisstable.SwissInt64Map()
        values = np.array([10, 20, 10, 30, 20], dtype=np.int64)
        uniques, labels = table.factorize(values)

        tm.assert_numpy_array_equal(uniques, np.array([10, 20, 30], dtype=np.intp))
        tm.assert_numpy_array_equal(labels, np.array([0, 1, 0, 2, 1], dtype=np.intp))

    def test_int64_single_value(self):
        """Test factorize with single unique value"""
        table = swisstable.SwissInt64Map()
        values = np.array([42, 42, 42], dtype=np.int64)
        uniques, labels = table.factorize(values)

        tm.assert_numpy_array_equal(uniques, np.array([42], dtype=np.intp))
        tm.assert_numpy_array_equal(labels, np.array([0, 0, 0], dtype=np.intp))

    def test_float64_with_nan_ignore(self):
        """Test factorize with NaN values (ignore_na=True)"""
        table = swisstable.SwissFloat64Map()
        values = np.array([1.0, np.nan, 2.0, np.nan, 1.0], dtype=np.float64)
        uniques, labels = table.factorize(values, na_sentinel=-1)

        tm.assert_numpy_array_equal(uniques, np.array([1.0, 2.0], dtype=np.float64))
        tm.assert_numpy_array_equal(labels, np.array([0, -1, 1, -1, 0], dtype=np.intp))


class TestValueCount:
    """Test value_count_*() standalone functions"""

    def test_int64_basic(self):
        """Test basic value_count functionality"""
        values = np.array([1, 2, 1, 3, 2, 1], dtype=np.int64)
        keys, counts = swisstable.value_count_int64(values)

        # Order should be by first occurrence
        tm.assert_numpy_array_equal(keys, np.array([1, 2, 3], dtype=np.intp))
        tm.assert_numpy_array_equal(counts, np.array([3, 2, 1], dtype=np.intp))

    def test_int64_all_same(self):
        """Test value_count with all same values"""
        values = np.array([42, 42, 42, 42, 42], dtype=np.int64)
        keys, counts = swisstable.value_count_int64(values)

        tm.assert_numpy_array_equal(keys, np.array([42], dtype=np.intp))
        tm.assert_numpy_array_equal(counts, np.array([5], dtype=np.intp))

    def test_int64_all_unique(self):
        """Test value_count with all unique values"""
        values = np.array([1, 2, 3, 4, 5], dtype=np.int64)
        keys, counts = swisstable.value_count_int64(values)

        tm.assert_numpy_array_equal(keys, np.array([1, 2, 3, 4, 5], dtype=np.intp))
        tm.assert_numpy_array_equal(counts, np.array([1, 1, 1, 1, 1], dtype=np.intp))

    def test_float64_with_nan_dropna_true(self):
        """Test value_count with NaN values (dropna=True)"""
        values = np.array([1.0, np.nan, 1.0, np.nan, 2.0], dtype=np.float64)
        keys, counts = swisstable.value_count_float64(values, dropna=True)

        tm.assert_numpy_array_equal(keys, np.array([1.0, 2.0], dtype=np.float64))
        tm.assert_numpy_array_equal(counts, np.array([2, 1], dtype=np.intp))

    def test_float64_with_nan_dropna_false(self):
        """Test value_count with NaN values (dropna=False)"""
        values = np.array([1.0, np.nan, 1.0, np.nan, 2.0], dtype=np.float64)
        keys, counts = swisstable.value_count_float64(values, dropna=False)

        # NaN is appended at the end
        assert len(keys) == 3
        assert keys[0] == 1.0
        assert keys[1] == 2.0
        assert np.isnan(keys[2])
        tm.assert_numpy_array_equal(counts, np.array([2, 1, 2], dtype=np.intp))


class TestDuplicated:
    """Test duplicated_*() standalone functions"""

    def test_int64_keep_first(self):
        """Test duplicated with keep='first'"""
        values = np.array([1, 2, 1, 3, 2], dtype=np.int64)
        result = swisstable.duplicated_int64(values, keep="first")

        tm.assert_numpy_array_equal(result, np.array([False, False, True, False, True]))

    def test_int64_keep_last(self):
        """Test duplicated with keep='last'"""
        values = np.array([1, 2, 1, 3, 2], dtype=np.int64)
        result = swisstable.duplicated_int64(values, keep="last")

        tm.assert_numpy_array_equal(result, np.array([True, True, False, False, False]))

    def test_int64_keep_false(self):
        """Test duplicated with keep=False (mark all duplicates)"""
        values = np.array([1, 2, 1, 3, 2], dtype=np.int64)
        result = swisstable.duplicated_int64(values, keep=False)

        tm.assert_numpy_array_equal(result, np.array([True, True, True, False, True]))

    def test_int64_no_duplicates(self):
        """Test duplicated with no duplicates"""
        values = np.array([1, 2, 3, 4, 5], dtype=np.int64)
        result = swisstable.duplicated_int64(values, keep="first")

        expected = np.array([False, False, False, False, False])
        tm.assert_numpy_array_equal(result, expected)

    def test_int64_all_duplicates(self):
        """Test duplicated with all same values"""
        values = np.array([42, 42, 42], dtype=np.int64)
        result = swisstable.duplicated_int64(values, keep="first")

        tm.assert_numpy_array_equal(result, np.array([False, True, True]))

    def test_float64_with_nan_keep_first(self):
        """Test duplicated with NaN values (keep='first')"""
        values = np.array([1.0, np.nan, 1.0, np.nan], dtype=np.float64)
        result = swisstable.duplicated_float64(values, keep="first")

        tm.assert_numpy_array_equal(result, np.array([False, False, True, True]))

    def test_float64_with_nan_keep_false(self):
        """Test duplicated with NaN values (keep=False)"""
        values = np.array([1.0, np.nan, 1.0, np.nan], dtype=np.float64)
        result = swisstable.duplicated_float64(values, keep=False)

        tm.assert_numpy_array_equal(result, np.array([True, True, True, True]))

    def test_invalid_keep_value(self):
        """Test duplicated with invalid keep value"""
        values = np.array([1, 2, 3], dtype=np.int64)
        with pytest.raises(ValueError, match="keep must be either"):
            swisstable.duplicated_int64(values, keep="invalid")
