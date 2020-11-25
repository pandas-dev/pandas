from contextlib import contextmanager
import tracemalloc

import numpy as np
import pytest

from pandas._libs import hashtable as ht

import pandas._testing as tm


@contextmanager
def activated_tracemalloc():
    tracemalloc.start()
    try:
        yield
    finally:
        tracemalloc.stop()


def get_allocated_khash_memory():
    snapshot = tracemalloc.take_snapshot()
    snapshot = snapshot.filter_traces(
        (tracemalloc.DomainFilter(True, ht.get_hashtable_trace_domain()),)
    )
    return sum(map(lambda x: x.size, snapshot.traces))


@pytest.mark.parametrize(
    "table_type, dtype",
    [
        (ht.PyObjectHashTable, np.object_),
        (ht.Int64HashTable, np.int64),
        (ht.UInt64HashTable, np.uint64),
        (ht.Float64HashTable, np.float64),
        (ht.Int32HashTable, np.int32),
        (ht.UInt32HashTable, np.uint32),
        (ht.Float32HashTable, np.float32),
        (ht.Int16HashTable, np.int16),
        (ht.UInt16HashTable, np.uint16),
        (ht.Int8HashTable, np.int8),
        (ht.UInt8HashTable, np.uint8),
    ],
)
class TestHashTable:
    def test_get_set_contains_len(self, table_type, dtype):
        index = 5
        table = table_type(55)
        assert len(table) == 0
        assert index not in table

        table.set_item(index, 42)
        assert len(table) == 1
        assert index in table
        assert table.get_item(index) == 42

        table.set_item(index + 1, 41)
        assert index in table
        assert index + 1 in table
        assert len(table) == 2
        assert table.get_item(index) == 42
        assert table.get_item(index + 1) == 41

        table.set_item(index, 21)
        assert index in table
        assert index + 1 in table
        assert len(table) == 2
        assert table.get_item(index) == 21
        assert table.get_item(index + 1) == 41
        assert index + 2 not in table

        with pytest.raises(KeyError) as excinfo:
            table.get_item(index + 2)
        assert str(index + 2) in str(excinfo.value)

    def test_map(self, table_type, dtype):
        # PyObjectHashTable has no map-method
        if table_type != ht.PyObjectHashTable:
            N = 77
            table = table_type()
            keys = np.arange(N).astype(dtype)
            vals = np.arange(N).astype(np.int64) + N
            table.map(keys, vals)
            for i in range(N):
                assert table.get_item(keys[i]) == i + N

    def test_map_locations(self, table_type, dtype):
        N = 8
        table = table_type()
        keys = (np.arange(N) + N).astype(dtype)
        table.map_locations(keys)
        for i in range(N):
            assert table.get_item(keys[i]) == i

    def test_lookup(self, table_type, dtype):
        N = 3
        table = table_type()
        keys = (np.arange(N) + N).astype(dtype)
        table.map_locations(keys)
        result = table.lookup(keys)
        expected = np.arange(N)
        tm.assert_numpy_array_equal(result.astype(np.int64), expected.astype(np.int64))

    def test_lookup_wrong(self, table_type, dtype):
        if dtype in (np.int8, np.uint8):
            N = 100
        else:
            N = 512
        table = table_type()
        keys = (np.arange(N) + N).astype(dtype)
        table.map_locations(keys)
        wrong_keys = np.arange(N).astype(dtype)
        result = table.lookup(wrong_keys)
        assert np.all(result == -1)

    def test_unique(self, table_type, dtype):
        if dtype in (np.int8, np.uint8):
            N = 88
        else:
            N = 1000
        table = table_type()
        expected = (np.arange(N) + N).astype(dtype)
        keys = np.repeat(expected, 5)
        unique = table.unique(keys)
        tm.assert_numpy_array_equal(unique, expected)

    def test_tracemalloc_works(self, table_type, dtype):
        if dtype in (np.int8, np.uint8):
            N = 256
        else:
            N = 30000
        keys = np.arange(N).astype(dtype)
        with activated_tracemalloc():
            table = table_type()
            table.map_locations(keys)
            used = get_allocated_khash_memory()
            my_size = table.sizeof()
            assert used == my_size
            del table
            assert get_allocated_khash_memory() == 0

    def test_tracemalloc_for_empty(self, table_type, dtype):
        with activated_tracemalloc():
            table = table_type()
            used = get_allocated_khash_memory()
            my_size = table.sizeof()
            assert used == my_size
            del table
            assert get_allocated_khash_memory() == 0


def test_tracemalloc_works_for_StringHashTable():
    N = 1000
    keys = np.arange(N).astype(np.compat.unicode).astype(np.object_)
    with activated_tracemalloc():
        table = ht.StringHashTable()
        table.map_locations(keys)
        used = get_allocated_khash_memory()
        my_size = table.sizeof()
        assert used == my_size
        del table
        assert get_allocated_khash_memory() == 0


def test_tracemalloc_for_empty_StringHashTable():
    with activated_tracemalloc():
        table = ht.StringHashTable()
        used = get_allocated_khash_memory()
        my_size = table.sizeof()
        assert used == my_size
        del table
        assert get_allocated_khash_memory() == 0


@pytest.mark.parametrize(
    "table_type, dtype",
    [
        (ht.Float64HashTable, np.float64),
        (ht.Float32HashTable, np.float32),
    ],
)
class TestHashTableWithNans:
    def test_get_set_contains_len(self, table_type, dtype):
        index = float("nan")
        table = table_type()
        assert index not in table

        table.set_item(index, 42)
        assert len(table) == 1
        assert index in table
        assert table.get_item(index) == 42

        table.set_item(index, 41)
        assert len(table) == 1
        assert index in table
        assert table.get_item(index) == 41

    def test_map(self, table_type, dtype):
        N = 332
        table = table_type()
        keys = np.full(N, np.nan, dtype=dtype)
        vals = (np.arange(N) + N).astype(np.int64)
        table.map(keys, vals)
        assert len(table) == 1
        assert table.get_item(np.nan) == 2 * N - 1

    def test_map_locations(self, table_type, dtype):
        N = 10
        table = table_type()
        keys = np.full(N, np.nan, dtype=dtype)
        table.map_locations(keys)
        assert len(table) == 1
        assert table.get_item(np.nan) == N - 1

    def test_unique(self, table_type, dtype):
        N = 1020
        table = table_type()
        keys = np.full(N, np.nan, dtype=dtype)
        unique = table.unique(keys)
        assert np.all(np.isnan(unique)) and len(unique) == 1


def get_ht_function(fun_name, type_suffix):
    return getattr(ht, fun_name + "_" + type_suffix)


@pytest.mark.parametrize(
    "dtype, type_suffix",
    [
        (np.object_, "object"),
        (np.int64, "int64"),
        (np.uint64, "uint64"),
        (np.float64, "float64"),
        (np.int32, "int32"),
        (np.uint32, "uint32"),
        (np.float32, "float32"),
        (np.int16, "int16"),
        (np.uint16, "uint16"),
        (np.int8, "int8"),
        (np.uint8, "uint8"),
    ],
)
class TestHelpFunctions:
    def test_value_count(self, dtype, type_suffix):
        N = 43
        value_count = get_ht_function("value_count", type_suffix)
        expected = (np.arange(N) + N).astype(dtype)
        values = np.repeat(expected, 5)
        keys, counts = value_count(values, False)
        tm.assert_numpy_array_equal(np.sort(keys), expected)
        assert np.all(counts == 5)

    def test_duplicated_first(self, dtype, type_suffix):
        N = 100
        duplicated = get_ht_function("duplicated", type_suffix)
        values = np.repeat(np.arange(N).astype(dtype), 5)
        result = duplicated(values)
        expected = np.ones_like(values, dtype=np.bool_)
        expected[::5] = False
        tm.assert_numpy_array_equal(result, expected)

    def test_ismember_yes(self, dtype, type_suffix):
        N = 127
        ismember = get_ht_function("ismember", type_suffix)
        arr = np.arange(N).astype(dtype)
        values = np.arange(N).astype(dtype)
        result = ismember(arr, values)
        expected = np.ones_like(values, dtype=np.bool_)
        tm.assert_numpy_array_equal(result, expected)

    def test_ismember_no(self, dtype, type_suffix):
        N = 17
        ismember = get_ht_function("ismember", type_suffix)
        arr = np.arange(N).astype(dtype)
        values = (np.arange(N) + N).astype(dtype)
        result = ismember(arr, values)
        expected = np.zeros_like(values, dtype=np.bool_)
        tm.assert_numpy_array_equal(result, expected)

    def test_mode(self, dtype, type_suffix):
        if dtype in (np.int8, np.uint8):
            N = 53
        else:
            N = 11111
        mode = get_ht_function("mode", type_suffix)
        values = np.repeat(np.arange(N).astype(dtype), 5)
        values[0] = 42
        result = mode(values, False)
        assert result == 42


@pytest.mark.parametrize(
    "dtype, type_suffix",
    [
        (np.float64, "float64"),
        (np.float32, "float32"),
    ],
)
class TestHelpFunctionsWithNans:
    def test_value_count(self, dtype, type_suffix):
        value_count = get_ht_function("value_count", type_suffix)
        values = np.array([np.nan, np.nan, np.nan], dtype=dtype)
        keys, counts = value_count(values, True)
        assert len(keys) == 0
        keys, counts = value_count(values, False)
        assert len(keys) == 1 and np.all(np.isnan(keys))
        assert counts[0] == 3

    def test_duplicated_first(self, dtype, type_suffix):
        duplicated = get_ht_function("duplicated", type_suffix)
        values = np.array([np.nan, np.nan, np.nan], dtype=dtype)
        result = duplicated(values)
        expected = np.array([False, True, True])
        tm.assert_numpy_array_equal(result, expected)

    def test_ismember_yes(self, dtype, type_suffix):
        ismember = get_ht_function("ismember", type_suffix)
        arr = np.array([np.nan, np.nan, np.nan], dtype=dtype)
        values = np.array([np.nan, np.nan], dtype=dtype)
        result = ismember(arr, values)
        expected = np.array([True, True, True], dtype=np.bool_)
        tm.assert_numpy_array_equal(result, expected)

    def test_ismember_no(self, dtype, type_suffix):
        ismember = get_ht_function("ismember", type_suffix)
        arr = np.array([np.nan, np.nan, np.nan], dtype=dtype)
        values = np.array([1], dtype=dtype)
        result = ismember(arr, values)
        expected = np.array([False, False, False], dtype=np.bool_)
        tm.assert_numpy_array_equal(result, expected)

    def test_mode(self, dtype, type_suffix):
        mode = get_ht_function("mode", type_suffix)
        values = np.array([42, np.nan, np.nan, np.nan], dtype=dtype)
        assert mode(values, True) == 42
        assert np.isnan(mode(values, False))
