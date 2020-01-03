import numpy as np

from pandas._libs import algos as libalgos, index as libindex

import pandas._testing as tm


class TestNumericEngine:
    def test_is_monotonic(self, numeric_indexing_engine_type_and_dtype):
        engine_type, dtype = numeric_indexing_engine_type_and_dtype
        num = 1000
        arr = np.array([1] * num + [2] * num + [3] * num, dtype=dtype)

        # monotonic increasing
        engine = engine_type(lambda: arr, len(arr))
        assert engine.is_monotonic_increasing is True
        assert engine.is_monotonic_decreasing is False

        # monotonic decreasing
        engine = engine_type(lambda: arr[::-1], len(arr))
        assert engine.is_monotonic_increasing is False
        assert engine.is_monotonic_decreasing is True

        # neither monotonic increasing or decreasing
        arr = np.array([1] * num + [2] * num + [1] * num, dtype=dtype)
        engine = engine_type(lambda: arr[::-1], len(arr))
        assert engine.is_monotonic_increasing is False
        assert engine.is_monotonic_decreasing is False

    def test_is_unique(self, numeric_indexing_engine_type_and_dtype):
        engine_type, dtype = numeric_indexing_engine_type_and_dtype

        # unique
        arr = np.array([1, 3, 2], dtype=dtype)
        engine = engine_type(lambda: arr, len(arr))
        assert engine.is_unique is True

        # not unique
        arr = np.array([1, 2, 1], dtype=dtype)
        engine = engine_type(lambda: arr, len(arr))
        assert engine.is_unique is False

    def test_get_loc(self, numeric_indexing_engine_type_and_dtype):
        engine_type, dtype = numeric_indexing_engine_type_and_dtype

        # unique
        arr = np.array([1, 2, 3], dtype=dtype)
        engine = engine_type(lambda: arr, len(arr))
        assert engine.get_loc(2) == 1

        # monotonic
        num = 1000
        arr = np.array([1] * num + [2] * num + [3] * num, dtype=dtype)
        engine = engine_type(lambda: arr, len(arr))
        assert engine.get_loc(2) == slice(1000, 2000)

        # not monotonic
        arr = np.array([1, 2, 3] * num, dtype=dtype)
        engine = engine_type(lambda: arr, len(arr))
        expected = np.array([False, True, False] * num, dtype=bool)
        result = engine.get_loc(2)
        assert (result == expected).all()

    def test_get_backfill_indexer(self, numeric_indexing_engine_type_and_dtype):
        engine_type, dtype = numeric_indexing_engine_type_and_dtype

        arr = np.array([1, 5, 10], dtype=dtype)
        engine = engine_type(lambda: arr, len(arr))

        new = np.arange(12, dtype=dtype)
        result = engine.get_backfill_indexer(new)

        expected = libalgos.backfill(arr, new)
        tm.assert_numpy_array_equal(result, expected)

    def test_get_pad_indexer(self, numeric_indexing_engine_type_and_dtype):
        engine_type, dtype = numeric_indexing_engine_type_and_dtype

        arr = np.array([1, 5, 10], dtype=dtype)
        engine = engine_type(lambda: arr, len(arr))

        new = np.arange(12, dtype=dtype)
        result = engine.get_pad_indexer(new)

        expected = libalgos.pad(arr, new)
        tm.assert_numpy_array_equal(result, expected)


class TestObjectEngine:
    engine_type = libindex.ObjectEngine
    dtype = np.object_
    values = list("abc")

    def test_is_monotonic(self):

        num = 1000
        arr = np.array(["a"] * num + ["a"] * num + ["c"] * num, dtype=self.dtype)

        # monotonic increasing
        engine = self.engine_type(lambda: arr, len(arr))
        assert engine.is_monotonic_increasing is True
        assert engine.is_monotonic_decreasing is False

        # monotonic decreasing
        engine = self.engine_type(lambda: arr[::-1], len(arr))
        assert engine.is_monotonic_increasing is False
        assert engine.is_monotonic_decreasing is True

        # neither monotonic increasing or decreasing
        arr = np.array(["a"] * num + ["b"] * num + ["a"] * num, dtype=self.dtype)
        engine = self.engine_type(lambda: arr[::-1], len(arr))
        assert engine.is_monotonic_increasing is False
        assert engine.is_monotonic_decreasing is False

    def test_is_unique(self):
        # unique
        arr = np.array(self.values, dtype=self.dtype)
        engine = self.engine_type(lambda: arr, len(arr))
        assert engine.is_unique is True

        # not unique
        arr = np.array(["a", "b", "a"], dtype=self.dtype)
        engine = self.engine_type(lambda: arr, len(arr))
        assert engine.is_unique is False

    def test_get_loc(self):
        # unique
        arr = np.array(self.values, dtype=self.dtype)
        engine = self.engine_type(lambda: arr, len(arr))
        assert engine.get_loc("b") == 1

        # monotonic
        num = 1000
        arr = np.array(["a"] * num + ["b"] * num + ["c"] * num, dtype=self.dtype)
        engine = self.engine_type(lambda: arr, len(arr))
        assert engine.get_loc("b") == slice(1000, 2000)

        # not monotonic
        arr = np.array(self.values * num, dtype=self.dtype)
        engine = self.engine_type(lambda: arr, len(arr))
        expected = np.array([False, True, False] * num, dtype=bool)
        result = engine.get_loc("b")
        assert (result == expected).all()

    def test_get_backfill_indexer(self):
        arr = np.array(["a", "e", "j"], dtype=self.dtype)
        engine = self.engine_type(lambda: arr, len(arr))

        new = np.array(list("abcdefghij"), dtype=self.dtype)
        result = engine.get_backfill_indexer(new)

        expected = libalgos.backfill["object"](arr, new)
        tm.assert_numpy_array_equal(result, expected)

    def test_get_pad_indexer(self):
        arr = np.array(["a", "e", "j"], dtype=self.dtype)
        engine = self.engine_type(lambda: arr, len(arr))

        new = np.array(list("abcdefghij"), dtype=self.dtype)
        result = engine.get_pad_indexer(new)

        expected = libalgos.pad["object"](arr, new)
        tm.assert_numpy_array_equal(result, expected)
