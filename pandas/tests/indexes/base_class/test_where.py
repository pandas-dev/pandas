import numpy as np

import pandas as pd
import pandas._testing as tm


class TestWhere:
    def test_where_intlike_str_doesnt_cast_ints(self):
        idx = pd.Index(range(3))
        mask = np.array([True, False, True])
        res = idx.where(mask, "2")
        expected = pd.Index([0, "2", 2])
        tm.assert_index_equal(res, expected)

    def test_where_preserves_object_dtype(self):
        # https://github.com/pandas-dev/pandas/pull/65653
        idx = pd.Index(list("abc"), dtype="object")
        result = idx.where(np.array([True, False, True]), other="x")
        expected = pd.Index(["a", "x", "c"], dtype="object")
        tm.assert_index_equal(result, expected)

    def test_where_index_with_index(self):
        # GH 65685
        idx = pd.Index(range(48))
        mask = np.ones(48, dtype=bool)

        result = idx.where(mask, idx)
        expected = idx.copy()

        tm.assert_index_equal(result, expected, exact=False)

    def test_where_index_with_ea_index(self):
        # GH 65685
        # Make BOTH indexes an Extension Array (Int64)
        idx = pd.Index([1, 2, 3], dtype="Int64")
        mask = np.array([False, False, False])

        other = pd.Index([4, 5, 6], dtype="Int64")

        # If mask is all False, result should be exactly 'other'
        result = idx.where(mask, other)
        expected = pd.Index([4, 5, 6], dtype="Int64")

        tm.assert_index_equal(result, expected)
