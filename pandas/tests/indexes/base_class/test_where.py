import numpy as np

from pandas import Index
import pandas._testing as tm


class TestWhere:
    def test_where_intlike_str_doesnt_cast_ints(self):
        idx = Index(range(3))
        mask = np.array([True, False, True])
        res = idx.where(mask, "2")
        expected = Index([0, "2", 2])
        tm.assert_index_equal(res, expected)

    def test_where_preserves_object_dtype(self):
        # https://github.com/pandas-dev/pandas/pull/65653
        idx = Index(list("abc"), dtype="object")
        result = idx.where(np.array([True, False, True]), other="x")
        expected = Index(["a", "x", "c"], dtype="object")
        tm.assert_index_equal(result, expected)

    def test_where_index_with_index(self):
        # GH 65685
        idx = Index(range(48))
        mask = np.ones(48, dtype=bool)

        result = idx.where(mask, idx)
        expected = idx.copy()

        tm.assert_index_equal(result, expected, exact=False)
