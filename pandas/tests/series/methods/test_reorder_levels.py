import numpy as np

from pandas import MultiIndex, Series
import pandas._testing as tm


class TestReorderLevels:
    def test_reorder_levels(self):
        index = MultiIndex(
            levels=[["bar"], ["one", "two", "three"], [0, 1]],
            codes=[[0, 0, 0, 0, 0, 0], [0, 1, 2, 0, 1, 2], [0, 1, 0, 1, 0, 1]],
            names=["L0", "L1", "L2"],
        )
        s = Series(np.arange(6), index=index)

        # no change, position
        result = s.reorder_levels([0, 1, 2])
        tm.assert_series_equal(s, result)

        # no change, labels
        result = s.reorder_levels(["L0", "L1", "L2"])
        tm.assert_series_equal(s, result)

        # rotate, position
        result = s.reorder_levels([1, 2, 0])
        e_idx = MultiIndex(
            levels=[["one", "two", "three"], [0, 1], ["bar"]],
            codes=[[0, 1, 2, 0, 1, 2], [0, 1, 0, 1, 0, 1], [0, 0, 0, 0, 0, 0]],
            names=["L1", "L2", "L0"],
        )
        expected = Series(np.arange(6), index=e_idx)
        tm.assert_series_equal(result, expected)
