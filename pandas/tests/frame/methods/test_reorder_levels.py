import numpy as np

from pandas import DataFrame, MultiIndex
import pandas._testing as tm


class TestReorderLevels:
    def test_reorder_levels(self):
        index = MultiIndex(
            levels=[["bar"], ["one", "two", "three"], [0, 1]],
            codes=[[0, 0, 0, 0, 0, 0], [0, 1, 2, 0, 1, 2], [0, 1, 0, 1, 0, 1]],
            names=["L0", "L1", "L2"],
        )
        df = DataFrame({"A": np.arange(6), "B": np.arange(6)}, index=index)

        # no change, position
        result = df.reorder_levels([0, 1, 2])
        tm.assert_frame_equal(df, result)

        # no change, labels
        result = df.reorder_levels(["L0", "L1", "L2"])
        tm.assert_frame_equal(df, result)

        # rotate, position
        result = df.reorder_levels([1, 2, 0])
        e_idx = MultiIndex(
            levels=[["one", "two", "three"], [0, 1], ["bar"]],
            codes=[[0, 1, 2, 0, 1, 2], [0, 1, 0, 1, 0, 1], [0, 0, 0, 0, 0, 0]],
            names=["L1", "L2", "L0"],
        )
        expected = DataFrame({"A": np.arange(6), "B": np.arange(6)}, index=e_idx)
        tm.assert_frame_equal(result, expected)

        result = df.reorder_levels([0, 0, 0])
        e_idx = MultiIndex(
            levels=[["bar"], ["bar"], ["bar"]],
            codes=[[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]],
            names=["L0", "L0", "L0"],
        )
        expected = DataFrame({"A": np.arange(6), "B": np.arange(6)}, index=e_idx)
        tm.assert_frame_equal(result, expected)

        result = df.reorder_levels(["L0", "L0", "L0"])
        tm.assert_frame_equal(result, expected)
