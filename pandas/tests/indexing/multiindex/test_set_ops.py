from numpy.random import randn

from pandas import DataFrame, MultiIndex, Series
import pandas._testing as tm


class TestMultiIndexSetOps:
    def test_multiindex_symmetric_difference(self):
        # GH 13490
        idx = MultiIndex.from_product([["a", "b"], ["A", "B"]], names=["a", "b"])
        result = idx ^ idx
        assert result.names == idx.names

        idx2 = idx.copy().rename(["A", "B"])
        result = idx ^ idx2
        assert result.names == [None, None]

    def test_mixed_depth_insert(self):
        arrays = [
            ["a", "top", "top", "routine1", "routine1", "routine2"],
            ["", "OD", "OD", "result1", "result2", "result1"],
            ["", "wx", "wy", "", "", ""],
        ]

        tuples = sorted(zip(*arrays))
        index = MultiIndex.from_tuples(tuples)
        df = DataFrame(randn(4, 6), columns=index)

        result = df.copy()
        expected = df.copy()
        result["b"] = [1, 2, 3, 4]
        expected["b", "", ""] = [1, 2, 3, 4]
        tm.assert_frame_equal(result, expected)

    def test_dataframe_insert_column_all_na(self):
        # GH #1534
        mix = MultiIndex.from_tuples([("1a", "2a"), ("1a", "2b"), ("1a", "2c")])
        df = DataFrame([[1, 2], [3, 4], [5, 6]], index=mix)
        s = Series({(1, 1): 1, (1, 2): 2})
        df["new"] = s
        assert df["new"].isna().all()
