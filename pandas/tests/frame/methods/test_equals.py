from pandas import DataFrame
import pandas._testing as tm


class TestEquals:
    def test_dataframe_not_equal(self):
        # see GH#28839
        df1 = DataFrame({"a": [1, 2], "b": ["s", "d"]})
        df2 = DataFrame({"a": ["s", "d"], "b": [1, 2]})
        assert df1.equals(df2) is False

    def test_equals_different_blocks(self):
        # GH#9330
        df0 = DataFrame({"A": ["x", "y"], "B": [1, 2], "C": ["w", "z"]})
        df1 = df0.reset_index()[["A", "B", "C"]]
        # this assert verifies that the above operations have
        # induced a block rearrangement
        assert df0._mgr.blocks[0].dtype != df1._mgr.blocks[0].dtype

        # do the real tests
        tm.assert_frame_equal(df0, df1)
        assert df0.equals(df1)
        assert df1.equals(df0)
