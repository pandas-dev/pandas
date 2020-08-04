from pandas import DataFrame, Index, MultiIndex
import pandas._testing as tm


class TestDropLevel:
    def test_droplevel(self):
        # GH#20342
        df = DataFrame([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
        df = df.set_index([0, 1]).rename_axis(["a", "b"])
        df.columns = MultiIndex.from_tuples(
            [("c", "e"), ("d", "f")], names=["level_1", "level_2"]
        )

        # test that dropping of a level in index works
        expected = df.reset_index("a", drop=True)
        result = df.droplevel("a", axis="index")
        tm.assert_frame_equal(result, expected)

        # test that dropping of a level in columns works
        expected = df.copy()
        expected.columns = Index(["c", "d"], name="level_1")
        result = df.droplevel("level_2", axis="columns")
        tm.assert_frame_equal(result, expected)
