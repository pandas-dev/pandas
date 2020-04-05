from pandas import DataFrame
import pandas._testing as tm


class TestCopy:
    def test_cache_on_copy(self):
        df = DataFrame({"a": [1]})

        df["x"] = [0]
        df["a"]

        df.copy()

        df["a"].values[0] = -1

        tm.assert_frame_equal(df, DataFrame({"a": [-1], "x": [0]}))

        df["y"] = [0]

        assert df["a"].values[0] == -1
        tm.assert_frame_equal(df, DataFrame({"a": [-1], "x": [0], "y": [0]}))
