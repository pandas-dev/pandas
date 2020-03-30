from pandas import DataFrame, Series
import pandas._testing as tm


class TestDataFramePop:
    def test_pop(self, float_frame):
        float_frame.columns.name = "baz"

        float_frame.pop("A")
        assert "A" not in float_frame

        float_frame["foo"] = "bar"
        float_frame.pop("foo")
        assert "foo" not in float_frame
        assert float_frame.columns.name == "baz"

        # gh-10912: inplace ops cause caching issue
        a = DataFrame([[1, 2, 3], [4, 5, 6]], columns=["A", "B", "C"], index=["X", "Y"])
        b = a.pop("B")
        b += 1

        # original frame
        expected = DataFrame([[1, 3], [4, 6]], columns=["A", "C"], index=["X", "Y"])
        tm.assert_frame_equal(a, expected)

        # result
        expected = Series([2, 5], index=["X", "Y"], name="B") + 1
        tm.assert_series_equal(b, expected)

    def test_pop_non_unique_cols(self):
        df = DataFrame({0: [0, 1], 1: [0, 1], 2: [4, 5]})
        df.columns = ["a", "b", "a"]

        res = df.pop("a")
        assert type(res) == DataFrame
        assert len(res) == 2
        assert len(df.columns) == 1
        assert "b" in df.columns
        assert "a" not in df.columns
        assert len(df.index) == 2
