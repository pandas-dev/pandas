from pandas import (
    DataFrame,
    Series,
)
import pandas._testing as tm


class TestDowncast:
    def test_downcast(self):
        df = DataFrame({"a": [1.0, 2.0], "b": 1.5, "c": 2})
        result = df.downcast("int8")
        expected = DataFrame(
            {
                "a": Series([1, 2], dtype="int8"),
                "b": 1.5,
                "c": Series([2, 2], dtype="int8"),
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_downcast_dict(self):
        df = DataFrame({"a": [1.0, 2.0], "b": 1.5, "c": 2, "d": 1.0})
        df_orig = df.copy()
        result = df.downcast({"a": "int8", "b": "int64", "c": "int32"})
        expected = DataFrame(
            {
                "a": Series([1, 2], dtype="int8"),
                "b": 1.5,
                "c": Series([2, 2], dtype="int32"),
                "d": 1.0,
            }
        )
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(df, df_orig)

    def test_downcast_infer(self):
        df = DataFrame({"a": [1.0, 2.0], "b": 1.5, "c": 2.0})
        result = df.downcast("infer")
        expected = DataFrame({"a": [1, 2], "b": 1.5, "c": 2})
        tm.assert_frame_equal(result, expected)
