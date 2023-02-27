from pandas import DataFrame
import pandas._testing as tm


class TestDowncast:
    def test_downcast(self):
        df = DataFrame({"a": [1.0, 2.0], "b": 1.5, "c": 2.0, "d": "a"})
        result = df.downcast()
        expected = DataFrame({"a": [1, 2], "b": 1.5, "c": 2, "d": "a"})
        tm.assert_frame_equal(result, expected)
