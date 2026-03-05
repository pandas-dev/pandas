from pandas import (
    NA,
    Series,
)
import pandas._testing as tm


class TestCombine:
    def test_combine_scalar(self):
        # GH#21248
        ser = Series([i * 10 for i in range(5)])
        result = ser.combine(3, lambda x, y: x + y)
        expected = Series([i * 10 + 3 for i in range(5)])
        tm.assert_series_equal(result, expected)

        result = ser.combine(22, lambda x, y: min(x, y))
        expected = Series([min(i * 10, 22) for i in range(5)])
        tm.assert_series_equal(result, expected)

    def test_combine_series(self):
        # GH#31899
        s1 = Series([91, NA, 94], dtype="Int8")
        s2 = Series([91, NA, 11], dtype="Int8")
        result = s1.combine(s2, lambda x, y: x + y)
        expected = Series([-74, NA, 105], dtype="Int8")  # dtype should be preserved
        tm.assert_series_equal(result, expected)
