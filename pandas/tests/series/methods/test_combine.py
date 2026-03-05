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

    def test_combine_fill_value_with_mixed_index(self):
        # GH#31142
        # Ensure combine with fill_value works correctly with mixed-type indices
        
        a = Series([1, 2, 3], index=["a", "b", "c"])
        b = Series([10, 20, 30], index=[0, "e", "f"])
        
        result = b.combine(a, lambda x, y: x + y, fill_value=0)
        
        # Index 0 should use fill_value for 'a' value (0 + 10 = 10, not 1 + 10 = 11)
        # The bug was that it incorrectly added b[0] with a['a'] when they shouldn't align
        expected = Series([10, 1, 2, 3, 20, 30], index=[0, "a", "b", "c", "e", "f"])
        tm.assert_series_equal(result, expected)
        
        # Test with integer index that doesn't match
        b2 = Series([10, 20, 30], index=[4, "e", "f"])
        result2 = b2.combine(a, lambda x, y: x + y, fill_value=0)
        expected2 = Series([10, 1, 2, 3, 20, 30], index=[4, "a", "b", "c", "e", "f"])
        tm.assert_series_equal(result2, expected2)
