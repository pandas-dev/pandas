from pandas.errors import Pandas4Warning

from pandas import (
    NA,
    Series,
)
import pandas._testing as tm


class TestCombine:
    def test_combine_scalar_deprecated(self):
        # GH#62918 - combining with a non-Series is deprecated
        ser = Series([i * 10 for i in range(5)])

        msg = "Series.combine with a non-Series argument is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = ser.combine(3, lambda x, y: x + y)
        expected = Series([i * 10 + 3 for i in range(5)])
        tm.assert_series_equal(result, expected)

        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = ser.combine(22, lambda x, y: min(x, y))
        expected = Series([min(i * 10, 22) for i in range(5)])
        tm.assert_series_equal(result, expected)

    def test_combine_fill_value_with_integer_index_key(self):
        # GH#31142 - fill_value should not use positional indexing
        a = Series([1, 2, 3], index=["a", "b", "c"])
        b = Series([10, 20, 30], index=[0, "e", "f"])
        result = b.combine(a, lambda x, y: x + y, fill_value=0)
        expected = Series(
            [10, 1, 2, 3, 20, 30],
            index=[0, "a", "b", "c", "e", "f"],
        )
        tm.assert_series_equal(result, expected)

    def test_combine_series(self):
        # GH#31899
        s1 = Series([91, NA, 94], dtype="Int8")
        s2 = Series([91, NA, 11], dtype="Int8")
        result = s1.combine(s2, lambda x, y: x + y)
        expected = Series([-74, NA, 105], dtype="Int8")  # dtype should be preserved
        tm.assert_series_equal(result, expected)
