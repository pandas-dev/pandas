import pandas as pd
import pandas._testing as tm


class TestCombine:
    def test_combine_scalar(self):
        # GH#21248
        ser = pd.Series([i * 10 for i in range(5)])
        result = ser.combine(3, lambda x, y: x + y)
        expected = pd.Series([i * 10 + 3 for i in range(5)])
        tm.assert_series_equal(result, expected)

        result = ser.combine(22, lambda x, y: min(x, y))
        expected = pd.Series([min(i * 10, 22) for i in range(5)])
        tm.assert_series_equal(result, expected)

    def test_combine_fill_value_with_integer_index_key(self):
        # GH#31142 - fill_value should not use positional indexing
        a = pd.Series([1, 2, 3], index=["a", "b", "c"])
        b = pd.Series([10, 20, 30], index=[0, "e", "f"])
        result = b.combine(a, lambda x, y: x + y, fill_value=0)
        expected = pd.Series(
            [10, 1, 2, 3, 20, 30],
            index=[0, "a", "b", "c", "e", "f"],
        )
        tm.assert_series_equal(result, expected)

    def test_combine_series(self):
        # GH#31899
        s1 = pd.Series([91, pd.NA, 94], dtype="Int8")
        s2 = pd.Series([91, pd.NA, 11], dtype="Int8")
        result = s1.combine(s2, lambda x, y: x + y)
        expected = pd.Series(
            [-74, pd.NA, 105], dtype="Int8"
        )  # dtype should be preserved
        tm.assert_series_equal(result, expected)

    def test_combine_non_unique_index(self):
        # https://github.com/pandas-dev/pandas/pull/67446
        # combine follows union rather than join semantics: occurrences of a
        # duplicated label are paired in order, not matched cartesian-style.
        # The second "a" has no partner in `right`, so it uses fill_value
        # and the result is [1 + 10, 2 + 0, 3 + 20].
        left = pd.Series([1, 2, 3], index=["a", "a", "b"])
        right = pd.Series([10, 20], index=["a", "b"])
        result = left.combine(right, lambda x, y: x + y, fill_value=0)
        expected = pd.Series([11, 2, 23], index=["a", "a", "b"])
        tm.assert_series_equal(result, expected)

        result = right.combine(left, lambda x, y: x + y, fill_value=0)
        tm.assert_series_equal(result, expected)

    def test_combine_non_unique_index_equal(self):
        # https://github.com/pandas-dev/pandas/pull/67446
        left = pd.Series([1, 2, 3], index=["a", "a", "b"])
        right = pd.Series([10, 20, 30], index=["a", "a", "b"])
        result = left.combine(right, lambda x, y: x + y)
        expected = pd.Series([11, 22, 33], index=["a", "a", "b"])
        tm.assert_series_equal(result, expected)

    def test_combine_non_unique_index_both_sides(self):
        # https://github.com/pandas-dev/pandas/pull/67446
        left = pd.Series([1, 2, 3, 4], index=["b", "a", "a", "a"])
        right = pd.Series([10, 20, 30, 40], index=["a", "a", "b", "c"])
        result = left.combine(right, lambda x, y: x + y, fill_value=0)
        expected = pd.Series([12, 23, 4, 31, 40], index=["a", "a", "a", "b", "c"])
        tm.assert_series_equal(result, expected)

    def test_combine_tz_naive_and_tz_aware(self):
        # https://github.com/pandas-dev/pandas/pull/67446
        left = pd.Series([1, 2], index=pd.date_range("2020", periods=2))
        right = pd.Series([10, 20], index=pd.date_range("2020", periods=2, tz="UTC"))
        result = left.combine(right, lambda x, y: x + y, fill_value=0)
        expected = pd.Series(
            [1, 2, 10, 20], index=left.index.union(right.index), dtype="int64"
        )
        tm.assert_series_equal(result, expected)

    def test_combine_multiindex_different_names(self):
        # https://github.com/pandas-dev/pandas/pull/67446
        mi1 = pd.MultiIndex.from_tuples([("a", 1), ("b", 2)], names=["x", "y"])
        mi2 = pd.MultiIndex.from_tuples([("a", 1), ("c", 3)], names=["p", "q"])
        left = pd.Series([1, 2], index=mi1)
        right = pd.Series([10, 20], index=mi2)
        result = left.combine(right, lambda x, y: x + y, fill_value=0)
        expected = pd.Series(
            [11, 2, 20], index=pd.MultiIndex.from_tuples([("a", 1), ("b", 2), ("c", 3)])
        )
        tm.assert_series_equal(result, expected)

    def test_combine_overlapping_interval_index(self):
        # https://github.com/pandas-dev/pandas/pull/67446
        # is_unique is True for an overlapping IntervalIndex, but get_indexer
        # is not supported; combine must use exact matching of intervals
        ii1 = pd.IntervalIndex.from_tuples([(0, 2), (1, 3)])
        ii2 = pd.IntervalIndex.from_tuples([(0, 2), (4, 5)])
        left = pd.Series([1, 2], index=ii1)
        right = pd.Series([10, 20], index=ii2)
        result = left.combine(right, lambda x, y: x + y, fill_value=0)
        expected = pd.Series(
            [11, 2, 20], index=pd.IntervalIndex.from_tuples([(0, 2), (1, 3), (4, 5)])
        )
        tm.assert_series_equal(result, expected)

    def test_combine_overlapping_interval_index_non_unique(self):
        # https://github.com/pandas-dev/pandas/pull/67446
        ii1 = pd.IntervalIndex.from_tuples([(0, 2), (1, 3), (0, 2)])
        ii2 = pd.IntervalIndex.from_tuples([(0, 2), (4, 5)])
        left = pd.Series([1, 2, 3], index=ii1)
        right = pd.Series([10, 20], index=ii2)
        result = left.combine(right, lambda x, y: x + y, fill_value=0)
        expected = pd.Series(
            [11, 3, 2, 20],
            index=pd.IntervalIndex.from_tuples([(0, 2), (0, 2), (1, 3), (4, 5)]),
        )
        tm.assert_series_equal(result, expected)
