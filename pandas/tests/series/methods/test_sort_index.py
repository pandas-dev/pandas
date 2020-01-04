import random

import numpy as np
import pytest

from pandas import IntervalIndex, MultiIndex, Series
import pandas._testing as tm


class TestSeriesSortIndex:
    def test_sort_index(self, datetime_series):
        rindex = list(datetime_series.index)
        random.shuffle(rindex)

        random_order = datetime_series.reindex(rindex)
        sorted_series = random_order.sort_index()
        tm.assert_series_equal(sorted_series, datetime_series)

        # descending
        sorted_series = random_order.sort_index(ascending=False)
        tm.assert_series_equal(
            sorted_series, datetime_series.reindex(datetime_series.index[::-1])
        )

        # compat on level
        sorted_series = random_order.sort_index(level=0)
        tm.assert_series_equal(sorted_series, datetime_series)

        # compat on axis
        sorted_series = random_order.sort_index(axis=0)
        tm.assert_series_equal(sorted_series, datetime_series)

        msg = "No axis named 1 for object type <class 'pandas.core.series.Series'>"
        with pytest.raises(ValueError, match=msg):
            random_order.sort_values(axis=1)

        sorted_series = random_order.sort_index(level=0, axis=0)
        tm.assert_series_equal(sorted_series, datetime_series)

        with pytest.raises(ValueError, match=msg):
            random_order.sort_index(level=0, axis=1)

    def test_sort_index_inplace(self, datetime_series):

        # For GH#11402
        rindex = list(datetime_series.index)
        random.shuffle(rindex)

        # descending
        random_order = datetime_series.reindex(rindex)
        result = random_order.sort_index(ascending=False, inplace=True)

        assert result is None
        tm.assert_series_equal(
            random_order, datetime_series.reindex(datetime_series.index[::-1])
        )

        # ascending
        random_order = datetime_series.reindex(rindex)
        result = random_order.sort_index(ascending=True, inplace=True)

        assert result is None
        tm.assert_series_equal(random_order, datetime_series)

    def test_sort_index_level(self):
        mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list("ABC"))
        s = Series([1, 2], mi)
        backwards = s.iloc[[1, 0]]

        res = s.sort_index(level="A")
        tm.assert_series_equal(backwards, res)

        res = s.sort_index(level=["A", "B"])
        tm.assert_series_equal(backwards, res)

        res = s.sort_index(level="A", sort_remaining=False)
        tm.assert_series_equal(s, res)

        res = s.sort_index(level=["A", "B"], sort_remaining=False)
        tm.assert_series_equal(s, res)

    @pytest.mark.parametrize("level", ["A", 0])  # GH#21052
    def test_sort_index_multiindex(self, level):

        mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list("ABC"))
        s = Series([1, 2], mi)
        backwards = s.iloc[[1, 0]]

        # implicit sort_remaining=True
        res = s.sort_index(level=level)
        tm.assert_series_equal(backwards, res)

        # GH#13496
        # sort has no effect without remaining lvls
        res = s.sort_index(level=level, sort_remaining=False)
        tm.assert_series_equal(s, res)

    def test_sort_index_kind(self):
        # GH#14444 & GH#13589:  Add support for sort algo choosing
        series = Series(index=[3, 2, 1, 4, 3], dtype=object)
        expected_series = Series(index=[1, 2, 3, 3, 4], dtype=object)

        index_sorted_series = series.sort_index(kind="mergesort")
        tm.assert_series_equal(expected_series, index_sorted_series)

        index_sorted_series = series.sort_index(kind="quicksort")
        tm.assert_series_equal(expected_series, index_sorted_series)

        index_sorted_series = series.sort_index(kind="heapsort")
        tm.assert_series_equal(expected_series, index_sorted_series)

    def test_sort_index_na_position(self):
        series = Series(index=[3, 2, 1, 4, 3, np.nan], dtype=object)
        expected_series_first = Series(index=[np.nan, 1, 2, 3, 3, 4], dtype=object)

        index_sorted_series = series.sort_index(na_position="first")
        tm.assert_series_equal(expected_series_first, index_sorted_series)

        expected_series_last = Series(index=[1, 2, 3, 3, 4, np.nan], dtype=object)

        index_sorted_series = series.sort_index(na_position="last")
        tm.assert_series_equal(expected_series_last, index_sorted_series)

    def test_sort_index_intervals(self):
        s = Series(
            [np.nan, 1, 2, 3], IntervalIndex.from_arrays([0, 1, 2, 3], [1, 2, 3, 4])
        )

        result = s.sort_index()
        expected = s
        tm.assert_series_equal(result, expected)

        result = s.sort_index(ascending=False)
        expected = Series(
            [3, 2, 1, np.nan], IntervalIndex.from_arrays([3, 2, 1, 0], [4, 3, 2, 1])
        )
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("inplace", [True, False])
    @pytest.mark.parametrize(
        "original_list, sorted_list, ascending, ignore_index, output_index",
        [
            ([2, 3, 6, 1], [2, 3, 6, 1], True, True, [0, 1, 2, 3]),
            ([2, 3, 6, 1], [2, 3, 6, 1], True, False, [0, 1, 2, 3]),
            ([2, 3, 6, 1], [1, 6, 3, 2], False, True, [0, 1, 2, 3]),
            ([2, 3, 6, 1], [1, 6, 3, 2], False, False, [3, 2, 1, 0]),
        ],
    )
    def test_sort_index_ignore_index(
        self, inplace, original_list, sorted_list, ascending, ignore_index, output_index
    ):
        # GH 30114
        ser = Series(original_list)
        expected = Series(sorted_list, index=output_index)
        kwargs = {
            "ascending": ascending,
            "ignore_index": ignore_index,
            "inplace": inplace,
        }

        if inplace:
            result_ser = ser.copy()
            result_ser.sort_index(**kwargs)
        else:
            result_ser = ser.sort_index(**kwargs)

        tm.assert_series_equal(result_ser, expected)
        tm.assert_series_equal(ser, Series(original_list))
