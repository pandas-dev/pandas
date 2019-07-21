import numpy as np
import pytest

from pandas import DataFrame, Series, date_range
from pandas.tests.window.common import Base
import pandas.util.testing as tm


# Tests for GH11704
class TestExpanding(Base):
    @pytest.mark.parametrize(
        "dataframe,expected,window",
        [
            (
                DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}),
                [({"A": [1, 2, 3], "B": [4, 5, 6]}, [0, 1, 2])],
                3,
            ),
            (
                DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}),
                [
                    ({"A": [1, 2], "B": [4, 5]}, [0, 1]),
                    ({"A": [1, 2, 3], "B": [4, 5, 6]}, [0, 1, 2]),
                ],
                2,
            ),
            (
                DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}),
                [
                    ({"A": [1], "B": [4]}, [0]),
                    ({"A": [1, 2], "B": [4, 5]}, [0, 1]),
                    ({"A": [1, 2, 3], "B": [4, 5, 6]}, [0, 1, 2]),
                ],
                1,
            ),
            (DataFrame({"A": [1], "B": [4]}), [], 1337),
            (DataFrame(), [({}, [])], 1337),
        ],
    )
    def test_iterator_expanding_dataframe(self, dataframe, expected, window):
        expected = [DataFrame(values, index=index) for (values, index) in expected]

        for (expected, actual) in zip(
            expected, dataframe.expanding(min_periods=window)
        ):
            tm.assert_frame_equal(actual, expected)

    @pytest.mark.parametrize(
        "series,expected,window",
        [
            (Series([1, 2, 3]), [([1, 2, 3], [0, 1, 2])], 3),
            (Series([1, 2, 3]), [([1, 2], [0, 1]), ([1, 2, 3], [0, 1, 2])], 2),
            (
                Series([1, 2, 3]),
                [([1], [0]), ([1, 2], [0, 1]), ([1, 2, 3], [0, 1, 2])],
                1,
            ),
            (Series([1, 2]), [([1, 2], [0, 1])], 1337),
            (Series([]), [], 1337),
        ],
    )
    def test_iterator_expanding_series(self, series, expected, window):
        expected = [Series(values, index=index) for (values, index) in expected]

        for (expected, actual) in zip(expected, series.expanding(min_periods=window)):
            tm.assert_series_equal(actual, expected)

    @pytest.mark.parametrize(
        "dataframe,expected,window",
        [
            (
                DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}),
                [({"A": [1, 2, 3], "B": [4, 5, 6]}, [0, 1, 2])],
                3,
            ),
            (
                DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}),
                [
                    ({"A": [1, 2], "B": [4, 5]}, [0, 1]),
                    ({"A": [2, 3], "B": [5, 6]}, [1, 2]),
                ],
                2,
            ),
            (
                DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}),
                [
                    ({"A": [1], "B": [4]}, [0]),
                    ({"A": [2], "B": [5]}, [1]),
                    ({"A": [3], "B": [6]}, [2]),
                ],
                1,
            ),
            (DataFrame({"A": [1], "B": [4]}), [], 1337),
            (DataFrame(), [({}, [])], 1337),
        ],
    )
    def test_iterator_rolling_dataframe(self, dataframe, expected, window):
        expected = [DataFrame(values, index=index) for (values, index) in expected]

        for (expected, actual) in zip(expected, dataframe.rolling(window)):
            tm.assert_frame_equal(actual, expected)

    @pytest.mark.parametrize(
        "series,expected,window",
        [
            (Series([1, 2, 3]), [([1, 2, 3], [0, 1, 2])], 3),
            (Series([1, 2, 3]), [([1, 2], [0, 1]), ([2, 3], [1, 2])], 2),
            (Series([1, 2, 3]), [([1], [0]), ([2], [1]), ([3], [2])], 1),
            (Series([1, 2]), [([1, 2], [0, 1])], 1337),
            (Series([]), [], 1337),
        ],
    )
    def test_iterator_rolling_series(self, series, expected, window):
        expected = [Series(values, index=index) for (values, index) in expected]

        for (expected, actual) in zip(expected, series.rolling(window)):
            tm.assert_series_equal(actual, expected)

    @pytest.mark.parametrize(
        "expected,window,minp",
        [
            (
                [
                    ({"A": [1, 2], "B": [4, 5]}, [0, 1]),
                    ({"A": [1, 2, 3], "B": [4, 5, 6]}, [0, 1, 2]),
                ],
                3,
                2,
            ),
            ([], 3, 4),
        ],
    )
    def test_iterator_rolling_dataframe_min_periods(self, expected, window, minp):
        dataframe = DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})

        expected = [DataFrame(values, index=index) for (values, index) in expected]

        for (expected, actual) in zip(
            expected, dataframe.rolling(window, min_periods=minp)
        ):
            tm.assert_frame_equal(actual, expected)

    @pytest.mark.parametrize(
        "expected,window,minp",
        [([([1, 2], [0, 1]), ([1, 2, 3], [0, 1, 2])], 3, 2), ([], 3, 4)],
    )
    def test_iterator_rolling_series_min_periods(self, expected, window, minp):
        series = Series([1, 2, 3])

        expected = [Series(values, index=index) for (values, index) in expected]

        for (expected, actual) in zip(
            expected, series.rolling(window, min_periods=minp)
        ):
            tm.assert_series_equal(actual, expected)

    @pytest.mark.parametrize(
        "expected,window,minp",
        [
            ([([1.0, np.nan, 3.0], [0, 1, 2]), ([np.nan, 3.0, 4.0], [1, 2, 3])], 3, 2),
            ([([1.0, np.nan, 3.0, 4.0], [0, 1, 2, 3])], 4, 3),
            ([], 4, 4),
        ],
    )
    def test_iterator_rolling_series_min_periods_nan(self, expected, window, minp):
        series = Series([1, np.nan, 3, 4])

        expected = [Series(values, index=index) for (values, index) in expected]

        for (expected, actual) in zip(
            expected, series.rolling(window, min_periods=minp)
        ):
            tm.assert_series_equal(actual, expected)

    @pytest.mark.parametrize(
        "expected,window,minp",
        [
            ([({"A": [1, np.nan, 3], "B": [np.nan, 5, 6]}, [0, 1, 2])], 3, 2),
            ([], 2, 2),
            ([], 3, 4),
        ],
    )
    def test_iterator_rolling_dataframe_min_periods_nan(self, expected, window, minp):
        dataframe = DataFrame({"A": [1, np.nan, 3], "B": [np.nan, 5, 6]})

        expected = [DataFrame(values, index=index) for (values, index) in expected]

        for (expected, actual) in zip(
            expected, dataframe.rolling(window, min_periods=minp)
        ):
            tm.assert_frame_equal(actual, expected)

    @pytest.mark.parametrize(
        "expected,window",
        [
            ([([0], [0]), ([1], [1]), ([2], [2]), ([3], [3]), ([4], [4])], "1s"),
            (
                [
                    ([0], [0]),
                    ([0, 1], [0, 1]),
                    ([1, 2], [1, 2]),
                    ([2, 3], [2, 3]),
                    ([3, 4], [3, 4]),
                ],
                "2S",
            ),
            (
                [
                    ([0], [0]),
                    ([0, 1], [0, 1]),
                    ([0, 1, 2], [0, 1, 2]),
                    ([1, 2, 3], [1, 2, 3]),
                    ([2, 3, 4], [2, 3, 4]),
                ],
                "3S",
            ),
        ],
    )
    def test_iterator_rolling_series_time(self, expected, window):
        series = Series(
            range(5), index=date_range(start="2016-01-01 09:30:00", periods=5, freq="s")
        )

        expected = [
            Series(values, index=series.index[index]) for (values, index) in expected
        ]

        for (expected, actual) in zip(expected, series.rolling(window)):
            tm.assert_series_equal(actual, expected)

    @pytest.mark.parametrize(
        "expected,window,minp",
        [
            ([], "1s", 2),
            (
                [
                    ([0, 1], [0, 1]),
                    ([1, 2], [1, 2]),
                    ([2, 3], [2, 3]),
                    ([3, 4], [3, 4]),
                ],
                "2S",
                2,
            ),
            (
                [
                    ([0], [0]),
                    ([0, 1], [0, 1]),
                    ([0, 1, 2], [0, 1, 2]),
                    ([1, 2, 3], [1, 2, 3]),
                    ([2, 3, 4], [2, 3, 4]),
                ],
                "3S",
                1,
            ),
        ],
    )
    def test_iterator_rolling_series_time_min_periods(self, expected, window, minp):
        series = Series(
            range(5), index=date_range(start="2016-01-01 09:30:00", periods=5, freq="s")
        )

        expected = [
            Series(values, index=series.index[index]) for (values, index) in expected
        ]

        for (expected, actual) in zip(
            expected, series.rolling(window, min_periods=minp)
        ):
            tm.assert_series_equal(actual, expected)
