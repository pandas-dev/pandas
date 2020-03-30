from datetime import datetime

import numpy as np

from pandas import Series
import pandas._testing as tm


def test_reindex_like(datetime_series):
    other = datetime_series[::2]
    tm.assert_series_equal(
        datetime_series.reindex(other.index), datetime_series.reindex_like(other)
    )

    # GH#7179
    day1 = datetime(2013, 3, 5)
    day2 = datetime(2013, 5, 5)
    day3 = datetime(2014, 3, 5)

    series1 = Series([5, None, None], [day1, day2, day3])
    series2 = Series([None, None], [day1, day3])

    result = series1.reindex_like(series2, method="pad")
    expected = Series([5, np.nan], index=[day1, day3])
    tm.assert_series_equal(result, expected)


def test_reindex_like_nearest():
    s = Series(np.arange(10, dtype="int64"))
    target = [0.1, 0.9, 1.5, 2.0]
    result = s.reindex(target, method="nearest")
    expected = Series(np.around(target).astype("int64"), target)
    tm.assert_series_equal(expected, result)

    result = s.reindex_like(result, method="nearest")
    tm.assert_series_equal(expected, result)

    result = s.reindex_like(result, method="nearest", tolerance=1)
    tm.assert_series_equal(expected, result)
    result = s.reindex_like(result, method="nearest", tolerance=[1, 2, 3, 4])
    tm.assert_series_equal(expected, result)

    result = s.reindex(target, method="nearest", tolerance=0.2)
    expected = Series([0, 1, np.nan, 2], target)
    tm.assert_series_equal(expected, result)

    result = s.reindex(target, method="nearest", tolerance=[0.3, 0.01, 0.4, 3])
    expected = Series([0, np.nan, np.nan, 2], target)
    tm.assert_series_equal(expected, result)
