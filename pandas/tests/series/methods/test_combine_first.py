import numpy as np

from pandas import Period, Series, date_range, period_range
import pandas._testing as tm


class TestCombineFirst:
    def test_combine_first_period_datetime(self):
        # GH#3367
        didx = date_range(start="1950-01-31", end="1950-07-31", freq="M")
        pidx = period_range(start=Period("1950-1"), end=Period("1950-7"), freq="M")
        # check to be consistent with DatetimeIndex
        for idx in [didx, pidx]:
            a = Series([1, np.nan, np.nan, 4, 5, np.nan, 7], index=idx)
            b = Series([9, 9, 9, 9, 9, 9, 9], index=idx)

            result = a.combine_first(b)
            expected = Series([1, 9, 9, 4, 5, 9, 7], index=idx, dtype=np.float64)
            tm.assert_series_equal(result, expected)
