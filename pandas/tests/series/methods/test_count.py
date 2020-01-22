import numpy as np

import pandas as pd
from pandas import Categorical, MultiIndex, Series
import pandas._testing as tm


class TestSeriesCount:
    def test_count(self, datetime_series):
        assert datetime_series.count() == len(datetime_series)

        datetime_series[::2] = np.NaN

        assert datetime_series.count() == np.isfinite(datetime_series).sum()

        mi = MultiIndex.from_arrays([list("aabbcc"), [1, 2, 2, np.nan, 1, 2]])
        ts = Series(np.arange(len(mi)), index=mi)

        left = ts.count(level=1)
        right = Series([2, 3, 1], index=[1, 2, np.nan])
        tm.assert_series_equal(left, right)

        ts.iloc[[0, 3, 5]] = np.nan
        tm.assert_series_equal(ts.count(level=1), right - 1)

        # GH#29478
        with pd.option_context("use_inf_as_na", True):
            assert pd.Series([pd.Timestamp("1990/1/1")]).count() == 1

    def test_count_categorical(self):

        ser = Series(
            Categorical(
                [np.nan, 1, 2, np.nan], categories=[5, 4, 3, 2, 1], ordered=True
            )
        )
        result = ser.count()
        assert result == 2
