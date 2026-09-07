import numpy as np

import pandas as pd


class TestSeriesCount:
    def test_count(self, datetime_series):
        assert datetime_series.count() == len(datetime_series)

        datetime_series[::2] = np.nan

        assert datetime_series.count() == np.isfinite(datetime_series).sum()

    def test_count_categorical(self):
        ser = pd.Series(
            pd.Categorical(
                [np.nan, 1, 2, np.nan], categories=[5, 4, 3, 2, 1], ordered=True
            )
        )
        result = ser.count()
        assert result == 2
