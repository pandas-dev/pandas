import numpy as np

import pandas as pd
from pandas import (
    Categorical,
    Series,
)


class TestSeriesCount:
    def test_count(self, datetime_series):
        assert datetime_series.count() == len(datetime_series)

        datetime_series[::2] = np.NaN

        assert datetime_series.count() == np.isfinite(datetime_series).sum()

        # GH#29478
        with pd.option_context("use_inf_as_na", True):
            assert Series([pd.Timestamp("1990/1/1")]).count() == 1

    def test_count_categorical(self):
        ser = Series(
            Categorical(
                [np.nan, 1, 2, np.nan], categories=[5, 4, 3, 2, 1], ordered=True
            )
        )
        result = ser.count()
        assert result == 2
