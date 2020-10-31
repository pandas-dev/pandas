import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Series, period_range


class TestSeriesPeriod:
    def setup_method(self, method):
        self.series = Series(period_range("2000-01-01", periods=10, freq="D"))

    # ---------------------------------------------------------------------
    # NaT support

    @pytest.mark.xfail(reason="PeriodDtype Series not supported yet")
    def test_NaT_scalar(self):
        series = Series([0, 1000, 2000, pd._libs.iNaT], dtype="period[D]")

        val = series[3]
        assert pd.isna(val)

        series[2] = val
        assert pd.isna(series[2])

    def test_intercept_astype_object(self):
        expected = self.series.astype("object")

        df = DataFrame({"a": self.series, "b": np.random.randn(len(self.series))})

        result = df.values.squeeze()
        assert (result[:, 0] == expected.values).all()

        df = DataFrame({"a": self.series, "b": ["foo"] * len(self.series)})

        result = df.values.squeeze()
        assert (result[:, 0] == expected.values).all()
