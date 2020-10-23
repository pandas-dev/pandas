import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Series, period_range
import pandas._testing as tm


class TestSeriesPeriod:
    def setup_method(self, method):
        self.series = Series(period_range("2000-01-01", periods=10, freq="D"))

    def test_isna(self):
        # GH 13737
        s = Series([pd.Period("2011-01", freq="M"), pd.Period("NaT", freq="M")])
        tm.assert_series_equal(s.isna(), Series([False, True]))
        tm.assert_series_equal(s.notna(), Series([True, False]))

    # ---------------------------------------------------------------------
    # NaT support

    @pytest.mark.xfail(reason="PeriodDtype Series not supported yet")
    def test_NaT_scalar(self):
        series = Series([0, 1000, 2000, pd._libs.iNaT], dtype="period[D]")

        val = series[3]
        assert pd.isna(val)

        series[2] = val
        assert pd.isna(series[2])

    def test_NaT_cast(self):
        result = Series([np.nan]).astype("period[D]")
        expected = Series([pd.NaT], dtype="period[D]")
        tm.assert_series_equal(result, expected)

    def test_intercept_astype_object(self):
        expected = self.series.astype("object")

        df = DataFrame({"a": self.series, "b": np.random.randn(len(self.series))})

        result = df.values.squeeze()
        assert (result[:, 0] == expected.values).all()

        df = DataFrame({"a": self.series, "b": ["foo"] * len(self.series)})

        result = df.values.squeeze()
        assert (result[:, 0] == expected.values).all()
