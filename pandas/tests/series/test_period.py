import numpy as np

from pandas import DataFrame, Series, period_range


class TestSeriesPeriod:

    # ---------------------------------------------------------------------
    # NaT support

    def test_intercept_astype_object(self):
        series = Series(period_range("2000-01-01", periods=10, freq="D"))

        expected = series.astype("object")

        df = DataFrame({"a": series, "b": np.random.randn(len(series))})

        result = df.values.squeeze()
        assert (result[:, 0] == expected.values).all()

        df = DataFrame({"a": series, "b": ["foo"] * len(series)})

        result = df.values.squeeze()
        assert (result[:, 0] == expected.values).all()
