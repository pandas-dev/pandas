import numpy as np

from pandas import DataFrame, period_range
import pandas._testing as tm


class TestPeriodIndex:
    def test_as_frame_columns(self):
        rng = period_range("1/1/2000", periods=5)
        df = DataFrame(np.random.randn(10, 5), columns=rng)

        ts = df[rng[0]]
        tm.assert_series_equal(ts, df.iloc[:, 0])

        # GH # 1211
        repr(df)

        ts = df["1/1/2000"]
        tm.assert_series_equal(ts, df.iloc[:, 0])
