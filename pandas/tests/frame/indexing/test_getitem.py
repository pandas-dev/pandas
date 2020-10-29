import numpy as np
import pytest

from pandas import DataFrame, MultiIndex, period_range
import pandas._testing as tm


class TestGetitem:
    def test_getitem_unused_level_raises(self):
        # GH#20410
        mi = MultiIndex(
            levels=[["a_lot", "onlyone", "notevenone"], [1970, ""]],
            codes=[[1, 0], [1, 0]],
        )
        df = DataFrame(-1, index=range(3), columns=mi)

        with pytest.raises(KeyError, match="notevenone"):
            df["notevenone"]

    def test_getitem_periodindex(self):
        rng = period_range("1/1/2000", periods=5)
        df = DataFrame(np.random.randn(10, 5), columns=rng)

        ts = df[rng[0]]
        tm.assert_series_equal(ts, df.iloc[:, 0])

        # GH#1211; smoketest unrelated to the rest of this test
        repr(df)

        ts = df["1/1/2000"]
        tm.assert_series_equal(ts, df.iloc[:, 0])
