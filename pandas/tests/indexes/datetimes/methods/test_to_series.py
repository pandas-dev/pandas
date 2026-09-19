import numpy as np

import pandas as pd
import pandas._testing as tm


class TestToSeries:
    def test_to_series(self):
        naive = pd.DatetimeIndex(["2013-1-1 13:00", "2013-1-2 14:00"], name="B")
        idx = naive.tz_localize("US/Pacific")

        expected = pd.Series(np.array(idx.tolist(), dtype="object"), name="B")
        result = idx.to_series(index=range(2))
        assert expected.dtype == idx.dtype
        tm.assert_series_equal(result, expected)
