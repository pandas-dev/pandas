import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm

from pandas.tseries.offsets import BDay


class TestTruncate:
    def test_truncate(self, datetime_series):
        offset = BDay()

        ts = datetime_series[::3]

        start, end = datetime_series.index[3], datetime_series.index[6]
        start_missing, end_missing = datetime_series.index[2], datetime_series.index[7]

        # neither specified
        truncated = ts.truncate()
        tm.assert_series_equal(truncated, ts)

        # both specified
        expected = ts[1:3]

        truncated = ts.truncate(start, end)
        tm.assert_series_equal(truncated, expected)

        truncated = ts.truncate(start_missing, end_missing)
        tm.assert_series_equal(truncated, expected)

        # start specified
        expected = ts[1:]

        truncated = ts.truncate(before=start)
        tm.assert_series_equal(truncated, expected)

        truncated = ts.truncate(before=start_missing)
        tm.assert_series_equal(truncated, expected)

        # end specified
        expected = ts[:3]

        truncated = ts.truncate(after=end)
        tm.assert_series_equal(truncated, expected)

        truncated = ts.truncate(after=end_missing)
        tm.assert_series_equal(truncated, expected)

        # corner case, empty series returned
        truncated = ts.truncate(after=datetime_series.index[0] - offset)
        assert len(truncated) == 0

        truncated = ts.truncate(before=datetime_series.index[-1] + offset)
        assert len(truncated) == 0

        msg = "Truncate: 1999-12-31 00:00:00 must be after 2000-02-14 00:00:00"
        with pytest.raises(ValueError, match=msg):
            ts.truncate(
                before=datetime_series.index[-1] + offset,
                after=datetime_series.index[0] - offset,
            )

    def test_truncate_nonsortedindex(self):
        # GH#17935

        s = pd.Series(["a", "b", "c", "d", "e"], index=[5, 3, 2, 9, 0])
        msg = "truncate requires a sorted index"

        with pytest.raises(ValueError, match=msg):
            s.truncate(before=3, after=9)

        rng = pd.date_range("2011-01-01", "2012-01-01", freq="W")
        ts = pd.Series(np.random.randn(len(rng)), index=rng)
        msg = "truncate requires a sorted index"

        with pytest.raises(ValueError, match=msg):
            ts.sort_values(ascending=False).truncate(before="2011-11", after="2011-12")
