import numpy as np
import pytest

from pandas import Series, timedelta_range
import pandas._testing as tm


class TestSlicing:
    def test_partial_slice(self):
        rng = timedelta_range("1 day 10:11:12", freq="h", periods=500)
        s = Series(np.arange(len(rng)), index=rng)

        result = s["5 day":"6 day"]
        expected = s.iloc[86:134]
        tm.assert_series_equal(result, expected)

        result = s["5 day":]
        expected = s.iloc[86:]
        tm.assert_series_equal(result, expected)

        result = s[:"6 day"]
        expected = s.iloc[:134]
        tm.assert_series_equal(result, expected)

        result = s["6 days, 23:11:12"]
        assert result == s.iloc[133]

        msg = r"^Timedelta\('50 days 00:00:00'\)$"
        with pytest.raises(KeyError, match=msg):
            s["50 days"]

    def test_partial_slice_high_reso(self):

        # higher reso
        rng = timedelta_range("1 day 10:11:12", freq="us", periods=2000)
        s = Series(np.arange(len(rng)), index=rng)

        result = s["1 day 10:11:12":]
        expected = s.iloc[0:]
        tm.assert_series_equal(result, expected)

        result = s["1 day 10:11:12.001":]
        expected = s.iloc[1000:]
        tm.assert_series_equal(result, expected)

        result = s["1 days, 10:11:12.001001"]
        assert result == s.iloc[1001]
