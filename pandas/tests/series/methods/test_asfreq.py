import numpy as np
import pytest

from pandas import DataFrame, Series, date_range, period_range
import pandas._testing as tm


class TestAsFreq:
    # TODO: de-duplicate/parametrize or move DataFrame test
    def test_asfreq_ts(self):
        index = period_range(freq="A", start="1/1/2001", end="12/31/2010")
        ts = Series(np.random.randn(len(index)), index=index)
        df = DataFrame(np.random.randn(len(index), 3), index=index)

        result = ts.asfreq("D", how="end")
        df_result = df.asfreq("D", how="end")
        exp_index = index.asfreq("D", how="end")
        assert len(result) == len(ts)
        tm.assert_index_equal(result.index, exp_index)
        tm.assert_index_equal(df_result.index, exp_index)

        result = ts.asfreq("D", how="start")
        assert len(result) == len(ts)
        tm.assert_index_equal(result.index, index.asfreq("D", how="start"))

    @pytest.mark.parametrize("tz", ["US/Eastern", "dateutil/US/Eastern"])
    def test_tz_aware_asfreq(self, tz):
        dr = date_range("2011-12-01", "2012-07-20", freq="D", tz=tz)

        ser = Series(np.random.randn(len(dr)), index=dr)

        # it works!
        ser.asfreq("T")
