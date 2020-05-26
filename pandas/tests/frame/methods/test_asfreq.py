from datetime import datetime

import numpy as np

from pandas import DataFrame, DatetimeIndex, Series, date_range
import pandas._testing as tm

from pandas.tseries import offsets


class TestAsFreq:
    def test_asfreq(self, datetime_frame):
        offset_monthly = datetime_frame.asfreq(offsets.BMonthEnd())
        rule_monthly = datetime_frame.asfreq("BM")

        tm.assert_almost_equal(offset_monthly["A"], rule_monthly["A"])

        filled = rule_monthly.asfreq("B", method="pad")  # noqa
        # TODO: actually check that this worked.

        # don't forget!
        filled_dep = rule_monthly.asfreq("B", method="pad")  # noqa

        # test does not blow up on length-0 DataFrame
        zero_length = datetime_frame.reindex([])
        result = zero_length.asfreq("BM")
        assert result is not zero_length

    def test_asfreq_datetimeindex(self):
        df = DataFrame(
            {"A": [1, 2, 3]},
            index=[datetime(2011, 11, 1), datetime(2011, 11, 2), datetime(2011, 11, 3)],
        )
        df = df.asfreq("B")
        assert isinstance(df.index, DatetimeIndex)

        ts = df["A"].asfreq("B")
        assert isinstance(ts.index, DatetimeIndex)

    def test_asfreq_fillvalue(self):
        # test for fill value during upsampling, related to issue 3715

        # setup
        rng = date_range("1/1/2016", periods=10, freq="2S")
        ts = Series(np.arange(len(rng)), index=rng)
        df = DataFrame({"one": ts})

        # insert pre-existing missing value
        df.loc["2016-01-01 00:00:08", "one"] = None

        actual_df = df.asfreq(freq="1S", fill_value=9.0)
        expected_df = df.asfreq(freq="1S").fillna(9.0)
        expected_df.loc["2016-01-01 00:00:08", "one"] = None
        tm.assert_frame_equal(expected_df, actual_df)

        expected_series = ts.asfreq(freq="1S").fillna(9.0)
        actual_series = ts.asfreq(freq="1S", fill_value=9.0)
        tm.assert_series_equal(expected_series, actual_series)
