from datetime import datetime

import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
from pandas import (
    Series,
    date_range,
)
import pandas._testing as tm


class TestTruncate:
    def test_truncate_datetimeindex_tz(self):
        # GH 9243
        idx = date_range("4/1/2005", "4/30/2005", freq="D", tz="US/Pacific")
        s = Series(range(len(idx)), index=idx)
        with pytest.raises(TypeError, match="Cannot compare tz-naive"):
            # GH#36148 as of 2.0 we require tzawareness compat
            s.truncate(datetime(2005, 4, 2), datetime(2005, 4, 4))

        lb = idx[1]
        ub = idx[3]
        result = s.truncate(lb.to_pydatetime(), ub.to_pydatetime())
        expected = Series([1, 2, 3], index=idx[1:4])
        tm.assert_series_equal(result, expected)

    def test_truncate_periodindex(self):
        # GH 17717
        idx1 = pd.PeriodIndex(
            [pd.Period("2017-09-02"), pd.Period("2017-09-02"), pd.Period("2017-09-03")]
        )
        series1 = Series([1, 2, 3], index=idx1)
        result1 = series1.truncate(after="2017-09-02")

        expected_idx1 = pd.PeriodIndex(
            [pd.Period("2017-09-02"), pd.Period("2017-09-02")]
        )
        tm.assert_series_equal(result1, Series([1, 2], index=expected_idx1))

        idx2 = pd.PeriodIndex(
            [pd.Period("2017-09-03"), pd.Period("2017-09-02"), pd.Period("2017-09-03")]
        )
        series2 = Series([1, 2, 3], index=idx2)
        result2 = series2.sort_index().truncate(after="2017-09-02")

        expected_idx2 = pd.PeriodIndex([pd.Period("2017-09-02")])
        tm.assert_series_equal(result2, Series([2], index=expected_idx2))

    def test_truncate_one_element_series(self):
        # GH 35544
        series = Series([0.1], index=pd.DatetimeIndex(["2020-08-04"]))
        before = pd.Timestamp("2020-08-02")
        after = pd.Timestamp("2020-08-04")

        result = series.truncate(before=before, after=after)

        # the input Series and the expected Series are the same
        tm.assert_series_equal(result, series)

    def test_truncate_index_only_one_unique_value(self):
        # GH 42365
        obj = Series(0, index=date_range("2021-06-30", "2021-06-30")).repeat(5)

        truncated = obj.truncate("2021-06-28", "2021-07-01")

        tm.assert_series_equal(truncated, obj)


@pytest.mark.parametrize("freq", ["Q-DEC", "Q-FEB"])
def test_truncate_periodindex_quarterly_string_no_deprecation_warning(freq):
    # GH#50907 truncate parses its bounds with to_datetime internally; on a
    # PeriodIndex the user is slicing by Period, so that must not warn
    pi = pd.period_range("2000Q1", periods=12, freq=freq)
    ser = Series(range(len(pi)), index=pi)

    with tm.assert_produces_warning(None):
        result = ser.truncate("2001Q1", "2001Q3")
    assert len(result) == 3


def test_truncate_datetimeindex_quarterly_string_still_warns():
    # GH#50907 on a DatetimeIndex the user did write a quarterly string, so the
    # suppression above must not reach this case
    dti = date_range("2001-01-01", periods=500, freq="D")
    ser = Series(range(len(dti)), index=dti)

    with tm.assert_produces_warning(
        Pandas4Warning, match="quarterly string is deprecated"
    ):
        result = ser.truncate("2001Q1", "2001Q3")
    assert len(result) == 182
