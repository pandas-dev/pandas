"""
Tests for DataFrame timezone-related methods
"""
import numpy as np
import pytest

from pandas.core.dtypes.dtypes import DatetimeTZDtype

from pandas import DataFrame, Series
import pandas._testing as tm
from pandas.core.indexes.datetimes import date_range


class TestDataFrameTimezones:
    @pytest.mark.parametrize("tz", ["US/Eastern", "dateutil/US/Eastern"])
    def test_frame_no_datetime64_dtype(self, tz):
        # after GH#7822
        # these retain the timezones on dict construction
        dr = date_range("2011/1/1", "2012/1/1", freq="W-FRI")
        dr_tz = dr.tz_localize(tz)
        df = DataFrame({"A": "foo", "B": dr_tz}, index=dr)
        tz_expected = DatetimeTZDtype("ns", dr_tz.tzinfo)
        assert df["B"].dtype == tz_expected

        # GH#2810 (with timezones)
        datetimes_naive = [ts.to_pydatetime() for ts in dr]
        datetimes_with_tz = [ts.to_pydatetime() for ts in dr_tz]
        df = DataFrame({"dr": dr})
        df["dr_tz"] = dr_tz
        df["datetimes_naive"] = datetimes_naive
        df["datetimes_with_tz"] = datetimes_with_tz
        result = df.dtypes
        expected = Series(
            [
                np.dtype("datetime64[ns]"),
                DatetimeTZDtype(tz=tz),
                np.dtype("datetime64[ns]"),
                DatetimeTZDtype(tz=tz),
            ],
            index=["dr", "dr_tz", "datetimes_naive", "datetimes_with_tz"],
        )
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("tz", ["US/Eastern", "dateutil/US/Eastern"])
    def test_frame_reset_index(self, tz):
        dr = date_range("2012-06-02", periods=10, tz=tz)
        df = DataFrame(np.random.randn(len(dr)), dr)
        roundtripped = df.reset_index().set_index("index")
        xp = df.index.tz
        rs = roundtripped.index.tz
        assert xp == rs

    @pytest.mark.parametrize("tz", [None, "America/New_York"])
    def test_boolean_compare_transpose_tzindex_with_dst(self, tz):
        # GH 19970
        idx = date_range("20161101", "20161130", freq="4H", tz=tz)
        df = DataFrame({"a": range(len(idx)), "b": range(len(idx))}, index=idx)
        result = df.T == df.T
        expected = DataFrame(True, index=list("ab"), columns=idx)
        tm.assert_frame_equal(result, expected)
