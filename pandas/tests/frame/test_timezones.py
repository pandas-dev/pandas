"""
Tests for DataFrame timezone-related methods
"""
import numpy as np
import pytest

from pandas.core.dtypes.dtypes import DatetimeTZDtype

import pandas as pd
from pandas import DataFrame, Series
import pandas._testing as tm
from pandas.core.indexes.datetimes import date_range


class TestDataFrameTimezones:
    def test_frame_join_tzaware(self):
        test1 = DataFrame(
            np.zeros((6, 3)),
            index=date_range(
                "2012-11-15 00:00:00", periods=6, freq="100L", tz="US/Central"
            ),
        )
        test2 = DataFrame(
            np.zeros((3, 3)),
            index=date_range(
                "2012-11-15 00:00:00", periods=3, freq="250L", tz="US/Central"
            ),
            columns=range(3, 6),
        )

        result = test1.join(test2, how="outer")
        ex_index = test1.index.union(test2.index)

        tm.assert_index_equal(result.index, ex_index)
        assert result.index.tz.zone == "US/Central"

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

    @pytest.mark.parametrize("copy", [True, False])
    @pytest.mark.parametrize(
        "method, tz", [["tz_localize", None], ["tz_convert", "Europe/Berlin"]]
    )
    def test_tz_localize_convert_copy_inplace_mutate(self, copy, method, tz):
        # GH 6326
        result = DataFrame(
            np.arange(0, 5), index=date_range("20131027", periods=5, freq="1H", tz=tz)
        )
        getattr(result, method)("UTC", copy=copy)
        expected = DataFrame(
            np.arange(0, 5), index=date_range("20131027", periods=5, freq="1H", tz=tz)
        )
        tm.assert_frame_equal(result, expected)

    def test_constructor_data_aware_dtype_naive(self, tz_aware_fixture):
        # GH 25843
        tz = tz_aware_fixture
        result = DataFrame({"d": [pd.Timestamp("2019", tz=tz)]}, dtype="datetime64[ns]")
        expected = DataFrame({"d": [pd.Timestamp("2019")]})
        tm.assert_frame_equal(result, expected)
