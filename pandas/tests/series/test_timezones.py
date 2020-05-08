"""
Tests for Series timezone-related methods
"""
from datetime import datetime

from dateutil.tz import tzoffset
import numpy as np
import pytest

from pandas import Series
import pandas._testing as tm
from pandas.core.indexes.datetimes import date_range


class TestSeriesTimezones:
    def test_dateutil_tzoffset_support(self):
        values = [188.5, 328.25]
        tzinfo = tzoffset(None, 7200)
        index = [
            datetime(2012, 5, 11, 11, tzinfo=tzinfo),
            datetime(2012, 5, 11, 12, tzinfo=tzinfo),
        ]
        series = Series(data=values, index=index)

        assert series.index.tz == tzinfo

        # it works! #2443
        repr(series.index[0])

    @pytest.mark.parametrize("copy", [True, False])
    @pytest.mark.parametrize(
        "method, tz", [["tz_localize", None], ["tz_convert", "Europe/Berlin"]]
    )
    def test_tz_localize_convert_copy_inplace_mutate(self, copy, method, tz):
        # GH 6326
        result = Series(
            np.arange(0, 5), index=date_range("20131027", periods=5, freq="1H", tz=tz)
        )
        getattr(result, method)("UTC", copy=copy)
        expected = Series(
            np.arange(0, 5), index=date_range("20131027", periods=5, freq="1H", tz=tz)
        )
        tm.assert_series_equal(result, expected)
