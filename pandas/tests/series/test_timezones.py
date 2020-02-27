"""
Tests for Series timezone-related methods
"""
from datetime import datetime

from dateutil.tz import tzoffset
import numpy as np
import pytest
import pytz

from pandas._libs.tslibs import conversion, timezones

from pandas import Series, Timestamp
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

    @pytest.mark.parametrize("tz", ["US/Eastern", "dateutil/US/Eastern"])
    def test_string_index_alias_tz_aware(self, tz):
        rng = date_range("1/1/2000", periods=10, tz=tz)
        ser = Series(np.random.randn(len(rng)), index=rng)

        result = ser["1/3/2000"]
        tm.assert_almost_equal(result, ser[2])

    def test_series_align_aware(self):
        idx1 = date_range("2001", periods=5, freq="H", tz="US/Eastern")
        ser = Series(np.random.randn(len(idx1)), index=idx1)
        ser_central = ser.tz_convert("US/Central")
        # # different timezones convert to UTC

        new1, new2 = ser.align(ser_central)
        assert new1.index.tz == pytz.UTC
        assert new2.index.tz == pytz.UTC

    @pytest.mark.parametrize("tzstr", ["Europe/Berlin", "dateutil/Europe/Berlin"])
    def test_getitem_pydatetime_tz(self, tzstr):
        tz = timezones.maybe_get_tz(tzstr)

        index = date_range(
            start="2012-12-24 16:00", end="2012-12-24 18:00", freq="H", tz=tzstr
        )
        ts = Series(index=index, data=index.hour)
        time_pandas = Timestamp("2012-12-24 17:00", tz=tzstr)

        dt = datetime(2012, 12, 24, 17, 0)
        time_datetime = conversion.localize_pydatetime(dt, tz)
        assert ts[time_pandas] == ts[time_datetime]

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

    def test_constructor_data_aware_dtype_naive(self, tz_aware_fixture):
        # GH 25843
        tz = tz_aware_fixture
        result = Series([Timestamp("2019", tz=tz)], dtype="datetime64[ns]")
        expected = Series([Timestamp("2019")])
        tm.assert_series_equal(result, expected)
