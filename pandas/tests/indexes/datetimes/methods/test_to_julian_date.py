from zoneinfo import ZoneInfo

import numpy as np
import pytest

from pandas import (
    Index,
    Timestamp,
    date_range,
    to_datetime,
)
import pandas._testing as tm


class TestDateTimeIndexToJulianDate:
    def test_1700(self):
        dr = date_range(start=Timestamp("1710-10-01"), periods=5, freq="D")
        r1 = Index([x.to_julian_date() for x in dr])
        r2 = dr.to_julian_date()
        assert isinstance(r2, Index) and r2.dtype == np.float64
        tm.assert_index_equal(r1, r2)

    def test_2000(self):
        dr = date_range(start=Timestamp("2000-02-27"), periods=5, freq="D")
        r1 = Index([x.to_julian_date() for x in dr])
        r2 = dr.to_julian_date()
        assert isinstance(r2, Index) and r2.dtype == np.float64
        tm.assert_index_equal(r1, r2)

    def test_hour(self):
        dr = date_range(start=Timestamp("2000-02-27"), periods=5, freq="h")
        r1 = Index([x.to_julian_date() for x in dr])
        r2 = dr.to_julian_date()
        assert isinstance(r2, Index) and r2.dtype == np.float64
        tm.assert_index_equal(r1, r2)

    def test_minute(self):
        dr = date_range(start=Timestamp("2000-02-27"), periods=5, freq="min")
        r1 = Index([x.to_julian_date() for x in dr])
        r2 = dr.to_julian_date()
        assert isinstance(r2, Index) and r2.dtype == np.float64
        tm.assert_index_equal(r1, r2)

    def test_second(self):
        dr = date_range(start=Timestamp("2000-02-27"), periods=5, freq="s")
        r1 = Index([x.to_julian_date() for x in dr])
        r2 = dr.to_julian_date()
        assert isinstance(r2, Index) and r2.dtype == np.float64
        tm.assert_index_equal(r1, r2)

    @pytest.mark.parametrize(
        "tz, expected",
        [
            (None, 2400000.5),
            (ZoneInfo("UTC"), 2400000.5),
            (ZoneInfo("US/Pacific"), 2400000.5 + (8 / 24)),
            (ZoneInfo("Europe/London"), 2400000.5 - (1 / 24)),
        ],
    )
    def test_to_julian_date_with_timezones_single_element(self, tz, expected):
        # GH54763: Timestamp.to_julian_date() must consider timezone
        dates = to_datetime(["1858-11-17T00:00:00.0"])
        if tz:
            dates = dates.tz_localize(tz)
        result = Index(dates.to_julian_date())
        expected = Index([expected])
        tm.assert_almost_equal(result, expected, rtol=1e-6, atol=1e-6)

    @pytest.mark.parametrize(
        "tz, offset",
        [
            (None, 0),
            (ZoneInfo("UTC"), 0),
            (ZoneInfo("US/Pacific"), 8),
            (ZoneInfo("Europe/London"), -1),
        ],
    )
    def test_to_julian_date_with_timezones_multiple_elements(self, tz, offset):
        # GH54763: Timestamp.to_julian_date() must consider timezone
        dates = to_datetime(
            [
                "1858-11-17T00:00:00",
                "1858-11-17T12:00:00",
                "2000-01-01T00:00:00",
                "2000-01-01T12:00:00",
                "2000-01-01T12:00:00",
            ]
        )
        if tz:
            dates = dates.tz_localize(tz)
        result = Index(dates.to_julian_date())
        expected = Index(
            [
                2400000.5 + (offset / 24),
                2400001.0 + (offset / 24),
                2451544.5 + (offset / 24),
                2451545.0 + (offset / 24),
                2451545.0 + (offset / 24),
            ]
        )
        tm.assert_almost_equal(result, expected, rtol=1e-6, atol=1e-6)
