from pandas import Timestamp

import zoneinfo

import pandas._testing as tm
import pytest


class TestTimestampToJulianDate:
    def test_compare_1700(self):
        ts = Timestamp("1700-06-23")
        res = ts.to_julian_date()
        assert res == 2_342_145.5

    def test_compare_2000(self):
        ts = Timestamp("2000-04-12")
        res = ts.to_julian_date()
        assert res == 2_451_646.5

    def test_compare_2100(self):
        ts = Timestamp("2100-08-12")
        res = ts.to_julian_date()
        assert res == 2_488_292.5

    def test_compare_hour01(self):
        ts = Timestamp("2000-08-12T01:00:00")
        res = ts.to_julian_date()
        assert res == 2_451_768.5416666666666666

    def test_compare_hour13(self):
        ts = Timestamp("2000-08-12T13:00:00")
        res = ts.to_julian_date()
        assert res == 2_451_769.0416666666666666

    @pytest.mark.parametrize("tz, expected", [
        (None, 2400000.5),
        (zoneinfo.ZoneInfo("UTC"), 2400000.5),
        (zoneinfo.ZoneInfo("US/Pacific"), 2400000.5 + (8 / 24)),
        (zoneinfo.ZoneInfo("Europe/London"), 2400000.5 - (1 / 24))
    ])
    def test_to_julian_date_with_timezones(self, tz, expected):
        # GH54763: Timestamp.to_julian_date() must consider timezone
        ts = Timestamp('1858-11-17T00:00:00.0')
        if tz:
            ts.tz_localize(tz)
        result = ts.to_julian_date()
        tm.assert_almost_equal(result, expected, rtol=1e-6, atol=1e-6)
