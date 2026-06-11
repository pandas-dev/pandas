import pytest

from pandas import Timestamp


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

    @pytest.mark.parametrize(
        "tz, utc_hour",
        [
            ("UTC", 13),
            ("US/Pacific", 20),  # PDT, UTC-7
            ("Europe/London", 12),  # BST, UTC+1
        ],
    )
    def test_tz_aware_uses_utc_instant(self, tz, utc_hour):
        # GH#54763 to_julian_date should reflect the UTC instant, not wall clock
        ts = Timestamp("2000-08-12T13:00:00", tz=tz)
        # JD 2451769.0 == 2000-08-12 12:00 UTC (JD epoch is noon)
        expected = 2_451_769.0 + (utc_hour - 12) / 24
        assert ts.to_julian_date() == expected

    def test_tz_aware_matches_utc_conversion(self):
        # GH#54763 result must agree with manually converting to UTC first
        ts = Timestamp("2020-03-14T15:32:52", tz="US/Eastern")
        expected = ts.tz_convert("UTC").tz_localize(None).to_julian_date()
        assert ts.to_julian_date() == expected
