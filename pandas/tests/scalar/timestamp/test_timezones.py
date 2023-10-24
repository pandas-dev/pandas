"""
Tests for Timestamp timezone-related methods
"""
from datetime import (
    date,
    datetime,
    timezone,
)

from dateutil.tz import (
    gettz,
    tzoffset,
)
import pytest
import pytz

from pandas._libs.tslibs import timezones

from pandas import Timestamp

try:
    from zoneinfo import ZoneInfo
except ImportError:
    # Cannot assign to a type
    ZoneInfo = None  # type: ignore[misc, assignment]


class TestTimestampTZOperations:
    # ------------------------------------------------------------------
    # Timestamp.__init__ with tz str or tzinfo

    def test_timestamp_constructor_tz_utc(self):
        utc_stamp = Timestamp("3/11/2012 05:00", tz="utc")
        assert utc_stamp.tzinfo is timezone.utc
        assert utc_stamp.hour == 5

        utc_stamp = Timestamp("3/11/2012 05:00").tz_localize("utc")
        assert utc_stamp.hour == 5

    def test_timestamp_to_datetime_tzoffset(self):
        tzinfo = tzoffset(None, 7200)
        expected = Timestamp("3/11/2012 04:00", tz=tzinfo)
        result = Timestamp(expected.to_pydatetime())
        assert expected == result

    def test_timestamp_constructor_near_dst_boundary(self):
        # GH#11481 & GH#15777
        # Naive string timestamps were being localized incorrectly
        # with tz_convert_from_utc_single instead of tz_localize_to_utc

        for tz in ["Europe/Brussels", "Europe/Prague"]:
            result = Timestamp("2015-10-25 01:00", tz=tz)
            expected = Timestamp("2015-10-25 01:00").tz_localize(tz)
            assert result == expected

            msg = "Cannot infer dst time from 2015-10-25 02:00:00"
            with pytest.raises(pytz.AmbiguousTimeError, match=msg):
                Timestamp("2015-10-25 02:00", tz=tz)

        result = Timestamp("2017-03-26 01:00", tz="Europe/Paris")
        expected = Timestamp("2017-03-26 01:00").tz_localize("Europe/Paris")
        assert result == expected

        msg = "2017-03-26 02:00"
        with pytest.raises(pytz.NonExistentTimeError, match=msg):
            Timestamp("2017-03-26 02:00", tz="Europe/Paris")

        # GH#11708
        naive = Timestamp("2015-11-18 10:00:00")
        result = naive.tz_localize("UTC").tz_convert("Asia/Kolkata")
        expected = Timestamp("2015-11-18 15:30:00+0530", tz="Asia/Kolkata")
        assert result == expected

        # GH#15823
        result = Timestamp("2017-03-26 00:00", tz="Europe/Paris")
        expected = Timestamp("2017-03-26 00:00:00+0100", tz="Europe/Paris")
        assert result == expected

        result = Timestamp("2017-03-26 01:00", tz="Europe/Paris")
        expected = Timestamp("2017-03-26 01:00:00+0100", tz="Europe/Paris")
        assert result == expected

        msg = "2017-03-26 02:00"
        with pytest.raises(pytz.NonExistentTimeError, match=msg):
            Timestamp("2017-03-26 02:00", tz="Europe/Paris")

        result = Timestamp("2017-03-26 02:00:00+0100", tz="Europe/Paris")
        naive = Timestamp(result.as_unit("ns")._value)
        expected = naive.tz_localize("UTC").tz_convert("Europe/Paris")
        assert result == expected

        result = Timestamp("2017-03-26 03:00", tz="Europe/Paris")
        expected = Timestamp("2017-03-26 03:00:00+0200", tz="Europe/Paris")
        assert result == expected

    @pytest.mark.parametrize(
        "tz",
        [
            pytz.timezone("US/Eastern"),
            gettz("US/Eastern"),
            "US/Eastern",
            "dateutil/US/Eastern",
        ],
    )
    def test_timestamp_constructed_by_date_and_tz(self, tz):
        # GH#2993, Timestamp cannot be constructed by datetime.date
        # and tz correctly

        result = Timestamp(date(2012, 3, 11), tz=tz)

        expected = Timestamp("3/11/2012", tz=tz)
        assert result.hour == expected.hour
        assert result == expected

    def test_timestamp_timetz_equivalent_with_datetime_tz(self, tz_naive_fixture):
        # GH21358
        tz = timezones.maybe_get_tz(tz_naive_fixture)

        stamp = Timestamp("2018-06-04 10:20:30", tz=tz)
        _datetime = datetime(2018, 6, 4, hour=10, minute=20, second=30, tzinfo=tz)

        result = stamp.timetz()
        expected = _datetime.timetz()

        assert result == expected
