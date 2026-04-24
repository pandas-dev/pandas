# Arithmetic tests specific to DatetimeIndex are generally about `freq`
#  retention or inference.  Other arithmetic tests belong in
#  tests/arithmetic/test_datetime64.py

from pandas import (
    Timedelta,
    TimedeltaIndex,
    Timestamp,
    date_range,
    offsets,
    timedelta_range,
)
import pandas._testing as tm


class TestDatetimeIndexArithmetic:
    def test_add_timedelta_preserves_freq(self):
        # GH#37295 should hold for any DTI with freq=None or Tick freq
        # GH#41943 as of pandas3 Tick does not include Day
        tz = "Canada/Eastern"
        dti = date_range(
            start=Timestamp("2019-03-26 00:00:00-0400", tz=tz),
            end=Timestamp("2020-10-17 00:00:00-0400", tz=tz),
            freq="12h",
        )
        result = dti + Timedelta(days=1)
        assert result.freq == dti.freq

    def test_sub_datetime_preserves_freq(self, tz_naive_fixture):
        # GH#48818
        # In pandas3 "D" preserves time-of-day across DST transitions, so
        #  is not preserved by subtraction.  Ticks offsets like "24h"
        #  are still preserved
        dti = date_range(
            "2016-01-01",
            periods=12,
            tz=tz_naive_fixture,
            freq=offsets.Hour(24),
        )

        res = dti - dti[0]
        expected = timedelta_range("0 Days", "11 Days", freq=offsets.Hour(24))
        tm.assert_index_equal(res, expected)
        assert res.freq == expected.freq

    def test_sub_datetime_preserves_freq_across_dst(self):
        # GH#48818
        ts = Timestamp("2016-03-11", tz="US/Pacific")
        dti = date_range(ts, periods=4, unit="ns")

        res = dti - dti[0]
        expected = TimedeltaIndex(
            [
                Timedelta(days=0),
                Timedelta(days=1),
                Timedelta(days=2),
                Timedelta(days=2, hours=23),
            ],
            dtype="m8[ns]",
        )
        tm.assert_index_equal(res, expected)
        assert res.freq == expected.freq

    def test_add_dti_day(self):
        # GH#35388
        dti = date_range("2020-03-28", periods=4, freq="D", tz="Europe/Berlin")
        result = (dti + dti.freq)[:-1]
        expected = dti[1:]
        tm.assert_index_equal(result, expected)

    def test_sub_timestamp_preserves_day_freq(self):
        # GH#62094
        dti = date_range("2021-01-01", periods=5, freq="D")
        ts = Timestamp("2020-01-01")

        result = dti - ts

        # The one crucial assertion:
        assert isinstance(result, TimedeltaIndex)
        assert result.freq == dti.freq
