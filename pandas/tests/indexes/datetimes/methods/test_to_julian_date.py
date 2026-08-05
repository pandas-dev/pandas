import numpy as np
import pytest

from pandas import (
    Index,
    Timestamp,
    date_range,
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

    @pytest.mark.parametrize("tz", ["UTC", "US/Pacific", "Europe/London"])
    def test_tz_aware_uses_utc_instant(self, tz):
        # GH#54763 to_julian_date should reflect the UTC instant, not wall clock
        dr = date_range(start="2000-02-27", periods=5, freq="h", tz=tz)
        result = dr.to_julian_date()
        expected = Index([ts.to_julian_date() for ts in dr])
        tm.assert_index_equal(result, expected)
        # also matches converting to naive UTC up front
        tm.assert_index_equal(
            result, dr.tz_convert("UTC").tz_localize(None).to_julian_date()
        )
