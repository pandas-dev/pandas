"""
Tests for DatetimeIndex methods behaving like their Timestamp counterparts
"""
from datetime import datetime

import numpy as np
import pytest

from pandas._libs.tslibs import OutOfBoundsDatetime

import pandas as pd
from pandas import (
    DatetimeIndex,
    Timestamp,
    date_range,
)
import pandas._testing as tm


class TestDatetimeIndexOps:
    def test_dti_time(self):
        rng = date_range("1/1/2000", freq="12min", periods=10)
        result = pd.Index(rng).time
        expected = [t.time() for t in rng]
        assert (result == expected).all()

    def test_dti_date(self):
        rng = date_range("1/1/2000", freq="12h", periods=10)
        result = pd.Index(rng).date
        expected = [t.date() for t in rng]
        assert (result == expected).all()

    @pytest.mark.parametrize("data", [["1400-01-01"], [datetime(1400, 1, 1)]])
    def test_dti_date_out_of_range(self, data):
        # GH#1475
        msg = (
            "^Out of bounds nanosecond timestamp: "
            "1400-01-01( 00:00:00)?, at position 0$"
        )
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            DatetimeIndex(data)

    @pytest.mark.parametrize(
        "field",
        [
            "dayofweek",
            "day_of_week",
            "dayofyear",
            "day_of_year",
            "quarter",
            "days_in_month",
            "is_month_start",
            "is_month_end",
            "is_quarter_start",
            "is_quarter_end",
            "is_year_start",
            "is_year_end",
        ],
    )
    def test_dti_timestamp_fields(self, field):
        # extra fields from DatetimeIndex like quarter and week
        idx = tm.makeDateIndex(100)
        expected = getattr(idx, field)[-1]

        result = getattr(Timestamp(idx[-1]), field)
        assert result == expected

    def test_dti_timestamp_isocalendar_fields(self):
        idx = tm.makeDateIndex(100)
        expected = tuple(idx.isocalendar().iloc[-1].to_list())
        result = idx[-1].isocalendar()
        assert result == expected

    # ----------------------------------------------------------------
    # DatetimeIndex.normalize

    def test_normalize(self):
        rng = date_range("1/1/2000 9:30", periods=10, freq="D")

        result = rng.normalize()
        expected = date_range("1/1/2000", periods=10, freq="D")
        tm.assert_index_equal(result, expected)

        arr_ns = np.array([1380585623454345752, 1380585612343234312]).astype(
            "datetime64[ns]"
        )
        rng_ns = DatetimeIndex(arr_ns)
        rng_ns_normalized = rng_ns.normalize()

        arr_ns = np.array([1380585600000000000, 1380585600000000000]).astype(
            "datetime64[ns]"
        )
        expected = DatetimeIndex(arr_ns)
        tm.assert_index_equal(rng_ns_normalized, expected)

        assert result.is_normalized
        assert not rng.is_normalized

    def test_normalize_nat(self):
        dti = DatetimeIndex([pd.NaT, Timestamp("2018-01-01 01:00:00")])
        result = dti.normalize()
        expected = DatetimeIndex([pd.NaT, Timestamp("2018-01-01")])
        tm.assert_index_equal(result, expected)


class TestDateTimeIndexToJulianDate:
    def test_1700(self):
        dr = date_range(start=Timestamp("1710-10-01"), periods=5, freq="D")
        r1 = pd.Index([x.to_julian_date() for x in dr])
        r2 = dr.to_julian_date()
        assert isinstance(r2, pd.Index) and r2.dtype == np.float64
        tm.assert_index_equal(r1, r2)

    def test_2000(self):
        dr = date_range(start=Timestamp("2000-02-27"), periods=5, freq="D")
        r1 = pd.Index([x.to_julian_date() for x in dr])
        r2 = dr.to_julian_date()
        assert isinstance(r2, pd.Index) and r2.dtype == np.float64
        tm.assert_index_equal(r1, r2)

    def test_hour(self):
        dr = date_range(start=Timestamp("2000-02-27"), periods=5, freq="h")
        r1 = pd.Index([x.to_julian_date() for x in dr])
        r2 = dr.to_julian_date()
        assert isinstance(r2, pd.Index) and r2.dtype == np.float64
        tm.assert_index_equal(r1, r2)

    def test_minute(self):
        dr = date_range(start=Timestamp("2000-02-27"), periods=5, freq="min")
        r1 = pd.Index([x.to_julian_date() for x in dr])
        r2 = dr.to_julian_date()
        assert isinstance(r2, pd.Index) and r2.dtype == np.float64
        tm.assert_index_equal(r1, r2)

    def test_second(self):
        dr = date_range(start=Timestamp("2000-02-27"), periods=5, freq="s")
        r1 = pd.Index([x.to_julian_date() for x in dr])
        r2 = dr.to_julian_date()
        assert isinstance(r2, pd.Index) and r2.dtype == np.float64
        tm.assert_index_equal(r1, r2)
