"""
Tests for DatetimeIndex methods behaving like their Timestamp counterparts
"""
from datetime import datetime

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
