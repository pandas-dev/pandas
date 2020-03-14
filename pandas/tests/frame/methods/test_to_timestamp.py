from datetime import timedelta

import numpy as np
import pytest

from pandas import (
    DataFrame,
    DatetimeIndex,
    Timedelta,
    date_range,
    period_range,
    to_datetime,
)
import pandas._testing as tm


class TestToTimestamp:
    def test_frame_to_time_stamp(self):
        K = 5
        index = period_range(freq="A", start="1/1/2001", end="12/1/2009")
        df = DataFrame(np.random.randn(len(index), K), index=index)
        df["mix"] = "a"

        exp_index = date_range("1/1/2001", end="12/31/2009", freq="A-DEC")
        exp_index = exp_index + Timedelta(1, "D") - Timedelta(1, "ns")
        result = df.to_timestamp("D", "end")
        tm.assert_index_equal(result.index, exp_index)
        tm.assert_numpy_array_equal(result.values, df.values)

        exp_index = date_range("1/1/2001", end="1/1/2009", freq="AS-JAN")
        result = df.to_timestamp("D", "start")
        tm.assert_index_equal(result.index, exp_index)

        def _get_with_delta(delta, freq="A-DEC"):
            return date_range(
                to_datetime("1/1/2001") + delta,
                to_datetime("12/31/2009") + delta,
                freq=freq,
            )

        delta = timedelta(hours=23)
        result = df.to_timestamp("H", "end")
        exp_index = _get_with_delta(delta)
        exp_index = exp_index + Timedelta(1, "h") - Timedelta(1, "ns")
        tm.assert_index_equal(result.index, exp_index)

        delta = timedelta(hours=23, minutes=59)
        result = df.to_timestamp("T", "end")
        exp_index = _get_with_delta(delta)
        exp_index = exp_index + Timedelta(1, "m") - Timedelta(1, "ns")
        tm.assert_index_equal(result.index, exp_index)

        result = df.to_timestamp("S", "end")
        delta = timedelta(hours=23, minutes=59, seconds=59)
        exp_index = _get_with_delta(delta)
        exp_index = exp_index + Timedelta(1, "s") - Timedelta(1, "ns")
        tm.assert_index_equal(result.index, exp_index)

        # columns
        df = df.T

        exp_index = date_range("1/1/2001", end="12/31/2009", freq="A-DEC")
        exp_index = exp_index + Timedelta(1, "D") - Timedelta(1, "ns")
        result = df.to_timestamp("D", "end", axis=1)
        tm.assert_index_equal(result.columns, exp_index)
        tm.assert_numpy_array_equal(result.values, df.values)

        exp_index = date_range("1/1/2001", end="1/1/2009", freq="AS-JAN")
        result = df.to_timestamp("D", "start", axis=1)
        tm.assert_index_equal(result.columns, exp_index)

        delta = timedelta(hours=23)
        result = df.to_timestamp("H", "end", axis=1)
        exp_index = _get_with_delta(delta)
        exp_index = exp_index + Timedelta(1, "h") - Timedelta(1, "ns")
        tm.assert_index_equal(result.columns, exp_index)

        delta = timedelta(hours=23, minutes=59)
        result = df.to_timestamp("T", "end", axis=1)
        exp_index = _get_with_delta(delta)
        exp_index = exp_index + Timedelta(1, "m") - Timedelta(1, "ns")
        tm.assert_index_equal(result.columns, exp_index)

        result = df.to_timestamp("S", "end", axis=1)
        delta = timedelta(hours=23, minutes=59, seconds=59)
        exp_index = _get_with_delta(delta)
        exp_index = exp_index + Timedelta(1, "s") - Timedelta(1, "ns")
        tm.assert_index_equal(result.columns, exp_index)

        # invalid axis
        with pytest.raises(ValueError, match="axis"):
            df.to_timestamp(axis=2)

        result1 = df.to_timestamp("5t", axis=1)
        result2 = df.to_timestamp("t", axis=1)
        expected = date_range("2001-01-01", "2009-01-01", freq="AS")
        assert isinstance(result1.columns, DatetimeIndex)
        assert isinstance(result2.columns, DatetimeIndex)
        tm.assert_numpy_array_equal(result1.columns.asi8, expected.asi8)
        tm.assert_numpy_array_equal(result2.columns.asi8, expected.asi8)
        # PeriodIndex.to_timestamp always use 'infer'
        assert result1.columns.freqstr == "AS-JAN"
        assert result2.columns.freqstr == "AS-JAN"
