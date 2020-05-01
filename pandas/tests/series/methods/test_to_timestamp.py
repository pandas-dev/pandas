from datetime import timedelta

from pandas import Series, Timedelta, date_range, period_range, to_datetime
import pandas._testing as tm


class TestToTimestamp:
    def test_to_timestamp(self):
        index = period_range(freq="A", start="1/1/2001", end="12/1/2009")
        series = Series(1, index=index, name="foo")

        exp_index = date_range("1/1/2001", end="12/31/2009", freq="A-DEC")
        result = series.to_timestamp(how="end")
        exp_index = exp_index + Timedelta(1, "D") - Timedelta(1, "ns")
        tm.assert_index_equal(result.index, exp_index)
        assert result.name == "foo"

        exp_index = date_range("1/1/2001", end="1/1/2009", freq="AS-JAN")
        result = series.to_timestamp(how="start")
        tm.assert_index_equal(result.index, exp_index)

        def _get_with_delta(delta, freq="A-DEC"):
            return date_range(
                to_datetime("1/1/2001") + delta,
                to_datetime("12/31/2009") + delta,
                freq=freq,
            )

        delta = timedelta(hours=23)
        result = series.to_timestamp("H", "end")
        exp_index = _get_with_delta(delta)
        exp_index = exp_index + Timedelta(1, "h") - Timedelta(1, "ns")
        tm.assert_index_equal(result.index, exp_index)

        delta = timedelta(hours=23, minutes=59)
        result = series.to_timestamp("T", "end")
        exp_index = _get_with_delta(delta)
        exp_index = exp_index + Timedelta(1, "m") - Timedelta(1, "ns")
        tm.assert_index_equal(result.index, exp_index)

        result = series.to_timestamp("S", "end")
        delta = timedelta(hours=23, minutes=59, seconds=59)
        exp_index = _get_with_delta(delta)
        exp_index = exp_index + Timedelta(1, "s") - Timedelta(1, "ns")
        tm.assert_index_equal(result.index, exp_index)

        index = period_range(freq="H", start="1/1/2001", end="1/2/2001")
        series = Series(1, index=index, name="foo")

        exp_index = date_range("1/1/2001 00:59:59", end="1/2/2001 00:59:59", freq="H")
        result = series.to_timestamp(how="end")
        exp_index = exp_index + Timedelta(1, "s") - Timedelta(1, "ns")
        tm.assert_index_equal(result.index, exp_index)
        assert result.name == "foo"
