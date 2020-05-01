"""Tests for PeriodIndex behaving like a vectorized Period scalar"""

from pandas import Timedelta, date_range, period_range
import pandas._testing as tm


class TestPeriodIndexOps:
    def test_start_time(self):
        # GH#17157
        index = period_range(freq="M", start="2016-01-01", end="2016-05-31")
        expected_index = date_range("2016-01-01", end="2016-05-31", freq="MS")
        tm.assert_index_equal(index.start_time, expected_index)

    def test_end_time(self):
        # GH#17157
        index = period_range(freq="M", start="2016-01-01", end="2016-05-31")
        expected_index = date_range("2016-01-01", end="2016-05-31", freq="M")
        expected_index += Timedelta(1, "D") - Timedelta(1, "ns")
        tm.assert_index_equal(index.end_time, expected_index)
