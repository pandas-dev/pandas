"""Tests for PeriodIndex behaving like a vectorized Period scalar"""

import numpy as np
import pytest

from pandas import (
    date_range,
    period_range,
)
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
        expected_index = date_range("2016-01-01", end="2016-05-31", freq="ME")
        expected_index += np.timedelta64(1, "D") - np.timedelta64(1, "us")
        tm.assert_index_equal(index.end_time, expected_index)

    @pytest.mark.filterwarnings(r"ignore:PeriodDtype\[B\] is deprecated:FutureWarning")
    @pytest.mark.filterwarnings(
        "ignore:Period with BDay freq is deprecated:FutureWarning"
    )
    def test_end_time_business_friday(self):
        # GH#34449
        pi = period_range("1990-01-05", freq="B", periods=1)
        result = pi.end_time

        dti = date_range("1990-01-05", freq="D", periods=1)._with_freq(None)
        expected = dti + np.timedelta64(1, "D") - np.timedelta64(1, "us")
        tm.assert_index_equal(result, expected)
