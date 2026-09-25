import numpy as np
import pytest

import pandas as pd


@pytest.fixture
def interval():
    return pd.Interval(0, 1)


class TestInterval:
    def test_properties(self, interval):
        assert interval.closed == "right"
        assert interval.left == 0
        assert interval.right == 1
        assert interval.mid == 0.5

    def test_hash(self, interval):
        # should not raise
        hash(interval)

    @pytest.mark.parametrize(
        "left, right, expected",
        [
            (0, 5, 5),
            (-2, 5.5, 7.5),
            (10, 10, 0),
            (10, np.inf, np.inf),
            (-np.inf, -5, np.inf),
            (-np.inf, np.inf, np.inf),
            (pd.Timedelta("0 days"), pd.Timedelta("5 days"), pd.Timedelta("5 days")),
            (pd.Timedelta("10 days"), pd.Timedelta("10 days"), pd.Timedelta("0 days")),
            (pd.Timedelta("1h10min"), pd.Timedelta("5h5min"), pd.Timedelta("3h55min")),
            (pd.Timedelta("5s"), pd.Timedelta("1h"), pd.Timedelta("59min55s")),
        ],
    )
    def test_length(self, left, right, expected):
        # GH 18789
        iv = pd.Interval(left, right)
        result = iv.length
        assert result == expected

    @pytest.mark.parametrize(
        "left, right, expected",
        [
            ("2017-01-01", "2017-01-06", "5 days"),
            ("2017-01-01", "2017-01-01 12:00:00", "12 hours"),
            ("2017-01-01 12:00", "2017-01-01 12:00:00", "0 days"),
            ("2017-01-01 12:01", "2017-01-05 17:31:00", "4 days 5 hours 30 min"),
        ],
    )
    @pytest.mark.parametrize("tz", (None, "UTC", "CET", "US/Eastern"))
    def test_length_timestamp(self, tz, left, right, expected):
        # GH 18789
        iv = pd.Interval(pd.Timestamp(left, tz=tz), pd.Timestamp(right, tz=tz))
        result = iv.length
        expected = pd.Timedelta(expected)
        assert result == expected

    @pytest.mark.parametrize(
        "left, right",
        [
            (0, 1),
            (pd.Timedelta("0 days"), pd.Timedelta("1 day")),
            (pd.Timestamp("2018-01-01"), pd.Timestamp("2018-01-02")),
            (
                pd.Timestamp("2018-01-01", tz="US/Eastern"),
                pd.Timestamp("2018-01-02", tz="US/Eastern"),
            ),
        ],
    )
    def test_is_empty(self, left, right, closed):
        # GH27219
        # non-empty always return False
        iv = pd.Interval(left, right, closed)
        assert iv.is_empty is False

        # same endpoint is empty except when closed='both' (contains one point)
        iv = pd.Interval(left, left, closed)
        result = iv.is_empty
        expected = closed != "both"
        assert result is expected
