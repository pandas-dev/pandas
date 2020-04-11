from datetime import timedelta

from numpy import timedelta64
import pytest

from pandas import Interval, Timedelta, Timestamp


@pytest.mark.parametrize("method", ["__add__", "__sub__"])
@pytest.mark.parametrize(
    "delta", [Timedelta(days=7), timedelta(7), timedelta64(7, "D")]
)
def test_timestamp_interval_add_subtract_timedelta(method, delta):
    # https://github.com/pandas-dev/pandas/issues/32023
    interval = Interval(
        Timestamp("2017-01-01 00:00:00"), Timestamp("2018-01-01 00:00:00")
    )

    result = getattr(interval, method)(delta)

    left = getattr(interval.left, method)(delta)
    right = getattr(interval.right, method)(delta)
    expected = Interval(left, right)

    assert result == expected


@pytest.mark.parametrize(
    "delta", [Timedelta(days=7), timedelta(7)]
)
def test_timedelta_add_timestamp_interval(delta):
    # https://github.com/pandas-dev/pandas/issues/32023
    interval = Interval(
        Timestamp("2017-01-01 00:00:00"), Timestamp("2018-01-01 00:00:00")
    )

    result = delta + interval

    left = interval.left + delta
    right = interval.right + delta
    expected = Interval(left, right)

    assert result == expected


@pytest.mark.parametrize("method", ["__add__", "__sub__"])
def test_numeric_interval_add_subtract_timedelta_raises(method):
    # https://github.com/pandas-dev/pandas/issues/32023
    interval = Interval(1, 2)
    delta = Timedelta(days=7)

    msg = "unsupported operand"
    with pytest.raises(TypeError, match=msg):
        getattr(interval, method)(delta)
