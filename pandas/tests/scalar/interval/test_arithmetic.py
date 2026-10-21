from datetime import timedelta

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestIntervalArithmetic:
    def test_interval_add(self, closed):
        interval = pd.Interval(0, 1, closed=closed)
        expected = pd.Interval(1, 2, closed=closed)

        result = interval + 1
        assert result == expected

        result = 1 + interval
        assert result == expected

        result = interval
        result += 1
        assert result == expected

        msg = r"unsupported operand type\(s\) for \+"
        with pytest.raises(TypeError, match=msg):
            interval + interval

        with pytest.raises(TypeError, match=msg):
            interval + "foo"

    def test_interval_sub(self, closed):
        interval = pd.Interval(0, 1, closed=closed)
        expected = pd.Interval(-1, 0, closed=closed)

        result = interval - 1
        assert result == expected

        result = interval
        result -= 1
        assert result == expected

        msg = r"unsupported operand type\(s\) for -"
        with pytest.raises(TypeError, match=msg):
            interval - interval

        with pytest.raises(TypeError, match=msg):
            interval - "foo"

    def test_interval_mult(self, closed):
        interval = pd.Interval(0, 1, closed=closed)
        expected = pd.Interval(0, 2, closed=closed)

        result = interval * 2
        assert result == expected

        result = 2 * interval
        assert result == expected

        result = interval
        result *= 2
        assert result == expected

        msg = r"unsupported operand type\(s\) for \*"
        with pytest.raises(TypeError, match=msg):
            interval * interval

        msg = r"can\'t multiply sequence by non-int"
        with pytest.raises(TypeError, match=msg):
            interval * "foo"

    def test_interval_div(self, closed):
        interval = pd.Interval(0, 1, closed=closed)
        expected = pd.Interval(0, 0.5, closed=closed)

        result = interval / 2.0
        assert result == expected

        result = interval
        result /= 2.0
        assert result == expected

        msg = r"unsupported operand type\(s\) for /"
        with pytest.raises(TypeError, match=msg):
            interval / interval

        with pytest.raises(TypeError, match=msg):
            interval / "foo"

    def test_interval_floordiv(self, closed):
        interval = pd.Interval(1, 2, closed=closed)
        expected = pd.Interval(0, 1, closed=closed)

        result = interval // 2
        assert result == expected

        result = interval
        result //= 2
        assert result == expected

        msg = r"unsupported operand type\(s\) for //"
        with pytest.raises(TypeError, match=msg):
            interval // interval

        with pytest.raises(TypeError, match=msg):
            interval // "foo"

    @pytest.mark.parametrize("method", ["__add__", "__sub__"])
    @pytest.mark.parametrize(
        "interval",
        [
            pd.Interval(
                pd.Timestamp("2017-01-01 00:00:00"), pd.Timestamp("2018-01-01 00:00:00")
            ),
            pd.Interval(pd.Timedelta(days=7), pd.Timedelta(days=14)),
        ],
    )
    @pytest.mark.parametrize(
        "delta", [pd.Timedelta(days=7), timedelta(7), np.timedelta64(7, "D")]
    )
    def test_time_interval_add_subtract_timedelta(self, interval, delta, method):
        # https://github.com/pandas-dev/pandas/issues/32023
        result = getattr(interval, method)(delta)
        left = getattr(interval.left, method)(delta)
        right = getattr(interval.right, method)(delta)
        expected = pd.Interval(left, right)

        assert result == expected

    @pytest.mark.parametrize("interval", [pd.Interval(1, 2), pd.Interval(1.0, 2.0)])
    @pytest.mark.parametrize(
        "delta", [pd.Timedelta(days=7), timedelta(7), np.timedelta64(7, "D")]
    )
    @pytest.mark.filterwarnings(
        "ignore:.*'generic' unit for NumPy timedelta:DeprecationWarning"
    )
    def test_numeric_interval_add_timedelta_raises(self, interval, delta):
        # https://github.com/pandas-dev/pandas/issues/32023
        msg = "|".join(
            [
                "unsupported operand",
                "cannot use operands",
                "Only numeric, Timestamp and Timedelta endpoints are allowed",
            ]
        )
        with pytest.raises((TypeError, ValueError), match=msg):
            interval + delta

        with pytest.raises((TypeError, ValueError), match=msg):
            delta + interval

    @pytest.mark.parametrize(
        "delta",
        [timedelta(0), np.timedelta64(0, "ns"), pd.Timedelta(0)],
    )
    def test_timedelta_add_timestamp_interval(self, delta):
        expected = pd.Interval(pd.Timestamp("2020-01-01"), pd.Timestamp("2020-02-01"))

        result = delta + expected
        assert result == expected

        result = expected + delta
        assert result == expected


class TestIntervalComparisons:
    def test_interval_equal(self):
        assert pd.Interval(0, 1) == pd.Interval(0, 1, closed="right")
        assert pd.Interval(0, 1) != pd.Interval(0, 1, closed="left")
        assert pd.Interval(0, 1) != 0

    def test_interval_comparison(self):
        msg = (
            "'<' not supported between instances of "
            "'pandas._libs.interval.Interval' and 'int'"
        )
        with pytest.raises(TypeError, match=msg):
            pd.Interval(0, 1) < 2

        assert pd.Interval(0, 1) < pd.Interval(1, 2)
        assert pd.Interval(0, 1) < pd.Interval(0, 2)
        assert pd.Interval(0, 1) < pd.Interval(0.5, 1.5)
        assert pd.Interval(0, 1) <= pd.Interval(0, 1)
        assert pd.Interval(0, 1) > pd.Interval(-1, 2)
        assert pd.Interval(0, 1) >= pd.Interval(0, 1)

    def test_equality_comparison_broadcasts_over_array(self):
        # https://github.com/pandas-dev/pandas/issues/35931
        interval = pd.Interval(0, 1)
        arr = np.array([interval, interval])
        result = interval == arr
        expected = np.array([True, True])
        tm.assert_numpy_array_equal(result, expected)
