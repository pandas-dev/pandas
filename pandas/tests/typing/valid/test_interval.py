# flake8: noqa: F841
# pyright: reportGeneralTypeIssues = true

import pandas as pd


def test_interval_init() -> None:
    i1: pd.Interval = pd.Interval(1, 2, closed="both")
    i2: pd.Interval = pd.Interval(1, right=2, closed="right")
    i3: pd.Interval = pd.Interval(left=1, right=2, closed="left")


def test_interval_arithmetic() -> None:
    i1: pd.Interval = pd.Interval(1, 2, closed="both")
    i2: pd.Interval = pd.Interval(1, right=2, closed="right")

    i3: pd.Interval = i1 + 1
    i4: pd.Interval = i1 - 1
    i5: pd.Interval = i1 * 2
    i6: pd.Interval = i1 / 2
    i7: pd.Interval = i1 // 2


def test_max_intervals() -> None:
    i1 = pd.Interval(
        pd.Timestamp("2000-01-01"), pd.Timestamp("2000-01-02"), closed="both"
    )
    i2 = pd.Interval(
        pd.Timestamp("2000-01-01T12:00:00"), pd.Timestamp("2000-01-02"), closed="both"
    )
    print(max(i1.left, i2.left))


def test_interval_length() -> None:
    i1 = pd.Interval(
        pd.Timestamp("2000-01-01"), pd.Timestamp("2000-01-02"), closed="both"
    )
    i1.length.total_seconds()
