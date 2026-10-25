import pandas as pd


def test_interval_repr():
    interval = pd.Interval(0, 1)
    assert repr(interval) == "Interval(0, 1, closed='right')"
    assert str(interval) == "(0, 1]"

    interval_left = pd.Interval(0, 1, closed="left")
    assert repr(interval_left) == "Interval(0, 1, closed='left')"
    assert str(interval_left) == "[0, 1)"
