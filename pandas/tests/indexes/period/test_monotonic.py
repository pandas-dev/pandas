import pandas as pd


def test_is_monotonic_increasing():
    # GH#17717
    p0 = pd.Period("2017-09-01")
    p1 = pd.Period("2017-09-02")
    p2 = pd.Period("2017-09-03")

    idx_inc0 = pd.PeriodIndex([p0, p1, p2])
    idx_inc1 = pd.PeriodIndex([p0, p1, p1])
    idx_dec0 = pd.PeriodIndex([p2, p1, p0])
    idx_dec1 = pd.PeriodIndex([p2, p1, p1])
    idx = pd.PeriodIndex([p1, p2, p0])

    assert idx_inc0.is_monotonic_increasing is True
    assert idx_inc1.is_monotonic_increasing is True
    assert idx_dec0.is_monotonic_increasing is False
    assert idx_dec1.is_monotonic_increasing is False
    assert idx.is_monotonic_increasing is False


def test_is_monotonic_decreasing():
    # GH#17717
    p0 = pd.Period("2017-09-01")
    p1 = pd.Period("2017-09-02")
    p2 = pd.Period("2017-09-03")

    idx_inc0 = pd.PeriodIndex([p0, p1, p2])
    idx_inc1 = pd.PeriodIndex([p0, p1, p1])
    idx_dec0 = pd.PeriodIndex([p2, p1, p0])
    idx_dec1 = pd.PeriodIndex([p2, p1, p1])
    idx = pd.PeriodIndex([p1, p2, p0])

    assert idx_inc0.is_monotonic_decreasing is False
    assert idx_inc1.is_monotonic_decreasing is False
    assert idx_dec0.is_monotonic_decreasing is True
    assert idx_dec1.is_monotonic_decreasing is True
    assert idx.is_monotonic_decreasing is False
