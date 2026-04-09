import pytest

from pandas._libs.tslibs import offsets

from pandas.tseries.frequencies import (
    is_subperiod,
    is_superperiod,
)


@pytest.mark.parametrize(
    "p1,p2,expected",
    [
        # Input validation.
        (offsets.MonthEnd(), None, False),
        (offsets.YearEnd(), None, False),
        (None, offsets.YearEnd(), False),
        (None, offsets.MonthEnd(), False),
        (None, None, False),
        (offsets.YearEnd(), offsets.MonthEnd(), True),
        (offsets.Hour(), offsets.Minute(), True),
        (offsets.Second(), offsets.Milli(), True),
        (offsets.Milli(), offsets.Micro(), True),
        (offsets.Micro(), offsets.Nano(), True),
    ],
)
def test_super_sub_symmetry(p1, p2, expected):
    assert is_superperiod(p1, p2) is expected
    assert is_subperiod(p2, p1) is expected


@pytest.mark.parametrize(
    "freq",
    ["D", "B", "C", "h", "min", "s", "ms", "us", "ns", "M", "BM", "W", "Q", "BQ", "Y"],
)
def test_is_subperiod_identity(freq):
    # GH#18553
    assert is_subperiod(freq, freq) is True
    assert is_superperiod(freq, freq) is True


@pytest.mark.parametrize(
    "source,target,expected",
    [
        # Annual with same/different month anchors
        ("Y-DEC", "Y-DEC", True),
        ("Y-MAR", "Y-DEC", False),
        # Quarterly with conforming/non-conforming months
        ("Q-DEC", "Q-DEC", True),
        ("Q-MAR", "Q-DEC", True),
        ("Q-FEB", "Q-DEC", False),
    ],
)
def test_is_subperiod_anchored_identity(source, target, expected):
    # GH#18553
    assert is_subperiod(source, target) is expected
    assert is_superperiod(target, source) is expected
