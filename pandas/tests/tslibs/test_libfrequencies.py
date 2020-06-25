import pytest

from pandas._libs.tslibs.parsing import get_rule_month

from pandas.tseries import offsets


@pytest.mark.parametrize(
    "obj,expected",
    [
        ("W", "DEC"),
        (offsets.Week(), "DEC"),
        ("D", "DEC"),
        (offsets.Day(), "DEC"),
        ("Q", "DEC"),
        (offsets.QuarterEnd(startingMonth=12), "DEC"),
        ("Q-JAN", "JAN"),
        (offsets.QuarterEnd(startingMonth=1), "JAN"),
        ("A-DEC", "DEC"),
        ("Y-DEC", "DEC"),
        (offsets.YearEnd(), "DEC"),
        ("A-MAY", "MAY"),
        ("Y-MAY", "MAY"),
        (offsets.YearEnd(month=5), "MAY"),
    ],
)
def test_get_rule_month(obj, expected):
    result = get_rule_month(obj)
    assert result == expected
