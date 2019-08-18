import pytest

from pandas._libs.tslibs.frequencies import (
    INVALID_FREQ_ERR_MSG,
    _period_str_to_code,
    get_rule_month,
    is_subperiod,
    is_superperiod,
)

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


@pytest.mark.parametrize(
    "obj,expected",
    [
        ("A", 1000),
        ("A-DEC", 1000),
        ("A-JAN", 1001),
        ("Y", 1000),
        ("Y-DEC", 1000),
        ("Y-JAN", 1001),
        ("Q", 2000),
        ("Q-DEC", 2000),
        ("Q-FEB", 2002),
        ("W", 4000),
        ("W-SUN", 4000),
        ("W-FRI", 4005),
        ("Min", 8000),
        ("ms", 10000),
        ("US", 11000),
        ("NS", 12000),
    ],
)
def test_period_str_to_code(obj, expected):
    assert _period_str_to_code(obj) == expected


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
    "freq,expected,aliases",
    [
        ("D", 6000, ["DAY", "DLY", "DAILY"]),
        ("M", 3000, ["MTH", "MONTH", "MONTHLY"]),
        ("N", 12000, ["NANOSECOND", "NANOSECONDLY"]),
        ("H", 7000, ["HR", "HOUR", "HRLY", "HOURLY"]),
        ("T", 8000, ["minute", "MINUTE", "MINUTELY"]),
        ("L", 10000, ["MILLISECOND", "MILLISECONDLY"]),
        ("U", 11000, ["MICROSECOND", "MICROSECONDLY"]),
        ("S", 9000, ["sec", "SEC", "SECOND", "SECONDLY"]),
        ("B", 5000, ["BUS", "BUSINESS", "BUSINESSLY", "WEEKDAY"]),
    ],
)
def test_assert_aliases_deprecated(freq, expected, aliases):
    assert isinstance(aliases, list)
    assert _period_str_to_code(freq) == expected

    for alias in aliases:
        with pytest.raises(ValueError, match=INVALID_FREQ_ERR_MSG):
            _period_str_to_code(alias)
