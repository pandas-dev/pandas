import pytest

from pandas._libs.tslibs import frequencies as libfrequencies, resolution
from pandas._libs.tslibs.frequencies import (
    FreqGroup,
    _period_code_map,
    get_freq,
    get_freq_code,
)

import pandas.tseries.offsets as offsets


@pytest.fixture(params=list(_period_code_map.items()))
def period_code_item(request):
    return request.param


@pytest.mark.parametrize(
    "freqstr,expected",
    [
        ("A", 1000),
        ("3A", 1000),
        ("-1A", 1000),
        ("Y", 1000),
        ("3Y", 1000),
        ("-1Y", 1000),
        ("W", 4000),
        ("W-MON", 4001),
        ("W-FRI", 4005),
    ],
)
def test_freq_code(freqstr, expected):
    assert get_freq(freqstr) == expected


def test_freq_code_match(period_code_item):
    freqstr, code = period_code_item
    assert get_freq(freqstr) == code


@pytest.mark.parametrize(
    "freqstr,expected",
    [
        ("A", 1000),
        ("3A", 1000),
        ("-1A", 1000),
        ("A-JAN", 1000),
        ("A-MAY", 1000),
        ("Y", 1000),
        ("3Y", 1000),
        ("-1Y", 1000),
        ("Y-JAN", 1000),
        ("Y-MAY", 1000),
        (offsets.YearEnd(), 1000),
        (offsets.YearEnd(month=1), 1000),
        (offsets.YearEnd(month=5), 1000),
        ("W", 4000),
        ("W-MON", 4000),
        ("W-FRI", 4000),
        (offsets.Week(), 4000),
        (offsets.Week(weekday=1), 4000),
        (offsets.Week(weekday=5), 4000),
        ("T", FreqGroup.FR_MIN),
    ],
)
def test_freq_group(freqstr, expected):
    assert resolution.get_freq_group(freqstr) == expected


def test_freq_group_match(period_code_item):
    freqstr, code = period_code_item

    str_group = resolution.get_freq_group(freqstr)
    code_group = resolution.get_freq_group(code)

    assert str_group == code_group == code // 1000 * 1000


@pytest.mark.parametrize(
    "freqstr,exp_freqstr",
    [("D", "D"), ("W", "D"), ("M", "D"), ("S", "S"), ("T", "S"), ("H", "S")],
)
def test_get_to_timestamp_base(freqstr, exp_freqstr):
    tsb = libfrequencies.get_to_timestamp_base

    assert tsb(get_freq_code(freqstr)[0]) == get_freq_code(exp_freqstr)[0]


_reso = resolution.Resolution


@pytest.mark.parametrize(
    "freqstr,expected",
    [
        ("A", "year"),
        ("Q", "quarter"),
        ("M", "month"),
        ("D", "day"),
        ("H", "hour"),
        ("T", "minute"),
        ("S", "second"),
        ("L", "millisecond"),
        ("U", "microsecond"),
        ("N", "nanosecond"),
    ],
)
def test_get_str_from_freq(freqstr, expected):
    assert _reso.get_str_from_freq(freqstr) == expected


@pytest.mark.parametrize("freq", ["A", "Q", "M", "D", "H", "T", "S", "L", "U", "N"])
def test_get_freq_roundtrip(freq):
    result = _reso.get_freq(_reso.get_str_from_freq(freq))
    assert freq == result


@pytest.mark.parametrize("freq", ["D", "H", "T", "S", "L", "U"])
def test_get_freq_roundtrip2(freq):
    result = _reso.get_freq(_reso.get_str(_reso.get_reso_from_freq(freq)))
    assert freq == result


@pytest.mark.parametrize(
    "args,expected",
    [
        ((1.5, "T"), (90, "S")),
        ((62.4, "T"), (3744, "S")),
        ((1.04, "H"), (3744, "S")),
        ((1, "D"), (1, "D")),
        ((0.342931, "H"), (1234551600, "U")),
        ((1.2345, "D"), (106660800, "L")),
    ],
)
def test_resolution_bumping(args, expected):
    # see gh-14378
    assert _reso.get_stride_from_decimal(*args) == expected


@pytest.mark.parametrize(
    "args",
    [
        (0.5, "N"),
        # Too much precision in the input can prevent.
        (0.3429324798798269273987982, "H"),
    ],
)
def test_cat(args):
    msg = "Could not convert to integer offset at any resolution"

    with pytest.raises(ValueError, match=msg):
        _reso.get_stride_from_decimal(*args)


@pytest.mark.parametrize(
    "freq_input,expected",
    [
        # Frequency string.
        ("A", (get_freq("A"), 1)),
        ("3D", (get_freq("D"), 3)),
        ("-2M", (get_freq("M"), -2)),
        # Tuple.
        (("D", 1), (get_freq("D"), 1)),
        (("A", 3), (get_freq("A"), 3)),
        (("M", -2), (get_freq("M"), -2)),
        ((5, "T"), (FreqGroup.FR_MIN, 5)),
        # Numeric Tuple.
        ((1000, 1), (1000, 1)),
        # Offsets.
        (offsets.Day(), (get_freq("D"), 1)),
        (offsets.Day(3), (get_freq("D"), 3)),
        (offsets.Day(-2), (get_freq("D"), -2)),
        (offsets.MonthEnd(), (get_freq("M"), 1)),
        (offsets.MonthEnd(3), (get_freq("M"), 3)),
        (offsets.MonthEnd(-2), (get_freq("M"), -2)),
        (offsets.Week(), (get_freq("W"), 1)),
        (offsets.Week(3), (get_freq("W"), 3)),
        (offsets.Week(-2), (get_freq("W"), -2)),
        (offsets.Hour(), (FreqGroup.FR_HR, 1)),
        # Monday is weekday=0.
        (offsets.Week(weekday=1), (get_freq("W-TUE"), 1)),
        (offsets.Week(3, weekday=0), (get_freq("W-MON"), 3)),
        (offsets.Week(-2, weekday=4), (get_freq("W-FRI"), -2)),
    ],
)
def test_get_freq_code(freq_input, expected):
    assert get_freq_code(freq_input) == expected


def test_get_code_invalid():
    with pytest.raises(ValueError, match="Invalid frequency"):
        get_freq_code((5, "baz"))
