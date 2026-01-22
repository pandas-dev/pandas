import pytest

from pandas import Timedelta

from pandas.tseries.frequencies import to_offset


def test_to_offset_timedelta_string_fallback():
    # This is the case you found: str(Timedelta)
    # Previously, this might have failed in the regex loop
    freq_str = "0 days 01:30:00"
    result = to_offset(freq_str)

    expected = to_offset("1h30min")
    assert result == expected


@pytest.mark.parametrize(
    "freq",
    [
        "1 days 02:00:00",
        "00:00:00.001",
        "1h 30min",  # spaces between components
    ],
)
def test_to_offset_consistent_with_timedelta(freq):
    # Ensure our fallback produces the same result as direct Timedelta conversion
    result = to_offset(freq)
    expected_td = Timedelta(freq)

    # Convert result to total seconds to compare magnitude
    assert result.nanos == expected_td.value


def test_to_offset_invalid_still_raises():
    # Ensure that if BOTH fail, we still get the Invalid Frequency error
    with pytest.raises(ValueError, match="Invalid frequency"):
        to_offset("not_a_time_at_all")
