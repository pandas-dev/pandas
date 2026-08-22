import pytest

from pandas import TimedeltaIndex


@pytest.mark.parametrize(
    "timedelta_val,expected",
    [
        ("1 day", "day"),
        ("1 hour", "hour"),
        ("1 minute", "minute"),
        ("1 second", "second"),
        ("1 millisecond", "millisecond"),
        ("1 microsecond", "microsecond"),
    ],
)
def test_tdi_resolution(timedelta_val, expected):
    # GH#65186
    tdi = TimedeltaIndex([timedelta_val])
    assert tdi.resolution == expected


def test_tdi_resolution_no_freq():
    # GH#65186 resolution should work even without a freq
    tdi = TimedeltaIndex(["1 day", "2 days", "4 days"])
    assert tdi.resolution == "day"

    tdi = TimedeltaIndex(["1 day 2 hours", "3 days"])
    assert tdi.resolution == "hour"
