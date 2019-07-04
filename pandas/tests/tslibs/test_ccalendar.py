from datetime import datetime

import numpy as np
import pytest

from pandas._libs.tslibs import ccalendar


@pytest.mark.parametrize(
    "date_tuple,expected",
    [
        ((2001, 3, 1), 60),
        ((2004, 3, 1), 61),
        ((1907, 12, 31), 365),  # End-of-year, non-leap year.
        ((2004, 12, 31), 366),  # End-of-year, leap year.
    ],
)
def test_get_day_of_year_numeric(date_tuple, expected):
    assert ccalendar.get_day_of_year(*date_tuple) == expected


def test_get_day_of_year_dt():
    dt = datetime.fromordinal(1 + np.random.randint(365 * 4000))
    result = ccalendar.get_day_of_year(dt.year, dt.month, dt.day)

    expected = (dt - dt.replace(month=1, day=1)).days + 1
    assert result == expected
