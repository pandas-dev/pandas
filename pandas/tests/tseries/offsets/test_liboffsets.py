# -*- coding: utf-8 -*-
"""
Tests for helper functions in the cython tslibs.offsets
"""
from datetime import datetime

import pytest

import pandas as pd

import pandas._libs.tslibs.offsets as liboffsets


def test_shift_month():
    dt = datetime(2017, 11, 15)

    assert liboffsets.shift_month(dt, 0, day_opt=None) == dt
    assert liboffsets.shift_month(dt, 0, day_opt=15) == dt

    assert liboffsets.shift_month(dt, 1,
                                  day_opt='start') == datetime(2017, 12, 1)

    assert liboffsets.shift_month(dt, -145,
                                  day_opt='end') == datetime(2005, 10, 31)

    with pytest.raises(ValueError):
        liboffsets.shift_month(dt, 3, day_opt='this should raise')


def test_get_day_of_month():
    # get_day_of_month is not directly exposed; we test it via roll_yearday
    dt = datetime(2017, 11, 15)

    with pytest.raises(ValueError):
        # To hit the raising case we need month == dt.month and n > 0
        liboffsets.roll_yearday(dt, n=3, month=11, day_opt='foo')


def test_roll_yearday():
    # Copied from doctest examples
    month = 3
    day_opt = 'start'              # `other` will be compared to March 1
    other = datetime(2017, 2, 10)  # before March 1
    assert liboffsets.roll_yearday(other, 2, month, day_opt) == 1
    assert liboffsets.roll_yearday(other, -7, month, day_opt) == -7
    assert liboffsets.roll_yearday(other, 0, month, day_opt) == 0

    other = pd.Timestamp('2014-03-15', tz='US/Eastern')  # after March 1
    assert liboffsets.roll_yearday(other, 2, month, day_opt) == 2
    assert liboffsets.roll_yearday(other, -7, month, day_opt) == -6
    assert liboffsets.roll_yearday(other, 0, month, day_opt) == 1

    month = 6
    day_opt = 'end'                # `other` will be compared to June 30
    other = datetime(1999, 6, 29)  # before June 30
    assert liboffsets.roll_yearday(other, 5, month, day_opt) == 4
    assert liboffsets.roll_yearday(other, -7, month, day_opt) == -7
    assert liboffsets.roll_yearday(other, 0, month, day_opt) == 0

    other = pd.Timestamp(2072, 8, 24, 6, 17, 18)  # after June 30
    assert liboffsets.roll_yearday(other, 5, month, day_opt) == 5
    assert liboffsets.roll_yearday(other, -7, month, day_opt) == -6
    assert liboffsets.roll_yearday(other, 0, month, day_opt) == 1
