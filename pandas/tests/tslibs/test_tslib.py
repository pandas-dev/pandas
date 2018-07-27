# -*- coding: utf-8 -*-
"""Tests for functions from pandas._libs.tslibs"""

from datetime import datetime, date

from pandas._libs import tslibs


def test_normalize_date():
    value = date(2012, 9, 7)

    result = tslibs.normalize_date(value)
    assert (result == datetime(2012, 9, 7))

    value = datetime(2012, 9, 7, 12)

    result = tslibs.normalize_date(value)
    assert (result == datetime(2012, 9, 7))

    value = datetime(2007, 10, 1, 1, 12, 5, 10)

    actual = tslibs.normalize_date(value)
    assert actual == datetime(2007, 10, 1)
