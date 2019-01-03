# -*- coding: utf-8 -*-
"""Tests for functions from pandas._libs.tslibs"""

from datetime import date, datetime

import pytest

from pandas._libs import tslibs


@pytest.mark.parametrize("value,expected", [
    (date(2012, 9, 7), datetime(2012, 9, 7)),
    (datetime(2012, 9, 7, 12), datetime(2012, 9, 7)),
    (datetime(2007, 10, 1, 1, 12, 5, 10), datetime(2007, 10, 1))
])
def test_normalize_date(value, expected):
    result = tslibs.normalize_date(value)
    assert result == expected
