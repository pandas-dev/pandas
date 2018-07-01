# -*- coding: utf-8 -*-
"""Tests for functions from pandas._libs.tslibs"""

from datetime import datetime, datetime

from pandas._libs import tslib


def test_normalize_date():
    value = date(2012, 9, 7)

    result = tslib.normalize_date(value)
    assert (result == datetime(2012, 9, 7))

    value = datetime(2012, 9, 7, 12)

    result = tslib.normalize_date(value)
    assert (result == datetime(2012, 9, 7))
