# -*- coding: utf-8 -*-
"""Tests for functions from pandas._libs.tslibs"""

from datetime import date, datetime

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


def test_rich_comparison_with_unsupported_type():
    # See https://github.com/pandas-dev/pandas/issues/24011

    class Inf:
        def __lt__(self, o):
            return False

        def __le__(self, o):
            return isinstance(o, Inf)

        def __gt__(self, o):
            return not isinstance(o, Inf)

        def __ge__(self, o):
            return True

        def __eq__(self, o):
            return isinstance(o, Inf)

    timestamp = tslibs.Timestamp('2018-11-30')

    # Comparison works if compared in *that* order, because
    # magic method is called on Inf
    assert Inf() > timestamp
    assert not (Inf() < timestamp)

    # ... but used to not work when magic method is called on Timestamp
    assert timestamp < Inf()
