# -*- coding: utf-8 -*-
from pandas import Timedelta


def test_repr():
    assert (repr(Timedelta(10, unit='d')) ==
            "Timedelta('10 days 00:00:00')")
    assert (repr(Timedelta(10, unit='s')) ==
            "Timedelta('0 days 00:00:10')")
    assert (repr(Timedelta(10, unit='ms')) ==
            "Timedelta('0 days 00:00:00.010000')")
    assert (repr(Timedelta(-10, unit='ms')) ==
            "Timedelta('-1 days +23:59:59.990000')")


def test_isoformat():
    td = Timedelta(days=6, minutes=50, seconds=3,
                   milliseconds=10, microseconds=10, nanoseconds=12)
    expected = 'P6DT0H50M3.010010012S'
    result = td.isoformat()
    assert result == expected

    td = Timedelta(days=4, hours=12, minutes=30, seconds=5)
    result = td.isoformat()
    expected = 'P4DT12H30M5S'
    assert result == expected

    td = Timedelta(nanoseconds=123)
    result = td.isoformat()
    expected = 'P0DT0H0M0.000000123S'
    assert result == expected

    # trim nano
    td = Timedelta(microseconds=10)
    result = td.isoformat()
    expected = 'P0DT0H0M0.00001S'
    assert result == expected

    # trim micro
    td = Timedelta(milliseconds=1)
    result = td.isoformat()
    expected = 'P0DT0H0M0.001S'
    assert result == expected

    # don't strip every 0
    result = Timedelta(minutes=1).isoformat()
    expected = 'P0DT0H1M0S'
    assert result == expected
