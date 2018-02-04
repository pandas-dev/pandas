# -*- coding: utf-8 -*-
"""Tests specific to tslibs.timezones"""
from datetime import datetime

import pytz
import dateutil.tz
import pytest
import numpy as np

import pandas.util.testing as tm

from pandas._libs import tslib
from pandas._libs.tslibs import timezones, conversion
from pandas import Timestamp, date_range


def test_maybe_get_tz_tzlocal(self):
    # GH#13583
    ts = Timestamp('2011-01-01', tz=dateutil.tz.tzlocal())
    assert ts.tz == dateutil.tz.tzlocal()
    assert "tz='tzlocal()')" in repr(ts)

    tz = timezones.maybe_get_tz('tzlocal()')
    assert tz == dateutil.tz.tzlocal()

    # get offset using normal datetime for test
    offset = dateutil.tz.tzlocal().utcoffset(datetime(2011, 1, 1))
    offset = offset.total_seconds() * 1000000000
    assert ts.value + offset == Timestamp('2011-01-01').value


class TestTimeZoneCacheKey(object):

    @pytest.mark.parametrize('tz_name', list(pytz.common_timezones))
    def test_cache_keys_are_distinct_for_pytz_vs_dateutil(self, tz_name):
        if tz_name == 'UTC':
            # skip utc as it's a special case in dateutil
            return
        tz_p = timezones.maybe_get_tz(tz_name)
        tz_d = timezones.maybe_get_tz('dateutil/' + tz_name)
        if tz_d is None:
            # skip timezones that dateutil doesn't know about.
            return
        assert (timezones._p_tz_cache_key(tz_p) !=
                timezones._p_tz_cache_key(tz_d))


class TestConversion(object):

    def test_tz_convert_empty(self):
        # Check empty array
        result = conversion.tz_convert(np.array([], dtype=np.int64),
                                       timezones.maybe_get_tz('US/Eastern'),
                                       timezones.maybe_get_tz('Asia/Tokyo'))
        expected = np.array([], dtype=np.int64)
        tm.assert_numpy_array_equal(result, expected)

    def test_tz_convert_nat(self):
        # Check all-NaT array
        result = conversion.tz_convert(np.array([tslib.iNaT], dtype=np.int64),
                                       timezones.maybe_get_tz('US/Eastern'),
                                       timezones.maybe_get_tz('Asia/Tokyo'))
        expected = np.array([tslib.iNaT], dtype=np.int64)
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize('tz', ['UTC', 'Asia/Tokyo',
                                    'US/Eastern', 'Europe/Moscow'])
    def test_tslib_tz_convert(self, tz):
        def compare_utc_to_local(tz_didx, utc_didx):
            f = lambda x: conversion.tz_convert_single(x, 'UTC', tz_didx.tz)
            result = conversion.tz_convert(tz_didx.asi8, 'UTC', tz_didx.tz)
            result_single = np.vectorize(f)(tz_didx.asi8)
            tm.assert_numpy_array_equal(result, result_single)

        def compare_local_to_utc(tz_didx, utc_didx):
            f = lambda x: conversion.tz_convert_single(x, tz_didx.tz, 'UTC')
            result = conversion.tz_convert(utc_didx.asi8, tz_didx.tz, 'UTC')
            result_single = np.vectorize(f)(utc_didx.asi8)
            tm.assert_numpy_array_equal(result, result_single)

        # US: 2014-03-09 - 2014-11-11
        # MOSCOW: 2014-10-26  /  2014-12-31
        tz_didx = date_range('2014-03-01', '2015-01-10', freq='H', tz=tz)
        utc_didx = date_range('2014-03-01', '2015-01-10', freq='H')
        compare_utc_to_local(tz_didx, utc_didx)
        # local tz to UTC can be differ in hourly (or higher) freqs because
        # of DST
        compare_local_to_utc(tz_didx, utc_didx)

        tz_didx = date_range('2000-01-01', '2020-01-01', freq='D', tz=tz)
        utc_didx = date_range('2000-01-01', '2020-01-01', freq='D')
        compare_utc_to_local(tz_didx, utc_didx)
        compare_local_to_utc(tz_didx, utc_didx)

        tz_didx = date_range('2000-01-01', '2100-01-01', freq='A', tz=tz)
        utc_didx = date_range('2000-01-01', '2100-01-01', freq='A')
        compare_utc_to_local(tz_didx, utc_didx)
        compare_local_to_utc(tz_didx, utc_didx)
