# -*- coding: utf-8 -*-
"""
Tests for Series timezone-related methods
"""
from datetime import datetime

import pytest
import pytz
from dateutil.tz import tzoffset
import numpy as np

import pandas.util.testing as tm
from pandas.compat import lrange
from pandas.core.indexes.datetimes import date_range
from pandas import Index, Series, Timestamp, DatetimeIndex


class TestSeriesTimezones(object):
    # ----------------------------------------------------------------
    # Series.append

    def test_ser_append_aware(self):
        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H',
                          tz='US/Eastern')
        rng2 = date_range('1/1/2011 02:00', periods=1, freq='H',
                          tz='US/Eastern')
        ts1 = Series([1], index=rng1)
        ts2 = Series([2], index=rng2)
        ts_result = ts1.append(ts2)

        exp_index = DatetimeIndex(['2011-01-01 01:00', '2011-01-01 02:00'],
                                  tz='US/Eastern')
        exp = Series([1, 2], index=exp_index)
        tm.assert_series_equal(ts_result, exp)
        assert ts_result.index.tz == rng1.tz

        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H', tz='UTC')
        rng2 = date_range('1/1/2011 02:00', periods=1, freq='H', tz='UTC')
        ts1 = Series([1], index=rng1)
        ts2 = Series([2], index=rng2)
        ts_result = ts1.append(ts2)

        exp_index = DatetimeIndex(['2011-01-01 01:00', '2011-01-01 02:00'],
                                  tz='UTC')
        exp = Series([1, 2], index=exp_index)
        tm.assert_series_equal(ts_result, exp)
        utc = rng1.tz
        assert utc == ts_result.index.tz

        # GH#7795
        # different tz coerces to object dtype, not UTC
        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H',
                          tz='US/Eastern')
        rng2 = date_range('1/1/2011 02:00', periods=1, freq='H',
                          tz='US/Central')
        ts1 = Series([1], index=rng1)
        ts2 = Series([2], index=rng2)
        ts_result = ts1.append(ts2)
        exp_index = Index([Timestamp('1/1/2011 01:00', tz='US/Eastern'),
                           Timestamp('1/1/2011 02:00', tz='US/Central')])
        exp = Series([1, 2], index=exp_index)
        tm.assert_series_equal(ts_result, exp)

    def test_ser_append_aware_naive(self):
        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H')
        rng2 = date_range('1/1/2011 02:00', periods=1, freq='H',
                          tz='US/Eastern')
        ser1 = Series(np.random.randn(len(rng1)), index=rng1)
        ser2 = Series(np.random.randn(len(rng2)), index=rng2)
        ts_result = ser1.append(ser2)

        assert ts_result.index.equals(
            ser1.index.astype(object).append(ser2.index.astype(object)))

        # mixed
        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H')
        rng2 = lrange(100)
        ser1 = Series(np.random.randn(len(rng1)), index=rng1)
        ser2 = Series(np.random.randn(len(rng2)), index=rng2)
        ts_result = ser1.append(ser2)
        assert ts_result.index.equals(
            ser1.index.astype(object).append(ser2.index))

    def test_ser_append_dst(self):
        rng1 = date_range('1/1/2016 01:00', periods=3, freq='H',
                          tz='US/Eastern')
        rng2 = date_range('8/1/2016 01:00', periods=3, freq='H',
                          tz='US/Eastern')
        ser1 = Series([1, 2, 3], index=rng1)
        ser2 = Series([10, 11, 12], index=rng2)
        ts_result = ser1.append(ser2)

        exp_index = DatetimeIndex(['2016-01-01 01:00', '2016-01-01 02:00',
                                   '2016-01-01 03:00', '2016-08-01 01:00',
                                   '2016-08-01 02:00', '2016-08-01 03:00'],
                                  tz='US/Eastern')
        exp = Series([1, 2, 3, 10, 11, 12], index=exp_index)
        tm.assert_series_equal(ts_result, exp)
        assert ts_result.index.tz == rng1.tz

    # ----------------------------------------------------------------
    # Series.tz_localize

    def test_ser_tz_localize(self):
        rng = date_range('1/1/2011', periods=100, freq='H')
        ser = Series(1, index=rng)

        result = ser.tz_localize('utc')
        assert result.index.tz.zone == 'UTC'

        # Can't localize if already tz-aware
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')
        ser = Series(1, index=rng)
        tm.assert_raises_regex(TypeError, 'Already tz-aware',
                               ser.tz_localize, 'US/Eastern')

    def test_ser_tz_localize_ambiguous_bool(self):
        # make sure that we are correctly accepting bool values as ambiguous

        # GH#14402
        naive = Timestamp('2015-11-01 01:00:03')
        expected0 = Timestamp('2015-11-01 01:00:03-0500', tz='US/Central')
        expected1 = Timestamp('2015-11-01 01:00:03-0600', tz='US/Central')

        with pytest.raises(pytz.AmbiguousTimeError):
            naive.tz_localize('US/Central')

        result = naive.tz_localize('US/Central', ambiguous=True)
        assert result == expected0

        result = naive.tz_localize('US/Central', ambiguous=False)
        assert result == expected1

        ser = Series([naive])
        expected0 = Series([expected0])
        expected1 = Series([expected1])

        with pytest.raises(pytz.AmbiguousTimeError):
            ser.dt.tz_localize('US/Central')

        result = ser.dt.tz_localize('US/Central', ambiguous=True)
        tm.assert_series_equal(result, expected0)

        result = ser.dt.tz_localize('US/Central', ambiguous=[True])
        tm.assert_series_equal(result, expected0)

        result = ser.dt.tz_localize('US/Central', ambiguous=False)
        tm.assert_series_equal(result, expected1)

        result = ser.dt.tz_localize('US/Central', ambiguous=[False])
        tm.assert_series_equal(result, expected1)

    # ----------------------------------------------------------------

    def test_ser_dateutil_tzoffset_support(self):
        values = [188.5, 328.25]
        tzinfo = tzoffset(None, 7200)
        index = [datetime(2012, 5, 11, 11, tzinfo=tzinfo),
                 datetime(2012, 5, 11, 12, tzinfo=tzinfo)]
        series = Series(data=values, index=index)

        assert series.index.tz == tzinfo

        # it works! GH#2443
        repr(series.index[0])

    def test_ser_tz_convert(self):
        rng = date_range('1/1/2011', periods=200, freq='D', tz='US/Eastern')
        ser = Series(1, index=rng)

        result = ser.tz_convert('Europe/Berlin')
        assert result.index.tz.zone == 'Europe/Berlin'

        # can't convert tz-naive
        rng = date_range('1/1/2011', periods=200, freq='D')
        ser = Series(1, index=rng)
        tm.assert_raises_regex(TypeError, "Cannot convert tz-naive",
                               ser.tz_convert, 'US/Eastern')

    # ----------------------------------------------------------------
    # Series.__add__

    def test_ser_add_tzcompat_raises(self):
        rng = date_range('1/1/2011', periods=10, freq='H')
        ser = Series(np.random.randn(len(rng)), index=rng)
        ser_utc = ser.tz_localize('utc')

        with pytest.raises(Exception):
            ser + ser_utc
        with pytest.raises(Exception):
            ser_utc + ser

    def test_ser_add_tz_mismatch_converts_to_utc(self):
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')

        perm = np.random.permutation(100)[:90]
        ser1 = Series(np.random.randn(90),
                      index=rng.take(perm).tz_convert('US/Eastern'))

        perm = np.random.permutation(100)[:90]
        ser2 = Series(np.random.randn(90),
                      index=rng.take(perm).tz_convert('Europe/Berlin'))

        result = ser1 + ser2

        uts1 = ser1.tz_convert('utc')
        uts2 = ser2.tz_convert('utc')
        expected = uts1 + uts2

        assert result.index.tz == pytz.UTC
        tm.assert_series_equal(result, expected)

    # TODO: test is redundant with test_ser_add_tz_mismatch_converts_to_utc
    def test_ser_add_tz_mismatch_returns_utc(self):
        rng = date_range('1/1/2011', periods=10, freq='H', tz='US/Eastern')
        ser = Series(np.random.randn(len(rng)), index=rng)

        ts_moscow = ser.tz_convert('Europe/Moscow')

        result = ser + ts_moscow
        assert result.index.tz is pytz.utc

        result = ts_moscow + ser
        assert result.index.tz is pytz.utc
