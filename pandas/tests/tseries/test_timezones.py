# pylint: disable-msg=E1101,W0612
import pytest

import pytz
import dateutil
import numpy as np

from dateutil.tz import tzoffset
from datetime import datetime

import pandas.util.testing as tm
from pandas.compat import lrange
from pandas.core.indexes.datetimes import date_range
from pandas._libs import tslib
from pandas._libs.tslibs import timezones, conversion
from pandas import Index, Series, Timestamp, DatetimeIndex
from pandas.util.testing import assert_series_equal


class TestTimeZoneSupportPytz(object):

    def tz(self, tz):
        # Construct a timezone object from a string. Overridden in subclass to
        # parameterize tests.
        return pytz.timezone(tz)

    def tzstr(self, tz):
        # Construct a timezone string from a string. Overridden in subclass to
        # parameterize tests.
        return tz

    def localize(self, tz, x):
        return tz.localize(x)

    def normalize(self, ts):
        tzinfo = ts.tzinfo
        return tzinfo.normalize(ts)

    def cmptz(self, tz1, tz2):
        # Compare two timezones. Overridden in subclass to parameterize
        # tests.
        return tz1.zone == tz2.zone

    def test_tz_localize_empty_series(self):
        # #2248

        ts = Series()

        ts2 = ts.tz_localize('utc')
        assert ts2.index.tz == pytz.utc

        ts2 = ts.tz_localize(self.tzstr('US/Eastern'))
        assert self.cmptz(ts2.index.tz, self.tz('US/Eastern'))

    def test_ambiguous_bool(self):
        # make sure that we are correctly accepting bool values as ambiguous

        # gh-14402
        t = Timestamp('2015-11-01 01:00:03')
        expected0 = Timestamp('2015-11-01 01:00:03-0500', tz='US/Central')
        expected1 = Timestamp('2015-11-01 01:00:03-0600', tz='US/Central')

        s = Series([t])
        expected0 = Series([expected0])
        expected1 = Series([expected1])

        def f():
            s.dt.tz_localize('US/Central')
        pytest.raises(pytz.AmbiguousTimeError, f)

        result = s.dt.tz_localize('US/Central', ambiguous=True)
        assert_series_equal(result, expected0)

        result = s.dt.tz_localize('US/Central', ambiguous=[True])
        assert_series_equal(result, expected0)

        result = s.dt.tz_localize('US/Central', ambiguous=False)
        assert_series_equal(result, expected1)

        result = s.dt.tz_localize('US/Central', ambiguous=[False])
        assert_series_equal(result, expected1)

    # test utility methods
    def test_infer_tz(self):
        eastern = self.tz('US/Eastern')
        utc = pytz.utc

        _start = datetime(2001, 1, 1)
        _end = datetime(2009, 1, 1)

        start = self.localize(eastern, _start)
        end = self.localize(eastern, _end)
        assert (timezones.infer_tzinfo(start, end) is
                self.localize(eastern, _start).tzinfo)
        assert (timezones.infer_tzinfo(start, None) is
                self.localize(eastern, _start).tzinfo)
        assert (timezones.infer_tzinfo(None, end) is
                self.localize(eastern, _end).tzinfo)

        start = utc.localize(_start)
        end = utc.localize(_end)
        assert (timezones.infer_tzinfo(start, end) is utc)

        end = self.localize(eastern, _end)
        pytest.raises(Exception, timezones.infer_tzinfo, start, end)
        pytest.raises(Exception, timezones.infer_tzinfo, end, start)

    def test_localized_at_time_between_time(self):
        from datetime import time

        rng = date_range('4/16/2012', '5/1/2012', freq='H')
        ts = Series(np.random.randn(len(rng)), index=rng)

        ts_local = ts.tz_localize(self.tzstr('US/Eastern'))

        result = ts_local.at_time(time(10, 0))
        expected = ts.at_time(time(10, 0)).tz_localize(self.tzstr(
            'US/Eastern'))
        assert_series_equal(result, expected)
        assert self.cmptz(result.index.tz, self.tz('US/Eastern'))

        t1, t2 = time(10, 0), time(11, 0)
        result = ts_local.between_time(t1, t2)
        expected = ts.between_time(t1,
                                   t2).tz_localize(self.tzstr('US/Eastern'))
        assert_series_equal(result, expected)
        assert self.cmptz(result.index.tz, self.tz('US/Eastern'))

    def test_string_index_alias_tz_aware(self):
        rng = date_range('1/1/2000', periods=10, tz=self.tzstr('US/Eastern'))
        ts = Series(np.random.randn(len(rng)), index=rng)

        result = ts['1/3/2000']
        tm.assert_almost_equal(result, ts[2])

    def test_tz_aware_asfreq(self):
        dr = date_range('2011-12-01', '2012-07-20', freq='D',
                        tz=self.tzstr('US/Eastern'))

        s = Series(np.random.randn(len(dr)), index=dr)

        # it works!
        s.asfreq('T')

    def test_dateutil_tzoffset_support(self):
        values = [188.5, 328.25]
        tzinfo = tzoffset(None, 7200)
        index = [datetime(2012, 5, 11, 11, tzinfo=tzinfo),
                 datetime(2012, 5, 11, 12, tzinfo=tzinfo)]
        series = Series(data=values, index=index)

        assert series.index.tz == tzinfo

        # it works! #2443
        repr(series.index[0])

    def test_getitem_pydatetime_tz(self):
        index = date_range(start='2012-12-24 16:00', end='2012-12-24 18:00',
                           freq='H', tz=self.tzstr('Europe/Berlin'))
        ts = Series(index=index, data=index.hour)
        time_pandas = Timestamp('2012-12-24 17:00',
                                tz=self.tzstr('Europe/Berlin'))
        time_datetime = self.localize(
            self.tz('Europe/Berlin'), datetime(2012, 12, 24, 17, 0))
        assert ts[time_pandas] == ts[time_datetime]

    def test_replace_across_dst(self):
        # GH#18319 check that 1) timezone is correctly normalized and
        # 2) that hour is not incorrectly changed by this normalization
        tz = self.tz('US/Eastern')

        ts_naive = Timestamp('2017-12-03 16:03:30')
        ts_aware = self.localize(tz, ts_naive)

        # Preliminary sanity-check
        assert ts_aware == self.normalize(ts_aware)

        # Replace across DST boundary
        ts2 = ts_aware.replace(month=6)

        # Check that `replace` preserves hour literal
        assert (ts2.hour, ts2.minute) == (ts_aware.hour, ts_aware.minute)

        # Check that post-replace object is appropriately normalized
        ts2b = self.normalize(ts2)
        assert ts2 == ts2b


class TestTimeZoneSupportDateutil(TestTimeZoneSupportPytz):

    def tz(self, tz):
        """
        Construct a dateutil timezone.
        Use tslib.maybe_get_tz so that we get the filename on the tz right
        on windows. See #7337.
        """
        return timezones.maybe_get_tz('dateutil/' + tz)

    def tzstr(self, tz):
        """ Construct a timezone string from a string. Overridden in subclass
        to parameterize tests. """
        return 'dateutil/' + tz

    def cmptz(self, tz1, tz2):
        """ Compare two timezones. Overridden in subclass to parameterize
        tests. """
        return tz1 == tz2

    def localize(self, tz, x):
        return x.replace(tzinfo=tz)

    def normalize(self, ts):
        # no-op for dateutil
        return ts

    def test_tzlocal(self):
        # GH 13583
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


class TestTimeZones(object):
    timezones = ['UTC', 'Asia/Tokyo', 'US/Eastern', 'dateutil/US/Pacific']

    def test_series_tz_localize(self):

        rng = date_range('1/1/2011', periods=100, freq='H')
        ts = Series(1, index=rng)

        result = ts.tz_localize('utc')
        assert result.index.tz.zone == 'UTC'

        # Can't localize if already tz-aware
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')
        ts = Series(1, index=rng)
        tm.assert_raises_regex(TypeError, 'Already tz-aware',
                               ts.tz_localize, 'US/Eastern')

    def test_series_tz_convert(self):
        rng = date_range('1/1/2011', periods=200, freq='D', tz='US/Eastern')
        ts = Series(1, index=rng)

        result = ts.tz_convert('Europe/Berlin')
        assert result.index.tz.zone == 'Europe/Berlin'

        # can't convert tz-naive
        rng = date_range('1/1/2011', periods=200, freq='D')
        ts = Series(1, index=rng)
        tm.assert_raises_regex(TypeError, "Cannot convert tz-naive",
                               ts.tz_convert, 'US/Eastern')

    def test_join_aware(self):
        rng = date_range('1/1/2011', periods=10, freq='H')
        ts = Series(np.random.randn(len(rng)), index=rng)

        ts_utc = ts.tz_localize('utc')

        pytest.raises(Exception, ts.__add__, ts_utc)
        pytest.raises(Exception, ts_utc.__add__, ts)

    def test_series_align_aware(self):
        idx1 = date_range('2001', periods=5, freq='H', tz='US/Eastern')
        ser = Series(np.random.randn(len(idx1)), index=idx1)
        ser_central = ser.tz_convert('US/Central')
        # # different timezones convert to UTC

        new1, new2 = ser.align(ser_central)
        assert new1.index.tz == pytz.UTC
        assert new2.index.tz == pytz.UTC

    def test_append_aware(self):
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
        assert_series_equal(ts_result, exp)
        assert ts_result.index.tz == rng1.tz

        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H', tz='UTC')
        rng2 = date_range('1/1/2011 02:00', periods=1, freq='H', tz='UTC')
        ts1 = Series([1], index=rng1)
        ts2 = Series([2], index=rng2)
        ts_result = ts1.append(ts2)

        exp_index = DatetimeIndex(['2011-01-01 01:00', '2011-01-01 02:00'],
                                  tz='UTC')
        exp = Series([1, 2], index=exp_index)
        assert_series_equal(ts_result, exp)
        utc = rng1.tz
        assert utc == ts_result.index.tz

        # GH 7795
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
        assert_series_equal(ts_result, exp)

    def test_append_dst(self):
        rng1 = date_range('1/1/2016 01:00', periods=3, freq='H',
                          tz='US/Eastern')
        rng2 = date_range('8/1/2016 01:00', periods=3, freq='H',
                          tz='US/Eastern')
        ts1 = Series([1, 2, 3], index=rng1)
        ts2 = Series([10, 11, 12], index=rng2)
        ts_result = ts1.append(ts2)

        exp_index = DatetimeIndex(['2016-01-01 01:00', '2016-01-01 02:00',
                                   '2016-01-01 03:00', '2016-08-01 01:00',
                                   '2016-08-01 02:00', '2016-08-01 03:00'],
                                  tz='US/Eastern')
        exp = Series([1, 2, 3, 10, 11, 12], index=exp_index)
        assert_series_equal(ts_result, exp)
        assert ts_result.index.tz == rng1.tz

    def test_append_aware_naive(self):
        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H')
        rng2 = date_range('1/1/2011 02:00', periods=1, freq='H',
                          tz='US/Eastern')
        ts1 = Series(np.random.randn(len(rng1)), index=rng1)
        ts2 = Series(np.random.randn(len(rng2)), index=rng2)
        ts_result = ts1.append(ts2)

        assert ts_result.index.equals(ts1.index.astype(object).append(
            ts2.index.astype(object)))

        # mixed
        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H')
        rng2 = lrange(100)
        ts1 = Series(np.random.randn(len(rng1)), index=rng1)
        ts2 = Series(np.random.randn(len(rng2)), index=rng2)
        ts_result = ts1.append(ts2)
        assert ts_result.index.equals(ts1.index.astype(object).append(
            ts2.index))

    def test_series_add_tz_mismatch_converts_to_utc(self):
        rng = date_range('1/1/2011', periods=10, freq='H', tz='US/Eastern')
        ts = Series(np.random.randn(len(rng)), index=rng)

        ts_moscow = ts.tz_convert('Europe/Moscow')

        result = ts + ts_moscow
        assert result.index.tz is pytz.utc

        result = ts_moscow + ts
        assert result.index.tz is pytz.utc

    def test_arith_utc_convert(self):
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')

        perm = np.random.permutation(100)[:90]
        ts1 = Series(np.random.randn(90),
                     index=rng.take(perm).tz_convert('US/Eastern'))

        perm = np.random.permutation(100)[:90]
        ts2 = Series(np.random.randn(90),
                     index=rng.take(perm).tz_convert('Europe/Berlin'))

        result = ts1 + ts2

        uts1 = ts1.tz_convert('utc')
        uts2 = ts2.tz_convert('utc')
        expected = uts1 + uts2

        assert result.index.tz == pytz.UTC
        assert_series_equal(result, expected)


class TestTslib(object):

    def test_tslib_tz_convert(self):
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

        for tz in ['UTC', 'Asia/Tokyo', 'US/Eastern', 'Europe/Moscow']:
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

        # Check empty array
        result = conversion.tz_convert(np.array([], dtype=np.int64),
                                       timezones.maybe_get_tz('US/Eastern'),
                                       timezones.maybe_get_tz('Asia/Tokyo'))
        tm.assert_numpy_array_equal(result, np.array([], dtype=np.int64))

        # Check all-NaT array
        result = conversion.tz_convert(np.array([tslib.iNaT], dtype=np.int64),
                                       timezones.maybe_get_tz('US/Eastern'),
                                       timezones.maybe_get_tz('Asia/Tokyo'))
        tm.assert_numpy_array_equal(result, np.array(
            [tslib.iNaT], dtype=np.int64))
