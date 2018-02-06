# pylint: disable-msg=E1101,W0612
import pytest

import pytz
import dateutil
import numpy as np

from datetime import datetime

import pandas.util.testing as tm
from pandas.core.indexes.datetimes import date_range
from pandas._libs import tslib
from pandas._libs.tslibs import timezones, conversion
from pandas import Timestamp


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
