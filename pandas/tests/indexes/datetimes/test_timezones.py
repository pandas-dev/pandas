# -*- coding: utf-8 -*-
"""Tests for DatetimeIndex timezone operations"""
from datetime import datetime, timedelta

import pytest
import numpy as np

import dateutil
from dateutil.tz import tzlocal, gettz
import pytz
from pytz import NonExistentTimeError
from distutils.version import LooseVersion

from pandas.compat import zip
import pandas.util.testing as tm
import pandas.util._test_decorators as td

import pandas as pd
from pandas import (Index, NaT, DatetimeIndex,
                    date_range, Timestamp)
from pandas.util.testing import set_timezone


class TestDatetimeIndexTimezones(object):

    # ----------------------------------------------------------------
    # DatetimeIndex.tz_convert

    def test_dti_tz_convert_hour_overflow_dst(self):
        # Regression test for GH#13306

        # sorted case US/Eastern -> UTC
        ts = ['2008-05-12 09:50:00',
              '2008-12-12 09:50:35',
              '2009-05-12 09:50:32']
        tt = DatetimeIndex(ts).tz_localize('US/Eastern')
        ut = tt.tz_convert('UTC')
        expected = Index([13, 14, 13])
        tm.assert_index_equal(ut.hour, expected)

        # sorted case UTC -> US/Eastern
        ts = ['2008-05-12 13:50:00',
              '2008-12-12 14:50:35',
              '2009-05-12 13:50:32']
        tt = DatetimeIndex(ts).tz_localize('UTC')
        ut = tt.tz_convert('US/Eastern')
        expected = Index([9, 9, 9])
        tm.assert_index_equal(ut.hour, expected)

        # unsorted case US/Eastern -> UTC
        ts = ['2008-05-12 09:50:00',
              '2008-12-12 09:50:35',
              '2008-05-12 09:50:32']
        tt = DatetimeIndex(ts).tz_localize('US/Eastern')
        ut = tt.tz_convert('UTC')
        expected = Index([13, 14, 13])
        tm.assert_index_equal(ut.hour, expected)

        # unsorted case UTC -> US/Eastern
        ts = ['2008-05-12 13:50:00',
              '2008-12-12 14:50:35',
              '2008-05-12 13:50:32']
        tt = DatetimeIndex(ts).tz_localize('UTC')
        ut = tt.tz_convert('US/Eastern')
        expected = Index([9, 9, 9])
        tm.assert_index_equal(ut.hour, expected)

    @pytest.mark.parametrize('tz', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_tz_convert_hour_overflow_dst_timestamps(self, tz):
        # GH#13306

        # sorted case US/Eastern -> UTC
        ts = [Timestamp('2008-05-12 09:50:00', tz=tz),
              Timestamp('2008-12-12 09:50:35', tz=tz),
              Timestamp('2009-05-12 09:50:32', tz=tz)]
        tt = DatetimeIndex(ts)
        ut = tt.tz_convert('UTC')
        expected = Index([13, 14, 13])
        tm.assert_index_equal(ut.hour, expected)

        # sorted case UTC -> US/Eastern
        ts = [Timestamp('2008-05-12 13:50:00', tz='UTC'),
              Timestamp('2008-12-12 14:50:35', tz='UTC'),
              Timestamp('2009-05-12 13:50:32', tz='UTC')]
        tt = DatetimeIndex(ts)
        ut = tt.tz_convert('US/Eastern')
        expected = Index([9, 9, 9])
        tm.assert_index_equal(ut.hour, expected)

        # unsorted case US/Eastern -> UTC
        ts = [Timestamp('2008-05-12 09:50:00', tz=tz),
              Timestamp('2008-12-12 09:50:35', tz=tz),
              Timestamp('2008-05-12 09:50:32', tz=tz)]
        tt = DatetimeIndex(ts)
        ut = tt.tz_convert('UTC')
        expected = Index([13, 14, 13])
        tm.assert_index_equal(ut.hour, expected)

        # unsorted case UTC -> US/Eastern
        ts = [Timestamp('2008-05-12 13:50:00', tz='UTC'),
              Timestamp('2008-12-12 14:50:35', tz='UTC'),
              Timestamp('2008-05-12 13:50:32', tz='UTC')]
        tt = DatetimeIndex(ts)
        ut = tt.tz_convert('US/Eastern')
        expected = Index([9, 9, 9])
        tm.assert_index_equal(ut.hour, expected)

    def test_dti_tz_convert_trans_pos_plus_1__bug(self):
        # Regression test for tslib.tz_convert(vals, tz1, tz2).
        # See https://github.com/pandas-dev/pandas/issues/4496 for details.
        for freq, n in [('H', 1), ('T', 60), ('S', 3600)]:
            idx = date_range(datetime(2011, 3, 26, 23),
                             datetime(2011, 3, 27, 1), freq=freq)
            idx = idx.tz_localize('UTC')
            idx = idx.tz_convert('Europe/Moscow')

            expected = np.repeat(np.array([3, 4, 5]), np.array([n, n, 1]))
            tm.assert_index_equal(idx.hour, Index(expected))

    def test_dti_tz_convert_dst(self):
        for freq, n in [('H', 1), ('T', 60), ('S', 3600)]:
            # Start DST
            idx = date_range('2014-03-08 23:00', '2014-03-09 09:00', freq=freq,
                             tz='UTC')
            idx = idx.tz_convert('US/Eastern')
            expected = np.repeat(np.array([18, 19, 20, 21, 22, 23,
                                           0, 1, 3, 4, 5]),
                                 np.array([n, n, n, n, n, n, n, n, n, n, 1]))
            tm.assert_index_equal(idx.hour, Index(expected))

            idx = date_range('2014-03-08 18:00', '2014-03-09 05:00', freq=freq,
                             tz='US/Eastern')
            idx = idx.tz_convert('UTC')
            expected = np.repeat(np.array([23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
                                 np.array([n, n, n, n, n, n, n, n, n, n, 1]))
            tm.assert_index_equal(idx.hour, Index(expected))

            # End DST
            idx = date_range('2014-11-01 23:00', '2014-11-02 09:00', freq=freq,
                             tz='UTC')
            idx = idx.tz_convert('US/Eastern')
            expected = np.repeat(np.array([19, 20, 21, 22, 23,
                                           0, 1, 1, 2, 3, 4]),
                                 np.array([n, n, n, n, n, n, n, n, n, n, 1]))
            tm.assert_index_equal(idx.hour, Index(expected))

            idx = date_range('2014-11-01 18:00', '2014-11-02 05:00', freq=freq,
                             tz='US/Eastern')
            idx = idx.tz_convert('UTC')
            expected = np.repeat(np.array([22, 23, 0, 1, 2, 3, 4, 5, 6,
                                           7, 8, 9, 10]),
                                 np.array([n, n, n, n, n, n, n, n, n,
                                           n, n, n, 1]))
            tm.assert_index_equal(idx.hour, Index(expected))

        # daily
        # Start DST
        idx = date_range('2014-03-08 00:00', '2014-03-09 00:00', freq='D',
                         tz='UTC')
        idx = idx.tz_convert('US/Eastern')
        tm.assert_index_equal(idx.hour, Index([19, 19]))

        idx = date_range('2014-03-08 00:00', '2014-03-09 00:00', freq='D',
                         tz='US/Eastern')
        idx = idx.tz_convert('UTC')
        tm.assert_index_equal(idx.hour, Index([5, 5]))

        # End DST
        idx = date_range('2014-11-01 00:00', '2014-11-02 00:00', freq='D',
                         tz='UTC')
        idx = idx.tz_convert('US/Eastern')
        tm.assert_index_equal(idx.hour, Index([20, 20]))

        idx = date_range('2014-11-01 00:00', '2014-11-02 000:00', freq='D',
                         tz='US/Eastern')
        idx = idx.tz_convert('UTC')
        tm.assert_index_equal(idx.hour, Index([4, 4]))

    def test_dti_tz_convert_tzlocal(self):
        # GH#13583
        # tz_convert doesn't affect to internal
        dti = date_range(start='2001-01-01', end='2001-03-01', tz='UTC')
        dti2 = dti.tz_convert(dateutil.tz.tzlocal())
        tm.assert_numpy_array_equal(dti2.asi8, dti.asi8)

        dti = date_range(start='2001-01-01', end='2001-03-01',
                         tz=dateutil.tz.tzlocal())
        dti2 = dti.tz_convert(None)
        tm.assert_numpy_array_equal(dti2.asi8, dti.asi8)

    @pytest.mark.parametrize('tz', ['UTC', 'Asia/Tokyo',
                                    'US/Eastern', 'dateutil/US/Pacific'])
    def test_dti_tz_convert_roundtrip(self, tz):
        idx1 = date_range(start='2014-01-01', end='2014-12-31', freq='M',
                          tz='UTC')
        exp1 = date_range(start='2014-01-01', end='2014-12-31', freq='M')

        idx2 = date_range(start='2014-01-01', end='2014-12-31', freq='D',
                          tz='UTC')
        exp2 = date_range(start='2014-01-01', end='2014-12-31', freq='D')

        idx3 = date_range(start='2014-01-01', end='2014-03-01', freq='H',
                          tz='UTC')
        exp3 = date_range(start='2014-01-01', end='2014-03-01', freq='H')

        idx4 = date_range(start='2014-08-01', end='2014-10-31', freq='T',
                          tz='UTC')
        exp4 = date_range(start='2014-08-01', end='2014-10-31', freq='T')

        for idx, expected in [(idx1, exp1), (idx2, exp2), (idx3, exp3),
                              (idx4, exp4)]:
            converted = idx.tz_convert(tz)
            reset = converted.tz_convert(None)
            tm.assert_index_equal(reset, expected)
            assert reset.tzinfo is None

            expected = converted.tz_convert('UTC').tz_localize(None)
            tm.assert_index_equal(reset, expected)

    @pytest.mark.parametrize('tz', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_tz_convert_unsorted(self, tz):
        dr = date_range('2012-03-09', freq='H', periods=100, tz='utc')
        dr = dr.tz_convert(tz)

        result = dr[::-1].hour
        exp = dr.hour[::-1]
        tm.assert_almost_equal(result, exp)

    # ----------------------------------------------------------------
    # DatetimeIndex.tz_localize

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern'),
                                    'US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_tz_localize_ambiguous(self, tz):
        dti = DatetimeIndex(start='11/6/2011 1:59', end='11/6/2011 2:00',
                            freq='L')
        with pytest.raises(pytz.AmbiguousTimeError):
            dti.tz_localize(tz)

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern'),
                                    'US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_tz_localize_nonexistent(self, tz):
        dti = DatetimeIndex(start='3/13/2011 1:59', end='3/13/2011 2:00',
                            freq='L')
        with pytest.raises(pytz.NonExistentTimeError):
            dti.tz_localize(tz)

    # TODO: De-duplicate with two tests above
    def test_dti_tz_localize_nonexistent_raise_coerce(self):
        # GH#13057
        times = ['2015-03-08 01:00', '2015-03-08 02:00', '2015-03-08 03:00']
        index = DatetimeIndex(times)
        tz = 'US/Eastern'

        with pytest.raises(NonExistentTimeError):
            index.tz_localize(tz=tz)

        with pytest.raises(NonExistentTimeError):
            index.tz_localize(tz=tz, errors='raise')

        result = index.tz_localize(tz=tz, errors='coerce')
        test_times = ['2015-03-08 01:00-05:00', 'NaT',
                      '2015-03-08 03:00-04:00']
        expected = DatetimeIndex(test_times)\
            .tz_localize('UTC').tz_convert('US/Eastern')
        tm.assert_index_equal(result, expected)

    def test_dti_tz_localize_tzlocal(self):
        # GH#13583
        offset = dateutil.tz.tzlocal().utcoffset(datetime(2011, 1, 1))
        offset = int(offset.total_seconds() * 1000000000)

        dti = date_range(start='2001-01-01', end='2001-03-01')
        dti2 = dti.tz_localize(dateutil.tz.tzlocal())
        tm.assert_numpy_array_equal(dti2.asi8 + offset, dti.asi8)

        dti = date_range(start='2001-01-01', end='2001-03-01',
                         tz=dateutil.tz.tzlocal())
        dti2 = dti.tz_localize(None)
        tm.assert_numpy_array_equal(dti2.asi8 - offset, dti.asi8)

    def test_dti_tz_localize_naive(self):
        rng = date_range('1/1/2011', periods=100, freq='H')

        conv = rng.tz_localize('US/Pacific')
        exp = date_range('1/1/2011', periods=100, freq='H', tz='US/Pacific')

        tm.assert_index_equal(conv, exp)

    @pytest.mark.parametrize('tz', ['UTC', 'Asia/Tokyo',
                                    'US/Eastern', 'dateutil/US/Pacific'])
    def test_dti_tz_localize_roundtrip(self, tz):
        idx1 = date_range(start='2014-01-01', end='2014-12-31', freq='M')
        idx2 = date_range(start='2014-01-01', end='2014-12-31', freq='D')
        idx3 = date_range(start='2014-01-01', end='2014-03-01', freq='H')
        idx4 = date_range(start='2014-08-01', end='2014-10-31', freq='T')
        for idx in [idx1, idx2, idx3, idx4]:
            localized = idx.tz_localize(tz)
            expected = date_range(start=idx[0], end=idx[-1], freq=idx.freq,
                                  tz=tz)
            tm.assert_index_equal(localized, expected)

            with pytest.raises(TypeError):
                localized.tz_localize(tz)

            reset = localized.tz_localize(None)
            tm.assert_index_equal(reset, idx)
            assert reset.tzinfo is None

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern'),
                                    'US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_tz_localize_utc_conversion(self, tz):
        # Localizing to time zone should:
        #  1) check for DST ambiguities
        #  2) convert to UTC

        rng = date_range('3/10/2012', '3/11/2012', freq='30T')

        converted = rng.tz_localize(tz)
        expected_naive = rng + pd.offsets.Hour(5)
        tm.assert_numpy_array_equal(converted.asi8, expected_naive.asi8)

        # DST ambiguity, this should fail
        rng = date_range('3/11/2012', '3/12/2012', freq='30T')
        # Is this really how it should fail??
        with pytest.raises(NonExistentTimeError):
            rng.tz_localize(tz)

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern'),
                                    'US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_tz_localize_ambiguous_nat(self, tz):
        times = ['11/06/2011 00:00', '11/06/2011 01:00', '11/06/2011 01:00',
                 '11/06/2011 02:00', '11/06/2011 03:00']
        naive = DatetimeIndex(times)
        localized = naive.tz_localize(tz, ambiguous='NaT')

        times = ['11/06/2011 00:00', np.NaN, np.NaN, '11/06/2011 02:00',
                 '11/06/2011 03:00']
        di_test = DatetimeIndex(times, tz='US/Eastern')

        # left dtype is datetime64[ns, US/Eastern]
        # right is datetime64[ns, tzfile('/usr/share/zoneinfo/US/Eastern')]
        tm.assert_numpy_array_equal(di_test.values, localized.values)

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern'),
                                    'US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_tz_localize_ambiguous_flags(self, tz):
        # November 6, 2011, fall back, repeat 2 AM hour

        # Pass in flags to determine right dst transition
        dr = date_range(datetime(2011, 11, 6, 0), periods=5,
                        freq=pd.offsets.Hour(), tz=tz)
        times = ['11/06/2011 00:00', '11/06/2011 01:00', '11/06/2011 01:00',
                 '11/06/2011 02:00', '11/06/2011 03:00']

        # Test tz_localize
        di = DatetimeIndex(times)
        is_dst = [1, 1, 0, 0, 0]
        localized = di.tz_localize(tz, ambiguous=is_dst)
        tm.assert_index_equal(dr, localized)
        tm.assert_index_equal(dr, DatetimeIndex(times, tz=tz,
                                                ambiguous=is_dst))

        localized = di.tz_localize(tz, ambiguous=np.array(is_dst))
        tm.assert_index_equal(dr, localized)

        localized = di.tz_localize(tz,
                                   ambiguous=np.array(is_dst).astype('bool'))
        tm.assert_index_equal(dr, localized)

        # Test constructor
        localized = DatetimeIndex(times, tz=tz, ambiguous=is_dst)
        tm.assert_index_equal(dr, localized)

        # Test duplicate times where infer_dst fails
        times += times
        di = DatetimeIndex(times)

        # When the sizes are incompatible, make sure error is raised
        pytest.raises(Exception, di.tz_localize, tz, ambiguous=is_dst)

        # When sizes are compatible and there are repeats ('infer' won't work)
        is_dst = np.hstack((is_dst, is_dst))
        localized = di.tz_localize(tz, ambiguous=is_dst)
        dr = dr.append(dr)
        tm.assert_index_equal(dr, localized)

        # When there is no dst transition, nothing special happens
        dr = date_range(datetime(2011, 6, 1, 0), periods=10,
                        freq=pd.offsets.Hour())
        is_dst = np.array([1] * 10)
        localized = dr.tz_localize(tz)
        localized_is_dst = dr.tz_localize(tz, ambiguous=is_dst)
        tm.assert_index_equal(localized, localized_is_dst)

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern')])
    def test_dti_tz_localize_ambiguous_infer_no_transition(self, tz):
        # When there is no dst transition, nothing special happens
        # When there is no dst transition, nothing special happens
        dr = date_range(datetime(2011, 6, 1, 0), periods=10,
                        freq=pd.offsets.Hour())
        localized = dr.tz_localize(tz)
        localized_infer = dr.tz_localize(tz, ambiguous='infer')
        tm.assert_index_equal(localized, localized_infer)
        with tm.assert_produces_warning(FutureWarning):
            localized_infer_old = dr.tz_localize(tz, infer_dst=True)
        tm.assert_index_equal(localized, localized_infer_old)

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern')])
    def test_dti_tz_localize_ambiguous_infer_repeated_hours(self, tz):
        # November 6, 2011, fall back, repeat 2 AM hour
        # With repeated hours, we can infer the transition
        dr = date_range(datetime(2011, 11, 6, 0), periods=5,
                        freq=pd.offsets.Hour(), tz=tz)
        times = ['11/06/2011 00:00', '11/06/2011 01:00', '11/06/2011 01:00',
                 '11/06/2011 02:00', '11/06/2011 03:00']
        naive = DatetimeIndex(times)
        localized = naive.tz_localize(tz, ambiguous='infer')
        tm.assert_index_equal(dr, localized)

        with tm.assert_produces_warning(FutureWarning):
            localized_old = naive.tz_localize(tz, infer_dst=True)

        tm.assert_index_equal(dr, localized_old)
        tm.assert_index_equal(dr, DatetimeIndex(times, tz=tz,
                                                ambiguous='infer'))

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern')])
    def test_dti_tz_localize_ambiguous_infer_no_repeated_hours(self, tz):
        # November 6, 2011, fall back, repeat 2 AM hour
        # With no repeated hours, we cannot infer the transition
        dr = date_range(datetime(2011, 11, 6, 0), periods=5,
                        freq=pd.offsets.Hour())
        with pytest.raises(pytz.AmbiguousTimeError):
            dr.tz_localize(tz)

    @pytest.mark.parametrize('tz', [pytz.timezone('US/Eastern'),
                                    gettz('US/Eastern')])
    def test_dti_tz_localize_ambiguous_times(self, tz):
        # March 13, 2011, spring forward, skip from 2 AM to 3 AM
        dr = date_range(datetime(2011, 3, 13, 1, 30), periods=3,
                        freq=pd.offsets.Hour())
        with pytest.raises(pytz.NonExistentTimeError):
            dr.tz_localize(tz)

        # after dst transition, it works
        dr = date_range(datetime(2011, 3, 13, 3, 30), periods=3,
                        freq=pd.offsets.Hour(), tz=tz)

        # November 6, 2011, fall back, repeat 2 AM hour
        dr = date_range(datetime(2011, 11, 6, 1, 30), periods=3,
                        freq=pd.offsets.Hour())
        with pytest.raises(pytz.AmbiguousTimeError):
            dr.tz_localize(tz)

        # UTC is OK
        dr = date_range(datetime(2011, 3, 13), periods=48,
                        freq=pd.offsets.Minute(30), tz=pytz.utc)

    # ----------------------------------------------------------------
    # DatetimeIndex.normalize

    def test_dti_normalize_tz(self):
        rng = date_range('1/1/2000 9:30', periods=10, freq='D',
                         tz='US/Eastern')

        result = rng.normalize()
        expected = date_range('1/1/2000', periods=10, freq='D',
                              tz='US/Eastern')
        tm.assert_index_equal(result, expected)

        assert result.is_normalized
        assert not rng.is_normalized

        rng = date_range('1/1/2000 9:30', periods=10, freq='D', tz='UTC')

        result = rng.normalize()
        expected = date_range('1/1/2000', periods=10, freq='D', tz='UTC')
        tm.assert_index_equal(result, expected)

        assert result.is_normalized
        assert not rng.is_normalized

        rng = date_range('1/1/2000 9:30', periods=10, freq='D', tz=tzlocal())
        result = rng.normalize()
        expected = date_range('1/1/2000', periods=10, freq='D', tz=tzlocal())
        tm.assert_index_equal(result, expected)

        assert result.is_normalized
        assert not rng.is_normalized

    @td.skip_if_windows
    @pytest.mark.parametrize('timezone', ['US/Pacific', 'US/Eastern', 'UTC',
                                          'Asia/Kolkata', 'Asia/Shanghai',
                                          'Australia/Canberra'])
    def test_dti_normalize_tz_local(self, timezone):
        # see GH#13459
        with set_timezone(timezone):
            rng = date_range('1/1/2000 9:30', periods=10, freq='D',
                             tz=tzlocal())

            result = rng.normalize()
            expected = date_range('1/1/2000', periods=10, freq='D',
                                  tz=tzlocal())
            tm.assert_index_equal(result, expected)

            assert result.is_normalized
            assert not rng.is_normalized

    # ----------------------------------------------------------------

    def test_dti_drop_dont_lose_tz(self):
        # GH#2621
        ind = date_range("2012-12-01", periods=10, tz="utc")
        ind = ind.drop(ind[-1])

        assert ind.tz is not None

    def test_dti_intersection_utc(self):
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')

        left = rng[10:90][::-1]
        right = rng[20:80][::-1]

        assert left.tz == rng.tz
        result = left.intersection(right)
        assert result.tz == left.tz

    # TODO: Break this into tests for tz_localize and tz_convert
    def test_dti_nat(self):
        # GH#5546
        dates = [NaT]
        idx = DatetimeIndex(dates)
        idx = idx.tz_localize('US/Pacific')
        tm.assert_index_equal(idx, DatetimeIndex(dates, tz='US/Pacific'))
        idx = idx.tz_convert('US/Eastern')
        tm.assert_index_equal(idx, DatetimeIndex(dates, tz='US/Eastern'))
        idx = idx.tz_convert('UTC')
        tm.assert_index_equal(idx, DatetimeIndex(dates, tz='UTC'))

        dates = ['2010-12-01 00:00', '2010-12-02 00:00', NaT]
        idx = DatetimeIndex(dates)
        idx = idx.tz_localize('US/Pacific')
        tm.assert_index_equal(idx, DatetimeIndex(dates, tz='US/Pacific'))
        idx = idx.tz_convert('US/Eastern')
        expected = ['2010-12-01 03:00', '2010-12-02 03:00', NaT]
        tm.assert_index_equal(idx, DatetimeIndex(expected, tz='US/Eastern'))

        idx = idx + pd.offsets.Hour(5)
        expected = ['2010-12-01 08:00', '2010-12-02 08:00', NaT]
        tm.assert_index_equal(idx, DatetimeIndex(expected, tz='US/Eastern'))
        idx = idx.tz_convert('US/Pacific')
        expected = ['2010-12-01 05:00', '2010-12-02 05:00', NaT]
        tm.assert_index_equal(idx, DatetimeIndex(expected, tz='US/Pacific'))

        idx = idx + np.timedelta64(3, 'h')
        expected = ['2010-12-01 08:00', '2010-12-02 08:00', NaT]
        tm.assert_index_equal(idx, DatetimeIndex(expected, tz='US/Pacific'))

        idx = idx.tz_convert('US/Eastern')
        expected = ['2010-12-01 11:00', '2010-12-02 11:00', NaT]
        tm.assert_index_equal(idx, DatetimeIndex(expected, tz='US/Eastern'))

    def test_dti_equals_with_tz(self):
        left = date_range('1/1/2011', periods=100, freq='H', tz='utc')
        right = date_range('1/1/2011', periods=100, freq='H', tz='US/Eastern')

        assert not left.equals(right)

    def test_dti_add_offset_tzaware(self):
        dates = date_range('2012-11-01', periods=3, tz='US/Pacific')
        offset = dates + pd.offsets.Hour(5)
        assert dates[0] + pd.offsets.Hour(5) == offset[0]

        # GH 6818
        for tz in ['UTC', 'US/Pacific', 'Asia/Tokyo']:
            dates = date_range('2010-11-01 00:00', periods=3, tz=tz, freq='H')
            expected = DatetimeIndex(['2010-11-01 05:00', '2010-11-01 06:00',
                                      '2010-11-01 07:00'], freq='H', tz=tz)

            offset = dates + pd.offsets.Hour(5)
            tm.assert_index_equal(offset, expected)
            offset = dates + np.timedelta64(5, 'h')
            tm.assert_index_equal(offset, expected)
            offset = dates + timedelta(hours=5)
            tm.assert_index_equal(offset, expected)

    def test_dti_union_tzaware(self):
        # non-overlapping
        rng = date_range("2012-11-15 00:00:00", periods=6, freq="H",
                         tz="US/Central")

        rng2 = date_range("2012-11-15 12:00:00", periods=6, freq="H",
                          tz="US/Eastern")

        result = rng.union(rng2)
        assert result.tz.zone == 'UTC'

    def test_dti_equality_different_timezones(self):
        utc_range = date_range('1/1/2000', periods=20, tz='UTC')
        eastern_range = utc_range.tz_convert('US/Eastern')
        berlin_range = utc_range.tz_convert('Europe/Berlin')

        for a, b, c in zip(utc_range, eastern_range, berlin_range):
            assert a == b
            assert b == c
            assert a == c

        assert (utc_range == eastern_range).all()
        assert (utc_range == berlin_range).all()
        assert (berlin_range == eastern_range).all()

    @pytest.mark.parametrize('tz', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_construction_ambiguous_endpoint(self, tz):
        # construction with an ambiguous end-point
        # GH#11626

        with pytest.raises(pytz.AmbiguousTimeError):
            date_range("2013-10-26 23:00", "2013-10-27 01:00",
                       tz="Europe/London", freq="H")

        times = date_range("2013-10-26 23:00", "2013-10-27 01:00", freq="H",
                           tz=tz, ambiguous='infer')
        assert times[0] == Timestamp('2013-10-26 23:00', tz=tz, freq="H")

        if str(tz).startswith('dateutil'):
            if LooseVersion(dateutil.__version__) < LooseVersion('2.6.0'):
                # see gh-14621
                assert times[-1] == Timestamp('2013-10-27 01:00:00+0000',
                                              tz=tz, freq="H")
            elif LooseVersion(dateutil.__version__) > LooseVersion('2.6.0'):
                # fixed ambiguous behavior
                assert times[-1] == Timestamp('2013-10-27 01:00:00+0100',
                                              tz=tz, freq="H")
        else:
            assert times[-1] == Timestamp('2013-10-27 01:00:00+0000',
                                          tz=tz, freq="H")

    @pytest.mark.parametrize('tz', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_construction_tz(self, tz):
        """ Test different DatetimeIndex constructions with timezone
        Follow-up of GH#4229
        """

        arr = ['11/10/2005 08:00:00', '11/10/2005 09:00:00']

        idx1 = pd.to_datetime(arr).tz_localize(tz)
        idx2 = DatetimeIndex(start="2005-11-10 08:00:00", freq='H', periods=2,
                             tz=tz)
        idx3 = DatetimeIndex(arr, tz=tz)
        idx4 = DatetimeIndex(np.array(arr), tz=tz)

        for other in [idx2, idx3, idx4]:
            tm.assert_index_equal(idx1, other)

    @pytest.mark.parametrize('tz', ['US/Eastern', 'dateutil/US/Eastern'])
    def test_dti_take_dont_lose_meta(self, tz):
        rng = date_range('1/1/2000', periods=20, tz=tz)

        result = rng.take(lrange(5))
        assert result.tz == rng.tz
        assert result.freq == rng.freq
