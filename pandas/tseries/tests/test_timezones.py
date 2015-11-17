# pylint: disable-msg=E1101,W0612
from datetime import datetime, timedelta, tzinfo, date
import nose

import numpy as np
import pytz

from pandas import (Index, Series, DataFrame, isnull, Timestamp)

from pandas import DatetimeIndex, to_datetime, NaT
from pandas import tslib

import pandas.core.datetools as datetools
import pandas.tseries.offsets as offsets
from pandas.tseries.index import bdate_range, date_range
import pandas.tseries.tools as tools
from pytz import NonExistentTimeError

import pandas.util.testing as tm
from pandas.core.dtypes import DatetimeTZDtype
from pandas.util.testing import assert_frame_equal
from pandas.compat import lrange, zip


try:
    import pytz
except ImportError:
    pass

try:
    import dateutil
except ImportError:
    pass


class FixedOffset(tzinfo):
    """Fixed offset in minutes east from UTC."""

    def __init__(self, offset, name):
        self.__offset = timedelta(minutes=offset)
        self.__name = name

    def utcoffset(self, dt):
        return self.__offset

    def tzname(self, dt):
        return self.__name

    def dst(self, dt):
        return timedelta(0)

fixed_off = FixedOffset(-420, '-07:00')
fixed_off_no_name = FixedOffset(-330, None)


class TestTimeZoneSupportPytz(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        tm._skip_if_no_pytz()

    def tz(self, tz):
        ''' Construct a timezone object from a string. Overridden in subclass to parameterize tests. '''
        return pytz.timezone(tz)

    def tzstr(self, tz):
        ''' Construct a timezone string from a string. Overridden in subclass to parameterize tests. '''
        return tz

    def localize(self, tz, x):
        return tz.localize(x)

    def cmptz(self, tz1, tz2):
        ''' Compare two timezones. Overridden in subclass to parameterize tests. '''
        return tz1.zone == tz2.zone

    def test_utc_to_local_no_modify(self):
        rng = date_range('3/11/2012', '3/12/2012', freq='H', tz='utc')
        rng_eastern = rng.tz_convert(self.tzstr('US/Eastern'))

        # Values are unmodified
        self.assertTrue(np.array_equal(rng.asi8, rng_eastern.asi8))

        self.assertTrue(self.cmptz(rng_eastern.tz, self.tz('US/Eastern')))

    def test_utc_to_local_no_modify_explicit(self):
        rng = date_range('3/11/2012', '3/12/2012', freq='H', tz='utc')
        rng_eastern = rng.tz_convert(self.tz('US/Eastern'))

        # Values are unmodified
        self.assert_numpy_array_equal(rng.asi8, rng_eastern.asi8)

        self.assertEqual(rng_eastern.tz, self.tz('US/Eastern'))


    def test_localize_utc_conversion(self):
        # Localizing to time zone should:
        #  1) check for DST ambiguities
        #  2) convert to UTC

        rng = date_range('3/10/2012', '3/11/2012', freq='30T')

        converted = rng.tz_localize(self.tzstr('US/Eastern'))
        expected_naive = rng + offsets.Hour(5)
        self.assert_numpy_array_equal(converted.asi8, expected_naive.asi8)

        # DST ambiguity, this should fail
        rng = date_range('3/11/2012', '3/12/2012', freq='30T')
        # Is this really how it should fail??
        self.assertRaises(NonExistentTimeError, rng.tz_localize, self.tzstr('US/Eastern'))

    def test_localize_utc_conversion_explicit(self):
        # Localizing to time zone should:
        #  1) check for DST ambiguities
        #  2) convert to UTC

        rng = date_range('3/10/2012', '3/11/2012', freq='30T')
        converted = rng.tz_localize(self.tz('US/Eastern'))
        expected_naive = rng + offsets.Hour(5)
        self.assertTrue(np.array_equal(converted.asi8, expected_naive.asi8))

        # DST ambiguity, this should fail
        rng = date_range('3/11/2012', '3/12/2012', freq='30T')
        # Is this really how it should fail??
        self.assertRaises(NonExistentTimeError, rng.tz_localize, self.tz('US/Eastern'))

    def test_timestamp_tz_localize(self):
        stamp = Timestamp('3/11/2012 04:00')

        result = stamp.tz_localize(self.tzstr('US/Eastern'))
        expected = Timestamp('3/11/2012 04:00', tz=self.tzstr('US/Eastern'))
        self.assertEqual(result.hour, expected.hour)
        self.assertEqual(result, expected)

    def test_timestamp_tz_localize_explicit(self):
        stamp = Timestamp('3/11/2012 04:00')

        result = stamp.tz_localize(self.tz('US/Eastern'))
        expected = Timestamp('3/11/2012 04:00', tz=self.tz('US/Eastern'))
        self.assertEqual(result.hour, expected.hour)
        self.assertEqual(result, expected)

    def test_timestamp_constructed_by_date_and_tz(self):
        # Fix Issue 2993, Timestamp cannot be constructed by datetime.date
        # and tz correctly

        result = Timestamp(date(2012, 3, 11), tz=self.tzstr('US/Eastern'))

        expected = Timestamp('3/11/2012', tz=self.tzstr('US/Eastern'))
        self.assertEqual(result.hour, expected.hour)
        self.assertEqual(result, expected)

    def test_timestamp_constructed_by_date_and_tz_explicit(self):
        # Fix Issue 2993, Timestamp cannot be constructed by datetime.date
        # and tz correctly

        result = Timestamp(date(2012, 3, 11), tz=self.tz('US/Eastern'))

        expected = Timestamp('3/11/2012', tz=self.tz('US/Eastern'))
        self.assertEqual(result.hour, expected.hour)
        self.assertEqual(result, expected)

    def test_timestamp_to_datetime_tzoffset(self):
        # tzoffset
        from dateutil.tz import tzoffset
        tzinfo = tzoffset(None, 7200)
        expected = Timestamp('3/11/2012 04:00', tz=tzinfo)
        result = Timestamp(expected.to_datetime())
        self.assertEqual(expected, result)

    def test_timedelta_push_over_dst_boundary(self):
        # #1389

        # 4 hours before DST transition
        stamp = Timestamp('3/10/2012 22:00', tz=self.tzstr('US/Eastern'))

        result = stamp + timedelta(hours=6)

        # spring forward, + "7" hours
        expected = Timestamp('3/11/2012 05:00', tz=self.tzstr('US/Eastern'))

        self.assertEqual(result, expected)

    def test_timedelta_push_over_dst_boundary_explicit(self):
        # #1389

        # 4 hours before DST transition
        stamp = Timestamp('3/10/2012 22:00', tz=self.tz('US/Eastern'))

        result = stamp + timedelta(hours=6)

        # spring forward, + "7" hours
        expected = Timestamp('3/11/2012 05:00', tz=self.tz('US/Eastern'))

        self.assertEqual(result, expected)

    def test_tz_localize_dti(self):
        from pandas.tseries.offsets import Hour

        dti = DatetimeIndex(start='1/1/2005', end='1/1/2005 0:00:30.256',
                            freq='L')
        dti2 = dti.tz_localize(self.tzstr('US/Eastern'))

        dti_utc = DatetimeIndex(start='1/1/2005 05:00',
                                end='1/1/2005 5:00:30.256', freq='L',
                                tz='utc')

        self.assert_numpy_array_equal(dti2.values, dti_utc.values)

        dti3 = dti2.tz_convert(self.tzstr('US/Pacific'))
        self.assert_numpy_array_equal(dti3.values, dti_utc.values)

        dti = DatetimeIndex(start='11/6/2011 1:59',
                            end='11/6/2011 2:00', freq='L')
        self.assertRaises(pytz.AmbiguousTimeError, dti.tz_localize,
                          self.tzstr('US/Eastern'))

        dti = DatetimeIndex(start='3/13/2011 1:59', end='3/13/2011 2:00',
                            freq='L')
        self.assertRaises(
            pytz.NonExistentTimeError, dti.tz_localize, self.tzstr('US/Eastern'))

    def test_tz_localize_empty_series(self):
        # #2248

        ts = Series()

        ts2 = ts.tz_localize('utc')
        self.assertTrue(ts2.index.tz == pytz.utc)

        ts2 = ts.tz_localize(self.tzstr('US/Eastern'))
        self.assertTrue(self.cmptz(ts2.index.tz, self.tz('US/Eastern')))

    def test_astimezone(self):
        utc = Timestamp('3/11/2012 22:00', tz='UTC')
        expected = utc.tz_convert(self.tzstr('US/Eastern'))
        result = utc.astimezone(self.tzstr('US/Eastern'))
        self.assertEqual(expected, result)
        tm.assertIsInstance(result, Timestamp)

    def test_create_with_tz(self):
        stamp = Timestamp('3/11/2012 05:00', tz=self.tzstr('US/Eastern'))
        self.assertEqual(stamp.hour, 5)

        rng = date_range(
            '3/11/2012 04:00', periods=10, freq='H', tz=self.tzstr('US/Eastern'))

        self.assertEqual(stamp, rng[1])

        utc_stamp = Timestamp('3/11/2012 05:00', tz='utc')
        self.assertIs(utc_stamp.tzinfo, pytz.utc)
        self.assertEqual(utc_stamp.hour, 5)

        stamp = Timestamp('3/11/2012 05:00').tz_localize('utc')
        self.assertEqual(utc_stamp.hour, 5)

    def test_create_with_fixed_tz(self):
        off = FixedOffset(420, '+07:00')
        start = datetime(2012, 3, 11, 5, 0, 0, tzinfo=off)
        end = datetime(2012, 6, 11, 5, 0, 0, tzinfo=off)
        rng = date_range(start=start, end=end)
        self.assertEqual(off, rng.tz)

        rng2 = date_range(start, periods=len(rng), tz=off)
        self.assertTrue(rng.equals(rng2))

        rng3 = date_range(
            '3/11/2012 05:00:00+07:00', '6/11/2012 05:00:00+07:00')
        self.assertTrue((rng.values == rng3.values).all())

    def test_create_with_fixedoffset_noname(self):
        off = fixed_off_no_name
        start = datetime(2012, 3, 11, 5, 0, 0, tzinfo=off)
        end = datetime(2012, 6, 11, 5, 0, 0, tzinfo=off)
        rng = date_range(start=start, end=end)
        self.assertEqual(off, rng.tz)

        idx = Index([start, end])
        self.assertEqual(off, idx.tz)

    def test_date_range_localize(self):
        rng = date_range(
            '3/11/2012 03:00', periods=15, freq='H', tz='US/Eastern')
        rng2 = DatetimeIndex(['3/11/2012 03:00', '3/11/2012 04:00'],
                             tz='US/Eastern')
        rng3 = date_range('3/11/2012 03:00', periods=15, freq='H')
        rng3 = rng3.tz_localize('US/Eastern')

        self.assertTrue(rng.equals(rng3))

        # DST transition time
        val = rng[0]
        exp = Timestamp('3/11/2012 03:00', tz='US/Eastern')

        self.assertEqual(val.hour, 3)
        self.assertEqual(exp.hour, 3)
        self.assertEqual(val, exp)  # same UTC value
        self.assertTrue(rng[:2].equals(rng2))

        # Right before the DST transition
        rng = date_range(
            '3/11/2012 00:00', periods=2, freq='H', tz='US/Eastern')
        rng2 = DatetimeIndex(['3/11/2012 00:00', '3/11/2012 01:00'],
                             tz='US/Eastern')
        self.assertTrue(rng.equals(rng2))
        exp = Timestamp('3/11/2012 00:00', tz='US/Eastern')
        self.assertEqual(exp.hour, 0)
        self.assertEqual(rng[0], exp)
        exp = Timestamp('3/11/2012 01:00', tz='US/Eastern')
        self.assertEqual(exp.hour, 1)
        self.assertEqual(rng[1], exp)

        rng = date_range('3/11/2012 00:00', periods=10, freq='H',
                         tz='US/Eastern')
        self.assertEqual(rng[2].hour, 3)

    def test_utc_box_timestamp_and_localize(self):
        rng = date_range('3/11/2012', '3/12/2012', freq='H', tz='utc')
        rng_eastern = rng.tz_convert(self.tzstr('US/Eastern'))

        tz = self.tz('US/Eastern')
        expected = rng[-1].astimezone(tz)

        stamp = rng_eastern[-1]
        self.assertEqual(stamp, expected)
        self.assertEqual(stamp.tzinfo, expected.tzinfo)

        # right tzinfo
        rng = date_range('3/13/2012', '3/14/2012', freq='H', tz='utc')
        rng_eastern = rng.tz_convert(self.tzstr('US/Eastern'))
        # test not valid for dateutil timezones.
        # self.assertIn('EDT', repr(rng_eastern[0].tzinfo))
        self.assertTrue('EDT' in repr(rng_eastern[0].tzinfo) or 'tzfile' in repr(rng_eastern[0].tzinfo))

    def test_timestamp_tz_convert(self):
        strdates = ['1/1/2012', '3/1/2012', '4/1/2012']
        idx = DatetimeIndex(strdates, tz=self.tzstr('US/Eastern'))

        conv = idx[0].tz_convert(self.tzstr('US/Pacific'))
        expected = idx.tz_convert(self.tzstr('US/Pacific'))[0]

        self.assertEqual(conv, expected)

    def test_pass_dates_localize_to_utc(self):
        strdates = ['1/1/2012', '3/1/2012', '4/1/2012']

        idx = DatetimeIndex(strdates)
        conv = idx.tz_localize(self.tzstr('US/Eastern'))

        fromdates = DatetimeIndex(strdates, tz=self.tzstr('US/Eastern'))

        self.assertEqual(conv.tz, fromdates.tz)
        self.assert_numpy_array_equal(conv.values, fromdates.values)

    def test_field_access_localize(self):
        strdates = ['1/1/2012', '3/1/2012', '4/1/2012']
        rng = DatetimeIndex(strdates, tz=self.tzstr('US/Eastern'))
        self.assertTrue((rng.hour == 0).all())

        # a more unusual time zone, #1946
        dr = date_range('2011-10-02 00:00', freq='h', periods=10,
                        tz=self.tzstr('America/Atikokan'))

        expected = np.arange(10)
        self.assert_numpy_array_equal(dr.hour, expected)

    def test_with_tz(self):
        tz = self.tz('US/Central')

        # just want it to work
        start = datetime(2011, 3, 12, tzinfo=pytz.utc)
        dr = bdate_range(start, periods=50, freq=datetools.Hour())
        self.assertIs(dr.tz, pytz.utc)

        # DateRange with naive datetimes
        dr = bdate_range('1/1/2005', '1/1/2009', tz=pytz.utc)
        dr = bdate_range('1/1/2005', '1/1/2009', tz=tz)

        # normalized
        central = dr.tz_convert(tz)
        self.assertIs(central.tz, tz)
        comp = self.localize(tz, central[0].to_pydatetime().replace(tzinfo=None)).tzinfo
        self.assertIs(central[0].tz, comp)

        # compare vs a localized tz
        comp = self.localize(tz, dr[0].to_pydatetime().replace(tzinfo=None)).tzinfo
        self.assertIs(central[0].tz, comp)

        # datetimes with tzinfo set
        dr = bdate_range(datetime(2005, 1, 1, tzinfo=pytz.utc),
                         '1/1/2009', tz=pytz.utc)

        self.assertRaises(Exception, bdate_range,
                          datetime(2005, 1, 1, tzinfo=pytz.utc),
                          '1/1/2009', tz=tz)

    def test_tz_localize(self):
        dr = bdate_range('1/1/2009', '1/1/2010')
        dr_utc = bdate_range('1/1/2009', '1/1/2010', tz=pytz.utc)
        localized = dr.tz_localize(pytz.utc)
        self.assert_numpy_array_equal(dr_utc, localized)

    def test_with_tz_ambiguous_times(self):
        tz = self.tz('US/Eastern')

        # March 13, 2011, spring forward, skip from 2 AM to 3 AM
        dr = date_range(datetime(2011, 3, 13, 1, 30), periods=3,
                        freq=datetools.Hour())
        self.assertRaises(pytz.NonExistentTimeError, dr.tz_localize, tz)

        # after dst transition, it works
        dr = date_range(datetime(2011, 3, 13, 3, 30), periods=3,
                        freq=datetools.Hour(), tz=tz)

        # November 6, 2011, fall back, repeat 2 AM hour
        dr = date_range(datetime(2011, 11, 6, 1, 30), periods=3,
                        freq=datetools.Hour())
        self.assertRaises(pytz.AmbiguousTimeError, dr.tz_localize, tz)

        # UTC is OK
        dr = date_range(datetime(2011, 3, 13), periods=48,
                        freq=datetools.Minute(30), tz=pytz.utc)

    def test_ambiguous_infer(self):
        # November 6, 2011, fall back, repeat 2 AM hour
        # With no repeated hours, we cannot infer the transition
        tz = self.tz('US/Eastern')
        dr = date_range(datetime(2011, 11, 6, 0), periods=5,
                        freq=datetools.Hour())
        self.assertRaises(pytz.AmbiguousTimeError, dr.tz_localize, tz)

        # With repeated hours, we can infer the transition
        dr = date_range(datetime(2011, 11, 6, 0), periods=5,
                        freq=datetools.Hour(), tz=tz)
        times = ['11/06/2011 00:00', '11/06/2011 01:00',
                 '11/06/2011 01:00', '11/06/2011 02:00',
                 '11/06/2011 03:00']
        di = DatetimeIndex(times)
        localized = di.tz_localize(tz, ambiguous='infer')
        self.assert_numpy_array_equal(dr, localized)
        with tm.assert_produces_warning(FutureWarning):
            localized_old = di.tz_localize(tz, infer_dst=True)
        self.assert_numpy_array_equal(dr, localized_old)
        self.assert_numpy_array_equal(dr, DatetimeIndex(times, tz=tz, ambiguous='infer'))

        # When there is no dst transition, nothing special happens
        dr = date_range(datetime(2011, 6, 1, 0), periods=10,
                        freq=datetools.Hour())
        localized = dr.tz_localize(tz)
        localized_infer = dr.tz_localize(tz, ambiguous='infer')
        self.assert_numpy_array_equal(localized, localized_infer)
        with tm.assert_produces_warning(FutureWarning):
            localized_infer_old = dr.tz_localize(tz, infer_dst=True)
        self.assert_numpy_array_equal(localized, localized_infer_old)

    def test_ambiguous_flags(self):
        # November 6, 2011, fall back, repeat 2 AM hour
        tz = self.tz('US/Eastern')

        # Pass in flags to determine right dst transition
        dr = date_range(datetime(2011, 11, 6, 0), periods=5,
                        freq=datetools.Hour(), tz=tz)
        times = ['11/06/2011 00:00', '11/06/2011 01:00',
                 '11/06/2011 01:00', '11/06/2011 02:00',
                 '11/06/2011 03:00']

        # Test tz_localize
        di = DatetimeIndex(times)
        is_dst = [1, 1, 0, 0, 0]
        localized = di.tz_localize(tz, ambiguous=is_dst)
        self.assert_numpy_array_equal(dr, localized)
        self.assert_numpy_array_equal(dr, DatetimeIndex(times, tz=tz, ambiguous=is_dst))

        localized = di.tz_localize(tz, ambiguous=np.array(is_dst))
        self.assert_numpy_array_equal(dr, localized)

        localized = di.tz_localize(tz, ambiguous=np.array(is_dst).astype('bool'))
        self.assert_numpy_array_equal(dr, localized)

        # Test constructor
        localized = DatetimeIndex(times, tz=tz, ambiguous=is_dst)
        self.assert_numpy_array_equal(dr, localized)

        # Test duplicate times where infer_dst fails
        times += times
        di = DatetimeIndex(times)

        # When the sizes are incompatible, make sure error is raised
        self.assertRaises(Exception, di.tz_localize, tz, ambiguous=is_dst)

        # When sizes are compatible and there are repeats ('infer' won't work)
        is_dst = np.hstack((is_dst, is_dst))
        localized = di.tz_localize(tz, ambiguous=is_dst)
        dr = dr.append(dr)
        self.assert_numpy_array_equal(dr, localized)

        # When there is no dst transition, nothing special happens
        dr = date_range(datetime(2011, 6, 1, 0), periods=10,
                        freq=datetools.Hour())
        is_dst = np.array([1] * 10)
        localized = dr.tz_localize(tz)
        localized_is_dst = dr.tz_localize(tz, ambiguous=is_dst)
        self.assert_numpy_array_equal(localized, localized_is_dst)

        # construction with an ambiguous end-point
        # GH 11626
        tz=self.tzstr("Europe/London")

        def f():
            date_range("2013-10-26 23:00", "2013-10-27 01:00",
                       tz="Europe/London",
                       freq="H")
            self.assertRaises(pytz.AmbiguousTimeError, f)
        times = date_range("2013-10-26 23:00", "2013-10-27 01:00",
                              freq="H",
                              tz=tz,
                              ambiguous='infer')
        self.assertEqual(times[0],Timestamp('2013-10-26 23:00',tz=tz))
        self.assertEqual(times[-1],Timestamp('2013-10-27 01:00',tz=tz))

    def test_ambiguous_nat(self):
        tz = self.tz('US/Eastern')
        times = ['11/06/2011 00:00', '11/06/2011 01:00',
                 '11/06/2011 01:00', '11/06/2011 02:00',
                 '11/06/2011 03:00']
        di = DatetimeIndex(times)
        localized = di.tz_localize(tz, ambiguous='NaT')

        times = ['11/06/2011 00:00', np.NaN,
                 np.NaN, '11/06/2011 02:00',
                 '11/06/2011 03:00']
        di_test = DatetimeIndex(times, tz='US/Eastern')
        self.assert_numpy_array_equal(di_test, localized)

    # test utility methods
    def test_infer_tz(self):
        eastern = self.tz('US/Eastern')
        utc = pytz.utc

        _start = datetime(2001, 1, 1)
        _end = datetime(2009, 1, 1)

        start = self.localize(eastern, _start)
        end = self.localize(eastern, _end)
        assert(tools._infer_tzinfo(start, end) is self.localize(eastern, _start).tzinfo)
        assert(tools._infer_tzinfo(start, None) is self.localize(eastern, _start).tzinfo)
        assert(tools._infer_tzinfo(None, end) is self.localize(eastern, _end).tzinfo)

        start = utc.localize(_start)
        end = utc.localize(_end)
        assert(tools._infer_tzinfo(start, end) is utc)

        end = self.localize(eastern, _end)
        self.assertRaises(Exception, tools._infer_tzinfo, start, end)
        self.assertRaises(Exception, tools._infer_tzinfo, end, start)

    def test_tz_string(self):
        result = date_range('1/1/2000', periods=10, tz=self.tzstr('US/Eastern'))
        expected = date_range('1/1/2000', periods=10,
                              tz=self.tz('US/Eastern'))

        self.assertTrue(result.equals(expected))

    def test_take_dont_lose_meta(self):
        tm._skip_if_no_pytz()
        rng = date_range('1/1/2000', periods=20, tz=self.tzstr('US/Eastern'))

        result = rng.take(lrange(5))
        self.assertEqual(result.tz, rng.tz)
        self.assertEqual(result.freq, rng.freq)

    def test_index_with_timezone_repr(self):
        rng = date_range('4/13/2010', '5/6/2010')

        rng_eastern = rng.tz_localize(self.tzstr('US/Eastern'))

        rng_repr = repr(rng_eastern)
        self.assertIn('2010-04-13 00:00:00', rng_repr)

    def test_index_astype_asobject_tzinfos(self):
        # #1345

        # dates around a dst transition
        rng = date_range('2/13/2010', '5/6/2010', tz=self.tzstr('US/Eastern'))

        objs = rng.asobject
        for i, x in enumerate(objs):
            exval = rng[i]
            self.assertEqual(x, exval)
            self.assertEqual(x.tzinfo, exval.tzinfo)

        objs = rng.astype(object)
        for i, x in enumerate(objs):
            exval = rng[i]
            self.assertEqual(x, exval)
            self.assertEqual(x.tzinfo, exval.tzinfo)

    def test_localized_at_time_between_time(self):
        from datetime import time

        rng = date_range('4/16/2012', '5/1/2012', freq='H')
        ts = Series(np.random.randn(len(rng)), index=rng)

        ts_local = ts.tz_localize(self.tzstr('US/Eastern'))

        result = ts_local.at_time(time(10, 0))
        expected = ts.at_time(time(10, 0)).tz_localize(self.tzstr('US/Eastern'))
        tm.assert_series_equal(result, expected)
        self.assertTrue(self.cmptz(result.index.tz, self.tz('US/Eastern')))

        t1, t2 = time(10, 0), time(11, 0)
        result = ts_local.between_time(t1, t2)
        expected = ts.between_time(t1, t2).tz_localize(self.tzstr('US/Eastern'))
        tm.assert_series_equal(result, expected)
        self.assertTrue(self.cmptz(result.index.tz, self.tz('US/Eastern')))

    def test_string_index_alias_tz_aware(self):
        rng = date_range('1/1/2000', periods=10, tz=self.tzstr('US/Eastern'))
        ts = Series(np.random.randn(len(rng)), index=rng)

        result = ts['1/3/2000']
        self.assertAlmostEqual(result, ts[2])

    def test_fixed_offset(self):
        dates = [datetime(2000, 1, 1, tzinfo=fixed_off),
                 datetime(2000, 1, 2, tzinfo=fixed_off),
                 datetime(2000, 1, 3, tzinfo=fixed_off)]
        result = to_datetime(dates)
        self.assertEqual(result.tz, fixed_off)

    def test_fixedtz_topydatetime(self):
        dates = np.array([datetime(2000, 1, 1, tzinfo=fixed_off),
                          datetime(2000, 1, 2, tzinfo=fixed_off),
                          datetime(2000, 1, 3, tzinfo=fixed_off)])
        result = to_datetime(dates).to_pydatetime()
        self.assert_numpy_array_equal(dates, result)
        result = to_datetime(dates)._mpl_repr()
        self.assert_numpy_array_equal(dates, result)

    def test_convert_tz_aware_datetime_datetime(self):
        # #1581

        tz = self.tz('US/Eastern')

        dates = [datetime(2000, 1, 1), datetime(2000, 1, 2),
                 datetime(2000, 1, 3)]

        dates_aware = [self.localize(tz, x) for x in dates]
        result = to_datetime(dates_aware)
        self.assertTrue(self.cmptz(result.tz, self.tz('US/Eastern')))

        converted = to_datetime(dates_aware, utc=True)
        ex_vals = [Timestamp(x).value for x in dates_aware]
        self.assert_numpy_array_equal(converted.asi8, ex_vals)
        self.assertIs(converted.tz, pytz.utc)

    def test_to_datetime_utc(self):
        from dateutil.parser import parse
        arr = np.array([parse('2012-06-13T01:39:00Z')], dtype=object)

        result = to_datetime(arr, utc=True)
        self.assertIs(result.tz, pytz.utc)

    def test_to_datetime_tzlocal(self):
        from dateutil.parser import parse
        from dateutil.tz import tzlocal
        dt = parse('2012-06-13T01:39:00Z')
        dt = dt.replace(tzinfo=tzlocal())

        arr = np.array([dt], dtype=object)

        result = to_datetime(arr, utc=True)
        self.assertIs(result.tz, pytz.utc)

        rng = date_range('2012-11-03 03:00', '2012-11-05 03:00', tz=tzlocal())
        arr = rng.to_pydatetime()
        result = to_datetime(arr, utc=True)
        self.assertIs(result.tz, pytz.utc)

    def test_frame_no_datetime64_dtype(self):

        # after 7822
        # these retain the timezones on dict construction

        dr = date_range('2011/1/1', '2012/1/1', freq='W-FRI')
        dr_tz = dr.tz_localize(self.tzstr('US/Eastern'))
        e = DataFrame({'A': 'foo', 'B': dr_tz}, index=dr)
        tz_expected = DatetimeTZDtype('ns',dr_tz.tzinfo)
        self.assertEqual(e['B'].dtype, tz_expected)

        # GH 2810 (with timezones)
        datetimes_naive   = [ ts.to_pydatetime() for ts in dr ]
        datetimes_with_tz = [ ts.to_pydatetime() for ts in dr_tz ]
        df = DataFrame({'dr' : dr, 'dr_tz' : dr_tz,
                        'datetimes_naive': datetimes_naive,
                        'datetimes_with_tz' : datetimes_with_tz })
        result = df.get_dtype_counts().sort_index()
        expected = Series({ 'datetime64[ns]' : 2, str(tz_expected) : 2 }).sort_index()
        tm.assert_series_equal(result, expected)

    def test_hongkong_tz_convert(self):
        # #1673
        dr = date_range(
            '2012-01-01', '2012-01-10', freq='D', tz='Hongkong')

        # it works!
        dr.hour

    def test_tz_convert_unsorted(self):
        dr = date_range('2012-03-09', freq='H', periods=100, tz='utc')
        dr = dr.tz_convert(self.tzstr('US/Eastern'))

        result = dr[::-1].hour
        exp = dr.hour[::-1]
        tm.assert_almost_equal(result, exp)

    def test_shift_localized(self):
        dr = date_range('2011/1/1', '2012/1/1', freq='W-FRI')
        dr_tz = dr.tz_localize(self.tzstr('US/Eastern'))

        result = dr_tz.shift(1, '10T')
        self.assertEqual(result.tz, dr_tz.tz)

    def test_tz_aware_asfreq(self):
        dr = date_range(
            '2011-12-01', '2012-07-20', freq='D', tz=self.tzstr('US/Eastern'))

        s = Series(np.random.randn(len(dr)), index=dr)

        # it works!
        s.asfreq('T')

    def test_static_tzinfo(self):
        # it works!
        index = DatetimeIndex([datetime(2012, 1, 1)], tz=self.tzstr('EST'))
        index.hour
        index[0]

    def test_tzaware_datetime_to_index(self):
        d = [datetime(2012, 8, 19, tzinfo=self.tz('US/Eastern'))]

        index = DatetimeIndex(d)
        self.assertTrue(self.cmptz(index.tz, self.tz('US/Eastern')))

    def test_date_range_span_dst_transition(self):
        # #1778

        # Standard -> Daylight Savings Time
        dr = date_range('03/06/2012 00:00', periods=200, freq='W-FRI',
                        tz='US/Eastern')

        self.assertTrue((dr.hour == 0).all())

        dr = date_range('2012-11-02', periods=10, tz=self.tzstr('US/Eastern'))
        self.assertTrue((dr.hour == 0).all())

    def test_convert_datetime_list(self):
        dr = date_range('2012-06-02', periods=10, tz=self.tzstr('US/Eastern'))

        dr2 = DatetimeIndex(list(dr), name='foo')
        self.assertTrue(dr.equals(dr2))
        self.assertEqual(dr.tz, dr2.tz)
        self.assertEqual(dr2.name, 'foo')

    def test_frame_from_records_utc(self):
        rec = {'datum': 1.5,
               'begin_time': datetime(2006, 4, 27, tzinfo=pytz.utc)}

        # it works
        DataFrame.from_records([rec], index='begin_time')

    def test_frame_reset_index(self):
        dr = date_range('2012-06-02', periods=10, tz=self.tzstr('US/Eastern'))
        df = DataFrame(np.random.randn(len(dr)), dr)
        roundtripped = df.reset_index().set_index('index')
        xp = df.index.tz
        rs = roundtripped.index.tz
        self.assertEqual(xp, rs)

    def test_dateutil_tzoffset_support(self):
        from dateutil.tz import tzoffset
        values = [188.5, 328.25]
        tzinfo = tzoffset(None, 7200)
        index = [datetime(2012, 5, 11, 11, tzinfo=tzinfo),
                 datetime(2012, 5, 11, 12, tzinfo=tzinfo)]
        series = Series(data=values, index=index)

        self.assertEqual(series.index.tz, tzinfo)

        # it works! #2443
        repr(series.index[0])

    def test_getitem_pydatetime_tz(self):
        index = date_range(start='2012-12-24 16:00',
                           end='2012-12-24 18:00', freq='H',
                           tz=self.tzstr('Europe/Berlin'))
        ts = Series(index=index, data=index.hour)
        time_pandas = Timestamp('2012-12-24 17:00', tz=self.tzstr('Europe/Berlin'))
        time_datetime = self.localize(self.tz('Europe/Berlin'), datetime(2012, 12, 24, 17, 0))
        self.assertEqual(ts[time_pandas], ts[time_datetime])

    def test_index_drop_dont_lose_tz(self):
        # #2621
        ind = date_range("2012-12-01", periods=10, tz="utc")
        ind = ind.drop(ind[-1])

        self.assertTrue(ind.tz is not None)

    def test_datetimeindex_tz(self):
        """ Test different DatetimeIndex constructions with timezone
        Follow-up of #4229
        """

        arr = ['11/10/2005 08:00:00', '11/10/2005 09:00:00']

        idx1 = to_datetime(arr).tz_localize(self.tzstr('US/Eastern'))
        idx2 = DatetimeIndex(start="2005-11-10 08:00:00", freq='H', periods=2, tz=self.tzstr('US/Eastern'))
        idx3 = DatetimeIndex(arr, tz=self.tzstr('US/Eastern'))
        idx4 = DatetimeIndex(np.array(arr), tz=self.tzstr('US/Eastern'))

        for other in [idx2, idx3, idx4]:
            self.assertTrue(idx1.equals(other))

    def test_datetimeindex_tz_nat(self):
        idx = to_datetime([Timestamp("2013-1-1", tz=self.tzstr('US/Eastern')), NaT])

        self.assertTrue(isnull(idx[1]))
        self.assertTrue(idx[0].tzinfo is not None)


class TestTimeZoneSupportDateutil(TestTimeZoneSupportPytz):
    _multiprocess_can_split_ = True

    def setUp(self):
        tm._skip_if_no_dateutil()

    def tz(self, tz):
        '''
        Construct a dateutil timezone.
        Use tslib.maybe_get_tz so that we get the filename on the tz right
        on windows. See #7337.
        '''
        return tslib.maybe_get_tz('dateutil/' + tz)

    def tzstr(self, tz):
        ''' Construct a timezone string from a string. Overridden in subclass to parameterize tests. '''
        return 'dateutil/' + tz

    def cmptz(self, tz1, tz2):
        ''' Compare two timezones. Overridden in subclass to parameterize tests. '''
        return tz1 == tz2

    def localize(self, tz, x):
        return x.replace(tzinfo=tz)

    def test_utc_with_system_utc(self):
        # Skipped on win32 due to dateutil bug
        tm._skip_if_windows()

        from pandas.tslib import maybe_get_tz

        # from system utc to real utc
        ts = Timestamp('2001-01-05 11:56', tz=maybe_get_tz('dateutil/UTC'))
        # check that the time hasn't changed.
        self.assertEqual(ts, ts.tz_convert(dateutil.tz.tzutc()))

        # from system utc to real utc
        ts = Timestamp('2001-01-05 11:56', tz=maybe_get_tz('dateutil/UTC'))
        # check that the time hasn't changed.
        self.assertEqual(ts, ts.tz_convert(dateutil.tz.tzutc()))

    def test_tslib_tz_convert_trans_pos_plus_1__bug(self):
        # Regression test for tslib.tz_convert(vals, tz1, tz2).
        # See https://github.com/pydata/pandas/issues/4496 for details.
        for freq, n in [('H', 1), ('T', 60), ('S', 3600)]:
            idx = date_range(datetime(2011, 3, 26, 23), datetime(2011, 3, 27, 1), freq=freq)
            idx = idx.tz_localize('UTC')
            idx = idx.tz_convert('Europe/Moscow')

            expected = np.repeat(np.array([3, 4, 5]), np.array([n, n, 1]))
            self.assert_numpy_array_equal(idx.hour, expected)

    def test_tslib_tz_convert_dst(self):
        for freq, n in [('H', 1), ('T', 60), ('S', 3600)]:
            # Start DST
            idx = date_range('2014-03-08 23:00', '2014-03-09 09:00', freq=freq, tz='UTC')
            idx = idx.tz_convert('US/Eastern')
            expected = np.repeat(np.array([18, 19, 20, 21, 22, 23, 0, 1, 3, 4, 5]),
                                 np.array([n, n, n, n, n, n, n, n, n, n, 1]))
            self.assert_numpy_array_equal(idx.hour, expected)

            idx = date_range('2014-03-08 18:00', '2014-03-09 05:00', freq=freq, tz='US/Eastern')
            idx = idx.tz_convert('UTC')
            expected = np.repeat(np.array([23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
                                 np.array([n, n, n, n, n, n, n, n, n, n, 1]))
            self.assert_numpy_array_equal(idx.hour, expected)

            # End DST
            idx = date_range('2014-11-01 23:00', '2014-11-02 09:00', freq=freq, tz='UTC')
            idx = idx.tz_convert('US/Eastern')
            expected = np.repeat(np.array([19, 20, 21, 22, 23, 0, 1, 1, 2, 3, 4]),
                                 np.array([n, n, n, n, n, n, n, n, n, n, 1]))
            self.assert_numpy_array_equal(idx.hour, expected)

            idx = date_range('2014-11-01 18:00', '2014-11-02 05:00', freq=freq, tz='US/Eastern')
            idx = idx.tz_convert('UTC')
            expected = np.repeat(np.array([22, 23, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
                                 np.array([n, n, n, n, n, n, n, n, n, n, n, n, 1]))
            self.assert_numpy_array_equal(idx.hour, expected)

        # daily
        # Start DST
        idx = date_range('2014-03-08 00:00', '2014-03-09 00:00', freq='D', tz='UTC')
        idx = idx.tz_convert('US/Eastern')
        self.assert_numpy_array_equal(idx.hour, np.array([19, 19]))

        idx = date_range('2014-03-08 00:00', '2014-03-09 00:00', freq='D', tz='US/Eastern')
        idx = idx.tz_convert('UTC')
        self.assert_numpy_array_equal(idx.hour, np.array([5, 5]))

        # End DST
        idx = date_range('2014-11-01 00:00', '2014-11-02 00:00', freq='D', tz='UTC')
        idx = idx.tz_convert('US/Eastern')
        self.assert_numpy_array_equal(idx.hour, np.array([20, 20]))

        idx = date_range('2014-11-01 00:00', '2014-11-02 000:00', freq='D', tz='US/Eastern')
        idx = idx.tz_convert('UTC')
        self.assert_numpy_array_equal(idx.hour, np.array([4, 4]))


class TestTimeZoneCacheKey(tm.TestCase):
    def test_cache_keys_are_distinct_for_pytz_vs_dateutil(self):
        tzs = pytz.common_timezones
        for tz_name in tzs:
            if tz_name == 'UTC':
                # skip utc as it's a special case in dateutil
                continue
            tz_p = tslib.maybe_get_tz(tz_name)
            tz_d = tslib.maybe_get_tz('dateutil/' + tz_name)
            if tz_d is None:
                # skip timezones that dateutil doesn't know about.
                continue
            self.assertNotEqual(tslib._p_tz_cache_key(tz_p), tslib._p_tz_cache_key(tz_d))


class TestTimeZones(tm.TestCase):
    _multiprocess_can_split_ = True
    timezones = ['UTC', 'Asia/Tokyo', 'US/Eastern', 'dateutil/US/Pacific']

    def setUp(self):
        tm._skip_if_no_pytz()

    def test_index_equals_with_tz(self):
        left = date_range('1/1/2011', periods=100, freq='H', tz='utc')
        right = date_range('1/1/2011', periods=100, freq='H',
                           tz='US/Eastern')

        self.assertFalse(left.equals(right))

    def test_tz_localize_naive(self):
        rng = date_range('1/1/2011', periods=100, freq='H')

        conv = rng.tz_localize('US/Pacific')
        exp = date_range('1/1/2011', periods=100, freq='H', tz='US/Pacific')

        self.assertTrue(conv.equals(exp))

    def test_tz_localize_roundtrip(self):
        for tz in self.timezones:
            idx1 = date_range(start='2014-01-01', end='2014-12-31', freq='M')
            idx2 = date_range(start='2014-01-01', end='2014-12-31', freq='D')
            idx3 = date_range(start='2014-01-01', end='2014-03-01', freq='H')
            idx4 = date_range(start='2014-08-01', end='2014-10-31', freq='T')
            for idx in [idx1, idx2, idx3, idx4]:
                localized = idx.tz_localize(tz)
                expected = date_range(start=idx[0], end=idx[-1], freq=idx.freq, tz=tz)
                tm.assert_index_equal(localized, expected)

                with tm.assertRaises(TypeError):
                    localized.tz_localize(tz)

                reset = localized.tz_localize(None)
                tm.assert_index_equal(reset, idx)
                self.assertTrue(reset.tzinfo is None)

    def test_series_frame_tz_localize(self):

        rng = date_range('1/1/2011', periods=100, freq='H')
        ts = Series(1, index=rng)

        result = ts.tz_localize('utc')
        self.assertEqual(result.index.tz.zone, 'UTC')

        df = DataFrame({'a': 1}, index=rng)
        result = df.tz_localize('utc')
        expected = DataFrame({'a': 1}, rng.tz_localize('UTC'))
        self.assertEqual(result.index.tz.zone, 'UTC')
        assert_frame_equal(result, expected)

        df = df.T
        result = df.tz_localize('utc', axis=1)
        self.assertEqual(result.columns.tz.zone, 'UTC')
        assert_frame_equal(result, expected.T)

        # Can't localize if already tz-aware
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')
        ts = Series(1, index=rng)
        tm.assertRaisesRegexp(TypeError, 'Already tz-aware', ts.tz_localize, 'US/Eastern')

    def test_series_frame_tz_convert(self):
        rng = date_range('1/1/2011', periods=200, freq='D',
                         tz='US/Eastern')
        ts = Series(1, index=rng)

        result = ts.tz_convert('Europe/Berlin')
        self.assertEqual(result.index.tz.zone, 'Europe/Berlin')

        df = DataFrame({'a': 1}, index=rng)
        result = df.tz_convert('Europe/Berlin')
        expected = DataFrame({'a': 1}, rng.tz_convert('Europe/Berlin'))
        self.assertEqual(result.index.tz.zone, 'Europe/Berlin')
        assert_frame_equal(result, expected)

        df = df.T
        result = df.tz_convert('Europe/Berlin', axis=1)
        self.assertEqual(result.columns.tz.zone, 'Europe/Berlin')
        assert_frame_equal(result, expected.T)

        # can't convert tz-naive
        rng = date_range('1/1/2011', periods=200, freq='D')
        ts = Series(1, index=rng)
        tm.assertRaisesRegexp(TypeError, "Cannot convert tz-naive", ts.tz_convert, 'US/Eastern')

    def test_tz_convert_roundtrip(self):
        for tz in self.timezones:
            idx1 = date_range(start='2014-01-01', end='2014-12-31', freq='M', tz='UTC')
            exp1 = date_range(start='2014-01-01', end='2014-12-31', freq='M')

            idx2 = date_range(start='2014-01-01', end='2014-12-31', freq='D', tz='UTC')
            exp2 = date_range(start='2014-01-01', end='2014-12-31', freq='D')

            idx3 = date_range(start='2014-01-01', end='2014-03-01', freq='H', tz='UTC')
            exp3 = date_range(start='2014-01-01', end='2014-03-01', freq='H')

            idx4 = date_range(start='2014-08-01', end='2014-10-31', freq='T', tz='UTC')
            exp4 = date_range(start='2014-08-01', end='2014-10-31', freq='T')


            for idx, expected in [(idx1, exp1), (idx2, exp2), (idx3, exp3), (idx4, exp4)]:
                converted = idx.tz_convert(tz)
                reset = converted.tz_convert(None)
                tm.assert_index_equal(reset, expected)
                self.assertTrue(reset.tzinfo is None)
                tm.assert_index_equal(reset, converted.tz_convert('UTC').tz_localize(None))


    def test_join_utc_convert(self):
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')

        left = rng.tz_convert('US/Eastern')
        right = rng.tz_convert('Europe/Berlin')

        for how in ['inner', 'outer', 'left', 'right']:
            result = left.join(left[:-5], how=how)
            tm.assertIsInstance(result, DatetimeIndex)
            self.assertEqual(result.tz, left.tz)

            result = left.join(right[:-5], how=how)
            tm.assertIsInstance(result, DatetimeIndex)
            self.assertEqual(result.tz.zone, 'UTC')

    def test_join_aware(self):
        rng = date_range('1/1/2011', periods=10, freq='H')
        ts = Series(np.random.randn(len(rng)), index=rng)

        ts_utc = ts.tz_localize('utc')

        self.assertRaises(Exception, ts.__add__, ts_utc)
        self.assertRaises(Exception, ts_utc.__add__, ts)

        test1 = DataFrame(np.zeros((6, 3)),
                          index=date_range("2012-11-15 00:00:00", periods=6,
                                           freq="100L", tz="US/Central"))
        test2 = DataFrame(np.zeros((3, 3)),
                          index=date_range("2012-11-15 00:00:00", periods=3,
                                           freq="250L", tz="US/Central"),
                          columns=lrange(3, 6))

        result = test1.join(test2, how='outer')
        ex_index = test1.index.union(test2.index)

        self.assertTrue(result.index.equals(ex_index))
        self.assertTrue(result.index.tz.zone == 'US/Central')

        # non-overlapping
        rng = date_range("2012-11-15 00:00:00", periods=6,
                         freq="H", tz="US/Central")

        rng2 = date_range("2012-11-15 12:00:00", periods=6,
                          freq="H", tz="US/Eastern")

        result = rng.union(rng2)
        self.assertTrue(result.tz.zone == 'UTC')

    def test_align_aware(self):
        idx1 = date_range('2001', periods=5, freq='H', tz='US/Eastern')
        idx2 = date_range('2001', periods=5, freq='2H', tz='US/Eastern')
        df1 = DataFrame(np.random.randn(len(idx1), 3), idx1)
        df2 = DataFrame(np.random.randn(len(idx2), 3), idx2)
        new1, new2 = df1.align(df2)
        self.assertEqual(df1.index.tz, new1.index.tz)
        self.assertEqual(df2.index.tz, new2.index.tz)

    def test_append_aware(self):
        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H',
                          tz='US/Eastern')
        rng2 = date_range('1/1/2011 02:00', periods=1, freq='H',
                          tz='US/Eastern')
        ts1 = Series(np.random.randn(len(rng1)), index=rng1)
        ts2 = Series(np.random.randn(len(rng2)), index=rng2)
        ts_result = ts1.append(ts2)
        self.assertEqual(ts_result.index.tz, rng1.tz)

        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H',
                          tz='UTC')
        rng2 = date_range('1/1/2011 02:00', periods=1, freq='H',
                          tz='UTC')
        ts1 = Series(np.random.randn(len(rng1)), index=rng1)
        ts2 = Series(np.random.randn(len(rng2)), index=rng2)
        ts_result = ts1.append(ts2)
        utc = rng1.tz
        self.assertEqual(utc, ts_result.index.tz)

        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H',
                          tz='US/Eastern')
        rng2 = date_range('1/1/2011 02:00', periods=1, freq='H',
                          tz='US/Central')
        ts1 = Series(np.random.randn(len(rng1)), index=rng1)
        ts2 = Series(np.random.randn(len(rng2)), index=rng2)
        ts_result = ts1.append(ts2)
        self.assertEqual(utc, ts_result.index.tz)

    def test_append_aware_naive(self):
        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H')
        rng2 = date_range('1/1/2011 02:00', periods=1, freq='H',
                          tz='US/Eastern')
        ts1 = Series(np.random.randn(len(rng1)), index=rng1)
        ts2 = Series(np.random.randn(len(rng2)), index=rng2)
        ts_result = ts1.append(ts2)
        self.assertTrue(ts_result.index.equals(
            ts1.index.asobject.append(ts2.index.asobject)))

        # mixed

        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H')
        rng2 = lrange(100)
        ts1 = Series(np.random.randn(len(rng1)), index=rng1)
        ts2 = Series(np.random.randn(len(rng2)), index=rng2)
        ts_result = ts1.append(ts2)
        self.assertTrue(ts_result.index.equals(
            ts1.index.asobject.append(ts2.index)))

    def test_equal_join_ensure_utc(self):
        rng = date_range('1/1/2011', periods=10, freq='H', tz='US/Eastern')
        ts = Series(np.random.randn(len(rng)), index=rng)

        ts_moscow = ts.tz_convert('Europe/Moscow')

        result = ts + ts_moscow
        self.assertIs(result.index.tz, pytz.utc)

        result = ts_moscow + ts
        self.assertIs(result.index.tz, pytz.utc)

        df = DataFrame({'a': ts})
        df_moscow = df.tz_convert('Europe/Moscow')
        result = df + df_moscow
        self.assertIs(result.index.tz, pytz.utc)

        result = df_moscow + df
        self.assertIs(result.index.tz, pytz.utc)

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

        self.assertEqual(result.index.tz, pytz.UTC)
        tm.assert_series_equal(result, expected)

    def test_intersection(self):
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')

        left = rng[10:90][::-1]
        right = rng[20:80][::-1]

        self.assertEqual(left.tz, rng.tz)
        result = left.intersection(right)
        self.assertEqual(result.tz, left.tz)

    def test_timestamp_equality_different_timezones(self):
        utc_range = date_range('1/1/2000', periods=20, tz='UTC')
        eastern_range = utc_range.tz_convert('US/Eastern')
        berlin_range = utc_range.tz_convert('Europe/Berlin')

        for a, b, c in zip(utc_range, eastern_range, berlin_range):
            self.assertEqual(a, b)
            self.assertEqual(b, c)
            self.assertEqual(a, c)

        self.assertTrue((utc_range == eastern_range).all())
        self.assertTrue((utc_range == berlin_range).all())
        self.assertTrue((berlin_range == eastern_range).all())

    def test_datetimeindex_tz(self):
        rng = date_range('03/12/2012 00:00', periods=10, freq='W-FRI',
                         tz='US/Eastern')
        rng2 = DatetimeIndex(data=rng, tz='US/Eastern')
        self.assertTrue(rng.equals(rng2))

    def test_normalize_tz(self):
        rng = date_range('1/1/2000 9:30', periods=10, freq='D',
                         tz='US/Eastern')

        result = rng.normalize()
        expected = date_range('1/1/2000', periods=10, freq='D',
                              tz='US/Eastern')
        self.assertTrue(result.equals(expected))

        self.assertTrue(result.is_normalized)
        self.assertFalse(rng.is_normalized)

        rng = date_range('1/1/2000 9:30', periods=10, freq='D',
                         tz='UTC')

        result = rng.normalize()
        expected = date_range('1/1/2000', periods=10, freq='D',
                              tz='UTC')
        self.assertTrue(result.equals(expected))

        self.assertTrue(result.is_normalized)
        self.assertFalse(rng.is_normalized)

        from dateutil.tz import tzlocal
        rng = date_range('1/1/2000 9:30', periods=10, freq='D',
                         tz=tzlocal())
        result = rng.normalize()
        expected = date_range('1/1/2000', periods=10, freq='D',
                              tz=tzlocal())
        self.assertTrue(result.equals(expected))

        self.assertTrue(result.is_normalized)
        self.assertFalse(rng.is_normalized)

    def test_tzaware_offset(self):
        dates = date_range('2012-11-01', periods=3, tz='US/Pacific')
        offset = dates + offsets.Hour(5)
        self.assertEqual(dates[0] + offsets.Hour(5), offset[0])

        # GH 6818
        for tz in ['UTC', 'US/Pacific', 'Asia/Tokyo']:
            dates = date_range('2010-11-01 00:00', periods=3, tz=tz, freq='H')
            expected = DatetimeIndex(['2010-11-01 05:00', '2010-11-01 06:00',
                                      '2010-11-01 07:00'], freq='H', tz=tz)

            offset = dates + offsets.Hour(5)
            self.assertTrue(offset.equals(expected))
            offset = dates + np.timedelta64(5, 'h')
            self.assertTrue(offset.equals(expected))
            offset = dates + timedelta(hours=5)
            self.assertTrue(offset.equals(expected))

    def test_nat(self):
        # GH 5546
        dates = [NaT]
        idx = DatetimeIndex(dates)
        idx = idx.tz_localize('US/Pacific')
        self.assertTrue(idx.equals(DatetimeIndex(dates, tz='US/Pacific')))
        idx = idx.tz_convert('US/Eastern')
        self.assertTrue(idx.equals(DatetimeIndex(dates, tz='US/Eastern')))
        idx = idx.tz_convert('UTC')
        self.assertTrue(idx.equals(DatetimeIndex(dates, tz='UTC')))

        dates = ['2010-12-01 00:00', '2010-12-02 00:00', NaT]
        idx = DatetimeIndex(dates)
        idx = idx.tz_localize('US/Pacific')
        self.assertTrue(idx.equals(DatetimeIndex(dates, tz='US/Pacific')))
        idx = idx.tz_convert('US/Eastern')
        expected = ['2010-12-01 03:00', '2010-12-02 03:00', NaT]
        self.assertTrue(idx.equals(DatetimeIndex(expected, tz='US/Eastern')))

        idx = idx + offsets.Hour(5)
        expected = ['2010-12-01 08:00', '2010-12-02 08:00', NaT]
        self.assertTrue(idx.equals(DatetimeIndex(expected, tz='US/Eastern')))
        idx = idx.tz_convert('US/Pacific')
        expected = ['2010-12-01 05:00', '2010-12-02 05:00', NaT]
        self.assertTrue(idx.equals(DatetimeIndex(expected, tz='US/Pacific')))

        idx = idx + np.timedelta64(3, 'h')
        expected = ['2010-12-01 08:00', '2010-12-02 08:00', NaT]
        self.assertTrue(idx.equals(DatetimeIndex(expected, tz='US/Pacific')))

        idx = idx.tz_convert('US/Eastern')
        expected = ['2010-12-01 11:00', '2010-12-02 11:00', NaT]
        self.assertTrue(idx.equals(DatetimeIndex(expected, tz='US/Eastern')))


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
