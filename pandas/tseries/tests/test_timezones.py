# pylint: disable-msg=E1101,W0612
from datetime import datetime, time, timedelta, tzinfo, date
import sys
import os
import unittest
import nose

import numpy as np
import pytz

from pandas import (Index, Series, TimeSeries, DataFrame, isnull,
                    date_range, Timestamp)

from pandas import DatetimeIndex, Int64Index, to_datetime

from pandas.core.daterange import DateRange
import pandas.core.datetools as datetools
import pandas.tseries.offsets as offsets
from pandas.tseries.index import bdate_range, date_range
import pandas.tseries.tools as tools
from pytz import NonExistentTimeError

from pandas.util.testing import assert_series_equal, assert_almost_equal, assertRaisesRegexp
import pandas.util.testing as tm

import pandas.lib as lib
import cPickle as pickle
import pandas.core.datetools as dt
from numpy.random import rand
from pandas.util.testing import assert_frame_equal
import pandas.util.py3compat as py3compat
from pandas.core.datetools import BDay
import pandas.core.common as com


def _skip_if_no_pytz():
    try:
        import pytz
    except ImportError:
        raise nose.SkipTest

try:
    import pytz
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


class TestTimeZoneSupport(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        _skip_if_no_pytz()

    def test_utc_to_local_no_modify(self):
        rng = date_range('3/11/2012', '3/12/2012', freq='H', tz='utc')
        rng_eastern = rng.tz_convert('US/Eastern')

        # Values are unmodified
        self.assert_(np.array_equal(rng.asi8, rng_eastern.asi8))

        self.assert_(rng_eastern.tz == pytz.timezone('US/Eastern'))


    def test_localize_utc_conversion(self):
        # Localizing to time zone should:
        #  1) check for DST ambiguities
        #  2) convert to UTC

        rng = date_range('3/10/2012', '3/11/2012', freq='30T')

        converted = rng.tz_localize('US/Eastern')
        expected_naive = rng + offsets.Hour(5)
        self.assert_(np.array_equal(converted.asi8, expected_naive.asi8))

        # DST ambiguity, this should fail
        rng = date_range('3/11/2012', '3/12/2012', freq='30T')
        # Is this really how it should fail??
        self.assertRaises(NonExistentTimeError, rng.tz_localize, 'US/Eastern')

    def test_timestamp_tz_localize(self):
        stamp = Timestamp('3/11/2012 04:00')

        result = stamp.tz_localize('US/Eastern')
        expected = Timestamp('3/11/2012 04:00', tz='US/Eastern')
        self.assertEquals(result.hour, expected.hour)
        self.assertEquals(result, expected)

    def test_timestamp_constructed_by_date_and_tz(self):
        """
        Fix Issue 2993, Timestamp cannot be constructed by datetime.date
        and tz correctly
        """

        result = Timestamp(date(2012, 3, 11), tz='US/Eastern')

        expected = Timestamp('3/11/2012', tz='US/Eastern')
        self.assertEquals(result.hour, expected.hour)
        self.assertEquals(result, expected)

    def test_timestamp_to_datetime_tzoffset(self):
        # tzoffset
        from dateutil.tz import tzoffset
        tzinfo = tzoffset(None, 7200)
        expected = Timestamp('3/11/2012 04:00', tz=tzinfo)
        result = Timestamp(expected.to_datetime())
        self.assertEquals(expected, result)

    def test_timedelta_push_over_dst_boundary(self):
        # #1389

        # 4 hours before DST transition
        stamp = Timestamp('3/10/2012 22:00', tz='US/Eastern')

        result = stamp + timedelta(hours=6)

        # spring forward, + "7" hours
        expected = Timestamp('3/11/2012 05:00', tz='US/Eastern')

        self.assertEquals(result, expected)

    def test_tz_localize_dti(self):
        from pandas.tseries.offsets import Hour

        dti = DatetimeIndex(start='1/1/2005', end='1/1/2005 0:00:30.256',
                            freq='L')
        dti2 = dti.tz_localize('US/Eastern')

        dti_utc = DatetimeIndex(start='1/1/2005 05:00',
                                end='1/1/2005 5:00:30.256', freq='L',
                                tz='utc')

        self.assert_(np.array_equal(dti2.values, dti_utc.values))

        dti3 = dti2.tz_convert('US/Pacific')
        self.assert_(np.array_equal(dti3.values, dti_utc.values))

        dti = DatetimeIndex(start='11/6/2011 1:59',
                            end='11/6/2011 2:00', freq='L')
        self.assertRaises(pytz.AmbiguousTimeError, dti.tz_localize,
                          'US/Eastern')

        dti = DatetimeIndex(start='3/13/2011 1:59', end='3/13/2011 2:00',
                            freq='L')
        self.assertRaises(
            pytz.NonExistentTimeError, dti.tz_localize, 'US/Eastern')

    def test_tz_localize_empty_series(self):
        # #2248

        ts = Series()

        ts2 = ts.tz_localize('utc')
        self.assertTrue(ts2.index.tz == pytz.utc)

        ts2 = ts.tz_localize('US/Eastern')
        self.assertTrue(ts2.index.tz == pytz.timezone('US/Eastern'))

    def test_astimezone(self):
        utc = Timestamp('3/11/2012 22:00', tz='UTC')
        expected = utc.tz_convert('US/Eastern')
        result = utc.astimezone('US/Eastern')
        self.assertEquals(expected, result)
        self.assert_(isinstance(result, Timestamp))

    def test_create_with_tz(self):
        stamp = Timestamp('3/11/2012 05:00', tz='US/Eastern')
        self.assertEquals(stamp.hour, 5)

        rng = date_range(
            '3/11/2012 04:00', periods=10, freq='H', tz='US/Eastern')

        self.assertEquals(stamp, rng[1])

        utc_stamp = Timestamp('3/11/2012 05:00', tz='utc')
        self.assert_(utc_stamp.tzinfo is pytz.utc)
        self.assertEquals(utc_stamp.hour, 5)

        stamp = Timestamp('3/11/2012 05:00').tz_localize('utc')
        self.assertEquals(utc_stamp.hour, 5)

    def test_create_with_fixed_tz(self):
        off = FixedOffset(420, '+07:00')
        start = datetime(2012, 3, 11, 5, 0, 0, tzinfo=off)
        end = datetime(2012, 6, 11, 5, 0, 0, tzinfo=off)
        rng = date_range(start=start, end=end)
        self.assertEqual(off, rng.tz)

        rng2 = date_range(start, periods=len(rng), tz=off)
        self.assert_(rng.equals(rng2))

        rng3 = date_range(
            '3/11/2012 05:00:00+07:00', '6/11/2012 05:00:00+07:00')
        self.assert_((rng.values == rng3.values).all())

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

        self.assert_(rng.equals(rng3))

        # DST transition time
        val = rng[0]
        exp = Timestamp('3/11/2012 03:00', tz='US/Eastern')

        self.assertEquals(val.hour, 3)
        self.assertEquals(exp.hour, 3)
        self.assertEquals(val, exp)  # same UTC value
        self.assert_(rng[:2].equals(rng2))

        # Right before the DST transition
        rng = date_range(
            '3/11/2012 00:00', periods=2, freq='H', tz='US/Eastern')
        rng2 = DatetimeIndex(['3/11/2012 00:00', '3/11/2012 01:00'],
                             tz='US/Eastern')
        self.assert_(rng.equals(rng2))
        exp = Timestamp('3/11/2012 00:00', tz='US/Eastern')
        self.assertEquals(exp.hour, 0)
        self.assertEquals(rng[0], exp)
        exp = Timestamp('3/11/2012 01:00', tz='US/Eastern')
        self.assertEquals(exp.hour, 1)
        self.assertEquals(rng[1], exp)

        rng = date_range('3/11/2012 00:00', periods=10, freq='H',
                         tz='US/Eastern')
        self.assert_(rng[2].hour == 3)

    def test_utc_box_timestamp_and_localize(self):
        rng = date_range('3/11/2012', '3/12/2012', freq='H', tz='utc')
        rng_eastern = rng.tz_convert('US/Eastern')

        tz = pytz.timezone('US/Eastern')
        expected = tz.normalize(rng[-1])

        stamp = rng_eastern[-1]
        self.assertEquals(stamp, expected)
        self.assertEquals(stamp.tzinfo, expected.tzinfo)

        # right tzinfo
        rng = date_range('3/13/2012', '3/14/2012', freq='H', tz='utc')
        rng_eastern = rng.tz_convert('US/Eastern')
        self.assert_('EDT' in repr(rng_eastern[0].tzinfo))

    def test_timestamp_tz_convert(self):
        strdates = ['1/1/2012', '3/1/2012', '4/1/2012']
        idx = DatetimeIndex(strdates, tz='US/Eastern')

        conv = idx[0].tz_convert('US/Pacific')
        expected = idx.tz_convert('US/Pacific')[0]

        self.assertEquals(conv, expected)

    def test_pass_dates_localize_to_utc(self):
        strdates = ['1/1/2012', '3/1/2012', '4/1/2012']

        idx = DatetimeIndex(strdates)
        conv = idx.tz_localize('US/Eastern')

        fromdates = DatetimeIndex(strdates, tz='US/Eastern')

        self.assert_(conv.tz == fromdates.tz)
        self.assert_(np.array_equal(conv.values, fromdates.values))

    def test_field_access_localize(self):
        strdates = ['1/1/2012', '3/1/2012', '4/1/2012']
        rng = DatetimeIndex(strdates, tz='US/Eastern')
        self.assert_((rng.hour == 0).all())

        # a more unusual time zone, #1946
        dr = date_range('2011-10-02 00:00', freq='h', periods=10,
                        tz='America/Atikokan')

        expected = np.arange(10)
        self.assert_(np.array_equal(dr.hour, expected))

    def test_with_tz(self):
        tz = pytz.timezone('US/Central')

        # just want it to work
        start = datetime(2011, 3, 12, tzinfo=pytz.utc)
        dr = bdate_range(start, periods=50, freq=datetools.Hour())
        self.assert_(dr.tz is pytz.utc)

        # DateRange with naive datetimes
        dr = bdate_range('1/1/2005', '1/1/2009', tz=pytz.utc)
        dr = bdate_range('1/1/2005', '1/1/2009', tz=tz)

        # normalized
        central = dr.tz_convert(tz)
        self.assert_(central.tz is tz)
        self.assert_(central[0].tz is tz)

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
        self.assert_(np.array_equal(dr_utc, localized))

    def test_with_tz_ambiguous_times(self):
        tz = pytz.timezone('US/Eastern')

        rng = bdate_range(datetime(2009, 1, 1), datetime(2010, 1, 1))

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

    # test utility methods
    def test_infer_tz(self):
        eastern = pytz.timezone('US/Eastern')
        utc = pytz.utc

        _start = datetime(2001, 1, 1)
        _end = datetime(2009, 1, 1)

        start = eastern.localize(_start)
        end = eastern.localize(_end)
        assert(tools._infer_tzinfo(start, end) is eastern)
        assert(tools._infer_tzinfo(start, None) is eastern)
        assert(tools._infer_tzinfo(None, end) is eastern)

        start = utc.localize(_start)
        end = utc.localize(_end)
        assert(tools._infer_tzinfo(start, end) is utc)

        end = eastern.localize(_end)
        self.assertRaises(Exception, tools._infer_tzinfo, start, end)
        self.assertRaises(Exception, tools._infer_tzinfo, end, start)

    def test_tz_string(self):
        result = date_range('1/1/2000', periods=10, tz='US/Eastern')
        expected = date_range('1/1/2000', periods=10,
                              tz=pytz.timezone('US/Eastern'))

        self.assert_(result.equals(expected))

    def test_take_dont_lose_meta(self):
        _skip_if_no_pytz()
        rng = date_range('1/1/2000', periods=20, tz='US/Eastern')

        result = rng.take(range(5))
        self.assert_(result.tz == rng.tz)
        self.assert_(result.freq == rng.freq)

    def test_index_with_timezone_repr(self):
        rng = date_range('4/13/2010', '5/6/2010')

        rng_eastern = rng.tz_localize('US/Eastern')

        rng_repr = repr(rng)
        self.assert_('2010-04-13 00:00:00' in rng_repr)

    def test_index_astype_asobject_tzinfos(self):
        # #1345

        # dates around a dst transition
        rng = date_range('2/13/2010', '5/6/2010', tz='US/Eastern')

        objs = rng.asobject
        for i, x in enumerate(objs):
            exval = rng[i]
            self.assertEquals(x, exval)
            self.assertEquals(x.tzinfo, exval.tzinfo)

        objs = rng.astype(object)
        for i, x in enumerate(objs):
            exval = rng[i]
            self.assertEquals(x, exval)
            self.assertEquals(x.tzinfo, exval.tzinfo)

    def test_localized_at_time_between_time(self):
        from datetime import time

        rng = date_range('4/16/2012', '5/1/2012', freq='H')
        ts = Series(np.random.randn(len(rng)), index=rng)

        ts_local = ts.tz_localize('US/Eastern')

        result = ts_local.at_time(time(10, 0))
        expected = ts.at_time(time(10, 0)).tz_localize('US/Eastern')
        assert_series_equal(result, expected)
        self.assert_(result.index.tz.zone == 'US/Eastern')

        t1, t2 = time(10, 0), time(11, 0)
        result = ts_local.between_time(t1, t2)
        expected = ts.between_time(t1, t2).tz_localize('US/Eastern')
        assert_series_equal(result, expected)
        self.assert_(result.index.tz.zone == 'US/Eastern')

    def test_string_index_alias_tz_aware(self):
        rng = date_range('1/1/2000', periods=10, tz='US/Eastern')
        ts = Series(np.random.randn(len(rng)), index=rng)

        result = ts['1/3/2000']
        self.assertAlmostEqual(result, ts[2])

    def test_fixed_offset(self):
        dates = [datetime(2000, 1, 1, tzinfo=fixed_off),
                 datetime(2000, 1, 2, tzinfo=fixed_off),
                 datetime(2000, 1, 3, tzinfo=fixed_off)]
        result = to_datetime(dates)
        self.assert_(result.tz == fixed_off)

    def test_fixedtz_topydatetime(self):
        dates = np.array([datetime(2000, 1, 1, tzinfo=fixed_off),
                          datetime(2000, 1, 2, tzinfo=fixed_off),
                          datetime(2000, 1, 3, tzinfo=fixed_off)])
        result = to_datetime(dates).to_pydatetime()
        self.assert_(np.array_equal(dates, result))
        result = to_datetime(dates)._mpl_repr()
        self.assert_(np.array_equal(dates, result))

    def test_convert_tz_aware_datetime_datetime(self):
        # #1581

        tz = pytz.timezone('US/Eastern')

        dates = [datetime(2000, 1, 1), datetime(2000, 1, 2),
                 datetime(2000, 1, 3)]

        dates_aware = [tz.localize(x) for x in dates]
        result = to_datetime(dates_aware)
        self.assert_(result.tz.zone == 'US/Eastern')

        converted = to_datetime(dates_aware, utc=True)
        ex_vals = [Timestamp(x).value for x in dates_aware]
        self.assert_(np.array_equal(converted.asi8, ex_vals))
        self.assert_(converted.tz is pytz.utc)

    def test_to_datetime_utc(self):
        from dateutil.parser import parse
        arr = np.array([parse('2012-06-13T01:39:00Z')], dtype=object)

        result = to_datetime(arr, utc=True)
        self.assert_(result.tz is pytz.utc)

    def test_to_datetime_tzlocal(self):
        from dateutil.parser import parse
        from dateutil.tz import tzlocal
        dt = parse('2012-06-13T01:39:00Z')
        dt = dt.replace(tzinfo=tzlocal())

        arr = np.array([dt], dtype=object)

        result = to_datetime(arr, utc=True)
        self.assert_(result.tz is pytz.utc)

        rng = date_range('2012-11-03 03:00', '2012-11-05 03:00', tz=tzlocal())
        arr = rng.to_pydatetime()
        result = to_datetime(arr, utc=True)
        self.assert_(result.tz is pytz.utc)

    def test_frame_no_datetime64_dtype(self):

        dr = date_range('2011/1/1', '2012/1/1', freq='W-FRI')
        dr_tz = dr.tz_localize('US/Eastern')
        e = DataFrame({'A': 'foo', 'B': dr_tz}, index=dr)
        self.assert_(e['B'].dtype == 'M8[ns]')

        # GH 2810 (with timezones)
        datetimes_naive   = [ ts.to_pydatetime() for ts in dr ]
        datetimes_with_tz = [ ts.to_pydatetime() for ts in dr_tz ]
        df = DataFrame({'dr' : dr, 'dr_tz' : dr_tz,
                        'datetimes_naive': datetimes_naive,
                        'datetimes_with_tz' : datetimes_with_tz })
        result = df.get_dtype_counts()
        expected = Series({ 'datetime64[ns]' : 3, 'object' : 1 })
        assert_series_equal(result, expected)

    def test_hongkong_tz_convert(self):
        # #1673
        dr = date_range(
            '2012-01-01', '2012-01-10', freq='D', tz='Hongkong')

        # it works!
        dr.hour

    def test_tz_convert_unsorted(self):
        dr = date_range('2012-03-09', freq='H', periods=100, tz='utc')
        dr = dr.tz_convert('US/Eastern')

        result = dr[::-1].hour
        exp = dr.hour[::-1]
        tm.assert_almost_equal(result, exp)

    def test_shift_localized(self):
        dr = date_range('2011/1/1', '2012/1/1', freq='W-FRI')
        dr_tz = dr.tz_localize('US/Eastern')

        result = dr_tz.shift(1, '10T')
        self.assert_(result.tz == dr_tz.tz)

    def test_tz_aware_asfreq(self):
        dr = date_range(
            '2011-12-01', '2012-07-20', freq='D', tz='US/Eastern')

        s = Series(np.random.randn(len(dr)), index=dr)

        # it works!
        s.asfreq('T')

    def test_static_tzinfo(self):
        # it works!
        index = DatetimeIndex([datetime(2012, 1, 1)], tz='EST')
        index.hour
        index[0]

    def test_tzaware_datetime_to_index(self):
        d = [datetime(2012, 8, 19, tzinfo=pytz.timezone('US/Eastern'))]

        index = DatetimeIndex(d)
        self.assert_(index.tz.zone == 'US/Eastern')

    def test_date_range_span_dst_transition(self):
        # #1778

        # Standard -> Daylight Savings Time
        dr = date_range('03/06/2012 00:00', periods=200, freq='W-FRI',
                        tz='US/Eastern')

        self.assert_((dr.hour == 0).all())

        dr = date_range('2012-11-02', periods=10, tz='US/Eastern')
        self.assert_((dr.hour == 0).all())

    def test_convert_datetime_list(self):
        dr = date_range('2012-06-02', periods=10, tz='US/Eastern')

        dr2 = DatetimeIndex(list(dr), name='foo')
        self.assert_(dr.equals(dr2))
        self.assert_(dr.tz == dr2.tz)
        self.assert_(dr2.name == 'foo')

    def test_frame_from_records_utc(self):
        rec = {'datum': 1.5,
               'begin_time': datetime(2006, 4, 27, tzinfo=pytz.utc)}

        # it works
        DataFrame.from_records([rec], index='begin_time')

    def test_frame_reset_index(self):
        dr = date_range('2012-06-02', periods=10, tz='US/Eastern')
        df = DataFrame(np.random.randn(len(dr)), dr)
        roundtripped = df.reset_index().set_index('index')
        xp = df.index.tz
        rs = roundtripped.index.tz
        self.assertEquals(xp, rs)

    def test_dateutil_tzoffset_support(self):
        from dateutil.tz import tzoffset
        values = [188.5, 328.25]
        tzinfo = tzoffset(None, 7200)
        index = [datetime(2012, 5, 11, 11, tzinfo=tzinfo),
                 datetime(2012, 5, 11, 12, tzinfo=tzinfo)]
        series = Series(data=values, index=index)

        self.assertEquals(series.index.tz, tzinfo)

        # it works! #2443
        repr(series.index[0])

    def test_getitem_pydatetime_tz(self):
        index = date_range(start='2012-12-24 16:00',
                           end='2012-12-24 18:00', freq='H',
                           tz='Europe/Berlin')
        ts = Series(index=index, data=index.hour)
        time_pandas = Timestamp('2012-12-24 17:00', tz='Europe/Berlin')
        time_datetime = datetime(2012, 12, 24, 17, 00,
                                 tzinfo=pytz.timezone('Europe/Berlin'))
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
        
        idx1 = to_datetime(arr).tz_localize('US/Eastern')
        idx2 = DatetimeIndex(start="2005-11-10 08:00:00", freq='H', periods=2, tz='US/Eastern')
        idx3 = DatetimeIndex(arr, tz='US/Eastern')
        idx4 = DatetimeIndex(np.array(arr), tz='US/Eastern')
        
        for other in [idx2, idx3, idx4]:
            self.assert_(idx1.equals(other))


class TestTimeZones(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        _skip_if_no_pytz()

    def test_index_equals_with_tz(self):
        left = date_range('1/1/2011', periods=100, freq='H', tz='utc')
        right = date_range('1/1/2011', periods=100, freq='H',
                           tz='US/Eastern')

        self.assert_(not left.equals(right))

    def test_tz_localize_naive(self):
        rng = date_range('1/1/2011', periods=100, freq='H')

        conv = rng.tz_localize('US/Pacific')
        exp = date_range('1/1/2011', periods=100, freq='H', tz='US/Pacific')

        self.assert_(conv.equals(exp))

    def test_series_frame_tz_localize(self):

        rng = date_range('1/1/2011', periods=100, freq='H')
        ts = Series(1, index=rng)

        result = ts.tz_localize('utc')
        self.assert_(result.index.tz.zone == 'UTC')

        df = DataFrame({'a': 1}, index=rng)
        result = df.tz_localize('utc')
        expected = DataFrame({'a': 1}, rng.tz_localize('UTC'))
        self.assert_(result.index.tz.zone == 'UTC')
        assert_frame_equal(result, expected)

        df = df.T
        result = df.tz_localize('utc', axis=1)
        self.assert_(result.columns.tz.zone == 'UTC')
        assert_frame_equal(result, expected.T)

        # Can't localize if already tz-aware
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')
        ts = Series(1, index=rng)
        assertRaisesRegexp(TypeError, 'Already tz-aware', ts.tz_localize, 'US/Eastern')

    def test_series_frame_tz_convert(self):
        rng = date_range('1/1/2011', periods=200, freq='D',
                         tz='US/Eastern')
        ts = Series(1, index=rng)

        result = ts.tz_convert('Europe/Berlin')
        self.assert_(result.index.tz.zone == 'Europe/Berlin')

        df = DataFrame({'a': 1}, index=rng)
        result = df.tz_convert('Europe/Berlin')
        expected = DataFrame({'a': 1}, rng.tz_convert('Europe/Berlin'))
        self.assert_(result.index.tz.zone == 'Europe/Berlin')
        assert_frame_equal(result, expected)

        df = df.T
        result = df.tz_convert('Europe/Berlin', axis=1)
        self.assert_(result.columns.tz.zone == 'Europe/Berlin')
        assert_frame_equal(result, expected.T)

        # can't convert tz-naive
        rng = date_range('1/1/2011', periods=200, freq='D')
        ts = Series(1, index=rng)
        assertRaisesRegexp(TypeError, "Cannot convert tz-naive", ts.tz_convert, 'US/Eastern')

    def test_join_utc_convert(self):
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')

        left = rng.tz_convert('US/Eastern')
        right = rng.tz_convert('Europe/Berlin')

        for how in ['inner', 'outer', 'left', 'right']:
            result = left.join(left[:-5], how=how)
            self.assert_(isinstance(result, DatetimeIndex))
            self.assert_(result.tz == left.tz)

            result = left.join(right[:-5], how=how)
            self.assert_(isinstance(result, DatetimeIndex))
            self.assert_(result.tz.zone == 'UTC')

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
                          columns=range(3, 6))

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
        self.assert_(ts_result.index.equals(
            ts1.index.asobject.append(ts2.index.asobject)))

        # mixed

        rng1 = date_range('1/1/2011 01:00', periods=1, freq='H')
        rng2 = range(100)
        ts1 = Series(np.random.randn(len(rng1)), index=rng1)
        ts2 = Series(np.random.randn(len(rng2)), index=rng2)
        ts_result = ts1.append(ts2)
        self.assert_(ts_result.index.equals(
            ts1.index.asobject.append(ts2.index)))

    def test_equal_join_ensure_utc(self):
        rng = date_range('1/1/2011', periods=10, freq='H', tz='US/Eastern')
        ts = Series(np.random.randn(len(rng)), index=rng)

        ts_moscow = ts.tz_convert('Europe/Moscow')

        result = ts + ts_moscow
        self.assert_(result.index.tz is pytz.utc)

        result = ts_moscow + ts
        self.assert_(result.index.tz is pytz.utc)

        df = DataFrame({'a': ts})
        df_moscow = df.tz_convert('Europe/Moscow')
        result = df + df_moscow
        self.assert_(result.index.tz is pytz.utc)

        result = df_moscow + df
        self.assert_(result.index.tz is pytz.utc)

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

        self.assert_(result.index.tz == pytz.UTC)
        assert_series_equal(result, expected)

    def test_intersection(self):
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')

        left = rng[10:90][::-1]
        right = rng[20:80][::-1]

        self.assert_(left.tz == rng.tz)
        result = left.intersection(right)
        self.assert_(result.tz == left.tz)

    def test_timestamp_equality_different_timezones(self):
        utc_range = date_range('1/1/2000', periods=20, tz='UTC')

        eastern_range = utc_range.tz_convert('US/Eastern')
        berlin_range = utc_range.tz_convert('Europe/Berlin')

        for a, b, c in zip(utc_range, eastern_range, berlin_range):
            self.assertEquals(a, b)
            self.assertEquals(b, c)
            self.assertEquals(a, c)

        self.assert_((utc_range == eastern_range).all())
        self.assert_((utc_range == berlin_range).all())
        self.assert_((berlin_range == eastern_range).all())

    def test_datetimeindex_tz(self):
        rng = date_range('03/12/2012 00:00', periods=10, freq='W-FRI',
                         tz='US/Eastern')
        rng2 = DatetimeIndex(data=rng, tz='US/Eastern')
        self.assert_(rng.equals(rng2))

    def test_normalize_tz(self):
        rng = date_range('1/1/2000 9:30', periods=10, freq='D',
                         tz='US/Eastern')

        result = rng.normalize()
        expected = date_range('1/1/2000', periods=10, freq='D',
                              tz='US/Eastern')
        self.assert_(result.equals(expected))

        self.assert_(result.is_normalized)
        self.assert_(not rng.is_normalized)

        rng = date_range('1/1/2000 9:30', periods=10, freq='D',
                         tz='UTC')

        result = rng.normalize()
        expected = date_range('1/1/2000', periods=10, freq='D',
                              tz='UTC')
        self.assert_(result.equals(expected))

        self.assert_(result.is_normalized)
        self.assert_(not rng.is_normalized)

        from dateutil.tz import tzlocal
        rng = date_range('1/1/2000 9:30', periods=10, freq='D',
                         tz=tzlocal())
        result = rng.normalize()
        expected = date_range('1/1/2000', periods=10, freq='D',
                              tz=tzlocal())
        self.assert_(result.equals(expected))

        self.assert_(result.is_normalized)
        self.assert_(not rng.is_normalized)

    def test_tzaware_offset(self):
        dates = date_range('2012-11-01', periods=3, tz='US/Pacific')
        offset = dates + offsets.Hour(5)
        self.assertEqual(dates[0] + offsets.Hour(5), offset[0])

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
