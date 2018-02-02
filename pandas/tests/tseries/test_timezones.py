# pylint: disable-msg=E1101,W0612
import pytest

import pytz
import numpy as np

from dateutil.parser import parse
from dateutil.tz import tzlocal
from datetime import datetime, timedelta, tzinfo, date

import pandas.util.testing as tm
import pandas.tseries.offsets as offsets
from pandas.compat import lrange
from pandas.core.indexes.datetimes import bdate_range, date_range
from pandas._libs.tslibs import timezones
from pandas import (Index, Series, isna, Timestamp, NaT,
                    DatetimeIndex, to_datetime)
from pandas.util.testing import assert_series_equal


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


class TestToDatetimeTimezones(object):
    def test_fixedtz_topydatetime(self):
        dates = np.array([datetime(2000, 1, 1, tzinfo=fixed_off),
                          datetime(2000, 1, 2, tzinfo=fixed_off),
                          datetime(2000, 1, 3, tzinfo=fixed_off)])
        result = to_datetime(dates).to_pydatetime()
        tm.assert_numpy_array_equal(dates, result)
        result = to_datetime(dates)._mpl_repr()
        tm.assert_numpy_array_equal(dates, result)

    def test_fixed_offset(self):
        dates = [datetime(2000, 1, 1, tzinfo=fixed_off),
                 datetime(2000, 1, 2, tzinfo=fixed_off),
                 datetime(2000, 1, 3, tzinfo=fixed_off)]
        result = to_datetime(dates)
        assert result.tz == fixed_off

    def test_to_datetime_utc(self):
        arr = np.array([parse('2012-06-13T01:39:00Z')], dtype=object)

        result = to_datetime(arr, utc=True)
        assert result.tz is pytz.utc

    def test_to_datetime_tzlocal(self):
        dt = parse('2012-06-13T01:39:00Z')
        dt = dt.replace(tzinfo=tzlocal())

        arr = np.array([dt], dtype=object)

        result = to_datetime(arr, utc=True)
        assert result.tz is pytz.utc

        rng = date_range('2012-11-03 03:00', '2012-11-05 03:00', tz=tzlocal())
        arr = rng.to_pydatetime()
        result = to_datetime(arr, utc=True)
        assert result.tz is pytz.utc


class TestDateRangeTimezones(object):
    def test_hongkong_tz_convert(self):
        # GH#1673
        dr = date_range('2012-01-01', '2012-01-10', freq='D', tz='Hongkong')

        # it works!
        dr.hour

    def test_bdate_range_tz_localize(self):
        dr = bdate_range('1/1/2009', '1/1/2010')
        dr_utc = bdate_range('1/1/2009', '1/1/2010', tz=pytz.utc)
        localized = dr.tz_localize(pytz.utc)
        tm.assert_index_equal(dr_utc, localized)

    def test_create_with_fixed_tz(self):
        off = FixedOffset(420, '+07:00')
        start = datetime(2012, 3, 11, 5, 0, 0, tzinfo=off)
        end = datetime(2012, 6, 11, 5, 0, 0, tzinfo=off)
        rng = date_range(start=start, end=end)
        assert off == rng.tz

        rng2 = date_range(start, periods=len(rng), tz=off)
        tm.assert_index_equal(rng, rng2)

        rng3 = date_range('3/11/2012 05:00:00+07:00',
                          '6/11/2012 05:00:00+07:00')
        assert (rng.values == rng3.values).all()

    def test_create_with_fixedoffset_noname(self):
        off = fixed_off_no_name
        start = datetime(2012, 3, 11, 5, 0, 0, tzinfo=off)
        end = datetime(2012, 6, 11, 5, 0, 0, tzinfo=off)
        rng = date_range(start=start, end=end)
        assert off == rng.tz

        idx = Index([start, end])
        assert off == idx.tz

    def test_date_range_localize(self):
        rng = date_range('3/11/2012 03:00', periods=15, freq='H',
                         tz='US/Eastern')
        rng2 = DatetimeIndex(['3/11/2012 03:00', '3/11/2012 04:00'],
                             tz='US/Eastern')
        rng3 = date_range('3/11/2012 03:00', periods=15, freq='H')
        rng3 = rng3.tz_localize('US/Eastern')

        tm.assert_index_equal(rng, rng3)

        # DST transition time
        val = rng[0]
        exp = Timestamp('3/11/2012 03:00', tz='US/Eastern')

        assert val.hour == 3
        assert exp.hour == 3
        assert val == exp  # same UTC value
        tm.assert_index_equal(rng[:2], rng2)

        # Right before the DST transition
        rng = date_range('3/11/2012 00:00', periods=2, freq='H',
                         tz='US/Eastern')
        rng2 = DatetimeIndex(['3/11/2012 00:00', '3/11/2012 01:00'],
                             tz='US/Eastern')
        tm.assert_index_equal(rng, rng2)
        exp = Timestamp('3/11/2012 00:00', tz='US/Eastern')
        assert exp.hour == 0
        assert rng[0] == exp
        exp = Timestamp('3/11/2012 01:00', tz='US/Eastern')
        assert exp.hour == 1
        assert rng[1] == exp

        rng = date_range('3/11/2012 00:00', periods=10, freq='H',
                         tz='US/Eastern')
        assert rng[2].hour == 3


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

    def test_utc_to_local_no_modify(self):
        rng = date_range('3/11/2012', '3/12/2012', freq='H', tz='utc')
        rng_eastern = rng.tz_convert(self.tzstr('US/Eastern'))

        # Values are unmodified
        tm.assert_numpy_array_equal(rng.asi8, rng_eastern.asi8)

        assert self.cmptz(rng_eastern.tz, self.tz('US/Eastern'))

    def test_utc_to_local_no_modify_explicit(self):
        rng = date_range('3/11/2012', '3/12/2012', freq='H', tz='utc')
        rng_eastern = rng.tz_convert(self.tz('US/Eastern'))

        # Values are unmodified
        tm.assert_numpy_array_equal(rng.asi8, rng_eastern.asi8)

        assert rng_eastern.tz == self.tz('US/Eastern')

    def test_timestamp_constructed_by_date_and_tz(self):
        # Fix Issue 2993, Timestamp cannot be constructed by datetime.date
        # and tz correctly

        result = Timestamp(date(2012, 3, 11), tz=self.tzstr('US/Eastern'))

        expected = Timestamp('3/11/2012', tz=self.tzstr('US/Eastern'))
        assert result.hour == expected.hour
        assert result == expected

    def test_timestamp_constructed_by_date_and_tz_explicit(self):
        # Fix Issue 2993, Timestamp cannot be constructed by datetime.date
        # and tz correctly

        result = Timestamp(date(2012, 3, 11), tz=self.tz('US/Eastern'))

        expected = Timestamp('3/11/2012', tz=self.tz('US/Eastern'))
        assert result.hour == expected.hour
        assert result == expected

    def test_timedelta_push_over_dst_boundary(self):
        # #1389

        # 4 hours before DST transition
        stamp = Timestamp('3/10/2012 22:00', tz=self.tzstr('US/Eastern'))

        result = stamp + timedelta(hours=6)

        # spring forward, + "7" hours
        expected = Timestamp('3/11/2012 05:00', tz=self.tzstr('US/Eastern'))

        assert result == expected

    def test_timedelta_push_over_dst_boundary_explicit(self):
        # #1389

        # 4 hours before DST transition
        stamp = Timestamp('3/10/2012 22:00', tz=self.tz('US/Eastern'))

        result = stamp + timedelta(hours=6)

        # spring forward, + "7" hours
        expected = Timestamp('3/11/2012 05:00', tz=self.tz('US/Eastern'))

        assert result == expected

    def test_dti_tz_localize(self):
        dti = DatetimeIndex(start='1/1/2005', end='1/1/2005 0:00:30.256',
                            freq='L')
        dti2 = dti.tz_localize(self.tzstr('US/Eastern'))

        dti_utc = DatetimeIndex(start='1/1/2005 05:00',
                                end='1/1/2005 5:00:30.256', freq='L', tz='utc')

        tm.assert_numpy_array_equal(dti2.values, dti_utc.values)

        dti3 = dti2.tz_convert(self.tzstr('US/Pacific'))
        tm.assert_numpy_array_equal(dti3.values, dti_utc.values)

    def test_tz_localize_empty_series(self):
        # #2248

        ts = Series()

        ts2 = ts.tz_localize('utc')
        assert ts2.index.tz == pytz.utc

        ts2 = ts.tz_localize(self.tzstr('US/Eastern'))
        assert self.cmptz(ts2.index.tz, self.tz('US/Eastern'))

    def test_create_with_tz(self):
        stamp = Timestamp('3/11/2012 05:00', tz=self.tzstr('US/Eastern'))
        assert stamp.hour == 5

        rng = date_range('3/11/2012 04:00', periods=10, freq='H',
                         tz=self.tzstr('US/Eastern'))

        assert stamp == rng[1]

        utc_stamp = Timestamp('3/11/2012 05:00', tz='utc')
        assert utc_stamp.tzinfo is pytz.utc
        assert utc_stamp.hour == 5

        utc_stamp = Timestamp('3/11/2012 05:00').tz_localize('utc')
        assert utc_stamp.hour == 5

    def test_utc_box_timestamp_and_localize(self):
        rng = date_range('3/11/2012', '3/12/2012', freq='H', tz='utc')
        rng_eastern = rng.tz_convert(self.tzstr('US/Eastern'))

        tz = self.tz('US/Eastern')
        expected = rng[-1].astimezone(tz)

        stamp = rng_eastern[-1]
        assert stamp == expected
        assert stamp.tzinfo == expected.tzinfo

        # right tzinfo
        rng = date_range('3/13/2012', '3/14/2012', freq='H', tz='utc')
        rng_eastern = rng.tz_convert(self.tzstr('US/Eastern'))
        # test not valid for dateutil timezones.
        # assert 'EDT' in repr(rng_eastern[0].tzinfo)
        assert ('EDT' in repr(rng_eastern[0].tzinfo) or
                'tzfile' in repr(rng_eastern[0].tzinfo))

    def test_timestamp_tz_convert(self):
        strdates = ['1/1/2012', '3/1/2012', '4/1/2012']
        idx = DatetimeIndex(strdates, tz=self.tzstr('US/Eastern'))

        conv = idx[0].tz_convert(self.tzstr('US/Pacific'))
        expected = idx.tz_convert(self.tzstr('US/Pacific'))[0]

        assert conv == expected

    def test_pass_dates_localize_to_utc(self):
        strdates = ['1/1/2012', '3/1/2012', '4/1/2012']

        idx = DatetimeIndex(strdates)
        conv = idx.tz_localize(self.tzstr('US/Eastern'))

        fromdates = DatetimeIndex(strdates, tz=self.tzstr('US/Eastern'))

        assert conv.tz == fromdates.tz
        tm.assert_numpy_array_equal(conv.values, fromdates.values)

    def test_field_access_localize(self):
        strdates = ['1/1/2012', '3/1/2012', '4/1/2012']
        rng = DatetimeIndex(strdates, tz=self.tzstr('US/Eastern'))
        assert (rng.hour == 0).all()

        # a more unusual time zone, #1946
        dr = date_range('2011-10-02 00:00', freq='h', periods=10,
                        tz=self.tzstr('America/Atikokan'))

        expected = Index(np.arange(10, dtype=np.int64))
        tm.assert_index_equal(dr.hour, expected)

    def test_with_tz(self):
        tz = self.tz('US/Central')

        # just want it to work
        start = datetime(2011, 3, 12, tzinfo=pytz.utc)
        dr = bdate_range(start, periods=50, freq=offsets.Hour())
        assert dr.tz is pytz.utc

        # DateRange with naive datetimes
        dr = bdate_range('1/1/2005', '1/1/2009', tz=pytz.utc)
        dr = bdate_range('1/1/2005', '1/1/2009', tz=tz)

        # normalized
        central = dr.tz_convert(tz)
        assert central.tz is tz
        comp = self.localize(tz, central[0].to_pydatetime().replace(
            tzinfo=None)).tzinfo
        assert central[0].tz is comp

        # compare vs a localized tz
        comp = self.localize(tz,
                             dr[0].to_pydatetime().replace(tzinfo=None)).tzinfo
        assert central[0].tz is comp

        # datetimes with tzinfo set
        dr = bdate_range(datetime(2005, 1, 1, tzinfo=pytz.utc),
                         datetime(2009, 1, 1, tzinfo=pytz.utc))

        pytest.raises(Exception, bdate_range,
                      datetime(2005, 1, 1, tzinfo=pytz.utc), '1/1/2009',
                      tz=tz)

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

    def test_tz_string(self):
        result = date_range('1/1/2000', periods=10,
                            tz=self.tzstr('US/Eastern'))
        expected = date_range('1/1/2000', periods=10, tz=self.tz('US/Eastern'))

        tm.assert_index_equal(result, expected)

    def test_dti_index_with_timezone_repr(self):
        rng = date_range('4/13/2010', '5/6/2010')

        rng_eastern = rng.tz_localize(self.tzstr('US/Eastern'))

        rng_repr = repr(rng_eastern)
        assert '2010-04-13 00:00:00' in rng_repr

    def test_dti_astype_asobject_tzinfos(self):
        # #1345

        # dates around a dst transition
        rng = date_range('2/13/2010', '5/6/2010', tz=self.tzstr('US/Eastern'))

        objs = rng.astype(object)
        for i, x in enumerate(objs):
            exval = rng[i]
            assert x == exval
            assert x.tzinfo == exval.tzinfo

        objs = rng.astype(object)
        for i, x in enumerate(objs):
            exval = rng[i]
            assert x == exval
            assert x.tzinfo == exval.tzinfo

    def test_dti_localized_at_time_between_time(self):
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

    def test_convert_tz_aware_datetime_datetime(self):
        # #1581

        tz = self.tz('US/Eastern')

        dates = [datetime(2000, 1, 1), datetime(2000, 1, 2),
                 datetime(2000, 1, 3)]

        dates_aware = [self.localize(tz, x) for x in dates]
        result = to_datetime(dates_aware)
        assert self.cmptz(result.tz, self.tz('US/Eastern'))

        converted = to_datetime(dates_aware, utc=True)
        ex_vals = np.array([Timestamp(x).value for x in dates_aware])
        tm.assert_numpy_array_equal(converted.asi8, ex_vals)
        assert converted.tz is pytz.utc

    def test_shift_localized(self):
        dr = date_range('2011/1/1', '2012/1/1', freq='W-FRI')
        dr_tz = dr.tz_localize(self.tzstr('US/Eastern'))

        result = dr_tz.shift(1, '10T')
        assert result.tz == dr_tz.tz

    def test_tz_aware_asfreq(self):
        dr = date_range('2011-12-01', '2012-07-20', freq='D',
                        tz=self.tzstr('US/Eastern'))

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
        assert self.cmptz(index.tz, self.tz('US/Eastern'))

    def test_date_range_span_dst_transition(self):
        # #1778

        # Standard -> Daylight Savings Time
        dr = date_range('03/06/2012 00:00', periods=200, freq='W-FRI',
                        tz='US/Eastern')

        assert (dr.hour == 0).all()

        dr = date_range('2012-11-02', periods=10, tz=self.tzstr('US/Eastern'))
        assert (dr.hour == 0).all()

    def test_convert_datetime_list(self):
        dr = date_range('2012-06-02', periods=10,
                        tz=self.tzstr('US/Eastern'), name='foo')
        dr2 = DatetimeIndex(list(dr), name='foo')
        tm.assert_index_equal(dr, dr2)
        assert dr.tz == dr2.tz
        assert dr2.name == 'foo'

    def test_getitem_pydatetime_tz(self):
        index = date_range(start='2012-12-24 16:00', end='2012-12-24 18:00',
                           freq='H', tz=self.tzstr('Europe/Berlin'))
        ts = Series(index=index, data=index.hour)
        time_pandas = Timestamp('2012-12-24 17:00',
                                tz=self.tzstr('Europe/Berlin'))
        time_datetime = self.localize(
            self.tz('Europe/Berlin'), datetime(2012, 12, 24, 17, 0))
        assert ts[time_pandas] == ts[time_datetime]

    def test_datetimeindex_tz_nat(self):
        idx = to_datetime([Timestamp("2013-1-1", tz=self.tzstr('US/Eastern')),
                           NaT])

        assert isna(idx[1])
        assert idx[0].tzinfo is not None

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


class TestTimeZones(object):

    def test_join_utc_convert(self):
        rng = date_range('1/1/2011', periods=100, freq='H', tz='utc')

        left = rng.tz_convert('US/Eastern')
        right = rng.tz_convert('Europe/Berlin')

        for how in ['inner', 'outer', 'left', 'right']:
            result = left.join(left[:-5], how=how)
            assert isinstance(result, DatetimeIndex)
            assert result.tz == left.tz

            result = left.join(right[:-5], how=how)
            assert isinstance(result, DatetimeIndex)
            assert result.tz.zone == 'UTC'

    def test_datetimeindex_tz(self):
        rng = date_range('03/12/2012 00:00', periods=10, freq='W-FRI',
                         tz='US/Eastern')
        rng2 = DatetimeIndex(data=rng, tz='US/Eastern')
        tm.assert_index_equal(rng, rng2)
