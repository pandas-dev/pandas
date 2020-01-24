""" test the scalar Timestamp """

import calendar
from datetime import datetime, timedelta
import locale
import unicodedata

import dateutil
from dateutil.tz import tzutc
import numpy as np
import pytest
import pytz
from pytz import timezone, utc

from pandas._libs.tslibs import conversion
from pandas._libs.tslibs.timezones import dateutil_gettz as gettz, get_timezone
import pandas.compat as compat
from pandas.compat.numpy import np_datetime64_compat
from pandas.errors import OutOfBoundsDatetime
import pandas.util._test_decorators as td

from pandas import NaT, Period, Timedelta, Timestamp
import pandas._testing as tm

from pandas.tseries import offsets


class TestTimestampProperties:
    def test_properties_business(self):
        ts = Timestamp("2017-10-01", freq="B")
        control = Timestamp("2017-10-01")
        assert ts.dayofweek == 6
        assert not ts.is_month_start  # not a weekday
        assert not ts.is_quarter_start  # not a weekday
        # Control case: non-business is month/qtr start
        assert control.is_month_start
        assert control.is_quarter_start

        ts = Timestamp("2017-09-30", freq="B")
        control = Timestamp("2017-09-30")
        assert ts.dayofweek == 5
        assert not ts.is_month_end  # not a weekday
        assert not ts.is_quarter_end  # not a weekday
        # Control case: non-business is month/qtr start
        assert control.is_month_end
        assert control.is_quarter_end

    def test_fields(self):
        def check(value, equal):
            # that we are int like
            assert isinstance(value, int)
            assert value == equal

        # GH 10050
        ts = Timestamp("2015-05-10 09:06:03.000100001")
        check(ts.year, 2015)
        check(ts.month, 5)
        check(ts.day, 10)
        check(ts.hour, 9)
        check(ts.minute, 6)
        check(ts.second, 3)
        msg = "'Timestamp' object has no attribute 'millisecond'"
        with pytest.raises(AttributeError, match=msg):
            ts.millisecond
        check(ts.microsecond, 100)
        check(ts.nanosecond, 1)
        check(ts.dayofweek, 6)
        check(ts.quarter, 2)
        check(ts.dayofyear, 130)
        check(ts.week, 19)
        check(ts.daysinmonth, 31)
        check(ts.daysinmonth, 31)

        # GH 13303
        ts = Timestamp("2014-12-31 23:59:00-05:00", tz="US/Eastern")
        check(ts.year, 2014)
        check(ts.month, 12)
        check(ts.day, 31)
        check(ts.hour, 23)
        check(ts.minute, 59)
        check(ts.second, 0)
        msg = "'Timestamp' object has no attribute 'millisecond'"
        with pytest.raises(AttributeError, match=msg):
            ts.millisecond
        check(ts.microsecond, 0)
        check(ts.nanosecond, 0)
        check(ts.dayofweek, 2)
        check(ts.quarter, 4)
        check(ts.dayofyear, 365)
        check(ts.week, 1)
        check(ts.daysinmonth, 31)

        ts = Timestamp("2014-01-01 00:00:00+01:00")
        starts = ["is_month_start", "is_quarter_start", "is_year_start"]
        for start in starts:
            assert getattr(ts, start)
        ts = Timestamp("2014-12-31 23:59:59+01:00")
        ends = ["is_month_end", "is_year_end", "is_quarter_end"]
        for end in ends:
            assert getattr(ts, end)

    # GH 12806
    @pytest.mark.parametrize(
        "data",
        [Timestamp("2017-08-28 23:00:00"), Timestamp("2017-08-28 23:00:00", tz="EST")],
    )
    @pytest.mark.parametrize(
        "time_locale", [None] if tm.get_locales() is None else [None] + tm.get_locales()
    )
    def test_names(self, data, time_locale):
        # GH 17354
        # Test .day_name(), .month_name
        if time_locale is None:
            expected_day = "Monday"
            expected_month = "August"
        else:
            with tm.set_locale(time_locale, locale.LC_TIME):
                expected_day = calendar.day_name[0].capitalize()
                expected_month = calendar.month_name[8].capitalize()

        result_day = data.day_name(time_locale)
        result_month = data.month_name(time_locale)

        # Work around https://github.com/pandas-dev/pandas/issues/22342
        # different normalizations
        expected_day = unicodedata.normalize("NFD", expected_day)
        expected_month = unicodedata.normalize("NFD", expected_month)

        result_day = unicodedata.normalize("NFD", result_day)
        result_month = unicodedata.normalize("NFD", result_month)

        assert result_day == expected_day
        assert result_month == expected_month

        # Test NaT
        nan_ts = Timestamp(NaT)
        assert np.isnan(nan_ts.day_name(time_locale))
        assert np.isnan(nan_ts.month_name(time_locale))

    def test_is_leap_year(self, tz_naive_fixture):
        tz = tz_naive_fixture
        # GH 13727
        dt = Timestamp("2000-01-01 00:00:00", tz=tz)
        assert dt.is_leap_year
        assert isinstance(dt.is_leap_year, bool)

        dt = Timestamp("1999-01-01 00:00:00", tz=tz)
        assert not dt.is_leap_year

        dt = Timestamp("2004-01-01 00:00:00", tz=tz)
        assert dt.is_leap_year

        dt = Timestamp("2100-01-01 00:00:00", tz=tz)
        assert not dt.is_leap_year

    def test_woy_boundary(self):
        # make sure weeks at year boundaries are correct
        d = datetime(2013, 12, 31)
        result = Timestamp(d).week
        expected = 1  # ISO standard
        assert result == expected

        d = datetime(2008, 12, 28)
        result = Timestamp(d).week
        expected = 52  # ISO standard
        assert result == expected

        d = datetime(2009, 12, 31)
        result = Timestamp(d).week
        expected = 53  # ISO standard
        assert result == expected

        d = datetime(2010, 1, 1)
        result = Timestamp(d).week
        expected = 53  # ISO standard
        assert result == expected

        d = datetime(2010, 1, 3)
        result = Timestamp(d).week
        expected = 53  # ISO standard
        assert result == expected

        result = np.array(
            [
                Timestamp(datetime(*args)).week
                for args in [(2000, 1, 1), (2000, 1, 2), (2005, 1, 1), (2005, 1, 2)]
            ]
        )
        assert (result == [52, 52, 53, 53]).all()

    def test_resolution(self):
        # GH#21336, GH#21365
        dt = Timestamp("2100-01-01 00:00:00")
        assert dt.resolution == Timedelta(nanoseconds=1)

        # Check that the attribute is available on the class, mirroring
        #  the stdlib datetime behavior
        assert Timestamp.resolution == Timedelta(nanoseconds=1)


class TestTimestampConstructors:
    def test_constructor(self):
        base_str = "2014-07-01 09:00"
        base_dt = datetime(2014, 7, 1, 9)
        base_expected = 1_404_205_200_000_000_000

        # confirm base representation is correct
        assert calendar.timegm(base_dt.timetuple()) * 1_000_000_000 == base_expected

        tests = [
            (base_str, base_dt, base_expected),
            (
                "2014-07-01 10:00",
                datetime(2014, 7, 1, 10),
                base_expected + 3600 * 1_000_000_000,
            ),
            (
                "2014-07-01 09:00:00.000008000",
                datetime(2014, 7, 1, 9, 0, 0, 8),
                base_expected + 8000,
            ),
            (
                "2014-07-01 09:00:00.000000005",
                Timestamp("2014-07-01 09:00:00.000000005"),
                base_expected + 5,
            ),
        ]

        timezones = [
            (None, 0),
            ("UTC", 0),
            (pytz.utc, 0),
            ("Asia/Tokyo", 9),
            ("US/Eastern", -4),
            ("dateutil/US/Pacific", -7),
            (pytz.FixedOffset(-180), -3),
            (dateutil.tz.tzoffset(None, 18000), 5),
        ]

        for date_str, date, expected in tests:
            for result in [Timestamp(date_str), Timestamp(date)]:
                # only with timestring
                assert result.value == expected
                assert conversion.pydt_to_i8(result) == expected

                # re-creation shouldn't affect to internal value
                result = Timestamp(result)
                assert result.value == expected
                assert conversion.pydt_to_i8(result) == expected

            # with timezone
            for tz, offset in timezones:
                for result in [Timestamp(date_str, tz=tz), Timestamp(date, tz=tz)]:
                    expected_tz = expected - offset * 3600 * 1_000_000_000
                    assert result.value == expected_tz
                    assert conversion.pydt_to_i8(result) == expected_tz

                    # should preserve tz
                    result = Timestamp(result)
                    assert result.value == expected_tz
                    assert conversion.pydt_to_i8(result) == expected_tz

                    # should convert to UTC
                    if tz is not None:
                        result = Timestamp(result).tz_convert("UTC")
                    else:
                        result = Timestamp(result, tz="UTC")
                    expected_utc = expected - offset * 3600 * 1_000_000_000
                    assert result.value == expected_utc
                    assert conversion.pydt_to_i8(result) == expected_utc

    def test_constructor_with_stringoffset(self):
        # GH 7833
        base_str = "2014-07-01 11:00:00+02:00"
        base_dt = datetime(2014, 7, 1, 9)
        base_expected = 1_404_205_200_000_000_000

        # confirm base representation is correct
        assert calendar.timegm(base_dt.timetuple()) * 1_000_000_000 == base_expected

        tests = [
            (base_str, base_expected),
            ("2014-07-01 12:00:00+02:00", base_expected + 3600 * 1_000_000_000),
            ("2014-07-01 11:00:00.000008000+02:00", base_expected + 8000),
            ("2014-07-01 11:00:00.000000005+02:00", base_expected + 5),
        ]

        timezones = [
            (None, 0),
            ("UTC", 0),
            (pytz.utc, 0),
            ("Asia/Tokyo", 9),
            ("US/Eastern", -4),
            ("dateutil/US/Pacific", -7),
            (pytz.FixedOffset(-180), -3),
            (dateutil.tz.tzoffset(None, 18000), 5),
        ]

        for date_str, expected in tests:
            for result in [Timestamp(date_str)]:
                # only with timestring
                assert result.value == expected
                assert conversion.pydt_to_i8(result) == expected

                # re-creation shouldn't affect to internal value
                result = Timestamp(result)
                assert result.value == expected
                assert conversion.pydt_to_i8(result) == expected

            # with timezone
            for tz, offset in timezones:
                result = Timestamp(date_str, tz=tz)
                expected_tz = expected
                assert result.value == expected_tz
                assert conversion.pydt_to_i8(result) == expected_tz

                # should preserve tz
                result = Timestamp(result)
                assert result.value == expected_tz
                assert conversion.pydt_to_i8(result) == expected_tz

                # should convert to UTC
                result = Timestamp(result).tz_convert("UTC")
                expected_utc = expected
                assert result.value == expected_utc
                assert conversion.pydt_to_i8(result) == expected_utc

        # This should be 2013-11-01 05:00 in UTC
        # converted to Chicago tz
        result = Timestamp("2013-11-01 00:00:00-0500", tz="America/Chicago")
        assert result.value == Timestamp("2013-11-01 05:00").value
        expected = "Timestamp('2013-11-01 00:00:00-0500', tz='America/Chicago')"  # noqa
        assert repr(result) == expected
        assert result == eval(repr(result))

        # This should be 2013-11-01 05:00 in UTC
        # converted to Tokyo tz (+09:00)
        result = Timestamp("2013-11-01 00:00:00-0500", tz="Asia/Tokyo")
        assert result.value == Timestamp("2013-11-01 05:00").value
        expected = "Timestamp('2013-11-01 14:00:00+0900', tz='Asia/Tokyo')"
        assert repr(result) == expected
        assert result == eval(repr(result))

        # GH11708
        # This should be 2015-11-18 10:00 in UTC
        # converted to Asia/Katmandu
        result = Timestamp("2015-11-18 15:45:00+05:45", tz="Asia/Katmandu")
        assert result.value == Timestamp("2015-11-18 10:00").value
        expected = "Timestamp('2015-11-18 15:45:00+0545', tz='Asia/Katmandu')"
        assert repr(result) == expected
        assert result == eval(repr(result))

        # This should be 2015-11-18 10:00 in UTC
        # converted to Asia/Kolkata
        result = Timestamp("2015-11-18 15:30:00+05:30", tz="Asia/Kolkata")
        assert result.value == Timestamp("2015-11-18 10:00").value
        expected = "Timestamp('2015-11-18 15:30:00+0530', tz='Asia/Kolkata')"
        assert repr(result) == expected
        assert result == eval(repr(result))

    def test_constructor_invalid(self):
        with pytest.raises(TypeError, match="Cannot convert input"):
            Timestamp(slice(2))
        with pytest.raises(ValueError, match="Cannot convert Period"):
            Timestamp(Period("1000-01-01"))

    def test_constructor_invalid_tz(self):
        # GH#17690
        with pytest.raises(TypeError, match="must be a datetime.tzinfo"):
            Timestamp("2017-10-22", tzinfo="US/Eastern")

        with pytest.raises(ValueError, match="at most one of"):
            Timestamp("2017-10-22", tzinfo=utc, tz="UTC")

        with pytest.raises(ValueError, match="Invalid frequency:"):
            # GH#5168
            # case where user tries to pass tz as an arg, not kwarg, gets
            # interpreted as a `freq`
            Timestamp("2012-01-01", "US/Pacific")

    def test_constructor_strptime(self):
        # GH25016
        # Test support for Timestamp.strptime
        fmt = "%Y%m%d-%H%M%S-%f%z"
        ts = "20190129-235348-000001+0000"
        with pytest.raises(NotImplementedError):
            Timestamp.strptime(ts, fmt)

    def test_constructor_tz_or_tzinfo(self):
        # GH#17943, GH#17690, GH#5168
        stamps = [
            Timestamp(year=2017, month=10, day=22, tz="UTC"),
            Timestamp(year=2017, month=10, day=22, tzinfo=utc),
            Timestamp(year=2017, month=10, day=22, tz=utc),
            Timestamp(datetime(2017, 10, 22), tzinfo=utc),
            Timestamp(datetime(2017, 10, 22), tz="UTC"),
            Timestamp(datetime(2017, 10, 22), tz=utc),
        ]
        assert all(ts == stamps[0] for ts in stamps)

    def test_constructor_positional(self):
        # see gh-10758
        with pytest.raises(TypeError):
            Timestamp(2000, 1)
        with pytest.raises(ValueError):
            Timestamp(2000, 0, 1)
        with pytest.raises(ValueError):
            Timestamp(2000, 13, 1)
        with pytest.raises(ValueError):
            Timestamp(2000, 1, 0)
        with pytest.raises(ValueError):
            Timestamp(2000, 1, 32)

        # see gh-11630
        assert repr(Timestamp(2015, 11, 12)) == repr(Timestamp("20151112"))
        assert repr(Timestamp(2015, 11, 12, 1, 2, 3, 999999)) == repr(
            Timestamp("2015-11-12 01:02:03.999999")
        )

    def test_constructor_keyword(self):
        # GH 10758
        with pytest.raises(TypeError):
            Timestamp(year=2000, month=1)
        with pytest.raises(ValueError):
            Timestamp(year=2000, month=0, day=1)
        with pytest.raises(ValueError):
            Timestamp(year=2000, month=13, day=1)
        with pytest.raises(ValueError):
            Timestamp(year=2000, month=1, day=0)
        with pytest.raises(ValueError):
            Timestamp(year=2000, month=1, day=32)

        assert repr(Timestamp(year=2015, month=11, day=12)) == repr(
            Timestamp("20151112")
        )

        assert repr(
            Timestamp(
                year=2015,
                month=11,
                day=12,
                hour=1,
                minute=2,
                second=3,
                microsecond=999999,
            )
        ) == repr(Timestamp("2015-11-12 01:02:03.999999"))

    def test_constructor_fromordinal(self):
        base = datetime(2000, 1, 1)

        ts = Timestamp.fromordinal(base.toordinal(), freq="D")
        assert base == ts
        assert ts.freq == "D"
        assert base.toordinal() == ts.toordinal()

        ts = Timestamp.fromordinal(base.toordinal(), tz="US/Eastern")
        assert Timestamp("2000-01-01", tz="US/Eastern") == ts
        assert base.toordinal() == ts.toordinal()

        # GH#3042
        dt = datetime(2011, 4, 16, 0, 0)
        ts = Timestamp.fromordinal(dt.toordinal())
        assert ts.to_pydatetime() == dt

        # with a tzinfo
        stamp = Timestamp("2011-4-16", tz="US/Eastern")
        dt_tz = stamp.to_pydatetime()
        ts = Timestamp.fromordinal(dt_tz.toordinal(), tz="US/Eastern")
        assert ts.to_pydatetime() == dt_tz

    @pytest.mark.parametrize(
        "result",
        [
            Timestamp(datetime(2000, 1, 2, 3, 4, 5, 6), nanosecond=1),
            Timestamp(
                year=2000,
                month=1,
                day=2,
                hour=3,
                minute=4,
                second=5,
                microsecond=6,
                nanosecond=1,
            ),
            Timestamp(
                year=2000,
                month=1,
                day=2,
                hour=3,
                minute=4,
                second=5,
                microsecond=6,
                nanosecond=1,
                tz="UTC",
            ),
            Timestamp(2000, 1, 2, 3, 4, 5, 6, 1, None),
            Timestamp(2000, 1, 2, 3, 4, 5, 6, 1, pytz.UTC),
        ],
    )
    def test_constructor_nanosecond(self, result):
        # GH 18898
        expected = Timestamp(datetime(2000, 1, 2, 3, 4, 5, 6), tz=result.tz)
        expected = expected + Timedelta(nanoseconds=1)
        assert result == expected

    @pytest.mark.parametrize("z", ["Z0", "Z00"])
    def test_constructor_invalid_Z0_isostring(self, z):
        # GH 8910
        with pytest.raises(ValueError):
            Timestamp("2014-11-02 01:00{}".format(z))

    @pytest.mark.parametrize(
        "arg",
        [
            "year",
            "month",
            "day",
            "hour",
            "minute",
            "second",
            "microsecond",
            "nanosecond",
        ],
    )
    def test_invalid_date_kwarg_with_string_input(self, arg):
        kwarg = {arg: 1}
        with pytest.raises(ValueError):
            Timestamp("2010-10-10 12:59:59.999999999", **kwarg)

    def test_out_of_bounds_integer_value(self):
        # GH#26651 check that we raise OutOfBoundsDatetime, not OverflowError
        with pytest.raises(OutOfBoundsDatetime):
            Timestamp(Timestamp.max.value * 2)
        with pytest.raises(OutOfBoundsDatetime):
            Timestamp(Timestamp.min.value * 2)

    def test_out_of_bounds_value(self):
        one_us = np.timedelta64(1).astype("timedelta64[us]")

        # By definition we can't go out of bounds in [ns], so we
        # convert the datetime64s to [us] so we can go out of bounds
        min_ts_us = np.datetime64(Timestamp.min).astype("M8[us]")
        max_ts_us = np.datetime64(Timestamp.max).astype("M8[us]")

        # No error for the min/max datetimes
        Timestamp(min_ts_us)
        Timestamp(max_ts_us)

        # One us less than the minimum is an error
        with pytest.raises(ValueError):
            Timestamp(min_ts_us - one_us)

        # One us more than the maximum is an error
        with pytest.raises(ValueError):
            Timestamp(max_ts_us + one_us)

    def test_out_of_bounds_string(self):
        with pytest.raises(ValueError):
            Timestamp("1676-01-01")
        with pytest.raises(ValueError):
            Timestamp("2263-01-01")

    def test_barely_out_of_bounds(self):
        # GH#19529
        # GH#19382 close enough to bounds that dropping nanos would result
        # in an in-bounds datetime
        with pytest.raises(OutOfBoundsDatetime):
            Timestamp("2262-04-11 23:47:16.854775808")

    def test_bounds_with_different_units(self):
        out_of_bounds_dates = ("1677-09-21", "2262-04-12")

        time_units = ("D", "h", "m", "s", "ms", "us")

        for date_string in out_of_bounds_dates:
            for unit in time_units:
                dt64 = np.datetime64(date_string, dtype="M8[{unit}]".format(unit=unit))
                with pytest.raises(ValueError):
                    Timestamp(dt64)

        in_bounds_dates = ("1677-09-23", "2262-04-11")

        for date_string in in_bounds_dates:
            for unit in time_units:
                dt64 = np.datetime64(date_string, dtype="M8[{unit}]".format(unit=unit))
                Timestamp(dt64)

    def test_min_valid(self):
        # Ensure that Timestamp.min is a valid Timestamp
        Timestamp(Timestamp.min)

    def test_max_valid(self):
        # Ensure that Timestamp.max is a valid Timestamp
        Timestamp(Timestamp.max)

    def test_now(self):
        # GH#9000
        ts_from_string = Timestamp("now")
        ts_from_method = Timestamp.now()
        ts_datetime = datetime.now()

        ts_from_string_tz = Timestamp("now", tz="US/Eastern")
        ts_from_method_tz = Timestamp.now(tz="US/Eastern")

        # Check that the delta between the times is less than 1s (arbitrarily
        # small)
        delta = Timedelta(seconds=1)
        assert abs(ts_from_method - ts_from_string) < delta
        assert abs(ts_datetime - ts_from_method) < delta
        assert abs(ts_from_method_tz - ts_from_string_tz) < delta
        assert (
            abs(
                ts_from_string_tz.tz_localize(None)
                - ts_from_method_tz.tz_localize(None)
            )
            < delta
        )

    def test_today(self):
        ts_from_string = Timestamp("today")
        ts_from_method = Timestamp.today()
        ts_datetime = datetime.today()

        ts_from_string_tz = Timestamp("today", tz="US/Eastern")
        ts_from_method_tz = Timestamp.today(tz="US/Eastern")

        # Check that the delta between the times is less than 1s (arbitrarily
        # small)
        delta = Timedelta(seconds=1)
        assert abs(ts_from_method - ts_from_string) < delta
        assert abs(ts_datetime - ts_from_method) < delta
        assert abs(ts_from_method_tz - ts_from_string_tz) < delta
        assert (
            abs(
                ts_from_string_tz.tz_localize(None)
                - ts_from_method_tz.tz_localize(None)
            )
            < delta
        )

    @pytest.mark.parametrize("tz", [None, pytz.timezone("US/Pacific")])
    def test_disallow_setting_tz(self, tz):
        # GH 3746
        ts = Timestamp("2010")
        with pytest.raises(AttributeError):
            ts.tz = tz

    @pytest.mark.parametrize("offset", ["+0300", "+0200"])
    def test_construct_timestamp_near_dst(self, offset):
        # GH 20854
        expected = Timestamp(
            "2016-10-30 03:00:00{}".format(offset), tz="Europe/Helsinki"
        )
        result = Timestamp(expected).tz_convert("Europe/Helsinki")
        assert result == expected

    @pytest.mark.parametrize(
        "arg", ["2013/01/01 00:00:00+09:00", "2013-01-01 00:00:00+09:00"]
    )
    def test_construct_with_different_string_format(self, arg):
        # GH 12064
        result = Timestamp(arg)
        expected = Timestamp(datetime(2013, 1, 1), tz=pytz.FixedOffset(540))
        assert result == expected

    def test_construct_timestamp_preserve_original_frequency(self):
        # GH 22311
        result = Timestamp(Timestamp("2010-08-08", freq="D")).freq
        expected = offsets.Day()
        assert result == expected

    def test_constructor_invalid_frequency(self):
        # GH 22311
        with pytest.raises(ValueError, match="Invalid frequency:"):
            Timestamp("2012-01-01", freq=[])

    @pytest.mark.parametrize("box", [datetime, Timestamp])
    def test_raise_tz_and_tzinfo_in_datetime_input(self, box):
        # GH 23579
        kwargs = {"year": 2018, "month": 1, "day": 1, "tzinfo": utc}
        with pytest.raises(ValueError, match="Cannot pass a datetime or Timestamp"):
            Timestamp(box(**kwargs), tz="US/Pacific")
        with pytest.raises(ValueError, match="Cannot pass a datetime or Timestamp"):
            Timestamp(box(**kwargs), tzinfo=pytz.timezone("US/Pacific"))

    def test_dont_convert_dateutil_utc_to_pytz_utc(self):
        result = Timestamp(datetime(2018, 1, 1), tz=tzutc())
        expected = Timestamp(datetime(2018, 1, 1)).tz_localize(tzutc())
        assert result == expected

    def test_constructor_subclassed_datetime(self):
        # GH 25851
        # ensure that subclassed datetime works for
        # Timestamp creation
        class SubDatetime(datetime):
            pass

        data = SubDatetime(2000, 1, 1)
        result = Timestamp(data)
        expected = Timestamp(2000, 1, 1)
        assert result == expected

    @pytest.mark.skipif(
        not compat.PY38,
        reason="datetime.fromisocalendar was added in Python version 3.8",
    )
    def test_constructor_fromisocalendar(self):
        # GH 30395
        expected_timestamp = Timestamp("2000-01-03 00:00:00")
        expected_stdlib = datetime.fromisocalendar(2000, 1, 1)
        result = Timestamp.fromisocalendar(2000, 1, 1)
        assert result == expected_timestamp
        assert result == expected_stdlib
        assert isinstance(result, Timestamp)


class TestTimestamp:
    def test_tz(self):
        tstr = "2014-02-01 09:00"
        ts = Timestamp(tstr)
        local = ts.tz_localize("Asia/Tokyo")
        assert local.hour == 9
        assert local == Timestamp(tstr, tz="Asia/Tokyo")
        conv = local.tz_convert("US/Eastern")
        assert conv == Timestamp("2014-01-31 19:00", tz="US/Eastern")
        assert conv.hour == 19

        # preserves nanosecond
        ts = Timestamp(tstr) + offsets.Nano(5)
        local = ts.tz_localize("Asia/Tokyo")
        assert local.hour == 9
        assert local.nanosecond == 5
        conv = local.tz_convert("US/Eastern")
        assert conv.nanosecond == 5
        assert conv.hour == 19

    def test_utc_z_designator(self):
        assert get_timezone(Timestamp("2014-11-02 01:00Z").tzinfo) is utc

    def test_asm8(self):
        np.random.seed(7_960_929)
        ns = [Timestamp.min.value, Timestamp.max.value, 1000]

        for n in ns:
            assert (
                Timestamp(n).asm8.view("i8") == np.datetime64(n, "ns").view("i8") == n
            )

        assert Timestamp("nat").asm8.view("i8") == np.datetime64("nat", "ns").view("i8")

    def test_class_ops_pytz(self):
        def compare(x, y):
            assert int((Timestamp(x).value - Timestamp(y).value) / 1e9) == 0

        compare(Timestamp.now(), datetime.now())
        compare(Timestamp.now("UTC"), datetime.now(timezone("UTC")))
        compare(Timestamp.utcnow(), datetime.utcnow())
        compare(Timestamp.today(), datetime.today())
        current_time = calendar.timegm(datetime.now().utctimetuple())
        compare(
            Timestamp.utcfromtimestamp(current_time),
            datetime.utcfromtimestamp(current_time),
        )
        compare(
            Timestamp.fromtimestamp(current_time), datetime.fromtimestamp(current_time)
        )

        date_component = datetime.utcnow()
        time_component = (date_component + timedelta(minutes=10)).time()
        compare(
            Timestamp.combine(date_component, time_component),
            datetime.combine(date_component, time_component),
        )

    def test_class_ops_dateutil(self):
        def compare(x, y):
            assert (
                int(
                    np.round(Timestamp(x).value / 1e9)
                    - np.round(Timestamp(y).value / 1e9)
                )
                == 0
            )

        compare(Timestamp.now(), datetime.now())
        compare(Timestamp.now("UTC"), datetime.now(tzutc()))
        compare(Timestamp.utcnow(), datetime.utcnow())
        compare(Timestamp.today(), datetime.today())
        current_time = calendar.timegm(datetime.now().utctimetuple())
        compare(
            Timestamp.utcfromtimestamp(current_time),
            datetime.utcfromtimestamp(current_time),
        )
        compare(
            Timestamp.fromtimestamp(current_time), datetime.fromtimestamp(current_time)
        )

        date_component = datetime.utcnow()
        time_component = (date_component + timedelta(minutes=10)).time()
        compare(
            Timestamp.combine(date_component, time_component),
            datetime.combine(date_component, time_component),
        )

    def test_basics_nanos(self):
        val = np.int64(946_684_800_000_000_000).view("M8[ns]")
        stamp = Timestamp(val.view("i8") + 500)
        assert stamp.year == 2000
        assert stamp.month == 1
        assert stamp.microsecond == 0
        assert stamp.nanosecond == 500

        # GH 14415
        val = np.iinfo(np.int64).min + 80_000_000_000_000
        stamp = Timestamp(val)
        assert stamp.year == 1677
        assert stamp.month == 9
        assert stamp.day == 21
        assert stamp.microsecond == 145224
        assert stamp.nanosecond == 192

    @pytest.mark.parametrize(
        "value, check_kwargs",
        [
            [946688461000000000, {}],
            [946688461000000000 / 1000, dict(unit="us")],
            [946688461000000000 / 1_000_000, dict(unit="ms")],
            [946688461000000000 / 1_000_000_000, dict(unit="s")],
            [10957, dict(unit="D", h=0)],
            [
                (946688461000000000 + 500000) / 1000000000,
                dict(unit="s", us=499, ns=964),
            ],
            [(946688461000000000 + 500000000) / 1000000000, dict(unit="s", us=500000)],
            [(946688461000000000 + 500000) / 1000000, dict(unit="ms", us=500)],
            [(946688461000000000 + 500000) / 1000, dict(unit="us", us=500)],
            [(946688461000000000 + 500000000) / 1000000, dict(unit="ms", us=500000)],
            [946688461000000000 / 1000.0 + 5, dict(unit="us", us=5)],
            [946688461000000000 / 1000.0 + 5000, dict(unit="us", us=5000)],
            [946688461000000000 / 1000000.0 + 0.5, dict(unit="ms", us=500)],
            [946688461000000000 / 1000000.0 + 0.005, dict(unit="ms", us=5, ns=5)],
            [946688461000000000 / 1000000000.0 + 0.5, dict(unit="s", us=500000)],
            [10957 + 0.5, dict(unit="D", h=12)],
        ],
    )
    def test_unit(self, value, check_kwargs):
        def check(value, unit=None, h=1, s=1, us=0, ns=0):
            stamp = Timestamp(value, unit=unit)
            assert stamp.year == 2000
            assert stamp.month == 1
            assert stamp.day == 1
            assert stamp.hour == h
            if unit != "D":
                assert stamp.minute == 1
                assert stamp.second == s
                assert stamp.microsecond == us
            else:
                assert stamp.minute == 0
                assert stamp.second == 0
                assert stamp.microsecond == 0
            assert stamp.nanosecond == ns

        check(value, **check_kwargs)

    def test_roundtrip(self):

        # test value to string and back conversions
        # further test accessors
        base = Timestamp("20140101 00:00:00")

        result = Timestamp(base.value + Timedelta("5ms").value)
        assert result == Timestamp(f"{base}.005000")
        assert result.microsecond == 5000

        result = Timestamp(base.value + Timedelta("5us").value)
        assert result == Timestamp(f"{base}.000005")
        assert result.microsecond == 5

        result = Timestamp(base.value + Timedelta("5ns").value)
        assert result == Timestamp(f"{base}.000000005")
        assert result.nanosecond == 5
        assert result.microsecond == 0

        result = Timestamp(base.value + Timedelta("6ms 5us").value)
        assert result == Timestamp(f"{base}.006005")
        assert result.microsecond == 5 + 6 * 1000

        result = Timestamp(base.value + Timedelta("200ms 5us").value)
        assert result == Timestamp(f"{base}.200005")
        assert result.microsecond == 5 + 200 * 1000

    def test_hash_equivalent(self):
        d = {datetime(2011, 1, 1): 5}
        stamp = Timestamp(datetime(2011, 1, 1))
        assert d[stamp] == 5

    def test_tz_conversion_freq(self, tz_naive_fixture):
        # GH25241
        t1 = Timestamp("2019-01-01 10:00", freq="H")
        assert t1.tz_localize(tz=tz_naive_fixture).freq == t1.freq
        t2 = Timestamp("2019-01-02 12:00", tz="UTC", freq="T")
        assert t2.tz_convert(tz="UTC").freq == t2.freq


class TestTimestampNsOperations:
    def test_nanosecond_string_parsing(self):
        ts = Timestamp("2013-05-01 07:15:45.123456789")
        # GH 7878
        expected_repr = "2013-05-01 07:15:45.123456789"
        expected_value = 1_367_392_545_123_456_789
        assert ts.value == expected_value
        assert expected_repr in repr(ts)

        ts = Timestamp("2013-05-01 07:15:45.123456789+09:00", tz="Asia/Tokyo")
        assert ts.value == expected_value - 9 * 3600 * 1_000_000_000
        assert expected_repr in repr(ts)

        ts = Timestamp("2013-05-01 07:15:45.123456789", tz="UTC")
        assert ts.value == expected_value
        assert expected_repr in repr(ts)

        ts = Timestamp("2013-05-01 07:15:45.123456789", tz="US/Eastern")
        assert ts.value == expected_value + 4 * 3600 * 1_000_000_000
        assert expected_repr in repr(ts)

        # GH 10041
        ts = Timestamp("20130501T071545.123456789")
        assert ts.value == expected_value
        assert expected_repr in repr(ts)

    def test_nanosecond_timestamp(self):
        # GH 7610
        expected = 1_293_840_000_000_000_005
        t = Timestamp("2011-01-01") + offsets.Nano(5)
        assert repr(t) == "Timestamp('2011-01-01 00:00:00.000000005')"
        assert t.value == expected
        assert t.nanosecond == 5

        t = Timestamp(t)
        assert repr(t) == "Timestamp('2011-01-01 00:00:00.000000005')"
        assert t.value == expected
        assert t.nanosecond == 5

        t = Timestamp(np_datetime64_compat("2011-01-01 00:00:00.000000005Z"))
        assert repr(t) == "Timestamp('2011-01-01 00:00:00.000000005')"
        assert t.value == expected
        assert t.nanosecond == 5

        expected = 1_293_840_000_000_000_010
        t = t + offsets.Nano(5)
        assert repr(t) == "Timestamp('2011-01-01 00:00:00.000000010')"
        assert t.value == expected
        assert t.nanosecond == 10

        t = Timestamp(t)
        assert repr(t) == "Timestamp('2011-01-01 00:00:00.000000010')"
        assert t.value == expected
        assert t.nanosecond == 10

        t = Timestamp(np_datetime64_compat("2011-01-01 00:00:00.000000010Z"))
        assert repr(t) == "Timestamp('2011-01-01 00:00:00.000000010')"
        assert t.value == expected
        assert t.nanosecond == 10


class TestTimestampToJulianDate:
    def test_compare_1700(self):
        r = Timestamp("1700-06-23").to_julian_date()
        assert r == 2_342_145.5

    def test_compare_2000(self):
        r = Timestamp("2000-04-12").to_julian_date()
        assert r == 2_451_646.5

    def test_compare_2100(self):
        r = Timestamp("2100-08-12").to_julian_date()
        assert r == 2_488_292.5

    def test_compare_hour01(self):
        r = Timestamp("2000-08-12T01:00:00").to_julian_date()
        assert r == 2_451_768.5416666666666666

    def test_compare_hour13(self):
        r = Timestamp("2000-08-12T13:00:00").to_julian_date()
        assert r == 2_451_769.0416666666666666


class TestTimestampConversion:
    def test_conversion(self):
        # GH#9255
        ts = Timestamp("2000-01-01")

        result = ts.to_pydatetime()
        expected = datetime(2000, 1, 1)
        assert result == expected
        assert type(result) == type(expected)

        result = ts.to_datetime64()
        expected = np.datetime64(ts.value, "ns")
        assert result == expected
        assert type(result) == type(expected)
        assert result.dtype == expected.dtype

    def test_to_pydatetime_nonzero_nano(self):
        ts = Timestamp("2011-01-01 9:00:00.123456789")

        # Warn the user of data loss (nanoseconds).
        with tm.assert_produces_warning(UserWarning, check_stacklevel=False):
            expected = datetime(2011, 1, 1, 9, 0, 0, 123456)
            result = ts.to_pydatetime()
            assert result == expected

    def test_timestamp_to_datetime(self):
        stamp = Timestamp("20090415", tz="US/Eastern", freq="D")
        dtval = stamp.to_pydatetime()
        assert stamp == dtval
        assert stamp.tzinfo == dtval.tzinfo

    def test_timestamp_to_datetime_dateutil(self):
        stamp = Timestamp("20090415", tz="dateutil/US/Eastern", freq="D")
        dtval = stamp.to_pydatetime()
        assert stamp == dtval
        assert stamp.tzinfo == dtval.tzinfo

    def test_timestamp_to_datetime_explicit_pytz(self):
        stamp = Timestamp("20090415", tz=pytz.timezone("US/Eastern"), freq="D")
        dtval = stamp.to_pydatetime()
        assert stamp == dtval
        assert stamp.tzinfo == dtval.tzinfo

    @td.skip_if_windows_python_3
    def test_timestamp_to_datetime_explicit_dateutil(self):
        stamp = Timestamp("20090415", tz=gettz("US/Eastern"), freq="D")
        dtval = stamp.to_pydatetime()
        assert stamp == dtval
        assert stamp.tzinfo == dtval.tzinfo

    def test_to_datetime_bijective(self):
        # Ensure that converting to datetime and back only loses precision
        # by going from nanoseconds to microseconds.
        exp_warning = None if Timestamp.max.nanosecond == 0 else UserWarning
        with tm.assert_produces_warning(exp_warning, check_stacklevel=False):
            assert (
                Timestamp(Timestamp.max.to_pydatetime()).value / 1000
                == Timestamp.max.value / 1000
            )

        exp_warning = None if Timestamp.min.nanosecond == 0 else UserWarning
        with tm.assert_produces_warning(exp_warning, check_stacklevel=False):
            assert (
                Timestamp(Timestamp.min.to_pydatetime()).value / 1000
                == Timestamp.min.value / 1000
            )

    def test_to_period_tz_warning(self):
        # GH#21333 make sure a warning is issued when timezone
        # info is lost
        ts = Timestamp("2009-04-15 16:17:18", tz="US/Eastern")
        with tm.assert_produces_warning(UserWarning):
            # warning that timezone info will be lost
            ts.to_period("D")

    def test_to_numpy_alias(self):
        # GH 24653: alias .to_numpy() for scalars
        ts = Timestamp(datetime.now())
        assert ts.to_datetime64() == ts.to_numpy()


class SubDatetime(datetime):
    pass


@pytest.mark.parametrize(
    "lh,rh",
    [
        (SubDatetime(2000, 1, 1), Timedelta(hours=1)),
        (Timedelta(hours=1), SubDatetime(2000, 1, 1)),
    ],
)
def test_dt_subclass_add_timedelta(lh, rh):
    # GH#25851
    # ensure that subclassed datetime works for
    # Timedelta operations
    result = lh + rh
    expected = SubDatetime(2000, 1, 1, 1)
    assert result == expected


def test_constructor_ambigous_dst():
    # GH 24329
    # Make sure that calling Timestamp constructor
    # on Timestamp created from ambiguous time
    # doesn't change Timestamp.value
    ts = Timestamp(1382835600000000000, tz="dateutil/Europe/London")
    expected = ts.value
    result = Timestamp(ts).value
    assert result == expected
