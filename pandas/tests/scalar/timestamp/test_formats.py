from datetime import (
    UTC,
    datetime,
    timedelta,
    timezone,
)
import pprint

import dateutil.tz
import numpy as np
import pytest

from pandas.compat import WASM

import pandas as pd

ts_no_ns = pd.Timestamp(
    year=2019,
    month=5,
    day=18,
    hour=15,
    minute=17,
    second=8,
    microsecond=132263,
)
ts_no_ns_year1 = pd.Timestamp(
    year=1,
    month=5,
    day=18,
    hour=15,
    minute=17,
    second=8,
    microsecond=132263,
)
ts_ns = pd.Timestamp(
    year=2019,
    month=5,
    day=18,
    hour=15,
    minute=17,
    second=8,
    microsecond=132263,
    nanosecond=123,
)
ts_ns_tz = pd.Timestamp(
    year=2019,
    month=5,
    day=18,
    hour=15,
    minute=17,
    second=8,
    microsecond=132263,
    nanosecond=123,
    tz="UTC",
)
ts_no_us = pd.Timestamp(
    year=2019,
    month=5,
    day=18,
    hour=15,
    minute=17,
    second=8,
    microsecond=0,
    nanosecond=123,
)


@pytest.mark.parametrize(
    "ts, timespec, expected_iso",
    [
        (ts_no_ns, "auto", "2019-05-18T15:17:08.132263"),
        (ts_no_ns, "seconds", "2019-05-18T15:17:08"),
        (ts_no_ns, "nanoseconds", "2019-05-18T15:17:08.132263000"),
        (ts_no_ns_year1, "seconds", "0001-05-18T15:17:08"),
        (ts_no_ns_year1, "nanoseconds", "0001-05-18T15:17:08.132263000"),
        (ts_ns, "auto", "2019-05-18T15:17:08.132263123"),
        (ts_ns, "hours", "2019-05-18T15"),
        (ts_ns, "minutes", "2019-05-18T15:17"),
        (ts_ns, "seconds", "2019-05-18T15:17:08"),
        (ts_ns, "milliseconds", "2019-05-18T15:17:08.132"),
        (ts_ns, "microseconds", "2019-05-18T15:17:08.132263"),
        (ts_ns, "nanoseconds", "2019-05-18T15:17:08.132263123"),
        (ts_ns_tz, "auto", "2019-05-18T15:17:08.132263123+00:00"),
        (ts_ns_tz, "hours", "2019-05-18T15+00:00"),
        (ts_ns_tz, "minutes", "2019-05-18T15:17+00:00"),
        (ts_ns_tz, "seconds", "2019-05-18T15:17:08+00:00"),
        (ts_ns_tz, "milliseconds", "2019-05-18T15:17:08.132+00:00"),
        (ts_ns_tz, "microseconds", "2019-05-18T15:17:08.132263+00:00"),
        (ts_ns_tz, "nanoseconds", "2019-05-18T15:17:08.132263123+00:00"),
        (ts_no_us, "auto", "2019-05-18T15:17:08.000000123"),
    ],
)
def test_isoformat(ts, timespec, expected_iso):
    assert ts.isoformat(timespec=timespec) == expected_iso


@pytest.mark.parametrize(
    "tz, expected_offset",
    [
        (timezone(timedelta(seconds=33539)), "+09:18:59"),
        (timezone(timedelta(seconds=-28378)), "-07:52:58"),
        (timezone(timedelta(seconds=33539, microseconds=7)), "+09:18:59.000007"),
        (timezone(timedelta(hours=9)), "+09:00"),
        (timezone(timedelta(hours=-8)), "-08:00"),
        ("UTC", "+00:00"),
    ],
)
@pytest.mark.parametrize(
    "microsecond, nanosecond, timespec, frac",
    [
        (0, 123, "auto", ".000000123"),
        (132263, 123, "auto", ".132263123"),
        (0, 0, "nanoseconds", ".000000000"),
        (132263, 0, "nanoseconds", ".132263000"),
    ],
)
def test_isoformat_sub_minute_utc_offset(
    tz, expected_offset, microsecond, nanosecond, timespec, frac
):
    # GH#66547 an offset that is not a whole number of minutes is rendered as
    #  "+HH:MM:SS", so the nanoseconds must not be spliced into it
    ts = pd.Timestamp(
        year=2019,
        month=5,
        day=18,
        hour=15,
        minute=17,
        second=8,
        microsecond=microsecond,
        nanosecond=nanosecond,
        tz=tz,
    )
    assert (
        ts.isoformat(timespec=timespec) == f"2019-05-18T15:17:08{frac}{expected_offset}"
    )


@pytest.mark.parametrize(
    "dt64, expected",
    [
        # years wider than 4 digits, and negative years, reachable at
        #  resolutions coarser than nanosecond
        (np.datetime64(2**63 - 1, "s"), "292277026596-12-04T15:30:07.000000000"),
        (np.datetime64(-(2**63) + 1, "s"), "-292277022657-01-27T08:29:53.000000000"),
        (np.datetime64(2**63 - 1, "ms"), "292278994-08-17T07:12:55.807000000"),
        (
            np.datetime64("10000-01-01T00:00:00.132263", "us"),
            "10000-01-01T00:00:00.132263000",
        ),
        (
            np.datetime64("-1000-01-01T00:00:00.132263", "us"),
            "-1000-01-01T00:00:00.132263000",
        ),
    ],
)
@pytest.mark.parametrize("sep", ["T", "-", "+"])
def test_isoformat_wide_year_no_date_splice(dt64, expected, sep):
    # GH#66547 locating the UTC offset must not mistake a date dash -- or an
    #  unusual separator -- for the offset's sign when the year is not 4 digits
    ts = pd.Timestamp(dt64)
    result = ts.isoformat(sep=sep, timespec="nanoseconds")
    assert result == expected.replace("T", sep)


@pytest.mark.parametrize(
    "tz, expected",
    [
        ("Asia/Tokyo", "1800-01-01T09:18:59.000000001+09:18:59"),
        ("US/Pacific", "1799-12-31T16:07:02.000000001-07:52:58"),
    ],
)
def test_isoformat_pre_standardization_offset(tz, expected):
    # GH#66547 zones carry a seconds component in their pre-standardization
    #  local-mean-time era
    ts = pd.Timestamp("1800-01-01 00:00:00.000000001", tz="UTC").tz_convert(tz)
    result = ts.isoformat()
    assert result == expected
    assert str(ts) == expected.replace("T", " ")
    # the offset must be rendered exactly as the stdlib renders it
    assert ts.to_pydatetime(warn=False).isoformat() == expected.replace(
        ".000000001", ""
    )


class TestTimestampRendering:
    @pytest.mark.parametrize(
        "tz", ["UTC", "Asia/Tokyo", "US/Eastern", "dateutil/America/Los_Angeles"]
    )
    @pytest.mark.parametrize("freq", ["D", "M", "S", "N"])
    @pytest.mark.parametrize(
        "date", ["2014-03-07", "2014-01-01 09:00", "2014-01-01 00:00:00.000000001"]
    )
    @pytest.mark.skipif(WASM, reason="tzset is not available on WASM")
    def test_repr(self, date, freq, tz):
        # avoid to match with timezone name
        freq_repr = f"'{freq}'"
        if tz.startswith("dateutil"):
            tz_repr = tz.replace("dateutil", "")
        else:
            tz_repr = tz

        date_only = pd.Timestamp(date)
        assert date in repr(date_only)
        assert tz_repr not in repr(date_only)
        assert freq_repr not in repr(date_only)
        assert date_only == eval(repr(date_only), dict(vars(pd)))

        date_tz = pd.Timestamp(date, tz=tz)
        assert date in repr(date_tz)
        assert tz_repr in repr(date_tz)
        assert freq_repr not in repr(date_tz)
        assert date_tz == eval(repr(date_tz), dict(vars(pd)))

    def test_repr_utcoffset(self):
        # This can cause the tz field to be populated, but it's redundant to
        # include this information in the date-string.
        date_with_utc_offset = pd.Timestamp("2014-03-13 00:00:00-0400")
        assert "2014-03-13 00:00:00-0400" in repr(date_with_utc_offset)
        assert "tzoffset" not in repr(date_with_utc_offset)
        assert "UTC-04:00" in repr(date_with_utc_offset)
        expr = repr(date_with_utc_offset)
        assert date_with_utc_offset == eval(expr, dict(vars(pd)))

    def test_timestamp_repr_pre1900(self):
        # pre-1900
        stamp = pd.Timestamp("1850-01-01", tz="US/Eastern")
        repr(stamp)

        iso8601 = "1850-01-01 01:23:45.012345"
        stamp = pd.Timestamp(iso8601, tz="US/Eastern")
        result = repr(stamp)
        assert iso8601 in result

    def test_pprint(self):
        # GH#12622
        nested_obj = {"foo": 1, "bar": [{"w": {"a": pd.Timestamp("2011-01-01")}}] * 10}
        result = pprint.pformat(nested_obj, width=50)
        expected = r"""{'bar': [{'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}},
         {'w': {'a': Timestamp('2011-01-01 00:00:00')}}],
 'foo': 1}"""
        assert result == expected

    def test_to_timestamp_repr_is_code(self):
        zs = [
            pd.Timestamp("99-04-17 00:00:00", tz="UTC"),
            pd.Timestamp("2001-04-17 00:00:00", tz="UTC"),
            pd.Timestamp("2001-04-17 00:00:00", tz="America/Los_Angeles"),
            pd.Timestamp("2001-04-17 00:00:00", tz=None),
        ]
        for z in zs:
            assert eval(repr(z), dict(vars(pd))) == z

    def test_repr_matches_pydatetime_no_tz(self):
        dt_date = datetime(2013, 1, 2)
        assert str(dt_date) == str(pd.Timestamp(dt_date))

        dt_datetime = datetime(2013, 1, 2, 12, 1, 3)
        assert str(dt_datetime) == str(pd.Timestamp(dt_datetime))

        dt_datetime_us = datetime(2013, 1, 2, 12, 1, 3, 45)
        assert str(dt_datetime_us) == str(pd.Timestamp(dt_datetime_us))

        ts_nanos_only = pd.Timestamp(200)
        assert str(ts_nanos_only) == "1970-01-01 00:00:00.000000200"

        ts_nanos_micros = pd.Timestamp(1200)
        assert str(ts_nanos_micros) == "1970-01-01 00:00:00.000001200"

    def test_repr_matches_pydatetime_tz_stdlib(self):
        dt_date = datetime(2013, 1, 2, tzinfo=UTC)
        assert str(dt_date) == str(pd.Timestamp(dt_date))

        dt_datetime = datetime(2013, 1, 2, 12, 1, 3, tzinfo=UTC)
        assert str(dt_datetime) == str(pd.Timestamp(dt_datetime))

        dt_datetime_us = datetime(2013, 1, 2, 12, 1, 3, 45, tzinfo=UTC)
        assert str(dt_datetime_us) == str(pd.Timestamp(dt_datetime_us))

    def test_repr_matches_pydatetime_tz_dateutil(self):
        utc = dateutil.tz.tzutc()

        dt_date = datetime(2013, 1, 2, tzinfo=utc)
        assert str(dt_date) == str(pd.Timestamp(dt_date))

        dt_datetime = datetime(2013, 1, 2, 12, 1, 3, tzinfo=utc)
        assert str(dt_datetime) == str(pd.Timestamp(dt_datetime))

        dt_datetime_us = datetime(2013, 1, 2, 12, 1, 3, 45, tzinfo=utc)
        assert str(dt_datetime_us) == str(pd.Timestamp(dt_datetime_us))


class TestTimestampStrftime:
    @pytest.mark.parametrize(
        "ts, fmt, expected",
        [
            # GH#29461 - %N for nanoseconds
            (ts_ns, "%Y-%m-%dT%H:%M:%S.%N", "2019-05-18T15:17:08.132263123"),
            (ts_no_ns, "%N", "132263000"),
            (ts_no_us, "%N", "000000123"),
            # %%N should produce literal %N, not nanoseconds
            (ts_ns, "%%N", "%N"),
            # standard directives still work
            (ts_ns, "%Y-%m-%d %X", "2019-05-18 15:17:08"),
            (ts_ns, "%f", "132263"),
        ],
    )
    def test_strftime(self, ts, fmt, expected):
        assert ts.strftime(fmt) == expected
