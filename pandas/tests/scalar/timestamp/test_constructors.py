import calendar
import copy
from datetime import (
    UTC,
    date,
    datetime,
    timedelta,
    timezone,
)
import re
import zoneinfo

import dateutil.tz
from dateutil.tz import (
    gettz,
    tzoffset,
    tzutc,
)
import numpy as np
import pytest

from pandas._libs.tslibs.dtypes import NpyDatetimeUnit
from pandas.compat import (
    PY313,
    PY314,
    PY315,
)
from pandas.errors import (
    OutOfBoundsDatetime,
    Pandas4Warning,
)

from pandas import (
    NA,
    NaT,
    Period,
    Timedelta,
    Timestamp,
)
import pandas._testing as tm


class TestTimestampConstructorUnitKeyword:
    @pytest.mark.parametrize("typ", [int, float])
    def test_constructor_int_float_with_YM_unit(self, typ):
        # GH#47266 avoid the conversions in cast_from_unit
        val = typ(150)

        ts = Timestamp(val, unit="Y")
        expected = Timestamp("2120-01-01")
        assert ts == expected

        ts = Timestamp(val, unit="M")
        expected = Timestamp("1982-07-01")
        assert ts == expected

    @pytest.mark.parametrize("typ", [int, float])
    def test_construct_from_int_float_with_unit_out_of_bound_raises(self, typ):
        # GH#50870  make sure we get an OutOfBoundsDatetime instead of OverflowError
        val = typ(150000000000000)

        msg = f"cannot convert input {int(val)} with the unit 'D'"
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            Timestamp(val, unit="D")

    def test_constructor_float_not_round_with_YM_unit_raises(self):
        # GH#47267 avoid the conversions in cast_from-unit

        msg = "Conversion of non-round float with unit=[MY] is ambiguous"
        with pytest.raises(ValueError, match=msg):
            Timestamp(150.5, unit="Y")

        with pytest.raises(ValueError, match=msg):
            Timestamp(150.5, unit="M")

    @pytest.mark.parametrize(
        "value, check_kwargs",
        [
            [946688461000000000, {}],
            [946688461000000000 / 1000, {"unit": "us"}],
            [946688461000000000 / 1_000_000, {"unit": "ms"}],
            [946688461000000000 / 1_000_000_000, {"unit": "s"}],
            [10957, {"unit": "D", "h": 0}],
            [
                (946688461000000000 + 500000) / 1000000000,
                {"unit": "s", "us": 499, "ns": 964},
            ],
            [
                (946688461000000000 + 500000000) / 1000000000,
                {"unit": "s", "us": 500000},
            ],
            [(946688461000000000 + 500000) / 1000000, {"unit": "ms", "us": 500}],
            [(946688461000000000 + 500000) / 1000, {"unit": "us", "us": 500}],
            [(946688461000000000 + 500000000) / 1000000, {"unit": "ms", "us": 500000}],
            [946688461000000000 / 1000.0 + 5, {"unit": "us", "us": 5}],
            [946688461000000000 / 1000.0 + 5000, {"unit": "us", "us": 5000}],
            [946688461000000000 / 1000000.0 + 0.5, {"unit": "ms", "us": 500}],
            [946688461000000000 / 1000000.0 + 0.005, {"unit": "ms", "us": 5, "ns": 5}],
            [946688461000000000 / 1000000000.0 + 0.5, {"unit": "s", "us": 500000}],
            [10957 + 0.5, {"unit": "D", "h": 12}],
        ],
    )
    def test_construct_with_unit(self, value, check_kwargs):
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


class TestTimestampConstructorFoldKeyword:
    def test_timestamp_constructor_invalid_fold_raise(self):
        # Test for GH#25057
        # Valid fold values are only [None, 0, 1]
        msg = "Valid values for the fold argument are None, 0, or 1."
        with pytest.raises(ValueError, match=msg):
            Timestamp(123, fold=2)

    def test_timestamp_constructor_pytz_fold_raise(self):
        # Test for GH#25057
        # pytz doesn't support fold. Check that we raise
        # if fold is passed with pytz
        pytz = pytest.importorskip("pytz")
        msg = "pytz timezones do not support fold. Please use dateutil timezones."
        tz = pytz.timezone("Europe/London")
        with pytest.raises(ValueError, match=msg):
            Timestamp(datetime(2019, 10, 27, 0, 30, 0, 0), tz=tz, fold=0)

    @pytest.mark.parametrize("fold", [0, 1])
    @pytest.mark.parametrize(
        "ts_input",
        [
            1572136200000000000,
            1572136200000000000.0,
            np.datetime64(1572136200000000000, "ns"),
            "2019-10-27 01:30:00+01:00",
            datetime(2019, 10, 27, 0, 30, 0, 0, tzinfo=UTC),
        ],
    )
    def test_timestamp_constructor_fold_conflict(self, ts_input, fold):
        # Test for GH#25057
        # Check that we raise on fold conflict
        msg = (
            "Cannot pass fold with possibly unambiguous input: int, float, "
            "numpy.datetime64, str, or timezone-aware datetime-like. "
            "Pass naive datetime-like or build Timestamp from components."
        )
        with pytest.raises(ValueError, match=msg):
            Timestamp(ts_input=ts_input, fold=fold)

    @pytest.mark.parametrize("tz", ["dateutil/Europe/London", None])
    @pytest.mark.parametrize("fold", [0, 1])
    def test_timestamp_constructor_retain_fold(self, tz, fold):
        # Test for GH#25057
        # Check that we retain fold
        ts = Timestamp(year=2019, month=10, day=27, hour=1, minute=30, tz=tz, fold=fold)
        result = ts.fold
        expected = fold
        assert result == expected

    def test_timestamp_constructor_string_tz_respects_fold(self):
        # GH#55932 fold must be respected when tz is given as a plain string
        #  (which now resolves to zoneinfo, not pytz)
        utc0 = Timestamp("2023-11-05T08:30:00Z")
        utc1 = Timestamp("2023-11-05T09:30:00Z")

        ts0 = Timestamp(
            year=2023, month=11, day=5, hour=1, minute=30, fold=0, tz="US/Pacific"
        )
        ts1 = Timestamp(
            year=2023, month=11, day=5, hour=1, minute=30, fold=1, tz="US/Pacific"
        )
        assert ts0 == utc0
        assert ts1 == utc1

        # Naive datetime + string tz + fold (second case in GH#55932 thread)
        naive = datetime(2022, 11, 6, 1, 6, 58)
        ts0 = Timestamp(naive, tz="America/New_York", fold=0)
        ts1 = Timestamp(naive, tz="America/New_York", fold=1)
        assert ts0 != ts1
        assert (ts1 - ts0) == Timedelta(hours=1)

    @pytest.mark.parametrize(
        "tz",
        [
            "dateutil/Europe/London",
            zoneinfo.ZoneInfo("Europe/London"),
        ],
    )
    @pytest.mark.parametrize(
        "ts_input,fold_out",
        [
            (1572136200000000000, 0),
            (1572139800000000000, 1),
            ("2019-10-27 01:30:00+01:00", 0),
            ("2019-10-27 01:30:00+00:00", 1),
            (datetime(2019, 10, 27, 1, 30, 0, 0, fold=0), 0),
            (datetime(2019, 10, 27, 1, 30, 0, 0, fold=1), 1),
        ],
    )
    def test_timestamp_constructor_infer_fold_from_value(self, tz, ts_input, fold_out):
        # Test for GH#25057
        # Check that we infer fold correctly based on timestamps since utc
        # or strings
        ts = Timestamp(ts_input, tz=tz)
        result = ts.fold
        expected = fold_out
        assert result == expected

    @pytest.mark.parametrize("tz", ["dateutil/Europe/London"])
    @pytest.mark.parametrize(
        "fold,value_out",
        [
            (0, 1572136200000000),
            (1, 1572139800000000),
        ],
    )
    def test_timestamp_constructor_adjust_value_for_fold(self, tz, fold, value_out):
        # Test for GH#25057
        # Check that we adjust value for fold correctly
        # based on timestamps since utc
        ts_input = datetime(2019, 10, 27, 1, 30)
        ts = Timestamp(ts_input, tz=tz, fold=fold)
        result = ts._value
        expected = value_out
        assert result == expected


class TestTimestampConstructorPositionalAndKeywordSupport:
    def test_constructor_positional(self):
        # see GH#10758
        msg = "'NoneType' object cannot be interpreted as an integer"
        with pytest.raises(TypeError, match=msg):
            Timestamp(2000, 1)

        msg = "month must be in 1..12"
        with pytest.raises(ValueError, match=msg):
            Timestamp(2000, 0, 1)
        with pytest.raises(ValueError, match=msg):
            Timestamp(2000, 13, 1)

        if PY314:
            msg = "must be in range 1..31 for month 1 in year 2000"
        else:
            msg = "day is out of range for month"
        with pytest.raises(ValueError, match=msg):
            Timestamp(2000, 1, 0)
        with pytest.raises(ValueError, match=msg):
            Timestamp(2000, 1, 32)

        # see gh-11630
        assert repr(Timestamp(2015, 11, 12)) == repr(Timestamp("20151112"))
        assert repr(Timestamp(2015, 11, 12, 1, 2, 3, 999999)) == repr(
            Timestamp("2015-11-12 01:02:03.999999")
        )

    def test_constructor_keyword(self):
        # GH#10758
        msg = "|".join(
            [
                r"datetime\(\) missing required argument 'day'",  # PY315
                "function missing required argument 'day'",
                "Required argument 'day'",
            ]
        )
        with pytest.raises(TypeError, match=msg):
            Timestamp(year=2000, month=1)

        msg = "month must be in 1..12"
        with pytest.raises(ValueError, match=msg):
            Timestamp(year=2000, month=0, day=1)
        with pytest.raises(ValueError, match=msg):
            Timestamp(year=2000, month=13, day=1)

        if PY314:
            msg = "must be in range 1..31 for month 1 in year 2000"
        else:
            msg = "day is out of range for month"
        with pytest.raises(ValueError, match=msg):
            Timestamp(year=2000, month=1, day=0)
        with pytest.raises(ValueError, match=msg):
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
        msg = "Cannot pass a date attribute keyword argument"
        with pytest.raises(ValueError, match=msg):
            Timestamp("2010-10-10 12:59:59.999999999", **kwarg)

    @pytest.mark.parametrize(
        "value",
        [
            datetime(2020, 12, 31),
            date(2020, 12, 31),
            Timestamp("2020-12-31"),
            np.datetime64("2020-12-31", "ns"),
            1609372800,
            1609372800.0,
        ],
    )
    @pytest.mark.parametrize(
        "arg", ["month", "day", "hour", "minute", "second", "microsecond"]
    )
    def test_invalid_date_kwarg_with_value_input(self, value, arg):
        # GH#31930 these used to be silently dropped
        msg = "Cannot pass both a value to convert and date attributes"
        with pytest.raises(ValueError, match=msg):
            Timestamp(value, **{arg: 1})

    @pytest.mark.parametrize(
        "value",
        [
            datetime(2020, 12, 31),
            date(2020, 12, 31),
            Timestamp("2020-12-31"),
            np.datetime64("2020-12-31", "ns"),
        ],
    )
    @pytest.mark.parametrize("year", [2020, np.int64(2020), np.uint32(2020)])
    def test_invalid_year_kwarg_with_value_input(self, value, year):
        # GH#31930 an integer `year` alongside a datetime-like value is the
        #  by-component form misfiring, not Timestamp(year, month, ...)
        msg = "Cannot pass both a value to convert and date attributes"
        with pytest.raises(ValueError, match=msg):
            Timestamp(value, year=year)
        with pytest.raises(ValueError, match=msg):
            Timestamp(value, year)

    def test_invalid_date_kwarg_names_the_offending_arguments(self):
        # GH#31930 the message lists what was actually passed, as Timedelta's
        #  equivalent check does
        msg = re.escape(
            "Cannot pass both a value to convert and date attributes, got "
            "['hour', 'minute']; 'tz' is keyword-only"
        )
        with pytest.raises(ValueError, match=msg):
            Timestamp(datetime(2020, 12, 31), hour=5, minute=3)

    @pytest.mark.parametrize("value", [datetime(2012, 1, 1), 1609372800])
    def test_tz_passed_positionally(self, value):
        # GH#31930, GH#5168 a tz passed as an arg rather than a kwarg lands on
        #  `year`; this is what the "'tz' is keyword-only" hint is for
        msg = re.escape(
            "Cannot pass both a value to convert and date attributes, got "
            "['year']; 'tz' is keyword-only"
        )
        with pytest.raises(ValueError, match=msg):
            Timestamp(value, "UTC")

    def test_components_without_year_report_value_form(self):
        # GH#31930 `2020` sits in the ts_input slot and is indistinguishable
        #  from an epoch value, so this reports the value-form error
        msg = re.escape(
            "Cannot pass both a value to convert and date attributes, got "
            "['month', 'day']; 'tz' is keyword-only"
        )
        with pytest.raises(ValueError, match=msg):
            Timestamp(2020, month=12, day=31)

    @pytest.mark.parametrize(
        "value", [[1, 2], b"2020-12-31", np.timedelta64(1, "s"), slice(2)]
    )
    def test_invalid_date_kwarg_with_unconvertible_input(self, value):
        # GH#31930 an input we cannot convert keeps complaining about the type
        msg = "Cannot convert input"
        with pytest.raises(TypeError, match=msg):
            Timestamp(value, hour=5)

    def test_float_year_positional_raises(self):
        # GH#31930 this is the by-component form with a bad year, not a value
        #  mixed with components
        msg = "'float' object cannot be interpreted as an integer"
        with pytest.raises(TypeError, match=msg):
            Timestamp(2020.0, 12, 31)

    @pytest.mark.parametrize(
        "value, msg",
        [
            (1609372800, "'NoneType' object cannot be interpreted as an integer"),
            (1609372800.0, "'float' object cannot be interpreted as an integer"),
        ],
    )
    def test_year_with_numeric_value_undetectable(self, value, msg):
        # GH#31930 Timestamp(epoch, year=2020) is read as the by-component
        #  form with month and day left out, so it fails the way
        #  Timestamp(2020, ...) does rather than reporting a mixed call
        with pytest.raises(TypeError, match=msg):
            Timestamp(value, year=2020)

    def test_int_value_with_year_kwarg_reads_as_components(self):
        # GH#31930 an int first argument is indistinguishable from a year, so
        #  this documented case is silently read as the by-component form
        assert Timestamp(2020, year=1, month=1) == Timestamp(2020, 1, 1)

    @pytest.mark.parametrize("value", [datetime(2020, 12, 31), Timestamp("2020-12-31")])
    def test_nanosecond_kwarg_allowed_with_pydatetime_input(self, value):
        # GH#31930 nanosecond is honored rather than dropped for pydatetime
        #  input, so it is not part of the ban
        result = Timestamp(value, nanosecond=5)
        assert result.nanosecond == 5

    @pytest.mark.parametrize(
        "value", [date(2020, 1, 1), np.datetime64("2020-01-01", "ns"), 1609372800000]
    )
    def test_nanosecond_kwarg_still_dropped_for_other_value_forms(self, value):
        # GH#31930 convert_to_tsobject only applies nanos on the pydatetime
        #  path, so nanosecond is silently ignored here.  Pinning the
        #  limitation the docstring and whatsnew call out; a separate bug.
        assert Timestamp(value, nanosecond=5).nanosecond == 0

    @pytest.mark.parametrize("kwargs", [{}, {"year": 2020}, {"year": 2020, "month": 1}])
    def test_constructor_missing_keyword(self, kwargs):
        # GH#31200

        # The exact error message of datetime() depends on its version
        if PY315:
            msg1 = (
                r"datetime\(\) missing required argument "
                r"'(year|month|day)' \(pos [123]\)"
            )
        else:
            msg1 = (
                r"function missing required argument '(year|month|day)' \(pos [123]\)"
            )
        msg2 = r"Required argument '(year|month|day)' \(pos [123]\) not found"
        msg = "|".join([msg1, msg2])

        with pytest.raises(TypeError, match=msg):
            Timestamp(**kwargs)

    def test_constructor_positional_with_tzinfo(self):
        # GH#31929
        ts = Timestamp(2020, 12, 31, tzinfo=UTC)
        expected = Timestamp("2020-12-31", tzinfo=UTC)
        assert ts == expected

    def test_constructor_positional_keyword_mixed_with_tzinfo(self):
        # GH#45307, GH#31930 nanosecond is keyword-only, so it does not
        #  compete for a positional slot
        ts = Timestamp(2020, 12, 31, tzinfo=UTC, nanosecond=4)
        expected = Timestamp("2020-12-31", tz=UTC) + Timedelta(nanoseconds=4)
        assert ts == expected

    @pytest.mark.parametrize("kwd", ["hour", "minute", "second", "microsecond"])
    def test_constructor_positional_keyword_gap_raises(self, kwd):
        # GH#31930 each parameter stands in for the *next* pydatetime field,
        #  so a keyword lands one field early; the gap it leaves in the
        #  positional arguments gives it away
        msg = "Cannot leave a gap in the positional date attributes"
        with pytest.raises(ValueError, match=msg):
            Timestamp(2020, 12, 31, **{kwd: 4})

    def test_constructor_positional_keyword_gap_raises_at_day(self):
        # GH#31930 the earliest gap, with `month` left unfilled
        msg = "Cannot leave a gap in the positional date attributes"
        with pytest.raises(ValueError, match=msg):
            Timestamp(2020, 12, day=5)

    @pytest.mark.parametrize(
        "args",
        [
            (2020, 1, 2, None, 5),
            (2020, 1, 2, None, None, None, None, UTC),
        ],
    )
    def test_constructor_positional_explicit_none_gap_raises(self, args):
        # GH#31930 an explicit None leaves the same gap a keyword would; these
        #  used to build a Timestamp with the later arguments dropped
        msg = "Cannot leave a gap in the positional date attributes"
        with pytest.raises(ValueError, match=msg):
            Timestamp(*args)

    def test_constructor_positional_keyword_no_gap(self):
        # GH#31930 a keyword that exactly fills the next free positional slot
        #  is indistinguishable from passing it positionally
        assert Timestamp(2020, 12, 31, day=5) == Timestamp(2020, 12, 31, 5)

    def test_constructor_positional_tzinfo_slot(self):
        # GH#31930 the 8th positional argument is pydatetime's tzinfo and is
        #  attached, not localized, just like datetime does it
        tz = zoneinfo.ZoneInfo("US/Eastern")
        result = Timestamp(2000, 1, 2, 3, 4, 5, 6, tz)
        assert result == Timestamp(datetime(2000, 1, 2, 3, 4, 5, 6, tz))
        assert result.tzinfo is tz

        # the 8th slot is `microsecond`, so the keyword spelling lands there too
        assert Timestamp(2000, 1, 2, 3, 4, 5, 6, microsecond=tz) == result

    def test_constructor_positional_tzinfo_slot_pytz(self):
        # GH#31930 attaching a pytz timezone gives LMT, matching pydatetime;
        #  localizing via tz=/tzinfo= would not
        pytz = pytest.importorskip("pytz")
        tz = pytz.timezone("US/Eastern")

        result = Timestamp(2000, 1, 2, 3, 4, 5, 6, tz)
        expected = datetime(2000, 1, 2, 3, 4, 5, 6, tz)
        assert result == expected
        assert result.utcoffset() == expected.utcoffset()

    def test_constructor_positional_tzinfo_slot_invalid(self):
        # GH#31930 same error as pydatetime gives for a non-tzinfo
        msg = "tzinfo argument must be None or of a tzinfo subclass"
        with pytest.raises(TypeError, match=msg):
            Timestamp(2000, 1, 2, 3, 4, 5, 6, 7)

    def test_constructor_positional_tzinfo_slot_with_tz_keyword(self):
        # GH#31930 the attached tzinfo conflicts with tz= just like it would
        #  for a tz-aware datetime
        msg = "Cannot pass a datetime or Timestamp with tzinfo with the tz parameter"
        with pytest.raises(ValueError, match=msg):
            Timestamp(2000, 1, 2, 3, 4, 5, 6, UTC, tz="UTC")

    def test_datetime_reconstruction_keeps_tzinfo(self):
        # GH#31930 CPython reconstructs datetime subclasses through the
        #  8-positional-argument form, e.g. in datetime.__replace__.
        #  Use a sub-microsecond, non-us-unit input so the assertions below
        #  pin the timezone rather than passing for unrelated reasons.
        tz = zoneinfo.ZoneInfo("US/Eastern")
        ts = Timestamp("2020-07-01 12:00:00.000000123", tz=tz)
        assert ts.unit == "ns"

        result = type(ts)(
            ts.year,
            ts.month,
            ts.day,
            ts.hour,
            ts.minute,
            ts.second,
            ts.microsecond,
            ts.tzinfo,
        )
        assert result.tz is tz
        assert result == ts.floor("us")

    @pytest.mark.skipif(not PY313, reason="copy.replace requires 3.13")
    def test_copy_replace_keeps_tzinfo(self):
        # GH#31930 copy.replace goes through datetime.__replace__, which
        #  reconstructs through the 8-argument form.  NB _Timestamp does not
        #  override __replace__, so the nanoseconds are still lost.
        tz = zoneinfo.ZoneInfo("US/Eastern")
        ts = Timestamp("2020-07-01 12:00:00.000000123", tz=tz)

        replaced = copy.replace(ts, day=2)
        assert replaced.tz is tz
        assert replaced == Timestamp("2020-07-02 12:00", tz=tz)
        assert replaced.nanosecond == 0

    @pytest.mark.skipif(not PY313, reason="copy.replace requires 3.13")
    def test_datetime_reconstruction_still_rejects_fold(self):
        # GH#31930 datetime.__replace__ passes fold alongside the 8 positional
        #  arguments, which the fold guard rejects for a non-datetime ts_input.
        #  Pinning the limitation the whatsnew calls out.
        ts = Timestamp(datetime(2020, 11, 1, 1, 30)).tz_localize(
            "US/Eastern", ambiguous=False
        )
        assert ts.fold == 1

        msg = "Cannot pass fold with possibly unambiguous input"
        with pytest.raises(ValueError, match=msg):
            copy.replace(ts, day=2)


class TestTimestampClassMethodConstructors:
    # Timestamp constructors other than __new__

    def test_utcnow_deprecated(self):
        # GH#56680
        msg = "Timestamp.utcnow is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            Timestamp.utcnow()  # noqa: TID251

    def test_utcfromtimestamp_deprecated(self):
        # GH#56680
        msg = "Timestamp.utcfromtimestamp is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            Timestamp.utcfromtimestamp(43)

    def test_constructor_strptime(self):
        # GH#25016
        # Test support for Timestamp.strptime
        fmt = "%Y%m%d-%H%M%S-%f%z"
        ts = "20190129-235348-000001+0000"
        msg = r"Timestamp.strptime\(\) is not implemented"
        with pytest.raises(NotImplementedError, match=msg):
            Timestamp.strptime(ts, fmt)

    def test_constructor_fromisocalendar(self):
        # GH#30395
        expected_timestamp = Timestamp("2000-01-03 00:00:00")
        expected_stdlib = datetime.fromisocalendar(2000, 1, 1)
        result = Timestamp.fromisocalendar(2000, 1, 1)
        assert result == expected_timestamp
        assert result == expected_stdlib
        assert isinstance(result, Timestamp)

    @pytest.mark.parametrize(
        "date_string",
        [
            "2023-11-05T08:30:00Z",
            "2023-11-05T08:30:00+0000",
            "2023-11-05T08:30:00+00:00",
            "2024-08-02T00:00:00-04:00",
        ],
    )
    def test_constructor_fromisoformat_preserves_tz(self, date_string):
        # GH#56398 fromisoformat used to drop tz info due to args offset
        #  in Timestamp.__new__ when the inherited datetime classmethod
        #  called cls(year, month, day, ..., tzinfo).
        result = Timestamp.fromisoformat(date_string)
        expected_stdlib = datetime.fromisoformat(date_string)
        assert isinstance(result, Timestamp)
        assert result.tzinfo is not None
        assert result == expected_stdlib
        assert result == Timestamp(date_string)

    def test_constructor_fromisoformat_naive(self):
        # GH#56398 the tz-naive case should remain tz-naive
        result = Timestamp.fromisoformat("2023-01-15T10:30:00")
        assert isinstance(result, Timestamp)
        assert result.tzinfo is None
        assert result == Timestamp("2023-01-15T10:30:00")

    def test_constructor_fromisoformat_roundtrip_isoformat(self):
        # GH#56398 fromisoformat(ts.isoformat()) should round-trip
        original = Timestamp.now("UTC")  # noqa: TID251
        result = Timestamp.fromisoformat(original.isoformat())
        assert result == original
        assert result.tzinfo is not None

    def test_constructor_fromordinal(self):
        base = datetime(2000, 1, 1)

        ts = Timestamp.fromordinal(base.toordinal())
        assert base == ts
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

    def test_now(self):
        # GH#9000
        ts_from_string = Timestamp("now")  # current-time-ok
        ts_from_method = Timestamp.now()  # noqa: TID251
        ts_datetime = datetime.now()  # noqa: TID251

        ts_from_string_tz = Timestamp("now", tz="US/Eastern")  # current-time-ok
        ts_from_method_tz = Timestamp.now(tz="US/Eastern")  # noqa: TID251

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
        ts_from_string = Timestamp("today")  # current-time-ok
        ts_from_method = Timestamp.today()  # noqa: TID251
        ts_datetime = datetime.today()  # noqa: TID251

        ts_from_string_tz = Timestamp("today", tz="US/Eastern")  # current-time-ok
        ts_from_method_tz = Timestamp.today(tz="US/Eastern")  # noqa: TID251

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


class TestTimestampResolutionInference:
    def test_construct_from_time_unit(self):
        # GH#54097 only passing a time component, no date
        ts = Timestamp("01:01:01.111")
        assert ts.unit == "us"

    def test_constructor_str_infer_reso(self):
        # non-iso8601 path

        # _parse_delimited_date path
        ts = Timestamp("01/30/2023")
        assert ts.unit == "us"

        # _parse_dateabbr_string path
        # GH#50907
        with tm.assert_produces_warning(
            Pandas4Warning, match="quarterly string is deprecated"
        ):
            ts = Timestamp("2015Q1")
        assert ts.unit == "us"

        # dateutil_parse path
        ts = Timestamp("2016-01-01 1:30:01 PM")
        assert ts.unit == "us"

        ts = Timestamp("2016 June 3 15:25:01.345")
        assert ts.unit == "us"

        ts = Timestamp("300-01-01")
        assert ts.unit == "us"

        ts = Timestamp("300 June 1:30:01.300")
        assert ts.unit == "us"

        # dateutil path -> don't drop trailing zeros
        ts = Timestamp("01-01-2013T00:00:00.000000000+0000")
        assert ts.unit == "ns"

        ts = Timestamp("2016/01/02 03:04:05.001000 UTC")
        assert ts.unit == "us"

        # higher-than-nanosecond -> we drop the trailing bits
        ts = Timestamp("01-01-2013T00:00:00.000000002100+0000")
        assert ts == Timestamp("01-01-2013T00:00:00.000000002+0000")
        assert ts.unit == "ns"

        # GH#56208 minute reso through the ISO8601 path with tz offset
        ts = Timestamp("2020-01-01 00:00+00:00")
        assert ts.unit == "us"

        ts = Timestamp("2020-01-01 00+00:00")
        assert ts.unit == "us"

    @pytest.mark.parametrize("method", ["now", "today"])
    def test_now_today_unit(self, method):
        # GH#55879
        ts_from_method = getattr(Timestamp, method)()
        ts_from_string = Timestamp(method)
        assert ts_from_method.unit == ts_from_string.unit == "us"


class TestTimestampConstructors:
    def test_disallow_dt64_with_weird_unit(self):
        # GH#25611
        dt64 = np.datetime64(1, "500m")
        msg = "np.datetime64 objects with units containing a multiplier"
        with pytest.raises(ValueError, match=msg):
            Timestamp(dt64)

    def test_weekday_but_no_day_raises(self):
        # GH#52659
        msg = "Parsing datetimes with weekday but no day information is not supported"
        with pytest.raises(ValueError, match=msg):
            Timestamp("2023 Sept Thu")

    def test_construct_from_string_invalid_raises(self):
        # dateutil (weirdly) parses "200622-12-31" as
        #  datetime(2022, 6, 20, 12, 0, tzinfo=tzoffset(None, -111600)
        #  which besides being mis-parsed, is a tzoffset that will cause
        #  str(ts) to raise ValueError.  Ensure we raise in the constructor
        #  instead.
        # see test_to_datetime_malformed_raise for analogous to_datetime test
        with pytest.raises(ValueError, match="gives an invalid tzoffset"):
            Timestamp("200622-12-31")

    def test_constructor_from_iso8601_str_with_offset_reso(self):
        # GH#49737
        ts = Timestamp("2016-01-01 04:05:06-01:00")
        assert ts.unit == "us"

        ts = Timestamp("2016-01-01 04:05:06.000-01:00")
        assert ts.unit == "us"

        ts = Timestamp("2016-01-01 04:05:06.000000-01:00")
        assert ts.unit == "us"

        ts = Timestamp("2016-01-01 04:05:06.000000001-01:00")
        assert ts.unit == "ns"

    def test_constructor_from_date_second_reso(self):
        # GH#49034 constructing from a pydate object gets lowest supported
        #  reso, i.e. seconds
        obj = date(2012, 9, 1)
        ts = Timestamp(obj)
        assert ts.unit == "s"

    def test_constructor_datetime64_with_tz(self):
        # GH#42288, GH#24559
        dt = np.datetime64("1970-01-01 05:00:00")
        tzstr = "UTC+05:00"

        # pre-2.0 this interpreted dt as a UTC time. in 2.0 this is treated
        #  as a wall-time, consistent with DatetimeIndex behavior
        ts = Timestamp(dt, tz=tzstr)

        alt = Timestamp(dt).tz_localize(tzstr)
        assert ts == alt
        assert ts.hour == 5

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
            (UTC, 0),
            ("Asia/Tokyo", 9),
            ("US/Eastern", -4),
            ("dateutil/US/Pacific", -7),
            (timezone(timedelta(hours=-3)), -3),
            (dateutil.tz.tzoffset(None, 18000), 5),
        ]

        for date_str, date_obj, expected in tests:
            for result in [Timestamp(date_str), Timestamp(date_obj)]:
                result = result.as_unit("ns")  # test originally written before non-nano
                # only with timestring
                assert result.as_unit("ns")._value == expected

                # re-creation shouldn't affect to internal value
                result = Timestamp(result)
                assert result.as_unit("ns")._value == expected

            # with timezone
            for tz, offset in timezones:
                for result in [Timestamp(date_str, tz=tz), Timestamp(date_obj, tz=tz)]:
                    result = result.as_unit(
                        "ns"
                    )  # test originally written before non-nano
                    expected_tz = expected - offset * 3600 * 1_000_000_000
                    assert result.as_unit("ns")._value == expected_tz

                    # should preserve tz
                    result = Timestamp(result)
                    assert result.as_unit("ns")._value == expected_tz

                    # should convert to UTC
                    if tz is not None:
                        result = Timestamp(result).tz_convert("UTC")
                    else:
                        result = Timestamp(result, tz="UTC")
                    expected_utc = expected - offset * 3600 * 1_000_000_000
                    assert result.as_unit("ns")._value == expected_utc

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
            ("UTC", 0),
            (UTC, 0),
            ("Asia/Tokyo", 9),
            ("US/Eastern", -4),
            ("dateutil/US/Pacific", -7),
            (timezone(timedelta(hours=-3)), -3),
            (dateutil.tz.tzoffset(None, 18000), 5),
        ]

        for date_str, expected in tests:
            for result in [Timestamp(date_str)]:
                # only with timestring
                assert result.as_unit("ns")._value == expected

                # re-creation shouldn't affect to internal value
                result = Timestamp(result)
                assert result.as_unit("ns")._value == expected

            # with timezone
            for tz, offset in timezones:
                result = Timestamp(date_str, tz=tz)
                expected_tz = expected
                assert result.as_unit("ns")._value == expected_tz

                # should preserve tz
                result = Timestamp(result)
                assert result.as_unit("ns")._value == expected_tz

                # should convert to UTC
                result = Timestamp(result).tz_convert("UTC")
                expected_utc = expected
                assert result.as_unit("ns")._value == expected_utc

        # This should be 2013-11-01 05:00 in UTC
        # converted to Chicago tz
        result = Timestamp("2013-11-01 00:00:00-0500", tz="America/Chicago")
        assert result._value == Timestamp("2013-11-01 05:00")._value
        expected = "Timestamp('2013-11-01 00:00:00-0500', tz='America/Chicago')"
        assert repr(result) == expected
        assert result == eval(repr(result))

        # This should be 2013-11-01 05:00 in UTC
        # converted to Tokyo tz (+09:00)
        result = Timestamp("2013-11-01 00:00:00-0500", tz="Asia/Tokyo")
        assert result._value == Timestamp("2013-11-01 05:00")._value
        expected = "Timestamp('2013-11-01 14:00:00+0900', tz='Asia/Tokyo')"
        assert repr(result) == expected
        assert result == eval(repr(result))

        # GH11708
        # This should be 2015-11-18 10:00 in UTC
        # converted to Asia/Katmandu
        result = Timestamp("2015-11-18 15:45:00+05:45", tz="Asia/Katmandu")
        assert result._value == Timestamp("2015-11-18 10:00")._value
        expected = "Timestamp('2015-11-18 15:45:00+0545', tz='Asia/Katmandu')"
        assert repr(result) == expected
        assert result == eval(repr(result))

        # This should be 2015-11-18 10:00 in UTC
        # converted to Asia/Kolkata
        result = Timestamp("2015-11-18 15:30:00+05:30", tz="Asia/Kolkata")
        assert result._value == Timestamp("2015-11-18 10:00")._value
        expected = "Timestamp('2015-11-18 15:30:00+0530', tz='Asia/Kolkata')"
        assert repr(result) == expected
        assert result == eval(repr(result))

    def test_constructor_invalid(self):
        msg = "Cannot convert input"
        with pytest.raises(TypeError, match=msg):
            Timestamp(slice(2))
        msg = "Cannot convert Period"
        with pytest.raises(ValueError, match=msg):
            Timestamp(Period("1000-01-01"))

    def test_constructor_invalid_tz(self):
        # GH#17690
        msg = (
            "Argument 'tzinfo' has incorrect type "
            r"\(expected datetime.tzinfo, got str\)"
        )
        with pytest.raises(TypeError, match=msg):
            Timestamp("2017-10-22", tzinfo="US/Eastern")

        msg = "at most one of"
        with pytest.raises(ValueError, match=msg):
            Timestamp("2017-10-22", tzinfo=UTC, tz="UTC")

        msg = "Cannot pass a date attribute keyword argument when passing a date string"
        with pytest.raises(ValueError, match=msg):
            # GH#5168
            # case where user tries to pass tz as an arg, not kwarg, gets
            # interpreted as `year`
            Timestamp("2012-01-01", "US/Pacific")

    def test_constructor_tz_or_tzinfo(self):
        # GH#17943, GH#17690, GH#5168
        stamps = [
            Timestamp(year=2017, month=10, day=22, tz="UTC"),
            Timestamp(year=2017, month=10, day=22, tzinfo=UTC),
            Timestamp(year=2017, month=10, day=22, tz=UTC),
            Timestamp(datetime(2017, 10, 22), tzinfo=UTC),
            Timestamp(datetime(2017, 10, 22), tz="UTC"),
            Timestamp(datetime(2017, 10, 22), tz=UTC),
        ]
        assert all(ts == stamps[0] for ts in stamps)

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
            Timestamp(2000, 1, 2, 3, 4, 5, 6, None, nanosecond=1),
            Timestamp(2000, 1, 2, 3, 4, 5, 6, tz=UTC, nanosecond=1),
        ],
    )
    def test_constructor_nanosecond(self, result):
        # GH 18898
        # As of 2.0 (GH 49416), nanosecond should not be accepted positionally
        expected = Timestamp(datetime(2000, 1, 2, 3, 4, 5, 6), tz=result.tz)
        expected = expected + Timedelta(nanoseconds=1)
        assert result == expected

    @pytest.mark.parametrize("z", ["Z0", "Z00"])
    def test_constructor_invalid_Z0_isostring(self, z):
        # GH 8910
        msg = f"Unknown datetime string format, unable to parse: 2014-11-02 01:00{z}"
        with pytest.raises(ValueError, match=msg):
            Timestamp(f"2014-11-02 01:00{z}")

    def test_out_of_bounds_integer_value(self):
        # GH#26651 check that we raise OutOfBoundsDatetime, not OverflowError
        msg = str(Timestamp.max._value * 2)
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            Timestamp(Timestamp.max._value * 2)
        msg = str(Timestamp.min._value * 2)
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            Timestamp(Timestamp.min._value * 2)

    def test_out_of_bounds_value(self):
        one_us = np.timedelta64(1, "ns").astype("timedelta64[us]")

        # By definition we can't go out of bounds in [ns], so we
        # convert the datetime64s to [us] so we can go out of bounds
        min_ts_us = np.datetime64(Timestamp.min).astype("M8[us]") + one_us
        max_ts_us = np.datetime64(Timestamp.max).astype("M8[us]")

        # No error for the min/max datetimes
        Timestamp(min_ts_us)
        Timestamp(max_ts_us)

        # We used to raise on these before supporting non-nano
        us_val = NpyDatetimeUnit.NPY_FR_us.value
        assert Timestamp(min_ts_us - one_us)._creso == us_val
        assert Timestamp(max_ts_us + one_us)._creso == us_val

        # https://github.com/numpy/numpy/issues/22346 for why
        #  we can't use the same construction as above with minute resolution

        # too_low, too_high are the _just_ outside the range of M8[s]
        too_low = np.datetime64("-292277022657-01-27T08:29", "m")
        too_high = np.datetime64("292277026596-12-04T15:31", "m")

        msg = "Out of bounds"
        # One us less than the minimum is an error
        with pytest.raises(ValueError, match=msg):
            Timestamp(too_low)

        # One us more than the maximum is an error
        with pytest.raises(ValueError, match=msg):
            Timestamp(too_high)

    def test_out_of_bounds_string(self):
        msg = "Cannot cast .* to unit='ns' without overflow"
        with pytest.raises(ValueError, match=msg):
            Timestamp("1676-01-01").as_unit("ns")
        with pytest.raises(ValueError, match=msg):
            Timestamp("2263-01-01").as_unit("ns")

        ts = Timestamp("2263-01-01")
        assert ts.unit == "us"

        ts = Timestamp("1676-01-01")
        assert ts.unit == "us"

    def test_barely_out_of_bounds(self):
        # GH#19529
        # GH#19382 close enough to bounds that dropping nanos would result
        # in an in-bounds datetime
        msg = "Out of bounds nanosecond timestamp: 2262-04-11 23:47:16"
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            Timestamp("2262-04-11 23:47:16.854775808")

    def test_bounds_with_different_units(self):
        out_of_bounds_dates = ("1677-09-21", "2262-04-12")

        time_units = ("D", "h", "m", "s", "ms", "us")

        for date_string in out_of_bounds_dates:
            for unit in time_units:
                dt64 = np.datetime64(date_string, unit)
                ts = Timestamp(dt64)
                if unit in ["s", "ms", "us"]:
                    # We can preserve the input unit
                    assert ts._value == dt64.view("i8")
                else:
                    # we chose the closest unit that we _do_ support
                    assert ts._creso == NpyDatetimeUnit.NPY_FR_s.value

        # With more extreme cases, we can't even fit inside second resolution
        info = np.iinfo(np.int64)
        msg = "Out of bounds second timestamp:"
        for value in [info.min + 1, info.max]:
            for unit in ["D", "h", "m"]:
                dt64 = np.datetime64(value, unit)
                with pytest.raises(OutOfBoundsDatetime, match=msg):
                    Timestamp(dt64)

        in_bounds_dates = ("1677-09-23", "2262-04-11")

        for date_string in in_bounds_dates:
            for unit in time_units:
                dt64 = np.datetime64(date_string, unit)
                Timestamp(dt64)

    @pytest.mark.parametrize("arg", ["001-01-01", "0001-01-01"])
    def test_out_of_bounds_string_consistency(self, arg):
        # GH 15829
        msg = "Cannot cast 0001-01-01 00:00:00 to unit='ns' without overflow"
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            Timestamp(arg).as_unit("ns")

        ts = Timestamp(arg)
        assert ts.unit == "us"
        assert ts.year == ts.month == ts.day == 1

    def test_min_valid(self):
        # Ensure that Timestamp.min is a valid Timestamp
        Timestamp(Timestamp.min)

    def test_max_valid(self):
        # Ensure that Timestamp.max is a valid Timestamp
        Timestamp(Timestamp.max)

    @pytest.mark.parametrize("offset", ["+0300", "+0200"])
    def test_construct_timestamp_near_dst(self, offset):
        # GH 20854
        expected = Timestamp(f"2016-10-30 03:00:00{offset}", tz="Europe/Helsinki")
        result = Timestamp(expected).tz_convert("Europe/Helsinki")
        assert result == expected

    @pytest.mark.parametrize(
        "arg", ["2013/01/01 00:00:00+09:00", "2013-01-01 00:00:00+09:00"]
    )
    def test_construct_with_different_string_format(self, arg):
        # GH 12064
        result = Timestamp(arg)
        expected = Timestamp(datetime(2013, 1, 1), tz=timezone(timedelta(hours=9)))
        assert result == expected

    def test_comma_separated_fractional_seconds(self):
        # GH#59256 - comma as decimal separator should parse correctly
        ts_comma = Timestamp("2024-06-17T18:57:43,567")
        ts_dot = Timestamp("2024-06-17T18:57:43.567")
        assert ts_comma == ts_dot
        assert ts_comma.unit == ts_dot.unit
        assert ts_comma.microsecond == 567000

    @pytest.mark.parametrize("box", [datetime, Timestamp])
    def test_raise_tz_and_tzinfo_in_datetime_input(self, box):
        # GH 23579
        kwargs = {"year": 2018, "month": 1, "day": 1, "tzinfo": UTC}
        msg = "Cannot pass a datetime or Timestamp"
        with pytest.raises(ValueError, match=msg):
            Timestamp(box(**kwargs), tz="US/Pacific")
        msg = "Cannot pass a datetime or Timestamp"
        with pytest.raises(ValueError, match=msg):
            Timestamp(box(**kwargs), tzinfo=zoneinfo.ZoneInfo("US/Pacific"))

    def test_dont_convert_dateutil_utc_to_default_utc(self):
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

    def test_timestamp_constructor_tz_utc(self):
        utc_stamp = Timestamp("3/11/2012 05:00", tz="utc")
        assert utc_stamp.tzinfo is UTC
        assert utc_stamp.hour == 5

        utc_stamp = Timestamp("3/11/2012 05:00").tz_localize("utc")
        assert utc_stamp.hour == 5

    def test_timestamp_to_datetime_tzoffset(self):
        tzinfo = tzoffset(None, 7200)
        expected = Timestamp("3/11/2012 04:00", tz=tzinfo)
        result = Timestamp(expected.to_pydatetime())
        assert expected == result

    def test_timestamp_constructor_near_dst_boundary(self):
        # GH#11481 & GH#15777
        # Naive string timestamps were being localized incorrectly
        # with tz_convert_from_utc_single instead of tz_localize_to_utc

        for tz in ["Europe/Brussels", "Europe/Prague"]:
            result = Timestamp("2015-10-25 01:00", tz=tz)
            expected = Timestamp("2015-10-25 01:00").tz_localize(tz)
            assert result == expected

            msg = "Cannot infer dst time from 2015-10-25 02:00:00"
            with pytest.raises(ValueError, match=msg):
                Timestamp("2015-10-25 02:00", tz=tz)

        result = Timestamp("2017-03-26 01:00", tz="Europe/Paris")
        expected = Timestamp("2017-03-26 01:00").tz_localize("Europe/Paris")
        assert result == expected

        msg = "2017-03-26 02:00"
        with pytest.raises(ValueError, match=msg):
            Timestamp("2017-03-26 02:00", tz="Europe/Paris")

        # GH#11708
        naive = Timestamp("2015-11-18 10:00:00")
        result = naive.tz_localize("UTC").tz_convert("Asia/Kolkata")
        expected = Timestamp("2015-11-18 15:30:00+0530", tz="Asia/Kolkata")
        assert result == expected

        # GH#15823
        result = Timestamp("2017-03-26 00:00", tz="Europe/Paris")
        expected = Timestamp("2017-03-26 00:00:00+0100", tz="Europe/Paris")
        assert result == expected

        result = Timestamp("2017-03-26 01:00", tz="Europe/Paris")
        expected = Timestamp("2017-03-26 01:00:00+0100", tz="Europe/Paris")
        assert result == expected

        msg = "2017-03-26 02:00"
        with pytest.raises(ValueError, match=msg):
            Timestamp("2017-03-26 02:00", tz="Europe/Paris")

        result = Timestamp("2017-03-26 02:00:00+0100", tz="Europe/Paris")
        naive = Timestamp(result.as_unit("ns")._value)
        expected = naive.tz_localize("UTC").tz_convert("Europe/Paris")
        assert result == expected

        result = Timestamp("2017-03-26 03:00", tz="Europe/Paris")
        expected = Timestamp("2017-03-26 03:00:00+0200", tz="Europe/Paris")
        assert result == expected

    @pytest.mark.parametrize(
        "tz",
        [
            "pytz/US/Eastern",
            gettz("US/Eastern"),
            "US/Eastern",
            "dateutil/US/Eastern",
        ],
    )
    def test_timestamp_constructed_by_date_and_tz(self, tz):
        # GH#2993, Timestamp cannot be constructed by datetime.date
        # and tz correctly
        if isinstance(tz, str) and tz.startswith("pytz/"):
            pytz = pytest.importorskip("pytz")
            tz = pytz.timezone(tz.removeprefix("pytz/"))
        result = Timestamp(date(2012, 3, 11), tz=tz)

        expected = Timestamp("3/11/2012", tz=tz)
        assert result.hour == expected.hour
        assert result == expected

    def test_explicit_tz_none(self):
        # GH#48688
        msg = "Passed data is timezone-aware, incompatible with 'tz=None'"
        with pytest.raises(ValueError, match=msg):
            Timestamp(datetime(2022, 1, 1, tzinfo=UTC), tz=None)

        with pytest.raises(ValueError, match=msg):
            Timestamp("2022-01-01 00:00:00", tzinfo=UTC, tz=None)

        with pytest.raises(ValueError, match=msg):
            Timestamp("2022-01-01 00:00:00-0400", tz=None)


def test_constructor_ambiguous_dst():
    # GH 24329
    # Make sure that calling Timestamp constructor
    # on Timestamp created from ambiguous time
    # doesn't change Timestamp.value
    ts = Timestamp(1382835600000000000, tz="dateutil/Europe/London")
    expected = ts._value
    result = Timestamp(ts)._value
    assert result == expected


@pytest.mark.parametrize("epoch", [1552211999999999872, 1552211999999999999])
def test_constructor_before_dst_switch(epoch):
    # GH 31043
    # Make sure that calling Timestamp constructor
    # on time just before DST switch doesn't lead to
    # nonexistent time or value change
    ts = Timestamp(epoch, tz="dateutil/America/Los_Angeles")
    result = ts.tz.dst(ts)
    expected = timedelta(seconds=0)
    assert Timestamp(ts)._value == epoch
    assert result == expected


def test_timestamp_constructor_identity():
    # Test for #30543
    expected = Timestamp("2017-01-01T12")
    result = Timestamp(expected)
    assert result is expected


@pytest.mark.parametrize("nano", [-1, 1000])
def test_timestamp_nano_range(nano):
    # GH 48255
    with pytest.raises(ValueError, match="nanosecond must be in 0..999"):
        Timestamp(year=2022, month=1, day=1, nanosecond=nano)


def test_non_nano_value():
    # https://github.com/pandas-dev/pandas/issues/49076
    msg = "The 'unit' keyword is only used when"
    with tm.assert_produces_warning(UserWarning, match=msg):
        result = Timestamp("1800-01-01", unit="s").value
    # `.value` shows nanoseconds, even though unit is 's'
    assert result == -5364662400000000000

    # out-of-nanoseconds-bounds `.value` raises informative message
    msg = (
        r"Cannot convert Timestamp to nanoseconds without overflow. "
        r"Use `.asm8.view\('i8'\)` to cast represent Timestamp in its "
        r"own unit \(here, us\).$"
    )
    ts = Timestamp("0300-01-01")
    assert ts.unit == "us"
    with pytest.raises(OverflowError, match=msg):
        ts.value
    # check that the suggested workaround actually works
    result = ts.asm8.view("i8")
    assert result == -52_700_112_000 * 10**6


def test_timestamp_constructor_np_str():
    # GH#48974 np.str_ should not break Timestamp construction
    result = Timestamp(np.str_("2023-01-01"))
    expected = Timestamp("2023-01-01")
    assert result == expected

    result = Timestamp(np.str_("2023-01-01 12:34:56"))
    expected = Timestamp("2023-01-01 12:34:56")
    assert result == expected


@pytest.mark.parametrize(
    "na_value", [None, np.nan, np.datetime64("NaT", "ns"), NaT, NA]
)
def test_timestamp_constructor_na_value(na_value):
    # GH45481
    result = Timestamp(na_value)
    expected = NaT
    assert result is expected


@pytest.mark.parametrize(
    "na_value", [None, np.nan, np.datetime64("NaT", "ns"), NaT, NA]
)
@pytest.mark.parametrize(
    "kwargs",
    [
        {"year": 2020},
        {"month": 12},
        {"day": 31},
        {"hour": 5},
        {"minute": 5},
        {"second": 5},
        {"microsecond": 5},
        {"nanosecond": 5},
        {"tz": "UTC"},
        {"year": 2020, "month": 12, "day": 31, "hour": 5},
    ],
)
def test_timestamp_constructor_na_value_with_components(na_value, kwargs):
    # GH#31930 null-like input is NaT for every combination of the remaining
    #  arguments; in particular `year` must not divert to the positional form
    assert Timestamp(na_value, **kwargs) is NaT


@pytest.mark.parametrize(
    "na_value", [None, np.nan, np.datetime64("NaT", "ns"), NaT, NA]
)
def test_timestamp_constructor_na_value_still_validates(na_value):
    # GH#31930 returning NaT for null-like input must not swallow arguments
    #  that are invalid on their own
    with pytest.raises(ValueError, match="nanosecond must be in 0..999"):
        Timestamp(na_value, nanosecond=5000)

    with pytest.raises(zoneinfo.ZoneInfoNotFoundError, match="Not/AZone"):
        Timestamp(na_value, tz="Not/AZone")


@pytest.mark.parametrize("tz", [None, "UTC"])
def test_timestamp_iso8601_offset_out_of_bounds(tz):
    # GH#65353 shifting to UTC takes this past Timestamp.max; it must raise
    # rather than silently wrap around to a bogus in-bounds datetime
    msg = "Out of bounds nanosecond timestamp"
    with pytest.raises(OutOfBoundsDatetime, match=msg):
        Timestamp("2262-04-11T20:00:00.000000000-11:00", tz=tz)


# GH#66510 iNaT == INT64_MIN, one below Timestamp.min. A value that renders or
#  shifts onto it is not NaT, but is indistinguishable from NaT once stored in a
#  datetime64 array, so it has to be rejected instead.
SENTINEL_UTC = "1677-09-21 00:12:43.145224192"
SENTINEL_PLUS_1H = "1677-09-21 01:12:43.145224192"
PLUS_1H = timezone(timedelta(hours=1))
MINUS_1H = timezone(timedelta(hours=-1))


def test_constructor_naive_datetime_renders_nat_sentinel():
    # GH#66510 dtstruct -> int64 rendering lands on iNaT
    msg = re.escape(f"Out of bounds nanosecond timestamp: {SENTINEL_UTC}")
    with pytest.raises(OutOfBoundsDatetime, match=msg):
        Timestamp(datetime(1677, 9, 21, 0, 12, 43, 145224), nanosecond=192)


def test_constructor_datetime64_tz_shift_hits_nat_sentinel():
    # GH#66510 np.datetime64 is treated as a wall time; localizing it lands on iNaT
    msg = re.escape(f"Out of bounds nanosecond timestamp: {SENTINEL_PLUS_1H}")
    with pytest.raises(OutOfBoundsDatetime, match=msg):
        Timestamp(np.datetime64(SENTINEL_PLUS_1H.replace(" ", "T")), tz=PLUS_1H)


def test_constructor_aware_datetime_offset_hits_nat_sentinel():
    # GH#66510 subtracting the UTC offset lands on iNaT
    msg = re.escape(f"Out of bounds nanosecond timestamp: {SENTINEL_PLUS_1H}")
    with pytest.raises(OutOfBoundsDatetime, match=msg):
        Timestamp(
            datetime(1677, 9, 21, 1, 12, 43, 145224, tzinfo=PLUS_1H), nanosecond=192
        )


@pytest.mark.parametrize(
    "arg, kwargs",
    [
        # the tz keyword shifts via tz_localize_to_utc_single
        (SENTINEL_PLUS_1H, {"tz": PLUS_1H}),
        # an embedded offset shifts via the checked_sub fastpath
        (f"{SENTINEL_PLUS_1H}+01:00", {}),
    ],
)
def test_constructor_str_shift_hits_nat_sentinel(arg, kwargs):
    # GH#66510 both shift paths in convert_str_to_tsobject
    msg = re.escape(f"Out of bounds nanosecond timestamp: {SENTINEL_PLUS_1H}")
    with pytest.raises(OutOfBoundsDatetime, match=msg):
        Timestamp(arg, **kwargs)


def test_constructor_aware_datetime_offset_overflow():
    # GH#66510 the offset subtraction used to raise a bare OverflowError
    msg = re.escape("Out of bounds nanosecond timestamp: 2262-04-11 23:47:16.854775807")
    tz = timezone(timedelta(hours=-23, minutes=-59))
    with pytest.raises(OutOfBoundsDatetime, match=msg):
        Timestamp(datetime(2262, 4, 11, 23, 47, 16, 854775, tzinfo=tz), nanosecond=807)


def test_constructor_offset_shifts_off_nat_sentinel():
    # GH#66510 the wall time renders onto iNaT but the shift to UTC moves it back
    #  in bounds, so the rejection has to happen after the shift, not before
    ts = Timestamp(
        datetime(1677, 9, 21, 0, 12, 43, 145224, tzinfo=MINUS_1H), nanosecond=192
    )
    assert ts._value == -(2**63) + 3600 * 10**9


def test_constructor_nat_sentinel_neighbours():
    # GH#66510 only the sentinel itself is out of bounds
    assert Timestamp("1677-09-21 00:12:43.145224193")._value == -(2**63) + 1
    assert Timestamp.min._value == -(2**63) + 1
    assert Timestamp(NaT) is NaT


@pytest.mark.parametrize(
    "date_str", ["2014Q2", "2014-Q2", "2Q2014", "2Q-2014", "2Q14", "14Q2"]
)
def test_constructor_quarterly_string_deprecated(date_str):
    # GH#50907
    with tm.assert_produces_warning(
        Pandas4Warning, match="quarterly string is deprecated"
    ):
        result = Timestamp(date_str)

    with tm.assert_produces_warning(None):
        expected = Period(date_str).to_timestamp()
    assert result == expected
