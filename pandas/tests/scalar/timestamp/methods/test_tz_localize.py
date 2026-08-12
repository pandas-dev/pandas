from datetime import (
    timedelta,
    timezone,
)
import re
import zoneinfo

import dateutil.tz
from dateutil.tz import gettz
import numpy as np
import pytest

from pandas._libs.tslibs.dtypes import NpyDatetimeUnit
from pandas.compat import WASM
from pandas.errors import OutOfBoundsDatetime
import pandas.util._test_decorators as td

from pandas import (
    NaT,
    Timedelta,
    Timestamp,
)
import pandas._testing as tm


def _distant_date_only_for_zoneinfo(year, tz):
    if int(year) >= 2200:
        if isinstance(tz, zoneinfo.ZoneInfo):
            return
        elif isinstance(tz, str) and not tz.startswith(("dateutil/", "pytz/")):
            # zoneinfo tz passed as string
            return
        else:
            pytest.skip("ZoneInfo required for year >= 2200")


class TestTimestampTZLocalize:
    def test_tz_localize_pushes_out_of_bounds(self):
        # GH#12677
        # tz_localize that pushes away from the boundary is OK
        pytz = pytest.importorskip("pytz")
        msg = (
            f"Converting {Timestamp.min.strftime('%Y-%m-%d %H:%M:%S')} "
            f"underflows past {Timestamp.min}"
        )
        pac = Timestamp.min.tz_localize(pytz.timezone("US/Pacific"))
        assert pac._value > Timestamp.min._value
        pac.tz_convert("Asia/Tokyo")  # tz_convert doesn't change value
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            Timestamp.min.tz_localize(pytz.timezone("Asia/Tokyo"))

        # tz_localize that pushes away from the boundary is OK
        msg = (
            f"Converting {Timestamp.max.strftime('%Y-%m-%d %H:%M:%S')} "
            f"overflows past {Timestamp.max}"
        )
        tokyo = Timestamp.max.tz_localize(pytz.timezone("Asia/Tokyo"))
        assert tokyo._value < Timestamp.max._value
        tokyo.tz_convert("US/Pacific")  # tz_convert doesn't change value
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            Timestamp.max.tz_localize(pytz.timezone("US/Pacific"))

    def test_tz_localize_shift_onto_nat_sentinel(self, unit):
        # GH#66550 Asia/Tokyo's pre-1888 LMT offset is +9:18:59, so this wall
        #  time shifts to exactly the NaT sentinel.  That is one below the
        #  minimum representable value, so it underflows; previously it was
        #  misreported as a nonexistent time.
        lmt = np.timedelta64(9 * 3600 + 18 * 60 + 59, "s")
        lmt_offset = int(lmt // np.timedelta64(1, unit))
        ts = Timestamp(np.datetime64(-(2**63) + lmt_offset, unit))
        with pytest.raises(OutOfBoundsDatetime, match="underflows past"):
            ts.tz_localize("Asia/Tokyo")
        with pytest.raises(OutOfBoundsDatetime, match="underflows past"):
            Timestamp(ts.to_datetime64(), tz="Asia/Tokyo")

    def test_tz_localize_fixed_offset_shift_onto_nat_sentinel(self):
        # GH#66550 same sentinel underflow via the scalar fixed-offset branch.  The
        #  message matches the Timestamp(..., tz=) constructor, which shifts through
        #  the same helper and already rejected this.
        ts = Timestamp(-(2**63) + 9 * 3600 * 10**9)
        msg = re.escape(f"Out of bounds nanosecond timestamp: {ts}")

        with pytest.raises(OutOfBoundsDatetime, match=msg):
            ts.tz_localize(timezone(timedelta(hours=9)))

    @td.skip_if_windows  # `tm.set_timezone` does not work on Windows
    @pytest.mark.skipif(WASM, reason="`tm.set_timezone` does not work on WASM")
    def test_tz_localize_tzlocal_shift_onto_nat_sentinel(self):
        # GH#66550 the tzlocal branch reaches the same guard as the fixed-offset
        #  one above; pin the ambient zone to +9 so the shift lands on the sentinel
        ts = Timestamp(-(2**63) + 9 * 3600 * 10**9)
        msg = re.escape(f"Out of bounds nanosecond timestamp: {ts}")

        with tm.set_timezone("Asia/Tokyo"):
            with pytest.raises(OutOfBoundsDatetime, match=msg):
                ts.tz_localize(dateutil.tz.tzlocal())

    def test_tz_localize_overflow_past_last_cached_transition(self):
        # GH#65733 wall times past the last cached DST transition are localized
        #  through the zoneinfo fallback.  A candidate UTC instant that overflows
        #  is out of bounds, but the unguarded subtraction wrapped instead, so the
        #  round-trip check failed and the wall time was reported as nonexistent.
        ts = Timestamp(Timestamp.max._value - 2 * 3600 * 10**9)

        with pytest.raises(OutOfBoundsDatetime, match="overflows past"):
            ts.tz_localize("America/New_York")
        # ...and with nonexistent="NaT" the misclassification silently won out
        with pytest.raises(OutOfBoundsDatetime, match="overflows past"):
            ts.tz_localize("America/New_York", nonexistent="NaT")

        # a wall time whose UTC instant still fits is unaffected
        ok = Timestamp(Timestamp.max._value - 5 * 3600 * 10**9)
        assert ok.tz_localize("America/New_York").utcoffset() == -timedelta(hours=4)

    @pytest.mark.parametrize(
        "tz",
        [zoneinfo.ZoneInfo("US/Central"), "dateutil/US/Central", "pytz/US/Central"],
    )
    def test_tz_localize_ambiguous_bool(self, unit, tz):
        # make sure that we are correctly accepting bool values as ambiguous
        # GH#14402
        if isinstance(tz, str) and tz.startswith("pytz/"):
            pytz = pytest.importorskip("pytz")
            tz = pytz.timezone(tz.removeprefix("pytz/"))
        ts = Timestamp("2015-11-01 01:00:03").as_unit(unit)
        expected0 = Timestamp("2015-11-01 01:00:03-0500", tz=tz)
        expected1 = Timestamp("2015-11-01 01:00:03-0600", tz=tz)

        msg = "Cannot infer dst time from 2015-11-01 01:00:03"
        with pytest.raises(ValueError, match=msg):
            ts.tz_localize(tz)

        result = ts.tz_localize(tz, ambiguous=True)
        assert result == expected0
        assert result._creso == getattr(NpyDatetimeUnit, f"NPY_FR_{unit}").value

        result = ts.tz_localize(tz, ambiguous=False)
        assert result == expected1
        assert result._creso == getattr(NpyDatetimeUnit, f"NPY_FR_{unit}").value

    @pytest.mark.parametrize("year", ["2014", "2200"])
    def test_tz_localize_ambiguous(self, year):
        ts = Timestamp(f"{year}-11-02 01:00").as_unit("s")
        ts_dst = ts.tz_localize("US/Eastern", ambiguous=True)
        ts_no_dst = ts.tz_localize("US/Eastern", ambiguous=False)

        assert ts_no_dst._value - ts_dst._value == 3600
        msg = re.escape(
            "'ambiguous' parameter must be one of: "
            "True, False, 'NaT', 'raise' (default)"
        )
        with pytest.raises(ValueError, match=msg):
            ts.tz_localize("US/Eastern", ambiguous="infer")

    def test_tz_localize_invalid(self):
        # GH#8025
        msg = "Cannot localize tz-aware Timestamp, use tz_convert for conversions"
        with pytest.raises(TypeError, match=msg):
            Timestamp("2011-01-01", tz="US/Eastern").tz_localize("Asia/Tokyo")

        msg = "Cannot convert tz-naive Timestamp, use tz_localize to localize"
        with pytest.raises(TypeError, match=msg):
            Timestamp("2011-01-01").tz_convert("Asia/Tokyo")

    @pytest.mark.parametrize("year", ["2015", "2201"])
    @pytest.mark.parametrize(
        "stamp, tz",
        [
            ("2015-03-08 02:00", "US/Eastern"),
            ("2015-03-08 02:30", "US/Pacific"),
            ("2015-03-29 02:00", "Europe/Paris"),
            ("2015-03-29 02:30", "Europe/Belgrade"),
        ],
    )
    def test_tz_localize_nonexistent(self, stamp, tz, year):
        # GH#13057
        stamp = stamp.replace("2015", year)
        ts = Timestamp(stamp)
        with pytest.raises(ValueError, match=stamp):
            ts.tz_localize(tz)
        # GH 22644
        with pytest.raises(ValueError, match=stamp):
            ts.tz_localize(tz, nonexistent="raise")
        assert ts.tz_localize(tz, nonexistent="NaT") is NaT

    @pytest.mark.parametrize(
        "stamp, tz, forward_expected, backward_expected",
        [
            (
                "2015-03-29 02:00:00",
                "Europe/Warsaw",
                "2015-03-29 03:00:00",
                "2015-03-29 01:59:59.999999",
            ),  # utc+1 -> utc+2
            (
                "2023-03-12 02:00:00",
                "America/Los_Angeles",
                "2023-03-12 03:00:00",
                "2023-03-12 01:59:59.999999",
            ),  # utc-8 -> utc-7
            (
                "2023-03-26 01:00:00",
                "Europe/London",
                "2023-03-26 02:00:00",
                "2023-03-26 00:59:59.999999",
            ),  # utc+0 -> utc+1
            (
                "2023-03-26 00:00:00",
                "Atlantic/Azores",
                "2023-03-26 01:00:00",
                "2023-03-25 23:59:59.999999",
            ),  # utc-1 -> utc+0
        ],
    )
    def test_tz_localize_nonexistent_shift(
        self, stamp, tz, forward_expected, backward_expected
    ):
        ts = Timestamp(stamp)
        forward_ts = ts.tz_localize(tz, nonexistent="shift_forward")
        assert forward_ts == Timestamp(forward_expected, tz=tz)
        backward_ts = ts.tz_localize(tz, nonexistent="shift_backward")
        assert backward_ts == Timestamp(backward_expected, tz=tz)

    def test_tz_localize_ambiguous_raise(self):
        # GH#13057
        ts = Timestamp("2015-11-1 01:00")
        msg = "Cannot infer dst time from 2015-11-01 01:00:00,"
        with pytest.raises(ValueError, match=msg):
            ts.tz_localize("US/Pacific", ambiguous="raise")

    def test_tz_localize_nonexistent_invalid_arg(self, warsaw):
        # GH 22644
        tz = warsaw
        ts = Timestamp("2015-03-29 02:00:00")
        msg = (
            "The nonexistent argument must be one of 'raise', 'NaT', "
            "'shift_forward', 'shift_backward' or a timedelta object"
        )
        with pytest.raises(ValueError, match=msg):
            ts.tz_localize(tz, nonexistent="foo")

    @pytest.mark.parametrize(
        "stamp",
        [
            "2014-02-01 09:00",
            "2014-07-08 09:00",
            "2014-11-01 17:00",
            "2014-11-05 00:00",
        ],
    )
    def test_tz_localize_roundtrip(self, stamp, tz_aware_fixture):
        tz = tz_aware_fixture
        ts = Timestamp(stamp)
        localized = ts.tz_localize(tz)
        assert localized == Timestamp(stamp, tz=tz)

        msg = "Cannot localize tz-aware Timestamp"
        with pytest.raises(TypeError, match=msg):
            localized.tz_localize(tz)

        reset = localized.tz_localize(None)
        assert reset == ts
        assert reset.tzinfo is None

    def test_tz_localize_ambiguous_compat(self):
        # validate that pytz and dateutil are compat for dst
        # when the transition happens
        pytz = pytest.importorskip("pytz")
        naive = Timestamp("2013-10-27 01:00:00")
        assert naive.unit == "us"

        pytz_zone = pytz.timezone("Europe/London")
        dateutil_zone = "dateutil/Europe/London"
        result_pytz = naive.tz_localize(pytz_zone, ambiguous=False)
        result_dateutil = naive.tz_localize(dateutil_zone, ambiguous=False)
        assert result_pytz._value == result_dateutil._value
        assert result_pytz._value == 1_382_835_600 * 10**6

        # fixed ambiguous behavior
        # see gh-14621, GH#45087
        assert result_pytz.to_pydatetime().tzname() == "GMT"
        assert result_dateutil.to_pydatetime().tzname() == "GMT"
        assert str(result_pytz) == str(result_dateutil)

        # 1 hour difference
        result_pytz = naive.tz_localize(pytz_zone, ambiguous=True)
        result_dateutil = naive.tz_localize(dateutil_zone, ambiguous=True)
        assert result_pytz._value == result_dateutil._value
        assert result_pytz._value == 1_382_832_000 * 10**6

        # see gh-14621
        assert str(result_pytz) == str(result_dateutil)
        assert (
            result_pytz.to_pydatetime().tzname()
            == result_dateutil.to_pydatetime().tzname()
        )

    @pytest.mark.parametrize(
        "tz",
        [
            "pytz/US/Eastern",
            gettz("US/Eastern"),
            zoneinfo.ZoneInfo("US/Eastern"),
            "dateutil/US/Eastern",
        ],
    )
    def test_timestamp_tz_localize(self, tz):
        if isinstance(tz, str) and tz.startswith("pytz/"):
            pytz = pytest.importorskip("pytz")
            tz = pytz.timezone(tz.removeprefix("pytz/"))
        stamp = Timestamp("3/11/2012 04:00")

        result = stamp.tz_localize(tz)
        expected = Timestamp("3/11/2012 04:00", tz=tz)
        assert result.hour == expected.hour
        assert result == expected

    @pytest.mark.parametrize("year", ["2018", "2204"])
    @pytest.mark.parametrize(
        "start_ts, tz, end_ts, shift",
        [
            ["2018-03-25 02:20:00", "Europe/Warsaw", "2018-03-25 03:00:00", "forward"],
            [
                "2018-03-25 02:20:00",
                "Europe/Warsaw",
                "2018-03-25 01:59:59.999999999",
                "backward",
            ],
            [
                "2018-03-25 02:20:00",
                "Europe/Warsaw",
                "2018-03-25 03:20:00",
                timedelta(hours=1),
            ],
            [
                "2018-03-25 02:20:00",
                "Europe/Warsaw",
                "2018-03-25 01:20:00",
                timedelta(hours=-1),
            ],
            ["2018-03-11 02:33:00", "US/Pacific", "2018-03-11 03:00:00", "forward"],
            [
                "2018-03-11 02:33:00",
                "US/Pacific",
                "2018-03-11 01:59:59.999999999",
                "backward",
            ],
            [
                "2018-03-11 02:33:00",
                "US/Pacific",
                "2018-03-11 03:33:00",
                timedelta(hours=1),
            ],
            [
                "2018-03-11 02:33:00",
                "US/Pacific",
                "2018-03-11 01:33:00",
                timedelta(hours=-1),
            ],
        ],
    )
    @pytest.mark.parametrize("tz_type", ["", "dateutil/"])
    def test_timestamp_tz_localize_nonexistent_shift(
        self, start_ts, tz, end_ts, shift, year, tz_type, unit
    ):
        # GH 8917, 24466
        if year == "2204" and tz_type == "dateutil/":
            pytest.skip("ZoneInfo required for year >= 2200")
        start_ts = start_ts.replace("2018", year)
        end_ts = end_ts.replace("2018", year)
        tz = tz_type + tz
        if isinstance(shift, str):
            shift = "shift_" + shift
        ts = Timestamp(start_ts).as_unit(unit)
        result = ts.tz_localize(tz, nonexistent=shift)
        expected = Timestamp(end_ts).tz_localize(tz)

        if unit == "us":
            assert result == expected.replace(nanosecond=0)
        elif unit == "ms":
            micros = expected.microsecond - expected.microsecond % 1000
            assert result == expected.replace(microsecond=micros, nanosecond=0)
        elif unit == "s":
            assert result == expected.replace(microsecond=0, nanosecond=0)
        else:
            assert result == expected
        assert result._creso == getattr(NpyDatetimeUnit, f"NPY_FR_{unit}").value

    @pytest.mark.parametrize("year", ["2015", "2201"])
    @pytest.mark.parametrize("offset", [-1, 1])
    def test_timestamp_tz_localize_nonexistent_shift_invalid(
        self, offset, year, warsaw
    ):
        # GH 8917, 24466
        tz = warsaw
        _distant_date_only_for_zoneinfo(year, tz)
        ts = Timestamp(f"{year}-03-29 02:20:00")
        msg = "The provided timedelta will relocalize on a nonexistent time"
        with pytest.raises(ValueError, match=msg):
            ts.tz_localize(tz, nonexistent=timedelta(seconds=offset))

    def test_timestamp_tz_localize_nonexistent_shift_overflow(self):
        # GH#66697
        ts = Timestamp("2011-03-13 02:30").as_unit("ns")

        with pytest.raises(OutOfBoundsDatetime, match="overflows past"):
            ts.tz_localize("US/Eastern", nonexistent=Timedelta.max)

    @pytest.mark.parametrize("year", ["2015", "2201"])
    def test_timestamp_tz_localize_nonexistent_NaT(self, year, warsaw, unit):
        # GH 8917
        tz = warsaw
        _distant_date_only_for_zoneinfo(year, tz)
        ts = Timestamp(f"{year}-03-29 02:20:00").as_unit(unit)
        result = ts.tz_localize(tz, nonexistent="NaT")
        assert result is NaT

    @pytest.mark.parametrize("year", ["2015", "2201"])
    def test_timestamp_tz_localize_nonexistent_raise(self, year, warsaw, unit):
        # GH 8917
        tz = warsaw
        _distant_date_only_for_zoneinfo(year, tz)
        ts = Timestamp(f"{year}-03-29 02:20:00").as_unit(unit)
        msg = f"{year}-03-29 02:20:00"
        with pytest.raises(ValueError, match=msg):
            ts.tz_localize(tz, nonexistent="raise")
        msg = (
            "The nonexistent argument must be one of 'raise', 'NaT', "
            "'shift_forward', 'shift_backward' or a timedelta object"
        )
        with pytest.raises(ValueError, match=msg):
            ts.tz_localize(tz, nonexistent="foo")
