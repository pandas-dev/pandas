from datetime import (
    UTC,
    datetime,
    timedelta,
    timezone,
)
from io import BytesIO
import struct
from zoneinfo import ZoneInfo

import dateutil.tz
from dateutil.tz import gettz
import numpy as np
import pytest

from pandas.compat import WASM
from pandas.errors import (
    OutOfBoundsDatetime,
    Pandas4Warning,
)
import pandas.util._test_decorators as td

from pandas import (
    DatetimeIndex,
    Timedelta,
    Timestamp,
    bdate_range,
    date_range,
    offsets,
    to_datetime,
)
import pandas._testing as tm


@pytest.fixture(params=["pytz/US/Eastern", gettz("US/Eastern"), ZoneInfo("US/Eastern")])
def tz(request):
    if isinstance(request.param, str) and request.param.startswith("pytz/"):
        pytz = pytest.importorskip("pytz")
        return pytz.timezone(request.param.removeprefix("pytz/"))
    return request.param


class TestTZLocalize:
    def test_tz_localize_shift_onto_nat_sentinel(self, unit):
        # GH#66550 Asia/Tokyo's pre-1888 LMT offset is +9:18:59, so this wall
        #  time shifts to exactly the NaT sentinel.  _get_utc_bounds handed that
        #  value to bisect_right_i8, which requires val >= tdata[0] (one above
        #  the sentinel), so the following deltas lookup read out of bounds and
        #  the result was misreported as a nonexistent time.
        lmt = np.timedelta64(9 * 3600 + 18 * 60 + 59, "s")
        lmt_offset = int(lmt // np.timedelta64(1, unit))
        dti = DatetimeIndex([np.datetime64(-(2**63) + lmt_offset, unit)])
        with pytest.raises(OutOfBoundsDatetime, match="underflows past"):
            dti.tz_localize("Asia/Tokyo")

    @pytest.mark.parametrize(
        "kwargs", [{}, {"nonexistent": "NaT"}, {"nonexistent": "shift_forward"}]
    )
    def test_tz_localize_overflow_past_last_cached_transition(self, kwargs):
        # GH#65733 the zoneinfo fallback used past the last cached DST transition
        #  wrapped instead of reporting an out-of-bounds UTC instant, so an
        #  ordinary wall time came back as NaT
        dti = DatetimeIndex([Timestamp(Timestamp.max._value - 2 * 3600 * 10**9)])

        with pytest.raises(OutOfBoundsDatetime, match="overflows past"):
            dti.tz_localize("America/New_York", **kwargs)

    @pytest.mark.parametrize(
        "tz", [timezone(timedelta(hours=9)), ZoneInfo("Etc/GMT-9")]
    )
    def test_tz_localize_fixed_offset_shift_onto_nat_sentinel(self, tz):
        # GH#66550 a shift landing exactly on the NaT sentinel is one below the
        #  minimum representable value, so it underflows.  The fixed-offset loop
        #  only checked for int64 wrap, so it stored the sentinel and the wall time
        #  read back as missing data.
        dti = DatetimeIndex(
            np.array([-(2**63) + 9 * 3600 * 10**9], dtype="i8").view("M8[ns]")
        )

        with pytest.raises(OutOfBoundsDatetime, match="underflows past"):
            dti.tz_localize(tz)

    @pytest.mark.parametrize(
        "tz", [timezone(timedelta(hours=-5)), ZoneInfo("Etc/GMT+5")]
    )
    def test_tz_localize_overflow_fixed_offset(self, tz):
        # GH#65733 the fixed-offset loop in tz_localize_to_utc shifted to UTC with a
        #  bare subtraction, so a wall time whose UTC instant is past Timestamp.max
        #  wrapped to a negative year instead of raising.  Needs >= 2 values: for a
        #  length-1 result DatetimeIndex.tz_localize boxes the lone value, and that
        #  incidentally catches the overflow downstream.
        val = Timestamp.max._value - 2 * 3600 * 10**9
        dti = DatetimeIndex(np.array([val, val - 10**9], dtype="i8").view("M8[ns]"))

        with pytest.raises(OutOfBoundsDatetime, match="overflows past"):
            dti.tz_localize(tz)

    @td.skip_if_windows  # `tm.set_timezone` does not work on Windows
    @pytest.mark.skipif(WASM, reason="`tm.set_timezone` does not work on WASM")
    def test_tz_localize_tzlocal_shift_onto_nat_sentinel(self):
        # GH#66550 the tzlocal loop needs the same sentinel guard as the
        #  fixed-offset one above.  Pin the ambient zone to +9 so the shift lands
        #  exactly on the sentinel.
        dti = DatetimeIndex(
            np.array([-(2**63) + 9 * 3600 * 10**9], dtype="i8").view("M8[ns]")
        )

        with tm.set_timezone("Asia/Tokyo"):
            with pytest.raises(OutOfBoundsDatetime, match="underflows past"):
                dti.tz_localize(dateutil.tz.tzlocal())

    @td.skip_if_windows  # `tm.set_timezone` does not work on Windows
    @pytest.mark.skipif(WASM, reason="`tm.set_timezone` does not work on WASM")
    @td.skip_if_32bit  # OverflowError inside tzlocal past 2038
    def test_tz_localize_overflow_tzlocal(self):
        # GH#65733 same bare-subtraction overflow in the tzlocal loop.  Pin the
        #  ambient zone: the shift only overflows west of UTC.
        val = Timestamp.max._value - 2 * 3600 * 10**9
        dti = DatetimeIndex(np.array([val, val - 10**9], dtype="i8").view("M8[ns]"))

        with tm.set_timezone("US/Eastern"):
            with pytest.raises(OutOfBoundsDatetime, match="overflows past"):
                dti.tz_localize(dateutil.tz.tzlocal())

    def test_tz_localize_invalidates_freq(self):
        # we only preserve freq in unambiguous cases

        # if localized to US/Eastern, this crosses a DST transition
        dti = date_range("2014-03-08 23:00", "2014-03-09 09:00", freq="h")
        assert dti.freq == "h"

        result = dti.tz_localize(None)  # no-op
        assert result.freq == "h"

        result = dti.tz_localize("UTC")  # unambiguous freq preservation
        assert result.freq == "h"

        result = dti.tz_localize("US/Eastern", nonexistent="shift_forward")
        assert result.freq is None
        assert result.inferred_freq is None  # i.e. we are not _too_ strict here

        # Case where we _can_ keep freq because we're length==1
        dti2 = dti[:1]
        result = dti2.tz_localize("US/Eastern")
        assert result.freq == "h"

    def test_tz_localize_utc_copies(self, utc_fixture):
        # GH#46460
        times = ["2015-03-08 01:00", "2015-03-08 02:00", "2015-03-08 03:00"]
        index = DatetimeIndex(times)

        res = index.tz_localize(utc_fixture)
        assert not tm.shares_memory(res, index)

        res2 = index._data.tz_localize(utc_fixture)
        assert not tm.shares_memory(index._data, res2)

    def test_dti_tz_localize_nonexistent_raise_coerce(self):
        # GH#13057
        times = ["2015-03-08 01:00", "2015-03-08 02:00", "2015-03-08 03:00"]
        index = DatetimeIndex(times)
        tz = "US/Eastern"
        with pytest.raises(ValueError, match="|".join(times)):
            index.tz_localize(tz=tz)

        with pytest.raises(ValueError, match="|".join(times)):
            index.tz_localize(tz=tz, nonexistent="raise")

        result = index.tz_localize(tz=tz, nonexistent="NaT")
        test_times = ["2015-03-08 01:00-05:00", "NaT", "2015-03-08 03:00-04:00"]
        dti = to_datetime(test_times, utc=True)
        expected = dti.tz_convert("US/Eastern")
        tm.assert_index_equal(result, expected)

    def test_dti_tz_localize_ambiguous_infer(self, tz):
        # November 6, 2011, fall back, repeat 2 AM hour
        # With no repeated hours, we cannot infer the transition
        dr = date_range(datetime(2011, 11, 6, 0), periods=5, freq=offsets.Hour())
        with pytest.raises(ValueError, match="Cannot infer dst time"):
            dr.tz_localize(tz)

    def test_dti_tz_localize_ambiguous_infer2(self, tz, unit):
        # With repeated hours, we can infer the transition
        dr = date_range(
            datetime(2011, 11, 6, 0), periods=5, freq=offsets.Hour(), tz=tz, unit=unit
        )
        times = [
            "11/06/2011 00:00",
            "11/06/2011 01:00",
            "11/06/2011 01:00",
            "11/06/2011 02:00",
            "11/06/2011 03:00",
        ]
        di = DatetimeIndex(times).as_unit(unit)
        result = di.tz_localize(tz, ambiguous="infer")
        expected = dr._with_freq(None)
        tm.assert_index_equal(result, expected)
        depr_msg = "The 'ambiguous' keyword in DatetimeIndex is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
            result2 = DatetimeIndex(times, tz=tz, ambiguous="infer").as_unit(unit)
        tm.assert_index_equal(result2, expected)

    def test_dti_tz_localize_ambiguous_infer3(self, tz):
        # When there is no dst transition, nothing special happens
        dr = date_range(datetime(2011, 6, 1, 0), periods=10, freq=offsets.Hour())
        localized = dr.tz_localize(tz)
        localized_infer = dr.tz_localize(tz, ambiguous="infer")
        tm.assert_index_equal(localized, localized_infer)

    def test_dti_tz_localize_ambiguous_times(self, tz):
        # March 13, 2011, spring forward, skip from 2 AM to 3 AM
        dr = date_range(datetime(2011, 3, 13, 1, 30), periods=3, freq=offsets.Hour())
        with pytest.raises(ValueError, match="2011-03-13 02:30:00"):
            dr.tz_localize(tz)

        # after dst transition, it works
        dr = date_range(
            datetime(2011, 3, 13, 3, 30), periods=3, freq=offsets.Hour(), tz=tz
        )

        # November 6, 2011, fall back, repeat 2 AM hour
        dr = date_range(datetime(2011, 11, 6, 1, 30), periods=3, freq=offsets.Hour())
        with pytest.raises(ValueError, match="Cannot infer dst time"):
            dr.tz_localize(tz)

        # UTC is OK
        dr = date_range(
            datetime(2011, 3, 13), periods=48, freq=offsets.Minute(30), tz=UTC
        )

    @pytest.mark.parametrize("tzstr", ["US/Eastern", "dateutil/US/Eastern"])
    def test_dti_tz_localize_pass_dates_to_utc(self, tzstr):
        strdates = ["1/1/2012", "3/1/2012", "4/1/2012"]

        idx = DatetimeIndex(strdates)
        conv = idx.tz_localize(tzstr)

        fromdates = DatetimeIndex(strdates, tz=tzstr)

        assert conv.tz == fromdates.tz
        tm.assert_numpy_array_equal(conv.asi8, fromdates.asi8)

    @pytest.mark.parametrize("prefix", ["", "dateutil/"])
    def test_dti_tz_localize(self, prefix):
        tzstr = prefix + "US/Eastern"
        dti = date_range(start="1/1/2005", end="1/1/2005 0:00:02.256", freq="ms")
        dti2 = dti.tz_localize(tzstr)

        dti_utc = date_range(
            start="1/1/2005 05:00", end="1/1/2005 5:00:02.256", freq="ms", tz="utc"
        )

        tm.assert_numpy_array_equal(dti2.to_numpy(), dti_utc.to_numpy())

        dti3 = dti2.tz_convert(prefix + "US/Pacific")
        tm.assert_numpy_array_equal(dti3.to_numpy(), dti_utc.to_numpy())

        dti = date_range(start="11/6/2011 1:59:59", end="11/6/2011 2:00", freq="ms")
        with pytest.raises(ValueError, match="Cannot infer dst time"):
            dti.tz_localize(tzstr)

        dti = date_range(start="3/13/2011 1:59:59", end="3/13/2011 2:00", freq="ms")
        with pytest.raises(ValueError, match="2011-03-13 02:00:00"):
            dti.tz_localize(tzstr)

    def test_dti_tz_localize_utc_conversion(self, tz):
        # Localizing to time zone should:
        #  1) check for DST ambiguities
        #  2) convert to UTC

        rng = date_range("3/10/2012", "3/11/2012", freq="30min")

        converted = rng.tz_localize(tz)
        expected_naive = rng + offsets.Hour(5)
        tm.assert_numpy_array_equal(converted.asi8, expected_naive.asi8)

        # DST ambiguity, this should fail
        rng = date_range("3/11/2012", "3/12/2012", freq="30min")
        # Is this really how it should fail??
        with pytest.raises(ValueError, match="2012-03-11 02:00:00"):
            rng.tz_localize(tz)

    def test_dti_tz_localize_roundtrip(self, tz_aware_fixture):
        # note: this tz tests that a tz-naive index can be localized
        # and de-localized successfully, when there are no DST transitions
        # in the range.
        idx = date_range(start="2014-06-01", end="2014-08-30", freq="15min")
        tz = tz_aware_fixture
        localized = idx.tz_localize(tz)
        # can't localize a tz-aware object
        with pytest.raises(
            TypeError, match="Already tz-aware, use tz_convert to convert"
        ):
            localized.tz_localize(tz)
        reset = localized.tz_localize(None)
        assert reset.tzinfo is None
        expected = idx._with_freq(None)
        tm.assert_index_equal(reset, expected)

    def test_dti_tz_localize_naive(self):
        rng = date_range("1/1/2011", periods=100, freq="h")

        conv = rng.tz_localize("US/Pacific")
        exp = date_range("1/1/2011", periods=100, freq="h", tz="US/Pacific")

        tm.assert_index_equal(conv, exp._with_freq(None))

    def test_dti_tz_localize_tzlocal(self):
        # GH#13583
        offset = dateutil.tz.tzlocal().utcoffset(datetime(2011, 1, 1))
        offset = int(offset.total_seconds() * 1000000000)

        dti = date_range(start="2001-01-01", end="2001-03-01", unit="ns")
        dti2 = dti.tz_localize(dateutil.tz.tzlocal())
        tm.assert_numpy_array_equal(dti2.asi8 + offset, dti.asi8)

        dti = date_range(
            start="2001-01-01", end="2001-03-01", tz=dateutil.tz.tzlocal(), unit="ns"
        )
        dti2 = dti.tz_localize(None)
        tm.assert_numpy_array_equal(dti2.asi8 - offset, dti.asi8)

    def test_dti_tz_localize_ambiguous_nat(self, tz):
        times = [
            "11/06/2011 00:00",
            "11/06/2011 01:00",
            "11/06/2011 01:00",
            "11/06/2011 02:00",
            "11/06/2011 03:00",
        ]
        di = DatetimeIndex(times)
        localized = di.tz_localize(tz, ambiguous="NaT")

        times = [
            "11/06/2011 00:00",
            np.nan,
            np.nan,
            "11/06/2011 02:00",
            "11/06/2011 03:00",
        ]
        di_test = DatetimeIndex(times, tz="US/Eastern")

        # left dtype is datetime64[ns, US/Eastern]
        # right is datetime64[ns, tzfile('/usr/share/zoneinfo/US/Eastern')]
        tm.assert_numpy_array_equal(di_test.asi8, localized.asi8)

    def test_dti_tz_localize_ambiguous_flags(self, tz, unit):
        # November 6, 2011, fall back, repeat 2 AM hour

        # Pass in flags to determine right dst transition
        dr = date_range(
            datetime(2011, 11, 6, 0), periods=5, freq=offsets.Hour(), tz=tz, unit=unit
        )
        times = [
            "11/06/2011 00:00",
            "11/06/2011 01:00",
            "11/06/2011 01:00",
            "11/06/2011 02:00",
            "11/06/2011 03:00",
        ]

        # Test tz_localize
        di = DatetimeIndex(times).as_unit(unit)
        is_dst = [1, 1, 0, 0, 0]
        localized = di.tz_localize(tz, ambiguous=is_dst)
        expected = dr._with_freq(None)
        tm.assert_index_equal(expected, localized, check_freq=False)

        depr_msg = "The 'ambiguous' keyword in DatetimeIndex is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
            result = DatetimeIndex(times, tz=tz, ambiguous=is_dst).as_unit(unit)
        tm.assert_index_equal(result, expected)

        localized = di.tz_localize(tz, ambiguous=np.array(is_dst))
        tm.assert_index_equal(dr, localized, check_freq=False)

        localized = di.tz_localize(tz, ambiguous=np.array(is_dst).astype("bool"))
        tm.assert_index_equal(dr, localized, check_freq=False)

        # Test constructor
        with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
            localized = DatetimeIndex(times, tz=tz, ambiguous=is_dst).as_unit(unit)
        tm.assert_index_equal(dr, localized, check_freq=False)

        # Test duplicate times where inferring the dst fails
        times += times
        di = DatetimeIndex(times).as_unit(unit)

        # When the sizes are incompatible, make sure error is raised
        msg = "Length of ambiguous bool-array must be the same size as vals"
        with pytest.raises(Exception, match=msg):
            di.tz_localize(tz, ambiguous=is_dst)

        # When sizes are compatible and there are repeats ('infer' won't work)
        is_dst = np.hstack((is_dst, is_dst))
        localized = di.tz_localize(tz, ambiguous=is_dst)
        dr = dr.append(dr)
        tm.assert_index_equal(dr, localized)

    def test_dti_tz_localize_ambiguous_flags2(self, tz):
        # When there is no dst transition, nothing special happens
        dr = date_range(datetime(2011, 6, 1, 0), periods=10, freq=offsets.Hour())
        is_dst = np.array([1] * 10)
        localized = dr.tz_localize(tz)
        localized_is_dst = dr.tz_localize(tz, ambiguous=is_dst)
        tm.assert_index_equal(localized, localized_is_dst)

    def test_dti_tz_localize_bdate_range(self):
        dr = bdate_range("1/1/2009", "1/1/2010")
        dr_utc = bdate_range("1/1/2009", "1/1/2010", tz=UTC)
        localized = dr.tz_localize(UTC)
        tm.assert_index_equal(dr_utc, localized)

    @pytest.mark.parametrize(
        "start_ts, tz, end_ts, shift",
        [
            ["2015-03-29 02:20:00", "Europe/Warsaw", "2015-03-29 03:00:00", "forward"],
            [
                "2015-03-29 02:20:00",
                "Europe/Warsaw",
                "2015-03-29 01:59:59.999999999",
                "backward",
            ],
            [
                "2015-03-29 02:20:00",
                "Europe/Warsaw",
                "2015-03-29 03:20:00",
                timedelta(hours=1),
            ],
            [
                "2015-03-29 02:20:00",
                "Europe/Warsaw",
                "2015-03-29 01:20:00",
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
            # GH#40705, GH#40915: shifting across the UTC+0 DST boundary
            ["2021-03-28 01:20:00", "Europe/London", "2021-03-28 02:00:00", "forward"],
            [
                "2021-03-28 01:20:00",
                "Europe/London",
                "2021-03-28 00:59:59.999999999",
                "backward",
            ],
            [
                "2021-03-28 01:20:00",
                "Europe/London",
                "2021-03-28 02:20:00",
                timedelta(hours=1),
            ],
            [
                "2021-03-28 01:20:00",
                "Europe/London",
                "2021-03-28 00:20:00",
                timedelta(hours=-1),
            ],
        ],
    )
    @pytest.mark.parametrize("tz_type", ["", "dateutil/"])
    def test_dti_tz_localize_nonexistent_shift(
        self, start_ts, tz, end_ts, shift, tz_type, unit
    ):
        # GH#8917
        tz = tz_type + tz
        if isinstance(shift, str):
            shift = "shift_" + shift
        dti = DatetimeIndex([Timestamp(start_ts)]).as_unit(unit)
        result = dti.tz_localize(tz, nonexistent=shift)
        expected = DatetimeIndex([Timestamp(end_ts)]).tz_localize(tz).as_unit(unit)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("offset", [-1, 1])
    def test_dti_tz_localize_nonexistent_shift_invalid(self, offset, warsaw):
        # GH#8917
        tz = warsaw
        dti = DatetimeIndex([Timestamp("2015-03-29 02:20:00")])
        msg = "The provided timedelta will relocalize on a nonexistent time"
        with pytest.raises(ValueError, match=msg):
            dti.tz_localize(tz, nonexistent=timedelta(seconds=offset))

    def test_dti_tz_localize_nonexistent_shift_overflow(self):
        # GH#66697
        dti = DatetimeIndex(["2011-03-13 02:30"]).as_unit("ns")

        with pytest.raises(OutOfBoundsDatetime, match="overflows past"):
            dti.tz_localize("US/Eastern", nonexistent=Timedelta.max)

    @pytest.mark.parametrize(
        "tz, start_ts",
        [
            # US/Eastern has a DST rule past the last cached transition and
            #  America/Sao_Paulo does not, so the two reach the offset that
            #  applies by different routes
            ("US/Eastern", "2011-03-13 02:30"),
            ("America/Sao_Paulo", "2018-11-04 00:30"),
        ],
    )
    def test_dti_tz_localize_nonexistent_shift_utc_overflow(self, tz, start_ts):
        # GH#66697
        ts = Timestamp(start_ts).as_unit("ns")
        dti = DatetimeIndex([ts, ts])
        shift = Timestamp.max - ts - Timedelta(hours=1)

        with pytest.raises(OutOfBoundsDatetime, match="overflows past"):
            dti.tz_localize(tz, nonexistent=shift)


def test_dti_tz_localize_nonexistent_shift_past_last_transition():
    # GH#66550 a shift landing past the last cached DST transition indexed the
    #  deltas array out of bounds, so the same call gave a different answer
    #  each time
    tz = gettz("America/Boise")
    dti = DatetimeIndex([Timestamp("2015-03-08 02:30")])
    shift = Timedelta(days=30_000)

    result = dti.tz_localize(tz, nonexistent=shift)

    expected = DatetimeIndex([dti[0] + shift]).tz_localize(tz)
    tm.assert_index_equal(result, expected)
    tm.assert_index_equal(result, dti.tz_localize(tz, nonexistent=shift))


@pytest.mark.parametrize(
    "tz, start_ts, shift_hours",
    [
        # east of UTC, shifting backward over the transition
        ("Europe/Warsaw", "2011-03-27 02:30", -2),
        ("Europe/Warsaw", "2011-03-27 02:30", -24),
        ("Europe/Berlin", "2011-03-27 02:30", -5),
        ("Asia/Tehran", "2011-03-22 00:30", -5),
        ("Pacific/Auckland", "2011-09-25 02:30", -3),
        ("Australia/Sydney", "2011-10-02 02:30", -3),
        # west of UTC, shifting forward over it
        ("US/Eastern", "2011-03-13 02:30", 5),
        ("US/Eastern", "2011-03-13 02:30", 24),
        ("America/Santiago", "2011-08-21 00:30", 5),
        # GH#40705 fixed this zone for a one-hour shift only
        ("Europe/London", "2021-03-28 01:20", 2),
    ],
)
def test_dti_tz_localize_nonexistent_timedelta_shift_across_transition(
    tz, start_ts, shift_hours, unit
):
    # GH#66820 the offset to apply was picked by bisecting the shifted *wall*
    #  time against the *UTC* transition instants, so any shift large enough to
    #  clear the transition took the offset from the wrong side of it and the
    #  result came back off by the DST delta.  One-hour shifts happened to land
    #  right, which is all the tests covered.
    # two elements so the vectorized loop is what is under test: at length 1
    #  tz_localize boxes element 0 to decide whether freq survives, which can
    #  raise on its own and mask what tz_localize_to_utc returned
    shift = Timedelta(hours=shift_hours)
    dti = DatetimeIndex([Timestamp(start_ts)] * 2).as_unit(unit)

    result = dti.tz_localize(tz, nonexistent=shift)

    expected = (
        DatetimeIndex([Timestamp(start_ts) + shift] * 2).tz_localize(tz).as_unit(unit)
    )
    tm.assert_index_equal(result, expected)


def test_tz_localize_nonexistent_timedelta_shift_scalar_matches_index():
    # GH#66820 the scalar path routes through the same block
    ts = Timestamp("2011-03-27 02:30")
    shift = Timedelta(hours=-2)

    result = ts.tz_localize("Europe/Warsaw", nonexistent=shift)

    assert result == Timestamp("2011-03-27 00:30", tz="Europe/Warsaw")
    dti = DatetimeIndex([ts] * 2).tz_localize("Europe/Warsaw", nonexistent=shift)
    assert (dti == result).all()


def test_dti_tz_localize_nonexistent_timedelta_shift_onto_nat_sentinel():
    # GH#66697 shifting a nonexistent time to a wall time whose UTC instant is
    #  exactly the NaT sentinel, one below Timestamp.min, used to come back as a
    #  missing value.  Asia/Tokyo's +9 offset is what makes the two line up.
    dti = DatetimeIndex([Timestamp("1948-05-02 00:30")]).as_unit("ns")
    jst = 9 * 3600 * 10**9
    shift = Timedelta(Timestamp.min._value - 1 + jst - dti.asi8[0], "ns")

    with pytest.raises(OutOfBoundsDatetime, match="underflows past"):
        dti.tz_localize("Asia/Tokyo", nonexistent=shift)


def _make_tzfile_ending_in_spring_forward(filename):
    """
    Build a UTC+11/UTC+12 zone whose final transition is a spring-forward.

    Real zones cannot be used here: a zone that ends this way exists in the IANA
    database (e.g. Asia/Anadyr), but most builds of the tz files append a
    sentinel transition at the 32-bit epoch rollover, which reuses the offset
    already in effect and so hides the case entirely.
    """
    abbrevs = b"+11\x00+12\x00"
    # (utc offset in seconds, isdst, index into abbrevs)
    ttinfos = [(39600, 0, 0), (43200, 0, 4)]
    transitions = [
        (int(datetime(2010, 10, 30, 15, tzinfo=UTC).timestamp()), 0),
        (int(datetime(2011, 3, 26, 15, tzinfo=UTC).timestamp()), 1),
    ]

    data = b"TZif" + b"\x00" + b"\x00" * 15
    data += struct.pack(">6l", 0, 0, 0, len(transitions), len(ttinfos), len(abbrevs))
    data += b"".join(struct.pack(">l", trans) for trans, _ in transitions)
    data += bytes(idx for _, idx in transitions)
    data += b"".join(struct.pack(">lBB", *ttinfo) for ttinfo in ttinfos)
    data += abbrevs
    return dateutil.tz.tzfile(BytesIO(data), filename=filename)


def test_dti_tz_localize_nonexistent_shift_at_last_transition():
    # GH#66550 when the last transition is a spring-forward, the delta in effect
    #  after it puts the shifted wall time back before that transition; the
    #  offset from before it is the one that applies
    tz = _make_tzfile_ending_in_spring_forward("/pandas-test/LastTransSpringForward")
    dti = DatetimeIndex([Timestamp("2011-03-27 02:30")])

    result = dti.tz_localize(tz, nonexistent="shift_backward")

    expected = DatetimeIndex([Timestamp("2011-03-27 01:59:59.999999")]).tz_localize(tz)
    tm.assert_index_equal(result, expected)
    assert result[0].utcoffset() == timedelta(hours=11)
