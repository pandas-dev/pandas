from datetime import (
    datetime,
    timedelta,
    timezone,
)
import subprocess
import sys
import textwrap
import zoneinfo

import dateutil.tz
import numpy as np
import pytest

from pandas._libs.tslibs import (
    conversion,
    timezones,
)
from pandas.compat import (
    IS64,
    is_platform_windows,
)

import pandas as pd
from pandas import (
    Index,
    Timedelta,
    Timestamp,
    date_range,
    to_datetime,
)
import pandas._testing as tm


@pytest.mark.single_cpu
def test_no_timezone_data():
    # https://github.com/pandas-dev/pandas/pull/63335
    # Test error message when timezone data is not available.
    msg = "'No time zone found with key Europe/Brussels'"
    code = textwrap.dedent(
        f"""\
        import sys, zoneinfo, pandas as pd
        sys.modules['tzdata'] = None
        zoneinfo.reset_tzpath(['/path/to/nowhere'])
        try:
            pd.to_datetime('2012-01-01').tz_localize('Europe/Brussels')
        except zoneinfo.ZoneInfoNotFoundError as err:
            assert str(err) == "{msg}"
        """
    )
    subprocess.check_call([sys.executable, "-c", code])


def test_is_utc(utc_fixture):
    tz = timezones.maybe_get_tz(utc_fixture)
    assert timezones.is_utc(tz)


def test_cache_keys_are_distinct_for_pytz_vs_dateutil():
    pytz = pytest.importorskip("pytz")
    for tz_name in pytz.common_timezones:
        tz_p = timezones.maybe_get_tz(tz_name)
        tz_d = timezones.maybe_get_tz("dateutil/" + tz_name)

    if tz_d is None:
        pytest.skip(tz_name + ": dateutil does not know about this one")

    if not (tz_name == "UTC" and is_platform_windows()):
        # they both end up as tzwin("UTC") on windows
        assert timezones._p_tz_cache_key(tz_p) != timezones._p_tz_cache_key(tz_d)


def test_tzlocal_repr():
    # see gh-13583
    ts = Timestamp("2011-01-01", tz=dateutil.tz.tzlocal())
    assert ts.tz == dateutil.tz.tzlocal()
    assert "tz='tzlocal()')" in repr(ts)


def test_tzlocal_maybe_get_tz():
    # see gh-13583
    tz = timezones.maybe_get_tz("tzlocal()")
    assert tz == dateutil.tz.tzlocal()


def test_tzlocal_offset():
    # see gh-13583
    #
    # Get offset using normal datetime for test.
    ts = Timestamp("2011-01-01", tz=dateutil.tz.tzlocal()).as_unit("s")

    offset = dateutil.tz.tzlocal().utcoffset(datetime(2011, 1, 1))
    offset = offset.total_seconds()

    assert ts._value + offset == Timestamp("2011-01-01").as_unit("s")._value


def test_tzlocal_is_not_utc():
    # even if the machine running the test is localized to UTC
    tz = dateutil.tz.tzlocal()
    assert not timezones.is_utc(tz)

    assert not timezones.tz_compare(tz, dateutil.tz.tzutc())


def test_tz_compare_utc(utc_fixture, utc_fixture2):
    tz = timezones.maybe_get_tz(utc_fixture)
    tz2 = timezones.maybe_get_tz(utc_fixture2)
    assert timezones.tz_compare(tz, tz2)


@pytest.fixture(
    params=[
        ("pytz/US/Eastern", lambda tz, x: tz.localize(x)),
        (dateutil.tz.gettz("US/Eastern"), lambda tz, x: x.replace(tzinfo=tz)),
    ]
)
def infer_setup(request):
    eastern, localize = request.param
    if isinstance(eastern, str) and eastern.startswith("pytz/"):
        pytz = pytest.importorskip("pytz")
        eastern = pytz.timezone(eastern.removeprefix("pytz/"))

    start_naive = datetime(2001, 1, 1)
    end_naive = datetime(2009, 1, 1)

    start = localize(eastern, start_naive)
    end = localize(eastern, end_naive)

    return eastern, localize, start, end, start_naive, end_naive


def test_infer_tz_compat(infer_setup):
    eastern, _, start, end, start_naive, end_naive = infer_setup

    assert (
        timezones.infer_tzinfo(start, end)
        is conversion.localize_pydatetime(start_naive, eastern).tzinfo
    )
    assert (
        timezones.infer_tzinfo(start, None)
        is conversion.localize_pydatetime(start_naive, eastern).tzinfo
    )
    assert (
        timezones.infer_tzinfo(None, end)
        is conversion.localize_pydatetime(end_naive, eastern).tzinfo
    )


def test_infer_tz_utc_localize(infer_setup):
    _, _, start, end, start_naive, end_naive = infer_setup
    utc = timezone.utc

    start = start_naive.astimezone(utc)
    end = end_naive.astimezone(utc)

    assert timezones.infer_tzinfo(start, end) is utc


@pytest.mark.parametrize("ordered", [True, False])
def test_infer_tz_mismatch(infer_setup, ordered):
    eastern, _, _, _, start_naive, end_naive = infer_setup
    msg = "Inputs must both have the same timezone"

    utc = timezone.utc
    start = start_naive.astimezone(utc)
    end = conversion.localize_pydatetime(end_naive, eastern)

    args = (start, end) if ordered else (end, start)

    with pytest.raises(AssertionError, match=msg):
        timezones.infer_tzinfo(*args)


def test_maybe_get_tz_invalid_types():
    with pytest.raises(TypeError, match="<class 'float'>"):
        timezones.maybe_get_tz(44.0)

    with pytest.raises(TypeError, match="<class 'module'>"):
        timezones.maybe_get_tz(pytest)

    msg = "<class 'pandas.Timestamp'>"
    with pytest.raises(TypeError, match=msg):
        timezones.maybe_get_tz(Timestamp("2021-01-01", tz="UTC"))


@pytest.mark.parametrize("tz_name", ["UTC", "GMT", "Etc/GMT+1", "Etc/GMT-5"])
def test_zoneinfo_fixed_offset(tz_name):
    # GH#64363
    zoneinfo = pytest.importorskip("zoneinfo")
    tz = zoneinfo.ZoneInfo(tz_name)
    assert timezones.is_fixed_offset(tz)


def test_zoneinfo_not_fixed_offset_with_historical_transition():
    # GH#64363
    zoneinfo = pytest.importorskip("zoneinfo")
    tz = zoneinfo.ZoneInfo("Africa/Lusaka")
    assert not timezones.is_fixed_offset(tz)


def test_maybe_get_tz_offset_only():
    # see gh-36004

    # timezone.utc
    tz = timezones.maybe_get_tz(timezone.utc)
    assert tz == timezone(timedelta(hours=0, minutes=0))

    # without UTC+- prefix
    tz = timezones.maybe_get_tz("+01:15")
    assert tz == timezone(timedelta(hours=1, minutes=15))

    tz = timezones.maybe_get_tz("-01:15")
    assert tz == timezone(-timedelta(hours=1, minutes=15))

    # with UTC+- prefix
    tz = timezones.maybe_get_tz("UTC+02:45")
    assert tz == timezone(timedelta(hours=2, minutes=45))

    tz = timezones.maybe_get_tz("UTC-02:45")
    assert tz == timezone(-timedelta(hours=2, minutes=45))


def test_zoneinfo_utc_to_local_post_2037():
    # GH#64363 - verify that ZoneInfo DST transitions after 2037
    # (generated from POSIX TZ string rules) produce correct local times.
    tz = zoneinfo.ZoneInfo("US/Pacific")
    utc_times = date_range("2040-07-01", periods=24, freq="h", tz="UTC")
    local = utc_times.tz_convert(tz)

    expected_hours = np.array(
        [
            datetime(2040, 7, 1, hour, tzinfo=timezone.utc).astimezone(tz).hour
            for hour in range(24)
        ],
        dtype=np.int32,
    )
    tm.assert_numpy_array_equal(local.hour.to_numpy(), expected_hours)


@pytest.mark.parametrize(
    "tz_name", ["America/New_York", "Australia/Melbourne", "Europe/Kaliningrad"]
)
def test_zoneinfo_utc_to_local_post_2100(tz_name):
    # GH#65712 - verify zoneinfo-generated future DST transitions are used
    # beyond year 2100 for both hemispheres.
    # Europe/Kaliningrad is a special case because it has no rule for future DST
    # transitions (and is also not a fixed offset timezone)
    tz = zoneinfo.ZoneInfo(tz_name)
    # fmt: off
    data = [
        "2038-01-01", "2038-07-01",
        "2100-01-01", "2100-07-01",
        "2200-01-01", "2200-07-01",
    ]
    # fmt: on
    utc_times = to_datetime(data, utc=True)
    local = utc_times.tz_convert(tz)

    expected = [
        datetime.fromisoformat(date).replace(tzinfo=timezone.utc).astimezone(tz)
        for date in data
    ]

    tm.assert_numpy_array_equal(local.to_pydatetime(), np.array(expected))
    assert [ts.utcoffset() for ts in local] == [dt.utcoffset() for dt in expected]


def test_zoneinfo_utc_to_local_far_future_seconds_resolution():
    # GH#65712 - for dates beyond cached transition data, we should fall back
    # to zoneinfo's API and preserve correct DST offsets.
    utc_times = to_datetime(
        np.array(["3000-01-01T00:00:00", "3000-07-01T00:00:00"], dtype="M8[s]"),
        utc=True,
    )
    local = utc_times.tz_convert("Europe/Brussels")

    tz = zoneinfo.ZoneInfo("Europe/Brussels")
    expected = [
        datetime(3000, 1, 1, tzinfo=timezone.utc).astimezone(tz),
        datetime(3000, 7, 1, tzinfo=timezone.utc).astimezone(tz),
    ]

    assert local.dtype == "datetime64[s, Europe/Brussels]"
    tm.assert_numpy_array_equal(local.to_pydatetime(), np.array(expected))


def test_zoneinfo_local_to_utc_far_future_seconds_resolution():
    # GH#65712 - localize should also preserve future DST rules when converting
    # from local wall times to UTC beyond cached transition data.
    local_times = to_datetime(
        np.array(["3000-01-01T00:00:00", "3000-07-01T00:00:00"], dtype="M8[s]")
    )
    localized = local_times.tz_localize("Europe/Brussels")
    result_utc = localized.tz_convert("UTC")

    tz = zoneinfo.ZoneInfo("Europe/Brussels")
    expected_utc = [
        datetime(3000, 1, 1, tzinfo=tz).astimezone(timezone.utc),
        datetime(3000, 7, 1, tzinfo=tz).astimezone(timezone.utc),
    ]

    assert result_utc.dtype == "datetime64[s, UTC]"
    tm.assert_numpy_array_equal(result_utc.to_pydatetime(), np.array(expected_utc))


@pytest.mark.skipif(
    not IS64,
    reason="stdlib datetime.fromtimestamp fails on 32-bit platforms with overflow",
)
@pytest.mark.parametrize(
    "tz_name",
    ["America/New_York", "America/Santiago", "Australia/Melbourne", "Europe/Brussels"],
)
def test_zoneinfo_boundary_at_last_cached_transition(tz_name):
    # GH#65733 - target the boundary where UTC is beyond the cached transition
    # range while local wall time can still compare below that UTC cutoff.
    tz = zoneinfo.ZoneInfo(tz_name)

    from zoneinfo._zoneinfo import ZoneInfo

    trans = datetime.fromtimestamp(max(ZoneInfo(tz_name)._tz_after.transitions(2099)))
    start = (trans - timedelta(days=1)).replace(tzinfo=timezone.utc)
    expected_utc = [(start + timedelta(minutes=30 * i)) for i in range(24 * 2 * 2)]
    expected_tz = [ts.astimezone(tz) for ts in expected_utc]
    expected_str = [str(ts) for ts in expected_tz]

    local = to_datetime([v[:-6] for v in expected_str])

    # ambiguous keyword only needed for timezones of northern hemisphere, ignored for
    # the others
    result = local.tz_localize(tz_name, ambiguous="infer")

    # verify local to UTC
    tm.assert_equal(
        result.tz_convert("UTC").to_pydatetime(), np.array(expected_utc, dtype="O")
    )
    # verify UTC to local by converting to pydatetime/str repr of underlying UTC value
    tm.assert_equal(result.astype(str), Index(expected_str))
    tm.assert_equal(result.to_pydatetime(), np.array(expected_tz, dtype="O"))


@pytest.mark.parametrize(
    "tz_name",
    ["America/New_York", "America/Santiago", "Australia/Melbourne", "Europe/Brussels"],
)
def test_zoneinfo_nonexistent_at_last_cached_transition(tz_name):
    # GH#65733 - target the boundary where UTC is beyond the cached transition
    # range while local wall time can still compare below that UTC cutoff.
    tz = zoneinfo.ZoneInfo(tz_name)

    from zoneinfo._zoneinfo import ZoneInfo

    start, end = ZoneInfo(tz_name)._tz_after.transitions(2099)
    if start > end:
        # to-dst nonexistent transition at the end of the year
        ts = Timestamp(start, unit="s") + Timedelta(minutes=30)

        with pytest.raises(ValueError, match="nonexistent time"):
            ts.tz_localize(tz_name)

        result = ts.tz_localize(tz_name, nonexistent="NaT")
        assert result is pd.NaT

        result_dst = ts.tz_localize(tz_name, nonexistent="shift_forward")
        expected_dst = (ts.to_pydatetime() + timedelta(minutes=30)).replace(tzinfo=tz)
        assert result_dst.to_pydatetime() == expected_dst

        result_std = ts.tz_localize(tz_name, nonexistent="shift_backward")
        expected_std = (
            ts.to_pydatetime() - timedelta(minutes=30, microseconds=1)
        ).replace(tzinfo=tz)
        assert result_std.to_pydatetime() == expected_std

    else:
        # to-std ambiguous transition at the end of the year
        ts = Timestamp(end, unit="s") - Timedelta(minutes=30)

        with pytest.raises(ValueError, match="Cannot infer dst time"):
            ts.tz_localize(tz_name)

        result = ts.tz_localize(tz_name, ambiguous="NaT")
        assert result is pd.NaT

        result_dst = ts.tz_localize(tz_name, ambiguous=True)
        expected_dst = ts.to_pydatetime().replace(fold=0, tzinfo=tz)
        assert result_dst.to_pydatetime() == expected_dst

        result_dst = ts.tz_localize(tz_name, ambiguous=False)
        expected_dst = ts.to_pydatetime().replace(fold=1, tzinfo=tz)
        assert result_dst.to_pydatetime() == expected_dst


@pytest.mark.parametrize("key", ["US/Eastern", "Africa/Lusaka", "Asia/Qyzylorda"])
def test_zoneinfo_utc_to_local_pre_first_transition(key):
    # GH#64363 - verify that ZoneInfo offsets before the first historical
    # transition (LMT era) match what ZoneInfo itself returns. The "before"
    # offset must be utcoff[0] from the TZ data (matching CPython's C
    # implementation), NOT the first non-DST transition offset.
    tz = zoneinfo.ZoneInfo(key)
    ts = Timestamp("1850-01-01", tz="UTC").tz_convert(tz)

    expected = datetime(1850, 1, 1, tzinfo=timezone.utc).astimezone(tz)
    assert ts.minute == expected.minute


def test_zoneinfo_conversion_outside_range_stdlib():
    # GH#65733 - verify that datetimes outside the range of Python's standard
    # library (year > 9999) raises a proper error message
    ts = Timestamp(np.datetime64("10000-01-01T09:00:00", "us"))

    msg = "Localizing Timestamps which are outside the range of Python"
    with pytest.raises(NotImplementedError, match=msg):
        ts.tz_localize("Europe/Brussels")

    with pytest.raises(NotImplementedError, match=msg):
        ts = Timestamp(ts._value, unit="us", tz="Europe/Brussels")


def test_normalize_pytz_timezone():
    pytz = pytest.importorskip("pytz")

    from pandas.io._util import _normalize_pytz_timezone

    for tz, expected in [
        (pytz.UTC, timezone.utc),
        (pytz.FixedOffset(90), timezone(timedelta(minutes=90))),
        (pytz.timezone("America/New_York"), zoneinfo.ZoneInfo("America/New_York")),
        (pytz.timezone("Etc/GMT+1"), zoneinfo.ZoneInfo("Etc/GMT+1")),
    ]:
        result = _normalize_pytz_timezone(tz)
        assert result == expected
