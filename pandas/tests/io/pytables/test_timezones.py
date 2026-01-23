from datetime import (
    date,
    timedelta,
)

import numpy as np
import pytest

from pandas._libs.tslibs.timezones import maybe_get_tz
import pandas.util._test_decorators as td

import pandas as pd
from pandas import (
    DataFrame,
    DatetimeIndex,
    Series,
    Timestamp,
    date_range,
)
import pandas._testing as tm


def _compare_with_tz(a, b):
    tm.assert_frame_equal(a, b)

    # compare the zones on each element
    for c in a.columns:
        for i in a.index:
            a_e = a.loc[i, c]
            b_e = b.loc[i, c]
            if not (a_e == b_e and a_e.tz == b_e.tz):
                raise AssertionError(f"invalid tz comparison [{a_e}] [{b_e}]")


# use maybe_get_tz instead of dateutil.tz.gettz to handle the windows
# filename issues.
gettz_dateutil = lambda x: maybe_get_tz("dateutil/" + x)
gettz_pytz = lambda x: x


@pytest.mark.filterwarnings(
    "ignore:`alltrue` is deprecated as of NumPy 1.25.0:DeprecationWarning"
)
@pytest.mark.parametrize("gettz", [gettz_dateutil, gettz_pytz])
def test_append_with_timezones(temp_hdfstore, gettz):
    # as columns

    # Single-tzinfo, no DST transition
    df_est = DataFrame(
        {
            "A": [
                Timestamp("20130102 2:00:00", tz=gettz("US/Eastern")).as_unit("ns")
                + timedelta(hours=1) * i
                for i in range(5)
            ]
        }
    )

    # frame with all columns having same tzinfo, but different sides
    #  of DST transition
    df_crosses_dst = DataFrame(
        {
            "A": Timestamp("20130102", tz=gettz("US/Eastern")).as_unit("ns"),
            "B": Timestamp("20130603", tz=gettz("US/Eastern")).as_unit("ns"),
        },
        index=range(5),
    )

    df_mixed_tz = DataFrame(
        {
            "A": Timestamp("20130102", tz=gettz("US/Eastern")).as_unit("ns"),
            "B": Timestamp("20130102", tz=gettz("EET")).as_unit("ns"),
        },
        index=range(5),
    )

    df_different_tz = DataFrame(
        {
            "A": Timestamp("20130102", tz=gettz("US/Eastern")).as_unit("ns"),
            "B": Timestamp("20130102", tz=gettz("CET")).as_unit("ns"),
        },
        index=range(5),
    )

    temp_hdfstore.append("df_tz", df_est, data_columns=["A"])
    result = temp_hdfstore["df_tz"]
    _compare_with_tz(result, df_est)
    tm.assert_frame_equal(result, df_est)

    # select with tz aware
    expected = df_est[df_est.A >= df_est.A[3]]
    result = temp_hdfstore.select("df_tz", where="A>=df_est.A[3]")
    _compare_with_tz(result, expected)

    # ensure we include dates in DST and STD time here.
    temp_hdfstore.remove("df_tz")
    temp_hdfstore.append("df_tz", df_crosses_dst)
    result = temp_hdfstore["df_tz"]
    _compare_with_tz(result, df_crosses_dst)
    tm.assert_frame_equal(result, df_crosses_dst)

    msg = (
        r"invalid info for \[values_block_1\] for \[tz\], "
        r"existing_value \[(dateutil/.*)?(US/Eastern|America/New_York)\] "
        r"conflicts with new value \[(dateutil/.*)?EET\]"
    )
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.append("df_tz", df_mixed_tz)

    # this is ok
    temp_hdfstore.remove("df_tz")
    temp_hdfstore.append("df_tz", df_mixed_tz, data_columns=["A", "B"])
    result = temp_hdfstore["df_tz"]
    _compare_with_tz(result, df_mixed_tz)
    tm.assert_frame_equal(result, df_mixed_tz)

    # can't append with diff timezone
    msg = (
        r"invalid info for \[B\] for \[tz\], "
        r"existing_value \[(dateutil/.*)?EET\] "
        r"conflicts with new value \[(dateutil/.*)?CET\]"
    )
    with pytest.raises(ValueError, match=msg):
        temp_hdfstore.append("df_tz", df_different_tz)


@pytest.mark.parametrize("gettz", [gettz_dateutil, gettz_pytz])
def test_append_with_timezones_as_index(temp_hdfstore, gettz):
    # GH#4098 example

    dti = date_range("2000-1-1", periods=3, freq="h", tz=gettz("US/Eastern"))
    dti = dti._with_freq(None)  # freq doesn't round-trip

    df = DataFrame({"A": Series(range(3), index=dti)})

    temp_hdfstore.put("df", df)
    result = temp_hdfstore.select("df")
    tm.assert_frame_equal(result, df)

    temp_hdfstore.remove("df")
    temp_hdfstore.append("df", df)
    result = temp_hdfstore.select("df")
    tm.assert_frame_equal(result, df)


def test_roundtrip_tz_aware_index(temp_hdfstore, unit):
    # GH 17618
    ts = Timestamp("2000-01-01 01:00:00", tz="US/Eastern")
    dti = DatetimeIndex([ts]).as_unit(unit)
    df = DataFrame(data=[0], index=dti)

    temp_hdfstore.put("frame", df, format="fixed")
    recons = temp_hdfstore["frame"]
    tm.assert_frame_equal(recons, df)

    value = recons.index[0]._value
    denom = {"ns": 1, "us": 1000, "ms": 10**6, "s": 10**9}[unit]
    assert value == 946706400000000000 // denom


def test_store_index_name_with_tz(temp_hdfstore):
    # GH 13884
    df = DataFrame({"A": [1, 2]})
    df.index = DatetimeIndex([1234567890123456787, 1234567890123456788])
    df.index = df.index.tz_localize("UTC")
    df.index.name = "foo"

    temp_hdfstore.put("frame", df, format="table")
    recons = temp_hdfstore["frame"]
    tm.assert_frame_equal(recons, df)


def test_tseries_select_index_column(temp_hdfstore):
    # GH7777
    # selecting a UTC datetimeindex column did
    # not preserve UTC tzinfo set before storing

    # check that no tz still works
    rng = date_range("1/1/2000", "1/30/2000")
    frame = DataFrame(
        np.random.default_rng(2).standard_normal((len(rng), 4)), index=rng
    )

    temp_hdfstore.append("frame", frame)
    result = temp_hdfstore.select_column("frame", "index")
    assert rng.tz == DatetimeIndex(result.values).tz

    # check utc
    rng = date_range("1/1/2000", "1/30/2000", tz="UTC")
    frame = DataFrame(
        np.random.default_rng(2).standard_normal((len(rng), 4)), index=rng
    )

    temp_hdfstore.remove("frame")
    temp_hdfstore.append("frame", frame)
    result = temp_hdfstore.select_column("frame", "index")
    assert rng.tz == result.dt.tz

    # double check non-utc
    rng = date_range("1/1/2000", "1/30/2000", tz="US/Eastern")
    frame = DataFrame(
        np.random.default_rng(2).standard_normal((len(rng), 4)), index=rng
    )

    temp_hdfstore.remove("frame")
    temp_hdfstore.append("frame", frame)
    result = temp_hdfstore.select_column("frame", "index")
    assert rng.tz == result.dt.tz


def test_timezones_fixed_format_frame_non_empty(temp_hdfstore):
    # index
    rng = date_range("1/1/2000", "1/30/2000", tz="US/Eastern")
    rng = rng._with_freq(None)  # freq doesn't round-trip
    df = DataFrame(np.random.default_rng(2).standard_normal((len(rng), 4)), index=rng)
    temp_hdfstore["df"] = df
    result = temp_hdfstore["df"]
    tm.assert_frame_equal(result, df)


def test_timezones_fixed_format_frame_non_empty_as_data(temp_hdfstore):
    # GH11411
    rng = date_range("1/1/2000", "1/30/2000", tz="US/Eastern")
    rng = rng._with_freq(None)  # freq doesn't round-trip
    df = DataFrame(
        {
            "A": rng,
            "B": rng.tz_convert("UTC").tz_localize(None),
            "C": rng.tz_convert("CET"),
            "D": range(len(rng)),
        },
        index=rng,
    )
    temp_hdfstore["df"] = df
    result = temp_hdfstore["df"]
    tm.assert_frame_equal(result, df)


def test_timezones_fixed_format_empty(temp_hdfstore, tz_aware_fixture, frame_or_series):
    # GH 20594

    dtype = pd.DatetimeTZDtype(tz=tz_aware_fixture)

    obj = Series(dtype=dtype, name="A")
    if frame_or_series is DataFrame:
        obj = obj.to_frame()

    temp_hdfstore["obj"] = obj
    result = temp_hdfstore["obj"]
    tm.assert_equal(result, obj)


def test_timezones_fixed_format_series_nonempty(temp_hdfstore, tz_aware_fixture):
    # GH 20594

    dtype = pd.DatetimeTZDtype(tz=tz_aware_fixture)

    s = Series([0], dtype=dtype)
    temp_hdfstore["s"] = s
    result = temp_hdfstore["s"]
    tm.assert_series_equal(result, s)


def test_fixed_offset_tz(temp_hdfstore):
    rng = date_range("1/1/2000 00:00:00-07:00", "1/30/2000 00:00:00-07:00")
    frame = DataFrame(
        np.random.default_rng(2).standard_normal((len(rng), 4)), index=rng
    )

    temp_hdfstore["frame"] = frame
    recons = temp_hdfstore["frame"]
    tm.assert_index_equal(recons.index, rng)
    assert rng.tz == recons.index.tz


@td.skip_if_windows
def test_store_timezone(temp_hdfstore):
    # GH2852
    # issue storing datetime.date with a timezone as it resets when read
    # back in a new timezone

    # original method
    today = date(2013, 9, 10)
    df = DataFrame([1, 2, 3], index=[today, today, today])
    temp_hdfstore["obj1"] = df
    result = temp_hdfstore["obj1"]
    tm.assert_frame_equal(result, df)

    # with tz setting
    with tm.set_timezone("EST5EDT"):
        today = date(2013, 9, 10)
        df = DataFrame([1, 2, 3], index=[today, today, today])
        temp_hdfstore["obj2"] = df

    with tm.set_timezone("CST6CDT"):
        result = temp_hdfstore["obj2"]

    tm.assert_frame_equal(result, df)


def test_dst_transitions(temp_hdfstore):
    # make sure we are not failing on transitions
    times = date_range(
        "2013-10-26 23:00",
        "2013-10-27 01:00",
        tz="Europe/London",
        freq="h",
        ambiguous="infer",
    )
    times = times._with_freq(None)  # freq doesn't round-trip

    for i in [times, times + pd.Timedelta("10min")]:
        df = DataFrame({"A": range(len(i)), "B": i}, index=i)
        temp_hdfstore.append("df", df)
        result = temp_hdfstore.select("df")
        tm.assert_frame_equal(result, df)
        temp_hdfstore.remove("df")


@pytest.mark.filterwarnings(
    "ignore:`alltrue` is deprecated as of NumPy 1.25.0:DeprecationWarning"
)
def test_read_with_where_tz_aware_index(temp_hdfstore):
    # GH 11926
    periods = 10
    dts = date_range("20151201", periods=periods, freq="D", tz="UTC", unit="ns")
    mi = pd.MultiIndex.from_arrays([dts, range(periods)], names=["DATE", "NO"])
    expected = DataFrame({"MYCOL": 0}, index=mi)

    key = "mykey"
    with pd.HDFStore(temp_hdfstore) as store:
        store.append(key, expected, format="table", append=True)
    result = pd.read_hdf(temp_hdfstore, key, where="DATE > 20151130")
    tm.assert_frame_equal(result, expected)
