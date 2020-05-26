import datetime

import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
from pandas import DataFrame, DatetimeIndex, Series, Timestamp, date_range
import pandas._testing as tm
from pandas.tests.io.pytables.common import (
    _maybe_remove,
    ensure_clean_path,
    ensure_clean_store,
)


def _compare_with_tz(a, b):
    tm.assert_frame_equal(a, b)

    # compare the zones on each element
    for c in a.columns:
        for i in a.index:
            a_e = a.loc[i, c]
            b_e = b.loc[i, c]
            if not (a_e == b_e and a_e.tz == b_e.tz):
                raise AssertionError(f"invalid tz comparison [{a_e}] [{b_e}]")


def test_append_with_timezones_dateutil(setup_path):

    from datetime import timedelta

    # use maybe_get_tz instead of dateutil.tz.gettz to handle the windows
    # filename issues.
    from pandas._libs.tslibs.timezones import maybe_get_tz

    gettz = lambda x: maybe_get_tz("dateutil/" + x)

    # as columns
    with ensure_clean_store(setup_path) as store:

        _maybe_remove(store, "df_tz")
        df = DataFrame(
            dict(
                A=[
                    Timestamp("20130102 2:00:00", tz=gettz("US/Eastern"))
                    + timedelta(hours=1) * i
                    for i in range(5)
                ]
            )
        )

        store.append("df_tz", df, data_columns=["A"])
        result = store["df_tz"]
        _compare_with_tz(result, df)
        tm.assert_frame_equal(result, df)

        # select with tz aware
        expected = df[df.A >= df.A[3]]
        result = store.select("df_tz", where="A>=df.A[3]")
        _compare_with_tz(result, expected)

        # ensure we include dates in DST and STD time here.
        _maybe_remove(store, "df_tz")
        df = DataFrame(
            dict(
                A=Timestamp("20130102", tz=gettz("US/Eastern")),
                B=Timestamp("20130603", tz=gettz("US/Eastern")),
            ),
            index=range(5),
        )
        store.append("df_tz", df)
        result = store["df_tz"]
        _compare_with_tz(result, df)
        tm.assert_frame_equal(result, df)

        df = DataFrame(
            dict(
                A=Timestamp("20130102", tz=gettz("US/Eastern")),
                B=Timestamp("20130102", tz=gettz("EET")),
            ),
            index=range(5),
        )
        with pytest.raises(ValueError):
            store.append("df_tz", df)

        # this is ok
        _maybe_remove(store, "df_tz")
        store.append("df_tz", df, data_columns=["A", "B"])
        result = store["df_tz"]
        _compare_with_tz(result, df)
        tm.assert_frame_equal(result, df)

        # can't append with diff timezone
        df = DataFrame(
            dict(
                A=Timestamp("20130102", tz=gettz("US/Eastern")),
                B=Timestamp("20130102", tz=gettz("CET")),
            ),
            index=range(5),
        )
        with pytest.raises(ValueError):
            store.append("df_tz", df)

    # as index
    with ensure_clean_store(setup_path) as store:

        dti = date_range("2000-1-1", periods=3, freq="H", tz=gettz("US/Eastern"))
        dti = dti._with_freq(None)  # freq doesnt round-trip

        # GH 4098 example
        df = DataFrame(dict(A=Series(range(3), index=dti,)))

        _maybe_remove(store, "df")
        store.put("df", df)
        result = store.select("df")
        tm.assert_frame_equal(result, df)

        _maybe_remove(store, "df")
        store.append("df", df)
        result = store.select("df")
        tm.assert_frame_equal(result, df)


def test_append_with_timezones_pytz(setup_path):

    from datetime import timedelta

    # as columns
    with ensure_clean_store(setup_path) as store:

        _maybe_remove(store, "df_tz")
        df = DataFrame(
            dict(
                A=[
                    Timestamp("20130102 2:00:00", tz="US/Eastern")
                    + timedelta(hours=1) * i
                    for i in range(5)
                ]
            )
        )
        store.append("df_tz", df, data_columns=["A"])
        result = store["df_tz"]
        _compare_with_tz(result, df)
        tm.assert_frame_equal(result, df)

        # select with tz aware
        _compare_with_tz(store.select("df_tz", where="A>=df.A[3]"), df[df.A >= df.A[3]])

        _maybe_remove(store, "df_tz")
        # ensure we include dates in DST and STD time here.
        df = DataFrame(
            dict(
                A=Timestamp("20130102", tz="US/Eastern"),
                B=Timestamp("20130603", tz="US/Eastern"),
            ),
            index=range(5),
        )
        store.append("df_tz", df)
        result = store["df_tz"]
        _compare_with_tz(result, df)
        tm.assert_frame_equal(result, df)

        df = DataFrame(
            dict(
                A=Timestamp("20130102", tz="US/Eastern"),
                B=Timestamp("20130102", tz="EET"),
            ),
            index=range(5),
        )
        with pytest.raises(ValueError):
            store.append("df_tz", df)

        # this is ok
        _maybe_remove(store, "df_tz")
        store.append("df_tz", df, data_columns=["A", "B"])
        result = store["df_tz"]
        _compare_with_tz(result, df)
        tm.assert_frame_equal(result, df)

        # can't append with diff timezone
        df = DataFrame(
            dict(
                A=Timestamp("20130102", tz="US/Eastern"),
                B=Timestamp("20130102", tz="CET"),
            ),
            index=range(5),
        )
        with pytest.raises(ValueError):
            store.append("df_tz", df)

    # as index
    with ensure_clean_store(setup_path) as store:

        dti = date_range("2000-1-1", periods=3, freq="H", tz="US/Eastern")
        dti = dti._with_freq(None)  # freq doesnt round-trip

        # GH 4098 example
        df = DataFrame(dict(A=Series(range(3), index=dti,)))

        _maybe_remove(store, "df")
        store.put("df", df)
        result = store.select("df")
        tm.assert_frame_equal(result, df)

        _maybe_remove(store, "df")
        store.append("df", df)
        result = store.select("df")
        tm.assert_frame_equal(result, df)


def test_tseries_select_index_column(setup_path):
    # GH7777
    # selecting a UTC datetimeindex column did
    # not preserve UTC tzinfo set before storing

    # check that no tz still works
    rng = date_range("1/1/2000", "1/30/2000")
    frame = DataFrame(np.random.randn(len(rng), 4), index=rng)

    with ensure_clean_store(setup_path) as store:
        store.append("frame", frame)
        result = store.select_column("frame", "index")
        assert rng.tz == DatetimeIndex(result.values).tz

    # check utc
    rng = date_range("1/1/2000", "1/30/2000", tz="UTC")
    frame = DataFrame(np.random.randn(len(rng), 4), index=rng)

    with ensure_clean_store(setup_path) as store:
        store.append("frame", frame)
        result = store.select_column("frame", "index")
        assert rng.tz == result.dt.tz

    # double check non-utc
    rng = date_range("1/1/2000", "1/30/2000", tz="US/Eastern")
    frame = DataFrame(np.random.randn(len(rng), 4), index=rng)

    with ensure_clean_store(setup_path) as store:
        store.append("frame", frame)
        result = store.select_column("frame", "index")
        assert rng.tz == result.dt.tz


def test_timezones_fixed(setup_path):
    with ensure_clean_store(setup_path) as store:

        # index
        rng = date_range("1/1/2000", "1/30/2000", tz="US/Eastern")
        rng = rng._with_freq(None)  # freq doesnt round-trip
        df = DataFrame(np.random.randn(len(rng), 4), index=rng)
        store["df"] = df
        result = store["df"]
        tm.assert_frame_equal(result, df)

        # as data
        # GH11411
        _maybe_remove(store, "df")
        df = DataFrame(
            {
                "A": rng,
                "B": rng.tz_convert("UTC").tz_localize(None),
                "C": rng.tz_convert("CET"),
                "D": range(len(rng)),
            },
            index=rng,
        )
        store["df"] = df
        result = store["df"]
        tm.assert_frame_equal(result, df)


def test_fixed_offset_tz(setup_path):
    rng = date_range("1/1/2000 00:00:00-07:00", "1/30/2000 00:00:00-07:00")
    frame = DataFrame(np.random.randn(len(rng), 4), index=rng)

    with ensure_clean_store(setup_path) as store:
        store["frame"] = frame
        recons = store["frame"]
        tm.assert_index_equal(recons.index, rng)
        assert rng.tz == recons.index.tz


@td.skip_if_windows
def test_store_timezone(setup_path):
    # GH2852
    # issue storing datetime.date with a timezone as it resets when read
    # back in a new timezone

    # original method
    with ensure_clean_store(setup_path) as store:

        today = datetime.date(2013, 9, 10)
        df = DataFrame([1, 2, 3], index=[today, today, today])
        store["obj1"] = df
        result = store["obj1"]
        tm.assert_frame_equal(result, df)

    # with tz setting
    with ensure_clean_store(setup_path) as store:

        with tm.set_timezone("EST5EDT"):
            today = datetime.date(2013, 9, 10)
            df = DataFrame([1, 2, 3], index=[today, today, today])
            store["obj1"] = df

        with tm.set_timezone("CST6CDT"):
            result = store["obj1"]

        tm.assert_frame_equal(result, df)


def test_legacy_datetimetz_object(datapath, setup_path):
    # legacy from < 0.17.0
    # 8260
    expected = DataFrame(
        dict(
            A=Timestamp("20130102", tz="US/Eastern"), B=Timestamp("20130603", tz="CET")
        ),
        index=range(5),
    )
    with ensure_clean_store(
        datapath("io", "data", "legacy_hdf", "datetimetz_object.h5"), mode="r"
    ) as store:
        result = store["df"]
        tm.assert_frame_equal(result, expected)


def test_dst_transitions(setup_path):
    # make sure we are not failing on transitions
    with ensure_clean_store(setup_path) as store:
        times = pd.date_range(
            "2013-10-26 23:00",
            "2013-10-27 01:00",
            tz="Europe/London",
            freq="H",
            ambiguous="infer",
        )
        times = times._with_freq(None)  # freq doesnt round-trip

        for i in [times, times + pd.Timedelta("10min")]:
            _maybe_remove(store, "df")
            df = DataFrame({"A": range(len(i)), "B": i}, index=i)
            store.append("df", df)
            result = store.select("df")
            tm.assert_frame_equal(result, df)


def test_read_with_where_tz_aware_index(setup_path):
    # GH 11926
    periods = 10
    dts = pd.date_range("20151201", periods=periods, freq="D", tz="UTC")
    mi = pd.MultiIndex.from_arrays([dts, range(periods)], names=["DATE", "NO"])
    expected = pd.DataFrame({"MYCOL": 0}, index=mi)

    key = "mykey"
    with ensure_clean_path(setup_path) as path:
        with pd.HDFStore(path) as store:
            store.append(key, expected, format="table", append=True)
        result = pd.read_hdf(path, key, where="DATE > 20151130")
        tm.assert_frame_equal(result, expected)


def test_py2_created_with_datetimez(datapath, setup_path):
    # The test HDF5 file was created in Python 2, but could not be read in
    # Python 3.
    #
    # GH26443
    index = [pd.Timestamp("2019-01-01T18:00").tz_localize("America/New_York")]
    expected = DataFrame({"data": 123}, index=index)
    with ensure_clean_store(
        datapath("io", "data", "legacy_hdf", "gh26443.h5"), mode="r"
    ) as store:
        result = store["key"]
        tm.assert_frame_equal(result, expected)
