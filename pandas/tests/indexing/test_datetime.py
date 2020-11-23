import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Index, Series, Timestamp, date_range
import pandas._testing as tm


class TestDatetimeIndex:
    def test_indexing_with_datetime_tz(self):

        # GH#8260
        # support datetime64 with tz

        idx = Index(date_range("20130101", periods=3, tz="US/Eastern"), name="foo")
        dr = date_range("20130110", periods=3)
        df = DataFrame({"A": idx, "B": dr})
        df["C"] = idx
        df.iloc[1, 1] = pd.NaT
        df.iloc[1, 2] = pd.NaT

        # indexing
        result = df.iloc[1]
        expected = Series(
            [Timestamp("2013-01-02 00:00:00-0500", tz="US/Eastern"), pd.NaT, pd.NaT],
            index=list("ABC"),
            dtype="object",
            name=1,
        )
        tm.assert_series_equal(result, expected)
        result = df.loc[1]
        expected = Series(
            [Timestamp("2013-01-02 00:00:00-0500", tz="US/Eastern"), pd.NaT, pd.NaT],
            index=list("ABC"),
            dtype="object",
            name=1,
        )
        tm.assert_series_equal(result, expected)

        # indexing - fast_xs
        df = DataFrame({"a": date_range("2014-01-01", periods=10, tz="UTC")})
        result = df.iloc[5]
        expected = Series(
            [Timestamp("2014-01-06 00:00:00+0000", tz="UTC")], index=["a"], name=5
        )
        tm.assert_series_equal(result, expected)

        result = df.loc[5]
        tm.assert_series_equal(result, expected)

        # indexing - boolean
        result = df[df.a > df.a[3]]
        expected = df.iloc[4:]
        tm.assert_frame_equal(result, expected)

        # indexing - setting an element
        df = DataFrame(
            data=pd.to_datetime(["2015-03-30 20:12:32", "2015-03-12 00:11:11"]),
            columns=["time"],
        )
        df["new_col"] = ["new", "old"]
        df.time = df.set_index("time").index.tz_localize("UTC")
        v = df[df.new_col == "new"].set_index("time").index.tz_convert("US/Pacific")

        # trying to set a single element on a part of a different timezone
        # this converts to object
        df2 = df.copy()
        df2.loc[df2.new_col == "new", "time"] = v

        expected = Series([v[0], df.loc[1, "time"]], name="time")
        tm.assert_series_equal(df2.time, expected)

        v = df.loc[df.new_col == "new", "time"] + pd.Timedelta("1s")
        df.loc[df.new_col == "new", "time"] = v
        tm.assert_series_equal(df.loc[df.new_col == "new", "time"], v)

    def test_consistency_with_tz_aware_scalar(self):
        # xef gh-12938
        # various ways of indexing the same tz-aware scalar
        df = Series([Timestamp("2016-03-30 14:35:25", tz="Europe/Brussels")]).to_frame()

        df = pd.concat([df, df]).reset_index(drop=True)
        expected = Timestamp("2016-03-30 14:35:25+0200", tz="Europe/Brussels")

        result = df[0][0]
        assert result == expected

        result = df.iloc[0, 0]
        assert result == expected

        result = df.loc[0, 0]
        assert result == expected

        result = df.iat[0, 0]
        assert result == expected

        result = df.at[0, 0]
        assert result == expected

        result = df[0].loc[0]
        assert result == expected

        result = df[0].at[0]
        assert result == expected

    def test_indexing_with_datetimeindex_tz(self):

        # GH 12050
        # indexing on a series with a datetimeindex with tz
        index = date_range("2015-01-01", periods=2, tz="utc")

        ser = Series(range(2), index=index, dtype="int64")

        # list-like indexing

        for sel in (index, list(index)):
            # getitem
            result = ser[sel]
            expected = ser.copy()
            if sel is not index:
                expected.index = expected.index._with_freq(None)
            tm.assert_series_equal(result, expected)

            # setitem
            result = ser.copy()
            result[sel] = 1
            expected = Series(1, index=index)
            tm.assert_series_equal(result, expected)

            # .loc getitem
            result = ser.loc[sel]
            expected = ser.copy()
            if sel is not index:
                expected.index = expected.index._with_freq(None)
            tm.assert_series_equal(result, expected)

            # .loc setitem
            result = ser.copy()
            result.loc[sel] = 1
            expected = Series(1, index=index)
            tm.assert_series_equal(result, expected)

        # single element indexing

        # getitem
        assert ser[index[1]] == 1

        # setitem
        result = ser.copy()
        result[index[1]] = 5
        expected = Series([0, 5], index=index)
        tm.assert_series_equal(result, expected)

        # .loc getitem
        assert ser.loc[index[1]] == 1

        # .loc setitem
        result = ser.copy()
        result.loc[index[1]] = 5
        expected = Series([0, 5], index=index)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("to_period", [True, False])
    def test_loc_getitem_listlike_of_datetimelike_keys(self, to_period):
        # GH 11497

        idx = date_range("2011-01-01", "2011-01-02", freq="D", name="idx")
        if to_period:
            idx = idx.to_period("D")
        ser = Series([0.1, 0.2], index=idx, name="s")

        keys = [Timestamp("2011-01-01"), Timestamp("2011-01-02")]
        if to_period:
            keys = [x.to_period("D") for x in keys]
        result = ser.loc[keys]
        exp = Series([0.1, 0.2], index=idx, name="s")
        if not to_period:
            exp.index = exp.index._with_freq(None)
        tm.assert_series_equal(result, exp, check_index_type=True)

        keys = [
            Timestamp("2011-01-02"),
            Timestamp("2011-01-02"),
            Timestamp("2011-01-01"),
        ]
        if to_period:
            keys = [x.to_period("D") for x in keys]
        exp = Series(
            [0.2, 0.2, 0.1], index=Index(keys, name="idx", dtype=idx.dtype), name="s"
        )
        result = ser.loc[keys]
        tm.assert_series_equal(result, exp, check_index_type=True)

        keys = [
            Timestamp("2011-01-03"),
            Timestamp("2011-01-02"),
            Timestamp("2011-01-03"),
        ]
        if to_period:
            keys = [x.to_period("D") for x in keys]

        with pytest.raises(KeyError, match="with any missing labels"):
            ser.loc[keys]

    def test_nanosecond_getitem_setitem_with_tz(self):
        # GH 11679
        data = ["2016-06-28 08:30:00.123456789"]
        index = pd.DatetimeIndex(data, dtype="datetime64[ns, America/Chicago]")
        df = DataFrame({"a": [10]}, index=index)
        result = df.loc[df.index[0]]
        expected = Series(10, index=["a"], name=df.index[0])
        tm.assert_series_equal(result, expected)

        result = df.copy()
        result.loc[df.index[0], "a"] = -1
        expected = DataFrame(-1, index=index, columns=["a"])
        tm.assert_frame_equal(result, expected)

    def test_loc_setitem_with_existing_dst(self):
        # GH 18308
        start = Timestamp("2017-10-29 00:00:00+0200", tz="Europe/Madrid")
        end = Timestamp("2017-10-29 03:00:00+0100", tz="Europe/Madrid")
        ts = Timestamp("2016-10-10 03:00:00", tz="Europe/Madrid")
        idx = pd.date_range(start, end, closed="left", freq="H")
        result = DataFrame(index=idx, columns=["value"])
        result.loc[ts, "value"] = 12
        expected = DataFrame(
            [np.nan] * len(idx) + [12],
            index=idx.append(pd.DatetimeIndex([ts])),
            columns=["value"],
            dtype=object,
        )
        tm.assert_frame_equal(result, expected)
