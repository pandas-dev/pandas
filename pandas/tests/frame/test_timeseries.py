from datetime import datetime, time
from itertools import product

import numpy as np
import pytest
import pytz

import pandas as pd
from pandas import (
    DataFrame,
    DatetimeIndex,
    Index,
    MultiIndex,
    Series,
    date_range,
    period_range,
    to_datetime,
)
import pandas._testing as tm

import pandas.tseries.offsets as offsets


@pytest.fixture(params=product([True, False], [True, False]))
def close_open_fixture(request):
    return request.param


class TestDataFrameTimeSeriesMethods:
    def test_frame_ctor_datetime64_column(self):
        rng = date_range("1/1/2000 00:00:00", "1/1/2000 1:59:50", freq="10s")
        dates = np.asarray(rng)

        df = DataFrame({"A": np.random.randn(len(rng)), "B": dates})
        assert np.issubdtype(df["B"].dtype, np.dtype("M8[ns]"))

    def test_frame_append_datetime64_column(self):
        rng = date_range("1/1/2000 00:00:00", "1/1/2000 1:59:50", freq="10s")
        df = DataFrame(index=np.arange(len(rng)))

        df["A"] = rng
        assert np.issubdtype(df["A"].dtype, np.dtype("M8[ns]"))

    def test_frame_datetime64_pre1900_repr(self):
        df = DataFrame({"year": date_range("1/1/1700", periods=50, freq="A-DEC")})
        # it works!
        repr(df)

    def test_frame_append_datetime64_col_other_units(self):
        n = 100

        units = ["h", "m", "s", "ms", "D", "M", "Y"]

        ns_dtype = np.dtype("M8[ns]")

        for unit in units:
            dtype = np.dtype(f"M8[{unit}]")
            vals = np.arange(n, dtype=np.int64).view(dtype)

            df = DataFrame({"ints": np.arange(n)}, index=np.arange(n))
            df[unit] = vals

            ex_vals = to_datetime(vals.astype("O")).values

            assert df[unit].dtype == ns_dtype
            assert (df[unit].values == ex_vals).all()

        # Test insertion into existing datetime64 column
        df = DataFrame({"ints": np.arange(n)}, index=np.arange(n))
        df["dates"] = np.arange(n, dtype=np.int64).view(ns_dtype)

        for unit in units:
            dtype = np.dtype(f"M8[{unit}]")
            vals = np.arange(n, dtype=np.int64).view(dtype)

            tmp = df.copy()

            tmp["dates"] = vals
            ex_vals = to_datetime(vals.astype("O")).values

            assert (tmp["dates"].values == ex_vals).all()

    def test_asfreq(self, datetime_frame):
        offset_monthly = datetime_frame.asfreq(offsets.BMonthEnd())
        rule_monthly = datetime_frame.asfreq("BM")

        tm.assert_almost_equal(offset_monthly["A"], rule_monthly["A"])

        filled = rule_monthly.asfreq("B", method="pad")  # noqa
        # TODO: actually check that this worked.

        # don't forget!
        filled_dep = rule_monthly.asfreq("B", method="pad")  # noqa

        # test does not blow up on length-0 DataFrame
        zero_length = datetime_frame.reindex([])
        result = zero_length.asfreq("BM")
        assert result is not zero_length

    def test_asfreq_datetimeindex(self):
        df = DataFrame(
            {"A": [1, 2, 3]},
            index=[datetime(2011, 11, 1), datetime(2011, 11, 2), datetime(2011, 11, 3)],
        )
        df = df.asfreq("B")
        assert isinstance(df.index, DatetimeIndex)

        ts = df["A"].asfreq("B")
        assert isinstance(ts.index, DatetimeIndex)

    def test_asfreq_fillvalue(self):
        # test for fill value during upsampling, related to issue 3715

        # setup
        rng = pd.date_range("1/1/2016", periods=10, freq="2S")
        ts = pd.Series(np.arange(len(rng)), index=rng)
        df = pd.DataFrame({"one": ts})

        # insert pre-existing missing value
        df.loc["2016-01-01 00:00:08", "one"] = None

        actual_df = df.asfreq(freq="1S", fill_value=9.0)
        expected_df = df.asfreq(freq="1S").fillna(9.0)
        expected_df.loc["2016-01-01 00:00:08", "one"] = None
        tm.assert_frame_equal(expected_df, actual_df)

        expected_series = ts.asfreq(freq="1S").fillna(9.0)
        actual_series = ts.asfreq(freq="1S", fill_value=9.0)
        tm.assert_series_equal(expected_series, actual_series)

    @pytest.mark.parametrize(
        "data,idx,expected_first,expected_last",
        [
            ({"A": [1, 2, 3]}, [1, 1, 2], 1, 2),
            ({"A": [1, 2, 3]}, [1, 2, 2], 1, 2),
            ({"A": [1, 2, 3, 4]}, ["d", "d", "d", "d"], "d", "d"),
            ({"A": [1, np.nan, 3]}, [1, 1, 2], 1, 2),
            ({"A": [np.nan, np.nan, 3]}, [1, 1, 2], 2, 2),
            ({"A": [1, np.nan, 3]}, [1, 2, 2], 1, 2),
        ],
    )
    def test_first_last_valid(
        self, float_frame, data, idx, expected_first, expected_last
    ):
        N = len(float_frame.index)
        mat = np.random.randn(N)
        mat[:5] = np.nan
        mat[-5:] = np.nan

        frame = DataFrame({"foo": mat}, index=float_frame.index)
        index = frame.first_valid_index()

        assert index == frame.index[5]

        index = frame.last_valid_index()
        assert index == frame.index[-6]

        # GH12800
        empty = DataFrame()
        assert empty.last_valid_index() is None
        assert empty.first_valid_index() is None

        # GH17400: no valid entries
        frame[:] = np.nan
        assert frame.last_valid_index() is None
        assert frame.first_valid_index() is None

        # GH20499: its preserves freq with holes
        frame.index = date_range("20110101", periods=N, freq="B")
        frame.iloc[1] = 1
        frame.iloc[-2] = 1
        assert frame.first_valid_index() == frame.index[1]
        assert frame.last_valid_index() == frame.index[-2]
        assert frame.first_valid_index().freq == frame.index.freq
        assert frame.last_valid_index().freq == frame.index.freq

        # GH 21441
        df = DataFrame(data, index=idx)
        assert expected_first == df.first_valid_index()
        assert expected_last == df.last_valid_index()

    @pytest.mark.parametrize("klass", [Series, DataFrame])
    def test_first_valid_index_all_nan(self, klass):
        # GH#9752 Series/DataFrame should both return None, not raise
        obj = klass([np.nan])

        assert obj.first_valid_index() is None
        assert obj.iloc[:0].first_valid_index() is None

    def test_first_subset(self):
        ts = tm.makeTimeDataFrame(freq="12h")
        result = ts.first("10d")
        assert len(result) == 20

        ts = tm.makeTimeDataFrame(freq="D")
        result = ts.first("10d")
        assert len(result) == 10

        result = ts.first("3M")
        expected = ts[:"3/31/2000"]
        tm.assert_frame_equal(result, expected)

        result = ts.first("21D")
        expected = ts[:21]
        tm.assert_frame_equal(result, expected)

        result = ts[:0].first("3M")
        tm.assert_frame_equal(result, ts[:0])

    def test_first_raises(self):
        # GH20725
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(TypeError):  # index is not a DatetimeIndex
            df.first("1D")

    def test_last_subset(self):
        ts = tm.makeTimeDataFrame(freq="12h")
        result = ts.last("10d")
        assert len(result) == 20

        ts = tm.makeTimeDataFrame(nper=30, freq="D")
        result = ts.last("10d")
        assert len(result) == 10

        result = ts.last("21D")
        expected = ts["2000-01-10":]
        tm.assert_frame_equal(result, expected)

        result = ts.last("21D")
        expected = ts[-21:]
        tm.assert_frame_equal(result, expected)

        result = ts[:0].last("3M")
        tm.assert_frame_equal(result, ts[:0])

    def test_last_raises(self):
        # GH20725
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(TypeError):  # index is not a DatetimeIndex
            df.last("1D")

    def test_at_time(self):
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)
        rs = ts.at_time(rng[1])
        assert (rs.index.hour == rng[1].hour).all()
        assert (rs.index.minute == rng[1].minute).all()
        assert (rs.index.second == rng[1].second).all()

        result = ts.at_time("9:30")
        expected = ts.at_time(time(9, 30))
        tm.assert_frame_equal(result, expected)

        result = ts.loc[time(9, 30)]
        expected = ts.loc[(rng.hour == 9) & (rng.minute == 30)]

        tm.assert_frame_equal(result, expected)

        # midnight, everything
        rng = date_range("1/1/2000", "1/31/2000")
        ts = DataFrame(np.random.randn(len(rng), 3), index=rng)

        result = ts.at_time(time(0, 0))
        tm.assert_frame_equal(result, ts)

        # time doesn't exist
        rng = date_range("1/1/2012", freq="23Min", periods=384)
        ts = DataFrame(np.random.randn(len(rng), 2), rng)
        rs = ts.at_time("16:00")
        assert len(rs) == 0

    @pytest.mark.parametrize(
        "hour", ["1:00", "1:00AM", time(1), time(1, tzinfo=pytz.UTC)]
    )
    def test_at_time_errors(self, hour):
        # GH 24043
        dti = pd.date_range("2018", periods=3, freq="H")
        df = pd.DataFrame(list(range(len(dti))), index=dti)
        if getattr(hour, "tzinfo", None) is None:
            result = df.at_time(hour)
            expected = df.iloc[1:2]
            tm.assert_frame_equal(result, expected)
        else:
            with pytest.raises(ValueError, match="Index must be timezone"):
                df.at_time(hour)

    def test_at_time_tz(self):
        # GH 24043
        dti = pd.date_range("2018", periods=3, freq="H", tz="US/Pacific")
        df = pd.DataFrame(list(range(len(dti))), index=dti)
        result = df.at_time(time(4, tzinfo=pytz.timezone("US/Eastern")))
        expected = df.iloc[1:2]
        tm.assert_frame_equal(result, expected)

    def test_at_time_raises(self):
        # GH20725
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(TypeError):  # index is not a DatetimeIndex
            df.at_time("00:00")

    @pytest.mark.parametrize("axis", ["index", "columns", 0, 1])
    def test_at_time_axis(self, axis):
        # issue 8839
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        ts = DataFrame(np.random.randn(len(rng), len(rng)))
        ts.index, ts.columns = rng, rng

        indices = rng[(rng.hour == 9) & (rng.minute == 30) & (rng.second == 0)]

        if axis in ["index", 0]:
            expected = ts.loc[indices, :]
        elif axis in ["columns", 1]:
            expected = ts.loc[:, indices]

        result = ts.at_time("9:30", axis=axis)
        tm.assert_frame_equal(result, expected)

    def test_between_time(self, close_open_fixture):
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)
        stime = time(0, 0)
        etime = time(1, 0)
        inc_start, inc_end = close_open_fixture

        filtered = ts.between_time(stime, etime, inc_start, inc_end)
        exp_len = 13 * 4 + 1
        if not inc_start:
            exp_len -= 5
        if not inc_end:
            exp_len -= 4

        assert len(filtered) == exp_len
        for rs in filtered.index:
            t = rs.time()
            if inc_start:
                assert t >= stime
            else:
                assert t > stime

            if inc_end:
                assert t <= etime
            else:
                assert t < etime

        result = ts.between_time("00:00", "01:00")
        expected = ts.between_time(stime, etime)
        tm.assert_frame_equal(result, expected)

        # across midnight
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)
        stime = time(22, 0)
        etime = time(9, 0)

        filtered = ts.between_time(stime, etime, inc_start, inc_end)
        exp_len = (12 * 11 + 1) * 4 + 1
        if not inc_start:
            exp_len -= 4
        if not inc_end:
            exp_len -= 4

        assert len(filtered) == exp_len
        for rs in filtered.index:
            t = rs.time()
            if inc_start:
                assert (t >= stime) or (t <= etime)
            else:
                assert (t > stime) or (t <= etime)

            if inc_end:
                assert (t <= etime) or (t >= stime)
            else:
                assert (t < etime) or (t >= stime)

    def test_between_time_raises(self):
        # GH20725
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(TypeError):  # index is not a DatetimeIndex
            df.between_time(start_time="00:00", end_time="12:00")

    def test_between_time_axis(self, axis):
        # issue 8839
        rng = date_range("1/1/2000", periods=100, freq="10min")
        ts = DataFrame(np.random.randn(len(rng), len(rng)))
        stime, etime = ("08:00:00", "09:00:00")
        exp_len = 7

        if axis in ["index", 0]:
            ts.index = rng
            assert len(ts.between_time(stime, etime)) == exp_len
            assert len(ts.between_time(stime, etime, axis=0)) == exp_len

        if axis in ["columns", 1]:
            ts.columns = rng
            selected = ts.between_time(stime, etime, axis=1).columns
            assert len(selected) == exp_len

    def test_between_time_axis_raises(self, axis):
        # issue 8839
        rng = date_range("1/1/2000", periods=100, freq="10min")
        mask = np.arange(0, len(rng))
        rand_data = np.random.randn(len(rng), len(rng))
        ts = DataFrame(rand_data, index=rng, columns=rng)
        stime, etime = ("08:00:00", "09:00:00")

        msg = "Index must be DatetimeIndex"
        if axis in ["columns", 1]:
            ts.index = mask
            with pytest.raises(TypeError, match=msg):
                ts.between_time(stime, etime)
            with pytest.raises(TypeError, match=msg):
                ts.between_time(stime, etime, axis=0)

        if axis in ["index", 0]:
            ts.columns = mask
            with pytest.raises(TypeError, match=msg):
                ts.between_time(stime, etime, axis=1)

    def test_operation_on_NaT(self):
        # Both NaT and Timestamp are in DataFrame.
        df = pd.DataFrame({"foo": [pd.NaT, pd.NaT, pd.Timestamp("2012-05-01")]})

        res = df.min()
        exp = pd.Series([pd.Timestamp("2012-05-01")], index=["foo"])
        tm.assert_series_equal(res, exp)

        res = df.max()
        exp = pd.Series([pd.Timestamp("2012-05-01")], index=["foo"])
        tm.assert_series_equal(res, exp)

        # GH12941, only NaTs are in DataFrame.
        df = pd.DataFrame({"foo": [pd.NaT, pd.NaT]})

        res = df.min()
        exp = pd.Series([pd.NaT], index=["foo"])
        tm.assert_series_equal(res, exp)

        res = df.max()
        exp = pd.Series([pd.NaT], index=["foo"])
        tm.assert_series_equal(res, exp)

    def test_datetime_assignment_with_NaT_and_diff_time_units(self):
        # GH 7492
        data_ns = np.array([1, "nat"], dtype="datetime64[ns]")
        result = pd.Series(data_ns).to_frame()
        result["new"] = data_ns
        expected = pd.DataFrame(
            {0: [1, None], "new": [1, None]}, dtype="datetime64[ns]"
        )
        tm.assert_frame_equal(result, expected)
        # OutOfBoundsDatetime error shouldn't occur
        data_s = np.array([1, "nat"], dtype="datetime64[s]")
        result["new"] = data_s
        expected = pd.DataFrame(
            {0: [1, None], "new": [1e9, None]}, dtype="datetime64[ns]"
        )
        tm.assert_frame_equal(result, expected)

    def test_frame_to_period(self):
        K = 5

        dr = date_range("1/1/2000", "1/1/2001")
        pr = period_range("1/1/2000", "1/1/2001")
        df = DataFrame(np.random.randn(len(dr), K), index=dr)
        df["mix"] = "a"

        pts = df.to_period()
        exp = df.copy()
        exp.index = pr
        tm.assert_frame_equal(pts, exp)

        pts = df.to_period("M")
        tm.assert_index_equal(pts.index, exp.index.asfreq("M"))

        df = df.T
        pts = df.to_period(axis=1)
        exp = df.copy()
        exp.columns = pr
        tm.assert_frame_equal(pts, exp)

        pts = df.to_period("M", axis=1)
        tm.assert_index_equal(pts.columns, exp.columns.asfreq("M"))

        msg = "No axis named 2 for object type <class 'pandas.core.frame.DataFrame'>"
        with pytest.raises(ValueError, match=msg):
            df.to_period(axis=2)

    @pytest.mark.parametrize("fn", ["tz_localize", "tz_convert"])
    def test_tz_convert_and_localize(self, fn):
        l0 = date_range("20140701", periods=5, freq="D")
        l1 = date_range("20140701", periods=5, freq="D")

        int_idx = Index(range(5))

        if fn == "tz_convert":
            l0 = l0.tz_localize("UTC")
            l1 = l1.tz_localize("UTC")

        for idx in [l0, l1]:

            l0_expected = getattr(idx, fn)("US/Pacific")
            l1_expected = getattr(idx, fn)("US/Pacific")

            df1 = DataFrame(np.ones(5), index=l0)
            df1 = getattr(df1, fn)("US/Pacific")
            tm.assert_index_equal(df1.index, l0_expected)

            # MultiIndex
            # GH7846
            df2 = DataFrame(np.ones(5), MultiIndex.from_arrays([l0, l1]))

            df3 = getattr(df2, fn)("US/Pacific", level=0)
            assert not df3.index.levels[0].equals(l0)
            tm.assert_index_equal(df3.index.levels[0], l0_expected)
            tm.assert_index_equal(df3.index.levels[1], l1)
            assert not df3.index.levels[1].equals(l1_expected)

            df3 = getattr(df2, fn)("US/Pacific", level=1)
            tm.assert_index_equal(df3.index.levels[0], l0)
            assert not df3.index.levels[0].equals(l0_expected)
            tm.assert_index_equal(df3.index.levels[1], l1_expected)
            assert not df3.index.levels[1].equals(l1)

            df4 = DataFrame(np.ones(5), MultiIndex.from_arrays([int_idx, l0]))

            # TODO: untested
            df5 = getattr(df4, fn)("US/Pacific", level=1)  # noqa

            tm.assert_index_equal(df3.index.levels[0], l0)
            assert not df3.index.levels[0].equals(l0_expected)
            tm.assert_index_equal(df3.index.levels[1], l1_expected)
            assert not df3.index.levels[1].equals(l1)

        # Bad Inputs

        # Not DatetimeIndex / PeriodIndex
        with pytest.raises(TypeError, match="DatetimeIndex"):
            df = DataFrame(index=int_idx)
            df = getattr(df, fn)("US/Pacific")

        # Not DatetimeIndex / PeriodIndex
        with pytest.raises(TypeError, match="DatetimeIndex"):
            df = DataFrame(np.ones(5), MultiIndex.from_arrays([int_idx, l0]))
            df = getattr(df, fn)("US/Pacific", level=0)

        # Invalid level
        with pytest.raises(ValueError, match="not valid"):
            df = DataFrame(index=l0)
            df = getattr(df, fn)("US/Pacific", level=1)
