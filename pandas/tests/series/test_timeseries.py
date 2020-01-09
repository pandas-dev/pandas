from datetime import datetime, time, timedelta
from io import StringIO
from itertools import product

import numpy as np
import pytest

from pandas._libs.tslib import iNaT
from pandas._libs.tslibs.np_datetime import OutOfBoundsDatetime
import pandas.util._test_decorators as td

import pandas as pd
from pandas import (
    DataFrame,
    DatetimeIndex,
    NaT,
    Series,
    Timestamp,
    concat,
    date_range,
    timedelta_range,
    to_datetime,
)
import pandas._testing as tm

from pandas.tseries.offsets import BDay, BMonthEnd


def _simple_ts(start, end, freq="D"):
    rng = date_range(start, end, freq=freq)
    return Series(np.random.randn(len(rng)), index=rng)


def assert_range_equal(left, right):
    assert left.equals(right)
    assert left.freq == right.freq
    assert left.tz == right.tz


class TestTimeSeries:
    def test_asfreq(self):
        ts = Series(
            [0.0, 1.0, 2.0],
            index=[
                datetime(2009, 10, 30),
                datetime(2009, 11, 30),
                datetime(2009, 12, 31),
            ],
        )

        daily_ts = ts.asfreq("B")
        monthly_ts = daily_ts.asfreq("BM")
        tm.assert_series_equal(monthly_ts, ts)

        daily_ts = ts.asfreq("B", method="pad")
        monthly_ts = daily_ts.asfreq("BM")
        tm.assert_series_equal(monthly_ts, ts)

        daily_ts = ts.asfreq(BDay())
        monthly_ts = daily_ts.asfreq(BMonthEnd())
        tm.assert_series_equal(monthly_ts, ts)

        result = ts[:0].asfreq("M")
        assert len(result) == 0
        assert result is not ts

        daily_ts = ts.asfreq("D", fill_value=-1)
        result = daily_ts.value_counts().sort_index()
        expected = Series([60, 1, 1, 1], index=[-1.0, 2.0, 1.0, 0.0]).sort_index()
        tm.assert_series_equal(result, expected)

    def test_asfreq_datetimeindex_empty_series(self):
        # GH 14320
        index = pd.DatetimeIndex(["2016-09-29 11:00"])
        expected = Series(index=index, dtype=object).asfreq("H")
        result = Series([3], index=index.copy()).asfreq("H")
        tm.assert_index_equal(expected.index, result.index)

    def test_autocorr(self, datetime_series):
        # Just run the function
        corr1 = datetime_series.autocorr()

        # Now run it with the lag parameter
        corr2 = datetime_series.autocorr(lag=1)

        # corr() with lag needs Series of at least length 2
        if len(datetime_series) <= 2:
            assert np.isnan(corr1)
            assert np.isnan(corr2)
        else:
            assert corr1 == corr2

        # Choose a random lag between 1 and length of Series - 2
        # and compare the result with the Series corr() function
        n = 1 + np.random.randint(max(1, len(datetime_series) - 2))
        corr1 = datetime_series.corr(datetime_series.shift(n))
        corr2 = datetime_series.autocorr(lag=n)

        # corr() with lag needs Series of at least length 2
        if len(datetime_series) <= 2:
            assert np.isnan(corr1)
            assert np.isnan(corr2)
        else:
            assert corr1 == corr2

    def test_first_last_valid(self, datetime_series):
        ts = datetime_series.copy()
        ts[:5] = np.NaN

        index = ts.first_valid_index()
        assert index == ts.index[5]

        ts[-5:] = np.NaN
        index = ts.last_valid_index()
        assert index == ts.index[-6]

        ts[:] = np.nan
        assert ts.last_valid_index() is None
        assert ts.first_valid_index() is None

        ser = Series([], index=[], dtype=object)
        assert ser.last_valid_index() is None
        assert ser.first_valid_index() is None

        # GH12800
        empty = Series(dtype=object)
        assert empty.last_valid_index() is None
        assert empty.first_valid_index() is None

        # GH20499: its preserves freq with holes
        ts.index = date_range("20110101", periods=len(ts), freq="B")
        ts.iloc[1] = 1
        ts.iloc[-2] = 1
        assert ts.first_valid_index() == ts.index[1]
        assert ts.last_valid_index() == ts.index[-2]
        assert ts.first_valid_index().freq == ts.index.freq
        assert ts.last_valid_index().freq == ts.index.freq

    def test_mpl_compat_hack(self, datetime_series):
        result = datetime_series[:, np.newaxis]
        expected = datetime_series.values[:, np.newaxis]
        tm.assert_almost_equal(result, expected)

    def test_timeseries_coercion(self):
        idx = tm.makeDateIndex(10000)
        ser = Series(np.random.randn(len(idx)), idx.astype(object))
        assert ser.index.is_all_dates
        assert isinstance(ser.index, DatetimeIndex)

    def test_contiguous_boolean_preserve_freq(self):
        rng = date_range("1/1/2000", "3/1/2000", freq="B")

        mask = np.zeros(len(rng), dtype=bool)
        mask[10:20] = True

        masked = rng[mask]
        expected = rng[10:20]
        assert expected.freq is not None
        assert_range_equal(masked, expected)

        mask[22] = True
        masked = rng[mask]
        assert masked.freq is None

    def test_to_datetime_unit(self):

        epoch = 1370745748
        s = Series([epoch + t for t in range(20)])
        result = to_datetime(s, unit="s")
        expected = Series(
            [Timestamp("2013-06-09 02:42:28") + timedelta(seconds=t) for t in range(20)]
        )
        tm.assert_series_equal(result, expected)

        s = Series([epoch + t for t in range(20)]).astype(float)
        result = to_datetime(s, unit="s")
        expected = Series(
            [Timestamp("2013-06-09 02:42:28") + timedelta(seconds=t) for t in range(20)]
        )
        tm.assert_series_equal(result, expected)

        s = Series([epoch + t for t in range(20)] + [iNaT])
        result = to_datetime(s, unit="s")
        expected = Series(
            [Timestamp("2013-06-09 02:42:28") + timedelta(seconds=t) for t in range(20)]
            + [NaT]
        )
        tm.assert_series_equal(result, expected)

        s = Series([epoch + t for t in range(20)] + [iNaT]).astype(float)
        result = to_datetime(s, unit="s")
        expected = Series(
            [Timestamp("2013-06-09 02:42:28") + timedelta(seconds=t) for t in range(20)]
            + [NaT]
        )
        tm.assert_series_equal(result, expected)

        # GH13834
        s = Series([epoch + t for t in np.arange(0, 2, 0.25)] + [iNaT]).astype(float)
        result = to_datetime(s, unit="s")
        expected = Series(
            [
                Timestamp("2013-06-09 02:42:28") + timedelta(seconds=t)
                for t in np.arange(0, 2, 0.25)
            ]
            + [NaT]
        )
        tm.assert_series_equal(result, expected)

        s = concat(
            [Series([epoch + t for t in range(20)]).astype(float), Series([np.nan])],
            ignore_index=True,
        )
        result = to_datetime(s, unit="s")
        expected = Series(
            [Timestamp("2013-06-09 02:42:28") + timedelta(seconds=t) for t in range(20)]
            + [NaT]
        )
        tm.assert_series_equal(result, expected)

        result = to_datetime([1, 2, "NaT", pd.NaT, np.nan], unit="D")
        expected = DatetimeIndex(
            [Timestamp("1970-01-02"), Timestamp("1970-01-03")] + ["NaT"] * 3
        )
        tm.assert_index_equal(result, expected)

        msg = "non convertible value foo with the unit 'D'"
        with pytest.raises(ValueError, match=msg):
            to_datetime([1, 2, "foo"], unit="D")
        msg = "cannot convert input 111111111 with the unit 'D'"
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            to_datetime([1, 2, 111111111], unit="D")

        # coerce we can process
        expected = DatetimeIndex(
            [Timestamp("1970-01-02"), Timestamp("1970-01-03")] + ["NaT"] * 1
        )
        result = to_datetime([1, 2, "foo"], unit="D", errors="coerce")
        tm.assert_index_equal(result, expected)

        result = to_datetime([1, 2, 111111111], unit="D", errors="coerce")
        tm.assert_index_equal(result, expected)

    def test_series_ctor_datetime64(self):
        rng = date_range("1/1/2000 00:00:00", "1/1/2000 1:59:50", freq="10s")
        dates = np.asarray(rng)

        series = Series(dates)
        assert np.issubdtype(series.dtype, np.dtype("M8[ns]"))

    def test_series_repr_nat(self):
        series = Series([0, 1000, 2000, iNaT], dtype="M8[ns]")

        result = repr(series)
        expected = (
            "0   1970-01-01 00:00:00.000000\n"
            "1   1970-01-01 00:00:00.000001\n"
            "2   1970-01-01 00:00:00.000002\n"
            "3                          NaT\n"
            "dtype: datetime64[ns]"
        )
        assert result == expected

    def test_asfreq_keep_index_name(self):
        # GH #9854
        index_name = "bar"
        index = pd.date_range("20130101", periods=20, name=index_name)
        df = pd.DataFrame(list(range(20)), columns=["foo"], index=index)

        assert index_name == df.index.name
        assert index_name == df.asfreq("10D").index.name

    def test_promote_datetime_date(self):
        rng = date_range("1/1/2000", periods=20)
        ts = Series(np.random.randn(20), index=rng)

        ts_slice = ts[5:]
        ts2 = ts_slice.copy()
        ts2.index = [x.date() for x in ts2.index]

        result = ts + ts2
        result2 = ts2 + ts
        expected = ts + ts[5:]
        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result2, expected)

        # test asfreq
        result = ts2.asfreq("4H", method="ffill")
        expected = ts[5:].asfreq("4H", method="ffill")
        tm.assert_series_equal(result, expected)

        result = rng.get_indexer(ts2.index)
        expected = rng.get_indexer(ts_slice.index)
        tm.assert_numpy_array_equal(result, expected)

    def test_asfreq_normalize(self):
        rng = date_range("1/1/2000 09:30", periods=20)
        norm = date_range("1/1/2000", periods=20)
        vals = np.random.randn(20)
        ts = Series(vals, index=rng)

        result = ts.asfreq("D", normalize=True)
        norm = date_range("1/1/2000", periods=20)
        expected = Series(vals, index=norm)

        tm.assert_series_equal(result, expected)

        vals = np.random.randn(20, 3)
        ts = DataFrame(vals, index=rng)

        result = ts.asfreq("D", normalize=True)
        expected = DataFrame(vals, index=norm)

        tm.assert_frame_equal(result, expected)

    def test_first_subset(self):
        ts = _simple_ts("1/1/2000", "1/1/2010", freq="12h")
        result = ts.first("10d")
        assert len(result) == 20

        ts = _simple_ts("1/1/2000", "1/1/2010")
        result = ts.first("10d")
        assert len(result) == 10

        result = ts.first("3M")
        expected = ts[:"3/31/2000"]
        tm.assert_series_equal(result, expected)

        result = ts.first("21D")
        expected = ts[:21]
        tm.assert_series_equal(result, expected)

        result = ts[:0].first("3M")
        tm.assert_series_equal(result, ts[:0])

    def test_first_raises(self):
        # GH20725
        ser = pd.Series("a b c".split())
        msg = "'first' only supports a DatetimeIndex index"
        with pytest.raises(TypeError, match=msg):
            ser.first("1D")

    def test_last_subset(self):
        ts = _simple_ts("1/1/2000", "1/1/2010", freq="12h")
        result = ts.last("10d")
        assert len(result) == 20

        ts = _simple_ts("1/1/2000", "1/1/2010")
        result = ts.last("10d")
        assert len(result) == 10

        result = ts.last("21D")
        expected = ts["12/12/2009":]
        tm.assert_series_equal(result, expected)

        result = ts.last("21D")
        expected = ts[-21:]
        tm.assert_series_equal(result, expected)

        result = ts[:0].last("3M")
        tm.assert_series_equal(result, ts[:0])

    def test_last_raises(self):
        # GH20725
        ser = pd.Series("a b c".split())
        msg = "'last' only supports a DatetimeIndex index"
        with pytest.raises(TypeError, match=msg):
            ser.last("1D")

    def test_format_pre_1900_dates(self):
        rng = date_range("1/1/1850", "1/1/1950", freq="A-DEC")
        rng.format()
        ts = Series(1, index=rng)
        repr(ts)

    def test_at_time(self):
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        ts = Series(np.random.randn(len(rng)), index=rng)
        rs = ts.at_time(rng[1])
        assert (rs.index.hour == rng[1].hour).all()
        assert (rs.index.minute == rng[1].minute).all()
        assert (rs.index.second == rng[1].second).all()

        result = ts.at_time("9:30")
        expected = ts.at_time(time(9, 30))
        tm.assert_series_equal(result, expected)

        df = DataFrame(np.random.randn(len(rng), 3), index=rng)

        result = ts[time(9, 30)]
        result_df = df.loc[time(9, 30)]
        expected = ts[(rng.hour == 9) & (rng.minute == 30)]
        exp_df = df[(rng.hour == 9) & (rng.minute == 30)]

        # FIXME: dont leave commented-out
        # expected.index = date_range('1/1/2000', '1/4/2000')

        tm.assert_series_equal(result, expected)
        tm.assert_frame_equal(result_df, exp_df)

        chunk = df.loc["1/4/2000":]
        result = chunk.loc[time(9, 30)]
        expected = result_df[-1:]
        tm.assert_frame_equal(result, expected)

        # midnight, everything
        rng = date_range("1/1/2000", "1/31/2000")
        ts = Series(np.random.randn(len(rng)), index=rng)

        result = ts.at_time(time(0, 0))
        tm.assert_series_equal(result, ts)

        # time doesn't exist
        rng = date_range("1/1/2012", freq="23Min", periods=384)
        ts = Series(np.random.randn(len(rng)), rng)
        rs = ts.at_time("16:00")
        assert len(rs) == 0

    def test_at_time_raises(self):
        # GH20725
        ser = pd.Series("a b c".split())
        msg = "Index must be DatetimeIndex"
        with pytest.raises(TypeError, match=msg):
            ser.at_time("00:00")

    def test_between(self):
        series = Series(date_range("1/1/2000", periods=10))
        left, right = series[[2, 7]]

        result = series.between(left, right)
        expected = (series >= left) & (series <= right)
        tm.assert_series_equal(result, expected)

    def test_between_time(self):
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        ts = Series(np.random.randn(len(rng)), index=rng)
        stime = time(0, 0)
        etime = time(1, 0)

        close_open = product([True, False], [True, False])
        for inc_start, inc_end in close_open:
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
        tm.assert_series_equal(result, expected)

        # across midnight
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        ts = Series(np.random.randn(len(rng)), index=rng)
        stime = time(22, 0)
        etime = time(9, 0)

        close_open = product([True, False], [True, False])
        for inc_start, inc_end in close_open:
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
        ser = pd.Series("a b c".split())
        msg = "Index must be DatetimeIndex"
        with pytest.raises(TypeError, match=msg):
            ser.between_time(start_time="00:00", end_time="12:00")

    def test_between_time_types(self):
        # GH11818
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        msg = r"Cannot convert arg \[datetime\.datetime\(2010, 1, 2, 1, 0\)\] to a time"
        with pytest.raises(ValueError, match=msg):
            rng.indexer_between_time(datetime(2010, 1, 2, 1), datetime(2010, 1, 2, 5))

        frame = DataFrame({"A": 0}, index=rng)
        with pytest.raises(ValueError, match=msg):
            frame.between_time(datetime(2010, 1, 2, 1), datetime(2010, 1, 2, 5))

        series = Series(0, index=rng)
        with pytest.raises(ValueError, match=msg):
            series.between_time(datetime(2010, 1, 2, 1), datetime(2010, 1, 2, 5))

    @td.skip_if_has_locale
    def test_between_time_formats(self):
        # GH11818
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)

        strings = [
            ("2:00", "2:30"),
            ("0200", "0230"),
            ("2:00am", "2:30am"),
            ("0200am", "0230am"),
            ("2:00:00", "2:30:00"),
            ("020000", "023000"),
            ("2:00:00am", "2:30:00am"),
            ("020000am", "023000am"),
        ]
        expected_length = 28

        for time_string in strings:
            assert len(ts.between_time(*time_string)) == expected_length

    def test_between_time_axis(self):
        # issue 8839
        rng = date_range("1/1/2000", periods=100, freq="10min")
        ts = Series(np.random.randn(len(rng)), index=rng)
        stime, etime = ("08:00:00", "09:00:00")
        expected_length = 7

        assert len(ts.between_time(stime, etime)) == expected_length
        assert len(ts.between_time(stime, etime, axis=0)) == expected_length
        msg = "No axis named 1 for object type <class 'pandas.core.series.Series'>"
        with pytest.raises(ValueError, match=msg):
            ts.between_time(stime, etime, axis=1)

    def test_to_period(self):
        from pandas.core.indexes.period import period_range

        ts = _simple_ts("1/1/2000", "1/1/2001")

        pts = ts.to_period()
        exp = ts.copy()
        exp.index = period_range("1/1/2000", "1/1/2001")
        tm.assert_series_equal(pts, exp)

        pts = ts.to_period("M")
        exp.index = exp.index.asfreq("M")
        tm.assert_index_equal(pts.index, exp.index.asfreq("M"))
        tm.assert_series_equal(pts, exp)

        # GH 7606 without freq
        idx = DatetimeIndex(["2011-01-01", "2011-01-02", "2011-01-03", "2011-01-04"])
        exp_idx = pd.PeriodIndex(
            ["2011-01-01", "2011-01-02", "2011-01-03", "2011-01-04"], freq="D"
        )

        s = Series(np.random.randn(4), index=idx)
        expected = s.copy()
        expected.index = exp_idx
        tm.assert_series_equal(s.to_period(), expected)

        df = DataFrame(np.random.randn(4, 4), index=idx, columns=idx)
        expected = df.copy()
        expected.index = exp_idx
        tm.assert_frame_equal(df.to_period(), expected)

        expected = df.copy()
        expected.columns = exp_idx
        tm.assert_frame_equal(df.to_period(axis=1), expected)

    def test_groupby_count_dateparseerror(self):
        dr = date_range(start="1/1/2012", freq="5min", periods=10)

        # BAD Example, datetimes first
        s = Series(np.arange(10), index=[dr, np.arange(10)])
        grouped = s.groupby(lambda x: x[1] % 2 == 0)
        result = grouped.count()

        s = Series(np.arange(10), index=[np.arange(10), dr])
        grouped = s.groupby(lambda x: x[0] % 2 == 0)
        expected = grouped.count()

        tm.assert_series_equal(result, expected)

    def test_to_csv_numpy_16_bug(self):
        frame = DataFrame({"a": date_range("1/1/2000", periods=10)})

        buf = StringIO()
        frame.to_csv(buf)

        result = buf.getvalue()
        assert "2000-01-01" in result

    def test_series_map_box_timedelta(self):
        # GH 11349
        s = Series(timedelta_range("1 day 1 s", periods=5, freq="h"))

        def f(x):
            return x.total_seconds()

        s.map(f)
        s.apply(f)
        DataFrame(s).applymap(f)

    def test_asfreq_resample_set_correct_freq(self):
        # GH5613
        # we test if .asfreq() and .resample() set the correct value for .freq
        df = pd.DataFrame(
            {"date": ["2012-01-01", "2012-01-02", "2012-01-03"], "col": [1, 2, 3]}
        )
        df = df.set_index(pd.to_datetime(df.date))

        # testing the settings before calling .asfreq() and .resample()
        assert df.index.freq is None
        assert df.index.inferred_freq == "D"

        # does .asfreq() set .freq correctly?
        assert df.asfreq("D").index.freq == "D"

        # does .resample() set .freq correctly?
        assert df.resample("D").asfreq().index.freq == "D"

    def test_pickle(self):

        # GH4606
        p = tm.round_trip_pickle(NaT)
        assert p is NaT

        idx = pd.to_datetime(["2013-01-01", NaT, "2014-01-06"])
        idx_p = tm.round_trip_pickle(idx)
        assert idx_p[0] == idx[0]
        assert idx_p[1] is NaT
        assert idx_p[2] == idx[2]

        # GH11002
        # don't infer freq
        idx = date_range("1750-1-1", "2050-1-1", freq="7D")
        idx_p = tm.round_trip_pickle(idx)
        tm.assert_index_equal(idx, idx_p)

    @pytest.mark.parametrize("tz", [None, "Asia/Tokyo", "US/Eastern"])
    def test_setops_preserve_freq(self, tz):
        rng = date_range("1/1/2000", "1/1/2002", name="idx", tz=tz)

        result = rng[:50].union(rng[50:100])
        assert result.name == rng.name
        assert result.freq == rng.freq
        assert result.tz == rng.tz

        result = rng[:50].union(rng[30:100])
        assert result.name == rng.name
        assert result.freq == rng.freq
        assert result.tz == rng.tz

        result = rng[:50].union(rng[60:100])
        assert result.name == rng.name
        assert result.freq is None
        assert result.tz == rng.tz

        result = rng[:50].intersection(rng[25:75])
        assert result.name == rng.name
        assert result.freqstr == "D"
        assert result.tz == rng.tz

        nofreq = DatetimeIndex(list(rng[25:75]), name="other")
        result = rng[:50].union(nofreq)
        assert result.name is None
        assert result.freq == rng.freq
        assert result.tz == rng.tz

        result = rng[:50].intersection(nofreq)
        assert result.name is None
        assert result.freq == rng.freq
        assert result.tz == rng.tz

    def test_from_M8_structured(self):
        dates = [(datetime(2012, 9, 9, 0, 0), datetime(2012, 9, 8, 15, 10))]
        arr = np.array(dates, dtype=[("Date", "M8[us]"), ("Forecasting", "M8[us]")])
        df = DataFrame(arr)

        assert df["Date"][0] == dates[0][0]
        assert df["Forecasting"][0] == dates[0][1]

        s = Series(arr["Date"])
        assert isinstance(s[0], Timestamp)
        assert s[0] == dates[0][0]

    def test_get_level_values_box(self):
        from pandas import MultiIndex

        dates = date_range("1/1/2000", periods=4)
        levels = [dates, [0, 1]]
        codes = [[0, 0, 1, 1, 2, 2, 3, 3], [0, 1, 0, 1, 0, 1, 0, 1]]

        index = MultiIndex(levels=levels, codes=codes)

        assert isinstance(index.get_level_values(0)[0], Timestamp)

    def test_view_tz(self):
        # GH#24024
        ser = pd.Series(pd.date_range("2000", periods=4, tz="US/Central"))
        result = ser.view("i8")
        expected = pd.Series(
            [
                946706400000000000,
                946792800000000000,
                946879200000000000,
                946965600000000000,
            ]
        )
        tm.assert_series_equal(result, expected)

    def test_asarray_tz_naive(self):
        # This shouldn't produce a warning.
        ser = pd.Series(pd.date_range("2000", periods=2))
        expected = np.array(["2000-01-01", "2000-01-02"], dtype="M8[ns]")
        result = np.asarray(ser)

        tm.assert_numpy_array_equal(result, expected)

        # optionally, object
        result = np.asarray(ser, dtype=object)

        expected = np.array([pd.Timestamp("2000-01-01"), pd.Timestamp("2000-01-02")])
        tm.assert_numpy_array_equal(result, expected)

    def test_asarray_tz_aware(self):
        tz = "US/Central"
        ser = pd.Series(pd.date_range("2000", periods=2, tz=tz))
        expected = np.array(["2000-01-01T06", "2000-01-02T06"], dtype="M8[ns]")
        result = np.asarray(ser, dtype="datetime64[ns]")

        tm.assert_numpy_array_equal(result, expected)

        # Old behavior with no warning
        result = np.asarray(ser, dtype="M8[ns]")

        tm.assert_numpy_array_equal(result, expected)

        # Future behavior with no warning
        expected = np.array(
            [pd.Timestamp("2000-01-01", tz=tz), pd.Timestamp("2000-01-02", tz=tz)]
        )
        result = np.asarray(ser, dtype=object)

        tm.assert_numpy_array_equal(result, expected)
