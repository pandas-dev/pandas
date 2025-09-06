from datetime import (
    datetime,
    timezone,
)
import re
import warnings
import zoneinfo

import dateutil
import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    Series,
    Timestamp,
)
import pandas._testing as tm
from pandas.core.indexes.datetimes import date_range
from pandas.core.indexes.period import (
    Period,
    PeriodIndex,
    period_range,
)
from pandas.core.resample import _get_period_range_edges

from pandas.tseries import offsets


@pytest.fixture
def simple_period_range_series():
    """
    Series with period range index and random data for test purposes.
    """

    def _simple_period_range_series(start, end, freq="D"):
        with warnings.catch_warnings():
            # suppress Period[B] deprecation warning
            msg = "|".join(["Period with BDay freq", r"PeriodDtype\[B\] is deprecated"])
            warnings.filterwarnings(
                "ignore",
                msg,
                category=FutureWarning,
            )
            rng = period_range(start, end, freq=freq)
        return Series(np.random.default_rng(2).standard_normal(len(rng)), index=rng)

    return _simple_period_range_series


class TestPeriodIndex:
    @pytest.mark.parametrize("freq", ["2D", "1h", "2h"])
    def test_asfreq(self, frame_or_series, freq):
        # GH 12884, 15944

        obj = frame_or_series(range(5), index=period_range("2020-01-01", periods=5))

        expected = obj.to_timestamp().resample(freq).asfreq()
        result = obj.to_timestamp().resample(freq).asfreq()
        tm.assert_almost_equal(result, expected)

    def test_asfreq_fill_value(self):
        # test for fill value during resampling, issue 3715

        index = period_range(datetime(2005, 1, 1), datetime(2005, 1, 10), freq="D")
        s = Series(range(len(index)), index=index)
        new_index = date_range(
            s.index[0].to_timestamp(how="start"),
            (s.index[-1]).to_timestamp(how="start"),
            freq="1h",
        )
        expected = s.to_timestamp().reindex(new_index, fill_value=4.0)
        result = s.to_timestamp().resample("1h").asfreq(fill_value=4.0)
        tm.assert_series_equal(result, expected)

        frame = s.to_frame("value")
        new_index = date_range(
            frame.index[0].to_timestamp(how="start"),
            (frame.index[-1]).to_timestamp(how="start"),
            freq="1h",
        )
        expected = frame.to_timestamp().reindex(new_index, fill_value=3.0)
        result = frame.to_timestamp().resample("1h").asfreq(fill_value=3.0)
        tm.assert_frame_equal(result, expected)

    def test_resample_basic(self):
        # GH3609
        s = Series(
            range(100),
            index=date_range("20130101", freq="s", periods=100, name="idx"),
            dtype="float",
        )
        s[10:30] = np.nan
        index = PeriodIndex(
            [Period("2013-01-01 00:00", "min"), Period("2013-01-01 00:01", "min")],
            name="idx",
        )
        expected = Series([34.5, 79.5], index=index)
        result = s.resample("min").mean().to_period()
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "tz",
        [
            zoneinfo.ZoneInfo("America/Los_Angeles"),
            dateutil.tz.gettz("America/Los_Angeles"),
        ],
    )
    def test_with_local_timezone(self, tz):
        # see gh-5430
        local_timezone = tz

        start = datetime(
            year=2013, month=11, day=1, hour=0, minute=0, tzinfo=timezone.utc
        )
        # 1 day later
        end = datetime(
            year=2013, month=11, day=2, hour=0, minute=0, tzinfo=timezone.utc
        )

        index = date_range(start, end, freq="h", name="idx")

        series = Series(1, index=index)
        series = series.tz_convert(local_timezone)
        msg = "Converting to PeriodArray/Index representation will drop timezone"
        with tm.assert_produces_warning(UserWarning, match=msg):
            result = series.resample("D").mean().to_period()

        # Create the expected series
        # Index is moved back a day with the timezone conversion from UTC to
        # Pacific
        expected_index = (
            period_range(start=start, end=end, freq="D", name="idx") - offsets.Day()
        )
        expected = Series(1.0, index=expected_index)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "tz",
        [
            zoneinfo.ZoneInfo("America/Los_Angeles"),
            dateutil.tz.gettz("America/Los_Angeles"),
        ],
    )
    def test_resample_with_tz(self, tz, unit):
        # GH 13238
        dti = date_range("2017-01-01", periods=48, freq="h", tz=tz, unit=unit)
        ser = Series(2, index=dti)
        result = ser.resample("D").mean()
        exp_dti = pd.DatetimeIndex(
            ["2017-01-01", "2017-01-02"], tz=tz, freq="D"
        ).as_unit(unit)
        expected = Series(
            2.0,
            index=exp_dti,
        )
        tm.assert_series_equal(result, expected)

    def test_resample_nonexistent_time_bin_edge(self):
        # GH 19375
        index = date_range("2017-03-12", "2017-03-12 1:45:00", freq="15min")
        s = Series(np.zeros(len(index)), index=index)
        expected = s.tz_localize("US/Pacific")
        expected.index = pd.DatetimeIndex(expected.index, freq="900s")
        result = expected.resample("900s").mean()
        tm.assert_series_equal(result, expected)

    def test_resample_nonexistent_time_bin_edge2(self):
        # GH 23742
        index = date_range(start="2017-10-10", end="2017-10-20", freq="1h")
        index = index.tz_localize("UTC").tz_convert("America/Sao_Paulo")
        df = DataFrame(data=list(range(len(index))), index=index)
        result = df.groupby(pd.Grouper(freq="1D")).count()
        expected = date_range(
            start="2017-10-09",
            end="2017-10-20",
            freq="D",
            tz="America/Sao_Paulo",
            nonexistent="shift_forward",
            inclusive="left",
        )
        tm.assert_index_equal(result.index, expected)

    def test_resample_ambiguous_time_bin_edge(self):
        # GH 10117
        idx = date_range(
            "2014-10-25 22:00:00",
            "2014-10-26 00:30:00",
            freq="30min",
            tz="Europe/London",
        )
        expected = Series(np.zeros(len(idx)), index=idx)
        result = expected.resample("30min").mean()
        tm.assert_series_equal(result, expected)

    def test_fill_method_and_how_upsample(self):
        # GH2073
        s = Series(
            np.arange(9, dtype="int64"),
            index=date_range("2010-01-01", periods=9, freq="QE"),
        )
        last = s.resample("ME").ffill()
        both = s.resample("ME").ffill().resample("ME").last().astype("int64")
        tm.assert_series_equal(last, both)

    def test_resample_irregular_sparse(self):
        dr = date_range(start="1/1/2012", freq="5min", periods=1000)
        s = Series(np.array(100), index=dr)
        # subset the data.
        subset = s[:"2012-01-04 06:55"]

        result = subset.resample("10min").apply(len)
        expected = s.resample("10min").apply(len).loc[result.index]
        tm.assert_series_equal(result, expected)

    def test_resample_weekly_all_na(self):
        rng = date_range("1/1/2000", periods=10, freq="W-WED")
        ts = Series(np.random.default_rng(2).standard_normal(len(rng)), index=rng)

        result = ts.resample("W-THU").asfreq()

        assert result.isna().all()

        result = ts.resample("W-THU").asfreq().ffill()[:-1]
        expected = ts.asfreq("W-THU").ffill()
        tm.assert_series_equal(result, expected)

    def test_resample_tz_localized(self, unit):
        dr = date_range(start="2012-4-13", end="2012-5-1", unit=unit)
        ts = Series(range(len(dr)), index=dr)

        ts_utc = ts.tz_localize("UTC")
        ts_local = ts_utc.tz_convert("America/Los_Angeles")

        result = ts_local.resample("W").mean()

        ts_local_naive = ts_local.copy()
        ts_local_naive.index = ts_local_naive.index.tz_localize(None)

        exp = ts_local_naive.resample("W").mean().tz_localize("America/Los_Angeles")
        exp.index = pd.DatetimeIndex(exp.index, freq="W")

        tm.assert_series_equal(result, exp)

        # it works
        result = ts_local.resample("D").mean()

    def test_resample_tz_localized2(self):
        # #2245
        idx = date_range(
            "2001-09-20 15:59", "2001-09-20 16:00", freq="min", tz="Australia/Sydney"
        )
        s = Series([1, 2], index=idx)

        result = s.resample("D", closed="right", label="right").mean()
        ex_index = date_range("2001-09-21", periods=1, freq="D", tz="Australia/Sydney")
        expected = Series([1.5], index=ex_index)

        tm.assert_series_equal(result, expected)

        # for good measure
        msg = "Converting to PeriodArray/Index representation will drop timezone "
        with tm.assert_produces_warning(UserWarning, match=msg):
            result = s.resample("D").mean().to_period()
        ex_index = period_range("2001-09-20", periods=1, freq="D")
        expected = Series([1.5], index=ex_index)
        tm.assert_series_equal(result, expected)

    def test_resample_tz_localized3(self):
        # GH 6397
        # comparing an offset that doesn't propagate tz's
        rng = date_range("1/1/2011", periods=20000, freq="h")
        rng = rng.tz_localize("EST")
        ts = DataFrame(index=rng)
        ts["first"] = np.random.default_rng(2).standard_normal(len(rng))
        ts["second"] = np.cumsum(np.random.default_rng(2).standard_normal(len(rng)))
        expected = DataFrame(
            {
                "first": ts.resample("YE").sum()["first"],
                "second": ts.resample("YE").mean()["second"],
            },
            columns=["first", "second"],
        )
        result = (
            ts.resample("YE")
            .agg({"first": "sum", "second": "mean"})
            .reindex(columns=["first", "second"])
        )
        tm.assert_frame_equal(result, expected)

    def test_closed_left_corner(self):
        # #1465
        s = Series(
            np.random.default_rng(2).standard_normal(21),
            index=date_range(start="1/1/2012 9:30", freq="1min", periods=21),
        )
        s.iloc[0] = np.nan

        result = s.resample("10min", closed="left", label="right").mean()
        exp = s[1:].resample("10min", closed="left", label="right").mean()
        tm.assert_series_equal(result, exp)

        result = s.resample("10min", closed="left", label="left").mean()
        exp = s[1:].resample("10min", closed="left", label="left").mean()

        ex_index = date_range(start="1/1/2012 9:30", freq="10min", periods=3)

        tm.assert_index_equal(result.index, ex_index)
        tm.assert_series_equal(result, exp)

    def test_resample_weekly_bug_1726(self):
        # 8/6/12 is a Monday
        ind = date_range(start="8/6/2012", end="8/26/2012", freq="D")
        n = len(ind)
        data = [[x] * 5 for x in range(n)]
        df = DataFrame(data, columns=["open", "high", "low", "close", "vol"], index=ind)

        # it works!
        df.resample("W-MON", closed="left", label="left").first()

    def test_resample_with_dst_time_change(self):
        # GH 15549
        index = (
            pd.DatetimeIndex([1457537600000000000, 1458059600000000000])
            .tz_localize("UTC")
            .tz_convert("America/Chicago")
        )
        df = DataFrame([1, 2], index=index)
        result = df.resample("12h", closed="right", label="right").last().ffill()

        expected_index_values = [
            "2016-03-09 12:00:00-06:00",
            "2016-03-10 00:00:00-06:00",
            "2016-03-10 12:00:00-06:00",
            "2016-03-11 00:00:00-06:00",
            "2016-03-11 12:00:00-06:00",
            "2016-03-12 00:00:00-06:00",
            "2016-03-12 12:00:00-06:00",
            "2016-03-13 00:00:00-06:00",
            "2016-03-13 13:00:00-05:00",
            "2016-03-14 01:00:00-05:00",
            "2016-03-14 13:00:00-05:00",
            "2016-03-15 01:00:00-05:00",
            "2016-03-15 13:00:00-05:00",
        ]
        index = (
            pd.to_datetime(expected_index_values, utc=True)
            .tz_convert("America/Chicago")
            .as_unit(index.unit)
        )
        index = pd.DatetimeIndex(index, freq="12h")
        expected = DataFrame(
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0],
            index=index,
        )
        tm.assert_frame_equal(result, expected)

    def test_resample_bms_2752(self):
        # GH2753
        timeseries = Series(
            index=pd.bdate_range("20000101", "20000201"), dtype=np.float64
        )
        res1 = timeseries.resample("BMS").mean()
        res2 = timeseries.resample("BMS").mean().resample("B").mean()
        assert res1.index[0] == Timestamp("20000103")
        assert res1.index[0] == res2.index[0]

    @pytest.mark.xfail(reason="Commented out for more than 3 years. Should this work?")
    def test_monthly_convention_span(self):
        rng = period_range("2000-01", periods=3, freq="ME")
        ts = Series(np.arange(3), index=rng)

        # hacky way to get same thing
        exp_index = period_range("2000-01-01", "2000-03-31", freq="D")
        expected = ts.asfreq("D", how="end").reindex(exp_index)
        expected = expected.fillna(method="bfill")

        result = ts.resample("D").mean()

        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "from_freq, to_freq", [("D", "ME"), ("QE", "YE"), ("ME", "QE"), ("D", "W")]
    )
    def test_default_right_closed_label(self, from_freq, to_freq):
        idx = date_range(start="8/15/2012", periods=100, freq=from_freq)
        df = DataFrame(np.random.default_rng(2).standard_normal((len(idx), 2)), idx)

        resampled = df.resample(to_freq).mean()
        tm.assert_frame_equal(
            resampled, df.resample(to_freq, closed="right", label="right").mean()
        )

    @pytest.mark.parametrize(
        "from_freq, to_freq",
        [("D", "MS"), ("QE", "YS"), ("ME", "QS"), ("h", "D"), ("min", "h")],
    )
    def test_default_left_closed_label(self, from_freq, to_freq):
        idx = date_range(start="8/15/2012", periods=100, freq=from_freq)
        df = DataFrame(np.random.default_rng(2).standard_normal((len(idx), 2)), idx)

        resampled = df.resample(to_freq).mean()
        tm.assert_frame_equal(
            resampled, df.resample(to_freq, closed="left", label="left").mean()
        )

    def test_evenly_divisible_with_no_extra_bins(self):
        # GH#4076
        # when the frequency is evenly divisible, sometimes extra bins

        df = DataFrame(
            np.random.default_rng(2).standard_normal((9, 3)),
            index=date_range("2000-1-1", periods=9),
        )
        result = df.resample("5D").mean()
        expected = pd.concat([df.iloc[0:5].mean(), df.iloc[5:].mean()], axis=1).T
        expected.index = pd.DatetimeIndex(
            [Timestamp("2000-1-1"), Timestamp("2000-1-6")], dtype="M8[ns]", freq="5D"
        )
        tm.assert_frame_equal(result, expected)

    def test_evenly_divisible_with_no_extra_bins2(self):
        index = date_range(start="2001-5-4", periods=28)
        df = DataFrame(
            [
                {
                    "REST_KEY": 1,
                    "DLY_TRN_QT": 80,
                    "DLY_SLS_AMT": 90,
                    "COOP_DLY_TRN_QT": 30,
                    "COOP_DLY_SLS_AMT": 20,
                }
            ]
            * 28
            + [
                {
                    "REST_KEY": 2,
                    "DLY_TRN_QT": 70,
                    "DLY_SLS_AMT": 10,
                    "COOP_DLY_TRN_QT": 50,
                    "COOP_DLY_SLS_AMT": 20,
                }
            ]
            * 28,
            index=index.append(index),
        ).sort_index()

        index = date_range("2001-5-4", periods=4, freq="7D")
        expected = DataFrame(
            [
                {
                    "REST_KEY": 14,
                    "DLY_TRN_QT": 14,
                    "DLY_SLS_AMT": 14,
                    "COOP_DLY_TRN_QT": 14,
                    "COOP_DLY_SLS_AMT": 14,
                }
            ]
            * 4,
            index=index,
        )
        result = df.resample("7D").count()
        tm.assert_frame_equal(result, expected)

        expected = DataFrame(
            [
                {
                    "REST_KEY": 21,
                    "DLY_TRN_QT": 1050,
                    "DLY_SLS_AMT": 700,
                    "COOP_DLY_TRN_QT": 560,
                    "COOP_DLY_SLS_AMT": 280,
                }
            ]
            * 4,
            index=index,
        )
        result = df.resample("7D").sum()
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "first,last,freq,freq_to_offset,exp_first,exp_last",
        [
            ("19910905", "19920406", "D", "D", "19910905", "19920406"),
            ("19910905 00:00", "19920406 06:00", "D", "D", "19910905", "19920406"),
            (
                "19910905 06:00",
                "19920406 06:00",
                "h",
                "h",
                "19910905 06:00",
                "19920406 06:00",
            ),
            ("19910906", "19920406", "M", "ME", "1991-09", "1992-04"),
            ("19910831", "19920430", "M", "ME", "1991-08", "1992-04"),
            ("1991-08", "1992-04", "M", "ME", "1991-08", "1992-04"),
        ],
    )
    def test_get_period_range_edges(
        self, first, last, freq, freq_to_offset, exp_first, exp_last
    ):
        first = Period(first)
        last = Period(last)

        exp_first = Period(exp_first, freq=freq)
        exp_last = Period(exp_last, freq=freq)

        freq = pd.tseries.frequencies.to_offset(freq_to_offset)
        result = _get_period_range_edges(first, last, freq)
        expected = (exp_first, exp_last)
        assert result == expected

    def test_resample_t_l_deprecated(self):
        # GH#52536
        msg_t = "Invalid frequency: T"
        msg_l = "Invalid frequency: L"

        with pytest.raises(ValueError, match=msg_l):
            period_range(
                "2020-01-01 00:00:00 00:00", "2020-01-01 00:00:00 00:01", freq="L"
            )
        rng_l = period_range(
            "2020-01-01 00:00:00 00:00", "2020-01-01 00:00:00 00:01", freq="ms"
        )
        ser = Series(np.arange(len(rng_l)), index=rng_l)

        with pytest.raises(ValueError, match=msg_t):
            ser.resample("T").mean()

    @pytest.mark.parametrize(
        "freq, freq_depr, freq_depr_res",
        [
            ("2Q", "2q", "2y"),
            ("2M", "2m", "2q"),
        ],
    )
    def test_resample_lowercase_frequency_raises(self, freq, freq_depr, freq_depr_res):
        msg = f"Invalid frequency: {freq_depr}"
        with pytest.raises(ValueError, match=msg):
            period_range("2020-01-01", "2020-08-01", freq=freq_depr)

        msg = f"Invalid frequency: {freq_depr_res}"
        rng = period_range("2020-01-01", "2020-08-01", freq=freq)
        ser = Series(np.arange(len(rng)), index=rng)
        with pytest.raises(ValueError, match=msg):
            ser.resample(freq_depr_res).mean()

    @pytest.mark.parametrize(
        "offset",
        [
            offsets.MonthBegin(),
            offsets.BYearBegin(2),
            offsets.BusinessHour(2),
        ],
    )
    def test_asfreq_invalid_period_offset(self, offset, frame_or_series):
        # GH#55785
        msg = re.escape(f"{offset} is not supported as period frequency")

        obj = frame_or_series(range(5), index=period_range("2020-01-01", periods=5))
        with pytest.raises(ValueError, match=msg):
            obj.asfreq(freq=offset)


@pytest.mark.parametrize(
    "freq",
    [
        ("2ME"),
        ("2QE"),
        ("2QE-FEB"),
        ("2YE"),
        ("2YE-MAR"),
        ("2me"),
        ("2qe"),
        ("2ye-mar"),
    ],
)
def test_resample_frequency_ME_QE_YE_raises(frame_or_series, freq):
    # GH#9586
    msg = f"{freq[1:]} is not supported as period frequency"

    obj = frame_or_series(range(5), index=period_range("2020-01-01", periods=5))
    msg = f"Invalid frequency: {freq}"
    with pytest.raises(ValueError, match=msg):
        obj.resample(freq)


@pytest.mark.parametrize("freq", ["2BME", "2CBME", "2SME", "2BQE-FEB", "2BYE-MAR"])
def test_resample_frequency_invalid_freq(frame_or_series, freq):
    # GH#9586
    msg = f"Invalid frequency: {freq}"

    obj = frame_or_series(range(5), index=period_range("2020-01-01", periods=5))
    with pytest.raises(ValueError, match=msg):
        obj.resample(freq)
