import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Period, Series, period_range
import pandas._testing as tm
from pandas.core.arrays import PeriodArray


class TestSeriesPeriod:
    def setup_method(self, method):
        self.series = Series(period_range("2000-01-01", periods=10, freq="D"))

    def test_auto_conversion(self):
        series = Series(list(period_range("2000-01-01", periods=10, freq="D")))
        assert series.dtype == "Period[D]"

        series = pd.Series(
            [pd.Period("2011-01-01", freq="D"), pd.Period("2011-02-01", freq="D")]
        )
        assert series.dtype == "Period[D]"

    def test_getitem(self):
        assert self.series[1] == pd.Period("2000-01-02", freq="D")

        result = self.series[[2, 4]]
        exp = pd.Series(
            [pd.Period("2000-01-03", freq="D"), pd.Period("2000-01-05", freq="D")],
            index=[2, 4],
            dtype="Period[D]",
        )
        tm.assert_series_equal(result, exp)
        assert result.dtype == "Period[D]"

    def test_isna(self):
        # GH 13737
        s = Series([pd.Period("2011-01", freq="M"), pd.Period("NaT", freq="M")])
        tm.assert_series_equal(s.isna(), Series([False, True]))
        tm.assert_series_equal(s.notna(), Series([True, False]))

    def test_fillna(self):
        # GH 13737
        s = Series([pd.Period("2011-01", freq="M"), pd.Period("NaT", freq="M")])

        res = s.fillna(pd.Period("2012-01", freq="M"))
        exp = Series([pd.Period("2011-01", freq="M"), pd.Period("2012-01", freq="M")])
        tm.assert_series_equal(res, exp)
        assert res.dtype == "Period[M]"

    def test_dropna(self):
        # GH 13737
        s = Series([pd.Period("2011-01", freq="M"), pd.Period("NaT", freq="M")])
        tm.assert_series_equal(s.dropna(), Series([pd.Period("2011-01", freq="M")]))

    def test_between(self):
        left, right = self.series[[2, 7]]
        result = self.series.between(left, right)
        expected = (self.series >= left) & (self.series <= right)
        tm.assert_series_equal(result, expected)

    # ---------------------------------------------------------------------
    # NaT support

    @pytest.mark.xfail(reason="PeriodDtype Series not supported yet")
    def test_NaT_scalar(self):
        series = Series([0, 1000, 2000, pd._libs.iNaT], dtype="period[D]")

        val = series[3]
        assert pd.isna(val)

        series[2] = val
        assert pd.isna(series[2])

    def test_NaT_cast(self):
        result = Series([np.nan]).astype("period[D]")
        expected = Series([pd.NaT], dtype="period[D]")
        tm.assert_series_equal(result, expected)

    def test_set_none(self):
        self.series[3] = None
        assert self.series[3] is pd.NaT

        self.series[3:5] = None
        assert self.series[4] is pd.NaT

    def test_set_nan(self):
        # Do we want to allow this?
        self.series[5] = np.nan
        assert self.series[5] is pd.NaT

        self.series[5:7] = np.nan
        assert self.series[6] is pd.NaT

    def test_intercept_astype_object(self):
        expected = self.series.astype("object")

        df = DataFrame({"a": self.series, "b": np.random.randn(len(self.series))})

        result = df.values.squeeze()
        assert (result[:, 0] == expected.values).all()

        df = DataFrame({"a": self.series, "b": ["foo"] * len(self.series)})

        result = df.values.squeeze()
        assert (result[:, 0] == expected.values).all()

    def test_align_series(self, join_type):
        rng = period_range("1/1/2000", "1/1/2010", freq="A")
        ts = Series(np.random.randn(len(rng)), index=rng)

        ts.align(ts[::2], join=join_type)

    def test_truncate(self):
        # GH 17717
        idx1 = pd.PeriodIndex(
            [pd.Period("2017-09-02"), pd.Period("2017-09-02"), pd.Period("2017-09-03")]
        )
        series1 = pd.Series([1, 2, 3], index=idx1)
        result1 = series1.truncate(after="2017-09-02")

        expected_idx1 = pd.PeriodIndex(
            [pd.Period("2017-09-02"), pd.Period("2017-09-02")]
        )
        tm.assert_series_equal(result1, pd.Series([1, 2], index=expected_idx1))

        idx2 = pd.PeriodIndex(
            [pd.Period("2017-09-03"), pd.Period("2017-09-02"), pd.Period("2017-09-03")]
        )
        series2 = pd.Series([1, 2, 3], index=idx2)
        result2 = series2.sort_index().truncate(after="2017-09-02")

        expected_idx2 = pd.PeriodIndex([pd.Period("2017-09-02")])
        tm.assert_series_equal(result2, pd.Series([2], index=expected_idx2))

    @pytest.mark.parametrize(
        "input_vals",
        [
            [Period("2016-01", freq="M"), Period("2016-02", freq="M")],
            [Period("2016-01-01", freq="D"), Period("2016-01-02", freq="D")],
            [
                Period("2016-01-01 00:00:00", freq="H"),
                Period("2016-01-01 01:00:00", freq="H"),
            ],
            [
                Period("2016-01-01 00:00:00", freq="M"),
                Period("2016-01-01 00:01:00", freq="M"),
            ],
            [
                Period("2016-01-01 00:00:00", freq="S"),
                Period("2016-01-01 00:00:01", freq="S"),
            ],
        ],
    )
    def test_end_time_timevalues(self, input_vals):
        # GH 17157
        # Check that the time part of the Period is adjusted by end_time
        # when using the dt accessor on a Series
        input_vals = PeriodArray._from_sequence(np.asarray(input_vals))

        s = Series(input_vals)
        result = s.dt.end_time
        expected = s.apply(lambda x: x.end_time)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("input_vals", [("2001"), ("NaT")])
    def test_to_period(self, input_vals):
        # GH 21205
        expected = Series([input_vals], dtype="Period[D]")
        result = Series([input_vals], dtype="datetime64[ns]").dt.to_period("D")
        tm.assert_series_equal(result, expected)
