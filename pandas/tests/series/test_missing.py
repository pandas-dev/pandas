from datetime import timedelta

import numpy as np
import pytest

from pandas._libs import iNaT

import pandas as pd
from pandas import (
    Categorical,
    DataFrame,
    Index,
    IntervalIndex,
    NaT,
    Series,
    Timestamp,
    date_range,
    isna,
)
import pandas._testing as tm


class TestSeriesMissingData:
    def test_categorical_nan_equality(self):
        cat = Series(Categorical(["a", "b", "c", np.nan]))
        exp = Series([True, True, True, False])
        res = cat == cat
        tm.assert_series_equal(res, exp)

    def test_categorical_nan_handling(self):

        # NaNs are represented as -1 in labels
        s = Series(Categorical(["a", "b", np.nan, "a"]))
        tm.assert_index_equal(s.cat.categories, Index(["a", "b"]))
        tm.assert_numpy_array_equal(
            s.values.codes, np.array([0, 1, -1, 0], dtype=np.int8)
        )

    def test_fillna_nat(self):
        series = Series([0, 1, 2, iNaT], dtype="M8[ns]")

        filled = series.fillna(method="pad")
        filled2 = series.fillna(value=series.values[2])

        expected = series.copy()
        expected.values[3] = expected.values[2]

        tm.assert_series_equal(filled, expected)
        tm.assert_series_equal(filled2, expected)

        df = DataFrame({"A": series})
        filled = df.fillna(method="pad")
        filled2 = df.fillna(value=series.values[2])
        expected = DataFrame({"A": expected})
        tm.assert_frame_equal(filled, expected)
        tm.assert_frame_equal(filled2, expected)

        series = Series([iNaT, 0, 1, 2], dtype="M8[ns]")

        filled = series.fillna(method="bfill")
        filled2 = series.fillna(value=series[1])

        expected = series.copy()
        expected[0] = expected[1]

        tm.assert_series_equal(filled, expected)
        tm.assert_series_equal(filled2, expected)

        df = DataFrame({"A": series})
        filled = df.fillna(method="bfill")
        filled2 = df.fillna(value=series[1])
        expected = DataFrame({"A": expected})
        tm.assert_frame_equal(filled, expected)
        tm.assert_frame_equal(filled2, expected)

    def test_isna_for_inf(self):
        s = Series(["a", np.inf, np.nan, pd.NA, 1.0])
        with pd.option_context("mode.use_inf_as_na", True):
            r = s.isna()
            dr = s.dropna()
        e = Series([False, True, True, True, False])
        de = Series(["a", 1.0], index=[0, 4])
        tm.assert_series_equal(r, e)
        tm.assert_series_equal(dr, de)

    def test_isnull_for_inf_deprecated(self):
        # gh-17115
        s = Series(["a", np.inf, np.nan, 1.0])
        with pd.option_context("mode.use_inf_as_null", True):
            r = s.isna()
            dr = s.dropna()

        e = Series([False, True, True, False])
        de = Series(["a", 1.0], index=[0, 3])
        tm.assert_series_equal(r, e)
        tm.assert_series_equal(dr, de)

    def test_timedelta64_nan(self):

        td = Series([timedelta(days=i) for i in range(10)])

        # nan ops on timedeltas
        td1 = td.copy()
        td1[0] = np.nan
        assert isna(td1[0])
        assert td1[0].value == iNaT
        td1[0] = td[0]
        assert not isna(td1[0])

        # GH#16674 iNaT is treated as an integer when given by the user
        td1[1] = iNaT
        assert not isna(td1[1])
        assert td1.dtype == np.object_
        assert td1[1] == iNaT
        td1[1] = td[1]
        assert not isna(td1[1])

        td1[2] = NaT
        assert isna(td1[2])
        assert td1[2].value == iNaT
        td1[2] = td[2]
        assert not isna(td1[2])

        # FIXME: don't leave commented-out
        # boolean setting
        # this doesn't work, not sure numpy even supports it
        # result = td[(td>np.timedelta64(timedelta(days=3))) &
        # td<np.timedelta64(timedelta(days=7)))] = np.nan
        # assert isna(result).sum() == 7

        # NumPy limitation =(

        # def test_logical_range_select(self):
        #     np.random.seed(12345)
        #     selector = -0.5 <= datetime_series <= 0.5
        #     expected = (datetime_series >= -0.5) & (datetime_series <= 0.5)
        #     tm.assert_series_equal(selector, expected)

    def test_dropna_empty(self):
        s = Series([], dtype=object)

        assert len(s.dropna()) == 0
        return_value = s.dropna(inplace=True)
        assert return_value is None
        assert len(s) == 0

        # invalid axis
        msg = "No axis named 1 for object type Series"
        with pytest.raises(ValueError, match=msg):
            s.dropna(axis=1)

    def test_datetime64_tz_dropna(self):
        # DatetimeBlock
        s = Series(
            [
                Timestamp("2011-01-01 10:00"),
                pd.NaT,
                Timestamp("2011-01-03 10:00"),
                pd.NaT,
            ]
        )
        result = s.dropna()
        expected = Series(
            [Timestamp("2011-01-01 10:00"), Timestamp("2011-01-03 10:00")], index=[0, 2]
        )
        tm.assert_series_equal(result, expected)

        # DatetimeBlockTZ
        idx = pd.DatetimeIndex(
            ["2011-01-01 10:00", pd.NaT, "2011-01-03 10:00", pd.NaT], tz="Asia/Tokyo"
        )
        s = Series(idx)
        assert s.dtype == "datetime64[ns, Asia/Tokyo]"
        result = s.dropna()
        expected = Series(
            [
                Timestamp("2011-01-01 10:00", tz="Asia/Tokyo"),
                Timestamp("2011-01-03 10:00", tz="Asia/Tokyo"),
            ],
            index=[0, 2],
        )
        assert result.dtype == "datetime64[ns, Asia/Tokyo]"
        tm.assert_series_equal(result, expected)

    def test_dropna_no_nan(self):
        for s in [Series([1, 2, 3], name="x"), Series([False, True, False], name="x")]:

            result = s.dropna()
            tm.assert_series_equal(result, s)
            assert result is not s

            s2 = s.copy()
            return_value = s2.dropna(inplace=True)
            assert return_value is None
            tm.assert_series_equal(s2, s)

    def test_dropna_intervals(self):
        s = Series(
            [np.nan, 1, 2, 3],
            IntervalIndex.from_arrays([np.nan, 0, 1, 2], [np.nan, 1, 2, 3]),
        )

        result = s.dropna()
        expected = s.iloc[1:]
        tm.assert_series_equal(result, expected)

    def test_valid(self, datetime_series):
        ts = datetime_series.copy()
        ts.index = ts.index._with_freq(None)
        ts[::2] = np.NaN

        result = ts.dropna()
        assert len(result) == ts.count()
        tm.assert_series_equal(result, ts[1::2])
        tm.assert_series_equal(result, ts[pd.notna(ts)])

    def test_isna(self):
        ser = Series([0, 5.4, 3, np.nan, -0.001])
        expected = Series([False, False, False, True, False])
        tm.assert_series_equal(ser.isna(), expected)

        ser = Series(["hi", "", np.nan])
        expected = Series([False, False, True])
        tm.assert_series_equal(ser.isna(), expected)

    def test_notna(self):
        ser = Series([0, 5.4, 3, np.nan, -0.001])
        expected = Series([True, True, True, False, True])
        tm.assert_series_equal(ser.notna(), expected)

        ser = Series(["hi", "", np.nan])
        expected = Series([True, True, False])
        tm.assert_series_equal(ser.notna(), expected)

    def test_pad_require_monotonicity(self):
        rng = date_range("1/1/2000", "3/1/2000", freq="B")

        # neither monotonic increasing or decreasing
        rng2 = rng[[1, 0, 2]]

        msg = "index must be monotonic increasing or decreasing"
        with pytest.raises(ValueError, match=msg):
            rng2.get_indexer(rng, method="pad")

    def test_dropna_preserve_name(self, datetime_series):
        datetime_series[:5] = np.nan
        result = datetime_series.dropna()
        assert result.name == datetime_series.name
        name = datetime_series.name
        ts = datetime_series.copy()
        return_value = ts.dropna(inplace=True)
        assert return_value is None
        assert ts.name == name
