from datetime import datetime, timedelta

import numpy as np
import pytest
import pytz

from pandas._libs import iNaT

import pandas as pd
from pandas import (
    Categorical,
    DataFrame,
    Index,
    IntervalIndex,
    NaT,
    Series,
    Timedelta,
    Timestamp,
    date_range,
    isna,
)
import pandas._testing as tm


class TestSeriesMissingData:
    def test_timedelta_fillna(self):
        # GH 3371
        s = Series(
            [
                Timestamp("20130101"),
                Timestamp("20130101"),
                Timestamp("20130102"),
                Timestamp("20130103 9:01:01"),
            ]
        )
        td = s.diff()

        # reg fillna
        result = td.fillna(Timedelta(seconds=0))
        expected = Series(
            [
                timedelta(0),
                timedelta(0),
                timedelta(1),
                timedelta(days=1, seconds=9 * 3600 + 60 + 1),
            ]
        )
        tm.assert_series_equal(result, expected)

        # interpreted as seconds, deprecated
        with pytest.raises(TypeError, match="Passing integers to fillna"):
            td.fillna(1)

        result = td.fillna(Timedelta(seconds=1))
        expected = Series(
            [
                timedelta(seconds=1),
                timedelta(0),
                timedelta(1),
                timedelta(days=1, seconds=9 * 3600 + 60 + 1),
            ]
        )
        tm.assert_series_equal(result, expected)

        result = td.fillna(timedelta(days=1, seconds=1))
        expected = Series(
            [
                timedelta(days=1, seconds=1),
                timedelta(0),
                timedelta(1),
                timedelta(days=1, seconds=9 * 3600 + 60 + 1),
            ]
        )
        tm.assert_series_equal(result, expected)

        result = td.fillna(np.timedelta64(int(1e9)))
        expected = Series(
            [
                timedelta(seconds=1),
                timedelta(0),
                timedelta(1),
                timedelta(days=1, seconds=9 * 3600 + 60 + 1),
            ]
        )
        tm.assert_series_equal(result, expected)

        result = td.fillna(NaT)
        expected = Series(
            [
                NaT,
                timedelta(0),
                timedelta(1),
                timedelta(days=1, seconds=9 * 3600 + 60 + 1),
            ],
            dtype="m8[ns]",
        )
        tm.assert_series_equal(result, expected)

        # ffill
        td[2] = np.nan
        result = td.ffill()
        expected = td.fillna(Timedelta(seconds=0))
        expected[0] = np.nan
        tm.assert_series_equal(result, expected)

        # bfill
        td[2] = np.nan
        result = td.bfill()
        expected = td.fillna(Timedelta(seconds=0))
        expected[2] = timedelta(days=1, seconds=9 * 3600 + 60 + 1)
        tm.assert_series_equal(result, expected)

    def test_datetime64_fillna(self):

        s = Series(
            [
                Timestamp("20130101"),
                Timestamp("20130101"),
                Timestamp("20130102"),
                Timestamp("20130103 9:01:01"),
            ]
        )
        s[2] = np.nan

        # ffill
        result = s.ffill()
        expected = Series(
            [
                Timestamp("20130101"),
                Timestamp("20130101"),
                Timestamp("20130101"),
                Timestamp("20130103 9:01:01"),
            ]
        )
        tm.assert_series_equal(result, expected)

        # bfill
        result = s.bfill()
        expected = Series(
            [
                Timestamp("20130101"),
                Timestamp("20130101"),
                Timestamp("20130103 9:01:01"),
                Timestamp("20130103 9:01:01"),
            ]
        )
        tm.assert_series_equal(result, expected)

        # GH 6587
        # make sure that we are treating as integer when filling
        # this also tests inference of a datetime-like with NaT's
        s = Series([pd.NaT, pd.NaT, "2013-08-05 15:30:00.000001"])
        expected = Series(
            [
                "2013-08-05 15:30:00.000001",
                "2013-08-05 15:30:00.000001",
                "2013-08-05 15:30:00.000001",
            ],
            dtype="M8[ns]",
        )
        result = s.fillna(method="backfill")
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("tz", ["US/Eastern", "Asia/Tokyo"])
    def test_datetime64_tz_fillna(self, tz):
        # DatetimeBlock
        s = Series(
            [
                Timestamp("2011-01-01 10:00"),
                pd.NaT,
                Timestamp("2011-01-03 10:00"),
                pd.NaT,
            ]
        )
        null_loc = pd.Series([False, True, False, True])

        result = s.fillna(pd.Timestamp("2011-01-02 10:00"))
        expected = Series(
            [
                Timestamp("2011-01-01 10:00"),
                Timestamp("2011-01-02 10:00"),
                Timestamp("2011-01-03 10:00"),
                Timestamp("2011-01-02 10:00"),
            ]
        )
        tm.assert_series_equal(expected, result)
        # check s is not changed
        tm.assert_series_equal(pd.isna(s), null_loc)

        result = s.fillna(pd.Timestamp("2011-01-02 10:00", tz=tz))
        expected = Series(
            [
                Timestamp("2011-01-01 10:00"),
                Timestamp("2011-01-02 10:00", tz=tz),
                Timestamp("2011-01-03 10:00"),
                Timestamp("2011-01-02 10:00", tz=tz),
            ]
        )
        tm.assert_series_equal(expected, result)
        tm.assert_series_equal(pd.isna(s), null_loc)

        result = s.fillna("AAA")
        expected = Series(
            [
                Timestamp("2011-01-01 10:00"),
                "AAA",
                Timestamp("2011-01-03 10:00"),
                "AAA",
            ],
            dtype=object,
        )
        tm.assert_series_equal(expected, result)
        tm.assert_series_equal(pd.isna(s), null_loc)

        result = s.fillna(
            {
                1: pd.Timestamp("2011-01-02 10:00", tz=tz),
                3: pd.Timestamp("2011-01-04 10:00"),
            }
        )
        expected = Series(
            [
                Timestamp("2011-01-01 10:00"),
                Timestamp("2011-01-02 10:00", tz=tz),
                Timestamp("2011-01-03 10:00"),
                Timestamp("2011-01-04 10:00"),
            ]
        )
        tm.assert_series_equal(expected, result)
        tm.assert_series_equal(pd.isna(s), null_loc)

        result = s.fillna(
            {1: pd.Timestamp("2011-01-02 10:00"), 3: pd.Timestamp("2011-01-04 10:00")}
        )
        expected = Series(
            [
                Timestamp("2011-01-01 10:00"),
                Timestamp("2011-01-02 10:00"),
                Timestamp("2011-01-03 10:00"),
                Timestamp("2011-01-04 10:00"),
            ]
        )
        tm.assert_series_equal(expected, result)
        tm.assert_series_equal(pd.isna(s), null_loc)

        # DatetimeBlockTZ
        idx = pd.DatetimeIndex(
            ["2011-01-01 10:00", pd.NaT, "2011-01-03 10:00", pd.NaT], tz=tz
        )
        s = pd.Series(idx)
        assert s.dtype == f"datetime64[ns, {tz}]"
        tm.assert_series_equal(pd.isna(s), null_loc)

        result = s.fillna(pd.Timestamp("2011-01-02 10:00"))
        expected = Series(
            [
                Timestamp("2011-01-01 10:00", tz=tz),
                Timestamp("2011-01-02 10:00"),
                Timestamp("2011-01-03 10:00", tz=tz),
                Timestamp("2011-01-02 10:00"),
            ]
        )
        tm.assert_series_equal(expected, result)
        tm.assert_series_equal(pd.isna(s), null_loc)

        result = s.fillna(pd.Timestamp("2011-01-02 10:00", tz=tz))
        idx = pd.DatetimeIndex(
            [
                "2011-01-01 10:00",
                "2011-01-02 10:00",
                "2011-01-03 10:00",
                "2011-01-02 10:00",
            ],
            tz=tz,
        )
        expected = Series(idx)
        tm.assert_series_equal(expected, result)
        tm.assert_series_equal(pd.isna(s), null_loc)

        result = s.fillna(pd.Timestamp("2011-01-02 10:00", tz=tz).to_pydatetime())
        idx = pd.DatetimeIndex(
            [
                "2011-01-01 10:00",
                "2011-01-02 10:00",
                "2011-01-03 10:00",
                "2011-01-02 10:00",
            ],
            tz=tz,
        )
        expected = Series(idx)
        tm.assert_series_equal(expected, result)
        tm.assert_series_equal(pd.isna(s), null_loc)

        result = s.fillna("AAA")
        expected = Series(
            [
                Timestamp("2011-01-01 10:00", tz=tz),
                "AAA",
                Timestamp("2011-01-03 10:00", tz=tz),
                "AAA",
            ],
            dtype=object,
        )
        tm.assert_series_equal(expected, result)
        tm.assert_series_equal(pd.isna(s), null_loc)

        result = s.fillna(
            {
                1: pd.Timestamp("2011-01-02 10:00", tz=tz),
                3: pd.Timestamp("2011-01-04 10:00"),
            }
        )
        expected = Series(
            [
                Timestamp("2011-01-01 10:00", tz=tz),
                Timestamp("2011-01-02 10:00", tz=tz),
                Timestamp("2011-01-03 10:00", tz=tz),
                Timestamp("2011-01-04 10:00"),
            ]
        )
        tm.assert_series_equal(expected, result)
        tm.assert_series_equal(pd.isna(s), null_loc)

        result = s.fillna(
            {
                1: pd.Timestamp("2011-01-02 10:00", tz=tz),
                3: pd.Timestamp("2011-01-04 10:00", tz=tz),
            }
        )
        expected = Series(
            [
                Timestamp("2011-01-01 10:00", tz=tz),
                Timestamp("2011-01-02 10:00", tz=tz),
                Timestamp("2011-01-03 10:00", tz=tz),
                Timestamp("2011-01-04 10:00", tz=tz),
            ]
        )
        tm.assert_series_equal(expected, result)
        tm.assert_series_equal(pd.isna(s), null_loc)

        # filling with a naive/other zone, coerce to object
        result = s.fillna(Timestamp("20130101"))
        expected = Series(
            [
                Timestamp("2011-01-01 10:00", tz=tz),
                Timestamp("2013-01-01"),
                Timestamp("2011-01-03 10:00", tz=tz),
                Timestamp("2013-01-01"),
            ]
        )
        tm.assert_series_equal(expected, result)
        tm.assert_series_equal(pd.isna(s), null_loc)

        result = s.fillna(Timestamp("20130101", tz="US/Pacific"))
        expected = Series(
            [
                Timestamp("2011-01-01 10:00", tz=tz),
                Timestamp("2013-01-01", tz="US/Pacific"),
                Timestamp("2011-01-03 10:00", tz=tz),
                Timestamp("2013-01-01", tz="US/Pacific"),
            ]
        )
        tm.assert_series_equal(expected, result)
        tm.assert_series_equal(pd.isna(s), null_loc)

    def test_fillna_dt64tz_with_method(self):
        # with timezone
        # GH 15855
        ser = pd.Series([pd.Timestamp("2012-11-11 00:00:00+01:00"), pd.NaT])
        exp = pd.Series(
            [
                pd.Timestamp("2012-11-11 00:00:00+01:00"),
                pd.Timestamp("2012-11-11 00:00:00+01:00"),
            ]
        )
        tm.assert_series_equal(ser.fillna(method="pad"), exp)

        ser = pd.Series([pd.NaT, pd.Timestamp("2012-11-11 00:00:00+01:00")])
        exp = pd.Series(
            [
                pd.Timestamp("2012-11-11 00:00:00+01:00"),
                pd.Timestamp("2012-11-11 00:00:00+01:00"),
            ]
        )
        tm.assert_series_equal(ser.fillna(method="bfill"), exp)

    def test_fillna_consistency(self):
        # GH 16402
        # fillna with a tz aware to a tz-naive, should result in object

        s = Series([Timestamp("20130101"), pd.NaT])

        result = s.fillna(Timestamp("20130101", tz="US/Eastern"))
        expected = Series(
            [Timestamp("20130101"), Timestamp("2013-01-01", tz="US/Eastern")],
            dtype="object",
        )
        tm.assert_series_equal(result, expected)

        # where (we ignore the errors=)
        result = s.where(
            [True, False], Timestamp("20130101", tz="US/Eastern"), errors="ignore"
        )
        tm.assert_series_equal(result, expected)

        result = s.where(
            [True, False], Timestamp("20130101", tz="US/Eastern"), errors="ignore"
        )
        tm.assert_series_equal(result, expected)

        # with a non-datetime
        result = s.fillna("foo")
        expected = Series([Timestamp("20130101"), "foo"])
        tm.assert_series_equal(result, expected)

        # assignment
        s2 = s.copy()
        s2[1] = "foo"
        tm.assert_series_equal(s2, expected)

    def test_datetime64tz_fillna_round_issue(self):
        # GH 14872

        data = pd.Series(
            [pd.NaT, pd.NaT, datetime(2016, 12, 12, 22, 24, 6, 100001, tzinfo=pytz.utc)]
        )

        filled = data.fillna(method="bfill")

        expected = pd.Series(
            [
                datetime(2016, 12, 12, 22, 24, 6, 100001, tzinfo=pytz.utc),
                datetime(2016, 12, 12, 22, 24, 6, 100001, tzinfo=pytz.utc),
                datetime(2016, 12, 12, 22, 24, 6, 100001, tzinfo=pytz.utc),
            ]
        )

        tm.assert_series_equal(filled, expected)

    def test_fillna_downcast(self):
        # GH 15277
        # infer int64 from float64
        s = pd.Series([1.0, np.nan])
        result = s.fillna(0, downcast="infer")
        expected = pd.Series([1, 0])
        tm.assert_series_equal(result, expected)

        # infer int64 from float64 when fillna value is a dict
        s = pd.Series([1.0, np.nan])
        result = s.fillna({1: 0}, downcast="infer")
        expected = pd.Series([1, 0])
        tm.assert_series_equal(result, expected)

    def test_fillna_int(self):
        s = Series(np.random.randint(-100, 100, 50))
        return_value = s.fillna(method="ffill", inplace=True)
        assert return_value is None
        tm.assert_series_equal(s.fillna(method="ffill", inplace=False), s)

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

    def test_fillna(self, datetime_series):
        ts = Series([0.0, 1.0, 2.0, 3.0, 4.0], index=tm.makeDateIndex(5))

        tm.assert_series_equal(ts, ts.fillna(method="ffill"))

        ts[2] = np.NaN

        exp = Series([0.0, 1.0, 1.0, 3.0, 4.0], index=ts.index)
        tm.assert_series_equal(ts.fillna(method="ffill"), exp)

        exp = Series([0.0, 1.0, 3.0, 3.0, 4.0], index=ts.index)
        tm.assert_series_equal(ts.fillna(method="backfill"), exp)

        exp = Series([0.0, 1.0, 5.0, 3.0, 4.0], index=ts.index)
        tm.assert_series_equal(ts.fillna(value=5), exp)

        msg = "Must specify a fill 'value' or 'method'"
        with pytest.raises(ValueError, match=msg):
            ts.fillna()

        msg = "Cannot specify both 'value' and 'method'"
        with pytest.raises(ValueError, match=msg):
            datetime_series.fillna(value=0, method="ffill")

        # GH 5703
        s1 = Series([np.nan])
        s2 = Series([1])
        result = s1.fillna(s2)
        expected = Series([1.0])
        tm.assert_series_equal(result, expected)
        result = s1.fillna({})
        tm.assert_series_equal(result, s1)
        result = s1.fillna(Series((), dtype=object))
        tm.assert_series_equal(result, s1)
        result = s2.fillna(s1)
        tm.assert_series_equal(result, s2)
        result = s1.fillna({0: 1})
        tm.assert_series_equal(result, expected)
        result = s1.fillna({1: 1})
        tm.assert_series_equal(result, Series([np.nan]))
        result = s1.fillna({0: 1, 1: 1})
        tm.assert_series_equal(result, expected)
        result = s1.fillna(Series({0: 1, 1: 1}))
        tm.assert_series_equal(result, expected)
        result = s1.fillna(Series({0: 1, 1: 1}, index=[4, 5]))
        tm.assert_series_equal(result, s1)

        s1 = Series([0, 1, 2], list("abc"))
        s2 = Series([0, np.nan, 2], list("bac"))
        result = s2.fillna(s1)
        expected = Series([0, 0, 2.0], list("bac"))
        tm.assert_series_equal(result, expected)

        # limit
        s = Series(np.nan, index=[0, 1, 2])
        result = s.fillna(999, limit=1)
        expected = Series([999, np.nan, np.nan], index=[0, 1, 2])
        tm.assert_series_equal(result, expected)

        result = s.fillna(999, limit=2)
        expected = Series([999, 999, np.nan], index=[0, 1, 2])
        tm.assert_series_equal(result, expected)

        # GH 9043
        # make sure a string representation of int/float values can be filled
        # correctly without raising errors or being converted
        vals = ["0", "1.5", "-0.3"]
        for val in vals:
            s = Series([0, 1, np.nan, np.nan, 4], dtype="float64")
            result = s.fillna(val)
            expected = Series([0, 1, val, val, 4], dtype="object")
            tm.assert_series_equal(result, expected)

    def test_fillna_bug(self):
        x = Series([np.nan, 1.0, np.nan, 3.0, np.nan], ["z", "a", "b", "c", "d"])
        filled = x.fillna(method="ffill")
        expected = Series([np.nan, 1.0, 1.0, 3.0, 3.0], x.index)
        tm.assert_series_equal(filled, expected)

        filled = x.fillna(method="bfill")
        expected = Series([1.0, 1.0, 3.0, 3.0, np.nan], x.index)
        tm.assert_series_equal(filled, expected)

    def test_fillna_invalid_method(self, datetime_series):
        try:
            datetime_series.fillna(method="ffil")
        except ValueError as inst:
            assert "ffil" in str(inst)

    def test_ffill(self):
        ts = Series([0.0, 1.0, 2.0, 3.0, 4.0], index=tm.makeDateIndex(5))
        ts[2] = np.NaN
        tm.assert_series_equal(ts.ffill(), ts.fillna(method="ffill"))

    def test_ffill_mixed_dtypes_without_missing_data(self):
        # GH14956
        series = pd.Series([datetime(2015, 1, 1, tzinfo=pytz.utc), 1])
        result = series.ffill()
        tm.assert_series_equal(series, result)

    def test_bfill(self):
        ts = Series([0.0, 1.0, 2.0, 3.0, 4.0], index=tm.makeDateIndex(5))
        ts[2] = np.NaN
        tm.assert_series_equal(ts.bfill(), ts.fillna(method="bfill"))

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
        s = pd.Series(idx)
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

    def test_pad_nan(self):
        x = Series(
            [np.nan, 1.0, np.nan, 3.0, np.nan], ["z", "a", "b", "c", "d"], dtype=float
        )

        return_value = x.fillna(method="pad", inplace=True)
        assert return_value is None

        expected = Series(
            [np.nan, 1.0, 1.0, 3.0, 3.0], ["z", "a", "b", "c", "d"], dtype=float
        )
        tm.assert_series_equal(x[1:], expected[1:])
        assert np.isnan(x[0]), np.isnan(expected[0])

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

    def test_series_fillna_limit(self):
        index = np.arange(10)
        s = Series(np.random.randn(10), index=index)

        result = s[:2].reindex(index)
        result = result.fillna(method="pad", limit=5)

        expected = s[:2].reindex(index).fillna(method="pad")
        expected[-3:] = np.nan
        tm.assert_series_equal(result, expected)

        result = s[-2:].reindex(index)
        result = result.fillna(method="bfill", limit=5)

        expected = s[-2:].reindex(index).fillna(method="backfill")
        expected[:3] = np.nan
        tm.assert_series_equal(result, expected)

    def test_series_pad_backfill_limit(self):
        index = np.arange(10)
        s = Series(np.random.randn(10), index=index)

        result = s[:2].reindex(index, method="pad", limit=5)

        expected = s[:2].reindex(index).fillna(method="pad")
        expected[-3:] = np.nan
        tm.assert_series_equal(result, expected)

        result = s[-2:].reindex(index, method="backfill", limit=5)

        expected = s[-2:].reindex(index).fillna(method="backfill")
        expected[:3] = np.nan
        tm.assert_series_equal(result, expected)
