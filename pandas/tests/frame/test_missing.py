import datetime

import dateutil
import numpy as np
import pytest

import pandas as pd
from pandas import Categorical, DataFrame, Series, Timestamp, date_range
import pandas._testing as tm
from pandas.tests.frame.common import _check_mixed_float


class TestDataFrameMissingData:
    def test_dropEmptyRows(self, float_frame):
        N = len(float_frame.index)
        mat = np.random.randn(N)
        mat[:5] = np.nan

        frame = DataFrame({"foo": mat}, index=float_frame.index)
        original = Series(mat, index=float_frame.index, name="foo")
        expected = original.dropna()
        inplace_frame1, inplace_frame2 = frame.copy(), frame.copy()

        smaller_frame = frame.dropna(how="all")
        # check that original was preserved
        tm.assert_series_equal(frame["foo"], original)
        return_value = inplace_frame1.dropna(how="all", inplace=True)
        tm.assert_series_equal(smaller_frame["foo"], expected)
        tm.assert_series_equal(inplace_frame1["foo"], expected)
        assert return_value is None

        smaller_frame = frame.dropna(how="all", subset=["foo"])
        return_value = inplace_frame2.dropna(how="all", subset=["foo"], inplace=True)
        tm.assert_series_equal(smaller_frame["foo"], expected)
        tm.assert_series_equal(inplace_frame2["foo"], expected)
        assert return_value is None

    def test_dropIncompleteRows(self, float_frame):
        N = len(float_frame.index)
        mat = np.random.randn(N)
        mat[:5] = np.nan

        frame = DataFrame({"foo": mat}, index=float_frame.index)
        frame["bar"] = 5
        original = Series(mat, index=float_frame.index, name="foo")
        inp_frame1, inp_frame2 = frame.copy(), frame.copy()

        smaller_frame = frame.dropna()
        tm.assert_series_equal(frame["foo"], original)
        return_value = inp_frame1.dropna(inplace=True)

        exp = Series(mat[5:], index=float_frame.index[5:], name="foo")
        tm.assert_series_equal(smaller_frame["foo"], exp)
        tm.assert_series_equal(inp_frame1["foo"], exp)
        assert return_value is None

        samesize_frame = frame.dropna(subset=["bar"])
        tm.assert_series_equal(frame["foo"], original)
        assert (frame["bar"] == 5).all()
        return_value = inp_frame2.dropna(subset=["bar"], inplace=True)
        tm.assert_index_equal(samesize_frame.index, float_frame.index)
        tm.assert_index_equal(inp_frame2.index, float_frame.index)
        assert return_value is None

    def test_dropna(self):
        df = DataFrame(np.random.randn(6, 4))
        df[2][:2] = np.nan

        dropped = df.dropna(axis=1)
        expected = df.loc[:, [0, 1, 3]]
        inp = df.copy()
        return_value = inp.dropna(axis=1, inplace=True)
        tm.assert_frame_equal(dropped, expected)
        tm.assert_frame_equal(inp, expected)
        assert return_value is None

        dropped = df.dropna(axis=0)
        expected = df.loc[list(range(2, 6))]
        inp = df.copy()
        return_value = inp.dropna(axis=0, inplace=True)
        tm.assert_frame_equal(dropped, expected)
        tm.assert_frame_equal(inp, expected)
        assert return_value is None

        # threshold
        dropped = df.dropna(axis=1, thresh=5)
        expected = df.loc[:, [0, 1, 3]]
        inp = df.copy()
        return_value = inp.dropna(axis=1, thresh=5, inplace=True)
        tm.assert_frame_equal(dropped, expected)
        tm.assert_frame_equal(inp, expected)
        assert return_value is None

        dropped = df.dropna(axis=0, thresh=4)
        expected = df.loc[range(2, 6)]
        inp = df.copy()
        return_value = inp.dropna(axis=0, thresh=4, inplace=True)
        tm.assert_frame_equal(dropped, expected)
        tm.assert_frame_equal(inp, expected)
        assert return_value is None

        dropped = df.dropna(axis=1, thresh=4)
        tm.assert_frame_equal(dropped, df)

        dropped = df.dropna(axis=1, thresh=3)
        tm.assert_frame_equal(dropped, df)

        # subset
        dropped = df.dropna(axis=0, subset=[0, 1, 3])
        inp = df.copy()
        return_value = inp.dropna(axis=0, subset=[0, 1, 3], inplace=True)
        tm.assert_frame_equal(dropped, df)
        tm.assert_frame_equal(inp, df)
        assert return_value is None

        # all
        dropped = df.dropna(axis=1, how="all")
        tm.assert_frame_equal(dropped, df)

        df[2] = np.nan
        dropped = df.dropna(axis=1, how="all")
        expected = df.loc[:, [0, 1, 3]]
        tm.assert_frame_equal(dropped, expected)

        # bad input
        msg = "No axis named 3 for object type DataFrame"
        with pytest.raises(ValueError, match=msg):
            df.dropna(axis=3)

    def test_drop_and_dropna_caching(self):
        # tst that cacher updates
        original = Series([1, 2, np.nan], name="A")
        expected = Series([1, 2], dtype=original.dtype, name="A")
        df = pd.DataFrame({"A": original.values.copy()})
        df2 = df.copy()
        df["A"].dropna()
        tm.assert_series_equal(df["A"], original)

        ser = df["A"]
        return_value = ser.dropna(inplace=True)
        tm.assert_series_equal(ser, expected)
        tm.assert_series_equal(df["A"], original)
        assert return_value is None

        df2["A"].drop([1])
        tm.assert_series_equal(df2["A"], original)

        ser = df2["A"]
        return_value = ser.drop([1], inplace=True)
        tm.assert_series_equal(ser, original.drop([1]))
        tm.assert_series_equal(df2["A"], original)
        assert return_value is None

    def test_dropna_corner(self, float_frame):
        # bad input
        msg = "invalid how option: foo"
        with pytest.raises(ValueError, match=msg):
            float_frame.dropna(how="foo")
        msg = "must specify how or thresh"
        with pytest.raises(TypeError, match=msg):
            float_frame.dropna(how=None)
        # non-existent column - 8303
        with pytest.raises(KeyError, match=r"^\['X'\]$"):
            float_frame.dropna(subset=["A", "X"])

    def test_dropna_multiple_axes(self):
        df = DataFrame(
            [
                [1, np.nan, 2, 3],
                [4, np.nan, 5, 6],
                [np.nan, np.nan, np.nan, np.nan],
                [7, np.nan, 8, 9],
            ]
        )

        # GH20987
        with pytest.raises(TypeError, match="supplying multiple axes"):
            df.dropna(how="all", axis=[0, 1])
        with pytest.raises(TypeError, match="supplying multiple axes"):
            df.dropna(how="all", axis=(0, 1))

        inp = df.copy()
        with pytest.raises(TypeError, match="supplying multiple axes"):
            inp.dropna(how="all", axis=(0, 1), inplace=True)

    def test_dropna_tz_aware_datetime(self):
        # GH13407
        df = DataFrame()
        dt1 = datetime.datetime(2015, 1, 1, tzinfo=dateutil.tz.tzutc())
        dt2 = datetime.datetime(2015, 2, 2, tzinfo=dateutil.tz.tzutc())
        df["Time"] = [dt1]
        result = df.dropna(axis=0)
        expected = DataFrame({"Time": [dt1]})
        tm.assert_frame_equal(result, expected)

        # Ex2
        df = DataFrame({"Time": [dt1, None, np.nan, dt2]})
        result = df.dropna(axis=0)
        expected = DataFrame([dt1, dt2], columns=["Time"], index=[0, 3])
        tm.assert_frame_equal(result, expected)

    def test_dropna_categorical_interval_index(self):
        # GH 25087
        ii = pd.IntervalIndex.from_breaks([0, 2.78, 3.14, 6.28])
        ci = pd.CategoricalIndex(ii)
        df = pd.DataFrame({"A": list("abc")}, index=ci)

        expected = df
        result = df.dropna()
        tm.assert_frame_equal(result, expected)

    def test_fillna_datetime(self, datetime_frame):
        tf = datetime_frame
        tf.loc[tf.index[:5], "A"] = np.nan
        tf.loc[tf.index[-5:], "A"] = np.nan

        zero_filled = datetime_frame.fillna(0)
        assert (zero_filled.loc[zero_filled.index[:5], "A"] == 0).all()

        padded = datetime_frame.fillna(method="pad")
        assert np.isnan(padded.loc[padded.index[:5], "A"]).all()
        assert (
            padded.loc[padded.index[-5:], "A"] == padded.loc[padded.index[-5], "A"]
        ).all()

        msg = "Must specify a fill 'value' or 'method'"
        with pytest.raises(ValueError, match=msg):
            datetime_frame.fillna()
        msg = "Cannot specify both 'value' and 'method'"
        with pytest.raises(ValueError, match=msg):
            datetime_frame.fillna(5, method="ffill")

    def test_fillna_mixed_type(self, float_string_frame):

        mf = float_string_frame
        mf.loc[mf.index[5:20], "foo"] = np.nan
        mf.loc[mf.index[-10:], "A"] = np.nan
        # TODO: make stronger assertion here, GH 25640
        mf.fillna(value=0)
        mf.fillna(method="pad")

    def test_fillna_mixed_float(self, mixed_float_frame):

        # mixed numeric (but no float16)
        mf = mixed_float_frame.reindex(columns=["A", "B", "D"])
        mf.loc[mf.index[-10:], "A"] = np.nan
        result = mf.fillna(value=0)
        _check_mixed_float(result, dtype=dict(C=None))

        result = mf.fillna(method="pad")
        _check_mixed_float(result, dtype=dict(C=None))

    def test_fillna_empty(self):
        # empty frame (GH #2778)
        df = DataFrame(columns=["x"])
        for m in ["pad", "backfill"]:
            df.x.fillna(method=m, inplace=True)
            df.x.fillna(method=m)

    def test_fillna_different_dtype(self):
        # with different dtype (GH#3386)
        df = DataFrame(
            [["a", "a", np.nan, "a"], ["b", "b", np.nan, "b"], ["c", "c", np.nan, "c"]]
        )

        result = df.fillna({2: "foo"})
        expected = DataFrame(
            [["a", "a", "foo", "a"], ["b", "b", "foo", "b"], ["c", "c", "foo", "c"]]
        )
        tm.assert_frame_equal(result, expected)

        return_value = df.fillna({2: "foo"}, inplace=True)
        tm.assert_frame_equal(df, expected)
        assert return_value is None

    def test_fillna_limit_and_value(self):
        # limit and value
        df = DataFrame(np.random.randn(10, 3))
        df.iloc[2:7, 0] = np.nan
        df.iloc[3:5, 2] = np.nan

        expected = df.copy()
        expected.iloc[2, 0] = 999
        expected.iloc[3, 2] = 999
        result = df.fillna(999, limit=1)
        tm.assert_frame_equal(result, expected)

    def test_fillna_datelike(self):
        # with datelike
        # GH#6344
        df = DataFrame(
            {
                "Date": [pd.NaT, Timestamp("2014-1-1")],
                "Date2": [Timestamp("2013-1-1"), pd.NaT],
            }
        )

        expected = df.copy()
        expected["Date"] = expected["Date"].fillna(df.loc[df.index[0], "Date2"])
        result = df.fillna(value={"Date": df["Date2"]})
        tm.assert_frame_equal(result, expected)

    def test_fillna_tzaware(self):
        # with timezone
        # GH#15855
        df = pd.DataFrame({"A": [pd.Timestamp("2012-11-11 00:00:00+01:00"), pd.NaT]})
        exp = pd.DataFrame(
            {
                "A": [
                    pd.Timestamp("2012-11-11 00:00:00+01:00"),
                    pd.Timestamp("2012-11-11 00:00:00+01:00"),
                ]
            }
        )
        tm.assert_frame_equal(df.fillna(method="pad"), exp)

        df = pd.DataFrame({"A": [pd.NaT, pd.Timestamp("2012-11-11 00:00:00+01:00")]})
        exp = pd.DataFrame(
            {
                "A": [
                    pd.Timestamp("2012-11-11 00:00:00+01:00"),
                    pd.Timestamp("2012-11-11 00:00:00+01:00"),
                ]
            }
        )
        tm.assert_frame_equal(df.fillna(method="bfill"), exp)

    def test_fillna_tzaware_different_column(self):
        # with timezone in another column
        # GH#15522
        df = pd.DataFrame(
            {
                "A": pd.date_range("20130101", periods=4, tz="US/Eastern"),
                "B": [1, 2, np.nan, np.nan],
            }
        )
        result = df.fillna(method="pad")
        expected = pd.DataFrame(
            {
                "A": pd.date_range("20130101", periods=4, tz="US/Eastern"),
                "B": [1.0, 2.0, 2.0, 2.0],
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_na_actions_categorical(self):

        cat = Categorical([1, 2, 3, np.nan], categories=[1, 2, 3])
        vals = ["a", "b", np.nan, "d"]
        df = DataFrame({"cats": cat, "vals": vals})
        cat2 = Categorical([1, 2, 3, 3], categories=[1, 2, 3])
        vals2 = ["a", "b", "b", "d"]
        df_exp_fill = DataFrame({"cats": cat2, "vals": vals2})
        cat3 = Categorical([1, 2, 3], categories=[1, 2, 3])
        vals3 = ["a", "b", np.nan]
        df_exp_drop_cats = DataFrame({"cats": cat3, "vals": vals3})
        cat4 = Categorical([1, 2], categories=[1, 2, 3])
        vals4 = ["a", "b"]
        df_exp_drop_all = DataFrame({"cats": cat4, "vals": vals4})

        # fillna
        res = df.fillna(value={"cats": 3, "vals": "b"})
        tm.assert_frame_equal(res, df_exp_fill)

        with pytest.raises(ValueError, match=("fill value must be in categories")):
            df.fillna(value={"cats": 4, "vals": "c"})

        res = df.fillna(method="pad")
        tm.assert_frame_equal(res, df_exp_fill)

        # dropna
        res = df.dropna(subset=["cats"])
        tm.assert_frame_equal(res, df_exp_drop_cats)

        res = df.dropna()
        tm.assert_frame_equal(res, df_exp_drop_all)

        # make sure that fillna takes missing values into account
        c = Categorical([np.nan, "b", np.nan], categories=["a", "b"])
        df = pd.DataFrame({"cats": c, "vals": [1, 2, 3]})

        cat_exp = Categorical(["a", "b", "a"], categories=["a", "b"])
        df_exp = DataFrame({"cats": cat_exp, "vals": [1, 2, 3]})

        res = df.fillna("a")
        tm.assert_frame_equal(res, df_exp)

    def test_fillna_categorical_nan(self):
        # GH 14021
        # np.nan should always be a valid filler
        cat = Categorical([np.nan, 2, np.nan])
        val = Categorical([np.nan, np.nan, np.nan])
        df = DataFrame({"cats": cat, "vals": val})

        # GH#32950 df.median() is poorly behaved because there is no
        #  Categorical.median
        median = Series({"cats": 2.0, "vals": np.nan})

        res = df.fillna(median)
        v_exp = [np.nan, np.nan, np.nan]
        df_exp = DataFrame({"cats": [2, 2, 2], "vals": v_exp}, dtype="category")
        tm.assert_frame_equal(res, df_exp)

        result = df.cats.fillna(np.nan)
        tm.assert_series_equal(result, df.cats)

        result = df.vals.fillna(np.nan)
        tm.assert_series_equal(result, df.vals)

        idx = pd.DatetimeIndex(
            ["2011-01-01 09:00", "2016-01-01 23:45", "2011-01-01 09:00", pd.NaT, pd.NaT]
        )
        df = DataFrame({"a": Categorical(idx)})
        tm.assert_frame_equal(df.fillna(value=pd.NaT), df)

        idx = pd.PeriodIndex(
            ["2011-01", "2011-01", "2011-01", pd.NaT, pd.NaT], freq="M"
        )
        df = DataFrame({"a": Categorical(idx)})
        tm.assert_frame_equal(df.fillna(value=pd.NaT), df)

        idx = pd.TimedeltaIndex(["1 days", "2 days", "1 days", pd.NaT, pd.NaT])
        df = DataFrame({"a": Categorical(idx)})
        tm.assert_frame_equal(df.fillna(value=pd.NaT), df)

    def test_fillna_downcast(self):
        # GH 15277
        # infer int64 from float64
        df = pd.DataFrame({"a": [1.0, np.nan]})
        result = df.fillna(0, downcast="infer")
        expected = pd.DataFrame({"a": [1, 0]})
        tm.assert_frame_equal(result, expected)

        # infer int64 from float64 when fillna value is a dict
        df = pd.DataFrame({"a": [1.0, np.nan]})
        result = df.fillna({"a": 0}, downcast="infer")
        expected = pd.DataFrame({"a": [1, 0]})
        tm.assert_frame_equal(result, expected)

    def test_fillna_dtype_conversion(self):
        # make sure that fillna on an empty frame works
        df = DataFrame(index=["A", "B", "C"], columns=[1, 2, 3, 4, 5])
        result = df.dtypes
        expected = Series([np.dtype("object")] * 5, index=[1, 2, 3, 4, 5])
        tm.assert_series_equal(result, expected)

        result = df.fillna(1)
        expected = DataFrame(1, index=["A", "B", "C"], columns=[1, 2, 3, 4, 5])
        tm.assert_frame_equal(result, expected)

        # empty block
        df = DataFrame(index=range(3), columns=["A", "B"], dtype="float64")
        result = df.fillna("nan")
        expected = DataFrame("nan", index=range(3), columns=["A", "B"])
        tm.assert_frame_equal(result, expected)

        # equiv of replace
        df = DataFrame(dict(A=[1, np.nan], B=[1.0, 2.0]))
        for v in ["", 1, np.nan, 1.0]:
            expected = df.replace(np.nan, v)
            result = df.fillna(v)
            tm.assert_frame_equal(result, expected)

    def test_fillna_datetime_columns(self):
        # GH 7095
        df = pd.DataFrame(
            {
                "A": [-1, -2, np.nan],
                "B": date_range("20130101", periods=3),
                "C": ["foo", "bar", None],
                "D": ["foo2", "bar2", None],
            },
            index=date_range("20130110", periods=3),
        )
        result = df.fillna("?")
        expected = pd.DataFrame(
            {
                "A": [-1, -2, "?"],
                "B": date_range("20130101", periods=3),
                "C": ["foo", "bar", "?"],
                "D": ["foo2", "bar2", "?"],
            },
            index=date_range("20130110", periods=3),
        )
        tm.assert_frame_equal(result, expected)

        df = pd.DataFrame(
            {
                "A": [-1, -2, np.nan],
                "B": [pd.Timestamp("2013-01-01"), pd.Timestamp("2013-01-02"), pd.NaT],
                "C": ["foo", "bar", None],
                "D": ["foo2", "bar2", None],
            },
            index=date_range("20130110", periods=3),
        )
        result = df.fillna("?")
        expected = pd.DataFrame(
            {
                "A": [-1, -2, "?"],
                "B": [pd.Timestamp("2013-01-01"), pd.Timestamp("2013-01-02"), "?"],
                "C": ["foo", "bar", "?"],
                "D": ["foo2", "bar2", "?"],
            },
            index=pd.date_range("20130110", periods=3),
        )
        tm.assert_frame_equal(result, expected)

    def test_ffill(self, datetime_frame):
        datetime_frame["A"][:5] = np.nan
        datetime_frame["A"][-5:] = np.nan

        tm.assert_frame_equal(
            datetime_frame.ffill(), datetime_frame.fillna(method="ffill")
        )

    def test_bfill(self, datetime_frame):
        datetime_frame["A"][:5] = np.nan
        datetime_frame["A"][-5:] = np.nan

        tm.assert_frame_equal(
            datetime_frame.bfill(), datetime_frame.fillna(method="bfill")
        )

    def test_frame_pad_backfill_limit(self):
        index = np.arange(10)
        df = DataFrame(np.random.randn(10, 4), index=index)

        result = df[:2].reindex(index, method="pad", limit=5)

        expected = df[:2].reindex(index).fillna(method="pad")
        expected.values[-3:] = np.nan
        tm.assert_frame_equal(result, expected)

        result = df[-2:].reindex(index, method="backfill", limit=5)

        expected = df[-2:].reindex(index).fillna(method="backfill")
        expected.values[:3] = np.nan
        tm.assert_frame_equal(result, expected)

    def test_frame_fillna_limit(self):
        index = np.arange(10)
        df = DataFrame(np.random.randn(10, 4), index=index)

        result = df[:2].reindex(index)
        result = result.fillna(method="pad", limit=5)

        expected = df[:2].reindex(index).fillna(method="pad")
        expected.values[-3:] = np.nan
        tm.assert_frame_equal(result, expected)

        result = df[-2:].reindex(index)
        result = result.fillna(method="backfill", limit=5)

        expected = df[-2:].reindex(index).fillna(method="backfill")
        expected.values[:3] = np.nan
        tm.assert_frame_equal(result, expected)

    def test_fillna_skip_certain_blocks(self):
        # don't try to fill boolean, int blocks

        df = DataFrame(np.random.randn(10, 4).astype(int))

        # it works!
        df.fillna(np.nan)

    @pytest.mark.parametrize("type", [int, float])
    def test_fillna_positive_limit(self, type):
        df = DataFrame(np.random.randn(10, 4)).astype(type)

        msg = "Limit must be greater than 0"
        with pytest.raises(ValueError, match=msg):
            df.fillna(0, limit=-5)

    @pytest.mark.parametrize("type", [int, float])
    def test_fillna_integer_limit(self, type):
        df = DataFrame(np.random.randn(10, 4)).astype(type)

        msg = "Limit must be an integer"
        with pytest.raises(ValueError, match=msg):
            df.fillna(0, limit=0.5)

    def test_fillna_inplace(self):
        df = DataFrame(np.random.randn(10, 4))
        df[1][:4] = np.nan
        df[3][-4:] = np.nan

        expected = df.fillna(value=0)
        assert expected is not df

        df.fillna(value=0, inplace=True)
        tm.assert_frame_equal(df, expected)

        expected = df.fillna(value={0: 0}, inplace=True)
        assert expected is None

        df[1][:4] = np.nan
        df[3][-4:] = np.nan
        expected = df.fillna(method="ffill")
        assert expected is not df

        df.fillna(method="ffill", inplace=True)
        tm.assert_frame_equal(df, expected)

    def test_fillna_dict_series(self):
        df = DataFrame(
            {
                "a": [np.nan, 1, 2, np.nan, np.nan],
                "b": [1, 2, 3, np.nan, np.nan],
                "c": [np.nan, 1, 2, 3, 4],
            }
        )

        result = df.fillna({"a": 0, "b": 5})

        expected = df.copy()
        expected["a"] = expected["a"].fillna(0)
        expected["b"] = expected["b"].fillna(5)
        tm.assert_frame_equal(result, expected)

        # it works
        result = df.fillna({"a": 0, "b": 5, "d": 7})

        # Series treated same as dict
        result = df.fillna(df.max())
        expected = df.fillna(df.max().to_dict())
        tm.assert_frame_equal(result, expected)

        # disable this for now
        with pytest.raises(NotImplementedError, match="column by column"):
            df.fillna(df.max(1), axis=1)

    def test_fillna_dataframe(self):
        # GH 8377
        df = DataFrame(
            {
                "a": [np.nan, 1, 2, np.nan, np.nan],
                "b": [1, 2, 3, np.nan, np.nan],
                "c": [np.nan, 1, 2, 3, 4],
            },
            index=list("VWXYZ"),
        )

        # df2 may have different index and columns
        df2 = DataFrame(
            {
                "a": [np.nan, 10, 20, 30, 40],
                "b": [50, 60, 70, 80, 90],
                "foo": ["bar"] * 5,
            },
            index=list("VWXuZ"),
        )

        result = df.fillna(df2)

        # only those columns and indices which are shared get filled
        expected = DataFrame(
            {
                "a": [np.nan, 1, 2, np.nan, 40],
                "b": [1, 2, 3, np.nan, 90],
                "c": [np.nan, 1, 2, 3, 4],
            },
            index=list("VWXYZ"),
        )

        tm.assert_frame_equal(result, expected)

    def test_fillna_columns(self):
        df = DataFrame(np.random.randn(10, 10))
        df.values[:, ::2] = np.nan

        result = df.fillna(method="ffill", axis=1)
        expected = df.T.fillna(method="pad").T
        tm.assert_frame_equal(result, expected)

        df.insert(6, "foo", 5)
        result = df.fillna(method="ffill", axis=1)
        expected = df.astype(float).fillna(method="ffill", axis=1)
        tm.assert_frame_equal(result, expected)

    def test_fillna_invalid_method(self, float_frame):
        with pytest.raises(ValueError, match="ffil"):
            float_frame.fillna(method="ffil")

    def test_fillna_invalid_value(self, float_frame):
        # list
        msg = '"value" parameter must be a scalar or dict, but you passed a "{}"'
        with pytest.raises(TypeError, match=msg.format("list")):
            float_frame.fillna([1, 2])
        # tuple
        with pytest.raises(TypeError, match=msg.format("tuple")):
            float_frame.fillna((1, 2))
        # frame with series
        msg = (
            '"value" parameter must be a scalar, dict or Series, but you '
            'passed a "DataFrame"'
        )
        with pytest.raises(TypeError, match=msg):
            float_frame.iloc[:, 0].fillna(float_frame)

    def test_fillna_col_reordering(self):
        cols = ["COL." + str(i) for i in range(5, 0, -1)]
        data = np.random.rand(20, 5)
        df = DataFrame(index=range(20), columns=cols, data=data)
        filled = df.fillna(method="ffill")
        assert df.columns.tolist() == filled.columns.tolist()

    def test_fill_corner(self, float_frame, float_string_frame):
        mf = float_string_frame
        mf.loc[mf.index[5:20], "foo"] = np.nan
        mf.loc[mf.index[-10:], "A"] = np.nan

        filled = float_string_frame.fillna(value=0)
        assert (filled.loc[filled.index[5:20], "foo"] == 0).all()
        del float_string_frame["foo"]

        empty_float = float_frame.reindex(columns=[])

        # TODO(wesm): unused?
        result = empty_float.fillna(value=0)  # noqa
