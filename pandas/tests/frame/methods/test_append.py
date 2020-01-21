import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Series, Timestamp
import pandas._testing as tm


class TestDataFrameAppend:
    def test_append_empty_list(self):
        # GH 28769
        df = DataFrame()
        result = df.append([])
        expected = df
        tm.assert_frame_equal(result, expected)
        assert result is not df

        df = DataFrame(np.random.randn(5, 4), columns=["foo", "bar", "baz", "qux"])
        result = df.append([])
        expected = df
        tm.assert_frame_equal(result, expected)
        assert result is not df  # .append() should return a new object

    def test_append_series_dict(self):
        df = DataFrame(np.random.randn(5, 4), columns=["foo", "bar", "baz", "qux"])

        series = df.loc[4]
        msg = "Indexes have overlapping values"
        with pytest.raises(ValueError, match=msg):
            df.append(series, verify_integrity=True)

        series.name = None
        msg = "Can only append a Series if ignore_index=True"
        with pytest.raises(TypeError, match=msg):
            df.append(series, verify_integrity=True)

        result = df.append(series[::-1], ignore_index=True)
        expected = df.append(
            DataFrame({0: series[::-1]}, index=df.columns).T, ignore_index=True
        )
        tm.assert_frame_equal(result, expected)

        # dict
        result = df.append(series.to_dict(), ignore_index=True)
        tm.assert_frame_equal(result, expected)

        result = df.append(series[::-1][:3], ignore_index=True)
        expected = df.append(
            DataFrame({0: series[::-1][:3]}).T, ignore_index=True, sort=True
        )
        tm.assert_frame_equal(result, expected.loc[:, result.columns])

        msg = "Can only append a dict if ignore_index=True"
        with pytest.raises(TypeError, match=msg):
            df.append(series.to_dict())

        # can append when name set
        row = df.loc[4]
        row.name = 5
        result = df.append(row)
        expected = df.append(df[-1:], ignore_index=True)
        tm.assert_frame_equal(result, expected)

    def test_append_list_of_series_dicts(self):
        df = DataFrame(np.random.randn(5, 4), columns=["foo", "bar", "baz", "qux"])

        dicts = [x.to_dict() for idx, x in df.iterrows()]

        result = df.append(dicts, ignore_index=True)
        expected = df.append(df, ignore_index=True)
        tm.assert_frame_equal(result, expected)

        # different columns
        dicts = [
            {"foo": 1, "bar": 2, "baz": 3, "peekaboo": 4},
            {"foo": 5, "bar": 6, "baz": 7, "peekaboo": 8},
        ]
        result = df.append(dicts, ignore_index=True, sort=True)
        expected = df.append(DataFrame(dicts), ignore_index=True, sort=True)
        tm.assert_frame_equal(result, expected)

    def test_append_missing_cols(self):
        # GH22252
        # exercise the conditional branch in append method where the data
        # to be appended is a list and does not contain all columns that are in
        # the target DataFrame
        df = DataFrame(np.random.randn(5, 4), columns=["foo", "bar", "baz", "qux"])

        dicts = [{"foo": 9}, {"bar": 10}]
        with tm.assert_produces_warning(None):
            result = df.append(dicts, ignore_index=True, sort=True)

        expected = df.append(DataFrame(dicts), ignore_index=True, sort=True)
        tm.assert_frame_equal(result, expected)

    def test_append_empty_dataframe(self):

        # Empty df append empty df
        df1 = DataFrame()
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        tm.assert_frame_equal(result, expected)

        # Non-empty df append empty df
        df1 = DataFrame(np.random.randn(5, 2))
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        tm.assert_frame_equal(result, expected)

        # Empty df with columns append empty df
        df1 = DataFrame(columns=["bar", "foo"])
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        tm.assert_frame_equal(result, expected)

        # Non-Empty df with columns append empty df
        df1 = DataFrame(np.random.randn(5, 2), columns=["bar", "foo"])
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        tm.assert_frame_equal(result, expected)

    def test_append_dtypes(self):

        # GH 5754
        # row appends of different dtypes (so need to do by-item)
        # can sometimes infer the correct type

        df1 = DataFrame({"bar": Timestamp("20130101")}, index=range(5))
        df2 = DataFrame()
        result = df1.append(df2)
        expected = df1.copy()
        tm.assert_frame_equal(result, expected)

        df1 = DataFrame({"bar": Timestamp("20130101")}, index=range(1))
        df2 = DataFrame({"bar": "foo"}, index=range(1, 2))
        result = df1.append(df2)
        expected = DataFrame({"bar": [Timestamp("20130101"), "foo"]})
        tm.assert_frame_equal(result, expected)

        df1 = DataFrame({"bar": Timestamp("20130101")}, index=range(1))
        df2 = DataFrame({"bar": np.nan}, index=range(1, 2))
        result = df1.append(df2)
        expected = DataFrame(
            {"bar": Series([Timestamp("20130101"), np.nan], dtype="M8[ns]")}
        )
        tm.assert_frame_equal(result, expected)

        df1 = DataFrame({"bar": Timestamp("20130101")}, index=range(1))
        df2 = DataFrame({"bar": np.nan}, index=range(1, 2), dtype=object)
        result = df1.append(df2)
        expected = DataFrame(
            {"bar": Series([Timestamp("20130101"), np.nan], dtype="M8[ns]")}
        )
        tm.assert_frame_equal(result, expected)

        df1 = DataFrame({"bar": np.nan}, index=range(1))
        df2 = DataFrame({"bar": Timestamp("20130101")}, index=range(1, 2))
        result = df1.append(df2)
        expected = DataFrame(
            {"bar": Series([np.nan, Timestamp("20130101")], dtype="M8[ns]")}
        )
        tm.assert_frame_equal(result, expected)

        df1 = DataFrame({"bar": Timestamp("20130101")}, index=range(1))
        df2 = DataFrame({"bar": 1}, index=range(1, 2), dtype=object)
        result = df1.append(df2)
        expected = DataFrame({"bar": Series([Timestamp("20130101"), 1])})
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "timestamp", ["2019-07-19 07:04:57+0100", "2019-07-19 07:04:57"]
    )
    def test_append_timestamps_aware_or_naive(self, tz_naive_fixture, timestamp):
        # GH 30238
        tz = tz_naive_fixture
        df = pd.DataFrame([pd.Timestamp(timestamp, tz=tz)])
        result = df.append(df.iloc[0]).iloc[-1]
        expected = pd.Series(pd.Timestamp(timestamp, tz=tz), name=0)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "data, dtype",
        [
            ([1], pd.Int64Dtype()),
            ([1], pd.CategoricalDtype()),
            ([pd.Interval(left=0, right=5)], pd.IntervalDtype()),
            ([pd.Period("2000-03", freq="M")], pd.PeriodDtype("M")),
            ([1], pd.SparseDtype()),
        ],
    )
    def test_other_dtypes(self, data, dtype):
        df = pd.DataFrame(data, dtype=dtype)
        result = df.append(df.iloc[0]).iloc[-1]
        expected = pd.Series(data, name=0, dtype=dtype)
        tm.assert_series_equal(result, expected)
