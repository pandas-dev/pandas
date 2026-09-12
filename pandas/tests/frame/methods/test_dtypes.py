from datetime import timedelta

import numpy as np
import pytest

from pandas.core.dtypes.dtypes import DatetimeTZDtype

import pandas as pd
import pandas._testing as tm


class TestDataFrameDataTypes:
    def test_empty_frame_dtypes(self):
        empty_df = pd.DataFrame()
        tm.assert_series_equal(empty_df.dtypes, pd.Series(dtype=object))

        nocols_df = pd.DataFrame(index=[1, 2, 3])
        tm.assert_series_equal(nocols_df.dtypes, pd.Series(dtype=object))

        norows_df = pd.DataFrame(columns=list("abc"))
        tm.assert_series_equal(norows_df.dtypes, pd.Series(object, index=list("abc")))

        norows_int_df = pd.DataFrame(columns=list("abc")).astype(np.int32)
        tm.assert_series_equal(
            norows_int_df.dtypes, pd.Series(np.dtype("int32"), index=list("abc"))
        )

        df = pd.DataFrame({"a": 1, "b": True, "c": 1.0}, index=[1, 2, 3])
        ex_dtypes = pd.Series({"a": np.int64, "b": np.bool_, "c": np.float64})
        tm.assert_series_equal(df.dtypes, ex_dtypes)

        # same but for empty slice of df
        tm.assert_series_equal(df[:0].dtypes, ex_dtypes)

    def test_datetime_with_tz_dtypes(self):
        tzframe = pd.DataFrame(
            {
                "A": pd.date_range("20130101", periods=3, unit="ns"),
                "B": pd.date_range("20130101", periods=3, tz="US/Eastern", unit="ns"),
                "C": pd.date_range("20130101", periods=3, tz="CET", unit="ns"),
            }
        )
        tzframe.iloc[1, 1] = pd.NaT
        tzframe.iloc[1, 2] = pd.NaT
        result = tzframe.dtypes.sort_index()
        expected = pd.Series(
            [
                np.dtype("datetime64[ns]"),
                DatetimeTZDtype("ns", "US/Eastern"),
                DatetimeTZDtype("ns", "CET"),
            ],
            ["A", "B", "C"],
        )

        tm.assert_series_equal(result, expected)

    def test_dtypes_are_correct_after_column_slice(self):
        # GH6525
        df = pd.DataFrame(index=range(5), columns=list("abc"), dtype=np.float64)
        tm.assert_series_equal(
            df.dtypes,
            pd.Series({"a": np.float64, "b": np.float64, "c": np.float64}),
        )
        tm.assert_series_equal(df.iloc[:, 2:].dtypes, pd.Series({"c": np.float64}))
        tm.assert_series_equal(
            df.dtypes,
            pd.Series({"a": np.float64, "b": np.float64, "c": np.float64}),
        )

    @pytest.mark.parametrize(
        "data",
        [pd.NA, True],
    )
    def test_dtypes_are_correct_after_groupby_last(self, data):
        # GH46409
        df = pd.DataFrame(
            {"id": [1, 2, 3, 4], "test": [True, pd.NA, data, False]}
        ).convert_dtypes()
        result = df.groupby("id").last().test
        expected = df.set_index("id").test
        assert result.dtype == pd.BooleanDtype()
        tm.assert_series_equal(expected, result)

    def test_dtypes_gh8722(self, float_string_frame):
        float_string_frame["bool"] = float_string_frame["A"] > 0
        result = float_string_frame.dtypes
        expected = pd.Series(
            {k: v.dtype for k, v in float_string_frame.items()}, index=result.index
        )
        tm.assert_series_equal(result, expected)

    def test_dtypes_timedeltas(self):
        df = pd.DataFrame(
            {
                "A": pd.Series(
                    pd.date_range("2012-1-1", periods=3, freq="D", unit="ns")
                ),
                "B": pd.Series([timedelta(days=i) for i in range(3)]),
            }
        )
        result = df.dtypes
        expected = pd.Series(
            [np.dtype("datetime64[ns]"), np.dtype("timedelta64[us]")], index=list("AB")
        )
        tm.assert_series_equal(result, expected)

        df["C"] = df["A"] + df["B"]
        result = df.dtypes
        expected = pd.Series(
            [
                np.dtype("datetime64[ns]"),
                np.dtype("timedelta64[us]"),
                np.dtype("datetime64[ns]"),
            ],
            index=list("ABC"),
        )
        tm.assert_series_equal(result, expected)

        # mixed int types
        df["D"] = 1
        result = df.dtypes
        expected = pd.Series(
            [
                np.dtype("datetime64[ns]"),
                np.dtype("timedelta64[us]"),
                np.dtype("datetime64[ns]"),
                np.dtype("int64"),
            ],
            index=list("ABCD"),
        )
        tm.assert_series_equal(result, expected)

    def test_frame_apply_np_array_return_type(self, using_infer_string):
        # GH 35517
        df = pd.DataFrame([["foo"]])
        result = df.apply(lambda col: np.array("bar"))
        expected = pd.Series(np.array("bar"))
        tm.assert_series_equal(result, expected)
