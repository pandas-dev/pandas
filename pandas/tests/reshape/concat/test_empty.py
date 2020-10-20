import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    Series,
    concat,
    date_range,
)
import pandas._testing as tm


@pytest.fixture(params=[True, False])
def sort(request):
    """Boolean sort keyword for concat and DataFrame.append."""
    return request.param


class TestEmptyConcat:
    def test_handle_empty_objects(self, sort):
        df = DataFrame(np.random.randn(10, 4), columns=list("abcd"))

        baz = df[:5].copy()
        baz["foo"] = "bar"
        empty = df[5:5]

        frames = [baz, empty, empty, df[5:]]
        concatted = concat(frames, axis=0, sort=sort)

        expected = df.reindex(columns=["a", "b", "c", "d", "foo"])
        expected["foo"] = expected["foo"].astype("O")
        expected.loc[0:4, "foo"] = "bar"

        tm.assert_frame_equal(concatted, expected)

        # empty as first element with time series
        # GH3259
        df = DataFrame(
            dict(A=range(10000)), index=date_range("20130101", periods=10000, freq="s")
        )
        empty = DataFrame()
        result = concat([df, empty], axis=1)
        tm.assert_frame_equal(result, df)
        result = concat([empty, df], axis=1)
        tm.assert_frame_equal(result, df)

        result = concat([df, empty])
        tm.assert_frame_equal(result, df)
        result = concat([empty, df])
        tm.assert_frame_equal(result, df)

    def test_concat_empty_series(self):
        # GH 11082
        s1 = pd.Series([1, 2, 3], name="x")
        s2 = pd.Series(name="y", dtype="float64")
        res = pd.concat([s1, s2], axis=1)
        exp = pd.DataFrame(
            {"x": [1, 2, 3], "y": [np.nan, np.nan, np.nan]},
            index=pd.Index([0, 1, 2], dtype="O"),
        )
        tm.assert_frame_equal(res, exp)

        s1 = pd.Series([1, 2, 3], name="x")
        s2 = pd.Series(name="y", dtype="float64")
        res = pd.concat([s1, s2], axis=0)
        # name will be reset
        exp = pd.Series([1, 2, 3])
        tm.assert_series_equal(res, exp)

        # empty Series with no name
        s1 = pd.Series([1, 2, 3], name="x")
        s2 = pd.Series(name=None, dtype="float64")
        res = pd.concat([s1, s2], axis=1)
        exp = pd.DataFrame(
            {"x": [1, 2, 3], 0: [np.nan, np.nan, np.nan]},
            columns=["x", 0],
            index=pd.Index([0, 1, 2], dtype="O"),
        )
        tm.assert_frame_equal(res, exp)

    def test_concat_empty_df_object_dtype(self):
        # GH 9149
        df_1 = pd.DataFrame({"Row": [0, 1, 1], "EmptyCol": np.nan, "NumberCol": [1, 2, 3]})
        df_2 = pd.DataFrame(columns=df_1.columns)
        result = pd.concat([df_1, df_2], axis=0)
        expected = df_1.astype(object)
        tm.assert_frame_equal(result, expected)

    def test_concat_inner_join_empty(self):
        # GH 15328
        df_empty = pd.DataFrame()
        df_a = pd.DataFrame({"a": [1, 2]}, index=[0, 1], dtype="int64")
        df_expected = pd.DataFrame({"a": []}, index=[], dtype="int64")

        for how, expected in [("inner", df_expected), ("outer", df_a)]:
            result = pd.concat([df_a, df_empty], axis=1, join=how)
            tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("tz", [None, "UTC"])
    @pytest.mark.parametrize("values", [[], [1, 2, 3]])
    def test_concat_empty_series_timelike(self, tz, values):
        # GH 18447

        first = Series([], dtype="M8[ns]").dt.tz_localize(tz)
        dtype = None if values else np.float64
        second = Series(values, dtype=dtype)

        expected = DataFrame(
            {
                0: pd.Series([pd.NaT] * len(values), dtype="M8[ns]").dt.tz_localize(tz),
                1: values,
            }
        )
        result = concat([first, second], axis=1)
        tm.assert_frame_equal(result, expected)

    def test_empty_dtype_coerce(self):

        # xref to #12411
        # xref to #12045
        # xref to #11594
        # see below

        # 10571
        df1 = DataFrame(data=[[1, None], [2, None]], columns=["a", "b"])
        df2 = DataFrame(data=[[3, None], [4, None]], columns=["a", "b"])
        result = concat([df1, df2])
        expected = df1.dtypes
        tm.assert_series_equal(result.dtypes, expected)
