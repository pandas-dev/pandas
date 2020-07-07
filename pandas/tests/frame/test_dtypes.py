from collections import OrderedDict
from datetime import timedelta

import numpy as np
import pytest

from pandas.core.dtypes.dtypes import DatetimeTZDtype

import pandas as pd
from pandas import DataFrame, Series, Timestamp, date_range, option_context
import pandas._testing as tm


def _check_cast(df, v):
    """
    Check if all dtypes of df are equal to v
    """
    assert all(s.dtype.name == v for _, s in df.items())


class TestDataFrameDataTypes:
    def test_concat_empty_dataframe_dtypes(self):
        df = DataFrame(columns=list("abc"))
        df["a"] = df["a"].astype(np.bool_)
        df["b"] = df["b"].astype(np.int32)
        df["c"] = df["c"].astype(np.float64)

        result = pd.concat([df, df])
        assert result["a"].dtype == np.bool_
        assert result["b"].dtype == np.int32
        assert result["c"].dtype == np.float64

        result = pd.concat([df, df.astype(np.float64)])
        assert result["a"].dtype == np.object_
        assert result["b"].dtype == np.float64
        assert result["c"].dtype == np.float64

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

        odict = OrderedDict
        df = pd.DataFrame(odict([("a", 1), ("b", True), ("c", 1.0)]), index=[1, 2, 3])
        ex_dtypes = pd.Series(
            odict([("a", np.int64), ("b", np.bool_), ("c", np.float64)])
        )
        tm.assert_series_equal(df.dtypes, ex_dtypes)

        # same but for empty slice of df
        tm.assert_series_equal(df[:0].dtypes, ex_dtypes)

    def test_datetime_with_tz_dtypes(self):
        tzframe = DataFrame(
            {
                "A": date_range("20130101", periods=3),
                "B": date_range("20130101", periods=3, tz="US/Eastern"),
                "C": date_range("20130101", periods=3, tz="CET"),
            }
        )
        tzframe.iloc[1, 1] = pd.NaT
        tzframe.iloc[1, 2] = pd.NaT
        result = tzframe.dtypes.sort_index()
        expected = Series(
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
        df = pd.DataFrame(index=range(5), columns=list("abc"), dtype=np.float_)
        odict = OrderedDict
        tm.assert_series_equal(
            df.dtypes,
            pd.Series(odict([("a", np.float_), ("b", np.float_), ("c", np.float_)])),
        )
        tm.assert_series_equal(
            df.iloc[:, 2:].dtypes, pd.Series(odict([("c", np.float_)]))
        )
        tm.assert_series_equal(
            df.dtypes,
            pd.Series(odict([("a", np.float_), ("b", np.float_), ("c", np.float_)])),
        )

    def test_dtypes_gh8722(self, float_string_frame):
        float_string_frame["bool"] = float_string_frame["A"] > 0
        result = float_string_frame.dtypes
        expected = Series(
            {k: v.dtype for k, v in float_string_frame.items()}, index=result.index
        )
        tm.assert_series_equal(result, expected)

        # compat, GH 8722
        with option_context("use_inf_as_na", True):
            df = DataFrame([[1]])
            result = df.dtypes
            tm.assert_series_equal(result, Series({0: np.dtype("int64")}))

    def test_singlerow_slice_categoricaldtype_gives_series(self):
        # GH29521
        df = pd.DataFrame({"x": pd.Categorical("a b c d e".split())})
        result = df.iloc[0]
        raw_cat = pd.Categorical(["a"], categories=["a", "b", "c", "d", "e"])
        expected = pd.Series(raw_cat, index=["x"], name=0, dtype="category")

        tm.assert_series_equal(result, expected)

    def test_timedeltas(self):
        df = DataFrame(
            dict(
                A=Series(date_range("2012-1-1", periods=3, freq="D")),
                B=Series([timedelta(days=i) for i in range(3)]),
            )
        )
        result = df.dtypes
        expected = Series(
            [np.dtype("datetime64[ns]"), np.dtype("timedelta64[ns]")], index=list("AB")
        )
        tm.assert_series_equal(result, expected)

        df["C"] = df["A"] + df["B"]
        result = df.dtypes
        expected = Series(
            [
                np.dtype("datetime64[ns]"),
                np.dtype("timedelta64[ns]"),
                np.dtype("datetime64[ns]"),
            ],
            index=list("ABC"),
        )
        tm.assert_series_equal(result, expected)

        # mixed int types
        df["D"] = 1
        result = df.dtypes
        expected = Series(
            [
                np.dtype("datetime64[ns]"),
                np.dtype("timedelta64[ns]"),
                np.dtype("datetime64[ns]"),
                np.dtype("int64"),
            ],
            index=list("ABCD"),
        )
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "input_vals",
        [
            ([1, 2]),
            (["1", "2"]),
            (list(pd.date_range("1/1/2011", periods=2, freq="H"))),
            (list(pd.date_range("1/1/2011", periods=2, freq="H", tz="US/Eastern"))),
            ([pd.Interval(left=0, right=5)]),
        ],
    )
    def test_constructor_list_str(self, input_vals, string_dtype):
        # GH 16605
        # Ensure that data elements are converted to strings when
        # dtype is str, 'str', or 'U'

        result = DataFrame({"A": input_vals}, dtype=string_dtype)
        expected = DataFrame({"A": input_vals}).astype({"A": string_dtype})
        tm.assert_frame_equal(result, expected)

    def test_constructor_list_str_na(self, string_dtype):

        result = DataFrame({"A": [1.0, 2.0, None]}, dtype=string_dtype)
        expected = DataFrame({"A": ["1.0", "2.0", None]}, dtype=object)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "data, expected",
        [
            # empty
            (DataFrame(), True),
            # multi-same
            (DataFrame({"A": [1, 2], "B": [1, 2]}), True),
            # multi-object
            (
                DataFrame(
                    {
                        "A": np.array([1, 2], dtype=object),
                        "B": np.array(["a", "b"], dtype=object),
                    }
                ),
                True,
            ),
            # multi-extension
            (
                DataFrame(
                    {"A": pd.Categorical(["a", "b"]), "B": pd.Categorical(["a", "b"])}
                ),
                True,
            ),
            # differ types
            (DataFrame({"A": [1, 2], "B": [1.0, 2.0]}), False),
            # differ sizes
            (
                DataFrame(
                    {
                        "A": np.array([1, 2], dtype=np.int32),
                        "B": np.array([1, 2], dtype=np.int64),
                    }
                ),
                False,
            ),
            # multi-extension differ
            (
                DataFrame(
                    {"A": pd.Categorical(["a", "b"]), "B": pd.Categorical(["b", "c"])}
                ),
                False,
            ),
        ],
    )
    def test_is_homogeneous_type(self, data, expected):
        assert data._is_homogeneous_type is expected

    def test_is_homogeneous_type_clears_cache(self):
        ser = pd.Series([1, 2, 3])
        df = ser.to_frame("A")
        df["B"] = ser

        assert len(df._mgr.blocks) == 2

        a = df["B"]  # caches lookup
        df._is_homogeneous_type  # _should_ clear cache
        assert len(df._mgr.blocks) == 1
        assert df["B"] is not a

    def test_asarray_homogenous(self):
        df = pd.DataFrame({"A": pd.Categorical([1, 2]), "B": pd.Categorical([1, 2])})
        result = np.asarray(df)
        # may change from object in the future
        expected = np.array([[1, 1], [2, 2]], dtype="object")
        tm.assert_numpy_array_equal(result, expected)

    def test_str_to_small_float_conversion_type(self):
        # GH 20388
        np.random.seed(13)
        col_data = [str(np.random.random() * 1e-12) for _ in range(5)]
        result = pd.DataFrame(col_data, columns=["A"])
        expected = pd.DataFrame(col_data, columns=["A"], dtype=object)
        tm.assert_frame_equal(result, expected)
        # change the dtype of the elements from object to float one by one
        result.loc[result.index, "A"] = [float(x) for x in col_data]
        expected = pd.DataFrame(col_data, columns=["A"], dtype=float)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "convert_integer, expected", [(False, np.dtype("int32")), (True, "Int32")]
    )
    def test_convert_dtypes(self, convert_integer, expected):
        # Specific types are tested in tests/series/test_dtypes.py
        # Just check that it works for DataFrame here
        df = pd.DataFrame(
            {
                "a": pd.Series([1, 2, 3], dtype=np.dtype("int32")),
                "b": pd.Series(["x", "y", "z"], dtype=np.dtype("O")),
            }
        )
        result = df.convert_dtypes(True, True, convert_integer, False)
        expected = pd.DataFrame(
            {
                "a": pd.Series([1, 2, 3], dtype=expected),
                "b": pd.Series(["x", "y", "z"], dtype="string"),
            }
        )
        tm.assert_frame_equal(result, expected)


class TestDataFrameDatetimeWithTZ:
    def test_interleave(self, timezone_frame):

        # interleave with object
        result = timezone_frame.assign(D="foo").values
        expected = np.array(
            [
                [
                    Timestamp("2013-01-01 00:00:00"),
                    Timestamp("2013-01-02 00:00:00"),
                    Timestamp("2013-01-03 00:00:00"),
                ],
                [
                    Timestamp("2013-01-01 00:00:00-0500", tz="US/Eastern"),
                    pd.NaT,
                    Timestamp("2013-01-03 00:00:00-0500", tz="US/Eastern"),
                ],
                [
                    Timestamp("2013-01-01 00:00:00+0100", tz="CET"),
                    pd.NaT,
                    Timestamp("2013-01-03 00:00:00+0100", tz="CET"),
                ],
                ["foo", "foo", "foo"],
            ],
            dtype=object,
        ).T
        tm.assert_numpy_array_equal(result, expected)

        # interleave with only datetime64[ns]
        result = timezone_frame.values
        expected = np.array(
            [
                [
                    Timestamp("2013-01-01 00:00:00"),
                    Timestamp("2013-01-02 00:00:00"),
                    Timestamp("2013-01-03 00:00:00"),
                ],
                [
                    Timestamp("2013-01-01 00:00:00-0500", tz="US/Eastern"),
                    pd.NaT,
                    Timestamp("2013-01-03 00:00:00-0500", tz="US/Eastern"),
                ],
                [
                    Timestamp("2013-01-01 00:00:00+0100", tz="CET"),
                    pd.NaT,
                    Timestamp("2013-01-03 00:00:00+0100", tz="CET"),
                ],
            ],
            dtype=object,
        ).T
        tm.assert_numpy_array_equal(result, expected)
