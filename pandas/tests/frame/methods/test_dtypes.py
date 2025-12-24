import numpy as np
import pytest

pa = pytest.importorskip("pyarrow")

from pandas.core.dtypes.dtypes import DatetimeTZDtype

import pandas as pd
from pandas import (
    DataFrame,
    Series,
    date_range,
)
import pandas._testing as tm


class TestDataFrameDataTypes:
    def test_empty_frame_dtypes(self):
        empty_df = DataFrame()
        tm.assert_series_equal(empty_df.dtypes, Series(dtype=object))

        nocols_df = DataFrame(index=[1, 2, 3])
        tm.assert_series_equal(nocols_df.dtypes, Series(dtype=object))

        norows_df = DataFrame(columns=list("abc"))
        tm.assert_series_equal(norows_df.dtypes, Series(object, index=list("abc")))

        norows_int_df = DataFrame(columns=list("abc")).astype(np.int32)
        tm.assert_series_equal(
            norows_int_df.dtypes, Series(np.dtype("int32"), index=list("abc"))
        )

        df = DataFrame({"a": 1, "b": True, "c": 1.0}, index=[1, 2, 3])
        ex_dtypes = Series({"a": np.int64, "b": np.bool_, "c": np.float64})
        tm.assert_series_equal(df.dtypes, ex_dtypes)

        # same but for empty slice of df
        tm.assert_series_equal(df[:0].dtypes, ex_dtypes)

    def test_datetime_with_tz_dtypes(self):
        tzframe = DataFrame(
            {
                "A": date_range("20130101", periods=3, unit="ns"),
                "B": date_range("20130101", periods=3, tz="US/Eastern", unit="ns"),
                "C": date_range("20130101", periods=3, tz="CET", unit="ns"),
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
        df = DataFrame(index=range(5), columns=list("abc"), dtype=np.float64)
        tm.assert_series_equal(
            df.dtypes,
            Series({"a": np.float64, "b": np.float64, "c": np.float64}),
        )
        tm.assert_series_equal(df.iloc[:, 2:].dtypes, Series({"c": np.float64}))
        tm.assert_series_equal(
            df.dtypes,
            Series({"a": np.float64, "b": np.float64, "c": np.float64}),
        )

    @pytest.mark.parametrize(
        "data",
        [pd.NA, True],
    )
    def test_dtypes_are_correct_after_groupby_last(self, data):
        # GH46409
        df = DataFrame(
            {"id": [1, 2, 3, 4], "test": [True, pd.NA, data, False]}
        ).convert_dtypes()
        result = df.groupby("id").last().test
        expected = df.set_index("id").test
        assert result.dtype == pd.BooleanDtype()
        tm.assert_series_equal(expected, result)

    def test_dtypes_gh8722(self, float_string_frame):
        float_string_frame["bool"] = float_string_frame["A"] > 0
        result = float_string_frame.dtypes
        expected = Series(
            {k: v.dtype for k, v in float_string_frame.items()}, index=result.index
        )
        tm.assert_series_equal(result, expected)

    def test_dtypes_timedeltas(self):
        df = DataFrame(
            {
                "A": Series(date_range("2012-1-1", periods=3, freq="D", unit="ns")),
                "B": Series([timedelta(days=i) for i in range(3)]),
            }
        )
        result = df.dtypes
        expected = Series(
            [np.dtype("datetime64[ns]"), np.dtype("timedelta64[us]")], index=list("AB")
        )
        tm.assert_series_equal(result, expected)

        df["C"] = df["A"] + df["B"]
        result = df.dtypes
        expected = Series(
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
        expected = Series(
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
        df = DataFrame([["foo"]])
        result = df.apply(lambda col: np.array("bar"))
        expected = Series(np.array("bar"))
        tm.assert_series_equal(result, expected)

    def test_apply_axis_columns_pyarrow_backend_does_not_preserve_arrowdtype(self):
    # GH#63317
        df = DataFrame(
        {
            "year": [2020, 2021, 2022],
            "make": ["Ford", "Toyota", "Honda"],
            "city08": [20, 30, 25],
            "highway08": [30, 40, 35],
        }
        ).convert_dtypes(dtype_backend="pyarrow")

        result = df.apply(
        lambda row: (row["city08"] + row["highway08"]) / 2,
        axis="columns",
        )

    # Row-wise UDF apply returns Python / NumPy scalars, so the result
    # dtype is inferred from those scalars and does not preserve the
    # pyarrow-backed dtype backend.
        assert result.dtype == "float64"
        assert not isinstance(result.dtype, pd.ArrowDtype)


    def test_apply_named_reduction_pyarrow_backend_preserved(self):
    # GH#63317 - contrast with a known reduction
        df = DataFrame(
        {
            "city08": [20, 30, 25],
            "highway08": [30, 40, 35],
        }
        ).convert_dtypes(dtype_backend="pyarrow")

        result = df.apply("mean", axis="columns")

    # For a named reduction, pandas knows the operation ahead of time
    # and can keep a pyarrow-backed result.
        assert isinstance(result.dtype, pd.ArrowDtype)
