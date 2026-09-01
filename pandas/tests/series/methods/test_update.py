import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm


class TestUpdate:
    def test_update(self):
        s = pd.Series([1.5, np.nan, 3.0, 4.0, np.nan])
        s2 = pd.Series([np.nan, 3.5, np.nan, 5.0])
        s.update(s2)

        expected = pd.Series([1.5, 3.5, 3.0, 5.0, np.nan])
        tm.assert_series_equal(s, expected)

        # GH 3217
        df = pd.DataFrame([{"a": 1}, {"a": 3, "b": 2}])
        df["c"] = np.nan
        # Cast to object to avoid upcast when setting "foo"
        df["c"] = df["c"].astype(object)
        df_orig = df.copy()

        with tm.raises_chained_assignment_error():
            df["c"].update(pd.Series(["foo"], index=[0]))
        expected = df_orig
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize(
        "other, dtype, expected, raises",
        [
            # other is int
            ([61, 63], "int32", pd.Series([10, 61, 12], dtype="int32"), False),
            ([61, 63], "int64", pd.Series([10, 61, 12]), False),
            ([61, 63], float, pd.Series([10.0, 61.0, 12.0]), False),
            ([61, 63], object, pd.Series([10, 61, 12], dtype=object), False),
            # other is float, but can be cast to int
            ([61.0, 63.0], "int32", pd.Series([10, 61, 12], dtype="int32"), False),
            ([61.0, 63.0], "int64", pd.Series([10, 61, 12]), False),
            ([61.0, 63.0], float, pd.Series([10.0, 61.0, 12.0]), False),
            ([61.0, 63.0], object, pd.Series([10, 61.0, 12], dtype=object), False),
            # others is float, cannot be cast to int
            ([61.1, 63.1], "int32", pd.Series([10.0, 61.1, 12.0]), True),
            ([61.1, 63.1], "int64", pd.Series([10.0, 61.1, 12.0]), True),
            ([61.1, 63.1], float, pd.Series([10.0, 61.1, 12.0]), False),
            ([61.1, 63.1], object, pd.Series([10, 61.1, 12], dtype=object), False),
            # other is object, cannot be cast
            ([(61,), (63,)], "int32", pd.Series([10, (61,), 12]), True),
            ([(61,), (63,)], "int64", pd.Series([10, (61,), 12]), True),
            ([(61,), (63,)], float, pd.Series([10.0, (61,), 12.0]), True),
            ([(61,), (63,)], object, pd.Series([10, (61,), 12]), False),
        ],
    )
    def test_update_dtypes(self, other, dtype, expected, raises):
        ser = pd.Series([10, 11, 12], dtype=dtype)
        other = pd.Series(other, index=[1, 3])
        if raises:
            with pytest.raises(TypeError, match="Invalid value"):
                ser.update(other)
        else:
            ser.update(other)
            tm.assert_series_equal(ser, expected)

    @pytest.mark.parametrize(
        "values, other, expected",
        [
            # update by key
            (
                {"a": 1, "b": 2, "c": 3, "d": 4},
                {"b": 5, "c": np.nan},
                {"a": 1, "b": 5, "c": 3, "d": 4},
            ),
            # update by position
            ([1, 2, 3, 4], [np.nan, 5, 1], [1, 5, 1, 4]),
        ],
    )
    def test_update_from_non_series(self, values, other, expected):
        # GH 33215
        series = pd.Series(values)
        series.update(other)
        expected = pd.Series(expected)
        tm.assert_series_equal(series, expected)

    @pytest.mark.parametrize(
        "data, other, expected, dtype",
        [
            (["a", None], [None, "b"], ["a", "b"], "string[python]"),
            pytest.param(
                ["a", None],
                [None, "b"],
                ["a", "b"],
                "string[pyarrow]",
                marks=td.skip_if_no("pyarrow"),
            ),
            ([1, None], [None, 2], [1, 2], "Int64"),
            ([True, None], [None, False], [True, False], "boolean"),
            (
                ["a", None],
                [None, "b"],
                ["a", "b"],
                pd.CategoricalDtype(categories=["a", "b"]),
            ),
            (
                [pd.Timestamp(year=2020, month=1, day=1, tz="Europe/London"), pd.NaT],
                [pd.NaT, pd.Timestamp(year=2020, month=1, day=1, tz="Europe/London")],
                [pd.Timestamp(year=2020, month=1, day=1, tz="Europe/London")] * 2,
                "datetime64[ns, Europe/London]",
            ),
        ],
    )
    def test_update_extension_array_series(self, data, other, expected, dtype):
        result = pd.Series(data, dtype=dtype)
        other = pd.Series(other, dtype=dtype)
        expected = pd.Series(expected, dtype=dtype)

        result.update(other)
        tm.assert_series_equal(result, expected)

    def test_update_with_categorical_type(self):
        # GH 25744
        dtype = pd.CategoricalDtype(["a", "b", "c", "d"])
        s1 = pd.Series(["a", "b", "c"], index=[1, 2, 3], dtype=dtype)
        s2 = pd.Series(["b", "a"], index=[1, 2], dtype=dtype)
        s1.update(s2)
        result = s1
        expected = pd.Series(["b", "a", "c"], index=[1, 2, 3], dtype=dtype)
        tm.assert_series_equal(result, expected)
