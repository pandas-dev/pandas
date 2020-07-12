import numpy as np
import pytest

from pandas import CategoricalDtype, DataFrame, NaT, Series, Timestamp
import pandas._testing as tm


class TestUpdate:
    def test_update(self):
        s = Series([1.5, np.nan, 3.0, 4.0, np.nan])
        s2 = Series([np.nan, 3.5, np.nan, 5.0])
        s.update(s2)

        expected = Series([1.5, 3.5, 3.0, 5.0, np.nan])
        tm.assert_series_equal(s, expected)

        # GH 3217
        df = DataFrame([{"a": 1}, {"a": 3, "b": 2}])
        df["c"] = np.nan

        df["c"].update(Series(["foo"], index=[0]))
        expected = DataFrame(
            [[1, np.nan, "foo"], [3, 2.0, np.nan]], columns=["a", "b", "c"]
        )
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize(
        "other, dtype, expected",
        [
            # other is int
            ([61, 63], "int32", Series([10, 61, 12], dtype="int32")),
            ([61, 63], "int64", Series([10, 61, 12])),
            ([61, 63], float, Series([10.0, 61.0, 12.0])),
            ([61, 63], object, Series([10, 61, 12], dtype=object)),
            # other is float, but can be cast to int
            ([61.0, 63.0], "int32", Series([10, 61, 12], dtype="int32")),
            ([61.0, 63.0], "int64", Series([10, 61, 12])),
            ([61.0, 63.0], float, Series([10.0, 61.0, 12.0])),
            ([61.0, 63.0], object, Series([10, 61.0, 12], dtype=object)),
            # others is float, cannot be cast to int
            ([61.1, 63.1], "int32", Series([10.0, 61.1, 12.0])),
            ([61.1, 63.1], "int64", Series([10.0, 61.1, 12.0])),
            ([61.1, 63.1], float, Series([10.0, 61.1, 12.0])),
            ([61.1, 63.1], object, Series([10, 61.1, 12], dtype=object)),
            # other is object, cannot be cast
            ([(61,), (63,)], "int32", Series([10, (61,), 12])),
            ([(61,), (63,)], "int64", Series([10, (61,), 12])),
            ([(61,), (63,)], float, Series([10.0, (61,), 12.0])),
            ([(61,), (63,)], object, Series([10, (61,), 12])),
        ],
    )
    def test_update_dtypes(self, other, dtype, expected):

        ser = Series([10, 11, 12], dtype=dtype)
        other = Series(other, index=[1, 3])
        ser.update(other)

        tm.assert_series_equal(ser, expected)

    @pytest.mark.parametrize(
        "series, other, expected",
        [
            # update by key
            (
                Series({"a": 1, "b": 2, "c": 3, "d": 4}),
                {"b": 5, "c": np.nan},
                Series({"a": 1, "b": 5, "c": 3, "d": 4}),
            ),
            # update by position
            (Series([1, 2, 3, 4]), [np.nan, 5, 1], Series([1, 5, 1, 4])),
        ],
    )
    def test_update_from_non_series(self, series, other, expected):
        # GH 33215
        series.update(other)
        tm.assert_series_equal(series, expected)

    @pytest.mark.parametrize(
        "result, target, expected",
        [
            (
                Series(["a", None], dtype="string"),
                Series([None, "b"], dtype="string"),
                Series(["a", "b"], dtype="string"),
            ),
            (
                Series([1, None], dtype="Int64"),
                Series([None, 2], dtype="Int64"),
                Series([1, 2], dtype="Int64"),
            ),
            (
                Series([True, None], dtype="boolean"),
                Series([None, False], dtype="boolean"),
                Series([True, False], dtype="boolean"),
            ),
            (
                Series(["a", None], dtype=CategoricalDtype(categories=["a", "b"])),
                Series([None, "b"], dtype=CategoricalDtype(categories=["a", "b"])),
                Series(["a", "b"], dtype=CategoricalDtype(categories=["a", "b"])),
            ),
            (
                Series([Timestamp(year=2020, month=1, day=1, tz="Europe/London"), NaT]),
                Series([NaT, Timestamp(year=2020, month=1, day=1, tz="Europe/London")]),
                Series([Timestamp(year=2020, month=1, day=1, tz="Europe/London")] * 2),
            ),
        ],
    )
    def test_update_extension_array_series(self, result, target, expected):
        result.update(target)
        tm.assert_series_equal(result, expected)
