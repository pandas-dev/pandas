import numpy as np
import pytest

from pandas import DataFrame, Series
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
