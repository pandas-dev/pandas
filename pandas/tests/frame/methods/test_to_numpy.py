import numpy as np
import pytest

from pandas import (
    DataFrame,
    NaT,
    Timestamp,
    date_range,
)
import pandas._testing as tm


class TestToNumpy:
    def test_to_numpy(self):
        df = DataFrame({"A": [1, 2], "B": [3, 4.5]})
        expected = np.array([[1, 3], [2, 4.5]])
        result = df.to_numpy()
        tm.assert_numpy_array_equal(result, expected)

    def test_to_numpy_dtype(self):
        df = DataFrame({"A": [1, 2], "B": [3, 4.5]})
        expected = np.array([[1, 3], [2, 4]], dtype="int64")
        result = df.to_numpy(dtype="int64")
        tm.assert_numpy_array_equal(result, expected)

    def test_to_numpy_copy(self):
        arr = np.random.default_rng(2).standard_normal((4, 3))
        df = DataFrame(arr)
        assert df.values.base is not arr
        assert df.to_numpy(copy=False).base is df.values.base
        assert df.to_numpy(copy=True).base is not arr

        # we still don't want a copy when na_value=np.nan is passed,
        #  and that can be respected because we are already numpy-float
        assert df.to_numpy(copy=False).base is df.values.base

    @pytest.mark.filterwarnings(
        "ignore:invalid value encountered in cast:RuntimeWarning"
    )
    def test_to_numpy_mixed_dtype_to_str(self):
        # https://github.com/pandas-dev/pandas/issues/35455
        df = DataFrame([[Timestamp("2020-01-01 00:00:00"), 100.0]])
        result = df.to_numpy(dtype=str)
        expected = np.array([["2020-01-01 00:00:00", "100.0"]], dtype=str)
        tm.assert_numpy_array_equal(result, expected)

    def test_to_numpy_datetime_with_na(self):
        # GH #53115
        dti = date_range("2016-01-01", periods=3, unit="ns")
        df = DataFrame(dti)
        df.iloc[0, 0] = NaT
        expected = np.array([[np.nan], [1.45169280e18], [1.45177920e18]])
        result = df.to_numpy(float, na_value=np.nan)
        tm.assert_numpy_array_equal(result, expected)

        df = DataFrame(
            {
                "a": [
                    Timestamp("1970-01-01").as_unit("s"),
                    Timestamp("1970-01-02").as_unit("s"),
                    NaT,
                ],
                "b": [
                    Timestamp("1970-01-01").as_unit("s"),
                    np.nan,
                    Timestamp("1970-01-02").as_unit("s"),
                ],
                "c": [
                    1,
                    np.nan,
                    2,
                ],
            }
        )
        expected = np.array(
            [
                [0.00e00, 0.00e00, 1.00e00],
                [8.64e04, np.nan, np.nan],
                [np.nan, 8.64e04, 2.00e00],
            ]
        )
        result = df.to_numpy(float, na_value=np.nan)
        tm.assert_numpy_array_equal(result, expected)
