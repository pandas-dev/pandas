"""
Tests for DataFrame cumulative operations

See also
--------
tests.series.test_cumulative
"""

import numpy as np
import pytest

from pandas import NA, DataFrame, Series
import pandas._testing as tm


class TestDataFrameCumulativeOps:
    # ---------------------------------------------------------------------
    # Cumulative Operations - cumsum, cummax, ...

    def test_cumsum_corner(self):
        dm = DataFrame(np.arange(20).reshape(4, 5), index=range(4), columns=range(5))
        # TODO(wesm): do something with this?
        result = dm.cumsum()  # noqa

    def test_cumsum(self, datetime_frame):
        datetime_frame.iloc[5:10, 0] = np.nan
        datetime_frame.iloc[10:15, 1] = np.nan
        datetime_frame.iloc[15:, 2] = np.nan

        # axis = 0
        cumsum = datetime_frame.cumsum()
        expected = datetime_frame.apply(Series.cumsum)
        tm.assert_frame_equal(cumsum, expected)

        # axis = 1
        cumsum = datetime_frame.cumsum(axis=1)
        expected = datetime_frame.apply(Series.cumsum, axis=1)
        tm.assert_frame_equal(cumsum, expected)

        # works
        df = DataFrame({"A": np.arange(20)}, index=np.arange(20))
        df.cumsum()

        # fix issue
        cumsum_xs = datetime_frame.cumsum(axis=1)
        assert np.shape(cumsum_xs) == np.shape(datetime_frame)

    def test_cumprod(self, datetime_frame):
        datetime_frame.iloc[5:10, 0] = np.nan
        datetime_frame.iloc[10:15, 1] = np.nan
        datetime_frame.iloc[15:, 2] = np.nan

        # axis = 0
        cumprod = datetime_frame.cumprod()
        expected = datetime_frame.apply(Series.cumprod)
        tm.assert_frame_equal(cumprod, expected)

        # axis = 1
        cumprod = datetime_frame.cumprod(axis=1)
        expected = datetime_frame.apply(Series.cumprod, axis=1)
        tm.assert_frame_equal(cumprod, expected)

        # fix issue
        cumprod_xs = datetime_frame.cumprod(axis=1)
        assert np.shape(cumprod_xs) == np.shape(datetime_frame)

        # ints
        df = datetime_frame.fillna(0).astype(int)
        df.cumprod(0)
        df.cumprod(1)

        # ints32
        df = datetime_frame.fillna(0).astype(np.int32)
        df.cumprod(0)
        df.cumprod(1)

    def test_cummin(self, datetime_frame):
        datetime_frame.iloc[5:10, 0] = np.nan
        datetime_frame.iloc[10:15, 1] = np.nan
        datetime_frame.iloc[15:, 2] = np.nan

        # axis = 0
        cummin = datetime_frame.cummin()
        expected = datetime_frame.apply(Series.cummin)
        tm.assert_frame_equal(cummin, expected)

        # axis = 1
        cummin = datetime_frame.cummin(axis=1)
        expected = datetime_frame.apply(Series.cummin, axis=1)
        tm.assert_frame_equal(cummin, expected)

        # it works
        df = DataFrame({"A": np.arange(20)}, index=np.arange(20))
        df.cummin()

        # fix issue
        cummin_xs = datetime_frame.cummin(axis=1)
        assert np.shape(cummin_xs) == np.shape(datetime_frame)

    def test_cummax(self, datetime_frame):
        datetime_frame.iloc[5:10, 0] = np.nan
        datetime_frame.iloc[10:15, 1] = np.nan
        datetime_frame.iloc[15:, 2] = np.nan

        # axis = 0
        cummax = datetime_frame.cummax()
        expected = datetime_frame.apply(Series.cummax)
        tm.assert_frame_equal(cummax, expected)

        # axis = 1
        cummax = datetime_frame.cummax(axis=1)
        expected = datetime_frame.apply(Series.cummax, axis=1)
        tm.assert_frame_equal(cummax, expected)

        # it works
        df = DataFrame({"A": np.arange(20)}, index=np.arange(20))
        df.cummax()

        # fix issue
        cummax_xs = datetime_frame.cummax(axis=1)
        assert np.shape(cummax_xs) == np.shape(datetime_frame)

    def test_cumulative_ops_preserve_dtypes(self):
        # GH#19296 dont incorrectly upcast to object
        df = DataFrame({"A": [1, 2, 3], "B": [1, 2, 3.0], "C": [True, False, False]})

        result = df.cumsum()

        expected = DataFrame(
            {
                "A": Series([1, 3, 6], dtype=np.int64),
                "B": Series([1, 3, 6], dtype=np.float64),
                "C": df["C"].cumsum(),
            }
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "func, exp",
        [
            ("cumsum", [2, NA, 7, 6, 6]),
            ("cumprod", [2, NA, 10, -10, 0]),
            ("cummin", [2, NA, 2, -1, -1]),
            ("cummax", [2, NA, 5, 5, 5]),
        ],
    )
    @pytest.mark.parametrize("dtype", ["Float64", "Int64"])
    def test_cummulative_ops_extension_dtype(self, frame_or_series, dtype, func, exp):
        # GH#39479
        obj = frame_or_series([2, np.nan, 5, -1, 0], dtype=dtype)
        result = getattr(obj, func)()
        expected = frame_or_series(exp, dtype=dtype)
        tm.assert_equal(result, expected)
