import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Series, Timestamp, date_range
import pandas._testing as tm


class TestDataFrameDiff:
    def test_diff(self, datetime_frame):
        the_diff = datetime_frame.diff(1)

        tm.assert_series_equal(
            the_diff["A"], datetime_frame["A"] - datetime_frame["A"].shift(1)
        )

        # int dtype
        a = 10_000_000_000_000_000
        b = a + 1
        s = Series([a, b])

        rs = DataFrame({"s": s}).diff()
        assert rs.s[1] == 1

        # mixed numeric
        tf = datetime_frame.astype("float32")
        the_diff = tf.diff(1)
        tm.assert_series_equal(the_diff["A"], tf["A"] - tf["A"].shift(1))

        # GH#10907
        df = pd.DataFrame({"y": pd.Series([2]), "z": pd.Series([3])})
        df.insert(0, "x", 1)
        result = df.diff(axis=1)
        expected = pd.DataFrame(
            {"x": np.nan, "y": pd.Series(1), "z": pd.Series(1)}
        ).astype("float64")
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("tz", [None, "UTC"])
    def test_diff_datetime_axis0(self, tz):
        # GH#18578
        df = DataFrame(
            {
                0: date_range("2010", freq="D", periods=2, tz=tz),
                1: date_range("2010", freq="D", periods=2, tz=tz),
            }
        )

        result = df.diff(axis=0)
        expected = DataFrame(
            {
                0: pd.TimedeltaIndex(["NaT", "1 days"]),
                1: pd.TimedeltaIndex(["NaT", "1 days"]),
            }
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("tz", [None, "UTC"])
    def test_diff_datetime_axis1(self, tz):
        # GH#18578
        df = DataFrame(
            {
                0: date_range("2010", freq="D", periods=2, tz=tz),
                1: date_range("2010", freq="D", periods=2, tz=tz),
            }
        )

        result = df.diff(axis=1)
        expected = DataFrame(
            {
                0: pd.TimedeltaIndex(["NaT", "NaT"]),
                1: pd.TimedeltaIndex(["0 days", "0 days"]),
            }
        )
        tm.assert_frame_equal(result, expected)

    def test_diff_timedelta(self):
        # GH#4533
        df = DataFrame(
            dict(
                time=[Timestamp("20130101 9:01"), Timestamp("20130101 9:02")],
                value=[1.0, 2.0],
            )
        )

        res = df.diff()
        exp = DataFrame(
            [[pd.NaT, np.nan], [pd.Timedelta("00:01:00"), 1]], columns=["time", "value"]
        )
        tm.assert_frame_equal(res, exp)

    def test_diff_mixed_dtype(self):
        df = DataFrame(np.random.randn(5, 3))
        df["A"] = np.array([1, 2, 3, 4, 5], dtype=object)

        result = df.diff()
        assert result[0].dtype == np.float64

    def test_diff_neg_n(self, datetime_frame):
        rs = datetime_frame.diff(-1)
        xp = datetime_frame - datetime_frame.shift(-1)
        tm.assert_frame_equal(rs, xp)

    def test_diff_float_n(self, datetime_frame):
        rs = datetime_frame.diff(1.0)
        xp = datetime_frame.diff(1)
        tm.assert_frame_equal(rs, xp)

    def test_diff_axis(self):
        # GH#9727
        df = DataFrame([[1.0, 2.0], [3.0, 4.0]])
        tm.assert_frame_equal(
            df.diff(axis=1), DataFrame([[np.nan, 1.0], [np.nan, 1.0]])
        )
        tm.assert_frame_equal(
            df.diff(axis=0), DataFrame([[np.nan, np.nan], [2.0, 2.0]])
        )

    @pytest.mark.xfail(
        reason="GH#32995 needs to operate column-wise or do inference",
        raises=AssertionError,
    )
    def test_diff_period(self):
        # GH#32995 Don't pass an incorrect axis
        #  TODO(EA2D): this bug wouldn't have happened with 2D EA
        pi = pd.date_range("2016-01-01", periods=3).to_period("D")
        df = pd.DataFrame({"A": pi})

        result = df.diff(1, axis=1)

        # TODO: should we make Block.diff do type inference?  or maybe algos.diff?
        expected = (df - pd.NaT).astype(object)
        tm.assert_frame_equal(result, expected)

    def test_diff_axis1_mixed_dtypes(self):
        # GH#32995 operate column-wise when we have mixed dtypes and axis=1
        df = pd.DataFrame({"A": range(3), "B": 2 * np.arange(3, dtype=np.float64)})

        expected = pd.DataFrame({"A": [np.nan, np.nan, np.nan], "B": df["B"] / 2})

        result = df.diff(axis=1)
        tm.assert_frame_equal(result, expected)

    def test_diff_axis1_mixed_dtypes_large_periods(self):
        # GH#32995 operate column-wise when we have mixed dtypes and axis=1
        df = pd.DataFrame({"A": range(3), "B": 2 * np.arange(3, dtype=np.float64)})

        expected = df * np.nan

        result = df.diff(axis=1, periods=3)
        tm.assert_frame_equal(result, expected)

    def test_diff_axis1_mixed_dtypes_negative_periods(self):
        # GH#32995 operate column-wise when we have mixed dtypes and axis=1
        df = pd.DataFrame({"A": range(3), "B": 2 * np.arange(3, dtype=np.float64)})

        expected = pd.DataFrame({"A": -1.0 * df["A"], "B": df["B"] * np.nan})

        result = df.diff(axis=1, periods=-1)
        tm.assert_frame_equal(result, expected)

    def test_diff_sparse(self):
        # GH#28813 .diff() should work for sparse dataframes as well
        sparse_df = pd.DataFrame([[0, 1], [1, 0]], dtype="Sparse[int]")

        result = sparse_df.diff()
        expected = pd.DataFrame(
            [[np.nan, np.nan], [1.0, -1.0]], dtype=pd.SparseDtype("float", 0.0)
        )

        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "axis,expected",
        [
            (
                0,
                pd.DataFrame(
                    {
                        "a": [np.nan, 0, 1, 0, np.nan, np.nan, np.nan, 0],
                        "b": [np.nan, 1, np.nan, np.nan, -2, 1, np.nan, np.nan],
                        "c": np.repeat(np.nan, 8),
                        "d": [np.nan, 3, 5, 7, 9, 11, 13, 15],
                    },
                    dtype="Int64",
                ),
            ),
            (
                1,
                pd.DataFrame(
                    {
                        "a": np.repeat(np.nan, 8),
                        "b": [0, 1, np.nan, 1, np.nan, np.nan, np.nan, 0],
                        "c": np.repeat(np.nan, 8),
                        "d": np.repeat(np.nan, 8),
                    },
                    dtype="Int64",
                ),
            ),
        ],
    )
    def test_diff_integer_na(self, axis, expected):
        # GH#24171 IntegerNA Support for DataFrame.diff()
        df = pd.DataFrame(
            {
                "a": np.repeat([0, 1, np.nan, 2], 2),
                "b": np.tile([0, 1, np.nan, 2], 2),
                "c": np.repeat(np.nan, 8),
                "d": np.arange(1, 9) ** 2,
            },
            dtype="Int64",
        )

        # Test case for default behaviour of diff
        result = df.diff(axis=axis)
        tm.assert_frame_equal(result, expected)

    def test_diff_readonly(self):
        # https://github.com/pandas-dev/pandas/issues/35559
        arr = np.random.randn(5, 2)
        arr.flags.writeable = False
        df = pd.DataFrame(arr)
        result = df.diff()
        expected = pd.DataFrame(np.array(df)).diff()
        tm.assert_frame_equal(result, expected)
