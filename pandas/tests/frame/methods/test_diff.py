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
        if tz is None:
            result = df.diff(axis=1)
            expected = DataFrame(
                {
                    0: pd.TimedeltaIndex(["NaT", "NaT"]),
                    1: pd.TimedeltaIndex(["0 days", "0 days"]),
                }
            )
            tm.assert_frame_equal(result, expected)
        else:
            with pytest.raises(NotImplementedError):
                result = df.diff(axis=1)

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
