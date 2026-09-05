import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestSeriesDiff:
    def test_diff_series_requires_integer(self):
        series = pd.Series(np.random.default_rng(2).standard_normal(2))
        with pytest.raises(ValueError, match="periods must be an integer"):
            series.diff(1.5)

    def test_diff_np(self):
        # TODO(__array_function__): could make np.diff return a Series
        #  matching ser.diff()

        ser = pd.Series(np.arange(5))

        res = np.diff(ser)
        expected = np.array([1, 1, 1, 1])
        tm.assert_numpy_array_equal(res, expected)

    def test_diff_int(self):
        # int dtype
        a = 10000000000000000
        b = a + 1
        ser = pd.Series([a, b])

        result = ser.diff()
        assert result[1] == 1

    def test_diff_tz(self):
        # Combined datetime diff, normal diff and boolean diff test
        ts = pd.Series(
            np.arange(10, dtype=np.float64),
            index=pd.date_range("2020-01-01", periods=10),
            name="ts",
        )
        ts.diff()

        # neg n
        result = ts.diff(-1)
        expected = ts - ts.shift(-1)
        tm.assert_series_equal(result, expected)

        # 0
        result = ts.diff(0)
        expected = ts - ts
        tm.assert_series_equal(result, expected)

    def test_diff_dt64(self):
        # datetime diff (GH#3100)
        ser = pd.Series(pd.date_range("20130102", periods=5))
        result = ser.diff()
        expected = ser - ser.shift(1)
        tm.assert_series_equal(result, expected)

        # timedelta diff
        result = result - result.shift(1)  # previous result
        expected = expected.diff()  # previously expected
        tm.assert_series_equal(result, expected)

    def test_diff_dt64tz(self):
        # with tz
        ser = pd.Series(
            pd.date_range("2000-01-01 09:00:00", periods=5, tz="US/Eastern"),
            name="foo",
        )
        result = ser.diff()
        expected = pd.Series(pd.TimedeltaIndex(["NaT"] + ["1 days"] * 4), name="foo")
        tm.assert_series_equal(result, expected)

    def test_diff_bool(self):
        # boolean series (test for fixing #17294)
        data = [False, True, True, False, False]
        output = [np.nan, True, False, True, False]
        ser = pd.Series(data)
        result = ser.diff()
        expected = pd.Series(output)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("dtype", ["uint8", "uint16", "uint32", "uint64"])
    def test_diff_unsigned_int_no_overflow(self, dtype):
        # GH#4899 unsigned subtraction used to wrap around instead of
        # producing a negative difference, e.g. uint32(0) - uint32(3)
        # silently became 4294967293 instead of -3.
        ser = pd.Series(np.array([3, 2, 0], dtype=dtype))
        result = ser.diff()
        expected = pd.Series([np.nan, -1, -2])
        tm.assert_series_equal(result, expected, check_dtype=False)
        assert (result.dropna() < 0).all()

    def test_diff_uint64_overflows_int64(self):
        # GH#4899 uint64 values that don't fit in int64 must not be
        # silently truncated; diff exactly using arbitrary-precision ints.
        ser = pd.Series(np.array([2**63 + 10, 2**63 + 5, 5], dtype="uint64"))
        result = ser.diff()
        expected = pd.Series([np.nan, -5, 5 - (2**63 + 5)], dtype=object)
        tm.assert_series_equal(result, expected)

    def test_diff_object_dtype(self):
        # object series
        ser = pd.Series([False, True, 5.0, np.nan, True, False])
        result = ser.diff()
        expected = ser - ser.shift(1)
        tm.assert_series_equal(result, expected)
