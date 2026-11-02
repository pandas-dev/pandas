"""
Tests for np.foo applied to Series, not necessarily ufuncs.
"""

import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm


class TestPtp:
    def test_ptp(self):
        # GH#21614
        N = 1000
        arr = np.random.default_rng(2).standard_normal(N)
        ser = pd.Series(arr)
        assert np.ptp(ser) == np.ptp(arr)

    def test_ptp_datetime(self):
        ser = pd.Series(pd.date_range("2020", periods=3))
        result = np.ptp(ser)
        assert result == pd.Timedelta(days=2)
        assert type(result) is pd.Timedelta

    def test_ptp_frame(self):
        df = pd.DataFrame({"a": [1.0, 5.0], "b": [2.0, 0.0]})
        result = np.ptp(df, axis=0)
        expected = pd.Series([4.0, 2.0], index=["a", "b"])
        tm.assert_series_equal(result, expected)


def test_numpy_unique(datetime_series):
    # it works!
    np.unique(datetime_series)


@pytest.mark.parametrize("index", [["a", "b", "c", "d", "e"], None])
def test_numpy_argwhere(index):
    # GH#35331

    s = pd.Series(range(5), index=index, dtype=np.int64)

    result = np.argwhere(s > 2).astype(np.int64)
    expected = np.array([[3], [4]], dtype=np.int64)

    tm.assert_numpy_array_equal(result, expected)


@td.skip_if_no("pyarrow")
def test_log_arrow_backed_missing_value(using_nan_is_na):
    # GH#56285, GH#62506
    ser = pd.Series([1, 2, None], dtype="float64[pyarrow]")
    result = np.log(ser)
    expected = np.log(pd.Series([1, 2, None], dtype="float64[pyarrow]"))
    tm.assert_series_equal(result, expected)
