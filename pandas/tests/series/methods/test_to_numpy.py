import numpy as np
import pytest

from pandas.errors import PerformanceWarning
import pandas.util._test_decorators as td

from pandas import (
    NA,
    DatetimeIndex,
    Series,
    Timedelta,
    date_range,
)
import pandas._testing as tm


@pytest.mark.parametrize("dtype", ["int64", "float64"])
def test_to_numpy_na_value(dtype):
    # GH#48951
    ser = Series([1, 2, NA, 4])
    result = ser.to_numpy(dtype=dtype, na_value=0)
    expected = np.array([1, 2, 0, 4], dtype=dtype)
    tm.assert_numpy_array_equal(result, expected)


def test_to_numpy_cast_before_setting_na():
    # GH#50600
    ser = Series([1])
    result = ser.to_numpy(dtype=np.float64, na_value=np.nan)
    expected = np.array([1.0])
    tm.assert_numpy_array_equal(result, expected)


def test_to_numpy_copy_false_returns_readonly_view():
    # GH#57431 - to_numpy(copy=False) should return a read-only view
    ser = Series([1.0, 2.0, 3.0])
    result = ser.to_numpy(copy=False)
    assert result.flags.writeable is False
    assert np.shares_memory(result, ser.to_numpy(copy=False))


@td.skip_if_no("pyarrow")
def test_to_numpy_arrow_dtype_given():
    # GH#57121
    ser = Series([1, NA], dtype="int64[pyarrow]")
    result = ser.to_numpy(dtype="float64")
    expected = np.array([1.0, np.nan])
    tm.assert_numpy_array_equal(result, expected)


def test_astype_ea_int_to_td_ts():
    # GH#57093
    ser = Series([1, None], dtype="Int64")
    result = ser.astype("m8[ns]")
    expected = Series([1, Timedelta("nat")], dtype="m8[ns]")
    tm.assert_series_equal(result, expected)

    result = ser.astype("M8[ns]")
    expected = Series([1, Timedelta("nat")], dtype="M8[ns]")
    tm.assert_series_equal(result, expected)


class TestToNumpyPerformanceWarning:
    # GH#62309
    def test_to_numpy_tz_aware_warns(self):
        # to_numpy() without dtype on tz-aware data should warn
        dti = date_range("2024-01-01", periods=5, freq="h", tz="UTC")
        ser = Series(dti)

        with tm.assert_produces_warning(PerformanceWarning):
            result = ser.to_numpy()
        assert result.dtype == np.dtype(object)

        with tm.assert_produces_warning(PerformanceWarning):
            result = DatetimeIndex(dti).to_numpy()
        assert result.dtype == np.dtype(object)

    def test_to_numpy_tz_aware_explicit_dtype_no_warn(self):
        # Specifying dtype explicitly should not warn
        dti = date_range("2024-01-01", periods=5, freq="h", tz="UTC")
        ser = Series(dti)

        with tm.assert_produces_warning(None):
            result = ser.to_numpy(dtype="datetime64[ns]")
        assert result.dtype == np.dtype("datetime64[ns]")

        with tm.assert_produces_warning(None):
            result = ser.to_numpy(dtype=object)
        assert result.dtype == np.dtype(object)

    def test_to_numpy_tz_naive_no_warn(self):
        # tz-naive datetime data should not warn
        dti = date_range("2024-01-01", periods=5, freq="h")
        ser = Series(dti)

        with tm.assert_produces_warning(None):
            result = ser.to_numpy()
        assert result.dtype.kind == "M"
