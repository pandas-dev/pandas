import numpy as np
import pytest

from pandas.errors import Pandas4Warning
import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize("dtype", ["int64", "float64"])
def test_to_numpy_na_value(dtype):
    # GH#48951
    ser = pd.Series([1, 2, pd.NA, 4])
    result = ser.to_numpy(dtype=dtype, na_value=0)
    expected = np.array([1, 2, 0, 4], dtype=dtype)
    tm.assert_numpy_array_equal(result, expected)


def test_to_numpy_cast_before_setting_na():
    # GH#50600
    ser = pd.Series([1])
    result = ser.to_numpy(dtype=np.float64, na_value=np.nan)
    expected = np.array([1.0])
    tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize("tz", [None, "UTC"])
def test_to_numpy_unitless_datetime64_deprecated(tz):
    # GH#59772
    ser = pd.Series(pd.NaT)
    if tz is not None:
        ser = ser.dt.tz_localize(tz)

    msg = "Using a unit-less 'datetime64' dtype in to_numpy is deprecated"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = ser.to_numpy("datetime64")

    assert result.dtype == np.dtype("datetime64[s]")
    assert np.isnat(result).all()


def test_to_numpy_unitless_datetime64_non_nat_deprecated():
    # GH#59772
    ser = pd.Series([pd.Timestamp("2020-01-01")])

    msg = "Using a unit-less 'datetime64' dtype in to_numpy is deprecated"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = ser.to_numpy("datetime64")

    expected = np.array(["2020-01-01"], dtype="datetime64[us]")
    tm.assert_numpy_array_equal(result, expected)


def test_to_numpy_copy_false_returns_readonly_view():
    # GH#57431 - to_numpy(copy=False) should return a read-only view
    ser = pd.Series([1.0, 2.0, 3.0])
    result = ser.to_numpy(copy=False)
    assert result.flags.writeable is False
    assert np.shares_memory(result, ser.to_numpy(copy=False))


@td.skip_if_no("pyarrow")
def test_to_numpy_arrow_dtype_given():
    # GH#57121
    ser = pd.Series([1, pd.NA], dtype="int64[pyarrow]")
    result = ser.to_numpy(dtype="float64")
    expected = np.array([1.0, np.nan])
    tm.assert_numpy_array_equal(result, expected)


def test_astype_ea_int_to_td_ts():
    # GH#57093
    ser = pd.Series([1, None], dtype="Int64")
    result = ser.astype("m8[ns]")
    expected = pd.Series([1, pd.Timedelta("nat")], dtype="m8[ns]")
    tm.assert_series_equal(result, expected)

    result = ser.astype("M8[ns]")
    expected = pd.Series([1, pd.Timedelta("nat")], dtype="M8[ns]")
    tm.assert_series_equal(result, expected)
