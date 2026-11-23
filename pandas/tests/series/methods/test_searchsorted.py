import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm
from pandas.api.types import is_scalar


class TestSeriesSearchSorted:
    def test_searchsorted(self):
        ser = pd.Series([1, 2, 3])

        result = ser.searchsorted(1, side="left")
        assert is_scalar(result)
        assert result == 0

        result = ser.searchsorted(1, side="right")
        assert is_scalar(result)
        assert result == 1

    def test_searchsorted_numeric_dtypes_scalar(self):
        ser = pd.Series([1, 2, 90, 1000, 3e9])
        res = ser.searchsorted(30)
        assert is_scalar(res)
        assert res == 2

        res = ser.searchsorted([30])
        exp = np.array([2], dtype=np.intp)
        tm.assert_numpy_array_equal(res, exp)

    def test_searchsorted_numeric_dtypes_vector(self):
        ser = pd.Series([1, 2, 90, 1000, 3e9])
        res = ser.searchsorted([91, 2e6])
        exp = np.array([3, 4], dtype=np.intp)
        tm.assert_numpy_array_equal(res, exp)

    def test_searchsorted_datetime64_scalar(self):
        ser = pd.Series(pd.date_range("20120101", periods=10, freq="2D"))
        val = pd.Timestamp("20120102")
        res = ser.searchsorted(val)
        assert is_scalar(res)
        assert res == 1

    def test_searchsorted_datetime64_scalar_mixed_timezones(self):
        # GH 30086
        ser = pd.Series(pd.date_range("20120101", periods=10, freq="2D", tz="UTC"))
        val = pd.Timestamp("20120102", tz="America/New_York")
        res = ser.searchsorted(val)
        assert is_scalar(res)
        assert res == 1

    def test_searchsorted_datetime64_list(self):
        ser = pd.Series(pd.date_range("20120101", periods=10, freq="2D"))
        vals = [pd.Timestamp("20120102"), pd.Timestamp("20120104")]
        res = ser.searchsorted(vals)
        exp = np.array([1, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(res, exp)

    def test_searchsorted_sorter(self):
        # GH8490
        ser = pd.Series([3, 1, 2])
        res = ser.searchsorted([0, 3], sorter=np.argsort(ser))
        exp = np.array([0, 2], dtype=np.intp)
        tm.assert_numpy_array_equal(res, exp)

    def test_searchsorted_dataframe_fail(self):
        # GH#49620
        ser = pd.Series([1, 2, 3, 4, 5])
        vals = pd.DataFrame([[1, 2], [3, 4]])
        msg = "Value must be 1-D array-like or scalar, DataFrame is not supported"
        with pytest.raises(ValueError, match=msg):
            ser.searchsorted(vals)

    @pytest.mark.parametrize(
        "dtype",
        [
            "int64",
            "uint8",
            "float64",
            "complex128",
            "Int64",
            "Float64",
            pytest.param("int64[pyarrow]", marks=td.skip_if_no("pyarrow")),
        ],
    )
    @pytest.mark.parametrize("na_value", [pd.NA, None])
    @pytest.mark.parametrize("side", ["left", "right"])
    @pytest.mark.parametrize(
        "box",
        [
            list,
            tuple,
            lambda x: np.array(x, dtype=object),
            lambda x: pd.array(x, dtype="Int64"),
        ],
        ids=["list", "tuple", "object-ndarray", "Int64-array"],
    )
    def test_searchsorted_na_value_numeric(self, dtype, na_value, side, box):
        # GH#30944 missing values in the search value sort last, like np.nan
        ser = pd.Series([1, 2, 3], dtype=dtype)
        assert ser.searchsorted(na_value, side=side) == 3

        expected = np.array([1 if side == "left" else 2, 3], dtype=np.intp)
        result = ser.searchsorted(box([2, na_value]), side=side)
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize("dtype", ["int64", "Int64", "string"])
    @pytest.mark.parametrize("as_array", [False, True])
    def test_searchsorted_na_value_routes_and_sorter(self, dtype, as_array):
        # GH#30944
        values = [3, 1, 2] if dtype != "string" else ["c", "a", "b"]
        haystack = pd.Index(values, dtype=dtype)
        if as_array:
            haystack = haystack.array
        sorter = np.array([1, 2, 0])
        assert haystack.searchsorted(pd.NA, sorter=sorter) == 3
        result = haystack.searchsorted([values[2], pd.NA], sorter=sorter)
        tm.assert_numpy_array_equal(result, np.array([1, 3], dtype=np.intp))

    @pytest.mark.parametrize("side", ["left", "right"])
    @pytest.mark.parametrize("na_value", [np.nan, pd.NA])
    def test_searchsorted_na_value_float_with_nan(self, side, na_value):
        # GH#30944 pd.NA is placed where np.nan is
        ser = pd.Series([1.0, 2.0, np.nan])
        expected = np.array([1, 2] if side == "left" else [2, 3], dtype=np.intp)
        result = ser.searchsorted([2.0, na_value], side=side)
        tm.assert_numpy_array_equal(result, expected)
        assert ser.searchsorted(na_value, side=side) == expected[1]

    @pytest.mark.parametrize("side", ["left", "right"])
    def test_searchsorted_nan_value_pyarrow_float(self, side):
        # GH#30944 NaN that is not null is placed by numpy, not last
        pa = pytest.importorskip("pyarrow")
        ser = pd.Series(pd.arrays.ArrowExtensionArray(pa.array([1.0, 2.0, np.nan])))
        expected = np.array([1, 2] if side == "left" else [2, 3], dtype=np.intp)
        result = ser.searchsorted([2.0, np.nan], side=side)
        tm.assert_numpy_array_equal(result, expected)
        assert ser.searchsorted(np.nan, side=side) == expected[1]

    def test_searchsorted_na_value_distinguish_nan_and_na(self):
        # GH#30944 NaN is a value, so it is placed by numpy; NA is placed last
        with pd.option_context("future.distinguish_nan_and_na", True):
            arr = pd.array([1.0, 2.0, np.nan], dtype="Float64")
            assert arr.searchsorted(np.nan) == 2
            result = arr.searchsorted([np.nan, pd.NA])
            tm.assert_numpy_array_equal(result, np.array([2, 3], dtype=np.intp))

    def test_searchsorted_nat_numeric_raises(self):
        # GH#30944 NaT is not a missing value for a numeric dtype
        ser = pd.Series([1, 2, 3])
        msg = "'<' not supported between instances of 'int' and 'NaTType'"
        with pytest.raises(TypeError, match=msg):
            ser.searchsorted(pd.NaT)
        with pytest.raises(TypeError, match=msg):
            ser.searchsorted([2, pd.NaT])

    def test_searchsorted_na_value_bool(self):
        # GH#30944
        ser = pd.Series([False, True])
        assert ser.searchsorted(pd.NA) == 2
        result = ser.searchsorted([True, None])
        tm.assert_numpy_array_equal(result, np.array([1, 2], dtype=np.intp))

    @td.skip_if_no("pyarrow")
    @pytest.mark.parametrize(
        "dtype",
        [
            "timestamp[s][pyarrow]",
            "timestamp[ns][pyarrow]",
            "timestamp[s, tz=UTC][pyarrow]",
            "duration[ns][pyarrow]",
        ],
    )
    def test_searchsorted_na_value_pyarrow_temporal(self, dtype):
        # GH#30944 NaT was placed first
        if dtype.startswith("duration"):
            values = pd.to_timedelta(["1D", "2D"])
        else:
            values = pd.to_datetime(["2020", "2021"], utc="UTC" in dtype)
        ser = pd.Series(values).astype(dtype)
        assert ser.searchsorted(pd.NaT) == 2
        result = ser.searchsorted([values[1], pd.NaT, None])
        tm.assert_numpy_array_equal(result, np.array([1, 2, 2], dtype=np.intp))

    def test_searchsorted_na_value_string(self, string_dtype_no_object):
        # GH#30944
        ser = pd.Series(["a", "b", "c"], dtype=string_dtype_no_object)
        na_value = string_dtype_no_object.na_value
        assert ser.searchsorted(na_value) == 3

        expected = np.array([1, 3], dtype=np.intp)
        result = ser.searchsorted(pd.array(["b", na_value], dtype=ser.dtype))
        tm.assert_numpy_array_equal(result, expected)
        result = ser.searchsorted(["b", na_value])
        tm.assert_numpy_array_equal(result, expected)
        result = ser.searchsorted(["b", np.nan])
        tm.assert_numpy_array_equal(result, expected)
        result = ser.searchsorted([np.nan])
        tm.assert_numpy_array_equal(result, np.array([3], dtype=np.intp))

    def test_searchsorted_nullable_value_beyond_float_precision(self):
        # GH#63938 the non-missing values are not cast to float
        base = 2**53
        ser = pd.Series([base, base + 1, base + 2, base + 3])
        value = pd.array([base + 1, pd.NA], dtype="Int64")
        result = ser.searchsorted(value)
        tm.assert_numpy_array_equal(result, np.array([1, 4], dtype=np.intp))
