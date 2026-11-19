import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestToNumpy:
    def test_to_numpy(self):
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4.5]})
        expected = np.array([[1, 3], [2, 4.5]])
        result = df.to_numpy()
        tm.assert_numpy_array_equal(result, expected)

    def test_to_numpy_dtype(self):
        df = pd.DataFrame({"A": [1, 2], "B": [3, 4.5]})
        expected = np.array([[1, 3], [2, 4]], dtype="int64")
        result = df.to_numpy(dtype="int64")
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize(
        "data, dtype",
        [
            ({"a": [1, 2, 3]}, "Int64"),
            ({"a": [1, 2, 3], "b": [4, 5, 6]}, "Int64"),
            ({"a": [1, 2, 3]}, "int64"),
        ],
    )
    def test_to_numpy_na_value_no_missing(self, data, dtype):
        # GH#56233 assigning na_value is a no-op when there is no NA, even
        # when the na_value cannot be cast to the target dtype
        df = pd.DataFrame(data, dtype=dtype)
        result = df.to_numpy(na_value=float("nan"))
        expected = df.to_numpy()
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype, expected_dtype",
        [
            ("Int64", "int64"),
            ("Float64", "float64"),
            ("boolean", "bool"),
            ("category", "int64"),
            ("int64[pyarrow]", "int64"),
        ],
    )
    def test_to_numpy_multiple_extension_columns(self, dtype, expected_dtype):
        # GH#22791 match the dtype a single column converts to
        if "pyarrow" in dtype:
            pytest.importorskip("pyarrow")
        df = pd.DataFrame({"a": [1, 0], "b": [0, 1]}, dtype=dtype)
        expected = np.array([[1, 0], [0, 1]], dtype=expected_dtype)
        tm.assert_numpy_array_equal(df.to_numpy(), expected)
        tm.assert_numpy_array_equal(df.values, expected)
        tm.assert_numpy_array_equal(np.asarray(df), expected)
        tm.assert_numpy_array_equal(df.to_numpy(na_value=np.nan), expected)

    def test_to_numpy_extension_columns_with_na(self):
        # GH#22791
        df = pd.DataFrame(
            {
                "a": pd.array([1, None], dtype="Int64"),
                "b": pd.array([3, 4], dtype="Int32"),
                "c": [5.5, 6.5],
            }
        )
        expected = np.array([[1, 3, 5.5], [np.nan, 4, 6.5]])
        tm.assert_numpy_array_equal(df.to_numpy(), expected)

        expected = np.array([[1, 3, 5.5], [-1, 4, 6.5]])
        tm.assert_numpy_array_equal(df.to_numpy(na_value=-1), expected)

        # None does not fit in float64, so fall back to object
        expected = np.array([[1, 3, 5.5], [None, 4, 6.5]], dtype=object)
        tm.assert_numpy_array_equal(df.to_numpy(na_value=None), expected)

        # converted per column, not via float64
        expected = np.array([[1, 3, 5.5], [pd.NA, 4, 6.5]], dtype=object)
        tm.assert_numpy_array_equal(np.asarray(df, dtype=object), expected)

        expected = np.array([["1", "3", "5.5"], ["<NA>", "4", "6.5"]], dtype="U4")
        tm.assert_numpy_array_equal(np.asarray(df, dtype="U4"), expected)

        # like a Series, raise instead of filling in the NA
        for dtype in ["int64", "bool"]:
            msg = f"cannot convert to '{dtype}'-dtype NumPy array with missing values"
            with pytest.raises(ValueError, match=msg):
                np.asarray(df, dtype=dtype)

    def test_to_numpy_extension_columns_datetime_na_value(self):
        # GH#22791 na_value is only filled into a float result
        cat = pd.Categorical(pd.to_datetime(["2020-01-01", None]))
        df = pd.DataFrame({"a": cat, "b": cat})
        ts = pd.Timestamp("2020-01-01")
        expected = np.array([[ts, ts], [pd.NaT, pd.NaT]], dtype=object)
        tm.assert_numpy_array_equal(df.to_numpy(na_value=pd.NaT), expected)

    def test_to_numpy_extension_columns_object_result(self):
        # GH#22791 object results keep the boxing of the object conversion
        pa = pytest.importorskip("pyarrow")
        ts = pd.array(
            [pd.Timestamp("2020-01-01"), None], dtype="timestamp[ns][pyarrow]"
        )
        dates = pd.array(
            pa.array([1, None], type=pa.date32()), dtype=pd.ArrowDtype(pa.date32())
        )
        df = pd.DataFrame({"a": ts, "b": dates})
        expected = np.array(
            [[pd.Timestamp("2020-01-01"), dates[0]], [pd.NA, pd.NA]], dtype=object
        )
        tm.assert_numpy_array_equal(df.to_numpy(), expected)

    def test_to_numpy_copy(self):
        arr = np.random.default_rng(2).standard_normal((4, 3))
        df = pd.DataFrame(arr)
        assert df.values.base is not arr
        assert df.to_numpy(copy=False).base is df.values.base
        assert df.to_numpy(copy=True).base is not arr

        # we still don't want a copy when na_value=np.nan is passed,
        #  and that can be respected because we are already numpy-float
        assert df.to_numpy(copy=False).base is df.values.base

    def test_to_numpy_mixed_dtype_to_str(self):
        # https://github.com/pandas-dev/pandas/issues/35455
        df = pd.DataFrame([[pd.Timestamp("2020-01-01 00:00:00"), 100.0]])
        result = df.to_numpy(dtype=str)
        expected = np.array([["2020-01-01 00:00:00", "100.0"]], dtype=str)
        tm.assert_numpy_array_equal(result, expected)

    def test_to_numpy_datetime_with_na(self):
        # GH #53115
        dti = pd.date_range("2016-01-01", periods=3, unit="ns")
        df = pd.DataFrame(dti)
        df.iloc[0, 0] = pd.NaT
        expected = np.array([[np.nan], [1.45169280e18], [1.45177920e18]])
        result = df.to_numpy(float, na_value=np.nan)
        tm.assert_numpy_array_equal(result, expected)

        df = pd.DataFrame(
            {
                "a": [
                    pd.Timestamp("1970-01-01").as_unit("s"),
                    pd.Timestamp("1970-01-02").as_unit("s"),
                    pd.NaT,
                ],
                "b": [
                    pd.Timestamp("1970-01-01").as_unit("s"),
                    np.nan,
                    pd.Timestamp("1970-01-02").as_unit("s"),
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
