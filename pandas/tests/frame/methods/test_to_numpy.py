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

    def test_to_numpy_timedelta_str_and_cross_kind(self):
        # GH#44484 numpy's repr gave "86400000000 microseco"
        tdi = pd.to_timedelta(["1D", None]).as_unit("us")
        df = pd.DataFrame({"a": tdi, "b": tdi})
        result = df.to_numpy(str, na_value="not available")
        expected = np.array([["1 days"] * 2, ["not available"] * 2], dtype=str)
        tm.assert_numpy_array_equal(result, expected)

        with pytest.raises(TypeError, match="Cannot cast TimedeltaArray"):
            df.to_numpy("M8[us]")

    def test_to_numpy_timedelta_multi_block(self):
        # GH#44484 frames with more than one block go through _interleave
        tdi = pd.to_timedelta(["1D", None]).as_unit("us")
        dti = pd.to_datetime(["2012-01-01", None])
        df = pd.DataFrame({"a": tdi, "b": dti})
        with pytest.raises(TypeError, match="Cannot cast DatetimeArray"):
            df.to_numpy("m8[us]")

        df = pd.DataFrame({"a": tdi, "b": [1.5, 2.5]})
        result = df.to_numpy(np.dtypes.StringDType())
        expected = np.array(
            [["1 days", "1.5"], ["NaT", "2.5"]],
            dtype=np.dtypes.StringDType(),
        )
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
