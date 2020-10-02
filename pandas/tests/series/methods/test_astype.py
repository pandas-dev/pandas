import numpy as np
import pytest

from pandas import Interval, Series, Timestamp, date_range
import pandas._testing as tm


class TestAstype:
    def test_astype_dt64_to_str(self):
        # GH#10442 : testing astype(str) is correct for Series/DatetimeIndex
        dti = date_range("2012-01-01", periods=3)
        result = Series(dti).astype(str)
        expected = Series(["2012-01-01", "2012-01-02", "2012-01-03"], dtype=object)
        tm.assert_series_equal(result, expected)

    def test_astype_dt64tz_to_str(self):
        # GH#10442 : testing astype(str) is correct for Series/DatetimeIndex
        dti_tz = date_range("2012-01-01", periods=3, tz="US/Eastern")
        result = Series(dti_tz).astype(str)
        expected = Series(
            [
                "2012-01-01 00:00:00-05:00",
                "2012-01-02 00:00:00-05:00",
                "2012-01-03 00:00:00-05:00",
            ],
            dtype=object,
        )
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "values",
        [
            Series(["x", "y", "z"], dtype="string"),
            Series(["x", "y", "z"], dtype="category"),
            Series(3 * [Timestamp("2020-01-01", tz="UTC")]),
            Series(3 * [Interval(0, 1)]),
        ],
    )
    @pytest.mark.parametrize("errors", ["raise", "ignore"])
    def test_astype_ignores_errors_for_extension_dtypes(self, values, errors):
        # https://github.com/pandas-dev/pandas/issues/35471
        if errors == "ignore":
            expected = values
            result = values.astype(float, errors="ignore")
            tm.assert_series_equal(result, expected)
        else:
            msg = "(Cannot cast)|(could not convert)"
            with pytest.raises((ValueError, TypeError), match=msg):
                values.astype(float, errors=errors)

    @pytest.mark.parametrize("dtype", [np.float16, np.float32, np.float64])
    def test_astype_from_float_to_str(self, dtype):
        # https://github.com/pandas-dev/pandas/issues/36451
        s = Series([0.1], dtype=dtype)
        result = s.astype(str)
        expected = Series(["0.1"])
        tm.assert_series_equal(result, expected)
