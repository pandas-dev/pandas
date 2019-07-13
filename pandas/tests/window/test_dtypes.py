from itertools import product

import numpy as np
import pytest

from pandas import DataFrame, Series
from pandas.core.base import DataError
import pandas.util.testing as tm

# gh-12373 : rolling functions error on float32 data
# make sure rolling functions works for different dtypes
#
# NOTE that these are yielded tests and so _create_data
# is explicitly called.
#
# further note that we are only checking rolling for fully dtype
# compliance (though both expanding and ewm inherit)


class Dtype:
    window = 2

    funcs = {
        "count": lambda v: v.count(),
        "max": lambda v: v.max(),
        "min": lambda v: v.min(),
        "sum": lambda v: v.sum(),
        "mean": lambda v: v.mean(),
        "std": lambda v: v.std(),
        "var": lambda v: v.var(),
        "median": lambda v: v.median(),
    }

    def get_expects(self):
        expects = {
            "sr1": {
                "count": Series([1, 2, 2, 2, 2], dtype="float64"),
                "max": Series([np.nan, 1, 2, 3, 4], dtype="float64"),
                "min": Series([np.nan, 0, 1, 2, 3], dtype="float64"),
                "sum": Series([np.nan, 1, 3, 5, 7], dtype="float64"),
                "mean": Series([np.nan, 0.5, 1.5, 2.5, 3.5], dtype="float64"),
                "std": Series([np.nan] + [np.sqrt(0.5)] * 4, dtype="float64"),
                "var": Series([np.nan, 0.5, 0.5, 0.5, 0.5], dtype="float64"),
                "median": Series([np.nan, 0.5, 1.5, 2.5, 3.5], dtype="float64"),
            },
            "sr2": {
                "count": Series([1, 2, 2, 2, 2], dtype="float64"),
                "max": Series([np.nan, 10, 8, 6, 4], dtype="float64"),
                "min": Series([np.nan, 8, 6, 4, 2], dtype="float64"),
                "sum": Series([np.nan, 18, 14, 10, 6], dtype="float64"),
                "mean": Series([np.nan, 9, 7, 5, 3], dtype="float64"),
                "std": Series([np.nan] + [np.sqrt(2)] * 4, dtype="float64"),
                "var": Series([np.nan, 2, 2, 2, 2], dtype="float64"),
                "median": Series([np.nan, 9, 7, 5, 3], dtype="float64"),
            },
            "sr3": {
                "count": Series([1, 2, 2, 1, 1], dtype="float64"),
                "max": Series([np.nan, 1, 2, np.nan, np.nan], dtype="float64"),
                "min": Series([np.nan, 0, 1, np.nan, np.nan], dtype="float64"),
                "sum": Series([np.nan, 1, 3, np.nan, np.nan], dtype="float64"),
                "mean": Series([np.nan, 0.5, 1.5, np.nan, np.nan], dtype="float64"),
                "std": Series(
                    [np.nan] + [np.sqrt(0.5)] * 2 + [np.nan] * 2, dtype="float64"
                ),
                "var": Series([np.nan, 0.5, 0.5, np.nan, np.nan], dtype="float64"),
                "median": Series([np.nan, 0.5, 1.5, np.nan, np.nan], dtype="float64"),
            },
            "df": {
                "count": DataFrame(
                    {0: Series([1, 2, 2, 2, 2]), 1: Series([1, 2, 2, 2, 2])},
                    dtype="float64",
                ),
                "max": DataFrame(
                    {0: Series([np.nan, 2, 4, 6, 8]), 1: Series([np.nan, 3, 5, 7, 9])},
                    dtype="float64",
                ),
                "min": DataFrame(
                    {0: Series([np.nan, 0, 2, 4, 6]), 1: Series([np.nan, 1, 3, 5, 7])},
                    dtype="float64",
                ),
                "sum": DataFrame(
                    {
                        0: Series([np.nan, 2, 6, 10, 14]),
                        1: Series([np.nan, 4, 8, 12, 16]),
                    },
                    dtype="float64",
                ),
                "mean": DataFrame(
                    {0: Series([np.nan, 1, 3, 5, 7]), 1: Series([np.nan, 2, 4, 6, 8])},
                    dtype="float64",
                ),
                "std": DataFrame(
                    {
                        0: Series([np.nan] + [np.sqrt(2)] * 4),
                        1: Series([np.nan] + [np.sqrt(2)] * 4),
                    },
                    dtype="float64",
                ),
                "var": DataFrame(
                    {0: Series([np.nan, 2, 2, 2, 2]), 1: Series([np.nan, 2, 2, 2, 2])},
                    dtype="float64",
                ),
                "median": DataFrame(
                    {0: Series([np.nan, 1, 3, 5, 7]), 1: Series([np.nan, 2, 4, 6, 8])},
                    dtype="float64",
                ),
            },
        }
        return expects

    def _create_dtype_data(self, dtype):
        sr1 = Series(np.arange(5), dtype=dtype)
        sr2 = Series(np.arange(10, 0, -2), dtype=dtype)
        sr3 = sr1.copy()
        sr3[3] = np.NaN
        df = DataFrame(np.arange(10).reshape((5, 2)), dtype=dtype)

        data = {"sr1": sr1, "sr2": sr2, "sr3": sr3, "df": df}

        return data

    def _create_data(self):
        self.data = self._create_dtype_data(self.dtype)
        self.expects = self.get_expects()

    def test_dtypes(self):
        self._create_data()
        for f_name, d_name in product(self.funcs.keys(), self.data.keys()):

            f = self.funcs[f_name]
            d = self.data[d_name]
            exp = self.expects[d_name][f_name]
            self.check_dtypes(f, f_name, d, d_name, exp)

    def check_dtypes(self, f, f_name, d, d_name, exp):
        roll = d.rolling(window=self.window)
        result = f(roll)
        tm.assert_almost_equal(result, exp)


class TestDtype_object(Dtype):
    dtype = object


class Dtype_integer(Dtype):
    pass


class TestDtype_int8(Dtype_integer):
    dtype = np.int8


class TestDtype_int16(Dtype_integer):
    dtype = np.int16


class TestDtype_int32(Dtype_integer):
    dtype = np.int32


class TestDtype_int64(Dtype_integer):
    dtype = np.int64


class Dtype_uinteger(Dtype):
    pass


class TestDtype_uint8(Dtype_uinteger):
    dtype = np.uint8


class TestDtype_uint16(Dtype_uinteger):
    dtype = np.uint16


class TestDtype_uint32(Dtype_uinteger):
    dtype = np.uint32


class TestDtype_uint64(Dtype_uinteger):
    dtype = np.uint64


class Dtype_float(Dtype):
    pass


class TestDtype_float16(Dtype_float):
    dtype = np.float16


class TestDtype_float32(Dtype_float):
    dtype = np.float32


class TestDtype_float64(Dtype_float):
    dtype = np.float64


class TestDtype_category(Dtype):
    dtype = "category"
    include_df = False

    def _create_dtype_data(self, dtype):
        sr1 = Series(range(5), dtype=dtype)
        sr2 = Series(range(10, 0, -2), dtype=dtype)

        data = {"sr1": sr1, "sr2": sr2}

        return data


class DatetimeLike(Dtype):
    def check_dtypes(self, f, f_name, d, d_name, exp):

        roll = d.rolling(window=self.window)
        if f_name == "count":
            result = f(roll)
            tm.assert_almost_equal(result, exp)

        else:
            with pytest.raises(DataError):
                f(roll)


class TestDtype_timedelta(DatetimeLike):
    dtype = np.dtype("m8[ns]")


class TestDtype_datetime(DatetimeLike):
    dtype = np.dtype("M8[ns]")


class TestDtype_datetime64UTC(DatetimeLike):
    dtype = "datetime64[ns, UTC]"

    def _create_data(self):
        pytest.skip(
            "direct creation of extension dtype "
            "datetime64[ns, UTC] is not supported ATM"
        )
