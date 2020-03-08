import operator
import warnings

import numpy as np

import pandas as pd
from pandas import DataFrame, Series, Timestamp, date_range, to_timedelta
import pandas._testing as tm
from pandas.core.algorithms import checked_add_with_arr

from .pandas_vb_common import numeric_dtypes

try:
    import pandas.core.computation.expressions as expr
except ImportError:
    import pandas.computation.expressions as expr
try:
    import pandas.tseries.holiday
except ImportError:
    pass


class IntFrameWithScalar:
    params = [
        [np.float64, np.int64],
        [2, 3.0, np.int32(4), np.float64(5)],
        [
            operator.add,
            operator.sub,
            operator.mul,
            operator.truediv,
            operator.floordiv,
            operator.pow,
            operator.mod,
            operator.eq,
            operator.ne,
            operator.gt,
            operator.ge,
            operator.lt,
            operator.le,
        ],
    ]
    param_names = ["dtype", "scalar", "op"]

    def setup(self, dtype, scalar, op):
        arr = np.random.randn(20000, 100)
        self.df = DataFrame(arr.astype(dtype))

    def time_frame_op_with_scalar(self, dtype, scalar, op):
        op(self.df, scalar)


class MixedFrameWithSeriesAxis0:
    params = [
        [
            "eq",
            "ne",
            "lt",
            "le",
            "ge",
            "gt",
            "add",
            "sub",
            "div",
            "floordiv",
            "mul",
            "pow",
        ]
    ]
    param_names = ["opname"]

    def setup(self, opname):
        arr = np.arange(10 ** 6).reshape(100, -1)
        df = DataFrame(arr)
        df["C"] = 1.0
        self.df = df
        self.ser = df[0]

    def time_frame_op_with_series_axis0(self, opname):
        getattr(self.df, opname)(self.ser, axis=0)


class Ops:

    params = [[True, False], ["default", 1]]
    param_names = ["use_numexpr", "threads"]

    def setup(self, use_numexpr, threads):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))

        if threads != "default":
            expr.set_numexpr_threads(threads)
        if not use_numexpr:
            expr.set_use_numexpr(False)

    def time_frame_add(self, use_numexpr, threads):
        self.df + self.df2

    def time_frame_mult(self, use_numexpr, threads):
        self.df * self.df2

    def time_frame_multi_and(self, use_numexpr, threads):
        self.df[(self.df > 0) & (self.df2 > 0)]

    def time_frame_comparison(self, use_numexpr, threads):
        self.df > self.df2

    def teardown(self, use_numexpr, threads):
        expr.set_use_numexpr(True)
        expr.set_numexpr_threads()


class Ops2:
    def setup(self):
        N = 10 ** 3
        self.df = DataFrame(np.random.randn(N, N))
        self.df2 = DataFrame(np.random.randn(N, N))

        self.df_int = DataFrame(
            np.random.randint(
                np.iinfo(np.int16).min, np.iinfo(np.int16).max, size=(N, N)
            )
        )
        self.df2_int = DataFrame(
            np.random.randint(
                np.iinfo(np.int16).min, np.iinfo(np.int16).max, size=(N, N)
            )
        )

        self.s = Series(np.random.randn(N))

    # Division

    def time_frame_float_div(self):
        self.df // self.df2

    def time_frame_float_div_by_zero(self):
        self.df / 0

    def time_frame_float_floor_by_zero(self):
        self.df // 0

    def time_frame_int_div_by_zero(self):
        self.df_int / 0

    # Modulo

    def time_frame_int_mod(self):
        self.df_int % self.df2_int

    def time_frame_float_mod(self):
        self.df % self.df2

    # Dot product

    def time_frame_dot(self):
        self.df.dot(self.df2)

    def time_series_dot(self):
        self.s.dot(self.s)

    def time_frame_series_dot(self):
        self.df.dot(self.s)


class Timeseries:

    params = [None, "US/Eastern"]
    param_names = ["tz"]

    def setup(self, tz):
        N = 10 ** 6
        halfway = (N // 2) - 1
        self.s = Series(date_range("20010101", periods=N, freq="T", tz=tz))
        self.ts = self.s[halfway]

        self.s2 = Series(date_range("20010101", periods=N, freq="s", tz=tz))

    def time_series_timestamp_compare(self, tz):
        self.s <= self.ts

    def time_timestamp_series_compare(self, tz):
        self.ts >= self.s

    def time_timestamp_ops_diff(self, tz):
        self.s2.diff()

    def time_timestamp_ops_diff_with_shift(self, tz):
        self.s - self.s.shift()


class IrregularOps:
    def setup(self):
        N = 10 ** 5
        idx = date_range(start="1/1/2000", periods=N, freq="s")
        s = Series(np.random.randn(N), index=idx)
        self.left = s.sample(frac=1)
        self.right = s.sample(frac=1)

    def time_add(self):
        self.left + self.right


class TimedeltaOps:
    def setup(self):
        self.td = to_timedelta(np.arange(1000000))
        self.ts = Timestamp("2000")

    def time_add_td_ts(self):
        self.td + self.ts


class CategoricalComparisons:
    params = ["__lt__", "__le__", "__eq__", "__ne__", "__ge__", "__gt__"]
    param_names = ["op"]

    def setup(self, op):
        N = 10 ** 5
        self.cat = pd.Categorical(list("aabbcd") * N, ordered=True)

    def time_categorical_op(self, op):
        getattr(self.cat, op)("b")


class IndexArithmetic:

    params = ["float", "int"]
    param_names = ["dtype"]

    def setup(self, dtype):
        N = 10 ** 6
        indexes = {"int": "makeIntIndex", "float": "makeFloatIndex"}
        self.index = getattr(tm, indexes[dtype])(N)

    def time_add(self, dtype):
        self.index + 2

    def time_subtract(self, dtype):
        self.index - 2

    def time_multiply(self, dtype):
        self.index * 2

    def time_divide(self, dtype):
        self.index / 2

    def time_modulo(self, dtype):
        self.index % 2


class NumericInferOps:
    # from GH 7332
    params = numeric_dtypes
    param_names = ["dtype"]

    def setup(self, dtype):
        N = 5 * 10 ** 5
        self.df = DataFrame(
            {"A": np.arange(N).astype(dtype), "B": np.arange(N).astype(dtype)}
        )

    def time_add(self, dtype):
        self.df["A"] + self.df["B"]

    def time_subtract(self, dtype):
        self.df["A"] - self.df["B"]

    def time_multiply(self, dtype):
        self.df["A"] * self.df["B"]

    def time_divide(self, dtype):
        self.df["A"] / self.df["B"]

    def time_modulo(self, dtype):
        self.df["A"] % self.df["B"]


class DateInferOps:
    # from GH 7332
    def setup_cache(self):
        N = 5 * 10 ** 5
        df = DataFrame({"datetime64": np.arange(N).astype("datetime64[ms]")})
        df["timedelta"] = df["datetime64"] - df["datetime64"]
        return df

    def time_subtract_datetimes(self, df):
        df["datetime64"] - df["datetime64"]

    def time_timedelta_plus_datetime(self, df):
        df["timedelta"] + df["datetime64"]

    def time_add_timedeltas(self, df):
        df["timedelta"] + df["timedelta"]


class AddOverflowScalar:

    params = [1, -1, 0]
    param_names = ["scalar"]

    def setup(self, scalar):
        N = 10 ** 6
        self.arr = np.arange(N)

    def time_add_overflow_scalar(self, scalar):
        checked_add_with_arr(self.arr, scalar)


class AddOverflowArray:
    def setup(self):
        N = 10 ** 6
        self.arr = np.arange(N)
        self.arr_rev = np.arange(-N, 0)
        self.arr_mixed = np.array([1, -1]).repeat(N / 2)
        self.arr_nan_1 = np.random.choice([True, False], size=N)
        self.arr_nan_2 = np.random.choice([True, False], size=N)

    def time_add_overflow_arr_rev(self):
        checked_add_with_arr(self.arr, self.arr_rev)

    def time_add_overflow_arr_mask_nan(self):
        checked_add_with_arr(self.arr, self.arr_mixed, arr_mask=self.arr_nan_1)

    def time_add_overflow_b_mask_nan(self):
        checked_add_with_arr(self.arr, self.arr_mixed, b_mask=self.arr_nan_1)

    def time_add_overflow_both_arg_nan(self):
        checked_add_with_arr(
            self.arr, self.arr_mixed, arr_mask=self.arr_nan_1, b_mask=self.arr_nan_2
        )


hcal = pd.tseries.holiday.USFederalHolidayCalendar()
# These offsets currently raise a NotImplimentedError with .apply_index()
non_apply = [
    pd.offsets.Day(),
    pd.offsets.BYearEnd(),
    pd.offsets.BYearBegin(),
    pd.offsets.BQuarterEnd(),
    pd.offsets.BQuarterBegin(),
    pd.offsets.BMonthEnd(),
    pd.offsets.BMonthBegin(),
    pd.offsets.CustomBusinessDay(),
    pd.offsets.CustomBusinessDay(calendar=hcal),
    pd.offsets.CustomBusinessMonthBegin(calendar=hcal),
    pd.offsets.CustomBusinessMonthEnd(calendar=hcal),
    pd.offsets.CustomBusinessMonthEnd(calendar=hcal),
]
other_offsets = [
    pd.offsets.YearEnd(),
    pd.offsets.YearBegin(),
    pd.offsets.QuarterEnd(),
    pd.offsets.QuarterBegin(),
    pd.offsets.MonthEnd(),
    pd.offsets.MonthBegin(),
    pd.offsets.DateOffset(months=2, days=2),
    pd.offsets.BusinessDay(),
    pd.offsets.SemiMonthEnd(),
    pd.offsets.SemiMonthBegin(),
]
offsets = non_apply + other_offsets


class OffsetArrayArithmetic:

    params = offsets
    param_names = ["offset"]

    def setup(self, offset):
        N = 10000
        rng = pd.date_range(start="1/1/2000", periods=N, freq="T")
        self.rng = rng
        self.ser = pd.Series(rng)

    def time_add_series_offset(self, offset):
        with warnings.catch_warnings(record=True):
            self.ser + offset

    def time_add_dti_offset(self, offset):
        with warnings.catch_warnings(record=True):
            self.rng + offset


class ApplyIndex:
    params = other_offsets
    param_names = ["offset"]

    def setup(self, offset):
        N = 10000
        rng = pd.date_range(start="1/1/2000", periods=N, freq="T")
        self.rng = rng

    def time_apply_index(self, offset):
        offset.apply_index(self.rng)


from .pandas_vb_common import setup  # noqa: F401 isort:skip
