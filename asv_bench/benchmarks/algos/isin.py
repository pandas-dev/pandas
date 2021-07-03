import numpy as np

try:
    from pandas.compat import np_version_under1p20
except ImportError:
    from pandas.compat.numpy import _np_version_under1p20 as np_version_under1p20

from pandas import (
    Categorical,
    NaT,
    Series,
    date_range,
)

from ..pandas_vb_common import tm


class IsIn:

    params = [
        "int64",
        "uint64",
        "object",
        "Int64",
        "boolean",
        "bool",
        "datetime64[ns]",
        "category[object]",
        "category[int]",
        "str",
        "string[python]",
        "string[pyarrow]",
    ]
    param_names = ["dtype"]

    def setup(self, dtype):
        N = 10000

        self.mismatched = [NaT.to_datetime64()] * 2

        if dtype in ["boolean", "bool"]:
            self.series = Series(np.random.randint(0, 2, N)).astype(dtype)
            self.values = [True, False]

        elif dtype == "datetime64[ns]":
            # Note: values here is much larger than non-dt64ns cases

            # dti has length=115777
            dti = date_range(start="2015-10-26", end="2016-01-01", freq="50s")
            self.series = Series(dti)
            self.values = self.series._values[::3]
            self.mismatched = [1, 2]

        elif dtype in ["category[object]", "category[int]"]:
            # Note: sizes are different in this case than others
            n = 5 * 10 ** 5
            sample_size = 100

            arr = list(np.random.randint(0, n // 10, size=n))
            if dtype == "category[object]":
                arr = [f"s{i:04d}" for i in arr]

            self.values = np.random.choice(arr, sample_size)
            self.series = Series(arr).astype("category")

        elif dtype in ["str", "string[python]", "string[pyarrow]"]:
            try:
                self.series = Series(tm.makeStringIndex(N), dtype=dtype)
            except ImportError:
                raise NotImplementedError
            self.values = list(self.series[:2])

        else:
            self.series = Series(np.random.randint(1, 10, N)).astype(dtype)
            self.values = [1, 2]

        self.cat_values = Categorical(self.values)

    def time_isin(self, dtype):
        self.series.isin(self.values)

    def time_isin_categorical(self, dtype):
        self.series.isin(self.cat_values)

    def time_isin_empty(self, dtype):
        self.series.isin([])

    def time_isin_mismatched_dtype(self, dtype):
        self.series.isin(self.mismatched)


class IsinAlmostFullWithRandomInt:
    params = [
        [np.float64, np.int64, np.uint64, np.object_],
        range(10, 21),
        ["inside", "outside"],
    ]
    param_names = ["dtype", "exponent", "title"]

    def setup(self, dtype, exponent, title):
        M = 3 * 2 ** (exponent - 2)
        # 0.77-the maximal share of occupied buckets
        self.series = Series(np.random.randint(0, M, M)).astype(dtype)

        values = np.random.randint(0, M, M).astype(dtype)
        if title == "inside":
            self.values = values
        elif title == "outside":
            self.values = values + M
        else:
            raise ValueError(title)

    def time_isin(self, dtype, exponent, title):
        self.series.isin(self.values)


class IsinWithRandomFloat:
    params = [
        [np.float64, np.object_],
        [
            1_300,
            2_000,
            7_000,
            8_000,
            70_000,
            80_000,
            750_000,
            900_000,
        ],
        ["inside", "outside"],
    ]
    param_names = ["dtype", "size", "title"]

    def setup(self, dtype, size, title):
        self.values = np.random.rand(size)
        self.series = Series(self.values).astype(dtype)
        np.random.shuffle(self.values)

        if title == "outside":
            self.values = self.values + 0.1

    def time_isin(self, dtype, size, title):
        self.series.isin(self.values)


class IsinWithArangeSorted:
    params = [
        [np.float64, np.int64, np.uint64, np.object_],
        [
            1_000,
            2_000,
            8_000,
            100_000,
            1_000_000,
        ],
    ]
    param_names = ["dtype", "size"]

    def setup(self, dtype, size):
        self.series = Series(np.arange(size)).astype(dtype)
        self.values = np.arange(size).astype(dtype)

    def time_isin(self, dtype, size):
        self.series.isin(self.values)


class IsinWithArange:
    params = [
        [np.float64, np.int64, np.uint64, np.object_],
        [
            1_000,
            2_000,
            8_000,
        ],
        [-2, 0, 2],
    ]
    param_names = ["dtype", "M", "offset_factor"]

    def setup(self, dtype, M, offset_factor):
        offset = int(M * offset_factor)
        tmp = Series(np.random.randint(offset, M + offset, 10 ** 6))
        self.series = tmp.astype(dtype)
        self.values = np.arange(M).astype(dtype)

    def time_isin(self, dtype, M, offset_factor):
        self.series.isin(self.values)


class IsInFloat64:

    params = [
        [np.float64, "Float64"],
        ["many_different_values", "few_different_values", "only_nans_values"],
    ]
    param_names = ["dtype", "title"]

    def setup(self, dtype, title):
        N_many = 10 ** 5
        N_few = 10 ** 6
        self.series = Series([1, 2], dtype=dtype)

        if title == "many_different_values":
            # runtime is dominated by creation of the lookup-table
            self.values = np.arange(N_many, dtype=np.float64)
        elif title == "few_different_values":
            # runtime is dominated by creation of the lookup-table
            self.values = np.zeros(N_few, dtype=np.float64)
        elif title == "only_nans_values":
            # runtime is dominated by creation of the lookup-table
            self.values = np.full(N_few, np.nan, dtype=np.float64)
        else:
            raise ValueError(title)

    def time_isin(self, dtype, title):
        self.series.isin(self.values)


class IsInForObjects:
    """
    A subset of the cartesian product of cases have special motivations:

    "nans" x "nans"
        if nan-objects are different objects,
        this has the potential to trigger O(n^2) running time

    "short" x "long"
        running time dominated by the preprocessing

    "long" x "short"
        running time dominated by look-up

    "long" x "long"
        no dominating part

    "long_floats" x "long_floats"
        because of nans floats are special
        no dominating part

    """

    variants = ["nans", "short", "long", "long_floats"]

    params = [variants, variants]
    param_names = ["series_type", "vals_type"]

    def setup(self, series_type, vals_type):
        N_many = 10 ** 5

        if series_type == "nans":
            ser_vals = np.full(10 ** 4, np.nan)
        elif series_type == "short":
            ser_vals = np.arange(2)
        elif series_type == "long":
            ser_vals = np.arange(N_many)
        elif series_type == "long_floats":
            ser_vals = np.arange(N_many, dtype=np.float_)

        self.series = Series(ser_vals).astype(object)

        if vals_type == "nans":
            values = np.full(10 ** 4, np.nan)
        elif vals_type == "short":
            values = np.arange(2)
        elif vals_type == "long":
            values = np.arange(N_many)
        elif vals_type == "long_floats":
            values = np.arange(N_many, dtype=np.float_)

        self.values = values.astype(object)

    def time_isin(self, series_type, vals_type):
        self.series.isin(self.values)


class IsInLongSeriesLookUpDominates:
    params = [
        ["int64", "int32", "float64", "float32", "object", "Int64", "Float64"],
        [5, 1000],
        ["random_hits", "random_misses", "monotone_hits", "monotone_misses"],
    ]
    param_names = ["dtype", "MaxNumber", "series_type"]

    def setup(self, dtype, MaxNumber, series_type):
        N = 10 ** 7

        # https://github.com/pandas-dev/pandas/issues/39844
        if not np_version_under1p20 and dtype in ("Int64", "Float64"):
            raise NotImplementedError

        if series_type == "random_hits":
            array = np.random.randint(0, MaxNumber, N)
        if series_type == "random_misses":
            array = np.random.randint(0, MaxNumber, N) + MaxNumber
        if series_type == "monotone_hits":
            array = np.repeat(np.arange(MaxNumber), N // MaxNumber)
        if series_type == "monotone_misses":
            array = np.arange(N) + MaxNumber

        self.series = Series(array).astype(dtype)
        self.values = np.arange(MaxNumber).astype(dtype)

    def time_isin(self, dtypes, MaxNumber, series_type):
        self.series.isin(self.values)


class IsInLongSeriesValuesDominate:
    params = [
        ["int64", "int32", "float64", "float32", "object", "Int64", "Float64"],
        ["random", "monotone"],
    ]
    param_names = ["dtype", "series_type"]

    def setup(self, dtype, series_type):
        N = 10 ** 7

        # https://github.com/pandas-dev/pandas/issues/39844
        if not np_version_under1p20 and dtype in ("Int64", "Float64"):
            raise NotImplementedError

        if series_type == "random":
            vals = np.random.randint(0, 10 * N, N)
        if series_type == "monotone":
            vals = np.arange(N)

        self.values = vals.astype(dtype)
        M = 10 ** 6 + 1
        self.series = Series(np.arange(M)).astype(dtype)

    def time_isin(self, dtypes, series_type):
        self.series.isin(self.values)


class IsInWithLongTupples:
    def setup(self):
        t = tuple(range(1000))
        self.series = Series([t] * 1000)
        self.values = [t]

    def time_isin(self):
        self.series.isin(self.values)
