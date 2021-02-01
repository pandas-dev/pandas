from datetime import datetime

import numpy as np

from pandas import Categorical, NaT, Series, date_range

from .pandas_vb_common import tm


class SeriesConstructor:

    params = [None, "dict"]
    param_names = ["data"]

    def setup(self, data):
        self.idx = date_range(
            start=datetime(2015, 10, 26), end=datetime(2016, 1, 1), freq="50s"
        )
        dict_data = dict(zip(self.idx, range(len(self.idx))))
        self.data = None if data is None else dict_data

    def time_constructor(self, data):
        Series(data=self.data, index=self.idx)


class IsIn:

    params = ["int64", "uint64", "object", "Int64"]
    param_names = ["dtype"]

    def setup(self, dtype):
        N = 10000
        self.s = Series(np.random.randint(1, 10, N)).astype(dtype)
        self.values = [1, 2]

    def time_isin(self, dtypes):
        self.s.isin(self.values)


class IsInBoolean:

    params = ["boolean", "bool"]
    param_names = ["dtype"]

    def setup(self, dtype):
        N = 10000
        self.s = Series(np.random.randint(0, 2, N)).astype(dtype)
        self.values = [True, False]

    def time_isin(self, dtypes):
        self.s.isin(self.values)


class IsInDatetime64:
    def setup(self):
        dti = date_range(
            start=datetime(2015, 10, 26), end=datetime(2016, 1, 1), freq="50s"
        )
        self.ser = Series(dti)
        self.subset = self.ser._values[::3]
        self.cat_subset = Categorical(self.subset)

    def time_isin(self):
        self.ser.isin(self.subset)

    def time_isin_cat_values(self):
        self.ser.isin(self.cat_subset)

    def time_isin_mismatched_dtype(self):
        self.ser.isin([1, 2])

    def time_isin_empty(self):
        self.ser.isin([])


class IsInFloat64:

    params = [np.float64, "Float64"]
    param_names = ["dtype"]

    def setup(self, dtype):
        N_many = 10 ** 5
        N_few = 10 ** 6
        self.small = Series([1, 2], dtype=dtype)
        self.many_different_values = np.arange(N_many, dtype=np.float64)
        self.few_different_values = np.zeros(N_few, dtype=np.float64)
        self.only_nans_values = np.full(N_few, np.nan, dtype=np.float64)

    def time_isin_many_different(self, dtypes):
        # runtime is dominated by creation of the lookup-table
        self.small.isin(self.many_different_values)

    def time_isin_few_different(self, dtypes):
        # runtime is dominated by creation of the lookup-table
        self.small.isin(self.few_different_values)

    def time_isin_nan_values(self, dtypes):
        # runtime is dominated by creation of the lookup-table
        self.small.isin(self.few_different_values)


class IsInForObjects:
    def setup(self):
        self.s_nans = Series(np.full(10 ** 4, np.nan)).astype(object)
        self.vals_nans = np.full(10 ** 4, np.nan).astype(object)
        self.s_short = Series(np.arange(2)).astype(object)
        self.s_long = Series(np.arange(10 ** 5)).astype(object)
        self.vals_short = np.arange(2).astype(object)
        self.vals_long = np.arange(10 ** 5).astype(object)
        # because of nans floats are special:
        self.s_long_floats = Series(np.arange(10 ** 5, dtype=np.float)).astype(object)
        self.vals_long_floats = np.arange(10 ** 5, dtype=np.float).astype(object)

    def time_isin_nans(self):
        # if nan-objects are different objects,
        # this has the potential to trigger O(n^2) running time
        self.s_nans.isin(self.vals_nans)

    def time_isin_short_series_long_values(self):
        # running time dominated by the preprocessing
        self.s_short.isin(self.vals_long)

    def time_isin_long_series_short_values(self):
        # running time dominated by look-up
        self.s_long.isin(self.vals_short)

    def time_isin_long_series_long_values(self):
        # no dominating part
        self.s_long.isin(self.vals_long)

    def time_isin_long_series_long_values_floats(self):
        # no dominating part
        self.s_long_floats.isin(self.vals_long_floats)


class IsInLongSeriesLookUpDominates:
    params = [
        ["int64", "int32", "float64", "float32", "object", "Int64", "Float64"],
        [5, 1000],
        ["random_hits", "random_misses", "monotone_hits", "monotone_misses"],
    ]
    param_names = ["dtype", "MaxNumber", "series_type"]

    def setup(self, dtype, MaxNumber, series_type):
        N = 10 ** 7
        if series_type == "random_hits":
            np.random.seed(42)
            array = np.random.randint(0, MaxNumber, N)
        if series_type == "random_misses":
            np.random.seed(42)
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
        if series_type == "random":
            np.random.seed(42)
            vals = np.random.randint(0, 10 * N, N)
        if series_type == "monotone":
            vals = np.arange(N)
        self.values = vals.astype(dtype)
        M = 10 ** 6 + 1
        self.series = Series(np.arange(M)).astype(dtype)

    def time_isin(self, dtypes, series_type):
        self.series.isin(self.values)


class NSort:

    params = ["first", "last", "all"]
    param_names = ["keep"]

    def setup(self, keep):
        self.s = Series(np.random.randint(1, 10, 100000))

    def time_nlargest(self, keep):
        self.s.nlargest(3, keep=keep)

    def time_nsmallest(self, keep):
        self.s.nsmallest(3, keep=keep)


class Dropna:

    params = ["int", "datetime"]
    param_names = ["dtype"]

    def setup(self, dtype):
        N = 10 ** 6
        data = {
            "int": np.random.randint(1, 10, N),
            "datetime": date_range("2000-01-01", freq="S", periods=N),
        }
        self.s = Series(data[dtype])
        if dtype == "datetime":
            self.s[np.random.randint(1, N, 100)] = NaT

    def time_dropna(self, dtype):
        self.s.dropna()


class SearchSorted:

    goal_time = 0.2
    params = [
        "int8",
        "int16",
        "int32",
        "int64",
        "uint8",
        "uint16",
        "uint32",
        "uint64",
        "float16",
        "float32",
        "float64",
        "str",
    ]
    param_names = ["dtype"]

    def setup(self, dtype):
        N = 10 ** 5
        data = np.array([1] * N + [2] * N + [3] * N).astype(dtype)
        self.s = Series(data)

    def time_searchsorted(self, dtype):
        key = "2" if dtype == "str" else 2
        self.s.searchsorted(key)


class Map:

    params = (["dict", "Series", "lambda"], ["object", "category", "int"])
    param_names = "mapper"

    def setup(self, mapper, dtype):
        map_size = 1000
        map_data = Series(map_size - np.arange(map_size), dtype=dtype)

        # construct mapper
        if mapper == "Series":
            self.map_data = map_data
        elif mapper == "dict":
            self.map_data = map_data.to_dict()
        elif mapper == "lambda":
            map_dict = map_data.to_dict()
            self.map_data = lambda x: map_dict[x]
        else:
            raise NotImplementedError

        self.s = Series(np.random.randint(0, map_size, 10000), dtype=dtype)

    def time_map(self, mapper, *args, **kwargs):
        self.s.map(self.map_data)


class Clip:
    params = [50, 1000, 10 ** 5]
    param_names = ["n"]

    def setup(self, n):
        self.s = Series(np.random.randn(n))

    def time_clip(self, n):
        self.s.clip(0, 1)


class ValueCounts:

    params = [[10 ** 3, 10 ** 4, 10 ** 5], ["int", "uint", "float", "object"]]
    param_names = ["N", "dtype"]

    def setup(self, N, dtype):
        self.s = Series(np.random.randint(0, N, size=10 * N)).astype(dtype)

    def time_value_counts(self, N, dtype):
        self.s.value_counts()


class Mode:

    params = [[10 ** 3, 10 ** 4, 10 ** 5], ["int", "uint", "float", "object"]]
    param_names = ["N", "dtype"]

    def setup(self, N, dtype):
        np.random.seed(42)
        self.s = Series(np.random.randint(0, N, size=10 * N)).astype(dtype)

    def time_mode(self, N, dtype):
        self.s.mode()


class Dir:
    def setup(self):
        self.s = Series(index=tm.makeStringIndex(10000))

    def time_dir_strings(self):
        dir(self.s)


class SeriesGetattr:
    # https://github.com/pandas-dev/pandas/issues/19764
    def setup(self):
        self.s = Series(1, index=date_range("2012-01-01", freq="s", periods=10 ** 6))

    def time_series_datetimeindex_repr(self):
        getattr(self.s, "a", None)


class All:

    params = [[10 ** 3, 10 ** 6], ["fast", "slow"], ["bool", "boolean"]]
    param_names = ["N", "case", "dtype"]

    def setup(self, N, case, dtype):
        val = case != "fast"
        self.s = Series([val] * N, dtype=dtype)

    def time_all(self, N, case, dtype):
        self.s.all()


class Any:

    params = [[10 ** 3, 10 ** 6], ["fast", "slow"], ["bool", "boolean"]]
    param_names = ["N", "case", "dtype"]

    def setup(self, N, case, dtype):
        val = case == "fast"
        self.s = Series([val] * N, dtype=dtype)

    def time_any(self, N, case, dtype):
        self.s.any()


class NanOps:

    params = [
        [
            "var",
            "mean",
            "median",
            "max",
            "min",
            "sum",
            "std",
            "sem",
            "argmax",
            "skew",
            "kurt",
            "prod",
        ],
        [10 ** 3, 10 ** 6],
        ["int8", "int32", "int64", "float64", "Int64", "boolean"],
    ]
    param_names = ["func", "N", "dtype"]

    def setup(self, func, N, dtype):
        if func == "argmax" and dtype in {"Int64", "boolean"}:
            # Skip argmax for nullable int since this doesn't work yet (GH-24382)
            raise NotImplementedError
        self.s = Series([1] * N, dtype=dtype)
        self.func = getattr(self.s, func)

    def time_func(self, func, N, dtype):
        self.func()


class Rank:

    param_names = ["dtype"]
    params = [
        ["int", "uint", "float", "object"],
    ]

    def setup(self, dtype):
        self.s = Series(np.random.randint(0, 1000, size=100000), dtype=dtype)

    def time_rank(self, dtype):
        self.s.rank()


from .pandas_vb_common import setup  # noqa: F401 isort:skip
