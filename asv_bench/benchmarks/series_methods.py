from datetime import datetime

import numpy as np

from pandas import NaT, Series, date_range

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

    params = ["int64", "uint64", "object"]
    param_names = ["dtype"]

    def setup(self, dtype):
        self.s = Series(np.random.randint(1, 10, 100000)).astype(dtype)
        self.values = [1, 2]

    def time_isin(self, dtypes):
        self.s.isin(self.values)


class IsInFloat64:
    def setup(self):
        self.small = Series([1, 2], dtype=np.float64)
        self.many_different_values = np.arange(10 ** 6, dtype=np.float64)
        self.few_different_values = np.zeros(10 ** 7, dtype=np.float64)
        self.only_nans_values = np.full(10 ** 7, np.nan, dtype=np.float64)

    def time_isin_many_different(self):
        # runtime is dominated by creation of the lookup-table
        self.small.isin(self.many_different_values)

    def time_isin_few_different(self):
        # runtime is dominated by creation of the lookup-table
        self.small.isin(self.few_different_values)

    def time_isin_nan_values(self):
        # runtime is dominated by creation of the lookup-table
        self.small.isin(self.few_different_values)


class IsInForObjects:
    def setup(self):
        self.s_nans = Series(np.full(10 ** 4, np.nan)).astype(np.object)
        self.vals_nans = np.full(10 ** 4, np.nan).astype(np.object)
        self.s_short = Series(np.arange(2)).astype(np.object)
        self.s_long = Series(np.arange(10 ** 5)).astype(np.object)
        self.vals_short = np.arange(2).astype(np.object)
        self.vals_long = np.arange(10 ** 5).astype(np.object)
        # because of nans floats are special:
        self.s_long_floats = Series(np.arange(10 ** 5, dtype=np.float)).astype(
            np.object
        )
        self.vals_long_floats = np.arange(10 ** 5, dtype=np.float).astype(np.object)

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

    params = ["int", "uint", "float", "object"]
    param_names = ["dtype"]

    def setup(self, dtype):
        self.s = Series(np.random.randint(0, 1000, size=100000)).astype(dtype)

    def time_value_counts(self, dtype):
        self.s.value_counts()


class Dir:
    def setup(self):
        self.s = Series(index=tm.makeStringIndex(10000))

    def time_dir_strings(self):
        dir(self.s)


class SeriesGetattr:
    # https://github.com/pandas-dev/pandas/issues/19764
    def setup(self):
        self.s = Series(1, index=date_range("2012-01-01", freq="s", periods=int(1e6)))

    def time_series_datetimeindex_repr(self):
        getattr(self.s, "a", None)


class All:

    params = [[10 ** 3, 10 ** 6], ["fast", "slow"]]
    param_names = ["N", "case"]

    def setup(self, N, case):
        val = case != "fast"
        self.s = Series([val] * N)

    def time_all(self, N, case):
        self.s.all()


class Any:

    params = [[10 ** 3, 10 ** 6], ["fast", "slow"]]
    param_names = ["N", "case"]

    def setup(self, N, case):
        val = case == "fast"
        self.s = Series([val] * N)

    def time_any(self, N, case):
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
        ["int8", "int32", "int64", "float64"],
    ]
    param_names = ["func", "N", "dtype"]

    def setup(self, func, N, dtype):
        self.s = Series([1] * N, dtype=dtype)
        self.func = getattr(self.s, func)

    def time_func(self, func, N, dtype):
        self.func()


from .pandas_vb_common import setup  # noqa: F401 isort:skip
