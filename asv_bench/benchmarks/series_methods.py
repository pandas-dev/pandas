from datetime import datetime

import numpy as np

from pandas import (
    NA,
    Index,
    NaT,
    Series,
    date_range,
)


class SeriesConstructor:
    def setup(self):
        self.idx = date_range(
            start=datetime(2015, 10, 26), end=datetime(2016, 1, 1), freq="50s"
        )
        self.data = dict(zip(self.idx, range(len(self.idx))))
        self.array = np.array([1, 2, 3])
        self.idx2 = Index(["a", "b", "c"])

    def time_constructor_dict(self):
        Series(data=self.data, index=self.idx)

    def time_constructor_no_data(self):
        Series(data=None, index=self.idx)


class ToFrame:
    params = [["int64", "datetime64[ns]", "category", "Int64"], [None, "foo"]]
    param_names = ["dtype", "name"]

    def setup(self, dtype, name):
        arr = np.arange(10**5)
        ser = Series(arr, dtype=dtype)
        self.ser = ser

    def time_to_frame(self, dtype, name):
        self.ser.to_frame(name)


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
        N = 10**6
        data = {
            "int": np.random.randint(1, 10, N),
            "datetime": date_range("2000-01-01", freq="s", periods=N),
        }
        self.s = Series(data[dtype])
        if dtype == "datetime":
            self.s[np.random.randint(1, N, 100)] = NaT

    def time_dropna(self, dtype):
        self.s.dropna()


class Fillna:
    params = [
        [
            "datetime64[ns]",
            "float32",
            "float64",
            "Float64",
            "Int64",
            "int64[pyarrow]",
            "string",
            "string[pyarrow]",
        ],
    ]
    param_names = ["dtype"]

    def setup(self, dtype):
        N = 10**6
        if dtype == "datetime64[ns]":
            data = date_range("2000-01-01", freq="s", periods=N)
            na_value = NaT
        elif dtype in ("float64", "Float64"):
            data = np.random.randn(N)
            na_value = np.nan
        elif dtype in ("Int64", "int64[pyarrow]"):
            data = np.arange(N)
            na_value = NA
        elif dtype in ("string", "string[pyarrow]"):
            data = np.array([str(i) * 5 for i in range(N)], dtype=object)
            na_value = NA
        else:
            raise NotImplementedError
        fill_value = data[0]
        ser = Series(data, dtype=dtype)
        ser[::2] = na_value
        self.ser = ser
        self.fill_value = fill_value

    def time_fillna(self, dtype):
        self.ser.fillna(value=self.fill_value)

    def time_ffill(self, dtype):
        self.ser.ffill()

    def time_bfill(self, dtype):
        self.ser.bfill()


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
        N = 10**5
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
    params = [50, 1000, 10**5]
    param_names = ["n"]

    def setup(self, n):
        self.s = Series(np.random.randn(n))

    def time_clip(self, n):
        self.s.clip(0, 1)


class ClipDt:
    def setup(self):
        dr = date_range("20220101", periods=100_000, freq="s", tz="UTC")
        self.clipper_dt = dr[0:1_000].repeat(100)
        self.s = Series(dr)

    def time_clip(self):
        self.s.clip(upper=self.clipper_dt)


class ValueCounts:
    params = [[10**3, 10**4, 10**5], ["int", "uint", "float", "object"]]
    param_names = ["N", "dtype"]

    def setup(self, N, dtype):
        self.s = Series(np.random.randint(0, N, size=10 * N)).astype(dtype)

    def time_value_counts(self, N, dtype):
        self.s.value_counts()


class ValueCountsEA:
    params = [[10**3, 10**4, 10**5], [True, False]]
    param_names = ["N", "dropna"]

    def setup(self, N, dropna):
        self.s = Series(np.random.randint(0, N, size=10 * N), dtype="Int64")
        self.s.loc[1] = NA

    def time_value_counts(self, N, dropna):
        self.s.value_counts(dropna=dropna)


class ValueCountsObjectDropNAFalse:
    params = [10**3, 10**4, 10**5]
    param_names = ["N"]

    def setup(self, N):
        self.s = Series(np.random.randint(0, N, size=10 * N)).astype("object")

    def time_value_counts(self, N):
        self.s.value_counts(dropna=False)


class Mode:
    params = [[10**3, 10**4, 10**5], ["int", "uint", "float", "object"]]
    param_names = ["N", "dtype"]

    def setup(self, N, dtype):
        self.s = Series(np.random.randint(0, N, size=10 * N)).astype(dtype)

    def time_mode(self, N, dtype):
        self.s.mode()


class ModeObjectDropNAFalse:
    params = [10**3, 10**4, 10**5]
    param_names = ["N"]

    def setup(self, N):
        self.s = Series(np.random.randint(0, N, size=10 * N)).astype("object")

    def time_mode(self, N):
        self.s.mode(dropna=False)


class Dir:
    def setup(self):
        self.s = Series(index=Index([f"i-{i}" for i in range(10000)], dtype=object))

    def time_dir_strings(self):
        dir(self.s)


class SeriesGetattr:
    # https://github.com/pandas-dev/pandas/issues/19764
    def setup(self):
        self.s = Series(1, index=date_range("2012-01-01", freq="s", periods=10**6))

    def time_series_datetimeindex_repr(self):
        getattr(self.s, "a", None)


class All:
    params = [[10**3, 10**6], ["fast", "slow"], ["bool", "boolean"]]
    param_names = ["N", "case", "dtype"]

    def setup(self, N, case, dtype):
        val = case != "fast"
        self.s = Series([val] * N, dtype=dtype)

    def time_all(self, N, case, dtype):
        self.s.all()


class Any:
    params = [[10**3, 10**6], ["fast", "slow"], ["bool", "boolean"]]
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
        [10**3, 10**6],
        ["int8", "int32", "int64", "float64", "Int64", "boolean"],
    ]
    param_names = ["func", "N", "dtype"]

    def setup(self, func, N, dtype):
        if func == "argmax" and dtype in {"Int64", "boolean"}:
            # Skip argmax for nullable int since this doesn't work yet (GH-24382)
            raise NotImplementedError
        self.s = Series(np.ones(N), dtype=dtype)
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


class Iter:
    param_names = ["dtype"]
    params = [
        "bool",
        "boolean",
        "int64",
        "Int64",
        "float64",
        "Float64",
        "datetime64[ns]",
    ]

    def setup(self, dtype):
        N = 10**5
        if dtype in ["bool", "boolean"]:
            data = np.repeat([True, False], N // 2)
        elif dtype in ["int64", "Int64"]:
            data = np.arange(N)
        elif dtype in ["float64", "Float64"]:
            data = np.random.randn(N)
        elif dtype == "datetime64[ns]":
            data = date_range("2000-01-01", freq="s", periods=N)
        else:
            raise NotImplementedError

        self.s = Series(data, dtype=dtype)

    def time_iter(self, dtype):
        for v in self.s:
            pass


class ToNumpy:
    def setup(self):
        N = 1_000_000
        self.ser = Series(
            np.random.randn(
                N,
            )
        )

    def time_to_numpy(self):
        self.ser.to_numpy()

    def time_to_numpy_double_copy(self):
        self.ser.to_numpy(dtype="float64", copy=True)

    def time_to_numpy_copy(self):
        self.ser.to_numpy(copy=True)

    def time_to_numpy_float_with_nan(self):
        self.ser.to_numpy(dtype="float64", na_value=np.nan)


class Replace:
    param_names = ["num_to_replace"]
    params = [100, 1000]

    def setup(self, num_to_replace):
        N = 1_000_000
        self.arr = np.random.randn(N)
        self.arr1 = self.arr.copy()
        np.random.shuffle(self.arr1)
        self.ser = Series(self.arr)

        self.to_replace_list = np.random.choice(self.arr, num_to_replace)
        self.values_list = np.random.choice(self.arr1, num_to_replace)

        self.replace_dict = dict(zip(self.to_replace_list, self.values_list))

    def time_replace_dict(self, num_to_replace):
        self.ser.replace(self.replace_dict)

    def peakmem_replace_dict(self, num_to_replace):
        self.ser.replace(self.replace_dict)

    def time_replace_list(self, num_to_replace):
        self.ser.replace(self.to_replace_list, self.values_list)

    def peakmem_replace_list(self, num_to_replace):
        self.ser.replace(self.to_replace_list, self.values_list)


from .pandas_vb_common import setup  # noqa: F401 isort:skip
