from importlib import import_module

import numpy as np

import pandas as pd

from .pandas_vb_common import tm

for imp in ["pandas.util", "pandas.tools.hashing"]:
    try:
        hashing = import_module(imp)
        break
    except (ImportError, TypeError, ValueError):
        pass


class Factorize:

    params = [
        [True, False],
        [True, False],
        [
            "int",
            "uint",
            "float",
            "object",
            "datetime64[ns]",
            "datetime64[ns, tz]",
            "Int64",
            "boolean",
            "string[pyarrow]",
        ],
    ]
    param_names = ["unique", "sort", "dtype"]

    def setup(self, unique, sort, dtype):
        N = 10 ** 5
        string_index = tm.makeStringIndex(N)
        string_arrow = None
        if dtype == "string[pyarrow]":
            try:
                string_arrow = pd.array(string_index, dtype="string[pyarrow]")
            except ImportError:
                raise NotImplementedError

        data = {
            "int": pd.Int64Index(np.arange(N)),
            "uint": pd.UInt64Index(np.arange(N)),
            "float": pd.Float64Index(np.random.randn(N)),
            "object": string_index,
            "datetime64[ns]": pd.date_range("2011-01-01", freq="H", periods=N),
            "datetime64[ns, tz]": pd.date_range(
                "2011-01-01", freq="H", periods=N, tz="Asia/Tokyo"
            ),
            "Int64": pd.array(np.arange(N), dtype="Int64"),
            "boolean": pd.array(np.random.randint(0, 2, N), dtype="boolean"),
            "string[pyarrow]": string_arrow,
        }[dtype]
        if not unique:
            data = data.repeat(5)
        self.data = data

    def time_factorize(self, unique, sort, dtype):
        pd.factorize(self.data, sort=sort)


class Duplicated:

    params = [
        [True, False],
        ["first", "last", False],
        ["int", "uint", "float", "string", "datetime64[ns]", "datetime64[ns, tz]"],
    ]
    param_names = ["unique", "keep", "dtype"]

    def setup(self, unique, keep, dtype):
        N = 10 ** 5
        data = {
            "int": pd.Int64Index(np.arange(N)),
            "uint": pd.UInt64Index(np.arange(N)),
            "float": pd.Float64Index(np.random.randn(N)),
            "string": tm.makeStringIndex(N),
            "datetime64[ns]": pd.date_range("2011-01-01", freq="H", periods=N),
            "datetime64[ns, tz]": pd.date_range(
                "2011-01-01", freq="H", periods=N, tz="Asia/Tokyo"
            ),
        }[dtype]
        if not unique:
            data = data.repeat(5)
        self.idx = data
        # cache is_unique
        self.idx.is_unique

    def time_duplicated(self, unique, keep, dtype):
        self.idx.duplicated(keep=keep)


class Hashing:
    def setup_cache(self):
        N = 10 ** 5

        df = pd.DataFrame(
            {
                "strings": pd.Series(
                    tm.makeStringIndex(10000).take(np.random.randint(0, 10000, size=N))
                ),
                "floats": np.random.randn(N),
                "ints": np.arange(N),
                "dates": pd.date_range("20110101", freq="s", periods=N),
                "timedeltas": pd.timedelta_range("1 day", freq="s", periods=N),
            }
        )
        df["categories"] = df["strings"].astype("category")
        df.iloc[10:20] = np.nan
        return df

    def time_frame(self, df):
        hashing.hash_pandas_object(df)

    def time_series_int(self, df):
        hashing.hash_pandas_object(df["ints"])

    def time_series_string(self, df):
        hashing.hash_pandas_object(df["strings"])

    def time_series_float(self, df):
        hashing.hash_pandas_object(df["floats"])

    def time_series_categorical(self, df):
        hashing.hash_pandas_object(df["categories"])

    def time_series_timedeltas(self, df):
        hashing.hash_pandas_object(df["timedeltas"])

    def time_series_dates(self, df):
        hashing.hash_pandas_object(df["dates"])


class Quantile:
    params = [
        [0, 0.5, 1],
        ["linear", "nearest", "lower", "higher", "midpoint"],
        ["float", "int", "uint"],
    ]
    param_names = ["quantile", "interpolation", "dtype"]

    def setup(self, quantile, interpolation, dtype):
        N = 10 ** 5
        data = {
            "int": np.arange(N),
            "uint": np.arange(N).astype(np.uint64),
            "float": np.random.randn(N),
        }
        self.idx = pd.Series(data[dtype].repeat(5))

    def time_quantile(self, quantile, interpolation, dtype):
        self.idx.quantile(quantile, interpolation=interpolation)


class SortIntegerArray:
    params = [10 ** 3, 10 ** 5]

    def setup(self, N):
        data = np.arange(N, dtype=float)
        data[40] = np.nan
        self.array = pd.array(data, dtype="Int64")

    def time_argsort(self, N):
        self.array.argsort()


from .pandas_vb_common import setup  # noqa: F401 isort:skip
