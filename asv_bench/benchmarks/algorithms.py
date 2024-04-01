from importlib import import_module

import numpy as np

import pandas as pd

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
            "int64",
            "uint64",
            "float64",
            "object",
            "object_str",
            "datetime64[ns]",
            "datetime64[ns, tz]",
            "Int64",
            "boolean",
            "string[pyarrow]",
        ],
    ]
    param_names = ["unique", "sort", "dtype"]

    def setup(self, unique, sort, dtype):
        N = 10**5

        if dtype in ["int64", "uint64", "Int64", "object"]:
            data = pd.Index(np.arange(N), dtype=dtype)
        elif dtype == "float64":
            data = pd.Index(np.random.randn(N), dtype=dtype)
        elif dtype == "boolean":
            data = pd.array(np.random.randint(0, 2, N), dtype=dtype)
        elif dtype == "datetime64[ns]":
            data = pd.date_range("2011-01-01", freq="h", periods=N)
        elif dtype == "datetime64[ns, tz]":
            data = pd.date_range("2011-01-01", freq="h", periods=N, tz="Asia/Tokyo")
        elif dtype == "object_str":
            data = pd.Index([f"i-{i}" for i in range(N)], dtype=object)
        elif dtype == "string[pyarrow]":
            data = pd.array(
                pd.Index([f"i-{i}" for i in range(N)], dtype=object),
                dtype="string[pyarrow]",
            )
        else:
            raise NotImplementedError

        if not unique:
            data = data.repeat(5)
        self.data = data

    def time_factorize(self, unique, sort, dtype):
        pd.factorize(self.data, sort=sort)

    def peakmem_factorize(self, unique, sort, dtype):
        pd.factorize(self.data, sort=sort)


class Duplicated:
    params = [
        [True, False],
        ["first", "last", False],
        [
            "int64",
            "uint64",
            "float64",
            "string",
            "datetime64[ns]",
            "datetime64[ns, tz]",
            "timestamp[ms][pyarrow]",
            "duration[s][pyarrow]",
        ],
    ]
    param_names = ["unique", "keep", "dtype"]

    def setup(self, unique, keep, dtype):
        N = 10**5
        if dtype in ["int64", "uint64"]:
            data = pd.Index(np.arange(N), dtype=dtype)
        elif dtype == "float64":
            data = pd.Index(np.random.randn(N), dtype="float64")
        elif dtype == "string":
            data = pd.Index([f"i-{i}" for i in range(N)], dtype=object)
        elif dtype == "datetime64[ns]":
            data = pd.date_range("2011-01-01", freq="h", periods=N)
        elif dtype == "datetime64[ns, tz]":
            data = pd.date_range("2011-01-01", freq="h", periods=N, tz="Asia/Tokyo")
        elif dtype in ["timestamp[ms][pyarrow]", "duration[s][pyarrow]"]:
            data = pd.Index(np.arange(N), dtype=dtype)
        else:
            raise NotImplementedError
        if not unique:
            data = data.repeat(5)
        self.idx = data
        # cache is_unique
        self.idx.is_unique

    def time_duplicated(self, unique, keep, dtype):
        self.idx.duplicated(keep=keep)


class DuplicatedMaskedArray:
    params = [
        [True, False],
        ["first", "last", False],
        ["Int64", "Float64"],
    ]
    param_names = ["unique", "keep", "dtype"]

    def setup(self, unique, keep, dtype):
        N = 10**5
        data = pd.Series(np.arange(N), dtype=dtype)
        data[list(range(1, N, 100))] = pd.NA
        if not unique:
            data = data.repeat(5)
        self.ser = data
        # cache is_unique
        self.ser.is_unique

    def time_duplicated(self, unique, keep, dtype):
        self.ser.duplicated(keep=keep)


class Hashing:
    def setup_cache(self):
        N = 10**5

        df = pd.DataFrame(
            {
                "strings": pd.Series(
                    pd.Index([f"i-{i}" for i in range(10000)], dtype=object).take(
                        np.random.randint(0, 10000, size=N)
                    )
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
        ["float64", "int64", "uint64"],
    ]
    param_names = ["quantile", "interpolation", "dtype"]

    def setup(self, quantile, interpolation, dtype):
        N = 10**5
        if dtype in ["int64", "uint64"]:
            data = np.arange(N, dtype=dtype)
        elif dtype == "float64":
            data = np.random.randn(N)
        else:
            raise NotImplementedError
        self.ser = pd.Series(data.repeat(5))

    def time_quantile(self, quantile, interpolation, dtype):
        self.ser.quantile(quantile, interpolation=interpolation)


class SortIntegerArray:
    params = [10**3, 10**5]

    def setup(self, N):
        data = np.arange(N, dtype=float)
        data[40] = np.nan
        self.array = pd.array(data, dtype="Int64")

    def time_argsort(self, N):
        self.array.argsort()


from .pandas_vb_common import setup  # noqa: F401 isort:skip
