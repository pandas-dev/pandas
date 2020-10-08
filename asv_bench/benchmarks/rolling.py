import numpy as np

import pandas as pd


class Methods:

    params = (
        ["DataFrame", "Series"],
        [10, 1000],
        ["int", "float"],
        ["median", "mean", "max", "min", "std", "count", "skew", "kurt", "sum"],
    )
    param_names = ["constructor", "window", "dtype", "method"]

    def setup(self, constructor, window, dtype, method):
        N = 10 ** 5
        arr = (100 * np.random.random(N)).astype(dtype)
        self.roll = getattr(pd, constructor)(arr).rolling(window)

    def time_rolling(self, constructor, window, dtype, method):
        getattr(self.roll, method)()

    def peakmem_rolling(self, constructor, window, dtype, method):
        getattr(self.roll, method)()


class Apply:
    params = (
        ["DataFrame", "Series"],
        [3, 300],
        ["int", "float"],
        [sum, np.sum, lambda x: np.sum(x) + 5],
        [True, False],
    )
    param_names = ["constructor", "window", "dtype", "function", "raw"]

    def setup(self, constructor, window, dtype, function, raw):
        N = 10 ** 3
        arr = (100 * np.random.random(N)).astype(dtype)
        self.roll = getattr(pd, constructor)(arr).rolling(window)

    def time_rolling(self, constructor, window, dtype, function, raw):
        self.roll.apply(function, raw=raw)


class Engine:
    params = (
        ["DataFrame", "Series"],
        ["int", "float"],
        [np.sum, lambda x: np.sum(x) + 5],
        ["cython", "numba"],
    )
    param_names = ["constructor", "dtype", "function", "engine"]

    def setup(self, constructor, dtype, function, engine):
        N = 10 ** 3
        arr = (100 * np.random.random(N)).astype(dtype)
        self.data = getattr(pd, constructor)(arr)

    def time_rolling_apply(self, constructor, dtype, function, engine):
        self.data.rolling(10).apply(function, raw=True, engine=engine)

    def time_expanding_apply(self, constructor, dtype, function, engine):
        self.data.expanding().apply(function, raw=True, engine=engine)


class ExpandingMethods:

    params = (
        ["DataFrame", "Series"],
        ["int", "float"],
        ["median", "mean", "max", "min", "std", "count", "skew", "kurt", "sum"],
    )
    param_names = ["constructor", "window", "dtype", "method"]

    def setup(self, constructor, dtype, method):
        N = 10 ** 5
        arr = (100 * np.random.random(N)).astype(dtype)
        self.expanding = getattr(pd, constructor)(arr).expanding()

    def time_expanding(self, constructor, dtype, method):
        getattr(self.expanding, method)()


class EWMMethods:

    params = (["DataFrame", "Series"], [10, 1000], ["int", "float"], ["mean", "std"])
    param_names = ["constructor", "window", "dtype", "method"]

    def setup(self, constructor, window, dtype, method):
        N = 10 ** 5
        arr = (100 * np.random.random(N)).astype(dtype)
        self.ewm = getattr(pd, constructor)(arr).ewm(halflife=window)

    def time_ewm(self, constructor, window, dtype, method):
        getattr(self.ewm, method)()


class VariableWindowMethods(Methods):
    params = (
        ["DataFrame", "Series"],
        ["50s", "1h", "1d"],
        ["int", "float"],
        ["median", "mean", "max", "min", "std", "count", "skew", "kurt", "sum"],
    )
    param_names = ["constructor", "window", "dtype", "method"]

    def setup(self, constructor, window, dtype, method):
        N = 10 ** 5
        arr = (100 * np.random.random(N)).astype(dtype)
        index = pd.date_range("2017-01-01", periods=N, freq="5s")
        self.roll = getattr(pd, constructor)(arr, index=index).rolling(window)


class Pairwise:

    params = ([10, 1000, None], ["corr", "cov"], [True, False])
    param_names = ["window", "method", "pairwise"]

    def setup(self, window, method, pairwise):
        N = 10 ** 4
        arr = np.random.random(N)
        self.df = pd.DataFrame(arr)

    def time_pairwise(self, window, method, pairwise):
        if window is None:
            r = self.df.expanding()
        else:
            r = self.df.rolling(window=window)
        getattr(r, method)(self.df, pairwise=pairwise)


class Quantile:
    params = (
        ["DataFrame", "Series"],
        [10, 1000],
        ["int", "float"],
        [0, 0.5, 1],
        ["linear", "nearest", "lower", "higher", "midpoint"],
    )
    param_names = ["constructor", "window", "dtype", "percentile"]

    def setup(self, constructor, window, dtype, percentile, interpolation):
        N = 10 ** 5
        arr = np.random.random(N).astype(dtype)
        self.roll = getattr(pd, constructor)(arr).rolling(window)

    def time_quantile(self, constructor, window, dtype, percentile, interpolation):
        self.roll.quantile(percentile, interpolation=interpolation)


class PeakMemFixedWindowMinMax:

    params = ["min", "max"]

    def setup(self, operation):
        N = int(1e6)
        arr = np.random.random(N)
        self.roll = pd.Series(arr).rolling(2)

    def peakmem_fixed(self, operation):
        for x in range(5):
            getattr(self.roll, operation)()


class ForwardWindowMethods:
    params = (
        ["DataFrame", "Series"],
        [10, 1000],
        ["int", "float"],
        ["median", "mean", "max", "min", "kurt", "sum"],
    )
    param_names = ["constructor", "window_size", "dtype", "method"]

    def setup(self, constructor, window_size, dtype, method):
        N = 10 ** 5
        arr = np.random.random(N).astype(dtype)
        indexer = pd.api.indexers.FixedForwardWindowIndexer(window_size=window_size)
        self.roll = getattr(pd, constructor)(arr).rolling(window=indexer)

    def time_rolling(self, constructor, window_size, dtype, method):
        getattr(self.roll, method)()

    def peakmem_rolling(self, constructor, window_size, dtype, method):
        getattr(self.roll, method)()


class Groupby:

    params = ["sum", "median", "mean", "max", "min", "kurt", "sum"]

    def setup(self, method):
        N = 1000
        df = pd.DataFrame(
            {
                "A": [str(i) for i in range(N)] * 10,
                "B": list(range(N)) * 10,
                "C": pd.date_range(start="1900-01-01", freq="1min", periods=N * 10),
            }
        )
        self.groupby_roll_int = df.groupby("A").rolling(window=2)
        self.groupby_roll_offset = df.groupby("A").rolling(window="30s", on="C")

    def time_rolling_int(self, method):
        getattr(self.groupby_roll_int, method)()

    def time_rolling_offset(self, method):
        getattr(self.groupby_roll_offset, method)()


from .pandas_vb_common import setup  # noqa: F401 isort:skip
