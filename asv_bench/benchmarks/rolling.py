import numpy as np

import pandas as pd


class Methods:

    params = (
        ["DataFrame", "Series"],
        [10, 1000],
        ["int", "float"],
        ["median", "mean", "max", "min", "std", "count", "skew", "kurt", "sum"],
    )
    param_names = ["contructor", "window", "dtype", "method"]

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


class ExpandingMethods:

    params = (
        ["DataFrame", "Series"],
        ["int", "float"],
        ["median", "mean", "max", "min", "std", "count", "skew", "kurt", "sum"],
    )
    param_names = ["contructor", "window", "dtype", "method"]

    def setup(self, constructor, dtype, method):
        N = 10 ** 5
        arr = (100 * np.random.random(N)).astype(dtype)
        self.expanding = getattr(pd, constructor)(arr).expanding()

    def time_expanding(self, constructor, dtype, method):
        getattr(self.expanding, method)()


class EWMMethods:

    params = (["DataFrame", "Series"], [10, 1000], ["int", "float"], ["mean", "std"])
    param_names = ["contructor", "window", "dtype", "method"]

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
    param_names = ["contructor", "window", "dtype", "method"]

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


class PeakMemFixed:
    def setup(self):
        N = 10
        arr = 100 * np.random.random(N)
        self.roll = pd.Series(arr).rolling(10)

    def peakmem_fixed(self):
        # GH 25926
        # This is to detect memory leaks in rolling operations.
        # To save time this is only ran on one method.
        # 6000 iterations is enough for most types of leaks to be detected
        for x in range(6000):
            self.roll.max()


from .pandas_vb_common import setup  # noqa: F401 isort:skip
