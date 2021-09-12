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
        ["sum", "max", "min", "median", "mean"],
    )
    param_names = ["constructor", "dtype", "function", "engine", "method"]

    def setup(self, constructor, dtype, function, engine, method):
        N = 10 ** 3
        arr = (100 * np.random.random(N)).astype(dtype)
        self.data = getattr(pd, constructor)(arr)

    def time_rolling_apply(self, constructor, dtype, function, engine, method):
        self.data.rolling(10).apply(function, raw=True, engine=engine)

    def time_expanding_apply(self, constructor, dtype, function, engine, method):
        self.data.expanding().apply(function, raw=True, engine=engine)

    def time_rolling_methods(self, constructor, dtype, function, engine, method):
        getattr(self.data.rolling(10), method)(engine=engine)


class ExpandingMethods:

    params = (
        ["DataFrame", "Series"],
        ["int", "float"],
        ["median", "mean", "max", "min", "std", "count", "skew", "kurt", "sum"],
    )
    param_names = ["constructor", "window", "dtype", "method"]

    def setup(self, constructor, dtype, method):
        N = 10 ** 5
        N_groupby = 100
        arr = (100 * np.random.random(N)).astype(dtype)
        self.expanding = getattr(pd, constructor)(arr).expanding()
        self.expanding_groupby = (
            pd.DataFrame({"A": arr[:N_groupby], "B": range(N_groupby)})
            .groupby("B")
            .expanding()
        )

    def time_expanding(self, constructor, dtype, method):
        getattr(self.expanding, method)()

    def time_expanding_groupby(self, constructor, dtype, method):
        getattr(self.expanding_groupby, method)()


class EWMMethods:

    params = (["DataFrame", "Series"], [10, 1000], ["int", "float"], ["mean", "std"])
    param_names = ["constructor", "window", "dtype", "method"]

    def setup(self, constructor, window, dtype, method):
        N = 10 ** 5
        arr = (100 * np.random.random(N)).astype(dtype)
        times = pd.date_range("1900", periods=N, freq="23s")
        self.ewm = getattr(pd, constructor)(arr).ewm(halflife=window)
        self.ewm_times = getattr(pd, constructor)(arr).ewm(
            halflife="1 Day", times=times
        )

    def time_ewm(self, constructor, window, dtype, method):
        getattr(self.ewm, method)()

    def time_ewm_times(self, constructor, window, dtype, method):
        self.ewm_times.mean()


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
        n_groups = 20
        groups = [i for _ in range(N // n_groups) for i in range(n_groups)]
        arr = np.random.random(N)
        self.df = pd.DataFrame(arr)
        self.df_group = pd.DataFrame({"A": groups, "B": arr}).groupby("A")

    def time_pairwise(self, window, method, pairwise):
        if window is None:
            r = self.df.expanding()
        else:
            r = self.df.rolling(window=window)
        getattr(r, method)(self.df, pairwise=pairwise)

    def time_groupby(self, window, method, pairwise):
        if window is None:
            r = self.df_group.expanding()
        else:
            r = self.df_group.rolling(window=window)
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
        N = 10 ** 6
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


class GroupbyLargeGroups:
    # https://github.com/pandas-dev/pandas/issues/38038
    # specific example where the rolling operation on a larger dataframe
    # is relatively cheap (few but large groups), but creation of
    # MultiIndex of result can be expensive

    def setup(self):
        N = 100000
        self.df = pd.DataFrame({"A": [1, 2] * (N // 2), "B": np.random.randn(N)})

    def time_rolling_multiindex_creation(self):
        self.df.groupby("A").rolling(3).mean()


class GroupbyEWM:

    params = ["var", "std", "cov", "corr"]
    param_names = ["method"]

    def setup(self, method):
        df = pd.DataFrame({"A": range(50), "B": range(50)})
        self.gb_ewm = df.groupby("A").ewm(com=1.0)

    def time_groupby_method(self, method):
        getattr(self.gb_ewm, method)()


class GroupbyEWMEngine:

    params = ["cython", "numba"]
    param_names = ["engine"]

    def setup(self, engine):
        df = pd.DataFrame({"A": range(50), "B": range(50)})
        self.gb_ewm = df.groupby("A").ewm(com=1.0)

    def time_groupby_mean(self, engine):
        self.gb_ewm.mean(engine=engine)


def table_method_func(x):
    return np.sum(x, axis=0) + 1


class TableMethod:

    params = ["single", "table"]
    param_names = ["method"]

    def setup(self, method):
        self.df = pd.DataFrame(np.random.randn(10, 1000))

    def time_apply(self, method):
        self.df.rolling(2, method=method).apply(
            table_method_func, raw=True, engine="numba"
        )

    def time_ewm_mean(self, method):
        self.df.ewm(1, method=method).mean(engine="numba")


from .pandas_vb_common import setup  # noqa: F401 isort:skip
