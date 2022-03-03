import warnings

import numpy as np

import pandas as pd


class Methods:

    params = (
        ["DataFrame", "Series"],
        [("rolling", {"window": 10}), ("rolling", {"window": 1000}), ("expanding", {})],
        ["int", "float"],
        ["median", "mean", "max", "min", "std", "count", "skew", "kurt", "sum", "sem"],
    )
    param_names = ["constructor", "window_kwargs", "dtype", "method"]

    def setup(self, constructor, window_kwargs, dtype, method):
        N = 10**5
        window, kwargs = window_kwargs
        arr = (100 * np.random.random(N)).astype(dtype)
        obj = getattr(pd, constructor)(arr)
        self.window = getattr(obj, window)(**kwargs)

    def time_method(self, constructor, window_kwargs, dtype, method):
        getattr(self.window, method)()

    def peakmem_method(self, constructor, window_kwargs, dtype, method):
        getattr(self.window, method)()


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
        N = 10**3
        arr = (100 * np.random.random(N)).astype(dtype)
        self.roll = getattr(pd, constructor)(arr).rolling(window)

    def time_rolling(self, constructor, window, dtype, function, raw):
        self.roll.apply(function, raw=raw)


class NumbaEngineMethods:
    params = (
        ["DataFrame", "Series"],
        ["int", "float"],
        [("rolling", {"window": 10}), ("expanding", {})],
        ["sum", "max", "min", "median", "mean", "var", "std"],
        [True, False],
        [None, 100],
    )
    param_names = [
        "constructor",
        "dtype",
        "window_kwargs",
        "method",
        "parallel",
        "cols",
    ]

    def setup(self, constructor, dtype, window_kwargs, method, parallel, cols):
        N = 10**3
        window, kwargs = window_kwargs
        shape = (N, cols) if cols is not None and constructor != "Series" else N
        arr = (100 * np.random.random(shape)).astype(dtype)
        data = getattr(pd, constructor)(arr)

        # Warm the cache
        with warnings.catch_warnings(record=True):
            # Catch parallel=True not being applicable e.g. 1D data
            self.window = getattr(data, window)(**kwargs)
            getattr(self.window, method)(
                engine="numba", engine_kwargs={"parallel": parallel}
            )

    def test_method(self, constructor, dtype, window_kwargs, method, parallel, cols):
        with warnings.catch_warnings(record=True):
            getattr(self.window, method)(
                engine="numba", engine_kwargs={"parallel": parallel}
            )


class NumbaEngineApply:
    params = (
        ["DataFrame", "Series"],
        ["int", "float"],
        [("rolling", {"window": 10}), ("expanding", {})],
        [np.sum, lambda x: np.sum(x) + 5],
        [True, False],
        [None, 100],
    )
    param_names = [
        "constructor",
        "dtype",
        "window_kwargs",
        "function",
        "parallel",
        "cols",
    ]

    def setup(self, constructor, dtype, window_kwargs, function, parallel, cols):
        N = 10**3
        window, kwargs = window_kwargs
        shape = (N, cols) if cols is not None and constructor != "Series" else N
        arr = (100 * np.random.random(shape)).astype(dtype)
        data = getattr(pd, constructor)(arr)

        # Warm the cache
        with warnings.catch_warnings(record=True):
            # Catch parallel=True not being applicable e.g. 1D data
            self.window = getattr(data, window)(**kwargs)
            self.window.apply(
                function, raw=True, engine="numba", engine_kwargs={"parallel": parallel}
            )

    def test_method(self, constructor, dtype, window_kwargs, function, parallel, cols):
        with warnings.catch_warnings(record=True):
            self.window.apply(
                function, raw=True, engine="numba", engine_kwargs={"parallel": parallel}
            )


class EWMMethods:

    params = (
        ["DataFrame", "Series"],
        [
            ({"halflife": 10}, "mean"),
            ({"halflife": 10}, "std"),
            ({"halflife": 1000}, "mean"),
            ({"halflife": 1000}, "std"),
            (
                {
                    "halflife": "1 Day",
                    "times": pd.date_range("1900", periods=10**5, freq="23s"),
                },
                "mean",
            ),
        ],
        ["int", "float"],
    )
    param_names = ["constructor", "kwargs_method", "dtype"]

    def setup(self, constructor, kwargs_method, dtype):
        N = 10**5
        kwargs, method = kwargs_method
        arr = (100 * np.random.random(N)).astype(dtype)
        self.method = method
        self.ewm = getattr(pd, constructor)(arr).ewm(**kwargs)

    def time_ewm(self, constructor, kwargs_method, dtype):
        getattr(self.ewm, self.method)()


class VariableWindowMethods(Methods):
    params = (
        ["DataFrame", "Series"],
        ["50s", "1h", "1d"],
        ["int", "float"],
        ["median", "mean", "max", "min", "std", "count", "skew", "kurt", "sum", "sem"],
    )
    param_names = ["constructor", "window", "dtype", "method"]

    def setup(self, constructor, window, dtype, method):
        N = 10**5
        arr = (100 * np.random.random(N)).astype(dtype)
        index = pd.date_range("2017-01-01", periods=N, freq="5s")
        self.window = getattr(pd, constructor)(arr, index=index).rolling(window)


class Pairwise:

    params = (
        [({"window": 10}, "rolling"), ({"window": 1000}, "rolling"), ({}, "expanding")],
        ["corr", "cov"],
        [True, False],
    )
    param_names = ["window_kwargs", "method", "pairwise"]

    def setup(self, kwargs_window, method, pairwise):
        N = 10**4
        n_groups = 20
        kwargs, window = kwargs_window
        groups = [i for _ in range(N // n_groups) for i in range(n_groups)]
        arr = np.random.random(N)
        self.df = pd.DataFrame(arr)
        self.window = getattr(self.df, window)(**kwargs)
        self.window_group = getattr(
            pd.DataFrame({"A": groups, "B": arr}).groupby("A"), window
        )(**kwargs)

    def time_pairwise(self, kwargs_window, method, pairwise):
        getattr(self.window, method)(self.df, pairwise=pairwise)

    def time_groupby(self, kwargs_window, method, pairwise):
        getattr(self.window_group, method)(self.df, pairwise=pairwise)


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
        N = 10**5
        arr = np.random.random(N).astype(dtype)
        self.roll = getattr(pd, constructor)(arr).rolling(window)

    def time_quantile(self, constructor, window, dtype, percentile, interpolation):
        self.roll.quantile(percentile, interpolation=interpolation)


class Rank:
    params = (
        ["DataFrame", "Series"],
        [10, 1000],
        ["int", "float"],
        [True, False],
        [True, False],
        ["min", "max", "average"],
    )
    param_names = [
        "constructor",
        "window",
        "dtype",
        "percentile",
        "ascending",
        "method",
    ]

    def setup(self, constructor, window, dtype, percentile, ascending, method):
        N = 10**5
        arr = np.random.random(N).astype(dtype)
        self.roll = getattr(pd, constructor)(arr).rolling(window)

    def time_rank(self, constructor, window, dtype, percentile, ascending, method):
        self.roll.rank(pct=percentile, ascending=ascending, method=method)


class PeakMemFixedWindowMinMax:

    params = ["min", "max"]

    def setup(self, operation):
        N = 10**6
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
        N = 10**5
        arr = np.random.random(N).astype(dtype)
        indexer = pd.api.indexers.FixedForwardWindowIndexer(window_size=window_size)
        self.roll = getattr(pd, constructor)(arr).rolling(window=indexer)

    def time_rolling(self, constructor, window_size, dtype, method):
        getattr(self.roll, method)()

    def peakmem_rolling(self, constructor, window_size, dtype, method):
        getattr(self.roll, method)()


class Groupby:

    params = (
        ["sum", "median", "mean", "max", "min", "kurt", "sum"],
        [
            ("rolling", {"window": 2}),
            ("rolling", {"window": "30s", "on": "C"}),
            ("expanding", {}),
        ],
    )

    def setup(self, method, window_kwargs):
        N = 1000
        window, kwargs = window_kwargs
        df = pd.DataFrame(
            {
                "A": [str(i) for i in range(N)] * 10,
                "B": list(range(N)) * 10,
                "C": pd.date_range(start="1900-01-01", freq="1min", periods=N * 10),
            }
        )
        self.groupby_window = getattr(df.groupby("A"), window)(**kwargs)

    def time_method(self, method, window_kwargs):
        getattr(self.groupby_window, method)()


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
