import importlib
import sys

import matplotlib
import numpy as np
import pkg_resources

from pandas import (
    DataFrame,
    DatetimeIndex,
    Series,
    date_range,
)

try:
    from pandas.plotting import andrews_curves
except ImportError:
    from pandas.tools.plotting import andrews_curves

from pandas.plotting._core import _get_plot_backend

matplotlib.use("Agg")


class SeriesPlotting:
    params = [["line", "bar", "area", "barh", "hist", "kde", "pie"]]
    param_names = ["kind"]

    def setup(self, kind):
        if kind in ["bar", "barh", "pie"]:
            n = 100
        elif kind in ["kde"]:
            n = 10000
        else:
            n = 1000000

        self.s = Series(np.random.randn(n))
        if kind in ["area", "pie"]:
            self.s = self.s.abs()

    def time_series_plot(self, kind):
        self.s.plot(kind=kind)


class FramePlotting:
    params = [
        ["line", "bar", "area", "barh", "hist", "kde", "pie", "scatter", "hexbin"]
    ]
    param_names = ["kind"]

    def setup(self, kind):
        if kind in ["bar", "barh", "pie"]:
            n = 100
        elif kind in ["kde", "scatter", "hexbin"]:
            n = 10000
        else:
            n = 1000000

        self.x = Series(np.random.randn(n))
        self.y = Series(np.random.randn(n))
        if kind in ["area", "pie"]:
            self.x = self.x.abs()
            self.y = self.y.abs()
        self.df = DataFrame({"x": self.x, "y": self.y})

    def time_frame_plot(self, kind):
        self.df.plot(x="x", y="y", kind=kind)


class TimeseriesPlotting:
    def setup(self):
        N = 2000
        M = 5
        idx = date_range("1/1/1975", periods=N)
        self.df = DataFrame(np.random.randn(N, M), index=idx)

        idx_irregular = DatetimeIndex(
            np.concatenate((idx.values[0:10], idx.values[12:]))
        )
        self.df2 = DataFrame(
            np.random.randn(len(idx_irregular), M), index=idx_irregular
        )

    def time_plot_regular(self):
        self.df.plot()

    def time_plot_regular_compat(self):
        self.df.plot(x_compat=True)

    def time_plot_irregular(self):
        self.df2.plot()

    def time_plot_table(self):
        self.df.plot(table=True)


class Misc:
    def setup(self):
        N = 500
        M = 10
        self.df = DataFrame(np.random.randn(N, M))
        self.df["Name"] = ["A"] * N

    def time_plot_andrews_curves(self):
        andrews_curves(self.df, "Name")


class BackendLoading:
    repeat = 1
    number = 1
    warmup_time = 0

    def setup(self):
        dist = pkg_resources.get_distribution("pandas")
        spec = importlib.machinery.ModuleSpec("my_backend", None)
        mod = importlib.util.module_from_spec(spec)
        mod.plot = lambda *args, **kwargs: 1

        backends = pkg_resources.get_entry_map("pandas")
        my_entrypoint = pkg_resources.EntryPoint(
            "pandas_plotting_backend", mod.__name__, dist=dist
        )
        backends["pandas_plotting_backends"][mod.__name__] = my_entrypoint
        for i in range(10):
            backends["pandas_plotting_backends"][str(i)] = my_entrypoint
        sys.modules["my_backend"] = mod

    def time_get_plot_backend(self):
        _get_plot_backend("my_backend")


from .pandas_vb_common import setup  # noqa: F401 isort:skip
