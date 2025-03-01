import contextlib
import importlib.machinery
import importlib.util
import os
import pathlib
import sys
import tempfile
from unittest import mock

import matplotlib
import numpy as np

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
        mod = importlib.util.module_from_spec(
            importlib.machinery.ModuleSpec("pandas_dummy_backend", None)
        )
        mod.plot = lambda *args, **kwargs: 1

        with contextlib.ExitStack() as stack:
            stack.enter_context(
                mock.patch.dict(sys.modules, {"pandas_dummy_backend": mod})
            )
            tmp_path = pathlib.Path(stack.enter_context(tempfile.TemporaryDirectory()))

            sys.path.insert(0, os.fsdecode(tmp_path))
            stack.callback(sys.path.remove, os.fsdecode(tmp_path))

            dist_info = tmp_path / "my_backend-0.0.0.dist-info"
            dist_info.mkdir()
            (dist_info / "entry_points.txt").write_bytes(
                b"[pandas_plotting_backends]\n"
                b"my_ep_backend = pandas_dummy_backend\n"
                b"my_ep_backend0 = pandas_dummy_backend\n"
                b"my_ep_backend1 = pandas_dummy_backend\n"
                b"my_ep_backend2 = pandas_dummy_backend\n"
                b"my_ep_backend3 = pandas_dummy_backend\n"
                b"my_ep_backend4 = pandas_dummy_backend\n"
                b"my_ep_backend5 = pandas_dummy_backend\n"
                b"my_ep_backend6 = pandas_dummy_backend\n"
                b"my_ep_backend7 = pandas_dummy_backend\n"
                b"my_ep_backend8 = pandas_dummy_backend\n"
                b"my_ep_backend9 = pandas_dummy_backend\n"
            )
            self.stack = stack.pop_all()

    def teardown(self):
        self.stack.close()

    def time_get_plot_backend(self):
        # finds the first my_ep_backend
        _get_plot_backend("my_ep_backend")

    def time_get_plot_backend_fallback(self):
        # iterates through all the my_ep_backend[0-9] before falling back
        # to importlib.import_module
        _get_plot_backend("pandas_dummy_backend")


from .pandas_vb_common import setup  # noqa: F401 isort:skip
