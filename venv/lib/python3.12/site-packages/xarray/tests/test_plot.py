from __future__ import annotations

import contextlib
import inspect
import math
from collections.abc import Callable, Generator, Hashable
from copy import copy
from datetime import date, timedelta
from typing import Any, Literal, cast

import numpy as np
import pandas as pd
import pytest

import xarray as xr
import xarray.plot as xplt
from xarray import DataArray, Dataset
from xarray.namedarray.utils import module_available
from xarray.plot.dataarray_plot import _infer_interval_breaks
from xarray.plot.dataset_plot import _infer_meta_data
from xarray.plot.utils import (
    _assert_valid_xy,
    _build_discrete_cmap,
    _color_palette,
    _determine_cmap_params,
    _maybe_gca,
    get_axis,
    label_from_attrs,
)
from xarray.tests import (
    assert_array_equal,
    assert_equal,
    assert_no_warnings,
    requires_cartopy,
    requires_cftime,
    requires_dask,
    requires_matplotlib,
    requires_seaborn,
)

# this should not be imported to test if the automatic lazy import works
has_nc_time_axis = module_available("nc_time_axis")

# import mpl and change the backend before other mpl imports
try:
    import matplotlib as mpl
    import matplotlib.dates
    import matplotlib.pyplot as plt
    import mpl_toolkits
except ImportError:
    pass

with contextlib.suppress(ImportError):
    import cartopy


@contextlib.contextmanager
def figure_context(*args, **kwargs):
    """context manager which autocloses a figure (even if the test failed)"""

    try:
        yield None
    finally:
        plt.close("all")


@pytest.fixture(autouse=True)
def test_all_figures_closed():
    """meta-test to ensure all figures are closed at the end of a test

    Notes:  Scope is kept to module (only invoke this function once per test
    module) else tests cannot be run in parallel (locally). Disadvantage: only
    catches one open figure per run. May still give a false positive if tests
    are run in parallel.
    """
    yield None

    open_figs = len(plt.get_fignums())
    if open_figs:
        raise RuntimeError(
            f"tests did not close all figures ({open_figs} figures open)"
        )


@pytest.mark.flaky
@pytest.mark.skip(reason="maybe flaky")
def text_in_fig() -> set[str]:
    """
    Return the set of all text in the figure
    """
    return {t.get_text() for t in plt.gcf().findobj(mpl.text.Text)}


def find_possible_colorbars() -> list[mpl.collections.QuadMesh]:
    # nb. this function also matches meshes from pcolormesh
    return plt.gcf().findobj(mpl.collections.QuadMesh)


def substring_in_axes(substring: str, ax: mpl.axes.Axes) -> bool:
    """
    Return True if a substring is found anywhere in an axes
    """
    alltxt: set[str] = {t.get_text() for t in ax.findobj(mpl.text.Text)}
    return any(substring in txt for txt in alltxt)


def substring_not_in_axes(substring: str, ax: mpl.axes.Axes) -> bool:
    """
    Return True if a substring is not found anywhere in an axes
    """
    alltxt: set[str] = {t.get_text() for t in ax.findobj(mpl.text.Text)}
    check = [(substring not in txt) for txt in alltxt]
    return all(check)


def property_in_axes_text(
    property, property_str, target_txt, ax: mpl.axes.Axes
) -> bool:
    """
    Return True if the specified text in an axes
    has the property assigned to property_str
    """
    alltxt: list[mpl.text.Text] = ax.findobj(mpl.text.Text)
    return all(
        plt.getp(t, property) == property_str
        for t in alltxt
        if t.get_text() == target_txt
    )


def easy_array(shape: tuple[int, ...], start: float = 0, stop: float = 1) -> np.ndarray:
    """
    Make an array with desired shape using np.linspace

    shape is a tuple like (2, 3)
    """
    a = np.linspace(start, stop, num=math.prod(shape))
    return a.reshape(shape)


def get_colorbar_label(colorbar) -> str:
    if colorbar.orientation == "vertical":
        return colorbar.ax.get_ylabel()
    else:
        return colorbar.ax.get_xlabel()


@requires_matplotlib
class PlotTestCase:
    @pytest.fixture(autouse=True)
    def setup(self) -> Generator:
        yield
        # Remove all matplotlib figures
        plt.close("all")

    def pass_in_axis(self, plotmethod, subplot_kw=None) -> None:
        fig, axs = plt.subplots(ncols=2, subplot_kw=subplot_kw, squeeze=False)
        ax = axs[0, 0]
        plotmethod(ax=ax)
        assert ax.has_data()

    @pytest.mark.slow
    def imshow_called(self, plotmethod) -> bool:
        plotmethod()
        images = plt.gca().findobj(mpl.image.AxesImage)
        return len(images) > 0

    def contourf_called(self, plotmethod) -> bool:
        plotmethod()

        # Compatible with mpl before (PathCollection) and after (QuadContourSet) 3.8
        def matchfunc(x) -> bool:
            return isinstance(
                x, mpl.collections.PathCollection | mpl.contour.QuadContourSet
            )

        paths = plt.gca().findobj(matchfunc)
        return len(paths) > 0


class TestPlot(PlotTestCase):
    @pytest.fixture(autouse=True)
    def setup_array(self) -> None:
        self.darray = DataArray(easy_array((2, 3, 4)))

    def test_accessor(self) -> None:
        from xarray.plot.accessor import DataArrayPlotAccessor

        assert DataArray.plot is DataArrayPlotAccessor
        assert isinstance(self.darray.plot, DataArrayPlotAccessor)

    def test_label_from_attrs(self) -> None:
        da = self.darray.copy()
        assert "" == label_from_attrs(da)

        da.name = 0
        assert "0" == label_from_attrs(da)

        da.name = "a"
        da.attrs["units"] = "a_units"
        da.attrs["long_name"] = "a_long_name"
        da.attrs["standard_name"] = "a_standard_name"
        assert "a_long_name [a_units]" == label_from_attrs(da)

        da.attrs.pop("long_name")
        assert "a_standard_name [a_units]" == label_from_attrs(da)
        da.attrs.pop("units")
        assert "a_standard_name" == label_from_attrs(da)

        da.attrs["units"] = "a_units"
        da.attrs.pop("standard_name")
        assert "a [a_units]" == label_from_attrs(da)

        da.attrs.pop("units")
        assert "a" == label_from_attrs(da)

        # Latex strings can be longer without needing a new line:
        long_latex_name = r"$Ra_s = \mathrm{mean}(\epsilon_k) / \mu M^2_\infty$"
        da.attrs = dict(long_name=long_latex_name)
        assert label_from_attrs(da) == long_latex_name

    def test1d(self) -> None:
        self.darray[:, 0, 0].plot()  # type: ignore[call-arg]

        with pytest.raises(ValueError, match=r"x must be one of None, 'dim_0'"):
            self.darray[:, 0, 0].plot(x="dim_1")  # type: ignore[call-arg]

        with pytest.raises(TypeError, match=r"complex128"):
            (self.darray[:, 0, 0] + 1j).plot()  # type: ignore[call-arg]

    def test_1d_bool(self) -> None:
        xr.ones_like(self.darray[:, 0, 0], dtype=bool).plot()  # type: ignore[call-arg]

    def test_1d_x_y_kw(self) -> None:
        z = np.arange(10)
        da = DataArray(np.cos(z), dims=["z"], coords=[z], name="f")

        xy: list[list[str | None]] = [[None, None], [None, "z"], ["z", None]]

        f, axs = plt.subplots(3, 1, squeeze=False)
        for aa, (x, y) in enumerate(xy):
            da.plot(x=x, y=y, ax=axs.flat[aa])  # type: ignore[call-arg]

        with pytest.raises(ValueError, match=r"Cannot specify both"):
            da.plot(x="z", y="z")  # type: ignore[call-arg]

        error_msg = "must be one of None, 'z'"
        with pytest.raises(ValueError, match=rf"x {error_msg}"):
            da.plot(x="f")  # type: ignore[call-arg]

        with pytest.raises(ValueError, match=rf"y {error_msg}"):
            da.plot(y="f")  # type: ignore[call-arg]

    def test_multiindex_level_as_coord(self) -> None:
        da = xr.DataArray(
            np.arange(5),
            dims="x",
            coords=dict(a=("x", np.arange(5)), b=("x", np.arange(5, 10))),
        )
        da = da.set_index(x=["a", "b"])

        for x in ["a", "b"]:
            h = da.plot(x=x)[0]  # type: ignore[call-arg]
            assert_array_equal(h.get_xdata(), da[x].values)

        for y in ["a", "b"]:
            h = da.plot(y=y)[0]  # type: ignore[call-arg]
            assert_array_equal(h.get_ydata(), da[y].values)

    # Test for bug in GH issue #2725
    def test_infer_line_data(self) -> None:
        current = DataArray(
            name="I",
            data=np.array([5, 8]),
            dims=["t"],
            coords={
                "t": (["t"], np.array([0.1, 0.2])),
                "V": (["t"], np.array([100, 200])),
            },
        )

        # Plot current against voltage
        line = current.plot.line(x="V")[0]
        assert_array_equal(line.get_xdata(), current.coords["V"].values)

        # Plot current against time
        line = current.plot.line()[0]
        assert_array_equal(line.get_xdata(), current.coords["t"].values)

    def test_line_plot_along_1d_coord(self) -> None:
        # Test for bug in GH #3334
        x_coord = xr.DataArray(data=[0.1, 0.2], dims=["x"])
        t_coord = xr.DataArray(data=[10, 20], dims=["t"])

        da = xr.DataArray(
            data=np.array([[0, 1], [5, 9]]),
            dims=["x", "t"],
            coords={"x": x_coord, "time": t_coord},
        )

        line = da.plot(x="time", hue="x")[0]  # type: ignore[call-arg]
        assert_array_equal(line.get_xdata(), da.coords["time"].values)

        line = da.plot(y="time", hue="x")[0]  # type: ignore[call-arg]
        assert_array_equal(line.get_ydata(), da.coords["time"].values)

    def test_line_plot_wrong_hue(self) -> None:
        da = xr.DataArray(
            data=np.array([[0, 1], [5, 9]]),
            dims=["x", "t"],
        )

        with pytest.raises(ValueError, match="hue must be one of"):
            da.plot(x="t", hue="wrong_coord")  # type: ignore[call-arg]

    def test_2d_line(self) -> None:
        with pytest.raises(ValueError, match=r"hue"):
            self.darray[:, :, 0].plot.line()

        self.darray[:, :, 0].plot.line(hue="dim_1")
        self.darray[:, :, 0].plot.line(x="dim_1")
        self.darray[:, :, 0].plot.line(y="dim_1")
        self.darray[:, :, 0].plot.line(x="dim_0", hue="dim_1")
        self.darray[:, :, 0].plot.line(y="dim_0", hue="dim_1")

        with pytest.raises(ValueError, match=r"Cannot"):
            self.darray[:, :, 0].plot.line(x="dim_1", y="dim_0", hue="dim_1")

    def test_2d_line_accepts_legend_kw(self) -> None:
        self.darray[:, :, 0].plot.line(x="dim_0", add_legend=False)
        assert not plt.gca().get_legend()
        plt.cla()
        self.darray[:, :, 0].plot.line(x="dim_0", add_legend=True)
        legend = plt.gca().get_legend()
        assert legend is not None
        # check whether legend title is set
        assert legend.get_title().get_text() == "dim_1"

    def test_2d_line_accepts_x_kw(self) -> None:
        self.darray[:, :, 0].plot.line(x="dim_0")
        assert plt.gca().get_xlabel() == "dim_0"
        plt.cla()
        self.darray[:, :, 0].plot.line(x="dim_1")
        assert plt.gca().get_xlabel() == "dim_1"

    def test_2d_line_accepts_hue_kw(self) -> None:
        self.darray[:, :, 0].plot.line(hue="dim_0")
        legend = plt.gca().get_legend()
        assert legend is not None
        assert legend.get_title().get_text() == "dim_0"
        plt.cla()
        self.darray[:, :, 0].plot.line(hue="dim_1")
        legend = plt.gca().get_legend()
        assert legend is not None
        assert legend.get_title().get_text() == "dim_1"

    def test_2d_coords_line_plot(self) -> None:
        lon, lat = np.meshgrid(np.linspace(-20, 20, 5), np.linspace(0, 30, 4))
        lon += lat / 10
        lat += lon / 10
        da = xr.DataArray(
            np.arange(20).reshape(4, 5),
            dims=["y", "x"],
            coords={"lat": (("y", "x"), lat), "lon": (("y", "x"), lon)},
        )

        with figure_context():
            hdl = da.plot.line(x="lon", hue="x")
            assert len(hdl) == 5

        with figure_context():
            hdl = da.plot.line(x="lon", hue="y")
            assert len(hdl) == 4

        with pytest.raises(ValueError, match="For 2D inputs, hue must be a dimension"):
            da.plot.line(x="lon", hue="lat")

    def test_2d_coord_line_plot_coords_transpose_invariant(self) -> None:
        # checks for bug reported in GH #3933
        x = np.arange(10)
        y = np.arange(20)
        ds = xr.Dataset(coords={"x": x, "y": y})

        for z in [ds.y + ds.x, ds.x + ds.y]:
            ds = ds.assign_coords(z=z)
            ds["v"] = ds.x + ds.y
            ds["v"].plot.line(y="z", hue="x")

    def test_2d_before_squeeze(self) -> None:
        a = DataArray(easy_array((1, 5)))
        a.plot()  # type: ignore[call-arg]

    def test2d_uniform_calls_imshow(self) -> None:
        assert self.imshow_called(self.darray[:, :, 0].plot.imshow)

    @pytest.mark.slow
    def test2d_nonuniform_calls_contourf(self) -> None:
        a = self.darray[:, :, 0]
        a.coords["dim_1"] = [2, 1, 89]
        assert self.contourf_called(a.plot.contourf)

    def test2d_1d_2d_coordinates_contourf(self) -> None:
        sz = (20, 10)
        depth = easy_array(sz)
        a = DataArray(
            easy_array(sz),
            dims=["z", "time"],
            coords={"depth": (["z", "time"], depth), "time": np.linspace(0, 1, sz[1])},
        )

        a.plot.contourf(x="time", y="depth")
        a.plot.contourf(x="depth", y="time")

    def test2d_1d_2d_coordinates_pcolormesh(self) -> None:
        # Test with equal coordinates to catch bug from #5097
        sz = 10
        y2d, x2d = np.meshgrid(np.arange(sz), np.arange(sz))
        a = DataArray(
            easy_array((sz, sz)),
            dims=["x", "y"],
            coords={"x2d": (["x", "y"], x2d), "y2d": (["x", "y"], y2d)},
        )

        for x, y in [
            ("x", "y"),
            ("y", "x"),
            ("x2d", "y"),
            ("y", "x2d"),
            ("x", "y2d"),
            ("y2d", "x"),
            ("x2d", "y2d"),
            ("y2d", "x2d"),
        ]:
            p = a.plot.pcolormesh(x=x, y=y)
            v = p.get_paths()[0].vertices
            assert isinstance(v, np.ndarray)

            # Check all vertices are different, except last vertex which should be the
            # same as the first
            _, unique_counts = np.unique(v[:-1], axis=0, return_counts=True)
            assert np.all(unique_counts == 1)

    def test_str_coordinates_pcolormesh(self) -> None:
        # test for #6775
        x = DataArray(
            [[1, 2, 3], [4, 5, 6]],
            dims=("a", "b"),
            coords={"a": [1, 2], "b": ["a", "b", "c"]},
        )
        x.plot.pcolormesh()
        x.T.plot.pcolormesh()

    def test_contourf_cmap_set(self) -> None:
        a = DataArray(easy_array((4, 4)), dims=["z", "time"])

        cmap_expected = mpl.colormaps["viridis"]

        # use copy to ensure cmap is not changed by contourf()
        # Set vmin and vmax so that _build_discrete_colormap is called with
        # extend='both'. extend is passed to
        # mpl.colors.from_levels_and_colors(), which returns a result with
        # sensible under and over values if extend='both', but not if
        # extend='neither' (but if extend='neither' the under and over values
        # would not be used because the data would all be within the plotted
        # range)
        pl = a.plot.contourf(cmap=copy(cmap_expected), vmin=0.1, vmax=0.9)

        # check the set_bad color
        cmap = pl.cmap
        assert cmap is not None
        assert_array_equal(
            cmap(np.ma.masked_invalid([np.nan]))[0],
            cmap_expected(np.ma.masked_invalid([np.nan]))[0],
        )

        # check the set_under color
        assert cmap(-np.inf) == cmap_expected(-np.inf)

        # check the set_over color
        assert cmap(np.inf) == cmap_expected(np.inf)

    def test_contourf_cmap_set_with_bad_under_over(self) -> None:
        a = DataArray(easy_array((4, 4)), dims=["z", "time"])

        # make a copy here because we want a local cmap that we will modify.
        cmap_expected = copy(mpl.colormaps["viridis"])

        cmap_expected.set_bad("w")
        # check we actually changed the set_bad color
        assert np.all(
            cmap_expected(np.ma.masked_invalid([np.nan]))[0]
            != mpl.colormaps["viridis"](np.ma.masked_invalid([np.nan]))[0]
        )

        cmap_expected.set_under("r")
        # check we actually changed the set_under color
        assert cmap_expected(-np.inf) != mpl.colormaps["viridis"](-np.inf)

        cmap_expected.set_over("g")
        # check we actually changed the set_over color
        assert cmap_expected(np.inf) != mpl.colormaps["viridis"](-np.inf)

        # copy to ensure cmap is not changed by contourf()
        pl = a.plot.contourf(cmap=copy(cmap_expected))
        cmap = pl.cmap
        assert cmap is not None

        # check the set_bad color has been kept
        assert_array_equal(
            cmap(np.ma.masked_invalid([np.nan]))[0],
            cmap_expected(np.ma.masked_invalid([np.nan]))[0],
        )

        # check the set_under color has been kept
        assert cmap(-np.inf) == cmap_expected(-np.inf)

        # check the set_over color has been kept
        assert cmap(np.inf) == cmap_expected(np.inf)

    def test3d(self) -> None:
        self.darray.plot()  # type: ignore[call-arg]

    def test_can_pass_in_axis(self) -> None:
        self.pass_in_axis(self.darray.plot)

    def test__infer_interval_breaks(self) -> None:
        assert_array_equal([-0.5, 0.5, 1.5], _infer_interval_breaks([0, 1]))
        assert_array_equal(
            [-0.5, 0.5, 5.0, 9.5, 10.5], _infer_interval_breaks([0, 1, 9, 10])
        )
        assert_array_equal(
            pd.date_range("20000101", periods=4) - np.timedelta64(12, "h"),  # type: ignore[operator]
            _infer_interval_breaks(pd.date_range("20000101", periods=3)),
        )

        # make a bounded 2D array that we will center and re-infer
        xref, yref = np.meshgrid(np.arange(6), np.arange(5))
        cx = (xref[1:, 1:] + xref[:-1, :-1]) / 2
        cy = (yref[1:, 1:] + yref[:-1, :-1]) / 2
        x = _infer_interval_breaks(cx, axis=1)
        x = _infer_interval_breaks(x, axis=0)
        y = _infer_interval_breaks(cy, axis=1)
        y = _infer_interval_breaks(y, axis=0)
        np.testing.assert_allclose(xref, x)
        np.testing.assert_allclose(yref, y)

        # test that ValueError is raised for non-monotonic 1D inputs
        with pytest.raises(ValueError):
            _infer_interval_breaks(np.array([0, 2, 1]), check_monotonic=True)

    def test__infer_interval_breaks_logscale(self) -> None:
        """
        Check if interval breaks are defined in the logspace if scale="log"
        """
        # Check for 1d arrays
        x = np.logspace(-4, 3, 8)
        expected_interval_breaks = 10 ** np.linspace(-4.5, 3.5, 9)
        np.testing.assert_allclose(
            _infer_interval_breaks(x, scale="log"), expected_interval_breaks
        )

        # Check for 2d arrays
        x = np.logspace(-4, 3, 8)
        y = np.linspace(-5, 5, 11)
        x, y = np.meshgrid(x, y)
        expected_interval_breaks = np.vstack([10 ** np.linspace(-4.5, 3.5, 9)] * 12)
        x = _infer_interval_breaks(x, axis=1, scale="log")
        x = _infer_interval_breaks(x, axis=0, scale="log")
        np.testing.assert_allclose(x, expected_interval_breaks)

    def test__infer_interval_breaks_logscale_invalid_coords(self) -> None:
        """
        Check error is raised when passing non-positive coordinates with logscale
        """
        # Check if error is raised after a zero value in the array
        x = np.linspace(0, 5, 6)
        with pytest.raises(ValueError):
            _infer_interval_breaks(x, scale="log")
        # Check if error is raised after negative values in the array
        x = np.linspace(-5, 5, 11)
        with pytest.raises(ValueError):
            _infer_interval_breaks(x, scale="log")

    def test_geo_data(self) -> None:
        # Regression test for gh2250
        # Realistic coordinates taken from the example dataset
        lat = np.array(
            [
                [16.28, 18.48, 19.58, 19.54, 18.35],
                [28.07, 30.52, 31.73, 31.68, 30.37],
                [39.65, 42.27, 43.56, 43.51, 42.11],
                [50.52, 53.22, 54.55, 54.50, 53.06],
            ]
        )
        lon = np.array(
            [
                [-126.13, -113.69, -100.92, -88.04, -75.29],
                [-129.27, -115.62, -101.54, -87.32, -73.26],
                [-133.10, -118.00, -102.31, -86.42, -70.76],
                [-137.85, -120.99, -103.28, -85.28, -67.62],
            ]
        )
        data = np.hypot(lon, lat)
        da = DataArray(
            data,
            dims=("y", "x"),
            coords={"lon": (("y", "x"), lon), "lat": (("y", "x"), lat)},
        )
        da.plot(x="lon", y="lat")  # type: ignore[call-arg]
        ax = plt.gca()
        assert ax.has_data()
        da.plot(x="lat", y="lon")  # type: ignore[call-arg]
        ax = plt.gca()
        assert ax.has_data()

    def test_datetime_dimension(self) -> None:
        nrow = 3
        ncol = 4
        time = pd.date_range("2000-01-01", periods=nrow)
        a = DataArray(
            easy_array((nrow, ncol)), coords=[("time", time), ("y", range(ncol))]
        )
        a.plot()  # type: ignore[call-arg]
        ax = plt.gca()
        assert ax.has_data()

    def test_date_dimension(self) -> None:
        nrow = 3
        ncol = 4
        start = date(2000, 1, 1)
        time = [start + timedelta(days=i) for i in range(nrow)]
        a = DataArray(
            easy_array((nrow, ncol)), coords=[("time", time), ("y", range(ncol))]
        )
        a.plot()  # type: ignore[call-arg]
        ax = plt.gca()
        assert ax.has_data()

    @pytest.mark.slow
    @pytest.mark.filterwarnings("ignore:tight_layout cannot")
    def test_convenient_facetgrid(self) -> None:
        a = easy_array((10, 15, 4))
        d = DataArray(a, dims=["y", "x", "z"])
        d.coords["z"] = list("abcd")
        g = d.plot(x="x", y="y", col="z", col_wrap=2, cmap="cool")  # type: ignore[call-arg]

        assert_array_equal(g.axs.shape, [2, 2])
        for ax in g.axs.flat:
            assert ax.has_data()

        with pytest.raises(ValueError, match=r"[Ff]acet"):
            d.plot(x="x", y="y", col="z", ax=plt.gca())  # type: ignore[call-arg]

        with pytest.raises(ValueError, match=r"[Ff]acet"):
            d[0].plot(x="x", y="y", col="z", ax=plt.gca())  # type: ignore[call-arg]

    @pytest.mark.slow
    def test_subplot_kws(self) -> None:
        a = easy_array((10, 15, 4))
        d = DataArray(a, dims=["y", "x", "z"])
        d.coords["z"] = list("abcd")
        g = d.plot(  # type: ignore[call-arg]
            x="x",
            y="y",
            col="z",
            col_wrap=2,
            cmap="cool",
            subplot_kws=dict(facecolor="r"),
        )
        for ax in g.axs.flat:
            # mpl V2
            assert ax.get_facecolor()[0:3] == mpl.colors.to_rgb("r")

    @pytest.mark.slow
    def test_plot_size(self) -> None:
        self.darray[:, 0, 0].plot(figsize=(13, 5))  # type: ignore[call-arg]
        assert tuple(plt.gcf().get_size_inches()) == (13, 5)

        self.darray.plot(figsize=(13, 5))  # type: ignore[call-arg]
        assert tuple(plt.gcf().get_size_inches()) == (13, 5)

        self.darray.plot(size=5)  # type: ignore[call-arg]
        assert plt.gcf().get_size_inches()[1] == 5

        self.darray.plot(size=5, aspect=2)  # type: ignore[call-arg]
        assert tuple(plt.gcf().get_size_inches()) == (10, 5)

        with pytest.raises(ValueError, match=r"cannot provide both"):
            self.darray.plot(ax=plt.gca(), figsize=(3, 4))  # type: ignore[call-arg]

        with pytest.raises(ValueError, match=r"cannot provide both"):
            self.darray.plot(size=5, figsize=(3, 4))  # type: ignore[call-arg]

        with pytest.raises(ValueError, match=r"cannot provide both"):
            self.darray.plot(size=5, ax=plt.gca())  # type: ignore[call-arg]

        with pytest.raises(ValueError, match=r"cannot provide `aspect`"):
            self.darray.plot(aspect=1)  # type: ignore[call-arg]

    @pytest.mark.slow
    @pytest.mark.filterwarnings("ignore:tight_layout cannot")
    def test_convenient_facetgrid_4d(self) -> None:
        a = easy_array((10, 15, 2, 3))
        d = DataArray(a, dims=["y", "x", "columns", "rows"])
        g = d.plot(x="x", y="y", col="columns", row="rows")  # type: ignore[call-arg]

        assert_array_equal(g.axs.shape, [3, 2])
        for ax in g.axs.flat:
            assert ax.has_data()

        with pytest.raises(ValueError, match=r"[Ff]acet"):
            d.plot(x="x", y="y", col="columns", ax=plt.gca())  # type: ignore[call-arg]

    def test_coord_with_interval(self) -> None:
        """Test line plot with intervals."""
        bins = [-1, 0, 1, 2]
        self.darray.groupby_bins("dim_0", bins).mean(...).plot()  # type: ignore[call-arg]

    def test_coord_with_interval_x(self) -> None:
        """Test line plot with intervals explicitly on x axis."""
        bins = [-1, 0, 1, 2]
        self.darray.groupby_bins("dim_0", bins).mean(...).plot(x="dim_0_bins")  # type: ignore[call-arg]

    def test_coord_with_interval_y(self) -> None:
        """Test line plot with intervals explicitly on y axis."""
        bins = [-1, 0, 1, 2]
        self.darray.groupby_bins("dim_0", bins).mean(...).plot(y="dim_0_bins")  # type: ignore[call-arg]

    def test_coord_with_interval_xy(self) -> None:
        """Test line plot with intervals on both x and y axes."""
        bins = [-1, 0, 1, 2]
        self.darray.groupby_bins("dim_0", bins).mean(...).dim_0_bins.plot()

    @pytest.mark.parametrize("dim", ("x", "y"))
    def test_labels_with_units_with_interval(self, dim) -> None:
        """Test line plot with intervals and a units attribute."""
        bins = [-1, 0, 1, 2]
        arr = self.darray.groupby_bins("dim_0", bins).mean(...)
        arr.dim_0_bins.attrs["units"] = "m"

        (mappable,) = arr.plot(**{dim: "dim_0_bins"})  # type: ignore[arg-type]
        ax = mappable.figure.gca()
        actual = getattr(ax, f"get_{dim}label")()

        expected = "dim_0_bins_center [m]"
        assert actual == expected

    def test_multiplot_over_length_one_dim(self) -> None:
        a = easy_array((3, 1, 1, 1))
        d = DataArray(a, dims=("x", "col", "row", "hue"))
        d.plot(col="col")  # type: ignore[call-arg]
        d.plot(row="row")  # type: ignore[call-arg]
        d.plot(hue="hue")  # type: ignore[call-arg]


class TestPlot1D(PlotTestCase):
    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        d = [0, 1.1, 0, 2]
        self.darray = DataArray(d, coords={"period": range(len(d))}, dims="period")
        self.darray.period.attrs["units"] = "s"

    def test_xlabel_is_index_name(self) -> None:
        self.darray.plot()  # type: ignore[call-arg]
        assert "period [s]" == plt.gca().get_xlabel()

    def test_no_label_name_on_x_axis(self) -> None:
        self.darray.plot(y="period")  # type: ignore[call-arg]
        assert "" == plt.gca().get_xlabel()

    def test_no_label_name_on_y_axis(self) -> None:
        self.darray.plot()  # type: ignore[call-arg]
        assert "" == plt.gca().get_ylabel()

    def test_ylabel_is_data_name(self) -> None:
        self.darray.name = "temperature"
        self.darray.attrs["units"] = "degrees_Celsius"
        self.darray.plot()  # type: ignore[call-arg]
        assert "temperature [degrees_Celsius]" == plt.gca().get_ylabel()

    def test_xlabel_is_data_name(self) -> None:
        self.darray.name = "temperature"
        self.darray.attrs["units"] = "degrees_Celsius"
        self.darray.plot(y="period")  # type: ignore[call-arg]
        assert "temperature [degrees_Celsius]" == plt.gca().get_xlabel()

    def test_format_string(self) -> None:
        self.darray.plot.line("ro")

    def test_can_pass_in_axis(self) -> None:
        self.pass_in_axis(self.darray.plot.line)

    def test_nonnumeric_index(self) -> None:
        a = DataArray([1, 2, 3], {"letter": ["a", "b", "c"]}, dims="letter")
        a.plot.line()

    def test_primitive_returned(self) -> None:
        p = self.darray.plot.line()
        assert isinstance(p[0], mpl.lines.Line2D)

    @pytest.mark.slow
    def test_plot_nans(self) -> None:
        self.darray[1] = np.nan
        self.darray.plot.line()

    def test_dates_are_concise(self) -> None:
        import matplotlib.dates as mdates

        time = pd.date_range("2000-01-01", "2000-01-10")
        a = DataArray(np.arange(len(time)), [("t", time)])
        a.plot.line()

        ax = plt.gca()

        assert isinstance(ax.xaxis.get_major_locator(), mdates.AutoDateLocator)
        assert isinstance(ax.xaxis.get_major_formatter(), mdates.ConciseDateFormatter)

    def test_xyincrease_false_changes_axes(self) -> None:
        self.darray.plot.line(xincrease=False, yincrease=False)
        xlim = plt.gca().get_xlim()
        ylim = plt.gca().get_ylim()
        diffs = xlim[1] - xlim[0], ylim[1] - ylim[0]
        assert all(x < 0 for x in diffs)

    def test_slice_in_title(self) -> None:
        self.darray.coords["d"] = 10.009
        self.darray.plot.line()
        title = plt.gca().get_title()
        assert "d = 10.01" == title

    def test_slice_in_title_single_item_array(self) -> None:
        """Edge case for data of shape (1, N) or (N, 1)."""
        darray = self.darray.expand_dims({"d": np.array([10.009])})
        darray.plot.line(x="period")
        title = plt.gca().get_title()
        assert "d = [10.009]" == title


class TestPlotStep(PlotTestCase):
    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        self.darray = DataArray(easy_array((2, 3, 4)))

    def test_step(self) -> None:
        hdl = self.darray[0, 0].plot.step()
        assert "steps" in hdl[0].get_drawstyle()

    @pytest.mark.parametrize("where", ["pre", "post", "mid"])
    def test_step_with_where(self, where) -> None:
        hdl = self.darray[0, 0].plot.step(where=where)
        assert hdl[0].get_drawstyle() == f"steps-{where}"

    def test_step_with_hue(self) -> None:
        hdl = self.darray[0].plot.step(hue="dim_2")
        assert hdl[0].get_drawstyle() == "steps-pre"

    @pytest.mark.parametrize("where", ["pre", "post", "mid"])
    def test_step_with_hue_and_where(self, where) -> None:
        hdl = self.darray[0].plot.step(hue="dim_2", where=where)
        assert hdl[0].get_drawstyle() == f"steps-{where}"

    def test_drawstyle_steps(self) -> None:
        hdl = self.darray[0].plot(hue="dim_2", drawstyle="steps")  # type: ignore[call-arg]
        assert hdl[0].get_drawstyle() == "steps"

    @pytest.mark.parametrize("where", ["pre", "post", "mid"])
    def test_drawstyle_steps_with_where(self, where) -> None:
        hdl = self.darray[0].plot(hue="dim_2", drawstyle=f"steps-{where}")  # type: ignore[call-arg]
        assert hdl[0].get_drawstyle() == f"steps-{where}"

    def test_coord_with_interval_step(self) -> None:
        """Test step plot with intervals."""
        bins = [-1, 0, 1, 2]
        self.darray.groupby_bins("dim_0", bins).mean(...).plot.step()
        line = plt.gca().lines[0]
        assert isinstance(line, mpl.lines.Line2D)
        assert len(np.asarray(line.get_xdata())) == ((len(bins) - 1) * 2)

    def test_coord_with_interval_step_x(self) -> None:
        """Test step plot with intervals explicitly on x axis."""
        bins = [-1, 0, 1, 2]
        self.darray.groupby_bins("dim_0", bins).mean(...).plot.step(x="dim_0_bins")
        line = plt.gca().lines[0]
        assert isinstance(line, mpl.lines.Line2D)
        assert len(np.asarray(line.get_xdata())) == ((len(bins) - 1) * 2)

    def test_coord_with_interval_step_y(self) -> None:
        """Test step plot with intervals explicitly on y axis."""
        bins = [-1, 0, 1, 2]
        self.darray.groupby_bins("dim_0", bins).mean(...).plot.step(y="dim_0_bins")
        line = plt.gca().lines[0]
        assert isinstance(line, mpl.lines.Line2D)
        assert len(np.asarray(line.get_xdata())) == ((len(bins) - 1) * 2)

    def test_coord_with_interval_step_x_and_y_raises_valueeerror(self) -> None:
        """Test that step plot with intervals both on x and y axes raises an error."""
        arr = xr.DataArray(
            [pd.Interval(0, 1), pd.Interval(1, 2)],
            coords=[("x", [pd.Interval(0, 1), pd.Interval(1, 2)])],
        )
        with pytest.raises(TypeError, match="intervals against intervals"):
            arr.plot.step()


class TestPlotHistogram(PlotTestCase):
    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        self.darray = DataArray(easy_array((2, 3, 4)))

    def test_3d_array(self) -> None:
        self.darray.plot.hist()  # type: ignore[call-arg]

    def test_xlabel_uses_name(self) -> None:
        self.darray.name = "testpoints"
        self.darray.attrs["units"] = "testunits"
        self.darray.plot.hist()  # type: ignore[call-arg]
        assert "testpoints [testunits]" == plt.gca().get_xlabel()

    def test_title_is_histogram(self) -> None:
        self.darray.coords["d"] = 10
        self.darray.plot.hist()  # type: ignore[call-arg]
        assert "d = 10" == plt.gca().get_title()

    def test_can_pass_in_kwargs(self) -> None:
        nbins = 5
        self.darray.plot.hist(bins=nbins)  # type: ignore[call-arg]
        assert nbins == len(plt.gca().patches)

    def test_can_pass_in_axis(self) -> None:
        self.pass_in_axis(self.darray.plot.hist)

    def test_primitive_returned(self) -> None:
        n, bins, patches = self.darray.plot.hist()  # type: ignore[call-arg]
        assert isinstance(n, np.ndarray)
        assert isinstance(bins, np.ndarray)
        assert isinstance(patches, mpl.container.BarContainer)
        assert isinstance(patches[0], mpl.patches.Rectangle)

    @pytest.mark.slow
    def test_plot_nans(self) -> None:
        self.darray[0, 0, 0] = np.nan
        self.darray.plot.hist()  # type: ignore[call-arg]

    def test_hist_coord_with_interval(self) -> None:
        (
            self.darray.groupby_bins("dim_0", [-1, 0, 1, 2])  # type: ignore[call-arg]
            .mean(...)
            .plot.hist(range=(-1, 2))
        )


@requires_matplotlib
class TestDetermineCmapParams:
    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        self.data = np.linspace(0, 1, num=100)

    def test_robust(self) -> None:
        cmap_params = _determine_cmap_params(self.data, robust=True)
        assert cmap_params["vmin"] == np.percentile(self.data, 2)
        assert cmap_params["vmax"] == np.percentile(self.data, 98)
        assert cmap_params["cmap"] == "viridis"
        assert cmap_params["extend"] == "both"
        assert cmap_params["levels"] is None
        assert cmap_params["norm"] is None

    def test_center(self) -> None:
        cmap_params = _determine_cmap_params(self.data, center=0.5)
        assert cmap_params["vmax"] - 0.5 == 0.5 - cmap_params["vmin"]
        assert cmap_params["cmap"] == "RdBu_r"
        assert cmap_params["extend"] == "neither"
        assert cmap_params["levels"] is None
        assert cmap_params["norm"] is None

    def test_cmap_sequential_option(self) -> None:
        with xr.set_options(cmap_sequential="magma"):
            cmap_params = _determine_cmap_params(self.data)
            assert cmap_params["cmap"] == "magma"

    def test_cmap_sequential_explicit_option(self) -> None:
        with xr.set_options(cmap_sequential=mpl.colormaps["magma"]):
            cmap_params = _determine_cmap_params(self.data)
            assert cmap_params["cmap"] == mpl.colormaps["magma"]

    def test_cmap_divergent_option(self) -> None:
        with xr.set_options(cmap_divergent="magma"):
            cmap_params = _determine_cmap_params(self.data, center=0.5)
            assert cmap_params["cmap"] == "magma"

    def test_nan_inf_are_ignored(self) -> None:
        cmap_params1 = _determine_cmap_params(self.data)
        data = self.data
        data[50:55] = np.nan
        data[56:60] = np.inf
        cmap_params2 = _determine_cmap_params(data)
        assert cmap_params1["vmin"] == cmap_params2["vmin"]
        assert cmap_params1["vmax"] == cmap_params2["vmax"]

    @pytest.mark.slow
    def test_integer_levels(self) -> None:
        data = self.data + 1

        # default is to cover full data range but with no guarantee on Nlevels
        for level in np.arange(2, 10, dtype=int):
            cmap_params = _determine_cmap_params(data, levels=level)
            assert cmap_params["vmin"] is None
            assert cmap_params["vmax"] is None
            assert cmap_params["norm"].vmin == cmap_params["levels"][0]
            assert cmap_params["norm"].vmax == cmap_params["levels"][-1]
            assert cmap_params["extend"] == "neither"

        # with min max we are more strict
        cmap_params = _determine_cmap_params(
            data, levels=5, vmin=0, vmax=5, cmap="Blues"
        )
        assert cmap_params["vmin"] is None
        assert cmap_params["vmax"] is None
        assert cmap_params["norm"].vmin == 0
        assert cmap_params["norm"].vmax == 5
        assert cmap_params["norm"].vmin == cmap_params["levels"][0]
        assert cmap_params["norm"].vmax == cmap_params["levels"][-1]
        assert cmap_params["cmap"].name == "Blues"
        assert cmap_params["extend"] == "neither"
        assert cmap_params["cmap"].N == 4
        assert cmap_params["norm"].N == 5

        cmap_params = _determine_cmap_params(data, levels=5, vmin=0.5, vmax=1.5)
        assert cmap_params["cmap"].name == "viridis"
        assert cmap_params["extend"] == "max"

        cmap_params = _determine_cmap_params(data, levels=5, vmin=1.5)
        assert cmap_params["cmap"].name == "viridis"
        assert cmap_params["extend"] == "min"

        cmap_params = _determine_cmap_params(data, levels=5, vmin=1.3, vmax=1.5)
        assert cmap_params["cmap"].name == "viridis"
        assert cmap_params["extend"] == "both"

    def test_list_levels(self) -> None:
        data = self.data + 1

        orig_levels = [0, 1, 2, 3, 4, 5]
        # vmin and vmax should be ignored if levels are explicitly provided
        cmap_params = _determine_cmap_params(data, levels=orig_levels, vmin=0, vmax=3)
        assert cmap_params["vmin"] is None
        assert cmap_params["vmax"] is None
        assert cmap_params["norm"].vmin == 0
        assert cmap_params["norm"].vmax == 5
        assert cmap_params["cmap"].N == 5
        assert cmap_params["norm"].N == 6

        for wrap_levels in cast(
            list[Callable[[Any], dict[Any, Any]]], [list, np.array, pd.Index, DataArray]
        ):
            cmap_params = _determine_cmap_params(data, levels=wrap_levels(orig_levels))
            assert_array_equal(cmap_params["levels"], orig_levels)

    def test_divergentcontrol(self) -> None:
        neg = self.data - 0.1
        pos = self.data

        # Default with positive data will be a normal cmap
        cmap_params = _determine_cmap_params(pos)
        assert cmap_params["vmin"] == 0
        assert cmap_params["vmax"] == 1
        assert cmap_params["cmap"] == "viridis"

        # Default with negative data will be a divergent cmap
        cmap_params = _determine_cmap_params(neg)
        assert cmap_params["vmin"] == -0.9
        assert cmap_params["vmax"] == 0.9
        assert cmap_params["cmap"] == "RdBu_r"

        # Setting vmin or vmax should prevent this only if center is false
        cmap_params = _determine_cmap_params(neg, vmin=-0.1, center=False)
        assert cmap_params["vmin"] == -0.1
        assert cmap_params["vmax"] == 0.9
        assert cmap_params["cmap"] == "viridis"
        cmap_params = _determine_cmap_params(neg, vmax=0.5, center=False)
        assert cmap_params["vmin"] == -0.1
        assert cmap_params["vmax"] == 0.5
        assert cmap_params["cmap"] == "viridis"

        # Setting center=False too
        cmap_params = _determine_cmap_params(neg, center=False)
        assert cmap_params["vmin"] == -0.1
        assert cmap_params["vmax"] == 0.9
        assert cmap_params["cmap"] == "viridis"

        # However, I should still be able to set center and have a div cmap
        cmap_params = _determine_cmap_params(neg, center=0)
        assert cmap_params["vmin"] == -0.9
        assert cmap_params["vmax"] == 0.9
        assert cmap_params["cmap"] == "RdBu_r"

        # Setting vmin or vmax alone will force symmetric bounds around center
        cmap_params = _determine_cmap_params(neg, vmin=-0.1)
        assert cmap_params["vmin"] == -0.1
        assert cmap_params["vmax"] == 0.1
        assert cmap_params["cmap"] == "RdBu_r"
        cmap_params = _determine_cmap_params(neg, vmax=0.5)
        assert cmap_params["vmin"] == -0.5
        assert cmap_params["vmax"] == 0.5
        assert cmap_params["cmap"] == "RdBu_r"
        cmap_params = _determine_cmap_params(neg, vmax=0.6, center=0.1)
        assert cmap_params["vmin"] == -0.4
        assert cmap_params["vmax"] == 0.6
        assert cmap_params["cmap"] == "RdBu_r"

        # But this is only true if vmin or vmax are negative
        cmap_params = _determine_cmap_params(pos, vmin=-0.1)
        assert cmap_params["vmin"] == -0.1
        assert cmap_params["vmax"] == 0.1
        assert cmap_params["cmap"] == "RdBu_r"
        cmap_params = _determine_cmap_params(pos, vmin=0.1)
        assert cmap_params["vmin"] == 0.1
        assert cmap_params["vmax"] == 1
        assert cmap_params["cmap"] == "viridis"
        cmap_params = _determine_cmap_params(pos, vmax=0.5)
        assert cmap_params["vmin"] == 0
        assert cmap_params["vmax"] == 0.5
        assert cmap_params["cmap"] == "viridis"

        # If both vmin and vmax are provided, output is non-divergent
        cmap_params = _determine_cmap_params(neg, vmin=-0.2, vmax=0.6)
        assert cmap_params["vmin"] == -0.2
        assert cmap_params["vmax"] == 0.6
        assert cmap_params["cmap"] == "viridis"

        # regression test for GH3524
        # infer diverging colormap from divergent levels
        cmap_params = _determine_cmap_params(pos, levels=[-0.1, 0, 1])
        # specifying levels makes cmap a Colormap object
        assert cmap_params["cmap"].name == "RdBu_r"

    def test_norm_sets_vmin_vmax(self) -> None:
        vmin = self.data.min()
        vmax = self.data.max()

        for norm, extend, levels in zip(
            [
                mpl.colors.Normalize(),
                mpl.colors.Normalize(),
                mpl.colors.Normalize(vmin + 0.1, vmax - 0.1),
                mpl.colors.Normalize(None, vmax - 0.1),
                mpl.colors.Normalize(vmin + 0.1, None),
            ],
            ["neither", "neither", "both", "max", "min"],
            [7, None, None, None, None],
            strict=True,
        ):
            test_min = vmin if norm.vmin is None else norm.vmin
            test_max = vmax if norm.vmax is None else norm.vmax

            cmap_params = _determine_cmap_params(self.data, norm=norm, levels=levels)
            assert cmap_params["vmin"] is None
            assert cmap_params["vmax"] is None
            assert cmap_params["norm"].vmin == test_min
            assert cmap_params["norm"].vmax == test_max
            assert cmap_params["extend"] == extend
            assert cmap_params["norm"] == norm


@requires_matplotlib
class TestDiscreteColorMap:
    @pytest.fixture(autouse=True)
    def setUp(self):
        x = np.arange(start=0, stop=10, step=2)
        y = np.arange(start=9, stop=-7, step=-3)
        xy = np.dstack(np.meshgrid(x, y))
        distance = np.linalg.norm(xy, axis=2)
        self.darray = DataArray(distance, list(zip(("y", "x"), (y, x), strict=True)))
        self.data_min = distance.min()
        self.data_max = distance.max()
        yield
        # Remove all matplotlib figures
        plt.close("all")

    @pytest.mark.slow
    def test_recover_from_seaborn_jet_exception(self) -> None:
        pal = _color_palette("jet", 4)
        assert type(pal) is np.ndarray
        assert len(pal) == 4

    @pytest.mark.slow
    def test_build_discrete_cmap(self) -> None:
        for cmap, levels, extend, filled in [
            ("jet", [0, 1], "both", False),
            ("hot", [-4, 4], "max", True),
        ]:
            ncmap, cnorm = _build_discrete_cmap(cmap, levels, extend, filled)
            assert ncmap.N == len(levels) - 1
            assert len(ncmap.colors) == len(levels) - 1
            assert cnorm.N == len(levels)
            assert_array_equal(cnorm.boundaries, levels)
            assert max(levels) == cnorm.vmax
            assert min(levels) == cnorm.vmin
            if filled:
                assert ncmap.colorbar_extend == extend
            else:
                assert ncmap.colorbar_extend == "max"

    @pytest.mark.slow
    def test_discrete_colormap_list_of_levels(self) -> None:
        for extend, levels in [
            ("max", [-1, 2, 4, 8, 10]),
            ("both", [2, 5, 10, 11]),
            ("neither", [0, 5, 10, 15]),
            ("min", [2, 5, 10, 15]),
        ]:
            for kind in ["imshow", "pcolormesh", "contourf", "contour"]:
                primitive = getattr(self.darray.plot, kind)(levels=levels)
                assert_array_equal(levels, primitive.norm.boundaries)
                assert max(levels) == primitive.norm.vmax
                assert min(levels) == primitive.norm.vmin
                if kind != "contour":
                    assert extend == primitive.cmap.colorbar_extend
                else:
                    assert "max" == primitive.cmap.colorbar_extend
                assert len(levels) - 1 == len(primitive.cmap.colors)

    @pytest.mark.slow
    def test_discrete_colormap_int_levels(self) -> None:
        for extend, levels, vmin, vmax, cmap in [
            ("neither", 7, None, None, None),
            ("neither", 7, None, 20, mpl.colormaps["RdBu"]),
            ("both", 7, 4, 8, None),
            ("min", 10, 4, 15, None),
        ]:
            for kind in ["imshow", "pcolormesh", "contourf", "contour"]:
                primitive = getattr(self.darray.plot, kind)(
                    levels=levels, vmin=vmin, vmax=vmax, cmap=cmap
                )
                assert levels >= len(primitive.norm.boundaries) - 1
                if vmax is None:
                    assert primitive.norm.vmax >= self.data_max
                else:
                    assert primitive.norm.vmax >= vmax
                if vmin is None:
                    assert primitive.norm.vmin <= self.data_min
                else:
                    assert primitive.norm.vmin <= vmin
                if kind != "contour":
                    assert extend == primitive.cmap.colorbar_extend
                else:
                    assert "max" == primitive.cmap.colorbar_extend
                assert levels >= len(primitive.cmap.colors)

    def test_discrete_colormap_list_levels_and_vmin_or_vmax(self) -> None:
        levels = [0, 5, 10, 15]
        primitive = self.darray.plot(levels=levels, vmin=-3, vmax=20)  # type: ignore[call-arg]
        assert primitive.norm.vmax == max(levels)
        assert primitive.norm.vmin == min(levels)

    def test_discrete_colormap_provided_boundary_norm(self) -> None:
        norm = mpl.colors.BoundaryNorm([0, 5, 10, 15], 4)
        primitive = self.darray.plot.contourf(norm=norm)
        np.testing.assert_allclose(list(primitive.levels), norm.boundaries)

    def test_discrete_colormap_provided_boundary_norm_matching_cmap_levels(
        self,
    ) -> None:
        norm = mpl.colors.BoundaryNorm([0, 5, 10, 15], 4)
        primitive = self.darray.plot.contourf(norm=norm)
        cbar = primitive.colorbar
        assert cbar is not None
        assert cbar.norm.Ncmap == cbar.norm.N  # type: ignore[attr-defined] # Exists, debatable if public though.


class Common2dMixin:
    """
    Common tests for 2d plotting go here.

    These tests assume that a staticmethod for `self.plotfunc` exists.
    Should have the same name as the method.
    """

    darray: DataArray
    plotfunc: staticmethod
    pass_in_axis: Callable

    # Needs to be overridden in TestSurface for facet grid plots
    subplot_kws: dict[Any, Any] | None = None

    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        da = DataArray(
            easy_array((10, 15), start=-1),
            dims=["y", "x"],
            coords={"y": np.arange(10), "x": np.arange(15)},
        )
        # add 2d coords
        ds = da.to_dataset(name="testvar")
        x, y = np.meshgrid(da.x.values, da.y.values)
        ds["x2d"] = DataArray(x, dims=["y", "x"])
        ds["y2d"] = DataArray(y, dims=["y", "x"])
        ds = ds.set_coords(["x2d", "y2d"])
        # set darray and plot method
        self.darray: DataArray = ds.testvar

        # Add CF-compliant metadata
        self.darray.attrs["long_name"] = "a_long_name"
        self.darray.attrs["units"] = "a_units"
        self.darray.x.attrs["long_name"] = "x_long_name"
        self.darray.x.attrs["units"] = "x_units"
        self.darray.y.attrs["long_name"] = "y_long_name"
        self.darray.y.attrs["units"] = "y_units"

        self.plotmethod = getattr(self.darray.plot, self.plotfunc.__name__)

    def test_label_names(self) -> None:
        self.plotmethod()
        assert "x_long_name [x_units]" == plt.gca().get_xlabel()
        assert "y_long_name [y_units]" == plt.gca().get_ylabel()

    def test_1d_raises_valueerror(self) -> None:
        with pytest.raises(ValueError, match=r"DataArray must be 2d"):
            self.plotfunc(self.darray[0, :])

    def test_bool(self) -> None:
        xr.ones_like(self.darray, dtype=bool).plot()  # type: ignore[call-arg]

    def test_complex_raises_typeerror(self) -> None:
        with pytest.raises(TypeError, match=r"complex128"):
            (self.darray + 1j).plot()  # type: ignore[call-arg]

    def test_3d_raises_valueerror(self) -> None:
        a = DataArray(easy_array((2, 3, 4)))
        if self.plotfunc.__name__ == "imshow":
            pytest.skip()
        with pytest.raises(ValueError, match=r"DataArray must be 2d"):
            self.plotfunc(a)

    def test_nonnumeric_index(self) -> None:
        a = DataArray(easy_array((3, 2)), coords=[["a", "b", "c"], ["d", "e"]])
        if self.plotfunc.__name__ == "surface":
            # ax.plot_surface errors with nonnumerics:
            with pytest.raises(TypeError, match="not supported for the input types"):
                self.plotfunc(a)
        else:
            self.plotfunc(a)

    def test_multiindex_raises_typeerror(self) -> None:
        a = DataArray(
            easy_array((3, 2)),
            dims=("x", "y"),
            coords=dict(x=("x", [0, 1, 2]), a=("y", [0, 1]), b=("y", [2, 3])),
        )
        a = a.set_index(y=("a", "b"))
        with pytest.raises(TypeError, match=r"[Pp]lot"):
            self.plotfunc(a)

    def test_can_pass_in_axis(self) -> None:
        self.pass_in_axis(self.plotmethod)

    def test_xyincrease_defaults(self) -> None:
        # With default settings the axis must be ordered regardless
        # of the coords order.
        self.plotfunc(DataArray(easy_array((3, 2)), coords=[[1, 2, 3], [1, 2]]))
        bounds = plt.gca().get_ylim()
        assert bounds[0] < bounds[1]
        bounds = plt.gca().get_xlim()
        assert bounds[0] < bounds[1]
        # Inverted coords
        self.plotfunc(DataArray(easy_array((3, 2)), coords=[[3, 2, 1], [2, 1]]))
        bounds = plt.gca().get_ylim()
        assert bounds[0] < bounds[1]
        bounds = plt.gca().get_xlim()
        assert bounds[0] < bounds[1]

    def test_xyincrease_false_changes_axes(self) -> None:
        self.plotmethod(xincrease=False, yincrease=False)
        xlim = plt.gca().get_xlim()
        ylim = plt.gca().get_ylim()
        diffs = xlim[0] - 14, xlim[1] - 0, ylim[0] - 9, ylim[1] - 0
        assert all(abs(x) < 1 for x in diffs)

    def test_xyincrease_true_changes_axes(self) -> None:
        self.plotmethod(xincrease=True, yincrease=True)
        xlim = plt.gca().get_xlim()
        ylim = plt.gca().get_ylim()
        diffs = xlim[0] - 0, xlim[1] - 14, ylim[0] - 0, ylim[1] - 9
        assert all(abs(x) < 1 for x in diffs)

    def test_dates_are_concise(self) -> None:
        import matplotlib.dates as mdates

        time = pd.date_range("2000-01-01", "2000-01-10")
        a = DataArray(np.random.randn(2, len(time)), [("xx", [1, 2]), ("t", time)])
        self.plotfunc(a, x="t")

        ax = plt.gca()

        assert isinstance(ax.xaxis.get_major_locator(), mdates.AutoDateLocator)
        assert isinstance(ax.xaxis.get_major_formatter(), mdates.ConciseDateFormatter)

    def test_plot_nans(self) -> None:
        x1 = self.darray[:5]
        x2 = self.darray.copy()
        x2[5:] = np.nan

        clim1 = self.plotfunc(x1).get_clim()
        clim2 = self.plotfunc(x2).get_clim()
        assert clim1 == clim2

    @pytest.mark.filterwarnings("ignore::UserWarning")
    @pytest.mark.filterwarnings("ignore:invalid value encountered")
    def test_can_plot_all_nans(self) -> None:
        # regression test for issue #1780
        self.plotfunc(DataArray(np.full((2, 2), np.nan)))

    @pytest.mark.filterwarnings("ignore: Attempting to set")
    def test_can_plot_axis_size_one(self) -> None:
        if self.plotfunc.__name__ not in ("contour", "contourf"):
            self.plotfunc(DataArray(np.ones((1, 1))))

    def test_disallows_rgb_arg(self) -> None:
        with pytest.raises(ValueError):
            # Always invalid for most plots.  Invalid for imshow with 2D data.
            self.plotfunc(DataArray(np.ones((2, 2))), rgb="not None")

    def test_viridis_cmap(self) -> None:
        cmap_name = self.plotmethod(cmap="viridis").get_cmap().name
        assert "viridis" == cmap_name

    def test_default_cmap(self) -> None:
        cmap_name = self.plotmethod().get_cmap().name
        assert "RdBu_r" == cmap_name

        cmap_name = self.plotfunc(abs(self.darray)).get_cmap().name
        assert "viridis" == cmap_name

    @requires_seaborn
    def test_seaborn_palette_as_cmap(self) -> None:
        cmap_name = self.plotmethod(levels=2, cmap="husl").get_cmap().name
        assert "husl" == cmap_name

    def test_can_change_default_cmap(self) -> None:
        cmap_name = self.plotmethod(cmap="Blues").get_cmap().name
        assert "Blues" == cmap_name

    def test_diverging_color_limits(self) -> None:
        artist = self.plotmethod()
        vmin, vmax = artist.get_clim()
        assert round(abs(-vmin - vmax), 7) == 0

    def test_xy_strings(self) -> None:
        self.plotmethod(x="y", y="x")
        ax = plt.gca()
        assert "y_long_name [y_units]" == ax.get_xlabel()
        assert "x_long_name [x_units]" == ax.get_ylabel()

    def test_positional_coord_string(self) -> None:
        self.plotmethod(y="x")
        ax = plt.gca()
        assert "x_long_name [x_units]" == ax.get_ylabel()
        assert "y_long_name [y_units]" == ax.get_xlabel()

        self.plotmethod(x="x")
        ax = plt.gca()
        assert "x_long_name [x_units]" == ax.get_xlabel()
        assert "y_long_name [y_units]" == ax.get_ylabel()

    def test_bad_x_string_exception(self) -> None:
        with pytest.raises(ValueError, match=r"x and y cannot be equal."):
            self.plotmethod(x="y", y="y")

        error_msg = "must be one of None, 'x', 'x2d', 'y', 'y2d'"
        with pytest.raises(ValueError, match=rf"x {error_msg}"):
            self.plotmethod(x="not_a_real_dim", y="y")
        with pytest.raises(ValueError, match=rf"x {error_msg}"):
            self.plotmethod(x="not_a_real_dim")
        with pytest.raises(ValueError, match=rf"y {error_msg}"):
            self.plotmethod(y="not_a_real_dim")
        self.darray.coords["z"] = 100

    def test_coord_strings(self) -> None:
        # 1d coords (same as dims)
        assert {"x", "y"} == set(self.darray.dims)
        self.plotmethod(y="y", x="x")

    def test_non_linked_coords(self) -> None:
        # plot with coordinate names that are not dimensions
        self.darray.coords["newy"] = self.darray.y + 150
        # Normal case, without transpose
        self.plotfunc(self.darray, x="x", y="newy")
        ax = plt.gca()
        assert "x_long_name [x_units]" == ax.get_xlabel()
        assert "newy" == ax.get_ylabel()
        # ax limits might change between plotfuncs
        # simply ensure that these high coords were passed over
        assert np.min(ax.get_ylim()) > 100.0

    def test_non_linked_coords_transpose(self) -> None:
        # plot with coordinate names that are not dimensions,
        # and with transposed y and x axes
        # This used to raise an error with pcolormesh and contour
        # https://github.com/pydata/xarray/issues/788
        self.darray.coords["newy"] = self.darray.y + 150
        self.plotfunc(self.darray, x="newy", y="x")
        ax = plt.gca()
        assert "newy" == ax.get_xlabel()
        assert "x_long_name [x_units]" == ax.get_ylabel()
        # ax limits might change between plotfuncs
        # simply ensure that these high coords were passed over
        assert np.min(ax.get_xlim()) > 100.0

    def test_multiindex_level_as_coord(self) -> None:
        da = DataArray(
            easy_array((3, 2)),
            dims=("x", "y"),
            coords=dict(x=("x", [0, 1, 2]), a=("y", [0, 1]), b=("y", [2, 3])),
        )
        da = da.set_index(y=["a", "b"])

        for x, y in (("a", "x"), ("b", "x"), ("x", "a"), ("x", "b")):
            self.plotfunc(da, x=x, y=y)

            ax = plt.gca()
            assert x == ax.get_xlabel()
            assert y == ax.get_ylabel()

        with pytest.raises(ValueError, match=r"levels of the same MultiIndex"):
            self.plotfunc(da, x="a", y="b")

        with pytest.raises(ValueError, match=r"y must be one of None, 'a', 'b', 'x'"):
            self.plotfunc(da, x="a", y="y")

    def test_default_title(self) -> None:
        a = DataArray(easy_array((4, 3, 2)), dims=["a", "b", "c"])
        a.coords["c"] = [0, 1]
        a.coords["d"] = "foo"
        self.plotfunc(a.isel(c=1))
        title = plt.gca().get_title()
        assert title in {"c = 1, d = foo", "d = foo, c = 1"}

    def test_colorbar_default_label(self) -> None:
        self.plotmethod(add_colorbar=True)
        assert "a_long_name [a_units]" in text_in_fig()

    def test_no_labels(self) -> None:
        self.darray.name = "testvar"
        self.darray.attrs["units"] = "test_units"
        self.plotmethod(add_labels=False)
        alltxt = text_in_fig()
        for string in [
            "x_long_name [x_units]",
            "y_long_name [y_units]",
            "testvar [test_units]",
        ]:
            assert string not in alltxt

    def test_colorbar_kwargs(self) -> None:
        # replace label
        self.darray.attrs.pop("long_name")
        self.darray.attrs["units"] = "test_units"
        # check default colorbar label
        self.plotmethod(add_colorbar=True)
        alltxt = text_in_fig()
        assert "testvar [test_units]" in alltxt
        self.darray.attrs.pop("units")

        self.darray.name = "testvar"
        self.plotmethod(add_colorbar=True, cbar_kwargs={"label": "MyLabel"})
        alltxt = text_in_fig()
        assert "MyLabel" in alltxt
        assert "testvar" not in alltxt
        # you can use anything accepted by the dict constructor as well
        self.plotmethod(add_colorbar=True, cbar_kwargs=(("label", "MyLabel"),))
        alltxt = text_in_fig()
        assert "MyLabel" in alltxt
        assert "testvar" not in alltxt
        # change cbar ax
        fig, axs = plt.subplots(1, 2, squeeze=False)
        ax = axs[0, 0]
        cax = axs[0, 1]
        self.plotmethod(
            ax=ax, cbar_ax=cax, add_colorbar=True, cbar_kwargs={"label": "MyBar"}
        )
        assert ax.has_data()
        assert cax.has_data()
        alltxt = text_in_fig()
        assert "MyBar" in alltxt
        assert "testvar" not in alltxt
        # note that there are two ways to achieve this
        fig, axs = plt.subplots(1, 2, squeeze=False)
        ax = axs[0, 0]
        cax = axs[0, 1]
        self.plotmethod(
            ax=ax, add_colorbar=True, cbar_kwargs={"label": "MyBar", "cax": cax}
        )
        assert ax.has_data()
        assert cax.has_data()
        alltxt = text_in_fig()
        assert "MyBar" in alltxt
        assert "testvar" not in alltxt
        # see that no colorbar is respected
        self.plotmethod(add_colorbar=False)
        assert "testvar" not in text_in_fig()
        # check that error is raised
        pytest.raises(
            ValueError,
            self.plotmethod,
            add_colorbar=False,
            cbar_kwargs={"label": "label"},
        )

    def test_verbose_facetgrid(self) -> None:
        a = easy_array((10, 15, 3))
        d = DataArray(a, dims=["y", "x", "z"])
        g = xplt.FacetGrid(d, col="z", subplot_kws=self.subplot_kws)
        g.map_dataarray(self.plotfunc, "x", "y")
        for ax in g.axs.flat:
            assert ax.has_data()

    def test_2d_function_and_method_signature_same(self) -> None:
        func_sig = inspect.signature(self.plotfunc)
        method_sig = inspect.signature(self.plotmethod)
        for argname, param in method_sig.parameters.items():
            assert func_sig.parameters[argname] == param

    @pytest.mark.filterwarnings("ignore:tight_layout cannot")
    def test_convenient_facetgrid(self) -> None:
        a = easy_array((10, 15, 4))
        d = DataArray(a, dims=["y", "x", "z"])
        g = self.plotfunc(d, x="x", y="y", col="z", col_wrap=2)

        assert_array_equal(g.axs.shape, [2, 2])
        for (y, x), ax in np.ndenumerate(g.axs):
            assert ax.has_data()
            if x == 0:
                assert "y" == ax.get_ylabel()
            else:
                assert "" == ax.get_ylabel()
            if y == 1:
                assert "x" == ax.get_xlabel()
            else:
                assert "" == ax.get_xlabel()

        # Inferring labels
        g = self.plotfunc(d, col="z", col_wrap=2)
        assert_array_equal(g.axs.shape, [2, 2])
        for (y, x), ax in np.ndenumerate(g.axs):
            assert ax.has_data()
            if x == 0:
                assert "y" == ax.get_ylabel()
            else:
                assert "" == ax.get_ylabel()
            if y == 1:
                assert "x" == ax.get_xlabel()
            else:
                assert "" == ax.get_xlabel()

    @pytest.mark.filterwarnings("ignore:tight_layout cannot")
    def test_convenient_facetgrid_4d(self) -> None:
        a = easy_array((10, 15, 2, 3))
        d = DataArray(a, dims=["y", "x", "columns", "rows"])
        g = self.plotfunc(d, x="x", y="y", col="columns", row="rows")

        assert_array_equal(g.axs.shape, [3, 2])
        for ax in g.axs.flat:
            assert ax.has_data()

    @pytest.mark.filterwarnings("ignore:This figure includes")
    def test_facetgrid_map_only_appends_mappables(self) -> None:
        a = easy_array((10, 15, 2, 3))
        d = DataArray(a, dims=["y", "x", "columns", "rows"])
        g = self.plotfunc(d, x="x", y="y", col="columns", row="rows")

        expected = g._mappables

        g.map(lambda: plt.plot(1, 1))
        actual = g._mappables

        assert expected == actual

    def test_facetgrid_cmap(self) -> None:
        # Regression test for GH592
        data = np.random.random(size=(20, 25, 12)) + np.linspace(-3, 3, 12)
        d = DataArray(data, dims=["x", "y", "time"])
        fg = d.plot.pcolormesh(col="time")
        # check that all color limits are the same
        assert len({m.get_clim() for m in fg._mappables}) == 1
        # check that all colormaps are the same
        assert len({m.get_cmap().name for m in fg._mappables}) == 1

    def test_facetgrid_cbar_kwargs(self) -> None:
        a = easy_array((10, 15, 2, 3))
        d = DataArray(a, dims=["y", "x", "columns", "rows"])
        g = self.plotfunc(
            d,
            x="x",
            y="y",
            col="columns",
            row="rows",
            cbar_kwargs={"label": "test_label"},
        )

        # catch contour case
        if g.cbar is not None:
            assert get_colorbar_label(g.cbar) == "test_label"

    def test_facetgrid_no_cbar_ax(self) -> None:
        a = easy_array((10, 15, 2, 3))
        d = DataArray(a, dims=["y", "x", "columns", "rows"])
        with pytest.raises(ValueError):
            self.plotfunc(d, x="x", y="y", col="columns", row="rows", cbar_ax=1)

    def test_cmap_and_color_both(self) -> None:
        with pytest.raises(ValueError):
            self.plotmethod(colors="k", cmap="RdBu")

    def test_2d_coord_with_interval(self) -> None:
        for dim in self.darray.dims:
            gp = self.darray.groupby_bins(dim, range(15), restore_coord_dims=True).mean(
                [dim]
            )
            for kind in ["imshow", "pcolormesh", "contourf", "contour"]:
                getattr(gp.plot, kind)()

    def test_colormap_error_norm_and_vmin_vmax(self) -> None:
        norm = mpl.colors.LogNorm(0.1, 1e1)

        with pytest.raises(ValueError):
            self.darray.plot(norm=norm, vmin=2)  # type: ignore[call-arg]

        with pytest.raises(ValueError):
            self.darray.plot(norm=norm, vmax=2)  # type: ignore[call-arg]


@pytest.mark.slow
class TestContourf(Common2dMixin, PlotTestCase):
    plotfunc = staticmethod(xplt.contourf)

    @pytest.mark.slow
    def test_contourf_called(self) -> None:
        # Having both statements ensures the test works properly
        assert not self.contourf_called(self.darray.plot.imshow)
        assert self.contourf_called(self.darray.plot.contourf)

    def test_primitive_artist_returned(self) -> None:
        artist = self.plotmethod()
        assert isinstance(artist, mpl.contour.QuadContourSet)

    @pytest.mark.slow
    def test_extend(self) -> None:
        artist = self.plotmethod()
        assert artist.extend == "neither"

        self.darray[0, 0] = -100
        self.darray[-1, -1] = 100
        artist = self.plotmethod(robust=True)
        assert artist.extend == "both"

        self.darray[0, 0] = 0
        self.darray[-1, -1] = 0
        artist = self.plotmethod(vmin=-0, vmax=10)
        assert artist.extend == "min"

        artist = self.plotmethod(vmin=-10, vmax=0)
        assert artist.extend == "max"

    @pytest.mark.slow
    def test_2d_coord_names(self) -> None:
        self.plotmethod(x="x2d", y="y2d")
        # make sure labels came out ok
        ax = plt.gca()
        assert "x2d" == ax.get_xlabel()
        assert "y2d" == ax.get_ylabel()

    @pytest.mark.slow
    def test_levels(self) -> None:
        artist = self.plotmethod(levels=[-0.5, -0.4, 0.1])
        assert artist.extend == "both"

        artist = self.plotmethod(levels=3)
        assert artist.extend == "neither"

    def test_colormap_norm(self) -> None:
        # Using a norm should plot a nice colorbar and look consistent with pcolormesh.
        norm = mpl.colors.LogNorm(0.1, 1e1)

        with pytest.warns(UserWarning):
            artist = self.plotmethod(norm=norm, add_colorbar=True)

        actual = artist.colorbar.locator()
        expected = np.array([0.01, 0.1, 1.0, 10.0])

        np.testing.assert_allclose(actual, expected)


@pytest.mark.slow
class TestContour(Common2dMixin, PlotTestCase):
    plotfunc = staticmethod(xplt.contour)

    # matplotlib cmap.colors gives an rgbA ndarray
    # when seaborn is used, instead we get an rgb tuple
    @staticmethod
    def _color_as_tuple(c: Any) -> tuple[Any, Any, Any]:
        return c[0], c[1], c[2]

    def test_colors(self) -> None:
        # with single color, we don't want rgb array
        artist = self.plotmethod(colors="k")
        assert artist.cmap.colors[0] == "k"

        # 2 colors, will repeat every other tick:
        artist = self.plotmethod(colors=["k", "b"])
        assert artist.cmap.colors[:2] == ["k", "b"]

        # 4 colors, will repeat every 4th tick:
        artist = self.darray.plot.contour(
            levels=[-0.5, 0.0, 0.5, 1.0], colors=["k", "r", "w", "b"]
        )
        assert artist.cmap.colors[:5] == ["k", "r", "w", "b"]  # type: ignore[attr-defined,unused-ignore]

        # the last color is now under "over"
        assert self._color_as_tuple(artist.cmap.get_over()) == (0.0, 0.0, 1.0)

    def test_colors_np_levels(self) -> None:
        # https://github.com/pydata/xarray/issues/3284
        levels = np.array([-0.5, 0.0, 0.5, 1.0])
        artist = self.darray.plot.contour(levels=levels, colors=["k", "r", "w", "b"])
        cmap = artist.cmap
        assert isinstance(cmap, mpl.colors.ListedColormap)

        assert artist.cmap.colors[:5] == ["k", "r", "w", "b"]  # type: ignore[attr-defined,unused-ignore]

        # the last color is now under "over"
        assert self._color_as_tuple(cmap.get_over()) == (0.0, 0.0, 1.0)

    def test_cmap_and_color_both(self) -> None:
        with pytest.raises(ValueError):
            self.plotmethod(colors="k", cmap="RdBu")

    def list_of_colors_in_cmap_raises_error(self) -> None:
        with pytest.raises(ValueError, match=r"list of colors"):
            self.plotmethod(cmap=["k", "b"])

    @pytest.mark.slow
    def test_2d_coord_names(self) -> None:
        self.plotmethod(x="x2d", y="y2d")
        # make sure labels came out ok
        ax = plt.gca()
        assert "x2d" == ax.get_xlabel()
        assert "y2d" == ax.get_ylabel()

    def test_single_level(self) -> None:
        # this used to raise an error, but not anymore since
        # add_colorbar defaults to false
        self.plotmethod(levels=[0.1])
        self.plotmethod(levels=1)

    def test_colormap_norm(self) -> None:
        # Using a norm should plot a nice colorbar and look consistent with pcolormesh.
        norm = mpl.colors.LogNorm(0.1, 1e1)

        with pytest.warns(UserWarning):
            artist = self.plotmethod(norm=norm, add_colorbar=True)

        actual = artist.colorbar.locator()
        expected = np.array([0.01, 0.1, 1.0, 10.0])

        np.testing.assert_allclose(actual, expected)


class TestPcolormesh(Common2dMixin, PlotTestCase):
    plotfunc = staticmethod(xplt.pcolormesh)

    def test_primitive_artist_returned(self) -> None:
        artist = self.plotmethod()
        assert isinstance(artist, mpl.collections.QuadMesh)

    def test_everything_plotted(self) -> None:
        artist = self.plotmethod()
        assert artist.get_array().size == self.darray.size

    @pytest.mark.slow
    def test_2d_coord_names(self) -> None:
        self.plotmethod(x="x2d", y="y2d")
        # make sure labels came out ok
        ax = plt.gca()
        assert "x2d" == ax.get_xlabel()
        assert "y2d" == ax.get_ylabel()

    def test_dont_infer_interval_breaks_for_cartopy(self) -> None:
        # Regression for GH 781
        ax = plt.gca()
        # Simulate a Cartopy Axis
        ax.projection = True  # type: ignore[attr-defined]
        artist = self.plotmethod(x="x2d", y="y2d", ax=ax)
        assert isinstance(artist, mpl.collections.QuadMesh)
        # Let cartopy handle the axis limits and artist size
        arr = artist.get_array()
        assert arr is not None
        assert arr.size <= self.darray.size


class TestPcolormeshLogscale(PlotTestCase):
    """
    Test pcolormesh axes when x and y are in logscale
    """

    plotfunc = staticmethod(xplt.pcolormesh)

    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        self.boundaries = (-1, 9, -4, 3)
        shape = (8, 11)
        x = np.logspace(self.boundaries[0], self.boundaries[1], shape[1])
        y = np.logspace(self.boundaries[2], self.boundaries[3], shape[0])
        da = DataArray(
            easy_array(shape, start=-1),
            dims=["y", "x"],
            coords={"y": y, "x": x},
            name="testvar",
        )
        self.darray = da

    def test_interval_breaks_logspace(self) -> None:
        """
        Check if the outer vertices of the pcolormesh are the expected values

        Checks bugfix for #5333
        """
        artist = self.darray.plot.pcolormesh(xscale="log", yscale="log")

        # Grab the coordinates of the vertices of the Patches
        x_vertices = [p.vertices[:, 0] for p in artist.properties()["paths"]]
        y_vertices = [p.vertices[:, 1] for p in artist.properties()["paths"]]

        # Get the maximum and minimum values for each set of vertices
        xmin, xmax = np.min(x_vertices), np.max(x_vertices)
        ymin, ymax = np.min(y_vertices), np.max(y_vertices)

        # Check if they are equal to 10 to the power of the outer value of its
        # corresponding axis plus or minus the interval in the logspace
        log_interval = 0.5
        np.testing.assert_allclose(xmin, 10 ** (self.boundaries[0] - log_interval))
        np.testing.assert_allclose(xmax, 10 ** (self.boundaries[1] + log_interval))
        np.testing.assert_allclose(ymin, 10 ** (self.boundaries[2] - log_interval))
        np.testing.assert_allclose(ymax, 10 ** (self.boundaries[3] + log_interval))


@pytest.mark.slow
class TestImshow(Common2dMixin, PlotTestCase):
    plotfunc = staticmethod(xplt.imshow)

    @pytest.mark.xfail(
        reason=(
            "Failing inside matplotlib. Should probably be fixed upstream because "
            "other plot functions can handle it. "
            "Remove this test when it works, already in Common2dMixin"
        )
    )
    def test_dates_are_concise(self) -> None:
        import matplotlib.dates as mdates

        time = pd.date_range("2000-01-01", "2000-01-10")
        a = DataArray(np.random.randn(2, len(time)), [("xx", [1, 2]), ("t", time)])
        self.plotfunc(a, x="t")

        ax = plt.gca()

        assert isinstance(ax.xaxis.get_major_locator(), mdates.AutoDateLocator)
        assert isinstance(ax.xaxis.get_major_formatter(), mdates.ConciseDateFormatter)

    @pytest.mark.slow
    def test_imshow_called(self) -> None:
        # Having both statements ensures the test works properly
        assert not self.imshow_called(self.darray.plot.contourf)
        assert self.imshow_called(self.darray.plot.imshow)

    def test_xy_pixel_centered(self) -> None:
        self.darray.plot.imshow(yincrease=False)
        assert np.allclose([-0.5, 14.5], plt.gca().get_xlim())
        assert np.allclose([9.5, -0.5], plt.gca().get_ylim())

    def test_default_aspect_is_auto(self) -> None:
        self.darray.plot.imshow()
        assert "auto" == plt.gca().get_aspect()

    @pytest.mark.slow
    def test_cannot_change_mpl_aspect(self) -> None:
        with pytest.raises(ValueError, match=r"not available in xarray"):
            self.darray.plot.imshow(aspect="equal")

        # with numbers we fall back to fig control
        self.darray.plot.imshow(size=5, aspect=2)
        assert "auto" == plt.gca().get_aspect()
        assert tuple(plt.gcf().get_size_inches()) == (10, 5)

    @pytest.mark.slow
    def test_primitive_artist_returned(self) -> None:
        artist = self.plotmethod()
        assert isinstance(artist, mpl.image.AxesImage)

    @pytest.mark.slow
    @requires_seaborn
    def test_seaborn_palette_needs_levels(self) -> None:
        with pytest.raises(ValueError):
            self.plotmethod(cmap="husl")

    def test_2d_coord_names(self) -> None:
        with pytest.raises(ValueError, match=r"requires 1D coordinates"):
            self.plotmethod(x="x2d", y="y2d")

    def test_plot_rgb_image(self) -> None:
        DataArray(
            easy_array((10, 15, 3), start=0), dims=["y", "x", "band"]
        ).plot.imshow()
        assert 0 == len(find_possible_colorbars())

    def test_plot_rgb_image_explicit(self) -> None:
        DataArray(
            easy_array((10, 15, 3), start=0), dims=["y", "x", "band"]
        ).plot.imshow(y="y", x="x", rgb="band")
        assert 0 == len(find_possible_colorbars())

    def test_plot_rgb_faceted(self) -> None:
        DataArray(
            easy_array((2, 2, 10, 15, 3), start=0), dims=["a", "b", "y", "x", "band"]
        ).plot.imshow(row="a", col="b")
        assert 0 == len(find_possible_colorbars())

    def test_plot_rgba_image_transposed(self) -> None:
        # We can handle the color axis being in any position
        DataArray(
            easy_array((4, 10, 15), start=0), dims=["band", "y", "x"]
        ).plot.imshow()

    def test_warns_ambiguous_dim(self) -> None:
        arr = DataArray(easy_array((3, 3, 3)), dims=["y", "x", "band"])
        with pytest.warns(UserWarning):
            arr.plot.imshow()
        # but doesn't warn if dimensions specified
        arr.plot.imshow(rgb="band")
        arr.plot.imshow(x="x", y="y")

    def test_rgb_errors_too_many_dims(self) -> None:
        arr = DataArray(easy_array((3, 3, 3, 3)), dims=["y", "x", "z", "band"])
        with pytest.raises(ValueError):
            arr.plot.imshow(rgb="band")

    def test_rgb_errors_bad_dim_sizes(self) -> None:
        arr = DataArray(easy_array((5, 5, 5)), dims=["y", "x", "band"])
        with pytest.raises(ValueError):
            arr.plot.imshow(rgb="band")

    @pytest.mark.parametrize(
        ["vmin", "vmax", "robust"],
        [
            (-1, None, False),
            (None, 2, False),
            (-1, 1, False),
            (0, 0, False),
            (0, None, True),
            (None, -1, True),
        ],
    )
    def test_normalize_rgb_imshow(
        self, vmin: float | None, vmax: float | None, robust: bool
    ) -> None:
        da = DataArray(easy_array((5, 5, 3), start=-0.6, stop=1.4))
        arr = da.plot.imshow(vmin=vmin, vmax=vmax, robust=robust).get_array()
        assert arr is not None
        assert 0 <= arr.min() <= arr.max() <= 1

    def test_normalize_rgb_one_arg_error(self) -> None:
        da = DataArray(easy_array((5, 5, 3), start=-0.6, stop=1.4))
        # If passed one bound that implies all out of range, error:
        for vmin, vmax in ((None, -1), (2, None)):
            with pytest.raises(ValueError):
                da.plot.imshow(vmin=vmin, vmax=vmax)
        # If passed two that's just moving the range, *not* an error:
        for vmin2, vmax2 in ((-1.2, -1), (2, 2.1)):
            da.plot.imshow(vmin=vmin2, vmax=vmax2)

    @pytest.mark.parametrize("dtype", [np.uint8, np.int8, np.int16])
    def test_imshow_rgb_values_in_valid_range(self, dtype) -> None:
        da = DataArray(np.arange(75, dtype=dtype).reshape((5, 5, 3)))
        _, ax = plt.subplots()
        out = da.plot.imshow(ax=ax).get_array()
        assert out is not None
        actual_dtype = out.dtype
        assert actual_dtype is not None
        assert actual_dtype == np.uint8
        assert (out[..., :3] == da.values).all()  # Compare without added alpha
        assert (out[..., -1] == 255).all()  # Compare alpha

    @pytest.mark.filterwarnings("ignore:Several dimensions of this array")
    def test_regression_rgb_imshow_dim_size_one(self) -> None:
        # Regression: https://github.com/pydata/xarray/issues/1966
        da = DataArray(easy_array((1, 3, 3), start=0.0, stop=1.0))
        da.plot.imshow()

    def test_origin_overrides_xyincrease(self) -> None:
        da = DataArray(easy_array((3, 2)), coords=[[-2, 0, 2], [-1, 1]])
        with figure_context():
            da.plot.imshow(origin="upper")
            assert plt.xlim()[0] < 0
            assert plt.ylim()[1] < 0

        with figure_context():
            da.plot.imshow(origin="lower")
            assert plt.xlim()[0] < 0
            assert plt.ylim()[0] < 0


class TestSurface(Common2dMixin, PlotTestCase):
    plotfunc = staticmethod(xplt.surface)
    subplot_kws = {"projection": "3d"}

    @pytest.mark.xfail(
        reason=(
            "Failing inside matplotlib. Should probably be fixed upstream because "
            "other plot functions can handle it. "
            "Remove this test when it works, already in Common2dMixin"
        )
    )
    def test_dates_are_concise(self) -> None:
        import matplotlib.dates as mdates

        time = pd.date_range("2000-01-01", "2000-01-10")
        a = DataArray(np.random.randn(2, len(time)), [("xx", [1, 2]), ("t", time)])
        self.plotfunc(a, x="t")

        ax = plt.gca()

        assert isinstance(ax.xaxis.get_major_locator(), mdates.AutoDateLocator)
        assert isinstance(ax.xaxis.get_major_formatter(), mdates.ConciseDateFormatter)

    def test_primitive_artist_returned(self) -> None:
        artist = self.plotmethod()
        assert isinstance(artist, mpl_toolkits.mplot3d.art3d.Poly3DCollection)

    @pytest.mark.slow
    def test_2d_coord_names(self) -> None:
        self.plotmethod(x="x2d", y="y2d")
        # make sure labels came out ok
        ax = plt.gca()
        assert isinstance(ax, mpl_toolkits.mplot3d.axes3d.Axes3D)
        assert "x2d" == ax.get_xlabel()
        assert "y2d" == ax.get_ylabel()
        assert f"{self.darray.long_name} [{self.darray.units}]" == ax.get_zlabel()

    def test_xyincrease_false_changes_axes(self) -> None:
        # Does not make sense for surface plots
        pytest.skip("does not make sense for surface plots")

    def test_xyincrease_true_changes_axes(self) -> None:
        # Does not make sense for surface plots
        pytest.skip("does not make sense for surface plots")

    def test_can_pass_in_axis(self) -> None:
        self.pass_in_axis(self.plotmethod, subplot_kw={"projection": "3d"})

    def test_default_cmap(self) -> None:
        # Does not make sense for surface plots with default arguments
        pytest.skip("does not make sense for surface plots")

    def test_diverging_color_limits(self) -> None:
        # Does not make sense for surface plots with default arguments
        pytest.skip("does not make sense for surface plots")

    def test_colorbar_kwargs(self) -> None:
        # Does not make sense for surface plots with default arguments
        pytest.skip("does not make sense for surface plots")

    def test_cmap_and_color_both(self) -> None:
        # Does not make sense for surface plots with default arguments
        pytest.skip("does not make sense for surface plots")

    def test_seaborn_palette_as_cmap(self) -> None:
        # seaborn does not work with mpl_toolkits.mplot3d
        with pytest.raises(ValueError):
            super().test_seaborn_palette_as_cmap()

    # Need to modify this test for surface(), because all subplots should have labels,
    # not just left and bottom
    @pytest.mark.filterwarnings("ignore:tight_layout cannot")
    def test_convenient_facetgrid(self) -> None:
        a = easy_array((10, 15, 4))
        d = DataArray(a, dims=["y", "x", "z"])
        g = self.plotfunc(d, x="x", y="y", col="z", col_wrap=2)  # type: ignore[arg-type] # https://github.com/python/mypy/issues/15015

        assert_array_equal(g.axs.shape, [2, 2])
        for (_y, _x), ax in np.ndenumerate(g.axs):
            assert ax.has_data()
            assert "y" == ax.get_ylabel()
            assert "x" == ax.get_xlabel()

        # Inferring labels
        g = self.plotfunc(d, col="z", col_wrap=2)  # type: ignore[arg-type] # https://github.com/python/mypy/issues/15015
        assert_array_equal(g.axs.shape, [2, 2])
        for (_y, _x), ax in np.ndenumerate(g.axs):
            assert ax.has_data()
            assert "y" == ax.get_ylabel()
            assert "x" == ax.get_xlabel()

    def test_viridis_cmap(self) -> None:
        return super().test_viridis_cmap()

    def test_can_change_default_cmap(self) -> None:
        return super().test_can_change_default_cmap()

    def test_colorbar_default_label(self) -> None:
        return super().test_colorbar_default_label()

    def test_facetgrid_map_only_appends_mappables(self) -> None:
        return super().test_facetgrid_map_only_appends_mappables()


class TestFacetGrid(PlotTestCase):
    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        d = easy_array((10, 15, 3))
        self.darray = DataArray(d, dims=["y", "x", "z"], coords={"z": ["a", "b", "c"]})
        self.g = xplt.FacetGrid(self.darray, col="z")

    @pytest.mark.slow
    def test_no_args(self) -> None:
        self.g.map_dataarray(xplt.contourf, "x", "y")

        # Don't want colorbar labeled with 'None'
        alltxt = text_in_fig()
        assert "None" not in alltxt

        for ax in self.g.axs.flat:
            assert ax.has_data()

    @pytest.mark.slow
    def test_names_appear_somewhere(self) -> None:
        self.darray.name = "testvar"
        self.g.map_dataarray(xplt.contourf, "x", "y")
        for k, ax in zip("abc", self.g.axs.flat, strict=True):
            assert f"z = {k}" == ax.get_title()

        alltxt = text_in_fig()
        assert self.darray.name in alltxt
        for label in ["x", "y"]:
            assert label in alltxt

    @pytest.mark.slow
    def test_text_not_super_long(self) -> None:
        self.darray.coords["z"] = [100 * letter for letter in "abc"]
        g = xplt.FacetGrid(self.darray, col="z")
        g.map_dataarray(xplt.contour, "x", "y")
        alltxt = text_in_fig()
        maxlen = max(len(txt) for txt in alltxt)
        assert maxlen < 50

        t0 = g.axs[0, 0].get_title()
        assert t0.endswith("...")

    @pytest.mark.slow
    def test_colorbar(self) -> None:
        vmin = self.darray.values.min()
        vmax = self.darray.values.max()
        expected = np.array((vmin, vmax))

        self.g.map_dataarray(xplt.imshow, "x", "y")

        for image in plt.gcf().findobj(mpl.image.AxesImage):
            assert isinstance(image, mpl.image.AxesImage)
            clim = np.array(image.get_clim())
            assert np.allclose(expected, clim)

        assert 1 == len(find_possible_colorbars())

    def test_colorbar_scatter(self) -> None:
        ds = Dataset({"a": (("x", "y"), np.arange(4).reshape(2, 2))})
        fg: xplt.FacetGrid = ds.plot.scatter(x="a", y="a", row="x", hue="a")
        cbar = fg.cbar
        assert cbar is not None
        assert hasattr(cbar, "vmin")
        assert cbar.vmin == 0
        assert hasattr(cbar, "vmax")
        assert cbar.vmax == 3

    @pytest.mark.slow
    def test_empty_cell(self) -> None:
        g = xplt.FacetGrid(self.darray, col="z", col_wrap=2)
        g.map_dataarray(xplt.imshow, "x", "y")

        bottomright = g.axs[-1, -1]
        assert not bottomright.has_data()
        assert not bottomright.get_visible()

    @pytest.mark.slow
    def test_norow_nocol_error(self) -> None:
        with pytest.raises(ValueError, match=r"[Rr]ow"):
            xplt.FacetGrid(self.darray)

    @pytest.mark.slow
    def test_groups(self) -> None:
        self.g.map_dataarray(xplt.imshow, "x", "y")
        upperleft_dict = self.g.name_dicts[0, 0]
        upperleft_array = self.darray.loc[upperleft_dict]
        z0 = self.darray.isel(z=0)

        assert_equal(upperleft_array, z0)

    @pytest.mark.slow
    def test_float_index(self) -> None:
        self.darray.coords["z"] = [0.1, 0.2, 0.4]
        g = xplt.FacetGrid(self.darray, col="z")
        g.map_dataarray(xplt.imshow, "x", "y")

    @pytest.mark.slow
    def test_nonunique_index_error(self) -> None:
        self.darray.coords["z"] = [0.1, 0.2, 0.2]
        with pytest.raises(ValueError, match=r"[Uu]nique"):
            xplt.FacetGrid(self.darray, col="z")

    @pytest.mark.slow
    def test_robust(self) -> None:
        z = np.zeros((20, 20, 2))
        darray = DataArray(z, dims=["y", "x", "z"])
        darray[:, :, 1] = 1
        darray[2, 0, 0] = -1000
        darray[3, 0, 0] = 1000
        g = xplt.FacetGrid(darray, col="z")
        g.map_dataarray(xplt.imshow, "x", "y", robust=True)

        # Color limits should be 0, 1
        # The largest number displayed in the figure should be less than 21
        numbers = set()
        alltxt = text_in_fig()
        for txt in alltxt:
            with contextlib.suppress(ValueError):
                numbers.add(float(txt))
        largest = max(abs(x) for x in numbers)
        assert largest < 21

    @pytest.mark.slow
    def test_can_set_vmin_vmax(self) -> None:
        vmin, vmax = 50.0, 1000.0
        expected = np.array((vmin, vmax))
        self.g.map_dataarray(xplt.imshow, "x", "y", vmin=vmin, vmax=vmax)

        for image in plt.gcf().findobj(mpl.image.AxesImage):
            assert isinstance(image, mpl.image.AxesImage)
            clim = np.array(image.get_clim())
            assert np.allclose(expected, clim)

    @pytest.mark.slow
    def test_vmin_vmax_equal(self) -> None:
        # regression test for GH3734
        fg = self.g.map_dataarray(xplt.imshow, "x", "y", vmin=50, vmax=50)
        for mappable in fg._mappables:
            assert mappable.norm.vmin != mappable.norm.vmax

    @pytest.mark.slow
    @pytest.mark.filterwarnings("ignore")
    def test_can_set_norm(self) -> None:
        norm = mpl.colors.SymLogNorm(0.1)
        self.g.map_dataarray(xplt.imshow, "x", "y", norm=norm)
        for image in plt.gcf().findobj(mpl.image.AxesImage):
            assert isinstance(image, mpl.image.AxesImage)
            assert image.norm is norm

    @pytest.mark.slow
    def test_figure_size(self) -> None:
        assert_array_equal(self.g.fig.get_size_inches(), (10, 3))

        g = xplt.FacetGrid(self.darray, col="z", size=6)
        assert_array_equal(g.fig.get_size_inches(), (19, 6))

        g = self.darray.plot.imshow(col="z", size=6)
        assert_array_equal(g.fig.get_size_inches(), (19, 6))

        g = xplt.FacetGrid(self.darray, col="z", size=4, aspect=0.5)
        assert_array_equal(g.fig.get_size_inches(), (7, 4))

        g = xplt.FacetGrid(self.darray, col="z", figsize=(9, 4))
        assert_array_equal(g.fig.get_size_inches(), (9, 4))

        with pytest.raises(ValueError, match=r"cannot provide both"):
            g = xplt.plot(self.darray, row=2, col="z", figsize=(6, 4), size=6)

        with pytest.raises(ValueError, match=r"Can't use"):
            g = xplt.plot(self.darray, row=2, col="z", ax=plt.gca(), size=6)

    @pytest.mark.slow
    def test_num_ticks(self) -> None:
        nticks = 99
        maxticks = nticks + 1
        self.g.map_dataarray(xplt.imshow, "x", "y")
        self.g.set_ticks(max_xticks=nticks, max_yticks=nticks)

        for ax in self.g.axs.flat:
            xticks = len(ax.get_xticks())
            yticks = len(ax.get_yticks())
            assert xticks <= maxticks
            assert yticks <= maxticks
            assert xticks >= nticks / 2.0
            assert yticks >= nticks / 2.0

    @pytest.mark.slow
    def test_map(self) -> None:
        assert self.g._finalized is False
        self.g.map(plt.contourf, "x", "y", ...)
        assert self.g._finalized is True
        self.g.map(lambda: None)

    @pytest.mark.slow
    def test_map_dataset(self) -> None:
        g = xplt.FacetGrid(self.darray.to_dataset(name="foo"), col="z")
        g.map(plt.contourf, "x", "y", "foo")

        alltxt = text_in_fig()
        for label in ["x", "y"]:
            assert label in alltxt
        # everything has a label
        assert "None" not in alltxt

        # colorbar can't be inferred automatically
        assert "foo" not in alltxt
        assert 0 == len(find_possible_colorbars())

        g.add_colorbar(label="colors!")
        assert "colors!" in text_in_fig()
        assert 1 == len(find_possible_colorbars())

    @pytest.mark.slow
    def test_set_axis_labels(self) -> None:
        g = self.g.map_dataarray(xplt.contourf, "x", "y")
        g.set_axis_labels("longitude", "latitude")
        alltxt = text_in_fig()
        for label in ["longitude", "latitude"]:
            assert label in alltxt

    @pytest.mark.slow
    def test_facetgrid_colorbar(self) -> None:
        a = easy_array((10, 15, 4))
        d = DataArray(a, dims=["y", "x", "z"], name="foo")

        d.plot.imshow(x="x", y="y", col="z")
        assert 1 == len(find_possible_colorbars())

        d.plot.imshow(x="x", y="y", col="z", add_colorbar=True)
        assert 1 == len(find_possible_colorbars())

        d.plot.imshow(x="x", y="y", col="z", add_colorbar=False)
        assert 0 == len(find_possible_colorbars())

    @pytest.mark.slow
    def test_facetgrid_polar(self) -> None:
        # test if polar projection in FacetGrid does not raise an exception
        self.darray.plot.pcolormesh(
            col="z", subplot_kws=dict(projection="polar"), sharex=False, sharey=False
        )

    @pytest.mark.slow
    def test_units_appear_somewhere(self) -> None:
        # assign coordinates to all dims so we can test for units
        darray = self.darray.assign_coords(
            {"x": np.arange(self.darray.x.size), "y": np.arange(self.darray.y.size)}
        )

        darray.x.attrs["units"] = "x_unit"
        darray.y.attrs["units"] = "y_unit"

        g = xplt.FacetGrid(darray, col="z")

        g.map_dataarray(xplt.contourf, "x", "y")

        alltxt = text_in_fig()

        # unit should appear as e.g. 'x [x_unit]'
        for unit_name in ["x_unit", "y_unit"]:
            assert unit_name in "".join(alltxt)


@pytest.mark.filterwarnings("ignore:tight_layout cannot")
class TestFacetGrid4d(PlotTestCase):
    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        a = easy_array((10, 15, 3, 2))
        darray = DataArray(a, dims=["y", "x", "col", "row"])
        darray.coords["col"] = np.array(
            ["col" + str(x) for x in darray.coords["col"].values]
        )
        darray.coords["row"] = np.array(
            ["row" + str(x) for x in darray.coords["row"].values]
        )

        self.darray = darray

    def test_title_kwargs(self) -> None:
        g = xplt.FacetGrid(self.darray, col="col", row="row")
        g.set_titles(template="{value}", weight="bold")

        # Rightmost column titles should be bold
        for label, ax in zip(
            self.darray.coords["row"].values, g.axs[:, -1], strict=True
        ):
            assert property_in_axes_text("weight", "bold", label, ax)

        # Top row titles should be bold
        for label, ax in zip(
            self.darray.coords["col"].values, g.axs[0, :], strict=True
        ):
            assert property_in_axes_text("weight", "bold", label, ax)

    @pytest.mark.slow
    def test_default_labels(self) -> None:
        g = xplt.FacetGrid(self.darray, col="col", row="row")
        assert (2, 3) == g.axs.shape

        g.map_dataarray(xplt.imshow, "x", "y")

        # Rightmost column should be labeled
        for label, ax in zip(
            self.darray.coords["row"].values, g.axs[:, -1], strict=True
        ):
            assert substring_in_axes(label, ax)

        # Top row should be labeled
        for label, ax in zip(
            self.darray.coords["col"].values, g.axs[0, :], strict=True
        ):
            assert substring_in_axes(label, ax)

        # ensure that row & col labels can be changed
        g.set_titles("abc={value}")
        for label, ax in zip(
            self.darray.coords["row"].values, g.axs[:, -1], strict=True
        ):
            assert substring_in_axes(f"abc={label}", ax)
            # previous labels were "row=row0" etc.
            assert substring_not_in_axes("row=", ax)

        for label, ax in zip(
            self.darray.coords["col"].values, g.axs[0, :], strict=True
        ):
            assert substring_in_axes(f"abc={label}", ax)
            # previous labels were "col=row0" etc.
            assert substring_not_in_axes("col=", ax)


@pytest.mark.filterwarnings("ignore:tight_layout cannot")
class TestFacetedLinePlotsLegend(PlotTestCase):
    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        self.darray = xr.tutorial.scatter_example_dataset()

    def test_legend_labels(self) -> None:
        fg = self.darray.A.plot.line(col="x", row="w", hue="z")
        all_legend_labels = [t.get_text() for t in fg.figlegend.texts]
        # labels in legend should be ['0', '1', '2', '3']
        assert sorted(all_legend_labels) == ["0", "1", "2", "3"]


@pytest.mark.filterwarnings("ignore:tight_layout cannot")
class TestFacetedLinePlots(PlotTestCase):
    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        self.darray = DataArray(
            np.random.randn(10, 6, 3, 4),
            dims=["hue", "x", "col", "row"],
            coords=[range(10), range(6), range(3), ["A", "B", "C", "C++"]],
            name="Cornelius Ortega the 1st",
        )

        self.darray.hue.name = "huename"
        self.darray.hue.attrs["units"] = "hunits"
        self.darray.x.attrs["units"] = "xunits"
        self.darray.col.attrs["units"] = "colunits"
        self.darray.row.attrs["units"] = "rowunits"

    def test_facetgrid_shape(self) -> None:
        g = self.darray.plot(row="row", col="col", hue="hue")  # type: ignore[call-arg]
        assert g.axs.shape == (len(self.darray.row), len(self.darray.col))

        g = self.darray.plot(row="col", col="row", hue="hue")  # type: ignore[call-arg]
        assert g.axs.shape == (len(self.darray.col), len(self.darray.row))

    def test_unnamed_args(self) -> None:
        g = self.darray.plot.line("o--", row="row", col="col", hue="hue")
        lines = [
            q for q in g.axs.flat[0].get_children() if isinstance(q, mpl.lines.Line2D)
        ]
        # passing 'o--' as argument should set marker and linestyle
        assert lines[0].get_marker() == "o"
        assert lines[0].get_linestyle() == "--"

    def test_default_labels(self) -> None:
        g = self.darray.plot(row="row", col="col", hue="hue")  # type: ignore[call-arg]
        # Rightmost column should be labeled
        for label, ax in zip(
            self.darray.coords["row"].values, g.axs[:, -1], strict=True
        ):
            assert substring_in_axes(label, ax)

        # Top row should be labeled
        for label, ax in zip(
            self.darray.coords["col"].values, g.axs[0, :], strict=True
        ):
            assert substring_in_axes(str(label), ax)

        # Leftmost column should have array name
        for ax in g.axs[:, 0]:
            assert substring_in_axes(str(self.darray.name), ax)

    def test_test_empty_cell(self) -> None:
        g = (
            self.darray.isel(row=1)  # type: ignore[call-arg]
            .drop_vars("row")
            .plot(col="col", hue="hue", col_wrap=2)
        )
        bottomright = g.axs[-1, -1]
        assert not bottomright.has_data()
        assert not bottomright.get_visible()

    def test_set_axis_labels(self) -> None:
        g = self.darray.plot(row="row", col="col", hue="hue")  # type: ignore[call-arg]
        g.set_axis_labels("longitude", "latitude")
        alltxt = text_in_fig()

        assert "longitude" in alltxt
        assert "latitude" in alltxt

    def test_axes_in_faceted_plot(self) -> None:
        with pytest.raises(ValueError):
            self.darray.plot.line(row="row", col="col", x="x", ax=plt.axes())

    def test_figsize_and_size(self) -> None:
        with pytest.raises(ValueError):
            self.darray.plot.line(row="row", col="col", x="x", size=3, figsize=(4, 3))

    def test_wrong_num_of_dimensions(self) -> None:
        with pytest.raises(ValueError):
            self.darray.plot(row="row", hue="hue")  # type: ignore[call-arg]
            self.darray.plot.line(row="row", hue="hue")


@requires_matplotlib
class TestDatasetQuiverPlots(PlotTestCase):
    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        das = [
            DataArray(
                np.random.randn(3, 3, 4, 4),
                dims=["x", "y", "row", "col"],
                coords=[range(k) for k in [3, 3, 4, 4]],
            )
            for _ in [1, 2]
        ]
        ds = Dataset({"u": das[0], "v": das[1]})
        ds.x.attrs["units"] = "xunits"
        ds.y.attrs["units"] = "yunits"
        ds.col.attrs["units"] = "colunits"
        ds.row.attrs["units"] = "rowunits"
        ds.u.attrs["units"] = "uunits"
        ds.v.attrs["units"] = "vunits"
        ds["mag"] = np.hypot(ds.u, ds.v)
        self.ds = ds

    def test_quiver(self) -> None:
        with figure_context():
            hdl = self.ds.isel(row=0, col=0).plot.quiver(x="x", y="y", u="u", v="v")
            assert isinstance(hdl, mpl.quiver.Quiver)
        with pytest.raises(ValueError, match=r"specify x, y, u, v"):
            self.ds.isel(row=0, col=0).plot.quiver(x="x", y="y", u="u")

        with pytest.raises(ValueError, match=r"hue_style"):
            self.ds.isel(row=0, col=0).plot.quiver(
                x="x", y="y", u="u", v="v", hue="mag", hue_style="discrete"
            )

    def test_facetgrid(self) -> None:
        with figure_context():
            fg = self.ds.plot.quiver(
                x="x", y="y", u="u", v="v", row="row", col="col", scale=1, hue="mag"
            )
            for handle in fg._mappables:
                assert isinstance(handle, mpl.quiver.Quiver)
            assert fg.quiverkey is not None
            assert "uunits" in fg.quiverkey.text.get_text()

        with figure_context():
            fg = self.ds.plot.quiver(
                x="x",
                y="y",
                u="u",
                v="v",
                row="row",
                col="col",
                scale=1,
                hue="mag",
                add_guide=False,
            )
            assert fg.quiverkey is None
        with pytest.raises(ValueError, match=r"Please provide scale"):
            self.ds.plot.quiver(x="x", y="y", u="u", v="v", row="row", col="col")

    @pytest.mark.parametrize(
        "add_guide, hue_style, legend, colorbar",
        [
            (None, None, False, True),
            (False, None, False, False),
            (True, None, False, True),
            (True, "continuous", False, True),
        ],
    )
    def test_add_guide(self, add_guide, hue_style, legend, colorbar) -> None:
        meta_data = _infer_meta_data(
            self.ds,
            x="x",
            y="y",
            hue="mag",
            hue_style=hue_style,
            add_guide=add_guide,
            funcname="quiver",
        )
        assert meta_data["add_legend"] is legend
        assert meta_data["add_colorbar"] is colorbar


@requires_matplotlib
class TestDatasetStreamplotPlots(PlotTestCase):
    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        das = [
            DataArray(
                np.random.randn(3, 4, 2, 2),
                dims=["x", "y", "row", "col"],
                coords=[range(k) for k in [3, 4, 2, 2]],
            )
            for _ in [1, 2]
        ]
        ds = Dataset({"u": das[0], "v": das[1]})
        ds.x.attrs["units"] = "xunits"
        ds.y.attrs["units"] = "yunits"
        ds.col.attrs["units"] = "colunits"
        ds.row.attrs["units"] = "rowunits"
        ds.u.attrs["units"] = "uunits"
        ds.v.attrs["units"] = "vunits"
        ds["mag"] = np.hypot(ds.u, ds.v)
        self.ds = ds

    def test_streamline(self) -> None:
        with figure_context():
            hdl = self.ds.isel(row=0, col=0).plot.streamplot(x="x", y="y", u="u", v="v")
            assert isinstance(hdl, mpl.collections.LineCollection)
        with pytest.raises(ValueError, match=r"specify x, y, u, v"):
            self.ds.isel(row=0, col=0).plot.streamplot(x="x", y="y", u="u")

        with pytest.raises(ValueError, match=r"hue_style"):
            self.ds.isel(row=0, col=0).plot.streamplot(
                x="x", y="y", u="u", v="v", hue="mag", hue_style="discrete"
            )

    def test_facetgrid(self) -> None:
        with figure_context():
            fg = self.ds.plot.streamplot(
                x="x", y="y", u="u", v="v", row="row", col="col", hue="mag"
            )
            for handle in fg._mappables:
                assert isinstance(handle, mpl.collections.LineCollection)

        with figure_context():
            fg = self.ds.plot.streamplot(
                x="x",
                y="y",
                u="u",
                v="v",
                row="row",
                col="col",
                hue="mag",
                add_guide=False,
            )


@requires_matplotlib
class TestDatasetScatterPlots(PlotTestCase):
    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        das = [
            DataArray(
                np.random.randn(3, 3, 4, 4),
                dims=["x", "row", "col", "hue"],
                coords=[range(k) for k in [3, 3, 4, 4]],
            )
            for _ in [1, 2]
        ]
        ds = Dataset({"A": das[0], "B": das[1]})
        ds.hue.name = "huename"
        ds.hue.attrs["units"] = "hunits"
        ds.x.attrs["units"] = "xunits"
        ds.col.attrs["units"] = "colunits"
        ds.row.attrs["units"] = "rowunits"
        ds.A.attrs["units"] = "Aunits"
        ds.B.attrs["units"] = "Bunits"
        self.ds = ds

    def test_accessor(self) -> None:
        from xarray.plot.accessor import DatasetPlotAccessor

        assert Dataset.plot is DatasetPlotAccessor
        assert isinstance(self.ds.plot, DatasetPlotAccessor)

    @pytest.mark.parametrize(
        "add_guide, hue_style, legend, colorbar",
        [
            (None, None, False, True),
            (False, None, False, False),
            (True, None, False, True),
            (True, "continuous", False, True),
            (False, "discrete", False, False),
            (True, "discrete", True, False),
        ],
    )
    def test_add_guide(
        self,
        add_guide: bool | None,
        hue_style: Literal["continuous", "discrete"] | None,
        legend: bool,
        colorbar: bool,
    ) -> None:
        meta_data = _infer_meta_data(
            self.ds,
            x="A",
            y="B",
            hue="hue",
            hue_style=hue_style,
            add_guide=add_guide,
            funcname="scatter",
        )
        assert meta_data["add_legend"] is legend
        assert meta_data["add_colorbar"] is colorbar

    def test_facetgrid_shape(self) -> None:
        g = self.ds.plot.scatter(x="A", y="B", row="row", col="col")
        assert g.axs.shape == (len(self.ds.row), len(self.ds.col))

        g = self.ds.plot.scatter(x="A", y="B", row="col", col="row")
        assert g.axs.shape == (len(self.ds.col), len(self.ds.row))

    def test_default_labels(self) -> None:
        g = self.ds.plot.scatter(x="A", y="B", row="row", col="col", hue="hue")

        # Top row should be labeled
        for label, ax in zip(self.ds.coords["col"].values, g.axs[0, :], strict=True):
            assert substring_in_axes(str(label), ax)

        # Bottom row should have name of x array name and units
        for ax in g.axs[-1, :]:
            assert ax.get_xlabel() == "A [Aunits]"

        # Leftmost column should have name of y array name and units
        for ax in g.axs[:, 0]:
            assert ax.get_ylabel() == "B [Bunits]"

    def test_axes_in_faceted_plot(self) -> None:
        with pytest.raises(ValueError):
            self.ds.plot.scatter(x="A", y="B", row="row", ax=plt.axes())

    def test_figsize_and_size(self) -> None:
        with pytest.raises(ValueError):
            self.ds.plot.scatter(x="A", y="B", row="row", size=3, figsize=(4, 3))

    @pytest.mark.parametrize(
        "x, y, hue, add_legend, add_colorbar, error_type",
        [
            pytest.param(
                "A", "The Spanish Inquisition", None, None, None, KeyError, id="bad_y"
            ),
            pytest.param(
                "The Spanish Inquisition", "B", None, None, True, ValueError, id="bad_x"
            ),
        ],
    )
    def test_bad_args(
        self,
        x: Hashable,
        y: Hashable,
        hue: Hashable | None,
        add_legend: bool | None,
        add_colorbar: bool | None,
        error_type: type[Exception],
    ) -> None:
        with pytest.raises(error_type):
            self.ds.plot.scatter(
                x=x, y=y, hue=hue, add_legend=add_legend, add_colorbar=add_colorbar
            )

    def test_datetime_hue(self) -> None:
        ds2 = self.ds.copy()

        # TODO: Currently plots as categorical, should it behave as numerical?
        ds2["hue"] = pd.date_range("2000-1-1", periods=4)
        ds2.plot.scatter(x="A", y="B", hue="hue")

        ds2["hue"] = pd.timedelta_range("-1D", periods=4, freq="D")
        ds2.plot.scatter(x="A", y="B", hue="hue")

    def test_facetgrid_hue_style(self) -> None:
        ds2 = self.ds.copy()

        # Numbers plots as continuous:
        g = ds2.plot.scatter(x="A", y="B", row="row", col="col", hue="hue")
        assert isinstance(g._mappables[-1], mpl.collections.PathCollection)

        # Datetimes plots as categorical:
        # TODO: Currently plots as categorical, should it behave as numerical?
        ds2["hue"] = pd.date_range("2000-1-1", periods=4)
        g = ds2.plot.scatter(x="A", y="B", row="row", col="col", hue="hue")
        assert isinstance(g._mappables[-1], mpl.collections.PathCollection)

        # Strings plots as categorical:
        ds2["hue"] = ["a", "a", "b", "b"]
        g = ds2.plot.scatter(x="A", y="B", row="row", col="col", hue="hue")
        assert isinstance(g._mappables[-1], mpl.collections.PathCollection)

    @pytest.mark.parametrize(
        ["x", "y", "hue", "markersize"],
        [("A", "B", "x", "col"), ("x", "row", "A", "B")],
    )
    def test_scatter(
        self, x: Hashable, y: Hashable, hue: Hashable, markersize: Hashable
    ) -> None:
        self.ds.plot.scatter(x=x, y=y, hue=hue, markersize=markersize)

        with pytest.raises(ValueError, match=r"u, v"):
            self.ds.plot.scatter(x=x, y=y, u="col", v="row")

    def test_non_numeric_legend(self) -> None:
        ds2 = self.ds.copy()
        ds2["hue"] = ["a", "b", "c", "d"]
        pc = ds2.plot.scatter(x="A", y="B", markersize="hue")
        axes = pc.axes
        assert axes is not None
        # should make a discrete legend
        assert hasattr(axes, "legend_")
        assert axes.legend_ is not None

    def test_legend_labels(self) -> None:
        # regression test for #4126: incorrect legend labels
        ds2 = self.ds.copy()
        ds2["hue"] = ["a", "a", "b", "b"]
        pc = ds2.plot.scatter(x="A", y="B", markersize="hue")
        axes = pc.axes
        assert axes is not None
        legend = axes.get_legend()
        assert legend is not None
        actual = [t.get_text() for t in legend.texts]
        expected = ["hue", "a", "b"]
        assert actual == expected

    def test_legend_labels_facetgrid(self) -> None:
        ds2 = self.ds.copy()
        ds2["hue"] = ["d", "a", "c", "b"]
        g = ds2.plot.scatter(x="A", y="B", hue="hue", markersize="x", col="col")
        legend = g.figlegend
        assert legend is not None
        actual = tuple(t.get_text() for t in legend.texts)
        expected = (
            "x [xunits]",
            "$\\mathdefault{0}$",
            "$\\mathdefault{1}$",
            "$\\mathdefault{2}$",
        )
        assert actual == expected

    def test_add_legend_by_default(self) -> None:
        sc = self.ds.plot.scatter(x="A", y="B", hue="hue")
        fig = sc.figure
        assert fig is not None
        assert len(fig.axes) == 2


class TestDatetimePlot(PlotTestCase):
    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        """
        Create a DataArray with a time-axis that contains datetime objects.
        """
        month = np.arange(1, 13, 1)
        data = np.sin(2 * np.pi * month / 12.0)
        times = pd.date_range(start="2017-01-01", freq="MS", periods=12)
        darray = DataArray(data, dims=["time"], coords=[times])

        self.darray = darray

    def test_datetime_line_plot(self) -> None:
        # test if line plot raises no Exception
        self.darray.plot.line()

    def test_datetime_units(self) -> None:
        # test that matplotlib-native datetime works:
        fig, ax = plt.subplots()
        ax.plot(self.darray["time"], self.darray)

        # Make sure only mpl converters are used, use type() so only
        # mpl.dates.AutoDateLocator passes and no other subclasses:
        assert type(ax.xaxis.get_major_locator()) is mpl.dates.AutoDateLocator

    def test_datetime_plot1d(self) -> None:
        # Test that matplotlib-native datetime works:
        p = self.darray.plot.line()
        ax = p[0].axes

        # Make sure only mpl converters are used, use type() so only
        # mpl.dates.AutoDateLocator passes and no other subclasses:
        assert type(ax.xaxis.get_major_locator()) is mpl.dates.AutoDateLocator

    def test_datetime_plot2d(self) -> None:
        # Test that matplotlib-native datetime works:
        da = DataArray(
            np.arange(3 * 4).reshape(3, 4),
            dims=("x", "y"),
            coords={
                "x": [1, 2, 3],
                "y": [np.datetime64(f"2000-01-{x:02d}") for x in range(1, 5)],
            },
        )

        p = da.plot.pcolormesh()
        ax = p.axes
        assert ax is not None

        # Make sure only mpl converters are used, use type() so only
        # mpl.dates.AutoDateLocator passes and no other subclasses:
        assert type(ax.xaxis.get_major_locator()) is mpl.dates.AutoDateLocator


@pytest.mark.filterwarnings("ignore:setting an array element with a sequence")
@requires_cftime
@pytest.mark.skipif(not has_nc_time_axis, reason="nc_time_axis is not installed")
class TestCFDatetimePlot(PlotTestCase):
    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        """
        Create a DataArray with a time-axis that contains cftime.datetime
        objects.
        """
        # case for 1d array
        data = np.random.rand(4, 12)
        time = xr.date_range(
            start="2017", periods=12, freq="1ME", calendar="noleap", use_cftime=True
        )
        darray = DataArray(data, dims=["x", "time"])
        darray.coords["time"] = time

        self.darray = darray

    def test_cfdatetime_line_plot(self) -> None:
        self.darray.isel(x=0).plot.line()

    def test_cfdatetime_pcolormesh_plot(self) -> None:
        self.darray.plot.pcolormesh()

    def test_cfdatetime_contour_plot(self) -> None:
        self.darray.plot.contour()


@requires_cftime
@pytest.mark.skipif(has_nc_time_axis, reason="nc_time_axis is installed")
class TestNcAxisNotInstalled(PlotTestCase):
    @pytest.fixture(autouse=True)
    def setUp(self) -> None:
        """
        Create a DataArray with a time-axis that contains cftime.datetime
        objects.
        """
        month = np.arange(1, 13, 1)
        data = np.sin(2 * np.pi * month / 12.0)
        darray = DataArray(data, dims=["time"])
        darray.coords["time"] = xr.date_range(
            start="2017", periods=12, freq="1ME", calendar="noleap", use_cftime=True
        )

        self.darray = darray

    def test_ncaxis_notinstalled_line_plot(self) -> None:
        with pytest.raises(ImportError, match=r"optional `nc-time-axis`"):
            self.darray.plot.line()


@requires_matplotlib
class TestAxesKwargs:
    @pytest.fixture(params=[1, 2, 3])
    def data_array(self, request) -> DataArray:
        """
        Return a simple DataArray
        """
        dims = request.param
        if dims == 1:
            return DataArray(easy_array((10,)))
        elif dims == 2:
            return DataArray(easy_array((10, 3)))
        elif dims == 3:
            return DataArray(easy_array((10, 3, 2)))
        else:
            raise ValueError(f"No DataArray implemented for {dims=}.")

    @pytest.fixture(params=[1, 2])
    def data_array_logspaced(self, request) -> DataArray:
        """
        Return a simple DataArray with logspaced coordinates
        """
        dims = request.param
        if dims == 1:
            return DataArray(
                np.arange(7), dims=("x",), coords={"x": np.logspace(-3, 3, 7)}
            )
        elif dims == 2:
            return DataArray(
                np.arange(16).reshape(4, 4),
                dims=("y", "x"),
                coords={"x": np.logspace(-1, 2, 4), "y": np.logspace(-5, -1, 4)},
            )
        else:
            raise ValueError(f"No DataArray implemented for {dims=}.")

    @pytest.mark.parametrize("xincrease", [True, False])
    def test_xincrease_kwarg(self, data_array, xincrease) -> None:
        with figure_context():
            data_array.plot(xincrease=xincrease)
            assert plt.gca().xaxis_inverted() == (not xincrease)

    @pytest.mark.parametrize("yincrease", [True, False])
    def test_yincrease_kwarg(self, data_array, yincrease) -> None:
        with figure_context():
            data_array.plot(yincrease=yincrease)
            assert plt.gca().yaxis_inverted() == (not yincrease)

    @pytest.mark.parametrize("xscale", ["linear", "logit", "symlog"])
    def test_xscale_kwarg(self, data_array, xscale) -> None:
        with figure_context():
            data_array.plot(xscale=xscale)
            assert plt.gca().get_xscale() == xscale

    @pytest.mark.parametrize("yscale", ["linear", "logit", "symlog"])
    def test_yscale_kwarg(self, data_array, yscale) -> None:
        with figure_context():
            data_array.plot(yscale=yscale)
            assert plt.gca().get_yscale() == yscale

    def test_xscale_log_kwarg(self, data_array_logspaced) -> None:
        xscale = "log"
        with figure_context():
            data_array_logspaced.plot(xscale=xscale)
            assert plt.gca().get_xscale() == xscale

    def test_yscale_log_kwarg(self, data_array_logspaced) -> None:
        yscale = "log"
        with figure_context():
            data_array_logspaced.plot(yscale=yscale)
            assert plt.gca().get_yscale() == yscale

    def test_xlim_kwarg(self, data_array) -> None:
        with figure_context():
            expected = (0.0, 1000.0)
            data_array.plot(xlim=[0, 1000])
            assert plt.gca().get_xlim() == expected

    def test_ylim_kwarg(self, data_array) -> None:
        with figure_context():
            data_array.plot(ylim=[0, 1000])
            expected = (0.0, 1000.0)
            assert plt.gca().get_ylim() == expected

    def test_xticks_kwarg(self, data_array) -> None:
        with figure_context():
            data_array.plot(xticks=np.arange(5))
            expected = np.arange(5).tolist()
            assert_array_equal(plt.gca().get_xticks(), expected)

    def test_yticks_kwarg(self, data_array) -> None:
        with figure_context():
            data_array.plot(yticks=np.arange(5))
            expected = np.arange(5)
            assert_array_equal(plt.gca().get_yticks(), expected)


@requires_matplotlib
@pytest.mark.parametrize("plotfunc", ["pcolormesh", "contourf", "contour"])
def test_plot_transposed_nondim_coord(plotfunc) -> None:
    x = np.linspace(0, 10, 101)
    h = np.linspace(3, 7, 101)
    s = np.linspace(0, 1, 51)
    z = s[:, np.newaxis] * h[np.newaxis, :]
    da = xr.DataArray(
        np.sin(x) * np.cos(z),
        dims=["s", "x"],
        coords={"x": x, "s": s, "z": (("s", "x"), z), "zt": (("x", "s"), z.T)},
    )
    with figure_context():
        getattr(da.plot, plotfunc)(x="x", y="zt")
    with figure_context():
        getattr(da.plot, plotfunc)(x="zt", y="x")


@requires_matplotlib
@pytest.mark.parametrize("plotfunc", ["pcolormesh", "imshow"])
def test_plot_transposes_properly(plotfunc) -> None:
    # test that we aren't mistakenly transposing when the 2 dimensions have equal sizes.
    da = xr.DataArray([np.sin(2 * np.pi / 10 * np.arange(10))] * 10, dims=("y", "x"))
    with figure_context():
        hdl = getattr(da.plot, plotfunc)(x="x", y="y")
        # get_array doesn't work for contour, contourf. It returns the colormap intervals.
        # pcolormesh returns 1D array but imshow returns a 2D array so it is necessary
        # to ravel() on the LHS
        assert_array_equal(hdl.get_array().ravel(), da.to_masked_array().ravel())


@requires_matplotlib
def test_facetgrid_single_contour() -> None:
    # regression test for GH3569
    x, y = np.meshgrid(np.arange(12), np.arange(12))
    z = xr.DataArray(np.hypot(x, y))
    z2 = xr.DataArray(np.hypot(x, y) + 1)
    ds = xr.concat([z, z2], dim="time")
    ds["time"] = [0, 1]

    with figure_context():
        ds.plot.contour(col="time", levels=[4], colors=["k"])


@requires_matplotlib
def test_get_axis_raises() -> None:
    # test get_axis raises an error if trying to do invalid things

    # cannot provide both ax and figsize
    with pytest.raises(ValueError, match="both `figsize` and `ax`"):
        get_axis(figsize=[4, 4], size=None, aspect=None, ax="something")  # type: ignore[arg-type]

    # cannot provide both ax and size
    with pytest.raises(ValueError, match="both `size` and `ax`"):
        get_axis(figsize=None, size=200, aspect=4 / 3, ax="something")  # type: ignore[arg-type]

    # cannot provide both size and figsize
    with pytest.raises(ValueError, match="both `figsize` and `size`"):
        get_axis(figsize=[4, 4], size=200, aspect=None, ax=None)

    # cannot provide aspect and size
    with pytest.raises(ValueError, match="`aspect` argument without `size`"):
        get_axis(figsize=None, size=None, aspect=4 / 3, ax=None)

    # cannot provide axis and subplot_kws
    with pytest.raises(ValueError, match="cannot use subplot_kws with existing ax"):
        get_axis(figsize=None, size=None, aspect=None, ax=1, something_else=5)  # type: ignore[arg-type]


@requires_matplotlib
@pytest.mark.parametrize(
    ["figsize", "size", "aspect", "ax", "kwargs"],
    [
        pytest.param((3, 2), None, None, False, {}, id="figsize"),
        pytest.param(
            (3.5, 2.5), None, None, False, {"label": "test"}, id="figsize_kwargs"
        ),
        pytest.param(None, 5, None, False, {}, id="size"),
        pytest.param(None, 5.5, None, False, {"label": "test"}, id="size_kwargs"),
        pytest.param(None, 5, 1, False, {}, id="size+aspect"),
        pytest.param(None, 5, "auto", False, {}, id="auto_aspect"),
        pytest.param(None, 5, "equal", False, {}, id="equal_aspect"),
        pytest.param(None, None, None, True, {}, id="ax"),
        pytest.param(None, None, None, False, {}, id="default"),
        pytest.param(None, None, None, False, {"label": "test"}, id="default_kwargs"),
    ],
)
def test_get_axis(
    figsize: tuple[float, float] | None,
    size: float | None,
    aspect: float | None,
    ax: bool,
    kwargs: dict[str, Any],
) -> None:
    with figure_context():
        inp_ax = plt.axes() if ax else None
        out_ax = get_axis(
            figsize=figsize, size=size, aspect=aspect, ax=inp_ax, **kwargs
        )
        assert isinstance(out_ax, mpl.axes.Axes)


@requires_matplotlib
@requires_cartopy
@pytest.mark.parametrize(
    ["figsize", "size", "aspect"],
    [
        pytest.param((3, 2), None, None, id="figsize"),
        pytest.param(None, 5, None, id="size"),
        pytest.param(None, 5, 1, id="size+aspect"),
        pytest.param(None, None, None, id="default"),
    ],
)
def test_get_axis_cartopy(
    figsize: tuple[float, float] | None, size: float | None, aspect: float | None
) -> None:
    kwargs = {"projection": cartopy.crs.PlateCarree()}
    with figure_context():
        out_ax = get_axis(figsize=figsize, size=size, aspect=aspect, **kwargs)
        assert isinstance(out_ax, cartopy.mpl.geoaxes.GeoAxesSubplot)


@requires_matplotlib
def test_get_axis_current() -> None:
    with figure_context():
        _, ax = plt.subplots()
        out_ax = get_axis()
        assert ax is out_ax


@requires_matplotlib
def test_maybe_gca() -> None:
    with figure_context():
        ax = _maybe_gca(aspect=1)

        assert isinstance(ax, mpl.axes.Axes)
        assert ax.get_aspect() == 1

    with figure_context():
        # create figure without axes
        plt.figure()
        ax = _maybe_gca(aspect=1)

        assert isinstance(ax, mpl.axes.Axes)
        assert ax.get_aspect() == 1

    with figure_context():
        existing_axes = plt.axes()
        ax = _maybe_gca(aspect=1)

        # reuses the existing axes
        assert existing_axes == ax
        # kwargs are ignored when reusing axes
        assert ax.get_aspect() == "auto"


@requires_matplotlib
@pytest.mark.parametrize(
    "x, y, z, hue, markersize, row, col, add_legend, add_colorbar",
    [
        ("A", "B", None, None, None, None, None, None, None),
        ("B", "A", None, "w", None, None, None, True, None),
        ("A", "B", None, "y", "x", None, None, True, True),
        ("A", "B", "z", None, None, None, None, None, None),
        ("B", "A", "z", "w", None, None, None, True, None),
        ("A", "B", "z", "y", "x", None, None, True, True),
        ("A", "B", "z", "y", "x", "w", None, True, True),
    ],
)
def test_datarray_scatter(
    x, y, z, hue, markersize, row, col, add_legend, add_colorbar
) -> None:
    """Test datarray scatter. Merge with TestPlot1D eventually."""
    ds = xr.tutorial.scatter_example_dataset()

    extra_coords = [v for v in [x, hue, markersize] if v is not None]

    # Base coords:
    coords = dict(ds.coords)

    # Add extra coords to the DataArray:
    coords.update({v: ds[v] for v in extra_coords})

    darray = xr.DataArray(ds[y], coords=coords)

    with figure_context():
        darray.plot.scatter(
            x=x,
            z=z,
            hue=hue,
            markersize=markersize,
            add_legend=add_legend,
            add_colorbar=add_colorbar,
        )


@requires_dask
@requires_matplotlib
@pytest.mark.parametrize(
    "plotfunc",
    ["scatter"],
)
def test_dataarray_not_loading_inplace(plotfunc: str) -> None:
    ds = xr.tutorial.scatter_example_dataset()
    ds = ds.chunk()

    with figure_context():
        getattr(ds.A.plot, plotfunc)(x="x")

    from dask.array import Array

    assert isinstance(ds.A.data, Array)


@requires_matplotlib
def test_assert_valid_xy() -> None:
    ds = xr.tutorial.scatter_example_dataset()
    darray = ds.A

    # x is valid and should not error:
    _assert_valid_xy(darray=darray, xy="x", name="x")

    # None should be valid as well even though it isn't in the valid list:
    _assert_valid_xy(darray=darray, xy=None, name="x")

    # A hashable that is not valid should error:
    with pytest.raises(ValueError, match="x must be one of"):
        _assert_valid_xy(darray=darray, xy="error_now", name="x")


@requires_matplotlib
@pytest.mark.parametrize(
    "val", [pytest.param([], id="empty"), pytest.param(0, id="scalar")]
)
@pytest.mark.parametrize(
    "method",
    [
        "__call__",
        "line",
        "step",
        "contour",
        "contourf",
        "hist",
        "imshow",
        "pcolormesh",
        "scatter",
        "surface",
    ],
)
def test_plot_empty_raises(val: list | float, method: str) -> None:
    da = xr.DataArray(val)
    with pytest.raises(TypeError, match="No numeric data"):
        getattr(da.plot, method)()


@requires_matplotlib
def test_facetgrid_axes_raises_deprecation_warning() -> None:
    with pytest.warns(
        DeprecationWarning,
        match=(
            "self.axes is deprecated since 2022.11 in order to align with "
            "matplotlibs plt.subplots, use self.axs instead."
        ),
    ):
        with figure_context():
            ds = xr.tutorial.scatter_example_dataset()
            g = ds.plot.scatter(x="A", y="B", col="x")
            _ = g.axes


@requires_matplotlib
def test_plot1d_default_rcparams() -> None:
    import matplotlib as mpl

    ds = xr.tutorial.scatter_example_dataset(seed=42)

    with figure_context():
        # scatter markers should by default have white edgecolor to better
        # see overlapping markers:
        fig, ax = plt.subplots(1, 1)
        ds.plot.scatter(x="A", y="B", marker="o", ax=ax)
        actual: np.ndarray = mpl.colors.to_rgba_array("w")
        expected: np.ndarray = ax.collections[0].get_edgecolor()  # type: ignore[assignment]
        np.testing.assert_allclose(actual, expected)

        # Facetgrids should have the default value as well:
        fg = ds.plot.scatter(x="A", y="B", col="x", marker="o")
        ax = fg.axs.ravel()[0]
        actual = mpl.colors.to_rgba_array("w")
        expected = ax.collections[0].get_edgecolor()  # type: ignore[assignment,unused-ignore]
        np.testing.assert_allclose(actual, expected)

        # scatter should not emit any warnings when using unfilled markers:
        with assert_no_warnings():
            fig, ax = plt.subplots(1, 1)
            ds.plot.scatter(x="A", y="B", ax=ax, marker="x")

        # Prioritize edgecolor argument over default plot1d values:
        fig, ax = plt.subplots(1, 1)
        ds.plot.scatter(x="A", y="B", marker="o", ax=ax, edgecolor="k")
        actual = mpl.colors.to_rgba_array("k")
        expected = ax.collections[0].get_edgecolor()  # type: ignore[assignment]
        np.testing.assert_allclose(actual, expected)


@requires_matplotlib
def test_plot1d_filtered_nulls() -> None:
    ds = xr.tutorial.scatter_example_dataset(seed=42)
    y = ds.y.where(ds.y > 0.2)
    expected = y.notnull().sum().item()

    with figure_context():
        pc = y.plot.scatter()
        actual = pc.get_offsets().shape[0]

        assert expected == actual


@requires_matplotlib
def test_9155() -> None:
    # A test for types from issue #9155

    with figure_context():
        data = xr.DataArray([1, 2, 3], dims=["x"])
        fig, ax = plt.subplots(ncols=1, nrows=1)
        data.plot(ax=ax)  # type: ignore[call-arg]


@requires_matplotlib
def test_temp_dataarray() -> None:
    from xarray.plot.dataset_plot import _temp_dataarray

    x = np.arange(1, 4)
    y = np.arange(4, 6)
    var1 = np.arange(x.size * y.size).reshape((x.size, y.size))
    var2 = np.arange(x.size * y.size).reshape((x.size, y.size))
    ds = xr.Dataset(
        {
            "var1": (["x", "y"], var1),
            "var2": (["x", "y"], 2 * var2),
            "var3": (["x"], 3 * x),
        },
        coords={
            "x": x,
            "y": y,
            "model": np.arange(7),
        },
    )

    # No broadcasting:
    y_ = "var1"
    locals_ = {"x": "var2"}
    da = _temp_dataarray(ds, y_, locals_)
    assert da.shape == (3, 2)

    # Broadcast from 1 to 2dim:
    y_ = "var3"
    locals_ = {"x": "var1"}
    da = _temp_dataarray(ds, y_, locals_)
    assert da.shape == (3, 2)

    # Ignore non-valid coord kwargs:
    y_ = "var3"
    locals_ = dict(x="x", extend="var2")
    da = _temp_dataarray(ds, y_, locals_)
    assert da.shape == (3,)
