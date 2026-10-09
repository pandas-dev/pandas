"""Test cases for DataFrame.plot"""

import pickle
import re

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.tests.plotting.common import (
    _check_colors,
    _check_plot_works,
    _unpack_cycler,
)

mpl = pytest.importorskip("matplotlib")
plt = pytest.importorskip("matplotlib.pyplot")
cm = pytest.importorskip("matplotlib.cm")


@pytest.mark.parametrize("as_column", [False, True])
@pytest.mark.parametrize(
    "values",
    [
        [0, 0, 0],
        [1, 1, 1],
        [-5, -5, -5],
        [5, 5, 5],
        [1e16, 1e16, 1e16],
        [0, 1, 0],
        [np.nan, np.nan, np.nan],
        [1, 2, np.inf],
        [5, 5, np.inf],
        [5, np.nan, 5],
        [1e-300, 1e-300, 1e-300],
        pd.array([1, pd.NA, 1], dtype="Int64"),
        [5],
    ],
)
def test_scatter_colorbar_consistent_mapping(values, as_column):
    # GH 64980: adding a colorbar must preserve the scatter's color mapping.
    df = pd.DataFrame({"x": np.arange(len(values)) + 1, "c": values})
    cmap = mpl.colors.ListedColormap(["blue", "red"])
    with tm.assert_produces_warning(False):
        expected_fig, expected_ax = plt.subplots()
        expected = expected_ax.scatter(df["x"], df["x"], c=df["c"], cmap=cmap)
        expected_fig.canvas.draw()

        ax = df.plot.scatter(
            "x", "x", c="c" if as_column else df["c"], cmap=cmap, colorbar=True
        )
        ax.figure.canvas.draw()

    result = ax.collections[0]
    assert result.get_clim() == expected.get_clim()
    tm.assert_numpy_array_equal(result.get_facecolors(), expected.get_facecolors())
    assert len(ax.figure.axes) == 2
    assert result.colorbar.mappable is result
    assert result.colorbar.norm is result.norm
    if result.norm.vmin == result.norm.vmax:
        tm.assert_numpy_array_equal(
            result.colorbar.solids.get_facecolors(), np.array([cmap(0.0)])
        )
        tm.assert_numpy_array_equal(
            result.colorbar.get_ticks(), np.array([result.norm.vmin])
        )


@pytest.mark.parametrize("values", [[1, 2, 3], [5, 5, 5]])
@pytest.mark.parametrize(
    "limits",
    [
        {"vmin": 0, "vmax": 10},
        {"vmin": 0},
        {"vmax": 10},
        {"vmin": 0, "vmax": 0},
        {"vmin": 5, "vmax": 5},
    ],
)
def test_scatter_colorbar_explicit_limits(values, limits):
    # GH 64980: do not introduce a norm alongside vmin/vmax.
    df = pd.DataFrame({"x": [1, 2, 3], "c": values})
    ax = df.plot.scatter("x", "x", c="c", colorbar=True, **limits)
    ax.figure.canvas.draw()
    expected = (limits.get("vmin", min(values)), limits.get("vmax", max(values)))
    assert ax.collections[0].get_clim() == expected


@pytest.mark.parametrize(
    "norm_factory",
    [
        lambda: mpl.colors.Normalize(),
        lambda: mpl.colors.Normalize(1, 100),
        lambda: mpl.colors.Normalize(5, 5),
        lambda: mpl.colors.Normalize(0, 0),
        lambda: mpl.colors.LogNorm(),
        lambda: mpl.colors.LogNorm(1, 100),
        lambda: mpl.colors.LogNorm(5, 5),
        lambda: mpl.colors.LogNorm(0, 0),
        lambda: mpl.colors.CenteredNorm(),
        lambda: mpl.colors.SymLogNorm(linthresh=1),
        lambda: mpl.colors.PowerNorm(1),
        lambda: mpl.colors.AsinhNorm(),
        lambda: mpl.colors.TwoSlopeNorm(vcenter=0),
    ],
)
def test_scatter_colorbar_custom_norm(norm_factory):
    # GH 64980: preserve the user's norm and its mapping, including log scales.
    df = pd.DataFrame({"x": [1, 2, 3], "c": [5, 5, 5]})
    norm = norm_factory()
    limits = (norm.vmin, norm.vmax) if norm.scaled() else None
    expected_fig, expected_ax = plt.subplots()
    expected = expected_ax.scatter(
        df["x"], df["x"], c=df["c"], cmap="viridis", norm=norm_factory()
    )
    expected_fig.canvas.draw()

    ax = df.plot.scatter("x", "x", c="c", cmap="viridis", norm=norm, colorbar=True)
    ax.figure.canvas.draw()
    result = ax.collections[0]
    assert result.norm is norm
    assert result.get_clim() == expected.get_clim()
    tm.assert_numpy_array_equal(result.get_facecolors(), expected.get_facecolors())
    if limits is not None:
        assert (norm.vmin, norm.vmax) == limits


@pytest.mark.parametrize("value", [0, -5, np.nan, 1e-300])
def test_scatter_no_colorbar_log_norm(value):
    # Preserve matplotlib's handling of nonpositive, missing, and tiny values
    # in a scatter plot without a colorbar.
    df = pd.DataFrame({"x": [1, 2, 3], "c": [value] * 3})
    expected_fig, expected_ax = plt.subplots()
    expected = expected_ax.scatter(
        df["x"], df["x"], c=df["c"], cmap="viridis", norm=mpl.colors.LogNorm()
    )
    expected_fig.canvas.draw()

    ax = df.plot.scatter(
        "x", "x", c="c", cmap="viridis", norm=mpl.colors.LogNorm(), colorbar=False
    )
    ax.figure.canvas.draw()
    assert ax.collections[0].get_clim() == expected.get_clim()
    tm.assert_numpy_array_equal(
        ax.collections[0].get_facecolors(), expected.get_facecolors()
    )


@pytest.mark.parametrize("colorbar", [False, True])
def test_scatter_colorbar_categorical_mapping(colorbar):
    # GH 64980: do not change categorical BoundaryNorm or its labels.
    df = pd.DataFrame({"x": [1, 2, 3], "c": pd.Categorical(["a", "b", "a"])})
    ax = df.plot.scatter("x", "x", c="c", cmap="viridis", colorbar=colorbar)
    ax.figure.canvas.draw()
    scatter = ax.collections[0]
    assert isinstance(scatter.norm, mpl.colors.BoundaryNorm)
    tm.assert_numpy_array_equal(scatter.norm.boundaries, np.array([0.0, 1.0, 2.0]))
    expected = mpl.colormaps["viridis"](np.array([0, 255, 0]))
    tm.assert_numpy_array_equal(scatter.get_facecolors(), expected)
    if colorbar:
        assert type(scatter.colorbar) is mpl.colorbar.Colorbar
        labels = [label.get_text() for label in scatter.colorbar.ax.get_yticklabels()]
        assert labels == ["a", "b"]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"c": [0, 1, 0]},
        {"c": [0, 0, 0], "vmin": 0, "vmax": 10},
        {"c": [1, 5, 10], "norm": mpl.colors.LogNorm(1, 10)},
        {"c": [0, 1, 0], "norm": mpl.colors.BoundaryNorm([0, 1, 2], 256)},
        {"c": [0, 0, 0], "norm": mpl.colors.NoNorm(0, 0)},
        {"c": ["red", "blue", "red"]},
    ],
)
def test_scatter_colorbar_standard_class(kwargs):
    # Only a singular continuous norm needs pandas' custom colorbar.
    df = pd.DataFrame({"x": [1, 2, 3]})
    ax = df.plot.scatter("x", "x", colorbar=True, **kwargs)
    ax.figure.canvas.draw()
    scatter = ax.collections[0]
    assert type(scatter.colorbar) is mpl.colorbar.Colorbar
    assert scatter.colorbar.mappable is scatter


def _constant_scatter():
    df = pd.DataFrame({"x": [1, 2, 3], "c": [0, 0, 0]})
    ax = df.plot.scatter("x", "x", c="c", cmap="viridis")
    scatter = ax.collections[0]
    return ax, scatter, scatter.colorbar


def test_scatter_colorbar_constant_to_range():
    # GH 64980: widening the norm switches from a swatch to a gradient.
    ax, scatter, colorbar = _constant_scatter()
    assert colorbar.mappable is scatter
    assert colorbar.ax.get_ylabel() == "c"

    scatter.set_clim(0, 10)
    ax.figure.canvas.draw()
    assert colorbar.norm is scatter.norm
    assert colorbar.ax.get_ylim() == (0, 10)
    assert colorbar.ax.get_ylabel() == "c"

    scatter.set_clim(0, 0)
    ax.figure.canvas.draw()
    with scatter.norm.callbacks.blocked():
        scatter.set_clim(0, 10)
    colorbar.update_normal()
    ax.figure.canvas.draw()
    assert colorbar.norm is scatter.norm
    assert colorbar.ax.get_ylim() == (0, 10)
    assert colorbar.values is None


def test_scatter_colorbar_manual_update_range_to_constant():
    # GH 64980: a manual update_normal applies a blocked return to a constant.
    cmap = mpl.colors.ListedColormap(["blue", "red"])
    df = pd.DataFrame({"x": [1, 2, 3], "c": [0, 0, 0]})
    ax = df.plot.scatter("x", "x", c="c", cmap=cmap)
    scatter = ax.collections[0]
    colorbar = scatter.colorbar

    scatter.set_clim(0, 10)
    scatter.set_array(np.array([5, 5, 5]))
    with scatter.norm.callbacks.blocked():
        scatter.set_clim(5, 5)
    ax.figure.canvas.draw()
    colorbar.update_normal()
    ax.figure.canvas.draw()

    assert scatter.get_clim() == (5, 5)
    tm.assert_numpy_array_equal(scatter.get_facecolors(), np.array([cmap(0.0)] * 3))
    tm.assert_numpy_array_equal(colorbar.solids.get_facecolors(), np.array([cmap(0.0)]))
    assert colorbar.get_ticks().tolist() == [5]


def test_scatter_colorbar_pickle_manual_refresh():
    # GH 64980: manual update_normal still refreshes a colorbar after pickling.
    ax, *_ = _constant_scatter()
    ax.figure.canvas.draw()
    fig2 = pickle.loads(pickle.dumps(ax.figure))
    scatter2 = fig2.axes[0].collections[0]
    colorbar2 = scatter2.colorbar

    scatter2.set_clim(0, 10)
    colorbar2.update_normal()
    fig2.canvas.draw()
    assert colorbar2.ax.get_ylim() == (0, 10)

    scatter2.set_clim(5, 5)
    colorbar2.update_normal()
    fig2.canvas.draw()
    assert scatter2.get_clim() == (5, 5)
    assert colorbar2.get_ticks().tolist() == [5]


def test_scatter_colorbar_range_to_constant():
    # A later constant range must not be expanded by the callback either.
    ax, scatter, colorbar = _constant_scatter()
    scatter.set_clim(0, 10)
    scatter.set_clim(5, 5)
    scatter.set_array(np.array([5, 5, 5]))
    ax.figure.canvas.draw()
    assert scatter.get_clim() == (5, 5)
    assert colorbar.ax.get_ylabel() == "c"
    tm.assert_numpy_array_equal(
        colorbar.solids.get_facecolors(), np.array([scatter.cmap(0.0)])
    )
    assert colorbar.get_ticks().tolist() == [5]

    scatter.set_array(np.array([0, 1, 0]))
    scatter.autoscale()
    ax.figure.canvas.draw()
    assert scatter.get_clim() == (0, 1)
    assert colorbar.ax.get_ylim() == (0, 1)


def test_scatter_colorbar_cmap_keeps_formatter():
    ax, scatter, colorbar = _constant_scatter()
    colorbar.formatter = mpl.ticker.FuncFormatter(lambda value, pos: "constant")
    scatter.set_cmap("plasma")
    ax.figure.canvas.draw()
    assert colorbar.ax.get_yticklabels()[0].get_text() == "constant"
    tm.assert_numpy_array_equal(
        colorbar.solids.get_facecolors(), np.array([scatter.cmap(0.0)])
    )


def test_scatter_colorbar_set_norm_log():
    ax, scatter, colorbar = _constant_scatter()
    values = np.array([1, 10, 100])
    norm = mpl.colors.LogNorm(1, 100)
    scatter.set_array(values)
    scatter.set_norm(norm)
    scatter.set_cmap("plasma")
    ax.figure.canvas.draw()
    assert colorbar.ax.get_yscale() == "log"
    assert colorbar.norm is norm
    assert colorbar.cmap is scatter.cmap
    tm.assert_almost_equal(colorbar.ax.get_ylim(), (1, 100))
    tm.assert_numpy_array_equal(
        scatter.get_facecolors(), mpl.colormaps["plasma"](norm(values))
    )


def test_scatter_colorbar_remove():
    ax, scatter, colorbar = _constant_scatter()
    colorbar.remove()
    assert scatter.colorbar is None
    assert not scatter.callbacks.callbacks.get("changed")
    scatter.set_clim(0, 0)
    ax.figure.canvas.draw()
    assert scatter.get_clim() == (0, 0)


@pytest.mark.parametrize(
    "layout", [None, "constrained", "tight", "manual", "subfigure"]
)
def test_scatter_colorbar_layout(layout):
    # The custom colorbar must use Figure.colorbar's axes placement rules.
    df = pd.DataFrame({"x": [1, 2, 3], "c": [0, 0, 0]})
    positions = []
    for use_pandas in [False, True]:
        fig = plt.figure(layout=layout if layout in ["constrained", "tight"] else None)
        if layout == "manual":
            ax = fig.add_axes([0.15, 0.15, 0.7, 0.7])
        elif layout == "subfigure":
            ax = fig.subfigures(1, 1).subplots()
        else:
            ax = fig.subplots()
        if use_pandas:
            df.plot.scatter("x", "x", c="c", ax=ax)
        else:
            scatter = ax.scatter(df["x"], df["x"], c=df["c"])
            ax.set_xlabel("x")
            ax.set_ylabel("x")
            lower, upper = mpl.ticker.AutoLocator().nonsingular(0, 0)
            fig.colorbar(
                scatter,
                ax=ax,
                boundaries=[lower, upper],
                values=[0],
                ticks=[0],
                label="c",
            )
        colorbar = ax.collections[0].colorbar
        assert ax.get_figure(root=False).gca() is ax
        fig.canvas.draw()
        positions.append((ax.get_position().bounds, colorbar.ax.get_position().bounds))
        colorbar.remove()
        fig.canvas.draw()
    tm.assert_almost_equal(positions[0], positions[1])


@pytest.mark.parametrize("norm_type", [mpl.colors.Normalize, mpl.colors.LogNorm])
def test_scatter_colorbar_replaces_unscaled_norm(norm_type):
    # GH 64980: replacing the norm with an unscaled one must autoscale to the
    # data before the colorbar decides whether the range is constant.
    df = pd.DataFrame({"x": [1, 2, 3], "c": [5, 5, 5]})
    cmap = mpl.colors.ListedColormap(["blue", "red"])

    def plotted(norm, *, colorbar):
        ax = df.plot.scatter("x", "x", c="c", cmap=cmap, colorbar=colorbar)
        scatter = ax.collections[0]
        scatter.set_norm(norm)
        ax.figure.canvas.draw()
        return scatter

    norm = norm_type()
    expected = plotted(norm_type(), colorbar=False)
    result = plotted(norm, colorbar=True)
    assert result.norm is norm
    assert result.get_clim() == expected.get_clim()
    tm.assert_numpy_array_equal(result.get_facecolors(), expected.get_facecolors())

    result.norm.vmin = result.norm.vmax = None
    result.autoscale_None()
    result.axes.figure.canvas.draw()
    expected.norm.vmin = expected.norm.vmax = None
    expected.autoscale_None()
    expected.axes.figure.canvas.draw()
    assert result.norm is norm
    assert result.get_clim() == expected.get_clim()
    tm.assert_numpy_array_equal(result.get_facecolors(), expected.get_facecolors())


def test_scatter_colorbar_set_norm_without_array():
    # GH 64980: clearing the array must not make autoscale raise.
    df = pd.DataFrame({"x": [1, 2, 3], "c": [5, 5, 5]})
    cmap = mpl.colors.ListedColormap(["blue", "red"])

    expected_fig, expected_ax = plt.subplots()
    expected = expected_ax.scatter(df["x"], df["x"], c=df["c"], cmap=cmap)
    expected_fig.colorbar(expected, ax=expected_ax)
    expected.set_array(None)
    expected.set_norm(mpl.colors.Normalize())
    expected_fig.canvas.draw()

    ax = df.plot.scatter("x", "x", c="c", cmap=cmap)
    result = ax.collections[0]
    result.set_array(None)
    result.set_norm(mpl.colors.Normalize())
    ax.figure.canvas.draw()
    assert result.get_clim() == expected.get_clim()
    tm.assert_numpy_array_equal(result.get_facecolors(), expected.get_facecolors())


def test_scatter_no_colorbar_near_constant():
    # Numerically close but distinct values must not be expanded when no
    # colorbar is requested (a regression in an earlier GH 64980 proposal).
    values = [1e16, 1e16 + 2, 1e16]
    df = pd.DataFrame({"x": [1, 2, 3], "c": values})
    norm = mpl.colors.Normalize(min(values), max(values))
    ax = df.plot.scatter("x", "x", c="c", norm=norm, colorbar=False)
    ax.figure.canvas.draw()
    assert ax.collections[0].norm is norm
    assert (norm.vmin, norm.vmax) == (min(values), max(values))


@pytest.mark.parametrize("extend", ["neither", "min", "max", "both"])
def test_scatter_colorbar_constant_extensions(extend):
    # Colorbar extensions need their own boundaries and color values.
    cmap = mpl.colors.ListedColormap(["blue", "red"])
    cmap.colorbar_extend = extend
    df = pd.DataFrame({"x": [1, 2, 3], "c": [0, 0, 0]})
    ax = df.plot.scatter("x", "x", c="c", cmap=cmap)
    scatter = ax.collections[0]
    colorbar = scatter.colorbar
    for value in [0, 5]:
        scatter.set_clim(value, value)
        ax.figure.canvas.draw()
        assert colorbar.extend == extend
        assert scatter.get_clim() == (value, value)
        tm.assert_numpy_array_equal(
            colorbar.solids.get_facecolors(), np.array([cmap(0.0)])
        )


def _check_colors_box(bp, box_c, whiskers_c, medians_c, caps_c="k", fliers_c=None):
    if fliers_c is None:
        fliers_c = "k"
    _check_colors(bp["boxes"], linecolors=[box_c] * len(bp["boxes"]))
    _check_colors(bp["whiskers"], linecolors=[whiskers_c] * len(bp["whiskers"]))
    _check_colors(bp["medians"], linecolors=[medians_c] * len(bp["medians"]))
    _check_colors(bp["fliers"], linecolors=[fliers_c] * len(bp["fliers"]))
    _check_colors(bp["caps"], linecolors=[caps_c] * len(bp["caps"]))


class TestDataFrameColor:
    @pytest.mark.parametrize("color", list(range(10)))
    def test_mpl2_color_cycle_str(self, color):
        # GH 15516
        color = f"C{color}"
        df = pd.DataFrame(
            np.random.default_rng(2).standard_normal((10, 3)), columns=["a", "b", "c"]
        )
        _check_plot_works(df.plot, color=color)

    def test_color_single_series_list(self):
        # GH 3486
        df = pd.DataFrame({"A": [1, 2, 3]})
        _check_plot_works(df.plot, color=["red"])

    @pytest.mark.parametrize("color", [(1, 0, 0), (1, 0, 0, 0.5)])
    def test_rgb_tuple_color(self, color):
        # GH 16695
        df = pd.DataFrame({"x": [1, 2], "y": [3, 4]})
        _check_plot_works(df.plot, x="x", y="y", color=color)

    def test_color_empty_string(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((10, 2)))
        with pytest.raises(ValueError, match="Invalid color argument:"):
            df.plot(color="")

    def test_color_and_style_arguments(self):
        df = pd.DataFrame({"x": [1, 2], "y": [3, 4]})
        # passing both 'color' and 'style' arguments should be allowed
        # if there is no color symbol in the style strings:
        ax = df.plot(color=["red", "black"], style=["-", "--"])
        # check that the linestyles are correctly set:
        linestyle = [line.get_linestyle() for line in ax.lines]
        assert linestyle == ["-", "--"]
        # check that the colors are correctly set:
        color = [line.get_color() for line in ax.lines]
        assert color == ["red", "black"]
        # passing both 'color' and 'style' arguments should not be allowed
        # if there is a color symbol in the style strings:
        msg = (
            "Cannot pass 'style' string with a color symbol and 'color' keyword "
            "argument. Please use one or the other or pass 'style' without a color "
            "symbol"
        )
        with pytest.raises(ValueError, match=msg):
            df.plot(color=["red", "black"], style=["k-", "r--"])

    @pytest.mark.parametrize(
        "color, expected",
        [
            ("green", ["green"] * 4),
            (["yellow", "red", "green", "blue"], ["yellow", "red", "green", "blue"]),
        ],
    )
    def test_color_and_marker(self, color, expected):
        # GH 21003
        df = pd.DataFrame(np.random.default_rng(2).random((7, 4)))
        ax = df.plot(color=color, style="d--")
        # check colors
        result = [i.get_color() for i in ax.lines]
        assert result == expected
        # check markers and linestyles
        assert all(i.get_linestyle() == "--" for i in ax.lines)
        assert all(i.get_marker() == "d" for i in ax.lines)

    def test_color_and_style(self):
        color = {"g": "black", "h": "brown"}
        style = {"g": "-", "h": "--"}
        expected_color = ["black", "brown"]
        expected_style = ["-", "--"]
        df = pd.DataFrame({"g": [1, 2], "h": [2, 3]}, index=[1, 2])
        ax = df.plot.line(color=color, style=style)
        color = [i.get_color() for i in ax.lines]
        style = [i.get_linestyle() for i in ax.lines]
        assert color == expected_color
        assert style == expected_style

    def test_bar_colors(self):
        default_colors = _unpack_cycler(plt.rcParams)

        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        ax = df.plot.bar()
        _check_colors(ax.patches[::5], facecolors=default_colors[:5])

    def test_bar_colors_custom(self):
        custom_colors = "rgcby"
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        ax = df.plot.bar(color=custom_colors)
        _check_colors(ax.patches[::5], facecolors=custom_colors)

    @pytest.mark.parametrize("colormap", ["jet", cm.jet])
    def test_bar_colors_cmap(self, colormap):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))

        ax = df.plot.bar(colormap=colormap)
        rgba_colors = [cm.jet(n) for n in np.linspace(0, 1, 5)]
        _check_colors(ax.patches[::5], facecolors=rgba_colors)

    def test_bar_colors_single_col(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        ax = df.loc[:, [0]].plot.bar(color="DodgerBlue")
        _check_colors([ax.patches[0]], facecolors=["DodgerBlue"])

    def test_bar_colors_single_col_list(self):
        # GH#18006 color maps per column, so a single-column frame uses only
        # the first color of the list; use a Series or y= (see
        # test_bar_colors_single_col_y_list) for one color per bar
        df = pd.DataFrame({"Test": [1, 3, 5, 7]})
        colors = [(0.9, 0.9, 0.4), (0.8, 0.4, 0.6), (0.2, 0.7, 0.9), (0.4, 0.4, 0.5)]
        ax = df.plot.bar(color=colors)
        _check_colors(ax.patches, facecolors=[colors[0]] * 4)

    def test_bar_colors_single_col_y_list(self):
        # GH#18006 selecting a single column with y= colors each bar
        # individually, unlike the per-column behavior in
        # test_bar_colors_single_col_list
        df = pd.DataFrame({"Test": [1, 3, 5, 7]})
        colors = [(0.9, 0.9, 0.4), (0.8, 0.4, 0.6), (0.2, 0.7, 0.9), (0.4, 0.4, 0.5)]
        ax = df.plot.bar(y="Test", color=colors)
        _check_colors(ax.patches, facecolors=colors)

    def test_bar_colors_green(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        ax = df.plot(kind="bar", color="green")
        _check_colors(ax.patches[::5], facecolors=["green"] * 5)

    def test_bar_user_colors(self):
        df = pd.DataFrame(
            {"A": range(4), "B": range(1, 5), "color": ["red", "blue", "blue", "red"]}
        )
        # This should *only* work when `y` is specified, else
        # we use one color per column
        ax = df.plot.bar(y="A", color=df["color"])
        result = [p.get_facecolor() for p in ax.patches]
        expected = [
            (1.0, 0.0, 0.0, 1.0),
            (0.0, 0.0, 1.0, 1.0),
            (0.0, 0.0, 1.0, 1.0),
            (1.0, 0.0, 0.0, 1.0),
        ]
        assert result == expected

    def test_if_scatterplot_colorbar_affects_xaxis_visibility(self):
        # addressing issue #10611, to ensure colobar does not
        # interfere with x-axis label and ticklabels with
        # ipython inline backend.
        random_array = np.random.default_rng(2).random((10, 3))
        df = pd.DataFrame(random_array, columns=["A label", "B label", "C label"])

        ax1 = df.plot.scatter(x="A label", y="B label")
        ax2 = df.plot.scatter(x="A label", y="B label", c="C label")

        vis1 = [vis.get_visible() for vis in ax1.xaxis.get_minorticklabels()]
        vis2 = [vis.get_visible() for vis in ax2.xaxis.get_minorticklabels()]
        assert vis1 == vis2

        vis1 = [vis.get_visible() for vis in ax1.xaxis.get_majorticklabels()]
        vis2 = [vis.get_visible() for vis in ax2.xaxis.get_majorticklabels()]
        assert vis1 == vis2

        assert (
            ax1.xaxis.get_label().get_visible() == ax2.xaxis.get_label().get_visible()
        )

    def test_if_hexbin_xaxis_label_is_visible(self):
        # addressing issue #10678, to ensure colobar does not
        # interfere with x-axis label and ticklabels with
        # ipython inline backend.
        random_array = np.random.default_rng(2).random((10, 3))
        df = pd.DataFrame(random_array, columns=["A label", "B label", "C label"])

        ax = df.plot.hexbin("A label", "B label", gridsize=12)
        assert all(vis.get_visible() for vis in ax.xaxis.get_minorticklabels())
        assert all(vis.get_visible() for vis in ax.xaxis.get_majorticklabels())
        assert ax.xaxis.get_label().get_visible()

    def test_if_scatterplot_colorbars_are_next_to_parent_axes(self):
        random_array = np.random.default_rng(2).random((10, 3))
        df = pd.DataFrame(random_array, columns=["A label", "B label", "C label"])

        fig, axes = plt.subplots(1, 2)
        df.plot.scatter("A label", "B label", c="C label", ax=axes[0])
        df.plot.scatter("A label", "B label", c="C label", ax=axes[1])
        plt.tight_layout()

        points = np.array([ax.get_position().get_points() for ax in fig.axes])
        axes_x_coords = points[:, :, 0]
        parent_distance = axes_x_coords[1, :] - axes_x_coords[0, :]
        colorbar_distance = axes_x_coords[3, :] - axes_x_coords[2, :]
        assert np.isclose(parent_distance, colorbar_distance, atol=1e-7).all()

    @pytest.mark.parametrize("cmap", [None, "Greys"])
    def test_scatter_with_c_column_name_with_colors(self, cmap):
        # https://github.com/pandas-dev/pandas/issues/34316

        df = pd.DataFrame(
            [[5.1, 3.5], [4.9, 3.0], [7.0, 3.2], [6.4, 3.2], [5.9, 3.0]],
            columns=["length", "width"],
        )
        df["species"] = ["r", "r", "g", "g", "b"]
        if cmap is not None:
            msg = "No data for colormapping provided via 'c'"
            with tm.assert_produces_warning(
                UserWarning, check_stacklevel=False, match=msg
            ):
                ax = df.plot.scatter(x=0, y=1, cmap=cmap, c="species")
        else:
            ax = df.plot.scatter(x=0, y=1, c="species", cmap=cmap)

        assert len(np.unique(ax.collections[0].get_facecolor(), axis=0)) == 3  # r/g/b
        assert (
            np.unique(ax.collections[0].get_facecolor(), axis=0)
            == np.array(
                [[0.0, 0.0, 1.0, 1.0], [0.0, 0.5, 0.0, 1.0], [1.0, 0.0, 0.0, 1.0]]
            )  # r/g/b
        ).all()
        assert ax.collections[0].colorbar is None

    def test_scatter_with_c_column_name_without_colors(self):
        # Given
        colors = ["NY", "MD", "MA", "CA"]
        color_count = 4  # 4 unique colors

        # When
        df = pd.DataFrame(
            {
                "dataX": range(100),
                "dataY": range(100),
                "color": (colors[i % len(colors)] for i in range(100)),
            }
        )

        # Then
        ax = df.plot.scatter("dataX", "dataY", c="color")
        assert len(np.unique(ax.collections[0].get_facecolor(), axis=0)) == color_count

        # Given
        colors = ["r", "g", "not-a-color"]
        color_count = 3
        # Also, since not all are mpl-colors, points matching 'r' or 'g'
        # are not necessarily red or green

        # When
        df = pd.DataFrame(
            {
                "dataX": range(100),
                "dataY": range(100),
                "color": (colors[i % len(colors)] for i in range(100)),
            }
        )

        # Then
        ax = df.plot.scatter("dataX", "dataY", c="color")
        assert len(np.unique(ax.collections[0].get_facecolor(), axis=0)) == color_count

    def test_scatter_colors(self):
        df = pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3], "c": [1, 2, 3]})
        with pytest.raises(TypeError, match="Specify exactly one of `c` and `color`"):
            df.plot.scatter(x="a", y="b", c="c", color="green")

    def test_scatter_colors_not_raising_warnings(self):
        # GH-53908. Do not raise UserWarning: No data for colormapping
        # provided via 'c'. Parameters 'cmap' will be ignored
        df = pd.DataFrame({"x": [1, 2, 3], "y": [1, 2, 3]})
        with tm.assert_produces_warning(None):
            ax = df.plot.scatter(x="x", y="y", c="b")
            assert (
                len(np.unique(ax.collections[0].get_facecolor(), axis=0)) == 1
            )  # blue
            assert (
                np.unique(ax.collections[0].get_facecolor(), axis=0)
                == np.array([[0.0, 0.0, 1.0, 1.0]])
            ).all()  # blue

    def test_scatter_colors_default(self):
        df = pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3], "c": [1, 2, 3]})
        default_colors = _unpack_cycler(mpl.pyplot.rcParams)

        ax = df.plot.scatter(x="a", y="b", c="c")
        tm.assert_numpy_array_equal(
            ax.collections[0].get_facecolor()[0],
            np.array(mpl.colors.ColorConverter.to_rgba(default_colors[0])),
        )

    def test_scatter_colors_white(self):
        df = pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3], "c": [1, 2, 3]})
        ax = df.plot.scatter(x="a", y="b", color="white")
        tm.assert_numpy_array_equal(
            ax.collections[0].get_facecolor()[0],
            np.array([1, 1, 1, 1], dtype=np.float64),
        )

    def test_scatter_colorbar_different_cmap(self):
        # GH 33389
        df = pd.DataFrame({"x": [1, 2, 3], "y": [1, 3, 2], "c": [1, 2, 3]})
        df["x2"] = df["x"] + 1

        _, ax = plt.subplots()
        df.plot("x", "y", c="c", kind="scatter", cmap="cividis", ax=ax)
        df.plot("x2", "y", c="c", kind="scatter", cmap="magma", ax=ax)

        assert ax.collections[0].cmap.name == "cividis"
        assert ax.collections[1].cmap.name == "magma"

    def test_line_colors(self):
        custom_colors = "rgcby"
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))

        ax = df.plot(color=custom_colors)
        _check_colors(ax.get_lines(), linecolors=custom_colors)

        plt.close("all")

        ax2 = df.plot(color=custom_colors)
        lines2 = ax2.get_lines()

        for l1, l2 in zip(ax.get_lines(), lines2, strict=True):
            assert l1.get_color() == l2.get_color()

    @pytest.mark.parametrize("colormap", ["jet", cm.jet])
    def test_line_colors_cmap(self, colormap):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        ax = df.plot(colormap=colormap)
        rgba_colors = [cm.jet(n) for n in np.linspace(0, 1, len(df))]
        _check_colors(ax.get_lines(), linecolors=rgba_colors)

    def test_line_colors_single_col(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        # make color a list if plotting one column frame
        # handles cases like df.plot(color='DodgerBlue')
        ax = df.loc[:, [0]].plot(color="DodgerBlue")
        _check_colors(ax.lines, linecolors=["DodgerBlue"])

    def test_line_colors_single_color(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        ax = df.plot(color="red")
        _check_colors(ax.get_lines(), linecolors=["red"] * 5)

    def test_line_colors_hex(self):
        # GH 10299
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        custom_colors = ["#FF0000", "#0000FF", "#FFFF00", "#000000", "#FFFFFF"]
        ax = df.plot(color=custom_colors)
        _check_colors(ax.get_lines(), linecolors=custom_colors)

    def test_dont_modify_colors(self):
        colors = ["r", "g", "b"]
        pd.DataFrame(np.random.default_rng(2).random((10, 2))).plot(color=colors)
        assert len(colors) == 3

    def test_line_colors_and_styles_subplots(self):
        # GH 9894
        default_colors = _unpack_cycler(mpl.pyplot.rcParams)

        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))

        axes = df.plot(subplots=True)
        for ax, c in zip(axes, list(default_colors), strict=False):
            _check_colors(ax.get_lines(), linecolors=[c])

    @pytest.mark.parametrize("color", ["k", "green"])
    def test_line_colors_and_styles_subplots_single_color_str(self, color):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        axes = df.plot(subplots=True, color=color)
        for ax in axes:
            _check_colors(ax.get_lines(), linecolors=[color])

    @pytest.mark.parametrize("color", ["rgcby", list("rgcby")])
    def test_line_colors_and_styles_subplots_custom_colors(self, color):
        # GH 9894
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        axes = df.plot(color=color, subplots=True)
        for ax, c in zip(axes, list(color), strict=True):
            _check_colors(ax.get_lines(), linecolors=[c])

    def test_line_colors_and_styles_subplots_colormap_hex(self):
        # GH 9894
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        # GH 10299
        custom_colors = ["#FF0000", "#0000FF", "#FFFF00", "#000000", "#FFFFFF"]
        axes = df.plot(color=custom_colors, subplots=True)
        for ax, c in zip(axes, list(custom_colors), strict=True):
            _check_colors(ax.get_lines(), linecolors=[c])

    @pytest.mark.parametrize("cmap", ["jet", cm.jet])
    def test_line_colors_and_styles_subplots_colormap_subplot(self, cmap):
        # GH 9894
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        rgba_colors = [cm.jet(n) for n in np.linspace(0, 1, len(df))]
        axes = df.plot(colormap=cmap, subplots=True)
        for ax, c in zip(axes, rgba_colors, strict=True):
            _check_colors(ax.get_lines(), linecolors=[c])

    def test_line_colors_and_styles_subplots_single_col(self):
        # GH 9894
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        # make color a list if plotting one column frame
        # handles cases like df.plot(color='DodgerBlue')
        axes = df.loc[:, [0]].plot(color="DodgerBlue", subplots=True)
        _check_colors(axes[0].lines, linecolors=["DodgerBlue"])

    def test_line_colors_and_styles_subplots_single_char(self):
        # GH 9894
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        # single character style
        axes = df.plot(style="r", subplots=True)
        for ax in axes:
            _check_colors(ax.get_lines(), linecolors=["r"])

    def test_line_colors_and_styles_subplots_list_styles(self):
        # GH 9894
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        # list of styles
        styles = list("rgcby")
        axes = df.plot(style=styles, subplots=True)
        for ax, c in zip(axes, styles, strict=True):
            _check_colors(ax.get_lines(), linecolors=[c])

    def test_area_colors(self):
        custom_colors = "rgcby"
        df = pd.DataFrame(np.random.default_rng(2).random((5, 5)))

        ax = df.plot.area(color=custom_colors)
        _check_colors(ax.get_lines(), linecolors=custom_colors)
        poly = [
            o
            for o in ax.get_children()
            if isinstance(o, mpl.collections.PolyCollection)
        ]
        _check_colors(poly, facecolors=custom_colors)

        handles, _ = ax.get_legend_handles_labels()
        _check_colors(handles, facecolors=custom_colors)

        for h in handles:
            assert h.get_alpha() is None

    def test_area_colors_poly(self):
        df = pd.DataFrame(np.random.default_rng(2).random((5, 5)))
        ax = df.plot.area(colormap="jet")
        jet_colors = [mpl.cm.jet(n) for n in np.linspace(0, 1, len(df))]
        _check_colors(ax.get_lines(), linecolors=jet_colors)
        poly = [
            o
            for o in ax.get_children()
            if isinstance(o, mpl.collections.PolyCollection)
        ]
        _check_colors(poly, facecolors=jet_colors)

        handles, _ = ax.get_legend_handles_labels()
        _check_colors(handles, facecolors=jet_colors)
        for h in handles:
            assert h.get_alpha() is None

    def test_area_colors_stacked_false(self):
        df = pd.DataFrame(np.random.default_rng(2).random((5, 5)))
        jet_colors = [mpl.cm.jet(n) for n in np.linspace(0, 1, len(df))]
        # When stacked=False, alpha is set to 0.5
        ax = df.plot.area(colormap=mpl.cm.jet, stacked=False)
        _check_colors(ax.get_lines(), linecolors=jet_colors)
        poly = [
            o
            for o in ax.get_children()
            if isinstance(o, mpl.collections.PolyCollection)
        ]
        jet_with_alpha = [(c[0], c[1], c[2], 0.5) for c in jet_colors]
        _check_colors(poly, facecolors=jet_with_alpha)

        handles, _ = ax.get_legend_handles_labels()
        linecolors = jet_with_alpha
        _check_colors(handles[: len(jet_colors)], linecolors=linecolors)
        for h in handles:
            assert h.get_alpha() == 0.5

    def test_hist_colors(self):
        default_colors = _unpack_cycler(mpl.pyplot.rcParams)

        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        ax = df.plot.hist()
        _check_colors(ax.patches[::10], facecolors=default_colors[:5])

    def test_hist_colors_single_custom(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        custom_colors = "rgcby"
        ax = df.plot.hist(color=custom_colors)
        _check_colors(ax.patches[::10], facecolors=custom_colors)

    @pytest.mark.parametrize("colormap", ["jet", cm.jet])
    def test_hist_colors_cmap(self, colormap):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        ax = df.plot.hist(colormap=colormap)
        rgba_colors = [cm.jet(n) for n in np.linspace(0, 1, 5)]
        _check_colors(ax.patches[::10], facecolors=rgba_colors)

    def test_hist_colors_single_col(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        ax = df.loc[:, [0]].plot.hist(color="DodgerBlue")
        _check_colors([ax.patches[0]], facecolors=["DodgerBlue"])

    def test_hist_colors_single_color(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        ax = df.plot(kind="hist", color="green")
        _check_colors(ax.patches[::10], facecolors=["green"] * 5)

    def test_kde_colors(self):
        pytest.importorskip("scipy")
        custom_colors = "rgcby"
        df = pd.DataFrame(np.random.default_rng(2).random((5, 5)))

        ax = df.plot.kde(color=custom_colors)
        _check_colors(ax.get_lines(), linecolors=custom_colors)

    @pytest.mark.parametrize("colormap", ["jet", cm.jet])
    def test_kde_colors_cmap(self, colormap):
        pytest.importorskip("scipy")
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        ax = df.plot.kde(colormap=colormap)
        rgba_colors = [cm.jet(n) for n in np.linspace(0, 1, len(df))]
        _check_colors(ax.get_lines(), linecolors=rgba_colors)

    def test_kde_colors_and_styles_subplots(self):
        pytest.importorskip("scipy")
        default_colors = _unpack_cycler(mpl.pyplot.rcParams)

        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))

        axes = df.plot(kind="kde", subplots=True)
        for ax, c in zip(axes, list(default_colors), strict=False):
            _check_colors(ax.get_lines(), linecolors=[c])

    @pytest.mark.parametrize("colormap", ["k", "red"])
    def test_kde_colors_and_styles_subplots_single_col_str(self, colormap):
        pytest.importorskip("scipy")
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        axes = df.plot(kind="kde", color=colormap, subplots=True)
        for ax in axes:
            _check_colors(ax.get_lines(), linecolors=[colormap])

    def test_kde_colors_and_styles_subplots_custom_color(self):
        pytest.importorskip("scipy")
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        custom_colors = "rgcby"
        axes = df.plot(kind="kde", color=custom_colors, subplots=True)
        for ax, c in zip(axes, list(custom_colors), strict=True):
            _check_colors(ax.get_lines(), linecolors=[c])

    @pytest.mark.parametrize("colormap", ["jet", cm.jet])
    def test_kde_colors_and_styles_subplots_cmap(self, colormap):
        pytest.importorskip("scipy")
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        rgba_colors = [cm.jet(n) for n in np.linspace(0, 1, len(df))]
        axes = df.plot(kind="kde", colormap=colormap, subplots=True)
        for ax, c in zip(axes, rgba_colors, strict=True):
            _check_colors(ax.get_lines(), linecolors=[c])

    def test_kde_colors_and_styles_subplots_single_col(self):
        pytest.importorskip("scipy")
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        # make color a list if plotting one column frame
        # handles cases like df.plot(color='DodgerBlue')
        axes = df.loc[:, [0]].plot(kind="kde", color="DodgerBlue", subplots=True)
        _check_colors(axes[0].lines, linecolors=["DodgerBlue"])

    def test_kde_colors_and_styles_subplots_single_char(self):
        pytest.importorskip("scipy")
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        # list of styles
        # single character style
        axes = df.plot(kind="kde", style="r", subplots=True)
        for ax in axes:
            _check_colors(ax.get_lines(), linecolors=["r"])

    def test_kde_colors_and_styles_subplots_list(self):
        pytest.importorskip("scipy")
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        # list of styles
        styles = list("rgcby")
        axes = df.plot(kind="kde", style=styles, subplots=True)
        for ax, c in zip(axes, styles, strict=True):
            _check_colors(ax.get_lines(), linecolors=[c])

    def test_boxplot_colors(self):
        default_colors = _unpack_cycler(mpl.pyplot.rcParams)

        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        bp = df.plot.box(return_type="dict")
        _check_colors_box(
            bp,
            default_colors[0],
            default_colors[0],
            default_colors[2],
            default_colors[0],
        )

    def test_boxplot_colors_dict_colors(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        dict_colors = {
            "boxes": "#572923",
            "whiskers": "#982042",
            "medians": "#804823",
            "caps": "#123456",
        }
        bp = df.plot.box(color=dict_colors, sym="r+", return_type="dict")
        _check_colors_box(
            bp,
            dict_colors["boxes"],
            dict_colors["whiskers"],
            dict_colors["medians"],
            dict_colors["caps"],
            "r",
        )

    def test_boxplot_colors_default_color(self):
        default_colors = _unpack_cycler(mpl.pyplot.rcParams)
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        # partial colors
        dict_colors = {"whiskers": "c", "medians": "m"}
        bp = df.plot.box(color=dict_colors, return_type="dict")
        _check_colors_box(bp, default_colors[0], "c", "m", default_colors[0])

    @pytest.mark.parametrize("colormap", ["jet", cm.jet])
    def test_boxplot_colors_cmap(self, colormap):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        bp = df.plot.box(colormap=colormap, return_type="dict")
        jet_colors = [cm.jet(n) for n in np.linspace(0, 1, 3)]
        _check_colors_box(
            bp, jet_colors[0], jet_colors[0], jet_colors[2], jet_colors[0]
        )

    def test_boxplot_colors_single(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        # string color is applied to all artists except fliers
        bp = df.plot.box(color="DodgerBlue", return_type="dict")
        _check_colors_box(bp, "DodgerBlue", "DodgerBlue", "DodgerBlue", "DodgerBlue")

    def test_boxplot_colors_tuple(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        # tuple is also applied to all artists except fliers
        bp = df.plot.box(color=(0, 1, 0), sym="#123456", return_type="dict")
        _check_colors_box(bp, (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), "#123456")

    def test_boxplot_colors_invalid(self):
        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 5)))
        msg = re.escape(
            "color dict contains invalid key 'xxxx'. The key must be either "
            "['boxes', 'whiskers', 'medians', 'caps']"
        )
        with pytest.raises(ValueError, match=msg):
            # Color contains invalid key results in ValueError
            df.plot.box(color={"boxes": "red", "xxxx": "blue"})

    def test_default_color_cycle(self):
        import cycler

        colors = list("rgbk")
        plt.rcParams["axes.prop_cycle"] = cycler.cycler("color", colors)

        df = pd.DataFrame(np.random.default_rng(2).standard_normal((5, 3)))
        ax = df.plot()

        expected = _unpack_cycler(plt.rcParams)[:3]
        _check_colors(ax.get_lines(), linecolors=expected)

    def test_no_color_bar(self):
        df = pd.DataFrame(
            {
                "A": np.random.default_rng(2).uniform(size=20),
                "B": np.random.default_rng(2).uniform(size=20),
                "C": np.arange(20) + np.random.default_rng(2).uniform(size=20),
            }
        )
        ax = df.plot.hexbin(x="A", y="B", colorbar=None)
        assert ax.collections[0].colorbar is None

    def test_mixing_cmap_and_colormap_raises(self):
        df = pd.DataFrame(
            {
                "A": np.random.default_rng(2).uniform(size=20),
                "B": np.random.default_rng(2).uniform(size=20),
                "C": np.arange(20) + np.random.default_rng(2).uniform(size=20),
            }
        )
        msg = "Only specify one of `cmap` and `colormap`"
        with pytest.raises(TypeError, match=msg):
            df.plot.hexbin(x="A", y="B", cmap="YlGn", colormap="BuGn")

    def test_passed_bar_colors(self):
        color_tuples = [(0.9, 0, 0, 1), (0, 0.9, 0, 1), (0, 0, 0.9, 1)]
        colormap = mpl.colors.ListedColormap(color_tuples)
        barplot = pd.DataFrame([[1, 2, 3]]).plot(kind="bar", cmap=colormap)
        assert color_tuples == [c.get_facecolor() for c in barplot.patches]

    def test_rcParams_bar_colors(self):
        color_tuples = [(0.9, 0, 0, 1), (0, 0.9, 0, 1), (0, 0, 0.9, 1)]
        with mpl.rc_context(rc={"axes.prop_cycle": mpl.cycler("color", color_tuples)}):
            barplot = pd.DataFrame([[1, 2, 3]]).plot(kind="bar")
        assert color_tuples == [c.get_facecolor() for c in barplot.patches]

    def test_colors_of_columns_with_same_name(self):
        # ISSUE 11136 -> https://github.com/pandas-dev/pandas/issues/11136
        # Creating a DataFrame with duplicate column labels and testing colors of them.
        df = pd.DataFrame({"b": [0, 1, 0], "a": [1, 2, 3]})
        df1 = pd.DataFrame({"a": [2, 4, 6]})
        df_concat = pd.concat([df, df1], axis=1)
        result = df_concat.plot()
        legend = result.get_legend()
        handles = legend.legend_handles
        for legend, line in zip(handles, result.lines, strict=True):
            assert legend.get_color() == line.get_color()

    def test_invalid_colormap(self):
        df = pd.DataFrame(
            np.random.default_rng(2).standard_normal((3, 2)), columns=["A", "B"]
        )
        msg = "|".join(["is not a valid value", "is not a known colormap"])
        with pytest.raises((ValueError, KeyError), match=msg):
            df.plot(colormap="invalid_colormap")

    def test_dataframe_none_color(self):
        # GH51953
        df = pd.DataFrame([[1, 2, 3]])
        ax = df.plot(color=None)
        expected = _unpack_cycler(mpl.pyplot.rcParams)[:3]
        _check_colors(ax.get_lines(), linecolors=expected)
