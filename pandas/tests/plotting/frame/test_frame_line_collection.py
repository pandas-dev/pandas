"""
A line plot draws all its columns as one collection where it can (GH#61532).
These check when that path is taken and that it draws what the per-column path
would have drawn.
"""

import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    Series,
    date_range,
    timedelta_range,
)
import pandas._testing as tm
from pandas.tests.plotting.common import (
    _check_legend_labels,
    _check_line_colors,
    _get_line_collection,
    _get_line_xydata,
)

mpl = pytest.importorskip("matplotlib")
plt = pytest.importorskip("matplotlib.pyplot")

from pandas.plotting._matplotlib.core import LinePlot


def _frame(nrows=20, ncols=4):
    values = np.arange(nrows * ncols, dtype="float64").reshape(nrows, ncols)
    return DataFrame(values, columns=[f"c{i}" for i in range(ncols)])


def _draw(df, collection, monkeypatch, **kwargs):
    with monkeypatch.context() as ctx:
        if not collection:
            ctx.setattr(LinePlot, "_collection_data", lambda self: None)
        return df.plot(**kwargs)


def _state(ax):
    """What the drawn columns look like, whichever artist carried them."""
    return {
        "xydata": _get_line_xydata(ax),
        "xlim": ax.get_xlim(),
        "ylim": ax.get_ylim(),
        "xticks": ax.get_xticks().tolist(),
        "yticks": ax.get_yticks().tolist(),
    }


def _assert_same(left, right):
    assert len(left["xydata"]) == len(right["xydata"])
    for lft, rgt in zip(left["xydata"], right["xydata"], strict=True):
        tm.assert_numpy_array_equal(
            np.ma.filled(np.asarray(lft, dtype="float64"), np.nan),
            np.ma.filled(np.asarray(rgt, dtype="float64"), np.nan),
        )
    for key in ("xlim", "ylim", "xticks", "yticks"):
        assert left[key] == right[key], key


def test_collection_used_for_wide_short_frame():
    ax = _frame().plot()
    assert len(ax.collections) == 1
    # the columns are no longer individually addressable artists
    assert len(ax.get_lines()) == 0
    assert len(ax.collections[0].get_paths()) == 4


def test_collection_not_used_when_matplotlib_would_simplify():
    # matplotlib simplifies paths of 128 vertices or more, and a collection
    # would not get that, so the per-column draw stays
    assert len(_frame(nrows=127).plot().collections) == 1
    ax = _frame(nrows=128).plot()
    assert len(ax.collections) == 0
    assert len(ax.get_lines()) == 4


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"legend": True},
        {"legend": False},
        {"color": "red"},
        {"color": ["red", "green", "blue", "black"]},
        {"colormap": "viridis"},
        {"linewidth": 3},
        {"lw": 0.5},
        {"linestyle": "--"},
        {"ls": ":"},
        {"alpha": 0.4},
        {"logy": True},
        {"logx": True},
        {"use_index": False},
        {"xlim": (0, 5)},
        {"grid": True},
        {"title": "t", "xlabel": "xx", "ylabel": "yy"},
    ],
)
def test_collection_matches_per_column(kwargs, monkeypatch):
    df = _frame()
    collection = _draw(df, True, monkeypatch, **kwargs)
    per_column = _draw(df, False, monkeypatch, **kwargs)
    _assert_same(_state(collection), _state(per_column))


@pytest.mark.parametrize(
    "index",
    [
        None,
        pd.Index([3.5, 1.0, 2.25, 0.5, 9.0]),
        pd.Index(list("vwxyz")),
        pd.MultiIndex.from_product([["x", "y"], [1, 2, 3]])[:5],
    ],
)
def test_collection_matches_per_column_index(index, monkeypatch):
    df = _frame(nrows=5, ncols=3)
    if index is not None:
        df.index = index
    collection = _draw(df, True, monkeypatch)
    per_column = _draw(df, False, monkeypatch)
    _assert_same(_state(collection), _state(per_column))


def test_collection_keeps_nan_gaps(monkeypatch):
    df = DataFrame({"a": [1.0, np.nan, 3.0, 4.0], "b": [1.0, 2.0, np.nan, 4.0]})
    collection = _draw(df, True, monkeypatch)
    per_column = _draw(df, False, monkeypatch)
    _assert_same(_state(collection), _state(per_column))
    ydata = [xy[:, 1] for xy in _get_line_xydata(collection)]
    assert np.isnan(ydata[0][1])
    assert np.isnan(ydata[1][2])


def test_collection_colors_and_legend(monkeypatch):
    df = _frame()
    colors = ["red", "green", "blue", "black"]
    ax = df.plot(color=colors, legend=True)
    _check_line_colors(ax, colors)
    _check_legend_labels(ax, labels=list(df.columns))
    handles = ax.get_legend().legend_handles
    for handle, color in zip(handles, colors, strict=True):
        assert mpl.colors.to_rgba(handle.get_color()) == mpl.colors.to_rgba(color)


def test_collection_legend_avoids_the_data():
    # the collection has to stay visible to matplotlib's "best" placement, or
    # the legend lands on top of the lines
    df = DataFrame({f"c{j}": np.linspace(0, 10, 20) + j for j in range(4)})
    ax = df.plot(legend=True)
    ax.figure.canvas.draw()
    box = ax.get_legend().get_window_extent()
    # the data climbs to the right, so "best" must not choose the right side
    assert box.x0 < ax.get_window_extent().x0 + ax.get_window_extent().width / 2


@pytest.mark.parametrize("linestyle", ["-", "--"])
def test_collection_caps_and_joins_like_a_line(linestyle):
    # a collection defaults to the backend's cap/join rather than the lines.*
    # rcParams a Line2D uses, which shows at every vertex
    ax = _frame().plot(linestyle=linestyle)
    collection = ax.collections[0]
    line = mpl.lines.Line2D([], [], linestyle=linestyle)
    if linestyle == "-":
        assert collection.get_capstyle() == line.get_solid_capstyle()
        assert collection.get_joinstyle() == line.get_solid_joinstyle()
    else:
        assert collection.get_capstyle() == line.get_dash_capstyle()
        assert collection.get_joinstyle() == line.get_dash_joinstyle()


def test_collection_draws_above_the_grid():
    ax = _frame().plot(grid=True)
    assert ax.collections[0].get_zorder() == mpl.lines.Line2D.zorder


def test_collection_keeps_left_axis_visible_with_secondary_y():
    # GH#61532: the collection counts as plotted data, so adding a secondary-y
    # series must not hide the left axis
    _, ax = plt.subplots()
    _frame().plot(ax=ax)
    Series(np.arange(20.0), name="x").plot(legend=True, secondary_y=True, ax=ax)
    assert ax.get_yaxis().get_visible()
    assert ax.right_ax.get_yaxis().get_visible()


def test_secondary_y_keeps_left_axis_for_collection_plots():
    # the same fix applies to any plot whose data is a collection, not just to
    # wide line plots
    _, ax = plt.subplots()
    DataFrame({"x": np.arange(10.0), "y": np.arange(10.0)}).plot.scatter(
        x="x", y="y", ax=ax
    )
    Series(np.arange(10.0), name="s").plot(secondary_y=True, ax=ax)
    assert ax.get_yaxis().get_visible()


def test_collection_keeps_irregular_index_rotation(monkeypatch):
    # the tick rotation depends on nothing having been plotted yet, which a
    # collection has to report just as lines do
    df = _frame(nrows=5, ncols=3)
    df.index = timedelta_range("1 day", periods=5)
    collection = _draw(df, True, monkeypatch, x_compat=True)
    per_column = _draw(df, False, monkeypatch, x_compat=True)
    assert [t.get_rotation() for t in collection.get_xticklabels()] == [
        t.get_rotation() for t in per_column.get_xticklabels()
    ]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"subplots": True},
        {"stacked": True},
        {"secondary_y": True},
        {"secondary_y": ["c1"]},
        {"style": "--"},
        {"style": ["--", "-.", ":", "-"]},
        {"marker": "o"},
        {"yerr": 0.1},
        {"xerr": 0.1},
        {"kind": "area"},
        # a collection cannot reproduce these, and dropping them silently
        # would change the plot
        {"drawstyle": "steps-post"},
        {"zorder": 5},
    ],
)
def test_collection_not_used(kwargs):
    assert _frame().plot(**kwargs) is not None
    ax = plt.gcf().axes[0]
    assert _get_line_collection(ax) is None


def test_collection_known_limitations():
    # Documented in whatsnew: without per-column artists these matplotlib APIs
    # cannot see the drawn data. Pinned here so any change to them is a
    # deliberate one rather than a silent regression.
    ax = _frame().plot()

    # matplotlib does not support collections in relim
    ax.relim()
    assert not np.isfinite(ax.dataLim.get_points()).any()

    # so a bare ax.legend() finds no labelled artists
    handles, labels = ax.get_legend_handles_labels()
    assert handles == []
    assert labels == []


def test_collection_not_used_for_single_column():
    assert len(_frame(ncols=1).plot().collections) == 0
    assert len(Series(np.arange(20.0)).plot().collections) == 0


@pytest.mark.parametrize(
    "index",
    [
        date_range("2020-01-01", periods=5),
        pd.period_range("2020-01-01", periods=5, freq="D"),
    ],
)
def test_collection_not_used_for_ts_plot(index):
    df = _frame(nrows=5, ncols=3)
    df.index = index
    assert len(df.plot().collections) == 0


@pytest.mark.parametrize(
    "data",
    [
        {"a": date_range("2020-01-01", periods=3), "b": [1.0, 2.0, 3.0]},
        {"a": [1 + 2j, 2 + 1j, 3 + 0j], "b": [1.0, 2.0, 3.0]},
        {"a": [True, False, True], "b": [1.0, 2.0, 3.0]},
    ],
)
def test_collection_not_used_for_unstackable_dtypes(data):
    # these do not stack into one numeric array of y values
    df = DataFrame(data)
    warning = np.exceptions.ComplexWarning if df.dtypes.eq("complex128").any() else None
    with tm.assert_produces_warning(warning, check_stacklevel=False):
        ax = df.plot(include_bool=True)
    assert _get_line_collection(ax) is None


def test_collection_not_used_when_marker_comes_from_rcparams():
    # every column would carry a marker, which a collection cannot draw
    with mpl.rc_context({"lines.marker": "x"}):
        assert _get_line_collection(_frame().plot()) is None


@pytest.mark.parametrize(
    "prop_cycle",
    [
        {"color": ["r", "g"], "linestyle": ["-", "--"]},
        {"color": ["r", "g"], "marker": ["o", "s"]},
        {"linestyle": ["-", "--"]},
    ],
)
def test_collection_not_used_when_cycle_carries_more_than_color(prop_cycle):
    # the per-column draw would vary a property the collection cannot
    with mpl.rc_context({"axes.prop_cycle": plt.cycler(**prop_cycle)}):
        if "color" not in prop_cycle:
            pytest.skip("a colorless cycle is unsupported for other reasons")
        assert _get_line_collection(_frame().plot()) is None
    # a cycle we cannot read is skipped just the same, so pin that a color-only
    # one is still read as such rather than passing the check vacuously
    with mpl.rc_context({"axes.prop_cycle": plt.cycler(color=["r", "g"])}):
        assert _get_line_collection(_frame().plot()) is not None


def test_collection_follows_axes_prop_cycle(monkeypatch):
    # a cycle set on the axes rather than through rcParams still counts
    _, ax = plt.subplots()
    ax.set_prop_cycle(color=["r", "g"], linestyle=["-", "--"])
    assert _get_line_collection(_frame().plot(ax=ax)) is None

    _, ax = plt.subplots()
    ax.set_prop_cycle(color=["r", "g"])
    assert _get_line_collection(_frame().plot(ax=ax)) is not None


def test_collection_continues_axes_color_cycle(monkeypatch):
    # the colors have to pick up where the axes left off, as a per-column draw
    # without an explicit color would
    df = _frame(ncols=3)
    _, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    collection = _draw(df, True, monkeypatch, ax=ax)
    colors = collection.collections[0].get_edgecolor()

    _, ax2 = plt.subplots()
    ax2.plot([0, 1], [0, 1])
    per_column = _draw(df, False, monkeypatch, ax=ax2)
    # skip the line that was already on the axes
    expected = [
        mpl.colors.to_rgba(line.get_color()) for line in per_column.get_lines()[1:]
    ]
    assert [tuple(color) for color in colors] == expected
