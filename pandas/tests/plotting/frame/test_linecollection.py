"""
Wide numeric DataFrame line plots are drawn with a single LineCollection
instead of one Line2D per column, for speed (GH#61532).
"""

import numpy as np
import pytest

import pandas as pd

mpl = pytest.importorskip("matplotlib")
plt = pytest.importorskip("matplotlib.pyplot")
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba


def _n_linecollections(ax) -> int:
    return sum(isinstance(coll, LineCollection) for coll in ax.collections)


def test_wide_frame_uses_linecollection():
    # > 200 numeric columns -> single LineCollection, no per-column Line2D
    df = pd.DataFrame(np.random.default_rng(0).standard_normal((10, 201)))
    ax = df.plot(legend=False)
    assert _n_linecollections(ax) == 1
    assert len(ax.lines) == 0
    plt.close(ax.figure)


def test_narrow_frame_uses_line2d():
    # at the threshold we keep the per-column Line2D path
    df = pd.DataFrame(np.random.default_rng(0).standard_normal((10, 200)))
    ax = df.plot(legend=False)
    assert _n_linecollections(ax) == 0
    assert len(ax.lines) == 200
    plt.close(ax.figure)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"marker": "o"},
        {"style": "-"},
        {"stacked": True},
        {"secondary_y": [0]},
        {"subplots": True},
    ],
)
def test_excluded_cases_keep_line2d(kwargs):
    # cases a LineCollection cannot reproduce fall back to the per-column path
    # (all-positive data so the stacked case is valid)
    df = pd.DataFrame(np.random.default_rng(0).random((10, 201)))
    ax = df.plot(legend=False, **kwargs)
    assert all(_n_linecollections(sub) == 0 for sub in np.ravel(ax))
    plt.close("all")


def test_datetime_index_keeps_line2d():
    # datetime x-axis is handled by the time-series path, not LineCollection
    idx = pd.date_range("2020-01-01", periods=10)
    df = pd.DataFrame(np.random.default_rng(0).standard_normal((10, 201)), index=idx)
    ax = df.plot(legend=False)
    assert _n_linecollections(ax) == 0
    plt.close(ax.figure)


def test_linecollection_segments_and_colors():
    df = pd.DataFrame(np.random.default_rng(0).standard_normal((10, 201)))
    ax = df.plot()
    lc = next(c for c in ax.collections if isinstance(c, LineCollection))

    segments = lc.get_segments()
    assert len(segments) == 201
    assert all(seg.shape == (10, 2) for seg in segments)

    # colors follow the default property cycle, as with per-column plotting
    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    colors = lc.get_colors()
    assert tuple(colors[0]) == to_rgba(cycle[0])
    assert tuple(colors[1]) == to_rgba(cycle[1 % len(cycle)])

    # every column is still represented in the legend
    assert len(ax.get_legend().get_texts()) == 201
    plt.close(ax.figure)


def test_linecollection_forwards_line_styling():
    # alpha / linewidth / linestyle passed by the user are reproduced
    df = pd.DataFrame(np.random.default_rng(0).standard_normal((10, 201)))
    ax = df.plot(legend=False, alpha=0.4, linewidth=3, linestyle="--")
    lc = next(c for c in ax.collections if isinstance(c, LineCollection))
    assert lc.get_alpha() == 0.4
    assert all(lw == 3 for lw in lc.get_linewidth())
    # matplotlib normalizes "--" to its dashed on/off sequence
    assert lc.get_linestyle()[0] != (0.0, None)
    plt.close(ax.figure)


def test_linecollection_matches_line2d_limits():
    # axis limits are identical whether or not the fast path is taken
    df = pd.DataFrame(np.random.default_rng(1).standard_normal((30, 201)).cumsum(0))

    ax_fast = df.plot(legend=False)

    # force the per-column path on the same data
    from pandas.plotting._matplotlib.core import LinePlot

    orig = LinePlot._use_collection
    LinePlot._use_collection = lambda self: False
    try:
        ax_slow = df.plot(legend=False)
    finally:
        LinePlot._use_collection = orig

    assert ax_fast.get_xlim() == pytest.approx(ax_slow.get_xlim())
    assert ax_fast.get_ylim() == pytest.approx(ax_slow.get_ylim())
    plt.close("all")
