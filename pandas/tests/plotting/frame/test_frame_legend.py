from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import pytest

import pandas.util._test_decorators as td

from pandas import DataFrame


@td.skip_if_no_mpl
@pytest.mark.xfail(
    reason=(
        "Open bug in matplotlib "
        "https://github.com/matplotlib/matplotlib/issues/11357"
    )
)
def test_mixed_yerr():
    # https://github.com/pandas-dev/pandas/issues/39522
    df = DataFrame([{"x": 1, "a": 1, "b": 1}, {"x": 2, "a": 2, "b": 3}])

    ax = df.plot("x", "a", c="orange", yerr=0.1, label="orange")
    df.plot("x", "b", c="blue", yerr=None, ax=ax, label="blue")

    legend = ax.get_legend()
    result_handles = legend.legendHandles

    assert isinstance(result_handles[0], LineCollection)
    assert isinstance(result_handles[1], Line2D)


@td.skip_if_no_mpl
def test_legend_false():
    # https://github.com/pandas-dev/pandas/issues/40044
    df = DataFrame({"a": [1, 1], "b": [2, 3]})
    df2 = DataFrame({"d": [2.5, 2.5]})

    ax = df.plot(legend=True, color={"a": "blue", "b": "green"}, secondary_y="b")
    df2.plot(legend=True, color={"d": "red"}, ax=ax)
    legend = ax.get_legend()
    result = [handle.get_color() for handle in legend.legendHandles]
    expected = ["blue", "green", "red"]
    assert result == expected
