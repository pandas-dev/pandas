from matplotlib.container import ErrorbarContainer
from matplotlib.lines import Line2D

from pandas import DataFrame


def test_mixed_yerr():
    # https://github.com/pandas-dev/pandas/issues/39522

    df = DataFrame([{"x": 1, "a": 1, "b": 1}, {"x": 2, "a": 2, "b": 3}])

    ax = df.plot("x", "a", c="orange", yerr=0.1, label="orange")
    df.plot("x", "b", c="blue", yerr=None, ax=ax, label="blue")

    result_handles, result_labels = ax.get_legend_handles_labels()

    assert isinstance(result_handles[0], Line2D)
    assert isinstance(result_handles[1], ErrorbarContainer)

    expected_labels = ["blue", "orange"]
    assert result_labels == expected_labels


def test_all_have_yerr():
    # https://github.com/pandas-dev/pandas/issues/39522

    df = DataFrame([{"x": 1, "a": 1, "b": 1}, {"x": 2, "a": 2, "b": 3}])

    ax = df.plot("x", "a", c="orange", yerr=0.1, label="orange")
    df.plot("x", "b", c="blue", yerr=0.1, ax=ax, label="blue")

    result_handles, result_labels = ax.get_legend_handles_labels()

    assert isinstance(result_handles[0], ErrorbarContainer)
    assert isinstance(result_handles[1], ErrorbarContainer)

    expected_labels = ["orange", "blue"]
    assert result_labels == expected_labels


def test_none_have_yerr():
    # https://github.com/pandas-dev/pandas/issues/39522

    df = DataFrame([{"x": 1, "a": 1, "b": 1}, {"x": 2, "a": 2, "b": 3}])

    ax = df.plot("x", "a", c="orange", label="orange")
    df.plot("x", "b", c="blue", ax=ax, label="blue")

    result_handles, result_labels = ax.get_legend_handles_labels()

    assert isinstance(result_handles[0], Line2D)
    assert isinstance(result_handles[1], Line2D)

    expected_labels = ["orange", "blue"]
    assert result_labels == expected_labels


def test_legend_false():
    # https://github.com/pandas-dev/pandas/issues/40044

    df = DataFrame({"a": [1, 1], "b": [2, 2]})
    df2 = DataFrame({"d": [2.5, 2.5]})

    ax = df.plot(legend=True, color={"a": "blue", "b": "green"}, secondary_y="b")
    df2.plot(legend=False, color="red", ax=ax)

    handles, result_labels = ax.get_legend_handles_labels()
    result_colors = [i.get_color() for i in handles]
    expected_labels = ["a", "d"]
    expected_colors = ["blue", "red"]

    assert result_labels == expected_labels
    assert result_colors == expected_colors
