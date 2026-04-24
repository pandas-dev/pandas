"""Test cases for misc plot functions"""

import os

import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import (
    DataFrame,
    Index,
    Series,
    Timestamp,
    date_range,
    interval_range,
    period_range,
    plotting,
    read_csv,
)
import pandas._testing as tm
from pandas.tests.plotting.common import (
    _check_colors,
    _check_legend_labels,
    _check_plot_works,
    _check_text_labels,
    _check_ticks_props,
)

mpl = pytest.importorskip("matplotlib")
plt = pytest.importorskip("matplotlib.pyplot")
cm = pytest.importorskip("matplotlib.cm")

import re

from pandas.plotting._matplotlib.style import get_standard_colors


@pytest.fixture
def iris(datapath) -> DataFrame:
    """
    The iris dataset as a DataFrame.
    """
    return read_csv(datapath("io", "data", "csv", "iris.csv"))


@td.skip_if_installed("matplotlib")
def test_import_error_message():
    # GH-19810
    df = DataFrame({"A": [1, 2]})

    with pytest.raises(ImportError, match="matplotlib is required for plotting"):
        df.plot()


def test_get_accessor_args():
    func = plotting._core.PlotAccessor._get_call_args

    msg = "Called plot accessor for type list, expected Series or DataFrame"
    with pytest.raises(TypeError, match=msg):
        func(backend_name="", data=[], args=[], kwargs={})

    msg = "should not be called with positional arguments"
    with pytest.raises(TypeError, match=msg):
        func(backend_name="", data=Series(dtype=object), args=["line", None], kwargs={})

    x, y, kind, kwargs = func(
        backend_name="",
        data=DataFrame(),
        args=["x"],
        kwargs={"y": "y", "kind": "bar", "grid": False},
    )
    assert x == "x"
    assert y == "y"
    assert kind == "bar"
    assert kwargs == {"grid": False}

    x, y, kind, kwargs = func(
        backend_name="pandas.plotting._matplotlib",
        data=Series(dtype=object),
        args=[],
        kwargs={},
    )
    assert x is None
    assert y is None
    assert kind == "line"
    assert len(kwargs) == 24


@pytest.mark.parametrize("kind", plotting.PlotAccessor._all_kinds)
@pytest.mark.parametrize(
    "data", [DataFrame(np.arange(15).reshape(5, 3)), Series(range(5))]
)
@pytest.mark.parametrize(
    "index",
    [
        Index(range(5)),
        date_range("2020-01-01", periods=5),
        period_range("2020-01-01", periods=5),
    ],
)
def test_savefig(kind, data, index):
    fig, ax = plt.subplots()
    data.index = index
    kwargs = {}
    if kind in ["hexbin", "scatter", "pie"]:
        if isinstance(data, Series):
            pytest.skip(f"{kind} not supported with Series")
        kwargs = {"x": 0, "y": 1}
    data.plot(kind=kind, ax=ax, **kwargs)
    fig.savefig(os.devnull)


class TestSeriesPlots:
    def test_autocorrelation_plot(self):
        ser = Series(
            np.arange(10, dtype=np.float64),
            index=date_range("2020-01-01", periods=10),
            name="ts",
        )
        # Ensure no UserWarning when making plot
        with tm.assert_produces_warning(None):
            _check_plot_works(plotting.autocorrelation_plot, series=ser)
            _check_plot_works(plotting.autocorrelation_plot, series=ser.values)

            ax = plotting.autocorrelation_plot(ser, label="Test")
        _check_legend_labels(ax, labels=["Test"])

    @pytest.mark.parametrize("kwargs", [{}, {"lag": 5}])
    def test_lag_plot(self, kwargs):
        ser = Series(
            np.arange(10, dtype=np.float64),
            index=date_range("2020-01-01", periods=10),
            name="ts",
        )
        _check_plot_works(plotting.lag_plot, series=ser, **kwargs)

    def test_bootstrap_plot(self):
        ser = Series(
            np.arange(10, dtype=np.float64),
            index=date_range("2020-01-01", periods=10),
            name="ts",
        )
        _check_plot_works(plotting.bootstrap_plot, series=ser, size=10)


class TestDataFramePlots:
    @pytest.mark.parametrize("pass_axis", [False, True])
    def test_scatter_matrix_axis(self, pass_axis):
        pytest.importorskip("scipy")
        scatter_matrix = plotting.scatter_matrix

        ax = None
        if pass_axis:
            _, ax = mpl.pyplot.subplots(3, 3)

        df = DataFrame(np.random.default_rng(2).standard_normal((10, 3)))

        # we are plotting multiples on a sub-plot
        with tm.assert_produces_warning(UserWarning, check_stacklevel=False):
            axes = _check_plot_works(
                scatter_matrix,
                frame=df,
                range_padding=0.1,
                ax=ax,
            )
        axes0_labels = axes[0][0].yaxis.get_majorticklabels()
        # GH 5662
        expected = ["-2", "-1", "0"]
        _check_text_labels(axes0_labels, expected)
        _check_ticks_props(axes, xlabelsize=8, xrot=90, ylabelsize=8, yrot=0)

    @pytest.mark.parametrize("pass_axis", [False, True])
    def test_scatter_matrix_axis_smaller(self, pass_axis):
        pytest.importorskip("scipy")
        scatter_matrix = plotting.scatter_matrix

        ax = None
        if pass_axis:
            _, ax = mpl.pyplot.subplots(3, 3)

        df = DataFrame(np.random.default_rng(11).standard_normal((10, 3)))
        df[0] = (df[0] - 2) / 3

        # we are plotting multiples on a sub-plot
        with tm.assert_produces_warning(UserWarning, check_stacklevel=False):
            axes = _check_plot_works(
                scatter_matrix,
                frame=df,
                range_padding=0.1,
                ax=ax,
            )
        axes0_labels = axes[0][0].yaxis.get_majorticklabels()
        expected = ["-1.25", "-1.0", "-0.75", "-0.5"]
        _check_text_labels(axes0_labels, expected)
        _check_ticks_props(axes, xlabelsize=8, xrot=90, ylabelsize=8, yrot=0)

    @pytest.mark.slow
    def test_andrews_curves_no_warning(self, iris):
        # Ensure no UserWarning when making plot
        with tm.assert_produces_warning(None):
            _check_plot_works(plotting.andrews_curves, frame=iris, class_column="Name")

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "linecolors",
        [
            ("#556270", "#4ECDC4", "#C7F464"),
            ["dodgerblue", "aquamarine", "seagreen"],
        ],
    )
    @pytest.mark.parametrize(
        "df",
        [
            "iris",
            DataFrame(
                {
                    "A": np.random.default_rng(2).standard_normal(10),
                    "B": np.random.default_rng(2).standard_normal(10),
                    "C": np.random.default_rng(2).standard_normal(10),
                    "Name": ["A"] * 10,
                }
            ),
        ],
    )
    def test_andrews_curves_linecolors(self, request, df, linecolors):
        if isinstance(df, str):
            df = request.getfixturevalue(df)
        ax = _check_plot_works(
            plotting.andrews_curves, frame=df, class_column="Name", color=linecolors
        )
        _check_colors(
            ax.get_lines()[:10], linecolors=linecolors, mapping=df["Name"][:10]
        )

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "df",
        [
            "iris",
            DataFrame(
                {
                    "A": np.random.default_rng(2).standard_normal(10),
                    "B": np.random.default_rng(2).standard_normal(10),
                    "C": np.random.default_rng(2).standard_normal(10),
                    "Name": ["A"] * 10,
                }
            ),
        ],
    )
    def test_andrews_curves_cmap(self, request, df):
        if isinstance(df, str):
            df = request.getfixturevalue(df)
        cmaps = [cm.jet(n) for n in np.linspace(0, 1, df["Name"].nunique())]
        ax = _check_plot_works(
            plotting.andrews_curves, frame=df, class_column="Name", color=cmaps
        )
        _check_colors(ax.get_lines()[:10], linecolors=cmaps, mapping=df["Name"][:10])

    @pytest.mark.slow
    def test_andrews_curves_handle(self):
        colors = ["b", "g", "r"]
        df = DataFrame({"A": [1, 2, 3], "B": [1, 2, 3], "C": [1, 2, 3], "Name": colors})
        ax = plotting.andrews_curves(df, "Name", color=colors)
        handles, _ = ax.get_legend_handles_labels()
        _check_colors(handles, linecolors=colors)

    @pytest.mark.slow
    @pytest.mark.parametrize(
        "color",
        [("#556270", "#4ECDC4", "#C7F464"), ["dodgerblue", "aquamarine", "seagreen"]],
    )
    def test_parallel_coordinates_colors(self, iris, color):
        df = iris

        ax = _check_plot_works(
            plotting.parallel_coordinates, frame=df, class_column="Name", color=color
        )
        _check_colors(ax.get_lines()[:10], linecolors=color, mapping=df["Name"][:10])

    @pytest.mark.slow
    def test_parallel_coordinates_cmap(self, iris):
        df = iris

        ax = _check_plot_works(
            plotting.parallel_coordinates,
            frame=df,
            class_column="Name",
            colormap=cm.jet,
        )
        cmaps = [mpl.cm.jet(n) for n in np.linspace(0, 1, df["Name"].nunique())]
        _check_colors(ax.get_lines()[:10], linecolors=cmaps, mapping=df["Name"][:10])

    @pytest.mark.slow
    def test_parallel_coordinates_line_diff(self, iris):
        df = iris

        ax = _check_plot_works(
            plotting.parallel_coordinates, frame=df, class_column="Name"
        )
        nlines = len(ax.get_lines())
        nxticks = len(ax.xaxis.get_ticklabels())

        ax = _check_plot_works(
            plotting.parallel_coordinates, frame=df, class_column="Name", axvlines=False
        )
        assert len(ax.get_lines()) == (nlines - nxticks)

    @pytest.mark.slow
    def test_parallel_coordinates_handles(self, iris):
        df = iris
        colors = ["b", "g", "r"]
        df = DataFrame({"A": [1, 2, 3], "B": [1, 2, 3], "C": [1, 2, 3], "Name": colors})
        ax = plotting.parallel_coordinates(df, "Name", color=colors)
        handles, _ = ax.get_legend_handles_labels()
        _check_colors(handles, linecolors=colors)

    # not sure if this is indicative of a problem
    @pytest.mark.filterwarnings("ignore:Attempting to set:UserWarning")
    def test_parallel_coordinates_with_sorted_labels(self):
        # GH 15908
        df = DataFrame(
            {
                "feat": list(range(30)),
                "class": [2 for _ in range(10)]
                + [3 for _ in range(10)]
                + [1 for _ in range(10)],
            }
        )
        ax = plotting.parallel_coordinates(df, "class", sort_labels=True)
        polylines, labels = ax.get_legend_handles_labels()
        color_label_tuples = zip(
            [polyline.get_color() for polyline in polylines], labels, strict=True
        )
        ordered_color_label_tuples = sorted(color_label_tuples, key=lambda x: x[1])
        prev_next_tupels = zip(
            list(ordered_color_label_tuples[0:-1]),
            list(ordered_color_label_tuples[1:]),
            strict=True,
        )
        for prev, nxt in prev_next_tupels:
            # labels and colors are ordered strictly increasing
            assert prev[1] < nxt[1] and prev[0] < nxt[0]

    def test_radviz_no_warning(self, iris):
        # Ensure no UserWarning when making plot
        with tm.assert_produces_warning(None):
            _check_plot_works(plotting.radviz, frame=iris, class_column="Name")

    @pytest.mark.parametrize(
        "color",
        [("#556270", "#4ECDC4", "#C7F464"), ["dodgerblue", "aquamarine", "seagreen"]],
    )
    def test_radviz_color(self, iris, color):
        df = iris
        ax = _check_plot_works(
            plotting.radviz, frame=df, class_column="Name", color=color
        )
        # skip Circle drawn as ticks
        patches = [p for p in ax.patches[:20] if p.get_label() != ""]
        _check_colors(patches[:10], facecolors=color, mapping=df["Name"][:10])

    def test_radviz_color_cmap(self, iris):
        df = iris
        ax = _check_plot_works(
            plotting.radviz, frame=df, class_column="Name", colormap=cm.jet
        )
        cmaps = [mpl.cm.jet(n) for n in np.linspace(0, 1, df["Name"].nunique())]
        patches = [p for p in ax.patches[:20] if p.get_label() != ""]
        _check_colors(patches, facecolors=cmaps, mapping=df["Name"][:10])

    def test_radviz_colors_handles(self):
        colors = [[0.0, 0.0, 1.0, 1.0], [0.0, 0.5, 1.0, 1.0], [1.0, 0.0, 0.0, 1.0]]
        df = DataFrame(
            {"A": [1, 2, 3], "B": [2, 1, 3], "C": [3, 2, 1], "Name": ["b", "g", "r"]}
        )
        ax = plotting.radviz(df, "Name", color=colors)
        handles, _ = ax.get_legend_handles_labels()
        _check_colors(handles, facecolors=colors)

    def test_subplot_titles(self, iris):
        df = iris.drop("Name", axis=1).head()
        # Use the column names as the subplot titles
        title = list(df.columns)

        # Case len(title) == len(df)
        plot = df.plot(subplots=True, title=title)
        assert [p.get_title() for p in plot] == title

    def test_subplot_titles_too_much(self, iris):
        df = iris.drop("Name", axis=1).head()
        # Use the column names as the subplot titles
        title = list(df.columns)
        # Case len(title) > len(df)
        msg = (
            "The length of `title` must equal the number of columns if "
            "using `title` of type `list` and `subplots=True`"
        )
        with pytest.raises(ValueError, match=msg):
            df.plot(subplots=True, title=[*title, "kittens > puppies"])

    def test_subplot_titles_too_little(self, iris):
        df = iris.drop("Name", axis=1).head()
        # Use the column names as the subplot titles
        title = list(df.columns)
        msg = (
            "The length of `title` must equal the number of columns if "
            "using `title` of type `list` and `subplots=True`"
        )
        # Case len(title) < len(df)
        with pytest.raises(ValueError, match=msg):
            df.plot(subplots=True, title=title[:2])

    def test_subplot_titles_subplots_false(self, iris):
        df = iris.drop("Name", axis=1).head()
        # Use the column names as the subplot titles
        title = list(df.columns)
        # Case subplots=False and title is of type list
        msg = (
            "Using `title` of type `list` is not supported unless "
            "`subplots=True` is passed"
        )
        with pytest.raises(ValueError, match=msg):
            df.plot(subplots=False, title=title)

    def test_subplot_titles_numeric_square_layout(self, iris):
        df = iris.drop("Name", axis=1).head()
        # Use the column names as the subplot titles
        title = list(df.columns)
        # Case df with 3 numeric columns but layout of (2,2)
        plot = df.drop("SepalWidth", axis=1).plot(
            subplots=True, layout=(2, 2), title=title[:-1]
        )
        title_list = [ax.get_title() for sublist in plot for ax in sublist]
        assert title_list == [*title[:3], ""]

    def test_get_standard_colors_random_seed(self):
        # GH17525
        df = DataFrame(np.zeros((10, 10)))

        # Make sure that the random seed isn't reset by get_standard_colors
        plotting.parallel_coordinates(df, 0)
        rand1 = np.random.default_rng(None).random()
        plotting.parallel_coordinates(df, 0)
        rand2 = np.random.default_rng(None).random()
        assert rand1 != rand2

    def test_get_standard_colors_consistency(self):
        # GH17525
        # Make sure it produces the same colors every time it's called
        color1 = get_standard_colors(1, color_type="random")
        color2 = get_standard_colors(1, color_type="random")
        assert color1 == color2

    def test_get_standard_colors_default_num_colors(self):
        # Make sure the default color_types returns the specified amount
        color1 = get_standard_colors(1, color_type="default")
        color2 = get_standard_colors(9, color_type="default")
        color3 = get_standard_colors(20, color_type="default")
        assert len(color1) == 1
        assert len(color2) == 9
        assert len(color3) == 20

    def test_plot_single_color(self):
        # Example from #20585. All 3 bars should have the same color
        df = DataFrame(
            {
                "account-start": ["2017-02-03", "2017-03-03", "2017-01-01"],
                "client": ["Alice Anders", "Bob Baker", "Charlie Chaplin"],
                "balance": [-1432.32, 10.43, 30000.00],
                "db-id": [1234, 2424, 251],
                "proxy-id": [525, 1525, 2542],
                "rank": [52, 525, 32],
            }
        )
        ax = df.client.value_counts().plot.bar()
        colors = [rect.get_facecolor() for rect in ax.get_children()[0:3]]
        assert all(color == colors[0] for color in colors)

    def test_get_standard_colors_no_appending(self):
        # GH20726

        # Make sure not to add more colors so that matplotlib can cycle
        # correctly.
        color_before = mpl.cm.gnuplot(range(5))
        color_after = get_standard_colors(1, color=color_before)
        assert len(color_after) == len(color_before)

        df = DataFrame(
            np.random.default_rng(2).standard_normal((48, 4)), columns=list("ABCD")
        )

        color_list = mpl.cm.gnuplot(np.linspace(0, 1, 16))
        p = df.A.plot.bar(figsize=(16, 7), color=color_list)
        assert p.patches[1].get_facecolor() == p.patches[17].get_facecolor()

    @pytest.mark.parametrize("kind", ["bar", "line"])
    def test_dictionary_color(self, kind):
        # issue-8193
        # Test plot color dictionary format
        data_files = ["a", "b"]

        expected = [(0.5, 0.24, 0.6), (0.3, 0.7, 0.7)]

        df1 = DataFrame(np.random.default_rng(2).random((2, 2)), columns=data_files)
        dic_color = {"b": (0.3, 0.7, 0.7), "a": (0.5, 0.24, 0.6)}

        ax = df1.plot(kind=kind, color=dic_color)
        if kind == "bar":
            colors = [rect.get_facecolor()[0:-1] for rect in ax.get_children()[0:3:2]]
        else:
            colors = [rect.get_color() for rect in ax.get_lines()[0:2]]
        assert all(color == expected[index] for index, color in enumerate(colors))

    def test_bar_plot(self):
        # GH38947
        # Test bar plot with string and int index
        expected = [mpl.text.Text(0, 0, "0"), mpl.text.Text(1, 0, "Total")]

        df = DataFrame(
            {
                "a": [1, 2],
            },
            index=Index([0, "Total"]),
        )
        plot_bar = df.plot.bar()
        assert all(
            (a.get_text() == b.get_text())
            for a, b in zip(plot_bar.get_xticklabels(), expected, strict=True)
        )

    def test_barh_plot_labels_mixed_integer_string(self):
        # GH39126
        # Test barh plot with string and integer at the same column
        df = DataFrame([{"word": 1, "value": 0}, {"word": "knowledge", "value": 2}])
        plot_barh = df.plot.barh(x="word", legend=None)
        expected_yticklabels = [
            mpl.text.Text(0, 0, "1"),
            mpl.text.Text(0, 1, "knowledge"),
        ]
        assert all(
            actual.get_text() == expected.get_text()
            for actual, expected in zip(
                plot_barh.get_yticklabels(), expected_yticklabels, strict=True
            )
        )

    def test_has_externally_shared_axis_x_axis(self):
        # GH33819
        # Test _has_externally_shared_axis() works for x-axis
        func = plotting._matplotlib.tools._has_externally_shared_axis

        fig = mpl.pyplot.figure()
        plots = fig.subplots(2, 4)

        # Create *externally* shared axes for first and third columns
        plots[0][0] = fig.add_subplot(231, sharex=plots[1][0])
        plots[0][2] = fig.add_subplot(233, sharex=plots[1][2])

        # Create *internally* shared axes for second and third columns
        plots[0][1].twinx()
        plots[0][2].twinx()

        # First  column is only externally shared
        # Second column is only internally shared
        # Third  column is both
        # Fourth column is neither
        assert func(plots[0][0], "x")
        assert not func(plots[0][1], "x")
        assert func(plots[0][2], "x")
        assert not func(plots[0][3], "x")

    def test_has_externally_shared_axis_y_axis(self):
        # GH33819
        # Test _has_externally_shared_axis() works for y-axis
        func = plotting._matplotlib.tools._has_externally_shared_axis

        fig = mpl.pyplot.figure()
        plots = fig.subplots(4, 2)

        # Create *externally* shared axes for first and third rows
        plots[0][0] = fig.add_subplot(321, sharey=plots[0][1])
        plots[2][0] = fig.add_subplot(325, sharey=plots[2][1])

        # Create *internally* shared axes for second and third rows
        plots[1][0].twiny()
        plots[2][0].twiny()

        # First  row is only externally shared
        # Second row is only internally shared
        # Third  row is both
        # Fourth row is neither
        assert func(plots[0][0], "y")
        assert not func(plots[1][0], "y")
        assert func(plots[2][0], "y")
        assert not func(plots[3][0], "y")

    def test_has_externally_shared_axis_invalid_compare_axis(self):
        # GH33819
        # Test _has_externally_shared_axis() raises an exception when
        # passed an invalid value as compare_axis parameter
        func = plotting._matplotlib.tools._has_externally_shared_axis

        fig = mpl.pyplot.figure()
        plots = fig.subplots(4, 2)

        # Create arbitrary axes
        plots[0][0] = fig.add_subplot(321, sharey=plots[0][1])

        # Check that an invalid compare_axis value triggers the expected exception
        msg = "needs 'x' or 'y' as a second parameter"
        with pytest.raises(ValueError, match=msg):
            func(plots[0][0], "z")

    def test_externally_shared_axes(self):
        # Example from GH33819
        # Create data
        df = DataFrame(
            {
                "a": np.random.default_rng(2).standard_normal(10),
                "b": np.random.default_rng(2).standard_normal(10),
            }
        )

        # Create figure
        fig = mpl.pyplot.figure()
        plots = fig.subplots(2, 3)

        # Create *externally* shared axes
        plots[0][0] = fig.add_subplot(231, sharex=plots[1][0])
        # note: no plots[0][1] that's the twin only case
        plots[0][2] = fig.add_subplot(233, sharex=plots[1][2])

        # Create *internally* shared axes
        # note: no plots[0][0] that's the external only case
        twin_ax1 = plots[0][1].twinx()
        twin_ax2 = plots[0][2].twinx()

        # Plot data to primary axes
        df["a"].plot(ax=plots[0][0], title="External share only").set_xlabel(
            "this label should never be visible"
        )
        df["a"].plot(ax=plots[1][0])

        df["a"].plot(ax=plots[0][1], title="Internal share (twin) only").set_xlabel(
            "this label should always be visible"
        )
        df["a"].plot(ax=plots[1][1])

        df["a"].plot(ax=plots[0][2], title="Both").set_xlabel(
            "this label should never be visible"
        )
        df["a"].plot(ax=plots[1][2])

        # Plot data to twinned axes
        df["b"].plot(ax=twin_ax1, color="green")
        df["b"].plot(ax=twin_ax2, color="yellow")

        assert not plots[0][0].xaxis.get_label().get_visible()
        assert plots[0][1].xaxis.get_label().get_visible()
        assert not plots[0][2].xaxis.get_label().get_visible()

    def test_plot_bar_axis_units_timestamp_conversion(self):
        # GH 38736
        # Ensure string x-axis from the second plot will not be converted to datetime
        # due to axis data from first plot
        df = DataFrame(
            [1.0],
            index=[Timestamp("2022-02-22 22:22:22")],
        )
        _check_plot_works(df.plot)
        s = Series({"A": 1.0})
        _check_plot_works(s.plot.bar)

    def test_bar_plt_xaxis_intervalrange(self):
        # GH 38969
        # Ensure IntervalIndex x-axis produces a bar plot as expected
        expected = [mpl.text.Text(0, 0, "([0, 1],)"), mpl.text.Text(1, 0, "([1, 2],)")]
        s = Series(
            [1, 2],
            index=[interval_range(0, 2, closed="both")],
        )
        _check_plot_works(s.plot.bar)
        assert all(
            (a.get_text() == b.get_text())
            for a, b in zip(s.plot.bar().get_xticklabels(), expected, strict=True)
        )


@pytest.fixture
def df_bar_data():
    return np.random.default_rng(3).integers(0, 100, 5)


@pytest.fixture
def df_bar_df(df_bar_data) -> DataFrame:
    df_bar_df = DataFrame(
        {
            "A": df_bar_data,
            "B": df_bar_data[::-1],
            "C": df_bar_data[0],
            "D": df_bar_data[-1],
        }
    )
    return df_bar_df


def _df_bar_xyheight_from_ax_helper(df_bar_data, ax, subplot_division):
    subplot_data_df_list = []

    # get xy and height of squares representing data, separated by subplots
    for i in range(len(subplot_division)):
        subplot_data = np.array(
            [
                (x.get_x(), x.get_y(), x.get_height())
                for x in ax[i].findobj(plt.Rectangle)
                if x.get_height() in df_bar_data
            ]
        )
        subplot_data_df_list.append(
            DataFrame(data=subplot_data, columns=["x_coord", "y_coord", "height"])
        )

    return subplot_data_df_list


def _df_bar_subplot_checker(df_bar_data, df_bar_df, subplot_data_df, subplot_columns):
    subplot_sliced_by_source = [
        subplot_data_df.iloc[
            len(df_bar_data) * i : len(df_bar_data) * (i + 1)
        ].reset_index()
        for i in range(len(subplot_columns))
    ]

    if len(subplot_columns) == 1:
        expected_total_height = df_bar_df.loc[:, subplot_columns[0]]
    else:
        expected_total_height = df_bar_df.loc[:, subplot_columns].sum(axis=1)

    for i in range(len(subplot_columns)):
        sliced_df = subplot_sliced_by_source[i]
        if i == 0:
            # Checks that the bar chart starts y=0
            assert (sliced_df["y_coord"] == 0).all()
            height_iter = sliced_df["y_coord"].add(sliced_df["height"])
        else:
            height_iter = height_iter + sliced_df["height"]

        if i + 1 == len(subplot_columns):
            # Checks final height matches what is expected
            tm.assert_series_equal(
                height_iter, expected_total_height, check_names=False, check_dtype=False
            )
        else:
            # Checks each preceding bar ends where the next one starts
            next_start_coord = subplot_sliced_by_source[i + 1]["y_coord"]
            tm.assert_series_equal(
                height_iter, next_start_coord, check_names=False, check_dtype=False
            )


# GH Issue 61018
@pytest.mark.parametrize("columns_used", [["A", "B"], ["C", "D"], ["D", "A"]])
def test_bar_1_subplot_1_double_stacked(df_bar_data, df_bar_df, columns_used):
    df_bar_df_trimmed = df_bar_df[columns_used]
    subplot_division = [columns_used]
    ax = df_bar_df_trimmed.plot(subplots=subplot_division, kind="bar", stacked=True)
    subplot_data_df_list = _df_bar_xyheight_from_ax_helper(
        df_bar_data, ax, subplot_division
    )
    for i in range(len(subplot_data_df_list)):
        _df_bar_subplot_checker(
            df_bar_data, df_bar_df_trimmed, subplot_data_df_list[i], subplot_division[i]
        )


@pytest.mark.parametrize(
    "columns_used", [["A", "B", "C"], ["A", "C", "B"], ["D", "A", "C"]]
)
def test_bar_2_subplot_1_double_stacked(df_bar_data, df_bar_df, columns_used):
    df_bar_df_trimmed = df_bar_df[columns_used]
    subplot_division = [(columns_used[0], columns_used[1]), (columns_used[2],)]
    ax = df_bar_df_trimmed.plot(subplots=subplot_division, kind="bar", stacked=True)
    subplot_data_df_list = _df_bar_xyheight_from_ax_helper(
        df_bar_data, ax, subplot_division
    )
    for i in range(len(subplot_data_df_list)):
        _df_bar_subplot_checker(
            df_bar_data, df_bar_df_trimmed, subplot_data_df_list[i], subplot_division[i]
        )


@pytest.mark.parametrize(
    "subplot_division",
    [
        [("A", "B"), ("C", "D")],
        [("A", "D"), ("C", "B")],
        [("B", "C"), ("D", "A")],
        [("B", "D"), ("C", "A")],
    ],
)
def test_bar_2_subplot_2_double_stacked(df_bar_data, df_bar_df, subplot_division):
    ax = df_bar_df.plot(subplots=subplot_division, kind="bar", stacked=True)
    subplot_data_df_list = _df_bar_xyheight_from_ax_helper(
        df_bar_data, ax, subplot_division
    )
    for i in range(len(subplot_data_df_list)):
        _df_bar_subplot_checker(
            df_bar_data, df_bar_df, subplot_data_df_list[i], subplot_division[i]
        )


@pytest.mark.parametrize(
    "subplot_division",
    [[("A", "B", "C")], [("A", "D", "B")], [("C", "A", "D")], [("D", "C", "A")]],
)
def test_bar_2_subplots_1_triple_stacked(df_bar_data, df_bar_df, subplot_division):
    ax = df_bar_df.plot(subplots=subplot_division, kind="bar", stacked=True)
    subplot_data_df_list = _df_bar_xyheight_from_ax_helper(
        df_bar_data, ax, subplot_division
    )
    for i in range(len(subplot_data_df_list)):
        _df_bar_subplot_checker(
            df_bar_data, df_bar_df, subplot_data_df_list[i], subplot_division[i]
        )


def test_bar_subplots_stacking_bool(df_bar_data, df_bar_df):
    subplot_division = [("A"), ("B"), ("C"), ("D")]
    ax = df_bar_df.plot(subplots=True, kind="bar", stacked=True)
    subplot_data_df_list = _df_bar_xyheight_from_ax_helper(
        df_bar_data, ax, subplot_division
    )
    for i in range(len(subplot_data_df_list)):
        _df_bar_subplot_checker(
            df_bar_data, df_bar_df, subplot_data_df_list[i], subplot_division[i]
        )


def test_plot_bar_label_count_default():
    df = DataFrame(
        [(30, 10, 10, 10), (20, 20, 20, 20), (10, 30, 30, 10)], columns=list("ABCD")
    )
    df.plot(subplots=True, kind="bar", title=["A", "B", "C", "D"])


def test_plot_bar_label_count_expected_fail():
    df = DataFrame(
        [(30, 10, 10, 10), (20, 20, 20, 20), (10, 30, 30, 10)], columns=list("ABCD")
    )
    error_regex = re.escape(
        "The number of titles (4) must equal the number of subplots (3)."
    )
    with pytest.raises(ValueError, match=error_regex):
        df.plot(
            subplots=[("A", "B")],
            kind="bar",
            title=["A&B", "C", "D", "Extra Title"],
        )


def test_plot_bar_label_count_expected_success():
    df = DataFrame(
        [(30, 10, 10, 10), (20, 20, 20, 20), (10, 30, 30, 10)], columns=list("ABCD")
    )
    df.plot(subplots=[("A", "B", "D")], kind="bar", title=["A&B&D", "C"])
