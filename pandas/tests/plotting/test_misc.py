# coding: utf-8

""" Test cases for misc plot functions """

import numpy as np
from numpy import random
from numpy.random import randn
import pytest

from pandas.compat import lmap
import pandas.util._test_decorators as td

from pandas import DataFrame
from pandas.tests.plotting.common import TestPlotBase, _check_plot_works
import pandas.util.testing as tm

import pandas.plotting as plotting


@td.skip_if_mpl
def test_import_error_message():
    # GH-19810
    df = DataFrame({"A": [1, 2]})

    with pytest.raises(ImportError, match='matplotlib is required'):
        df.plot()


@td.skip_if_no_mpl
class TestSeriesPlots(TestPlotBase):

    def setup_method(self, method):
        TestPlotBase.setup_method(self, method)
        import matplotlib as mpl
        mpl.rcdefaults()

        self.ts = tm.makeTimeSeries()
        self.ts.name = 'ts'

    @pytest.mark.slow
    def test_autocorrelation_plot(self):
        from pandas.plotting import autocorrelation_plot
        _check_plot_works(autocorrelation_plot, series=self.ts)
        _check_plot_works(autocorrelation_plot, series=self.ts.values)

        ax = autocorrelation_plot(self.ts, label='Test')
        self._check_legend_labels(ax, labels=['Test'])

    @pytest.mark.slow
    def test_lag_plot(self):
        from pandas.plotting import lag_plot
        _check_plot_works(lag_plot, series=self.ts)
        _check_plot_works(lag_plot, series=self.ts, lag=5)

    @pytest.mark.slow
    def test_bootstrap_plot(self):
        from pandas.plotting import bootstrap_plot
        _check_plot_works(bootstrap_plot, series=self.ts, size=10)


@td.skip_if_no_mpl
class TestDataFramePlots(TestPlotBase):

    # This XPASSES when tested with mpl == 3.0.1
    @td.xfail_if_mpl_2_2
    @td.skip_if_no_scipy
    def test_scatter_matrix_axis(self):
        scatter_matrix = plotting.scatter_matrix

        with tm.RNGContext(42):
            df = DataFrame(randn(100, 3))

        # we are plotting multiples on a sub-plot
        with tm.assert_produces_warning(UserWarning):
            axes = _check_plot_works(scatter_matrix, filterwarnings='always',
                                     frame=df, range_padding=.1)
        axes0_labels = axes[0][0].yaxis.get_majorticklabels()

        # GH 5662
        expected = ['-2', '0', '2']
        self._check_text_labels(axes0_labels, expected)
        self._check_ticks_props(
            axes, xlabelsize=8, xrot=90, ylabelsize=8, yrot=0)

        df[0] = ((df[0] - 2) / 3)

        # we are plotting multiples on a sub-plot
        with tm.assert_produces_warning(UserWarning):
            axes = _check_plot_works(scatter_matrix, filterwarnings='always',
                                     frame=df, range_padding=.1)
        axes0_labels = axes[0][0].yaxis.get_majorticklabels()
        expected = ['-1.0', '-0.5', '0.0']
        self._check_text_labels(axes0_labels, expected)
        self._check_ticks_props(
            axes, xlabelsize=8, xrot=90, ylabelsize=8, yrot=0)

    @pytest.mark.slow
    def test_andrews_curves(self, iris):
        from pandas.plotting import andrews_curves
        from matplotlib import cm

        df = iris

        _check_plot_works(andrews_curves, frame=df, class_column='Name')

        rgba = ('#556270', '#4ECDC4', '#C7F464')
        ax = _check_plot_works(andrews_curves, frame=df,
                               class_column='Name', color=rgba)
        self._check_colors(
            ax.get_lines()[:10], linecolors=rgba, mapping=df['Name'][:10])

        cnames = ['dodgerblue', 'aquamarine', 'seagreen']
        ax = _check_plot_works(andrews_curves, frame=df,
                               class_column='Name', color=cnames)
        self._check_colors(
            ax.get_lines()[:10], linecolors=cnames, mapping=df['Name'][:10])

        ax = _check_plot_works(andrews_curves, frame=df,
                               class_column='Name', colormap=cm.jet)
        cmaps = lmap(cm.jet, np.linspace(0, 1, df['Name'].nunique()))
        self._check_colors(
            ax.get_lines()[:10], linecolors=cmaps, mapping=df['Name'][:10])

        length = 10
        df = DataFrame({"A": random.rand(length),
                        "B": random.rand(length),
                        "C": random.rand(length),
                        "Name": ["A"] * length})

        _check_plot_works(andrews_curves, frame=df, class_column='Name')

        rgba = ('#556270', '#4ECDC4', '#C7F464')
        ax = _check_plot_works(andrews_curves, frame=df,
                               class_column='Name', color=rgba)
        self._check_colors(
            ax.get_lines()[:10], linecolors=rgba, mapping=df['Name'][:10])

        cnames = ['dodgerblue', 'aquamarine', 'seagreen']
        ax = _check_plot_works(andrews_curves, frame=df,
                               class_column='Name', color=cnames)
        self._check_colors(
            ax.get_lines()[:10], linecolors=cnames, mapping=df['Name'][:10])

        ax = _check_plot_works(andrews_curves, frame=df,
                               class_column='Name', colormap=cm.jet)
        cmaps = lmap(cm.jet, np.linspace(0, 1, df['Name'].nunique()))
        self._check_colors(
            ax.get_lines()[:10], linecolors=cmaps, mapping=df['Name'][:10])

        colors = ['b', 'g', 'r']
        df = DataFrame({"A": [1, 2, 3],
                        "B": [1, 2, 3],
                        "C": [1, 2, 3],
                        "Name": colors})
        ax = andrews_curves(df, 'Name', color=colors)
        handles, labels = ax.get_legend_handles_labels()
        self._check_colors(handles, linecolors=colors)

        with tm.assert_produces_warning(FutureWarning):
            andrews_curves(data=df, class_column='Name')

    @pytest.mark.slow
    def test_parallel_coordinates(self, iris):
        from pandas.plotting import parallel_coordinates
        from matplotlib import cm

        df = iris

        ax = _check_plot_works(parallel_coordinates,
                               frame=df, class_column='Name')
        nlines = len(ax.get_lines())
        nxticks = len(ax.xaxis.get_ticklabels())

        rgba = ('#556270', '#4ECDC4', '#C7F464')
        ax = _check_plot_works(parallel_coordinates,
                               frame=df, class_column='Name', color=rgba)
        self._check_colors(
            ax.get_lines()[:10], linecolors=rgba, mapping=df['Name'][:10])

        cnames = ['dodgerblue', 'aquamarine', 'seagreen']
        ax = _check_plot_works(parallel_coordinates,
                               frame=df, class_column='Name', color=cnames)
        self._check_colors(
            ax.get_lines()[:10], linecolors=cnames, mapping=df['Name'][:10])

        ax = _check_plot_works(parallel_coordinates,
                               frame=df, class_column='Name', colormap=cm.jet)
        cmaps = lmap(cm.jet, np.linspace(0, 1, df['Name'].nunique()))
        self._check_colors(
            ax.get_lines()[:10], linecolors=cmaps, mapping=df['Name'][:10])

        ax = _check_plot_works(parallel_coordinates,
                               frame=df, class_column='Name', axvlines=False)
        assert len(ax.get_lines()) == (nlines - nxticks)

        colors = ['b', 'g', 'r']
        df = DataFrame({"A": [1, 2, 3],
                        "B": [1, 2, 3],
                        "C": [1, 2, 3],
                        "Name": colors})
        ax = parallel_coordinates(df, 'Name', color=colors)
        handles, labels = ax.get_legend_handles_labels()
        self._check_colors(handles, linecolors=colors)

        with tm.assert_produces_warning(FutureWarning):
            parallel_coordinates(data=df, class_column='Name')
        with tm.assert_produces_warning(FutureWarning):
            parallel_coordinates(df, 'Name', colors=colors)

    # not sure if this is indicative of a problem
    @pytest.mark.filterwarnings("ignore:Attempting to set:UserWarning")
    def test_parallel_coordinates_with_sorted_labels(self):
        """ For #15908 """
        from pandas.plotting import parallel_coordinates

        df = DataFrame({"feat": [i for i in range(30)],
                        "class": [2 for _ in range(10)] +
                                 [3 for _ in range(10)] +
                                 [1 for _ in range(10)]})
        ax = parallel_coordinates(df, 'class', sort_labels=True)
        polylines, labels = ax.get_legend_handles_labels()
        color_label_tuples = \
            zip([polyline.get_color() for polyline in polylines], labels)
        ordered_color_label_tuples = sorted(color_label_tuples,
                                            key=lambda x: x[1])
        prev_next_tupels = zip([i for i in ordered_color_label_tuples[0:-1]],
                               [i for i in ordered_color_label_tuples[1:]])
        for prev, nxt in prev_next_tupels:
            # labels and colors are ordered strictly increasing
            assert prev[1] < nxt[1] and prev[0] < nxt[0]

    @pytest.mark.slow
    def test_radviz(self, iris):
        from pandas.plotting import radviz
        from matplotlib import cm

        df = iris
        _check_plot_works(radviz, frame=df, class_column='Name')

        rgba = ('#556270', '#4ECDC4', '#C7F464')
        ax = _check_plot_works(
            radviz, frame=df, class_column='Name', color=rgba)
        # skip Circle drawn as ticks
        patches = [p for p in ax.patches[:20] if p.get_label() != '']
        self._check_colors(
            patches[:10], facecolors=rgba, mapping=df['Name'][:10])

        cnames = ['dodgerblue', 'aquamarine', 'seagreen']
        _check_plot_works(radviz, frame=df, class_column='Name', color=cnames)
        patches = [p for p in ax.patches[:20] if p.get_label() != '']
        self._check_colors(patches, facecolors=cnames, mapping=df['Name'][:10])

        _check_plot_works(radviz, frame=df,
                          class_column='Name', colormap=cm.jet)
        cmaps = lmap(cm.jet, np.linspace(0, 1, df['Name'].nunique()))
        patches = [p for p in ax.patches[:20] if p.get_label() != '']
        self._check_colors(patches, facecolors=cmaps, mapping=df['Name'][:10])

        colors = [[0., 0., 1., 1.],
                  [0., 0.5, 1., 1.],
                  [1., 0., 0., 1.]]
        df = DataFrame({"A": [1, 2, 3],
                        "B": [2, 1, 3],
                        "C": [3, 2, 1],
                        "Name": ['b', 'g', 'r']})
        ax = radviz(df, 'Name', color=colors)
        handles, labels = ax.get_legend_handles_labels()
        self._check_colors(handles, facecolors=colors)

    @pytest.mark.slow
    def test_subplot_titles(self, iris):
        df = iris.drop('Name', axis=1).head()
        # Use the column names as the subplot titles
        title = list(df.columns)

        # Case len(title) == len(df)
        plot = df.plot(subplots=True, title=title)
        assert [p.get_title() for p in plot] == title

        # Case len(title) > len(df)
        pytest.raises(ValueError, df.plot, subplots=True,
                      title=title + ["kittens > puppies"])

        # Case len(title) < len(df)
        pytest.raises(ValueError, df.plot, subplots=True, title=title[:2])

        # Case subplots=False and title is of type list
        pytest.raises(ValueError, df.plot, subplots=False, title=title)

        # Case df with 3 numeric columns but layout of (2,2)
        plot = df.drop('SepalWidth', axis=1).plot(subplots=True, layout=(2, 2),
                                                  title=title[:-1])
        title_list = [ax.get_title() for sublist in plot for ax in sublist]
        assert title_list == title[:3] + ['']

    def test_get_standard_colors_random_seed(self):
        # GH17525
        df = DataFrame(np.zeros((10, 10)))

        # Make sure that the random seed isn't reset by _get_standard_colors
        plotting.parallel_coordinates(df, 0)
        rand1 = random.random()
        plotting.parallel_coordinates(df, 0)
        rand2 = random.random()
        assert rand1 != rand2

        # Make sure it produces the same colors every time it's called
        from pandas.plotting._style import _get_standard_colors
        color1 = _get_standard_colors(1, color_type='random')
        color2 = _get_standard_colors(1, color_type='random')
        assert color1 == color2

    def test_get_standard_colors_default_num_colors(self):
        from pandas.plotting._style import _get_standard_colors

        # Make sure the default color_types returns the specified amount
        color1 = _get_standard_colors(1, color_type='default')
        color2 = _get_standard_colors(9, color_type='default')
        color3 = _get_standard_colors(20, color_type='default')
        assert len(color1) == 1
        assert len(color2) == 9
        assert len(color3) == 20

    def test_plot_single_color(self):
        # Example from #20585. All 3 bars should have the same color
        df = DataFrame({'account-start': ['2017-02-03', '2017-03-03',
                                          '2017-01-01'],
                        'client': ['Alice Anders', 'Bob Baker',
                                   'Charlie Chaplin'],
                        'balance': [-1432.32, 10.43, 30000.00],
                        'db-id': [1234, 2424, 251],
                        'proxy-id': [525, 1525, 2542],
                        'rank': [52, 525, 32],
                        })
        ax = df.client.value_counts().plot.bar()
        colors = lmap(lambda rect: rect.get_facecolor(),
                      ax.get_children()[0:3])
        assert all(color == colors[0] for color in colors)

    def test_get_standard_colors_no_appending(self):
        # GH20726

        # Make sure not to add more colors so that matplotlib can cycle
        # correctly.
        from matplotlib import cm
        color_before = cm.gnuplot(range(5))
        color_after = plotting._style._get_standard_colors(
            1, color=color_before)
        assert len(color_after) == len(color_before)

        df = DataFrame(np.random.randn(48, 4), columns=list("ABCD"))

        color_list = cm.gnuplot(np.linspace(0, 1, 16))
        p = df.A.plot.bar(figsize=(16, 7), color=color_list)
        assert (p.patches[1].get_facecolor()
                == p.patches[17].get_facecolor())
