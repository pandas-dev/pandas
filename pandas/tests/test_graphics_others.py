#!/usr/bin/env python
# coding: utf-8

import nose
import itertools
import os
import string
import warnings
from distutils.version import LooseVersion

from pandas import Series, DataFrame, MultiIndex
from pandas.compat import range, lmap, lzip
import pandas.util.testing as tm

import numpy as np
from numpy import random
from numpy.random import randn

from numpy.testing.decorators import slow
import pandas.tools.plotting as plotting

from pandas.tests.test_graphics import (TestPlotBase, _check_plot_works,
                                        curpath, _ok_for_gaussian_kde)


"""
These tests are for ``DataFrame.hist``, ``DataFrame.boxplot`` and
other miscellaneous plots.
`Dataframe.plot`` and ``Series.plot`` are tested in test_graphics.py
"""


def _skip_if_mpl_14_or_dev_boxplot():
    # GH 8382
    # Boxplot failures on 1.4 and 1.4.1
    # Don't need try / except since that's done at class level
    import matplotlib
    if str(matplotlib.__version__) >= LooseVersion('1.4'):
        raise nose.SkipTest("Matplotlib Regression in 1.4 and current dev.")


@tm.mplskip
class TestSeriesPlots(TestPlotBase):

    def setUp(self):
        TestPlotBase.setUp(self)
        import matplotlib as mpl
        mpl.rcdefaults()

        self.ts = tm.makeTimeSeries()
        self.ts.name = 'ts'

        self.series = tm.makeStringSeries()
        self.series.name = 'series'

        self.iseries = tm.makePeriodSeries()
        self.iseries.name = 'iseries'

    @slow
    def test_hist_legacy(self):
        _check_plot_works(self.ts.hist)
        _check_plot_works(self.ts.hist, grid=False)
        _check_plot_works(self.ts.hist, figsize=(8, 10))
        _check_plot_works(self.ts.hist, by=self.ts.index.month)
        _check_plot_works(self.ts.hist, by=self.ts.index.month, bins=5)

        fig, ax = self.plt.subplots(1, 1)
        _check_plot_works(self.ts.hist, ax=ax)
        _check_plot_works(self.ts.hist, ax=ax, figure=fig)
        _check_plot_works(self.ts.hist, figure=fig)
        tm.close()

        fig, (ax1, ax2) = self.plt.subplots(1, 2)
        _check_plot_works(self.ts.hist, figure=fig, ax=ax1)
        _check_plot_works(self.ts.hist, figure=fig, ax=ax2)

        with tm.assertRaises(ValueError):
            self.ts.hist(by=self.ts.index, figure=fig)

    @slow
    def test_hist_bins_legacy(self):
        df = DataFrame(np.random.randn(10, 2))
        ax = df.hist(bins=2)[0][0]
        self.assertEqual(len(ax.patches), 2)

    @slow
    def test_hist_layout(self):
        df = self.hist_df
        with tm.assertRaises(ValueError):
            df.height.hist(layout=(1, 1))

        with tm.assertRaises(ValueError):
            df.height.hist(layout=[1, 1])

    @slow
    def test_hist_layout_with_by(self):
        df = self.hist_df

        axes = _check_plot_works(df.height.hist, by=df.gender, layout=(2, 1))
        self._check_axes_shape(axes, axes_num=2, layout=(2, 1))

        axes = _check_plot_works(df.height.hist, by=df.gender, layout=(3, -1))
        self._check_axes_shape(axes, axes_num=2, layout=(3, 1))

        axes = _check_plot_works(df.height.hist, by=df.category, layout=(4, 1))
        self._check_axes_shape(axes, axes_num=4, layout=(4, 1))

        axes = _check_plot_works(
            df.height.hist, by=df.category, layout=(2, -1))
        self._check_axes_shape(axes, axes_num=4, layout=(2, 2))

        axes = _check_plot_works(
            df.height.hist, by=df.category, layout=(3, -1))
        self._check_axes_shape(axes, axes_num=4, layout=(3, 2))

        axes = _check_plot_works(
            df.height.hist, by=df.category, layout=(-1, 4))
        self._check_axes_shape(axes, axes_num=4, layout=(1, 4))

        axes = _check_plot_works(
            df.height.hist, by=df.classroom, layout=(2, 2))
        self._check_axes_shape(axes, axes_num=3, layout=(2, 2))

        axes = df.height.hist(by=df.category, layout=(4, 2), figsize=(12, 7))
        self._check_axes_shape(
            axes, axes_num=4, layout=(4, 2), figsize=(12, 7))

    @slow
    def test_hist_no_overlap(self):
        from matplotlib.pyplot import subplot, gcf
        x = Series(randn(2))
        y = Series(randn(2))
        subplot(121)
        x.hist()
        subplot(122)
        y.hist()
        fig = gcf()
        axes = fig.get_axes()
        self.assertEqual(len(axes), 2)

    @slow
    def test_hist_by_no_extra_plots(self):
        df = self.hist_df
        axes = df.height.hist(by=df.gender)  # noqa
        self.assertEqual(len(self.plt.get_fignums()), 1)

    @slow
    def test_plot_fails_when_ax_differs_from_figure(self):
        from pylab import figure
        fig1 = figure()
        fig2 = figure()
        ax1 = fig1.add_subplot(111)
        with tm.assertRaises(AssertionError):
            self.ts.hist(ax=ax1, figure=fig2)

    @slow
    def test_autocorrelation_plot(self):
        from pandas.tools.plotting import autocorrelation_plot
        _check_plot_works(autocorrelation_plot, series=self.ts)
        _check_plot_works(autocorrelation_plot, series=self.ts.values)

        ax = autocorrelation_plot(self.ts, label='Test')
        self._check_legend_labels(ax, labels=['Test'])

    @slow
    def test_lag_plot(self):
        from pandas.tools.plotting import lag_plot
        _check_plot_works(lag_plot, series=self.ts)
        _check_plot_works(lag_plot, series=self.ts, lag=5)

    @slow
    def test_bootstrap_plot(self):
        from pandas.tools.plotting import bootstrap_plot
        _check_plot_works(bootstrap_plot, series=self.ts, size=10)


@tm.mplskip
class TestDataFramePlots(TestPlotBase):

    def setUp(self):
        TestPlotBase.setUp(self)
        import matplotlib as mpl
        mpl.rcdefaults()

        self.tdf = tm.makeTimeDataFrame()
        self.hexbin_df = DataFrame({
            "A": np.random.uniform(size=20),
            "B": np.random.uniform(size=20),
            "C": np.arange(20) + np.random.uniform(size=20)})

        from pandas import read_csv
        path = os.path.join(curpath(), 'data', 'iris.csv')
        self.iris = read_csv(path)

    @slow
    def test_boxplot_legacy(self):
        df = DataFrame(randn(6, 4),
                       index=list(string.ascii_letters[:6]),
                       columns=['one', 'two', 'three', 'four'])
        df['indic'] = ['foo', 'bar'] * 3
        df['indic2'] = ['foo', 'bar', 'foo'] * 2

        _check_plot_works(df.boxplot, return_type='dict')
        _check_plot_works(df.boxplot, column=[
                          'one', 'two'], return_type='dict')
        _check_plot_works(df.boxplot, column=['one', 'two'], by='indic')
        _check_plot_works(df.boxplot, column='one', by=['indic', 'indic2'])
        _check_plot_works(df.boxplot, by='indic')
        _check_plot_works(df.boxplot, by=['indic', 'indic2'])
        _check_plot_works(plotting.boxplot, data=df['one'], return_type='dict')
        _check_plot_works(df.boxplot, notch=1, return_type='dict')
        _check_plot_works(df.boxplot, by='indic', notch=1)

        df = DataFrame(np.random.rand(10, 2), columns=['Col1', 'Col2'])
        df['X'] = Series(['A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B'])
        df['Y'] = Series(['A'] * 10)
        _check_plot_works(df.boxplot, by='X')

        # When ax is supplied and required number of axes is 1,
        # passed ax should be used:
        fig, ax = self.plt.subplots()
        axes = df.boxplot('Col1', by='X', ax=ax)
        self.assertIs(ax.get_axes(), axes)

        fig, ax = self.plt.subplots()
        axes = df.groupby('Y').boxplot(ax=ax, return_type='axes')
        self.assertIs(ax.get_axes(), axes['A'])

        # Multiple columns with an ax argument should use same figure
        fig, ax = self.plt.subplots()
        axes = df.boxplot(column=['Col1', 'Col2'],
                          by='X', ax=ax, return_type='axes')
        self.assertIs(axes['Col1'].get_figure(), fig)

        # When by is None, check that all relevant lines are present in the
        # dict
        fig, ax = self.plt.subplots()
        d = df.boxplot(ax=ax, return_type='dict')
        lines = list(itertools.chain.from_iterable(d.values()))
        self.assertEqual(len(ax.get_lines()), len(lines))

    @slow
    def test_boxplot_return_type_legacy(self):
        # API change in https://github.com/pydata/pandas/pull/7096
        import matplotlib as mpl  # noqa

        df = DataFrame(randn(6, 4),
                       index=list(string.ascii_letters[:6]),
                       columns=['one', 'two', 'three', 'four'])
        with tm.assertRaises(ValueError):
            df.boxplot(return_type='NOTATYPE')

        with tm.assert_produces_warning(FutureWarning):
            result = df.boxplot()
        # change to Axes in future
        self._check_box_return_type(result, 'dict')

        with tm.assert_produces_warning(False):
            result = df.boxplot(return_type='dict')
        self._check_box_return_type(result, 'dict')

        with tm.assert_produces_warning(False):
            result = df.boxplot(return_type='axes')
        self._check_box_return_type(result, 'axes')

        with tm.assert_produces_warning(False):
            result = df.boxplot(return_type='both')
        self._check_box_return_type(result, 'both')

    @slow
    def test_boxplot_axis_limits(self):

        def _check_ax_limits(col, ax):
            y_min, y_max = ax.get_ylim()
            self.assertTrue(y_min <= col.min())
            self.assertTrue(y_max >= col.max())

        df = self.hist_df.copy()
        df['age'] = np.random.randint(1, 20, df.shape[0])
        # One full row
        height_ax, weight_ax = df.boxplot(['height', 'weight'], by='category')
        _check_ax_limits(df['height'], height_ax)
        _check_ax_limits(df['weight'], weight_ax)
        self.assertEqual(weight_ax._sharey, height_ax)

        # Two rows, one partial
        p = df.boxplot(['height', 'weight', 'age'], by='category')
        height_ax, weight_ax, age_ax = p[0, 0], p[0, 1], p[1, 0]
        dummy_ax = p[1, 1]
        _check_ax_limits(df['height'], height_ax)
        _check_ax_limits(df['weight'], weight_ax)
        _check_ax_limits(df['age'], age_ax)
        self.assertEqual(weight_ax._sharey, height_ax)
        self.assertEqual(age_ax._sharey, height_ax)
        self.assertIsNone(dummy_ax._sharey)

    @slow
    def test_boxplot_empty_column(self):
        _skip_if_mpl_14_or_dev_boxplot()
        df = DataFrame(np.random.randn(20, 4))
        df.loc[:, 0] = np.nan
        _check_plot_works(df.boxplot, return_type='axes')

    @slow
    def test_hist_df_legacy(self):
        from matplotlib.patches import Rectangle
        _check_plot_works(self.hist_df.hist)

        # make sure layout is handled
        df = DataFrame(randn(100, 3))
        axes = _check_plot_works(df.hist, grid=False)
        self._check_axes_shape(axes, axes_num=3, layout=(2, 2))
        self.assertFalse(axes[1, 1].get_visible())

        df = DataFrame(randn(100, 1))
        _check_plot_works(df.hist)

        # make sure layout is handled
        df = DataFrame(randn(100, 6))
        axes = _check_plot_works(df.hist, layout=(4, 2))
        self._check_axes_shape(axes, axes_num=6, layout=(4, 2))

        # make sure sharex, sharey is handled
        _check_plot_works(df.hist, sharex=True, sharey=True)

        # handle figsize arg
        _check_plot_works(df.hist, figsize=(8, 10))

        # check bins argument
        _check_plot_works(df.hist, bins=5)

        # make sure xlabelsize and xrot are handled
        ser = df[0]
        xf, yf = 20, 18
        xrot, yrot = 30, 40
        axes = ser.hist(xlabelsize=xf, xrot=xrot, ylabelsize=yf, yrot=yrot)
        self._check_ticks_props(axes, xlabelsize=xf, xrot=xrot,
                                ylabelsize=yf, yrot=yrot)

        xf, yf = 20, 18
        xrot, yrot = 30, 40
        axes = df.hist(xlabelsize=xf, xrot=xrot, ylabelsize=yf, yrot=yrot)
        self._check_ticks_props(axes, xlabelsize=xf, xrot=xrot,
                                ylabelsize=yf, yrot=yrot)

        tm.close()
        # make sure kwargs to hist are handled
        ax = ser.hist(normed=True, cumulative=True, bins=4)
        # height of last bin (index 5) must be 1.0
        rects = [x for x in ax.get_children() if isinstance(x, Rectangle)]
        self.assertAlmostEqual(rects[-1].get_height(), 1.0)

        tm.close()
        ax = ser.hist(log=True)
        # scale of y must be 'log'
        self._check_ax_scales(ax, yaxis='log')

        tm.close()

        # propagate attr exception from matplotlib.Axes.hist
        with tm.assertRaises(AttributeError):
            ser.hist(foo='bar')

    @slow
    def test_hist_layout(self):
        df = DataFrame(randn(100, 3))

        layout_to_expected_size = (
            {'layout': None, 'expected_size': (2, 2)},  # default is 2x2
            {'layout': (2, 2), 'expected_size': (2, 2)},
            {'layout': (4, 1), 'expected_size': (4, 1)},
            {'layout': (1, 4), 'expected_size': (1, 4)},
            {'layout': (3, 3), 'expected_size': (3, 3)},
            {'layout': (-1, 4), 'expected_size': (1, 4)},
            {'layout': (4, -1), 'expected_size': (4, 1)},
            {'layout': (-1, 2), 'expected_size': (2, 2)},
            {'layout': (2, -1), 'expected_size': (2, 2)}
        )

        for layout_test in layout_to_expected_size:
            axes = df.hist(layout=layout_test['layout'])
            expected = layout_test['expected_size']
            self._check_axes_shape(axes, axes_num=3, layout=expected)

        # layout too small for all 4 plots
        with tm.assertRaises(ValueError):
            df.hist(layout=(1, 1))

        # invalid format for layout
        with tm.assertRaises(ValueError):
            df.hist(layout=(1,))
        with tm.assertRaises(ValueError):
            df.hist(layout=(-1, -1))

    @slow
    def test_scatter_plot_legacy(self):
        tm._skip_if_no_scipy()

        df = DataFrame(randn(100, 2))

        def scat(**kwds):
            return plotting.scatter_matrix(df, **kwds)

        _check_plot_works(scat)
        _check_plot_works(scat, marker='+')
        _check_plot_works(scat, vmin=0)
        if _ok_for_gaussian_kde('kde'):
            _check_plot_works(scat, diagonal='kde')
        if _ok_for_gaussian_kde('density'):
            _check_plot_works(scat, diagonal='density')
        _check_plot_works(scat, diagonal='hist')
        _check_plot_works(scat, range_padding=.1)

        def scat2(x, y, by=None, ax=None, figsize=None):
            return plotting.scatter_plot(df, x, y, by, ax, figsize=None)

        _check_plot_works(scat2, x=0, y=1)
        grouper = Series(np.repeat([1, 2, 3, 4, 5], 20), df.index)
        _check_plot_works(scat2, x=0, y=1, by=grouper)

    def test_scatter_matrix_axis(self):
        tm._skip_if_no_scipy()
        scatter_matrix = plotting.scatter_matrix

        with tm.RNGContext(42):
            df = DataFrame(randn(100, 3))

        # we are plotting multiples on a sub-plot
        with tm.assert_produces_warning(UserWarning):
            axes = _check_plot_works(scatter_matrix, filterwarnings='always',
                                     frame=df, range_padding=.1)
        axes0_labels = axes[0][0].yaxis.get_majorticklabels()

        # GH 5662
        expected = ['-2', '-1', '0', '1', '2']
        self._check_text_labels(axes0_labels, expected)
        self._check_ticks_props(
            axes, xlabelsize=8, xrot=90, ylabelsize=8, yrot=0)

        df[0] = ((df[0] - 2) / 3)

        # we are plotting multiples on a sub-plot
        with tm.assert_produces_warning(UserWarning):
            axes = _check_plot_works(scatter_matrix, filterwarnings='always',
                                     frame=df, range_padding=.1)
        axes0_labels = axes[0][0].yaxis.get_majorticklabels()
        expected = ['-1.2', '-1.0', '-0.8', '-0.6', '-0.4', '-0.2', '0.0']
        self._check_text_labels(axes0_labels, expected)
        self._check_ticks_props(
            axes, xlabelsize=8, xrot=90, ylabelsize=8, yrot=0)

    @slow
    def test_andrews_curves(self):
        from pandas.tools.plotting import andrews_curves
        from matplotlib import cm

        df = self.iris

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

    @slow
    def test_parallel_coordinates(self):
        from pandas.tools.plotting import parallel_coordinates
        from matplotlib import cm

        df = self.iris

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

    @slow
    def test_radviz(self):
        from pandas.tools.plotting import radviz
        from matplotlib import cm

        df = self.iris
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


@tm.mplskip
class TestDataFrameGroupByPlots(TestPlotBase):

    @slow
    def test_boxplot_legacy(self):
        grouped = self.hist_df.groupby(by='gender')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            axes = _check_plot_works(grouped.boxplot, return_type='axes')
        self._check_axes_shape(list(axes.values()), axes_num=2, layout=(1, 2))

        axes = _check_plot_works(grouped.boxplot, subplots=False,
                                 return_type='axes')
        self._check_axes_shape(axes, axes_num=1, layout=(1, 1))
        tuples = lzip(string.ascii_letters[:10], range(10))
        df = DataFrame(np.random.rand(10, 3),
                       index=MultiIndex.from_tuples(tuples))

        grouped = df.groupby(level=1)
        axes = _check_plot_works(grouped.boxplot, return_type='axes')
        self._check_axes_shape(list(axes.values()), axes_num=10, layout=(4, 3))

        axes = _check_plot_works(grouped.boxplot, subplots=False,
                                 return_type='axes')
        self._check_axes_shape(axes, axes_num=1, layout=(1, 1))

        grouped = df.unstack(level=1).groupby(level=0, axis=1)
        axes = _check_plot_works(grouped.boxplot, return_type='axes')
        self._check_axes_shape(list(axes.values()), axes_num=3, layout=(2, 2))

        axes = _check_plot_works(grouped.boxplot, subplots=False,
                                 return_type='axes')
        self._check_axes_shape(axes, axes_num=1, layout=(1, 1))

    @slow
    def test_grouped_plot_fignums(self):
        n = 10
        weight = Series(np.random.normal(166, 20, size=n))
        height = Series(np.random.normal(60, 10, size=n))
        with tm.RNGContext(42):
            gender = np.random.choice(['male', 'female'], size=n)
        df = DataFrame({'height': height, 'weight': weight, 'gender': gender})
        gb = df.groupby('gender')

        res = gb.plot()
        self.assertEqual(len(self.plt.get_fignums()), 2)
        self.assertEqual(len(res), 2)
        tm.close()

        res = gb.boxplot(return_type='axes')
        self.assertEqual(len(self.plt.get_fignums()), 1)
        self.assertEqual(len(res), 2)
        tm.close()

        # now works with GH 5610 as gender is excluded
        res = df.groupby('gender').hist()
        tm.close()

    @slow
    def test_grouped_hist_legacy(self):
        from matplotlib.patches import Rectangle

        df = DataFrame(randn(500, 2), columns=['A', 'B'])
        df['C'] = np.random.randint(0, 4, 500)
        df['D'] = ['X'] * 500

        axes = plotting.grouped_hist(df.A, by=df.C)
        self._check_axes_shape(axes, axes_num=4, layout=(2, 2))

        tm.close()
        axes = df.hist(by=df.C)
        self._check_axes_shape(axes, axes_num=4, layout=(2, 2))

        tm.close()
        # group by a key with single value
        axes = df.hist(by='D', rot=30)
        self._check_axes_shape(axes, axes_num=1, layout=(1, 1))
        self._check_ticks_props(axes, xrot=30)

        tm.close()
        # make sure kwargs to hist are handled
        xf, yf = 20, 18
        xrot, yrot = 30, 40
        axes = plotting.grouped_hist(df.A, by=df.C, normed=True,
                                     cumulative=True, bins=4,
                                     xlabelsize=xf, xrot=xrot,
                                     ylabelsize=yf, yrot=yrot)
        # height of last bin (index 5) must be 1.0
        for ax in axes.ravel():
            rects = [x for x in ax.get_children() if isinstance(x, Rectangle)]
            height = rects[-1].get_height()
            self.assertAlmostEqual(height, 1.0)
        self._check_ticks_props(axes, xlabelsize=xf, xrot=xrot,
                                ylabelsize=yf, yrot=yrot)

        tm.close()
        axes = plotting.grouped_hist(df.A, by=df.C, log=True)
        # scale of y must be 'log'
        self._check_ax_scales(axes, yaxis='log')

        tm.close()
        # propagate attr exception from matplotlib.Axes.hist
        with tm.assertRaises(AttributeError):
            plotting.grouped_hist(df.A, by=df.C, foo='bar')

        with tm.assert_produces_warning(FutureWarning):
            df.hist(by='C', figsize='default')

    @slow
    def test_grouped_hist_legacy2(self):
        n = 10
        weight = Series(np.random.normal(166, 20, size=n))
        height = Series(np.random.normal(60, 10, size=n))
        with tm.RNGContext(42):
            gender_int = np.random.choice([0, 1], size=n)
        df_int = DataFrame({'height': height, 'weight': weight,
                            'gender': gender_int})
        gb = df_int.groupby('gender')
        axes = gb.hist()
        self.assertEqual(len(axes), 2)
        self.assertEqual(len(self.plt.get_fignums()), 2)
        tm.close()

    @slow
    def test_grouped_box_return_type(self):
        df = self.hist_df

        # old style: return_type=None
        result = df.boxplot(by='gender')
        self.assertIsInstance(result, np.ndarray)
        self._check_box_return_type(
            result, None,
            expected_keys=['height', 'weight', 'category'])

        # now for groupby
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = df.groupby('gender').boxplot()
        self._check_box_return_type(
            result, 'dict', expected_keys=['Male', 'Female'])

        columns2 = 'X B C D A G Y N Q O'.split()
        df2 = DataFrame(random.randn(50, 10), columns=columns2)
        categories2 = 'A B C D E F G H I J'.split()
        df2['category'] = categories2 * 5

        for t in ['dict', 'axes', 'both']:
            returned = df.groupby('classroom').boxplot(return_type=t)
            self._check_box_return_type(
                returned, t, expected_keys=['A', 'B', 'C'])

            returned = df.boxplot(by='classroom', return_type=t)
            self._check_box_return_type(
                returned, t,
                expected_keys=['height', 'weight', 'category'])

            returned = df2.groupby('category').boxplot(return_type=t)
            self._check_box_return_type(returned, t, expected_keys=categories2)

            returned = df2.boxplot(by='category', return_type=t)
            self._check_box_return_type(returned, t, expected_keys=columns2)

    @slow
    def test_grouped_box_layout(self):
        df = self.hist_df

        self.assertRaises(ValueError, df.boxplot, column=['weight', 'height'],
                          by=df.gender, layout=(1, 1))
        self.assertRaises(ValueError, df.boxplot,
                          column=['height', 'weight', 'category'],
                          layout=(2, 1), return_type='dict')
        self.assertRaises(ValueError, df.boxplot, column=['weight', 'height'],
                          by=df.gender, layout=(-1, -1))

        box = _check_plot_works(df.groupby('gender').boxplot, column='height',
                                return_type='dict')
        self._check_axes_shape(self.plt.gcf().axes, axes_num=2, layout=(1, 2))

        box = _check_plot_works(df.groupby('category').boxplot,
                                column='height',
                                return_type='dict')
        self._check_axes_shape(self.plt.gcf().axes, axes_num=4, layout=(2, 2))

        # GH 6769
        box = _check_plot_works(df.groupby('classroom').boxplot,
                                column='height', return_type='dict')
        self._check_axes_shape(self.plt.gcf().axes, axes_num=3, layout=(2, 2))

        # GH 5897
        axes = df.boxplot(column=['height', 'weight', 'category'], by='gender',
                          return_type='axes')
        self._check_axes_shape(self.plt.gcf().axes, axes_num=3, layout=(2, 2))
        for ax in [axes['height']]:
            self._check_visible(ax.get_xticklabels(), visible=False)
            self._check_visible([ax.xaxis.get_label()], visible=False)
        for ax in [axes['weight'], axes['category']]:
            self._check_visible(ax.get_xticklabels())
            self._check_visible([ax.xaxis.get_label()])

        box = df.groupby('classroom').boxplot(
            column=['height', 'weight', 'category'], return_type='dict')
        self._check_axes_shape(self.plt.gcf().axes, axes_num=3, layout=(2, 2))

        box = _check_plot_works(df.groupby('category').boxplot,
                                column='height',
                                layout=(3, 2), return_type='dict')
        self._check_axes_shape(self.plt.gcf().axes, axes_num=4, layout=(3, 2))
        box = _check_plot_works(df.groupby('category').boxplot,
                                column='height',
                                layout=(3, -1), return_type='dict')
        self._check_axes_shape(self.plt.gcf().axes, axes_num=4, layout=(3, 2))

        box = df.boxplot(column=['height', 'weight', 'category'], by='gender',
                         layout=(4, 1))
        self._check_axes_shape(self.plt.gcf().axes, axes_num=3, layout=(4, 1))

        box = df.boxplot(column=['height', 'weight', 'category'], by='gender',
                         layout=(-1, 1))
        self._check_axes_shape(self.plt.gcf().axes, axes_num=3, layout=(3, 1))

        box = df.groupby('classroom').boxplot(
            column=['height', 'weight', 'category'], layout=(1, 4),
            return_type='dict')
        self._check_axes_shape(self.plt.gcf().axes, axes_num=3, layout=(1, 4))

        box = df.groupby('classroom').boxplot(  # noqa
            column=['height', 'weight', 'category'], layout=(1, -1),
            return_type='dict')
        self._check_axes_shape(self.plt.gcf().axes, axes_num=3, layout=(1, 3))

    @slow
    def test_grouped_box_multiple_axes(self):
        # GH 6970, GH 7069
        df = self.hist_df

        # check warning to ignore sharex / sharey
        # this check should be done in the first function which
        # passes multiple axes to plot, hist or boxplot
        # location should be changed if other test is added
        # which has earlier alphabetical order
        with tm.assert_produces_warning(UserWarning):
            fig, axes = self.plt.subplots(2, 2)
            df.groupby('category').boxplot(
                column='height', return_type='axes', ax=axes)
            self._check_axes_shape(self.plt.gcf().axes,
                                   axes_num=4, layout=(2, 2))

        fig, axes = self.plt.subplots(2, 3)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            returned = df.boxplot(column=['height', 'weight', 'category'],
                                  by='gender', return_type='axes', ax=axes[0])
        returned = np.array(list(returned.values()))
        self._check_axes_shape(returned, axes_num=3, layout=(1, 3))
        self.assert_numpy_array_equal(returned, axes[0])
        self.assertIs(returned[0].figure, fig)

        # draw on second row
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            returned = df.groupby('classroom').boxplot(
                column=['height', 'weight', 'category'],
                return_type='axes', ax=axes[1])
        returned = np.array(list(returned.values()))
        self._check_axes_shape(returned, axes_num=3, layout=(1, 3))
        self.assert_numpy_array_equal(returned, axes[1])
        self.assertIs(returned[0].figure, fig)

        with tm.assertRaises(ValueError):
            fig, axes = self.plt.subplots(2, 3)
            # pass different number of axes from required
            axes = df.groupby('classroom').boxplot(ax=axes)

    @slow
    def test_grouped_hist_layout(self):
        df = self.hist_df
        self.assertRaises(ValueError, df.hist, column='weight', by=df.gender,
                          layout=(1, 1))
        self.assertRaises(ValueError, df.hist, column='height', by=df.category,
                          layout=(1, 3))
        self.assertRaises(ValueError, df.hist, column='height', by=df.category,
                          layout=(-1, -1))

        axes = _check_plot_works(df.hist, column='height', by=df.gender,
                                 layout=(2, 1))
        self._check_axes_shape(axes, axes_num=2, layout=(2, 1))

        axes = _check_plot_works(df.hist, column='height', by=df.gender,
                                 layout=(2, -1))
        self._check_axes_shape(axes, axes_num=2, layout=(2, 1))

        axes = df.hist(column='height', by=df.category, layout=(4, 1))
        self._check_axes_shape(axes, axes_num=4, layout=(4, 1))

        axes = df.hist(column='height', by=df.category, layout=(-1, 1))
        self._check_axes_shape(axes, axes_num=4, layout=(4, 1))

        axes = df.hist(column='height', by=df.category,
                       layout=(4, 2), figsize=(12, 8))
        self._check_axes_shape(
            axes, axes_num=4, layout=(4, 2), figsize=(12, 8))
        tm.close()

        # GH 6769
        axes = _check_plot_works(
            df.hist, column='height', by='classroom', layout=(2, 2))
        self._check_axes_shape(axes, axes_num=3, layout=(2, 2))

        # without column
        axes = _check_plot_works(df.hist, by='classroom')
        self._check_axes_shape(axes, axes_num=3, layout=(2, 2))

        axes = df.hist(by='gender', layout=(3, 5))
        self._check_axes_shape(axes, axes_num=2, layout=(3, 5))

        axes = df.hist(column=['height', 'weight', 'category'])
        self._check_axes_shape(axes, axes_num=3, layout=(2, 2))

    @slow
    def test_grouped_hist_multiple_axes(self):
        # GH 6970, GH 7069
        df = self.hist_df

        fig, axes = self.plt.subplots(2, 3)
        returned = df.hist(column=['height', 'weight', 'category'], ax=axes[0])
        self._check_axes_shape(returned, axes_num=3, layout=(1, 3))
        self.assert_numpy_array_equal(returned, axes[0])
        self.assertIs(returned[0].figure, fig)
        returned = df.hist(by='classroom', ax=axes[1])
        self._check_axes_shape(returned, axes_num=3, layout=(1, 3))
        self.assert_numpy_array_equal(returned, axes[1])
        self.assertIs(returned[0].figure, fig)

        with tm.assertRaises(ValueError):
            fig, axes = self.plt.subplots(2, 3)
            # pass different number of axes from required
            axes = df.hist(column='height', ax=axes)

    @slow
    def test_axis_share_x(self):
        df = self.hist_df
        # GH4089
        ax1, ax2 = df.hist(column='height', by=df.gender, sharex=True)

        # share x
        self.assertTrue(ax1._shared_x_axes.joined(ax1, ax2))
        self.assertTrue(ax2._shared_x_axes.joined(ax1, ax2))

        # don't share y
        self.assertFalse(ax1._shared_y_axes.joined(ax1, ax2))
        self.assertFalse(ax2._shared_y_axes.joined(ax1, ax2))

    @slow
    def test_axis_share_y(self):
        df = self.hist_df
        ax1, ax2 = df.hist(column='height', by=df.gender, sharey=True)

        # share y
        self.assertTrue(ax1._shared_y_axes.joined(ax1, ax2))
        self.assertTrue(ax2._shared_y_axes.joined(ax1, ax2))

        # don't share x
        self.assertFalse(ax1._shared_x_axes.joined(ax1, ax2))
        self.assertFalse(ax2._shared_x_axes.joined(ax1, ax2))

    @slow
    def test_axis_share_xy(self):
        df = self.hist_df
        ax1, ax2 = df.hist(column='height', by=df.gender, sharex=True,
                           sharey=True)

        # share both x and y
        self.assertTrue(ax1._shared_x_axes.joined(ax1, ax2))
        self.assertTrue(ax2._shared_x_axes.joined(ax1, ax2))

        self.assertTrue(ax1._shared_y_axes.joined(ax1, ax2))
        self.assertTrue(ax2._shared_y_axes.joined(ax1, ax2))


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
