#!/usr/bin/env python
# coding: utf-8

import nose
import itertools
import os
import string
import warnings
from distutils.version import LooseVersion

from datetime import datetime, date

import pandas as pd
from pandas import (Series, DataFrame, MultiIndex, PeriodIndex, date_range,
                    bdate_range)
from pandas.compat import (range, lrange, StringIO, lmap, lzip, u, zip,
                           iteritems, OrderedDict, PY3)
from pandas.util.decorators import cache_readonly
import pandas.core.common as com
import pandas.util.testing as tm
from pandas.util.testing import ensure_clean
from pandas.core.config import set_option


import numpy as np
from numpy import random
from numpy.random import rand, randn

from numpy.testing import assert_allclose
from numpy.testing.decorators import slow
import pandas.tools.plotting as plotting


"""
These tests are for ``Dataframe.plot`` and ``Series.plot``.
Other plot methods such as ``.hist``, ``.boxplot`` and other miscellaneous
are tested in test_graphics_others.py
"""


def _skip_if_no_scipy_gaussian_kde():
    try:
        import scipy
        from scipy.stats import gaussian_kde
    except ImportError:
        raise nose.SkipTest("scipy version doesn't support gaussian_kde")


def _ok_for_gaussian_kde(kind):
    if kind in ['kde','density']:
        try:
            import scipy
            from scipy.stats import gaussian_kde
        except ImportError:
            return False
    return True


@tm.mplskip
class TestPlotBase(tm.TestCase):

    def setUp(self):

        import matplotlib as mpl
        mpl.rcdefaults()

        n = 100
        with tm.RNGContext(42):
            gender = tm.choice(['Male', 'Female'], size=n)
            classroom = tm.choice(['A', 'B', 'C'], size=n)

            self.hist_df = DataFrame({'gender': gender,
                                      'classroom': classroom,
                                      'height': random.normal(66, 4, size=n),
                                      'weight': random.normal(161, 32, size=n),
                                      'category': random.randint(4, size=n)})

        self.mpl_le_1_2_1 = plotting._mpl_le_1_2_1()
        self.mpl_ge_1_3_1 = plotting._mpl_ge_1_3_1()
        self.mpl_ge_1_4_0 = plotting._mpl_ge_1_4_0()
        self.mpl_ge_1_5_0 = plotting._mpl_ge_1_5_0()

        if self.mpl_ge_1_4_0:
            self.bp_n_objects = 7
        else:
            self.bp_n_objects = 8
        if self.mpl_ge_1_5_0:
            # 1.5 added PolyCollections to legend handler
            # so we have twice as many items.
            self.polycollection_factor = 2
        else:
            self.polycollection_factor = 1


    def tearDown(self):
        tm.close()

    @cache_readonly
    def plt(self):
        import matplotlib.pyplot as plt
        return plt

    @cache_readonly
    def colorconverter(self):
        import matplotlib.colors as colors
        return colors.colorConverter

    def _check_legend_labels(self, axes, labels=None, visible=True):
        """
        Check each axes has expected legend labels

        Parameters
        ----------
        axes : matplotlib Axes object, or its list-like
        labels : list-like
            expected legend labels
        visible : bool
            expected legend visibility. labels are checked only when visible is True
        """

        if visible and (labels is None):
            raise ValueError('labels must be specified when visible is True')
        axes = self._flatten_visible(axes)
        for ax in axes:
            if visible:
                self.assertTrue(ax.get_legend() is not None)
                self._check_text_labels(ax.get_legend().get_texts(), labels)
            else:
                self.assertTrue(ax.get_legend() is None)

    def _check_data(self, xp, rs):
        """
        Check each axes has identical lines

        Parameters
        ----------
        xp : matplotlib Axes object
        rs : matplotlib Axes object
        """
        xp_lines = xp.get_lines()
        rs_lines = rs.get_lines()

        def check_line(xpl, rsl):
            xpdata = xpl.get_xydata()
            rsdata = rsl.get_xydata()
            assert_allclose(xpdata, rsdata)

        self.assertEqual(len(xp_lines), len(rs_lines))
        [check_line(xpl, rsl) for xpl, rsl in zip(xp_lines, rs_lines)]
        tm.close()

    def _check_visible(self, collections, visible=True):
        """
        Check each artist is visible or not

        Parameters
        ----------
        collections : matplotlib Artist or its list-like
            target Artist or its list or collection
        visible : bool
            expected visibility
        """
        from matplotlib.collections import Collection
        if not isinstance(collections, Collection) and not com.is_list_like(collections):
            collections = [collections]

        for patch in collections:
            self.assertEqual(patch.get_visible(), visible)

    def _get_colors_mapped(self, series, colors):
        unique = series.unique()
        # unique and colors length can be differed
        # depending on slice value
        mapped = dict(zip(unique, colors))
        return [mapped[v] for v in series.values]

    def _check_colors(self, collections, linecolors=None, facecolors=None,
                      mapping=None):
        """
        Check each artist has expected line colors and face colors

        Parameters
        ----------
        collections : list-like
            list or collection of target artist
        linecolors : list-like which has the same length as collections
            list of expected line colors
        facecolors : list-like which has the same length as collections
            list of expected face colors
        mapping : Series
            Series used for color grouping key
            used for andrew_curves, parallel_coordinates, radviz test
        """

        from matplotlib.lines import Line2D
        from matplotlib.collections import Collection, PolyCollection
        conv = self.colorconverter
        if linecolors is not None:

            if mapping is not None:
                linecolors = self._get_colors_mapped(mapping, linecolors)
                linecolors = linecolors[:len(collections)]

            self.assertEqual(len(collections), len(linecolors))
            for patch, color in zip(collections, linecolors):
                if isinstance(patch, Line2D):
                    result = patch.get_color()
                    # Line2D may contains string color expression
                    result = conv.to_rgba(result)
                elif isinstance(patch, PolyCollection):
                    result = tuple(patch.get_edgecolor()[0])
                else:
                    result = patch.get_edgecolor()

                expected = conv.to_rgba(color)
                self.assertEqual(result, expected)

        if facecolors is not None:

            if mapping is not None:
                facecolors = self._get_colors_mapped(mapping, facecolors)
                facecolors = facecolors[:len(collections)]

            self.assertEqual(len(collections), len(facecolors))
            for patch, color in zip(collections, facecolors):
                if isinstance(patch, Collection):
                    # returned as list of np.array
                    result = patch.get_facecolor()[0]
                else:
                    result = patch.get_facecolor()

                if isinstance(result, np.ndarray):
                    result = tuple(result)

                expected = conv.to_rgba(color)
                self.assertEqual(result, expected)

    def _check_text_labels(self, texts, expected):
        """
        Check each text has expected labels

        Parameters
        ----------
        texts : matplotlib Text object, or its list-like
            target text, or its list
        expected : str or list-like which has the same length as texts
            expected text label, or its list
        """
        if not com.is_list_like(texts):
            self.assertEqual(texts.get_text(), expected)
        else:
            labels = [t.get_text() for t in texts]
            self.assertEqual(len(labels), len(expected))
            for l, e in zip(labels, expected):
                self.assertEqual(l, e)

    def _check_ticks_props(self, axes, xlabelsize=None, xrot=None,
                           ylabelsize=None, yrot=None):
        """
        Check each axes has expected tick properties

        Parameters
        ----------
        axes : matplotlib Axes object, or its list-like
        xlabelsize : number
            expected xticks font size
        xrot : number
            expected xticks rotation
        ylabelsize : number
            expected yticks font size
        yrot : number
            expected yticks rotation
        """
        from matplotlib.ticker import NullFormatter
        axes = self._flatten_visible(axes)
        for ax in axes:
            if xlabelsize or xrot:
                if isinstance(ax.xaxis.get_minor_formatter(), NullFormatter):
                    # If minor ticks has NullFormatter, rot / fontsize are not retained
                    labels = ax.get_xticklabels()
                else:
                    labels = ax.get_xticklabels() + ax.get_xticklabels(minor=True)

                for label in labels:
                    if xlabelsize is not None:
                        self.assertAlmostEqual(label.get_fontsize(), xlabelsize)
                    if xrot is not None:
                        self.assertAlmostEqual(label.get_rotation(), xrot)

            if ylabelsize or yrot:
                if isinstance(ax.yaxis.get_minor_formatter(), NullFormatter):
                    labels = ax.get_yticklabels()
                else:
                    labels = ax.get_yticklabels() + ax.get_yticklabels(minor=True)

                for label in labels:
                    if ylabelsize is not None:
                        self.assertAlmostEqual(label.get_fontsize(), ylabelsize)
                    if yrot is not None:
                        self.assertAlmostEqual(label.get_rotation(), yrot)

    def _check_ax_scales(self, axes, xaxis='linear', yaxis='linear'):
        """
        Check each axes has expected scales

        Parameters
        ----------
        axes : matplotlib Axes object, or its list-like
        xaxis : {'linear', 'log'}
            expected xaxis scale
        yaxis :  {'linear', 'log'}
            expected yaxis scale
        """
        axes = self._flatten_visible(axes)
        for ax in axes:
            self.assertEqual(ax.xaxis.get_scale(), xaxis)
            self.assertEqual(ax.yaxis.get_scale(), yaxis)

    def _check_axes_shape(self, axes, axes_num=None, layout=None, figsize=(8.0, 6.0)):
        """
        Check expected number of axes is drawn in expected layout

        Parameters
        ----------
        axes : matplotlib Axes object, or its list-like
        axes_num : number
            expected number of axes. Unnecessary axes should be set to invisible.
        layout :  tuple
            expected layout, (expected number of rows , columns)
        figsize : tuple
            expected figsize. default is matplotlib default
        """
        visible_axes = self._flatten_visible(axes)

        if axes_num is not None:
            self.assertEqual(len(visible_axes), axes_num)
            for ax in visible_axes:
                # check something drawn on visible axes
                self.assertTrue(len(ax.get_children()) > 0)

        if layout is not None:
            result = self._get_axes_layout(plotting._flatten(axes))
            self.assertEqual(result, layout)

        self.assert_numpy_array_equal(np.round(visible_axes[0].figure.get_size_inches()),
                                      np.array(figsize))

    def _get_axes_layout(self, axes):
        x_set = set()
        y_set = set()
        for ax in axes:
            # check axes coordinates to estimate layout
            points = ax.get_position().get_points()
            x_set.add(points[0][0])
            y_set.add(points[0][1])
        return (len(y_set), len(x_set))

    def _flatten_visible(self, axes):
        """
        Flatten axes, and filter only visible

        Parameters
        ----------
        axes : matplotlib Axes object, or its list-like

        """
        axes = plotting._flatten(axes)
        axes = [ax for ax in axes if ax.get_visible()]
        return axes

    def _check_has_errorbars(self, axes, xerr=0, yerr=0):
        """
        Check axes has expected number of errorbars

        Parameters
        ----------
        axes : matplotlib Axes object, or its list-like
        xerr : number
            expected number of x errorbar
        yerr : number
            expected number of y errorbar
        """
        axes = self._flatten_visible(axes)
        for ax in axes:
            containers = ax.containers
            xerr_count = 0
            yerr_count = 0
            for c in containers:
                has_xerr = getattr(c, 'has_xerr', False)
                has_yerr = getattr(c, 'has_yerr', False)
                if has_xerr:
                    xerr_count += 1
                if has_yerr:
                    yerr_count += 1
            self.assertEqual(xerr, xerr_count)
            self.assertEqual(yerr, yerr_count)

    def _check_box_return_type(self, returned, return_type, expected_keys=None,
                               check_ax_title=True):
        """
        Check box returned type is correct

        Parameters
        ----------
        returned : object to be tested, returned from boxplot
        return_type : str
            return_type passed to boxplot
        expected_keys : list-like, optional
            group labels in subplot case. If not passed,
            the function checks assuming boxplot uses single ax
        check_ax_title : bool
            Whether to check the ax.title is the same as expected_key
            Intended to be checked by calling from ``boxplot``.
            Normal ``plot`` doesn't attach ``ax.title``, it must be disabled.
        """
        from matplotlib.axes import Axes
        types = {'dict': dict, 'axes': Axes, 'both': tuple}
        if expected_keys is None:
            # should be fixed when the returning default is changed
            if return_type is None:
                return_type = 'dict'

            self.assertTrue(isinstance(returned, types[return_type]))
            if return_type == 'both':
                self.assertIsInstance(returned.ax, Axes)
                self.assertIsInstance(returned.lines, dict)
        else:
            # should be fixed when the returning default is changed
            if return_type is None:
                for r in self._flatten_visible(returned):
                    self.assertIsInstance(r, Axes)
                return

            self.assertTrue(isinstance(returned, OrderedDict))
            self.assertEqual(sorted(returned.keys()), sorted(expected_keys))
            for key, value in iteritems(returned):
                self.assertTrue(isinstance(value, types[return_type]))
                # check returned dict has correct mapping
                if return_type == 'axes':
                    if check_ax_title:
                        self.assertEqual(value.get_title(), key)
                elif return_type == 'both':
                    if check_ax_title:
                        self.assertEqual(value.ax.get_title(), key)
                    self.assertIsInstance(value.ax, Axes)
                    self.assertIsInstance(value.lines, dict)
                elif return_type == 'dict':
                    line = value['medians'][0]
                    if check_ax_title:
                        self.assertEqual(line.get_axes().get_title(), key)
                else:
                    raise AssertionError

    def _check_grid_settings(self, obj, kinds, kws={}):
        # Make sure plot defaults to rcParams['axes.grid'] setting, GH 9792

        import matplotlib as mpl

        def is_grid_on():
            xoff = all(not g.gridOn for g in self.plt.gca().xaxis.get_major_ticks())
            yoff = all(not g.gridOn for g in self.plt.gca().yaxis.get_major_ticks())
            return not(xoff and yoff)

        spndx=1
        for kind in kinds:
            if not _ok_for_gaussian_kde(kind):
                continue

            self.plt.subplot(1,4*len(kinds),spndx); spndx+=1
            mpl.rc('axes',grid=False)
            obj.plot(kind=kind, **kws)
            self.assertFalse(is_grid_on())

            self.plt.subplot(1,4*len(kinds),spndx); spndx+=1
            mpl.rc('axes',grid=True)
            obj.plot(kind=kind, grid=False, **kws)
            self.assertFalse(is_grid_on())

            if kind != 'pie':
                self.plt.subplot(1,4*len(kinds),spndx); spndx+=1
                mpl.rc('axes',grid=True)
                obj.plot(kind=kind, **kws)
                self.assertTrue(is_grid_on())

                self.plt.subplot(1,4*len(kinds),spndx); spndx+=1
                mpl.rc('axes',grid=False)
                obj.plot(kind=kind, grid=True, **kws)
                self.assertTrue(is_grid_on())

    def _maybe_unpack_cycler(self, rcParams, field='color'):
        """
        Compat layer for MPL 1.5 change to color cycle

        Before: plt.rcParams['axes.color_cycle'] -> ['b', 'g', 'r'...]
        After : plt.rcParams['axes.prop_cycle'] -> cycler(...)
        """
        if self.mpl_ge_1_5_0:
            cyl = rcParams['axes.prop_cycle']
            colors = [v[field] for v in cyl]
        else:
            colors = rcParams['axes.color_cycle']
        return colors


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
    def test_plot(self):
        _check_plot_works(self.ts.plot, label='foo')
        _check_plot_works(self.ts.plot, use_index=False)
        axes = _check_plot_works(self.ts.plot, rot=0)
        self._check_ticks_props(axes, xrot=0)

        ax = _check_plot_works(self.ts.plot, style='.', logy=True)
        self._check_ax_scales(ax, yaxis='log')

        ax = _check_plot_works(self.ts.plot, style='.', logx=True)
        self._check_ax_scales(ax, xaxis='log')

        ax = _check_plot_works(self.ts.plot, style='.', loglog=True)
        self._check_ax_scales(ax, xaxis='log', yaxis='log')

        _check_plot_works(self.ts[:10].plot.bar)
        _check_plot_works(self.ts.plot.area, stacked=False)
        _check_plot_works(self.iseries.plot)

        for kind in ['line', 'bar', 'barh', 'kde', 'hist', 'box']:
            if not _ok_for_gaussian_kde(kind):
                continue
            _check_plot_works(self.series[:5].plot, kind=kind)

        _check_plot_works(self.series[:10].plot.barh)
        ax = _check_plot_works(Series(randn(10)).plot.bar, color='black')
        self._check_colors([ax.patches[0]], facecolors=['black'])

        # GH 6951
        ax = _check_plot_works(self.ts.plot, subplots=True)
        self._check_axes_shape(ax, axes_num=1, layout=(1, 1))

        ax = _check_plot_works(self.ts.plot, subplots=True, layout=(-1, 1))
        self._check_axes_shape(ax, axes_num=1, layout=(1, 1))
        ax = _check_plot_works(self.ts.plot, subplots=True, layout=(1, -1))
        self._check_axes_shape(ax, axes_num=1, layout=(1, 1))

    @slow
    def test_plot_figsize_and_title(self):
        # figsize and title
        ax = self.series.plot(title='Test', figsize=(16, 8))
        self._check_text_labels(ax.title, 'Test')
        self._check_axes_shape(ax, axes_num=1, layout=(1, 1), figsize=(16, 8))

    def test_dont_modify_rcParams(self):
        # GH 8242
        if self.mpl_ge_1_5_0:
            key = 'axes.prop_cycle'
        else:
            key = 'axes.color_cycle'
        colors = self.plt.rcParams[key]
        Series([1, 2, 3]).plot()
        self.assertEqual(colors, self.plt.rcParams[key])

    def test_ts_line_lim(self):
        ax = self.ts.plot()
        xmin, xmax = ax.get_xlim()
        lines = ax.get_lines()
        self.assertEqual(xmin, lines[0].get_data(orig=False)[0][0])
        self.assertEqual(xmax, lines[0].get_data(orig=False)[0][-1])
        tm.close()

        ax = self.ts.plot(secondary_y=True)
        xmin, xmax = ax.get_xlim()
        lines = ax.get_lines()
        self.assertEqual(xmin, lines[0].get_data(orig=False)[0][0])
        self.assertEqual(xmax, lines[0].get_data(orig=False)[0][-1])

    def test_ts_area_lim(self):
        ax = self.ts.plot.area(stacked=False)
        xmin, xmax = ax.get_xlim()
        line = ax.get_lines()[0].get_data(orig=False)[0]
        self.assertEqual(xmin, line[0])
        self.assertEqual(xmax, line[-1])
        tm.close()

        # GH 7471
        ax = self.ts.plot.area(stacked=False, x_compat=True)
        xmin, xmax = ax.get_xlim()
        line = ax.get_lines()[0].get_data(orig=False)[0]
        self.assertEqual(xmin, line[0])
        self.assertEqual(xmax, line[-1])
        tm.close()

        tz_ts = self.ts.copy()
        tz_ts.index = tz_ts.tz_localize('GMT').tz_convert('CET')
        ax = tz_ts.plot.area(stacked=False, x_compat=True)
        xmin, xmax = ax.get_xlim()
        line = ax.get_lines()[0].get_data(orig=False)[0]
        self.assertEqual(xmin, line[0])
        self.assertEqual(xmax, line[-1])
        tm.close()

        ax = tz_ts.plot.area(stacked=False, secondary_y=True)
        xmin, xmax = ax.get_xlim()
        line = ax.get_lines()[0].get_data(orig=False)[0]
        self.assertEqual(xmin, line[0])
        self.assertEqual(xmax, line[-1])

    def test_label(self):
        s = Series([1, 2])
        ax = s.plot(label='LABEL', legend=True)
        self._check_legend_labels(ax, labels=['LABEL'])
        self.plt.close()
        ax = s.plot(legend=True)
        self._check_legend_labels(ax, labels=['None'])
        self.plt.close()
        # get name from index
        s.name = 'NAME'
        ax = s.plot(legend=True)
        self._check_legend_labels(ax, labels=['NAME'])
        self.plt.close()
        # override the default
        ax = s.plot(legend=True, label='LABEL')
        self._check_legend_labels(ax, labels=['LABEL'])
        self.plt.close()
        # Add lebel info, but don't draw
        ax = s.plot(legend=False, label='LABEL')
        self.assertEqual(ax.get_legend(), None)  # Hasn't been drawn
        ax.legend()  # draw it
        self._check_legend_labels(ax, labels=['LABEL'])

    def test_line_area_nan_series(self):
        values = [1, 2, np.nan, 3]
        s = Series(values)
        ts = Series(values, index=tm.makeDateIndex(k=4))

        for d in [s, ts]:
            ax = _check_plot_works(d.plot)
            masked = ax.lines[0].get_ydata()
            # remove nan for comparison purpose
            self.assert_numpy_array_equal(np.delete(masked.data, 2), np.array([1, 2, 3]))
            self.assert_numpy_array_equal(masked.mask, np.array([False, False, True, False]))

            expected = np.array([1, 2, 0, 3])
            ax = _check_plot_works(d.plot, stacked=True)
            self.assert_numpy_array_equal(ax.lines[0].get_ydata(), expected)
            ax = _check_plot_works(d.plot.area)
            self.assert_numpy_array_equal(ax.lines[0].get_ydata(), expected)
            ax = _check_plot_works(d.plot.area, stacked=False)
            self.assert_numpy_array_equal(ax.lines[0].get_ydata(), expected)

    def test_line_use_index_false(self):
        s = Series([1, 2, 3], index=['a', 'b', 'c'])
        s.index.name = 'The Index'
        ax = s.plot(use_index=False)
        label = ax.get_xlabel()
        self.assertEqual(label, '')
        ax2 = s.plot.bar(use_index=False)
        label2 = ax2.get_xlabel()
        self.assertEqual(label2, '')

    @slow
    def test_bar_log(self):
        expected = np.array([1., 10., 100., 1000.])

        if not self.mpl_le_1_2_1:
            expected = np.hstack((.1, expected, 1e4))

        ax = Series([200, 500]).plot.bar(log=True)
        tm.assert_numpy_array_equal(ax.yaxis.get_ticklocs(), expected)
        tm.close()

        ax = Series([200, 500]).plot.barh(log=True)
        tm.assert_numpy_array_equal(ax.xaxis.get_ticklocs(), expected)
        tm.close()

        # GH 9905
        expected = np.array([1.0e-03, 1.0e-02, 1.0e-01, 1.0e+00])

        if not self.mpl_le_1_2_1:
            expected = np.hstack((1.0e-04, expected, 1.0e+01))

        ax = Series([0.1, 0.01, 0.001]).plot(log=True, kind='bar')
        tm.assert_numpy_array_equal(ax.get_ylim(), (0.001, 0.10000000000000001))
        tm.assert_numpy_array_equal(ax.yaxis.get_ticklocs(), expected)
        tm.close()

        ax = Series([0.1, 0.01, 0.001]).plot(log=True, kind='barh')
        tm.assert_numpy_array_equal(ax.get_xlim(), (0.001, 0.10000000000000001))
        tm.assert_numpy_array_equal(ax.xaxis.get_ticklocs(), expected)

    @slow
    def test_bar_ignore_index(self):
        df = Series([1, 2, 3, 4], index=['a', 'b', 'c', 'd'])
        ax = df.plot.bar(use_index=False)
        self._check_text_labels(ax.get_xticklabels(), ['0', '1', '2', '3'])

    def test_rotation(self):
        df = DataFrame(randn(5, 5))
        # Default rot 0
        axes = df.plot()
        self._check_ticks_props(axes, xrot=0)

        axes = df.plot(rot=30)
        self._check_ticks_props(axes, xrot=30)

    def test_irregular_datetime(self):
        rng = date_range('1/1/2000', '3/1/2000')
        rng = rng[[0, 1, 2, 3, 5, 9, 10, 11, 12]]
        ser = Series(randn(len(rng)), rng)
        ax = ser.plot()
        xp = datetime(1999, 1, 1).toordinal()
        ax.set_xlim('1/1/1999', '1/1/2001')
        self.assertEqual(xp, ax.get_xlim()[0])

    @slow
    def test_pie_series(self):
        # if sum of values is less than 1.0, pie handle them as rate and draw semicircle.
        series = Series(np.random.randint(1, 5),
                        index=['a', 'b', 'c', 'd', 'e'], name='YLABEL')
        ax = _check_plot_works(series.plot.pie)
        self._check_text_labels(ax.texts, series.index)
        self.assertEqual(ax.get_ylabel(), 'YLABEL')

        # without wedge labels
        ax = _check_plot_works(series.plot.pie, labels=None)
        self._check_text_labels(ax.texts, [''] * 5)

        # with less colors than elements
        color_args = ['r', 'g', 'b']
        ax = _check_plot_works(series.plot.pie, colors=color_args)

        color_expected = ['r', 'g', 'b', 'r', 'g']
        self._check_colors(ax.patches, facecolors=color_expected)

        # with labels and colors
        labels = ['A', 'B', 'C', 'D', 'E']
        color_args = ['r', 'g', 'b', 'c', 'm']
        ax = _check_plot_works(series.plot.pie, labels=labels, colors=color_args)
        self._check_text_labels(ax.texts, labels)
        self._check_colors(ax.patches, facecolors=color_args)

        # with autopct and fontsize
        ax = _check_plot_works(series.plot.pie, colors=color_args,
                               autopct='%.2f', fontsize=7)
        pcts = ['{0:.2f}'.format(s * 100) for s in series.values / float(series.sum())]
        iters = [iter(series.index), iter(pcts)]
        expected_texts = list(next(it) for it in itertools.cycle(iters))
        self._check_text_labels(ax.texts, expected_texts)
        for t in ax.texts:
            self.assertEqual(t.get_fontsize(), 7)

        # includes negative value
        with tm.assertRaises(ValueError):
            series = Series([1, 2, 0, 4, -1], index=['a', 'b', 'c', 'd', 'e'])
            series.plot.pie()

        # includes nan
        series = Series([1, 2, np.nan, 4],
                        index=['a', 'b', 'c', 'd'], name='YLABEL')
        ax = _check_plot_works(series.plot.pie)
        self._check_text_labels(ax.texts, ['a', 'b', '', 'd'])

    def test_pie_nan(self):
        s = Series([1, np.nan, 1, 1])
        ax = s.plot.pie(legend=True)
        expected = ['0', '', '2', '3']
        result = [x.get_text() for x in ax.texts]
        self.assertEqual(result, expected)

    @slow
    def test_hist_df_kwargs(self):
        df = DataFrame(np.random.randn(10, 2))
        ax = df.plot.hist(bins=5)
        self.assertEqual(len(ax.patches), 10)

    @slow
    def test_hist_df_with_nonnumerics(self):
        # GH 9853
        with tm.RNGContext(1):
            df = DataFrame(np.random.randn(10, 4), columns=['A', 'B', 'C', 'D'])
        df['E'] = ['x', 'y'] * 5
        ax = df.plot.hist(bins=5)
        self.assertEqual(len(ax.patches), 20)

        ax = df.plot.hist() # bins=10
        self.assertEqual(len(ax.patches), 40)

    @slow
    def test_hist_legacy(self):
        _check_plot_works(self.ts.hist)
        _check_plot_works(self.ts.hist, grid=False)
        _check_plot_works(self.ts.hist, figsize=(8, 10))
        _check_plot_works(self.ts.hist, filterwarnings='ignore', by=self.ts.index.month)
        _check_plot_works(self.ts.hist, filterwarnings='ignore', by=self.ts.index.month, bins=5)

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

        axes = _check_plot_works(df.height.hist, filterwarnings='ignore',
                                 by=df.gender, layout=(2, 1))
        self._check_axes_shape(axes, axes_num=2, layout=(2, 1))

        axes = _check_plot_works(df.height.hist, filterwarnings='ignore',
                                 by=df.gender, layout=(3, -1))
        self._check_axes_shape(axes, axes_num=2, layout=(3, 1))

        axes = _check_plot_works(df.height.hist, filterwarnings='ignore',
                                 by=df.category, layout=(4, 1))
        self._check_axes_shape(axes, axes_num=4, layout=(4, 1))

        axes = _check_plot_works(df.height.hist, filterwarnings='ignore',
                                 by=df.category, layout=(2, -1))
        self._check_axes_shape(axes, axes_num=4, layout=(2, 2))

        axes = _check_plot_works(df.height.hist, filterwarnings='ignore',
                                 by=df.category, layout=(3, -1))
        self._check_axes_shape(axes, axes_num=4, layout=(3, 2))

        axes = _check_plot_works(df.height.hist, filterwarnings='ignore',
                                 by=df.category, layout=(-1, 4))
        self._check_axes_shape(axes, axes_num=4, layout=(1, 4))

        axes = _check_plot_works(df.height.hist, filterwarnings='ignore',
                                 by=df.classroom, layout=(2, 2))
        self._check_axes_shape(axes, axes_num=3, layout=(2, 2))

        axes = df.height.hist(by=df.category, layout=(4, 2), figsize=(12, 7))
        self._check_axes_shape(axes, axes_num=4, layout=(4, 2), figsize=(12, 7))

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
    def test_hist_secondary_legend(self):
        # GH 9610
        df = DataFrame(np.random.randn(30, 4), columns=list('abcd'))

        # primary -> secondary
        ax = df['a'].plot.hist(legend=True)
        df['b'].plot.hist(ax=ax, legend=True, secondary_y=True)
        # both legends are dran on left ax
        # left and right axis must be visible
        self._check_legend_labels(ax, labels=['a', 'b (right)'])
        self.assertTrue(ax.get_yaxis().get_visible())
        self.assertTrue(ax.right_ax.get_yaxis().get_visible())
        tm.close()

        # secondary -> secondary
        ax = df['a'].plot.hist(legend=True, secondary_y=True)
        df['b'].plot.hist(ax=ax, legend=True, secondary_y=True)
        # both legends are draw on left ax
        # left axis must be invisible, right axis must be visible
        self._check_legend_labels(ax.left_ax, labels=['a (right)', 'b (right)'])
        self.assertFalse(ax.left_ax.get_yaxis().get_visible())
        self.assertTrue(ax.get_yaxis().get_visible())
        tm.close()

        # secondary -> primary
        ax = df['a'].plot.hist(legend=True, secondary_y=True)
        # right axes is returned
        df['b'].plot.hist(ax=ax, legend=True)
        # both legends are draw on left ax
        # left and right axis must be visible
        self._check_legend_labels(ax.left_ax, labels=['a (right)', 'b'])
        self.assertTrue(ax.left_ax.get_yaxis().get_visible())
        self.assertTrue(ax.get_yaxis().get_visible())
        tm.close()

    @slow
    def test_df_series_secondary_legend(self):
        # GH 9779
        df = DataFrame(np.random.randn(30, 3), columns=list('abc'))
        s = Series(np.random.randn(30), name='x')

        # primary -> secondary (without passing ax)
        ax = df.plot()
        s.plot(legend=True, secondary_y=True)
        # both legends are dran on left ax
        # left and right axis must be visible
        self._check_legend_labels(ax, labels=['a', 'b', 'c', 'x (right)'])
        self.assertTrue(ax.get_yaxis().get_visible())
        self.assertTrue(ax.right_ax.get_yaxis().get_visible())
        tm.close()

        # primary -> secondary (with passing ax)
        ax = df.plot()
        s.plot(ax=ax, legend=True, secondary_y=True)
        # both legends are dran on left ax
        # left and right axis must be visible
        self._check_legend_labels(ax, labels=['a', 'b', 'c', 'x (right)'])
        self.assertTrue(ax.get_yaxis().get_visible())
        self.assertTrue(ax.right_ax.get_yaxis().get_visible())
        tm.close()

        # seconcary -> secondary (without passing ax)
        ax = df.plot(secondary_y=True)
        s.plot(legend=True, secondary_y=True)
        # both legends are dran on left ax
        # left axis must be invisible and right axis must be visible
        expected = ['a (right)', 'b (right)', 'c (right)', 'x (right)']
        self._check_legend_labels(ax.left_ax, labels=expected)
        self.assertFalse(ax.left_ax.get_yaxis().get_visible())
        self.assertTrue(ax.get_yaxis().get_visible())
        tm.close()

        # secondary -> secondary (with passing ax)
        ax = df.plot(secondary_y=True)
        s.plot(ax=ax, legend=True, secondary_y=True)
        # both legends are dran on left ax
        # left axis must be invisible and right axis must be visible
        expected = ['a (right)', 'b (right)', 'c (right)', 'x (right)']
        self._check_legend_labels(ax.left_ax, expected)
        self.assertFalse(ax.left_ax.get_yaxis().get_visible())
        self.assertTrue(ax.get_yaxis().get_visible())
        tm.close()

        # secondary -> secondary (with passing ax)
        ax = df.plot(secondary_y=True, mark_right=False)
        s.plot(ax=ax, legend=True, secondary_y=True)
        # both legends are dran on left ax
        # left axis must be invisible and right axis must be visible
        expected = ['a', 'b', 'c', 'x (right)']
        self._check_legend_labels(ax.left_ax, expected)
        self.assertFalse(ax.left_ax.get_yaxis().get_visible())
        self.assertTrue(ax.get_yaxis().get_visible())
        tm.close()

    @slow
    def test_plot_fails_with_dupe_color_and_style(self):
        x = Series(randn(2))
        with tm.assertRaises(ValueError):
            x.plot(style='k--', color='k')

    @slow
    def test_hist_kde(self):
        ax = self.ts.plot.hist(logy=True)
        self._check_ax_scales(ax, yaxis='log')
        xlabels = ax.get_xticklabels()
        # ticks are values, thus ticklabels are blank
        self._check_text_labels(xlabels, [''] * len(xlabels))
        ylabels = ax.get_yticklabels()
        self._check_text_labels(ylabels, [''] * len(ylabels))

        tm._skip_if_no_scipy()
        _skip_if_no_scipy_gaussian_kde()
        _check_plot_works(self.ts.plot.kde)
        _check_plot_works(self.ts.plot.density)
        ax = self.ts.plot.kde(logy=True)
        self._check_ax_scales(ax, yaxis='log')
        xlabels = ax.get_xticklabels()
        self._check_text_labels(xlabels, [''] * len(xlabels))
        ylabels = ax.get_yticklabels()
        self._check_text_labels(ylabels, [''] * len(ylabels))

    @slow
    def test_kde_kwargs(self):
        tm._skip_if_no_scipy()
        _skip_if_no_scipy_gaussian_kde()
        from numpy import linspace
        _check_plot_works(self.ts.plot.kde, bw_method=.5, ind=linspace(-100,100,20))
        _check_plot_works(self.ts.plot.density, bw_method=.5, ind=linspace(-100,100,20))
        ax = self.ts.plot.kde(logy=True, bw_method=.5, ind=linspace(-100,100,20))
        self._check_ax_scales(ax, yaxis='log')
        self._check_text_labels(ax.yaxis.get_label(), 'Density')

    @slow
    def test_kde_missing_vals(self):
        tm._skip_if_no_scipy()
        _skip_if_no_scipy_gaussian_kde()
        s = Series(np.random.uniform(size=50))
        s[0] = np.nan
        ax = _check_plot_works(s.plot.kde)

    @slow
    def test_hist_kwargs(self):
        ax = self.ts.plot.hist(bins=5)
        self.assertEqual(len(ax.patches), 5)
        self._check_text_labels(ax.yaxis.get_label(), 'Frequency')
        tm.close()

        if self.mpl_ge_1_3_1:
            ax = self.ts.plot.hist(orientation='horizontal')
            self._check_text_labels(ax.xaxis.get_label(), 'Frequency')
            tm.close()

            ax = self.ts.plot.hist(align='left', stacked=True)
            tm.close()

    @slow
    def test_hist_kde_color(self):
        ax = self.ts.plot.hist(logy=True, bins=10, color='b')
        self._check_ax_scales(ax, yaxis='log')
        self.assertEqual(len(ax.patches), 10)
        self._check_colors(ax.patches, facecolors=['b'] * 10)

        tm._skip_if_no_scipy()
        _skip_if_no_scipy_gaussian_kde()
        ax = self.ts.plot.kde(logy=True, color='r')
        self._check_ax_scales(ax, yaxis='log')
        lines = ax.get_lines()
        self.assertEqual(len(lines), 1)
        self._check_colors(lines, ['r'])

    @slow
    def test_boxplot_series(self):
        ax = self.ts.plot.box(logy=True)
        self._check_ax_scales(ax, yaxis='log')
        xlabels = ax.get_xticklabels()
        self._check_text_labels(xlabels, [self.ts.name])
        ylabels = ax.get_yticklabels()
        self._check_text_labels(ylabels, [''] * len(ylabels))

    @slow
    def test_kind_both_ways(self):
        s = Series(range(3))
        for kind in plotting._common_kinds + plotting._series_kinds:
            if not _ok_for_gaussian_kde(kind):
                continue
            s.plot(kind=kind)
            getattr(s.plot, kind)()

    @slow
    def test_invalid_plot_data(self):
        s = Series(list('abcd'))
        for kind in plotting._common_kinds:
            if not _ok_for_gaussian_kde(kind):
                continue
            with tm.assertRaises(TypeError):
                s.plot(kind=kind)

    @slow
    def test_valid_object_plot(self):
        s = Series(lrange(10), dtype=object)
        for kind in plotting._common_kinds:
            if not _ok_for_gaussian_kde(kind):
                continue
            _check_plot_works(s.plot, kind=kind)

    def test_partially_invalid_plot_data(self):
        s = Series(['a', 'b', 1.0, 2])
        for kind in plotting._common_kinds:
            if not _ok_for_gaussian_kde(kind):
                continue
            with tm.assertRaises(TypeError):
                s.plot(kind=kind)

    def test_invalid_kind(self):
        s = Series([1, 2])
        with tm.assertRaises(ValueError):
            s.plot(kind='aasdf')

    @slow
    def test_dup_datetime_index_plot(self):
        dr1 = date_range('1/1/2009', periods=4)
        dr2 = date_range('1/2/2009', periods=4)
        index = dr1.append(dr2)
        values = randn(index.size)
        s = Series(values, index=index)
        _check_plot_works(s.plot)

    @slow
    def test_errorbar_plot(self):

        s = Series(np.arange(10), name='x')
        s_err = np.random.randn(10)
        d_err =  DataFrame(randn(10, 2), index=s.index, columns=['x', 'y'])
        # test line and bar plots
        kinds = ['line', 'bar']
        for kind in kinds:
            ax = _check_plot_works(s.plot, yerr=Series(s_err), kind=kind)
            self._check_has_errorbars(ax, xerr=0, yerr=1)
            ax = _check_plot_works(s.plot, yerr=s_err, kind=kind)
            self._check_has_errorbars(ax, xerr=0, yerr=1)
            ax = _check_plot_works(s.plot, yerr=s_err.tolist(), kind=kind)
            self._check_has_errorbars(ax, xerr=0, yerr=1)
            ax = _check_plot_works(s.plot, yerr=d_err, kind=kind)
            self._check_has_errorbars(ax, xerr=0, yerr=1)
            ax = _check_plot_works(s.plot, xerr=0.2, yerr=0.2, kind=kind)
            self._check_has_errorbars(ax, xerr=1, yerr=1)

        ax = _check_plot_works(s.plot, xerr=s_err)
        self._check_has_errorbars(ax, xerr=1, yerr=0)

        # test time series plotting
        ix = date_range('1/1/2000', '1/1/2001', freq='M')
        ts = Series(np.arange(12), index=ix, name='x')
        ts_err = Series(np.random.randn(12), index=ix)
        td_err =  DataFrame(randn(12, 2), index=ix, columns=['x', 'y'])

        ax = _check_plot_works(ts.plot, yerr=ts_err)
        self._check_has_errorbars(ax, xerr=0, yerr=1)
        ax = _check_plot_works(ts.plot, yerr=td_err)
        self._check_has_errorbars(ax, xerr=0, yerr=1)

        # check incorrect lengths and types
        with tm.assertRaises(ValueError):
            s.plot(yerr=np.arange(11))

        s_err = ['zzz']*10
        # in mpl 1.5+ this is a TypeError
        with tm.assertRaises((ValueError, TypeError)):
            s.plot(yerr=s_err)

    def test_table(self):
        _check_plot_works(self.series.plot, table=True)
        _check_plot_works(self.series.plot, table=self.series)

    @slow
    def test_series_grid_settings(self):
        # Make sure plot defaults to rcParams['axes.grid'] setting, GH 9792
        self._check_grid_settings(Series([1,2,3]),
            plotting._series_kinds + plotting._common_kinds)

    @slow
    def test_standard_colors(self):
        for c in ['r', 'red', 'green', '#FF0000']:
            result = plotting._get_standard_colors(1, color=c)
            self.assertEqual(result, [c])

            result = plotting._get_standard_colors(1, color=[c])
            self.assertEqual(result, [c])

            result = plotting._get_standard_colors(3, color=c)
            self.assertEqual(result, [c] * 3)

            result = plotting._get_standard_colors(3, color=[c])
            self.assertEqual(result, [c] * 3)

    @slow
    def test_standard_colors_all(self):
        import matplotlib.colors as colors

        # multiple colors like mediumaquamarine
        for c in colors.cnames:
            result = plotting._get_standard_colors(num_colors=1, color=c)
            self.assertEqual(result, [c])

            result = plotting._get_standard_colors(num_colors=1, color=[c])
            self.assertEqual(result, [c])

            result = plotting._get_standard_colors(num_colors=3, color=c)
            self.assertEqual(result, [c] * 3)

            result = plotting._get_standard_colors(num_colors=3, color=[c])
            self.assertEqual(result, [c] * 3)

        # single letter colors like k
        for c in colors.ColorConverter.colors:
            result = plotting._get_standard_colors(num_colors=1, color=c)
            self.assertEqual(result, [c])

            result = plotting._get_standard_colors(num_colors=1, color=[c])
            self.assertEqual(result, [c])

            result = plotting._get_standard_colors(num_colors=3, color=c)
            self.assertEqual(result, [c] * 3)

            result = plotting._get_standard_colors(num_colors=3, color=[c])
            self.assertEqual(result, [c] * 3)

    def test_series_plot_color_kwargs(self):
        # GH1890
        ax = Series(np.arange(12) + 1).plot(color='green')
        self._check_colors(ax.get_lines(), linecolors=['green'])

    def test_time_series_plot_color_kwargs(self):
        # #1890
        ax = Series(np.arange(12) + 1, index=date_range(
            '1/1/2000', periods=12)).plot(color='green')
        self._check_colors(ax.get_lines(), linecolors=['green'])

    def test_time_series_plot_color_with_empty_kwargs(self):
        import matplotlib as mpl

        if self.mpl_ge_1_5_0:
            def_colors = self._maybe_unpack_cycler(mpl.rcParams)
        else:
            def_colors = mpl.rcParams['axes.color_cycle']
        index = date_range('1/1/2000', periods=12)
        s = Series(np.arange(1, 13), index=index)

        ncolors = 3

        for i in range(ncolors):
            ax = s.plot()
        self._check_colors(ax.get_lines(), linecolors=def_colors[:ncolors])

    def test_xticklabels(self):
        # GH11529
        s = Series(np.arange(10), index=['P%02d' % i for i in range(10)])
        ax = s.plot(xticks=[0,3,5,9])
        exp = ['P%02d' % i for i in [0,3,5,9]]
        self._check_text_labels(ax.get_xticklabels(), exp)


@tm.mplskip
class TestDataFramePlots(TestPlotBase):
    def setUp(self):
        TestPlotBase.setUp(self)
        import matplotlib as mpl
        mpl.rcdefaults()

        self.tdf = tm.makeTimeDataFrame()
        self.hexbin_df = DataFrame({"A": np.random.uniform(size=20),
                               "B": np.random.uniform(size=20),
                               "C": np.arange(20) + np.random.uniform(size=20)})

        from pandas import read_csv
        path = os.path.join(curpath(), 'data', 'iris.csv')
        self.iris = read_csv(path)

    @slow
    def test_plot(self):
        df = self.tdf
        _check_plot_works(df.plot, filterwarnings='ignore', grid=False)
        axes = _check_plot_works(df.plot, filterwarnings='ignore', subplots=True)
        self._check_axes_shape(axes, axes_num=4, layout=(4, 1))

        axes = _check_plot_works(df.plot, filterwarnings='ignore',
                                 subplots=True, layout=(-1, 2))
        self._check_axes_shape(axes, axes_num=4, layout=(2, 2))

        axes = _check_plot_works(df.plot, filterwarnings='ignore',
                                 subplots=True, use_index=False)
        self._check_axes_shape(axes, axes_num=4, layout=(4, 1))

        df = DataFrame({'x': [1, 2], 'y': [3, 4]})
        with tm.assertRaises(TypeError):
            df.plot.line(blarg=True)

        df = DataFrame(np.random.rand(10, 3),
                       index=list(string.ascii_letters[:10]))

        _check_plot_works(df.plot, use_index=True)
        _check_plot_works(df.plot, sort_columns=False)
        _check_plot_works(df.plot, yticks=[1, 5, 10])
        _check_plot_works(df.plot, xticks=[1, 5, 10])
        _check_plot_works(df.plot, ylim=(-100, 100), xlim=(-100, 100))

        _check_plot_works(df.plot, filterwarnings='ignore', subplots=True, title='blah')
        # We have to redo it here because _check_plot_works does two plots, once without an ax
        # kwarg and once with an ax kwarg and the new sharex behaviour does not remove the
        # visibility of the latter axis (as ax is present).
        # see: https://github.com/pydata/pandas/issues/9737
        axes = df.plot(subplots=True, title='blah')
        self._check_axes_shape(axes, axes_num=3, layout=(3, 1))
        #axes[0].figure.savefig("test.png")
        for ax in axes[:2]:
            self._check_visible(ax.xaxis)   # xaxis must be visible for grid
            self._check_visible(ax.get_xticklabels(), visible=False)
            self._check_visible(ax.get_xticklabels(minor=True), visible=False)
            self._check_visible([ax.xaxis.get_label()], visible=False)
        for ax in [axes[2]]:
            self._check_visible(ax.xaxis)
            self._check_visible(ax.get_xticklabels())
            self._check_visible([ax.xaxis.get_label()])
            self._check_ticks_props(ax, xrot=0)

        _check_plot_works(df.plot, title='blah')

        tuples = lzip(string.ascii_letters[:10], range(10))
        df = DataFrame(np.random.rand(10, 3),
                       index=MultiIndex.from_tuples(tuples))
        _check_plot_works(df.plot, use_index=True)

        # unicode
        index = MultiIndex.from_tuples([(u('\u03b1'), 0),
                                        (u('\u03b1'), 1),
                                        (u('\u03b2'), 2),
                                        (u('\u03b2'), 3),
                                        (u('\u03b3'), 4),
                                        (u('\u03b3'), 5),
                                        (u('\u03b4'), 6),
                                        (u('\u03b4'), 7)], names=['i0', 'i1'])
        columns = MultiIndex.from_tuples([('bar', u('\u0394')),
                                        ('bar', u('\u0395'))], names=['c0',
                                                                    'c1'])
        df = DataFrame(np.random.randint(0, 10, (8, 2)),
                       columns=columns,
                       index=index)
        _check_plot_works(df.plot, title=u('\u03A3'))

        # GH 6951
        # Test with single column
        df = DataFrame({'x': np.random.rand(10)})
        axes = _check_plot_works(df.plot.bar, subplots=True)
        self._check_axes_shape(axes, axes_num=1, layout=(1, 1))

        axes = _check_plot_works(df.plot.bar, subplots=True,
                                layout=(-1, 1))
        self._check_axes_shape(axes, axes_num=1, layout=(1, 1))
        # When ax is supplied and required number of axes is 1,
        # passed ax should be used:
        fig, ax = self.plt.subplots()
        axes = df.plot.bar(subplots=True, ax=ax)
        self.assertEqual(len(axes), 1)
        if self.mpl_ge_1_5_0:
            result = ax.axes
        else:
            result = ax.get_axes()  # deprecated
        self.assertIs(result, axes[0])

    def test_color_and_style_arguments(self):
        df = DataFrame({'x': [1, 2], 'y': [3, 4]})
        # passing both 'color' and 'style' arguments should be allowed
        # if there is no color symbol in the style strings:
        ax = df.plot(color = ['red', 'black'], style = ['-', '--'])
        # check that the linestyles are correctly set:
        linestyle = [line.get_linestyle() for line in ax.lines]
        self.assertEqual(linestyle, ['-', '--'])
        # check that the colors are correctly set:
        color = [line.get_color() for line in ax.lines]
        self.assertEqual(color, ['red', 'black'])
        # passing both 'color' and 'style' arguments should not be allowed
        # if there is a color symbol in the style strings:
        with tm.assertRaises(ValueError):
            df.plot(color = ['red', 'black'], style = ['k-', 'r--'])

    def test_nonnumeric_exclude(self):
        df = DataFrame({'A': ["x", "y", "z"], 'B': [1, 2, 3]})
        ax = df.plot()
        self.assertEqual(len(ax.get_lines()), 1)  # B was plotted

    @slow
    def test_implicit_label(self):
        df = DataFrame(randn(10, 3), columns=['a', 'b', 'c'])
        ax = df.plot(x='a', y='b')
        self._check_text_labels(ax.xaxis.get_label(), 'a')

    @slow
    def test_donot_overwrite_index_name(self):
        # GH 8494
        df = DataFrame(randn(2, 2), columns=['a', 'b'])
        df.index.name = 'NAME'
        df.plot(y='b', label='LABEL')
        self.assertEqual(df.index.name, 'NAME')

    @slow
    def test_plot_xy(self):
        # columns.inferred_type == 'string'
        df = self.tdf
        self._check_data(df.plot(x=0, y=1),
                         df.set_index('A')['B'].plot())
        self._check_data(df.plot(x=0), df.set_index('A').plot())
        self._check_data(df.plot(y=0), df.B.plot())
        self._check_data(df.plot(x='A', y='B'),
                         df.set_index('A').B.plot())
        self._check_data(df.plot(x='A'), df.set_index('A').plot())
        self._check_data(df.plot(y='B'), df.B.plot())

        # columns.inferred_type == 'integer'
        df.columns = lrange(1, len(df.columns) + 1)
        self._check_data(df.plot(x=1, y=2),
                         df.set_index(1)[2].plot())
        self._check_data(df.plot(x=1), df.set_index(1).plot())
        self._check_data(df.plot(y=1), df[1].plot())

        # figsize and title
        ax = df.plot(x=1, y=2, title='Test', figsize=(16, 8))
        self._check_text_labels(ax.title, 'Test')
        self._check_axes_shape(ax, axes_num=1, layout=(1, 1), figsize=(16., 8.))

        # columns.inferred_type == 'mixed'
        # TODO add MultiIndex test

    @slow
    def test_logscales(self):
        df = DataFrame({'a': np.arange(100)},
                       index=np.arange(100))
        ax = df.plot(logy=True)
        self._check_ax_scales(ax, yaxis='log')

        ax = df.plot(logx=True)
        self._check_ax_scales(ax, xaxis='log')

        ax = df.plot(loglog=True)
        self._check_ax_scales(ax, xaxis='log', yaxis='log')

    @slow
    def test_xcompat(self):
        import pandas as pd

        df = self.tdf
        ax = df.plot(x_compat=True)
        lines = ax.get_lines()
        self.assertNotIsInstance(lines[0].get_xdata(), PeriodIndex)

        tm.close()
        pd.plot_params['xaxis.compat'] = True
        ax = df.plot()
        lines = ax.get_lines()
        self.assertNotIsInstance(lines[0].get_xdata(), PeriodIndex)

        tm.close()
        pd.plot_params['x_compat'] = False
        ax = df.plot()
        lines = ax.get_lines()
        self.assertNotIsInstance(lines[0].get_xdata(), PeriodIndex)
        self.assertIsInstance(PeriodIndex(lines[0].get_xdata()), PeriodIndex)

        tm.close()
        # useful if you're plotting a bunch together
        with pd.plot_params.use('x_compat', True):
            ax = df.plot()
            lines = ax.get_lines()
            self.assertNotIsInstance(lines[0].get_xdata(), PeriodIndex)

        tm.close()
        ax = df.plot()
        lines = ax.get_lines()
        self.assertNotIsInstance(lines[0].get_xdata(), PeriodIndex)
        self.assertIsInstance(PeriodIndex(lines[0].get_xdata()), PeriodIndex)

    def test_period_compat(self):
        # GH 9012
        # period-array conversions
        df = DataFrame(
            np.random.rand(21, 2),
            index=bdate_range(datetime(2000, 1, 1), datetime(2000, 1, 31)),
            columns=['a', 'b'])

        df.plot()
        self.plt.axhline(y=0)
        tm.close()

    def test_unsorted_index(self):
        df = DataFrame({'y': np.arange(100)},
                       index=np.arange(99, -1, -1), dtype=np.int64)
        ax = df.plot()
        l = ax.get_lines()[0]
        rs = l.get_xydata()
        rs = Series(rs[:, 1], rs[:, 0], dtype=np.int64, name='y')
        tm.assert_series_equal(rs, df.y, check_index_type=False)
        tm.close()

        df.index = pd.Index(np.arange(99, -1, -1), dtype=np.float64)
        ax = df.plot()
        l = ax.get_lines()[0]
        rs = l.get_xydata()
        rs = Series(rs[:, 1], rs[:, 0], dtype=np.int64, name='y')
        tm.assert_series_equal(rs, df.y)

    @slow
    def test_subplots(self):
        df = DataFrame(np.random.rand(10, 3),
                       index=list(string.ascii_letters[:10]))

        for kind in ['bar', 'barh', 'line', 'area']:
            axes = df.plot(kind=kind, subplots=True, sharex=True, legend=True)
            self._check_axes_shape(axes, axes_num=3, layout=(3, 1))
            self.assertEqual(axes.shape, (3, ))

            for ax, column in zip(axes, df.columns):
                self._check_legend_labels(ax, labels=[com.pprint_thing(column)])

            for ax in axes[:-2]:
                self._check_visible(ax.xaxis)   # xaxis must be visible for grid
                self._check_visible(ax.get_xticklabels(), visible=False)
                self._check_visible(ax.get_xticklabels(minor=True), visible=False)
                self._check_visible(ax.xaxis.get_label(), visible=False)
                self._check_visible(ax.get_yticklabels())

            self._check_visible(axes[-1].xaxis)
            self._check_visible(axes[-1].get_xticklabels())
            self._check_visible(axes[-1].get_xticklabels(minor=True))
            self._check_visible(axes[-1].xaxis.get_label())
            self._check_visible(axes[-1].get_yticklabels())

            axes = df.plot(kind=kind, subplots=True, sharex=False)
            for ax in axes:
                self._check_visible(ax.xaxis)
                self._check_visible(ax.get_xticklabels())
                self._check_visible(ax.get_xticklabels(minor=True))
                self._check_visible(ax.xaxis.get_label())
                self._check_visible(ax.get_yticklabels())

            axes = df.plot(kind=kind, subplots=True, legend=False)
            for ax in axes:
                self.assertTrue(ax.get_legend() is None)

    @slow
    def test_subplots_timeseries(self):
        idx = date_range(start='2014-07-01', freq='M', periods=10)
        df = DataFrame(np.random.rand(10, 3), index=idx)

        for kind in ['line', 'area']:
            axes = df.plot(kind=kind, subplots=True, sharex=True)
            self._check_axes_shape(axes, axes_num=3, layout=(3, 1))

            for ax in axes[:-2]:
                # GH 7801
                self._check_visible(ax.xaxis)   # xaxis must be visible for grid
                self._check_visible(ax.get_xticklabels(), visible=False)
                self._check_visible(ax.get_xticklabels(minor=True), visible=False)
                self._check_visible(ax.xaxis.get_label(), visible=False)
                self._check_visible(ax.get_yticklabels())

            self._check_visible(axes[-1].xaxis)
            self._check_visible(axes[-1].get_xticklabels())
            self._check_visible(axes[-1].get_xticklabels(minor=True))
            self._check_visible(axes[-1].xaxis.get_label())
            self._check_visible(axes[-1].get_yticklabels())
            self._check_ticks_props(axes, xrot=0)

            axes = df.plot(kind=kind, subplots=True, sharex=False, rot=45, fontsize=7)
            for ax in axes:
                self._check_visible(ax.xaxis)
                self._check_visible(ax.get_xticklabels())
                self._check_visible(ax.get_xticklabels(minor=True))
                self._check_visible(ax.xaxis.get_label())
                self._check_visible(ax.get_yticklabels())
                self._check_ticks_props(ax, xlabelsize=7, xrot=45, ylabelsize=7)

    @slow
    def test_subplots_layout(self):
        # GH 6667
        df = DataFrame(np.random.rand(10, 3),
                       index=list(string.ascii_letters[:10]))

        axes = df.plot(subplots=True, layout=(2, 2))
        self._check_axes_shape(axes, axes_num=3, layout=(2, 2))
        self.assertEqual(axes.shape, (2, 2))

        axes = df.plot(subplots=True, layout=(-1, 2))
        self._check_axes_shape(axes, axes_num=3, layout=(2, 2))
        self.assertEqual(axes.shape, (2, 2))

        axes = df.plot(subplots=True, layout=(2, -1))
        self._check_axes_shape(axes, axes_num=3, layout=(2, 2))
        self.assertEqual(axes.shape, (2, 2))

        axes = df.plot(subplots=True, layout=(1, 4))
        self._check_axes_shape(axes, axes_num=3, layout=(1, 4))
        self.assertEqual(axes.shape, (1, 4))

        axes = df.plot(subplots=True, layout=(-1, 4))
        self._check_axes_shape(axes, axes_num=3, layout=(1, 4))
        self.assertEqual(axes.shape, (1, 4))

        axes = df.plot(subplots=True, layout=(4, -1))
        self._check_axes_shape(axes, axes_num=3, layout=(4, 1))
        self.assertEqual(axes.shape, (4, 1))

        with tm.assertRaises(ValueError):
            axes = df.plot(subplots=True, layout=(1, 1))
        with tm.assertRaises(ValueError):
            axes = df.plot(subplots=True, layout=(-1, -1))

        # single column
        df = DataFrame(np.random.rand(10, 1),
                       index=list(string.ascii_letters[:10]))
        axes = df.plot(subplots=True)
        self._check_axes_shape(axes, axes_num=1, layout=(1, 1))
        self.assertEqual(axes.shape, (1, ))

        axes = df.plot(subplots=True, layout=(3, 3))
        self._check_axes_shape(axes, axes_num=1, layout=(3, 3))
        self.assertEqual(axes.shape, (3, 3))

    @slow
    def test_subplots_warnings(self):
        # GH 9464
        warnings.simplefilter('error')
        try:
            df = DataFrame(np.random.randn(100, 4))
            df.plot(subplots=True, layout=(3, 2))

            df = DataFrame(np.random.randn(100, 4),
                           index=date_range('1/1/2000', periods=100))
            df.plot(subplots=True, layout=(3, 2))
        except Warning as w:
            self.fail(w)
        warnings.simplefilter('default')

    @slow
    def test_subplots_multiple_axes(self):
        # GH 5353, 6970, GH 7069
        fig, axes = self.plt.subplots(2, 3)
        df = DataFrame(np.random.rand(10, 3),
                       index=list(string.ascii_letters[:10]))

        returned = df.plot(subplots=True, ax=axes[0], sharex=False, sharey=False)
        self._check_axes_shape(returned, axes_num=3, layout=(1, 3))
        self.assertEqual(returned.shape, (3, ))
        self.assertIs(returned[0].figure, fig)
        # draw on second row
        returned = df.plot(subplots=True, ax=axes[1], sharex=False, sharey=False)
        self._check_axes_shape(returned, axes_num=3, layout=(1, 3))
        self.assertEqual(returned.shape, (3, ))
        self.assertIs(returned[0].figure, fig)
        self._check_axes_shape(axes, axes_num=6, layout=(2, 3))
        tm.close()

        with tm.assertRaises(ValueError):
            fig, axes = self.plt.subplots(2, 3)
            # pass different number of axes from required
            df.plot(subplots=True, ax=axes)

        # pass 2-dim axes and invalid layout
        # invalid lauout should not affect to input and return value
        # (show warning is tested in
        # TestDataFrameGroupByPlots.test_grouped_box_multiple_axes
        fig, axes = self.plt.subplots(2, 2)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            df = DataFrame(np.random.rand(10, 4),
                           index=list(string.ascii_letters[:10]))

            returned = df.plot(subplots=True, ax=axes, layout=(2, 1),
                               sharex=False, sharey=False)
            self._check_axes_shape(returned, axes_num=4, layout=(2, 2))
            self.assertEqual(returned.shape, (4, ))

            returned = df.plot(subplots=True, ax=axes, layout=(2, -1),
                               sharex=False, sharey=False)
            self._check_axes_shape(returned, axes_num=4, layout=(2, 2))
            self.assertEqual(returned.shape, (4, ))

            returned = df.plot(subplots=True, ax=axes, layout=(-1, 2),
                               sharex=False, sharey=False)
        self._check_axes_shape(returned, axes_num=4, layout=(2, 2))
        self.assertEqual(returned.shape, (4, ))

        # single column
        fig, axes = self.plt.subplots(1, 1)
        df = DataFrame(np.random.rand(10, 1),
                       index=list(string.ascii_letters[:10]))

        axes = df.plot(subplots=True, ax=[axes], sharex=False, sharey=False)
        self._check_axes_shape(axes, axes_num=1, layout=(1, 1))
        self.assertEqual(axes.shape, (1, ))

    def test_subplots_ts_share_axes(self):
        # GH 3964
        fig, axes = self.plt.subplots(3, 3, sharex=True, sharey=True)
        self.plt.subplots_adjust(left=0.05, right=0.95, hspace=0.3, wspace=0.3)
        df = DataFrame(np.random.randn(10, 9), index=date_range(start='2014-07-01', freq='M', periods=10))
        for i, ax in enumerate(axes.ravel()):
            df[i].plot(ax=ax, fontsize=5)

        #Rows other than bottom should not be visible
        for ax in axes[0:-1].ravel():
            self._check_visible(ax.get_xticklabels(), visible=False)

        #Bottom row should be visible
        for ax in axes[-1].ravel():
            self._check_visible(ax.get_xticklabels(), visible=True)

        #First column should be visible
        for ax in axes[[0, 1, 2], [0]].ravel():
            self._check_visible(ax.get_yticklabels(), visible=True)

        #Other columns should not be visible
        for ax in axes[[0, 1, 2], [1]].ravel():
            self._check_visible(ax.get_yticklabels(), visible=False)
        for ax in axes[[0, 1, 2], [2]].ravel():
            self._check_visible(ax.get_yticklabels(), visible=False)

    def test_subplots_sharex_axes_existing_axes(self):
        # GH 9158
        d = {'A': [1., 2., 3., 4.], 'B': [4., 3., 2., 1.], 'C': [5, 1, 3, 4]}
        df = DataFrame(d, index=date_range('2014 10 11', '2014 10 14'))

        axes = df[['A', 'B']].plot(subplots=True)
        df['C'].plot(ax=axes[0], secondary_y=True)

        self._check_visible(axes[0].get_xticklabels(), visible=False)
        self._check_visible(axes[1].get_xticklabels(), visible=True)
        for ax in axes.ravel():
            self._check_visible(ax.get_yticklabels(), visible=True)

    @slow
    def test_subplots_dup_columns(self):
        # GH 10962
        df = DataFrame(np.random.rand(5, 5), columns=list('aaaaa'))
        axes = df.plot(subplots=True)
        for ax in axes:
            self._check_legend_labels(ax, labels=['a'])
            self.assertEqual(len(ax.lines), 1)
        tm.close()

        axes = df.plot(subplots=True, secondary_y='a')
        for ax in axes:
            # (right) is only attached when subplots=False
            self._check_legend_labels(ax, labels=['a'])
            self.assertEqual(len(ax.lines), 1)
        tm.close()

        ax = df.plot(secondary_y='a')
        self._check_legend_labels(ax, labels=['a (right)'] * 5)
        self.assertEqual(len(ax.lines), 0)
        self.assertEqual(len(ax.right_ax.lines), 5)

    def test_negative_log(self):
        df = - DataFrame(rand(6, 4),
                       index=list(string.ascii_letters[:6]),
                       columns=['x', 'y', 'z', 'four'])

        with tm.assertRaises(ValueError):
            df.plot.area(logy=True)
        with tm.assertRaises(ValueError):
            df.plot.area(loglog=True)

    def _compare_stacked_y_cood(self, normal_lines, stacked_lines):
        base = np.zeros(len(normal_lines[0].get_data()[1]))
        for nl, sl in zip(normal_lines, stacked_lines):
            base += nl.get_data()[1] # get y coodinates
            sy = sl.get_data()[1]
            self.assert_numpy_array_equal(base, sy)

    def test_line_area_stacked(self):
        with tm.RNGContext(42):
            df = DataFrame(rand(6, 4),
                           columns=['w', 'x', 'y', 'z'])
            neg_df = - df
            # each column has either positive or negative value
            sep_df = DataFrame({'w': rand(6), 'x': rand(6),
                                'y': - rand(6), 'z': - rand(6)})
            # each column has positive-negative mixed value
            mixed_df = DataFrame(randn(6, 4), index=list(string.ascii_letters[:6]),
                                 columns=['w', 'x', 'y', 'z'])

            for kind in ['line', 'area']:
                ax1 = _check_plot_works(df.plot, kind=kind, stacked=False)
                ax2 = _check_plot_works(df.plot, kind=kind, stacked=True)
                self._compare_stacked_y_cood(ax1.lines, ax2.lines)

                ax1 = _check_plot_works(neg_df.plot, kind=kind, stacked=False)
                ax2 = _check_plot_works(neg_df.plot, kind=kind, stacked=True)
                self._compare_stacked_y_cood(ax1.lines, ax2.lines)

                ax1 = _check_plot_works(sep_df.plot, kind=kind, stacked=False)
                ax2 = _check_plot_works(sep_df.plot, kind=kind, stacked=True)
                self._compare_stacked_y_cood(ax1.lines[:2], ax2.lines[:2])
                self._compare_stacked_y_cood(ax1.lines[2:], ax2.lines[2:])

                _check_plot_works(mixed_df.plot, stacked=False)
                with tm.assertRaises(ValueError):
                    mixed_df.plot(stacked=True)

                _check_plot_works(df.plot, kind=kind, logx=True, stacked=True)

    def test_line_area_nan_df(self):
        values1 = [1, 2, np.nan, 3]
        values2 = [3, np.nan, 2, 1]
        df = DataFrame({'a': values1, 'b': values2})
        tdf = DataFrame({'a': values1, 'b': values2}, index=tm.makeDateIndex(k=4))

        for d in [df, tdf]:
            ax = _check_plot_works(d.plot)
            masked1 = ax.lines[0].get_ydata()
            masked2 = ax.lines[1].get_ydata()
            # remove nan for comparison purpose
            self.assert_numpy_array_equal(np.delete(masked1.data, 2), np.array([1, 2, 3]))
            self.assert_numpy_array_equal(np.delete(masked2.data, 1), np.array([3, 2, 1]))
            self.assert_numpy_array_equal(masked1.mask, np.array([False, False, True, False]))
            self.assert_numpy_array_equal(masked2.mask, np.array([False, True, False, False]))

            expected1 = np.array([1, 2, 0, 3])
            expected2 = np.array([3, 0, 2, 1])

            ax = _check_plot_works(d.plot, stacked=True)
            self.assert_numpy_array_equal(ax.lines[0].get_ydata(), expected1)
            self.assert_numpy_array_equal(ax.lines[1].get_ydata(), expected1 + expected2)

            ax = _check_plot_works(d.plot.area)
            self.assert_numpy_array_equal(ax.lines[0].get_ydata(), expected1)
            self.assert_numpy_array_equal(ax.lines[1].get_ydata(), expected1 + expected2)

            ax = _check_plot_works(d.plot.area, stacked=False)
            self.assert_numpy_array_equal(ax.lines[0].get_ydata(), expected1)
            self.assert_numpy_array_equal(ax.lines[1].get_ydata(), expected2)

    def test_line_lim(self):
        df = DataFrame(rand(6, 3), columns=['x', 'y', 'z'])
        ax = df.plot()
        xmin, xmax = ax.get_xlim()
        lines = ax.get_lines()
        self.assertEqual(xmin, lines[0].get_data()[0][0])
        self.assertEqual(xmax, lines[0].get_data()[0][-1])

        ax = df.plot(secondary_y=True)
        xmin, xmax = ax.get_xlim()
        lines = ax.get_lines()
        self.assertEqual(xmin, lines[0].get_data()[0][0])
        self.assertEqual(xmax, lines[0].get_data()[0][-1])

        axes = df.plot(secondary_y=True, subplots=True)
        self._check_axes_shape(axes, axes_num=3, layout=(3, 1))
        for ax in axes:
            self.assertTrue(hasattr(ax, 'left_ax'))
            self.assertFalse(hasattr(ax, 'right_ax'))
            xmin, xmax = ax.get_xlim()
            lines = ax.get_lines()
            self.assertEqual(xmin, lines[0].get_data()[0][0])
            self.assertEqual(xmax, lines[0].get_data()[0][-1])

    def test_area_lim(self):
        df = DataFrame(rand(6, 4),
                       columns=['x', 'y', 'z', 'four'])

        neg_df = - df
        for stacked in [True, False]:
            ax = _check_plot_works(df.plot.area, stacked=stacked)
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            lines = ax.get_lines()
            self.assertEqual(xmin, lines[0].get_data()[0][0])
            self.assertEqual(xmax, lines[0].get_data()[0][-1])
            self.assertEqual(ymin, 0)

            ax = _check_plot_works(neg_df.plot.area, stacked=stacked)
            ymin, ymax = ax.get_ylim()
            self.assertEqual(ymax, 0)

    @slow
    def test_bar_colors(self):
        import matplotlib.pyplot as plt
        default_colors = self._maybe_unpack_cycler(plt.rcParams)

        df = DataFrame(randn(5, 5))
        ax = df.plot.bar()
        self._check_colors(ax.patches[::5], facecolors=default_colors[:5])
        tm.close()

        custom_colors = 'rgcby'
        ax = df.plot.bar(color=custom_colors)
        self._check_colors(ax.patches[::5], facecolors=custom_colors)
        tm.close()

        from matplotlib import cm
        # Test str -> colormap functionality
        ax = df.plot.bar(colormap='jet')
        rgba_colors = lmap(cm.jet, np.linspace(0, 1, 5))
        self._check_colors(ax.patches[::5], facecolors=rgba_colors)
        tm.close()

        # Test colormap functionality
        ax = df.plot.bar(colormap=cm.jet)
        rgba_colors = lmap(cm.jet, np.linspace(0, 1, 5))
        self._check_colors(ax.patches[::5], facecolors=rgba_colors)
        tm.close()

        ax = df.ix[:, [0]].plot.bar(color='DodgerBlue')
        self._check_colors([ax.patches[0]], facecolors=['DodgerBlue'])
        tm.close()

        ax = df.plot(kind='bar', color='green')
        self._check_colors(ax.patches[::5], facecolors=['green'] * 5)
        tm.close()

    @slow
    def test_bar_linewidth(self):
        df = DataFrame(randn(5, 5))

        # regular
        ax = df.plot.bar(linewidth=2)
        for r in ax.patches:
            self.assertEqual(r.get_linewidth(), 2)

        # stacked
        ax = df.plot.bar(stacked=True, linewidth=2)
        for r in ax.patches:
            self.assertEqual(r.get_linewidth(), 2)

        # subplots
        axes = df.plot.bar(linewidth=2, subplots=True)
        self._check_axes_shape(axes, axes_num=5, layout=(5, 1))
        for ax in axes:
            for r in ax.patches:
                self.assertEqual(r.get_linewidth(), 2)

    @slow
    def test_bar_barwidth(self):
        df = DataFrame(randn(5, 5))

        width = 0.9

        # regular
        ax = df.plot.bar(width=width)
        for r in ax.patches:
            self.assertEqual(r.get_width(), width / len(df.columns))

        # stacked
        ax = df.plot.bar(stacked=True, width=width)
        for r in ax.patches:
            self.assertEqual(r.get_width(), width)

        # horizontal regular
        ax = df.plot.barh(width=width)
        for r in ax.patches:
            self.assertEqual(r.get_height(), width / len(df.columns))

        # horizontal stacked
        ax = df.plot.barh(stacked=True, width=width)
        for r in ax.patches:
            self.assertEqual(r.get_height(), width)

        # subplots
        axes = df.plot.bar(width=width, subplots=True)
        for ax in axes:
            for r in ax.patches:
                self.assertEqual(r.get_width(), width)

        # horizontal subplots
        axes = df.plot.barh(width=width, subplots=True)
        for ax in axes:
            for r in ax.patches:
                self.assertEqual(r.get_height(), width)

    @slow
    def test_bar_barwidth_position(self):
        df = DataFrame(randn(5, 5))
        self._check_bar_alignment(df, kind='bar', stacked=False, width=0.9, position=0.2)
        self._check_bar_alignment(df, kind='bar', stacked=True, width=0.9, position=0.2)
        self._check_bar_alignment(df, kind='barh', stacked=False, width=0.9, position=0.2)
        self._check_bar_alignment(df, kind='barh', stacked=True, width=0.9, position=0.2)
        self._check_bar_alignment(df, kind='bar', subplots=True, width=0.9, position=0.2)
        self._check_bar_alignment(df, kind='barh', subplots=True, width=0.9, position=0.2)

    @slow
    def test_bar_bottom_left(self):
        df = DataFrame(rand(5, 5))
        ax = df.plot.bar(stacked=False, bottom=1)
        result = [p.get_y() for p in ax.patches]
        self.assertEqual(result, [1] * 25)

        ax = df.plot.bar(stacked=True, bottom=[-1, -2, -3, -4, -5])
        result = [p.get_y() for p in ax.patches[:5]]
        self.assertEqual(result, [-1, -2, -3, -4, -5])

        ax = df.plot.barh(stacked=False, left=np.array([1, 1, 1, 1, 1]))
        result = [p.get_x() for p in ax.patches]
        self.assertEqual(result, [1] * 25)

        ax = df.plot.barh(stacked=True, left=[1, 2, 3, 4, 5])
        result = [p.get_x() for p in ax.patches[:5]]
        self.assertEqual(result, [1, 2, 3, 4, 5])

        axes = df.plot.bar(subplots=True, bottom=-1)
        for ax in axes:
            result = [p.get_y() for p in ax.patches]
            self.assertEqual(result, [-1] * 5)

        axes = df.plot.barh(subplots=True, left=np.array([1, 1, 1, 1, 1]))
        for ax in axes:
            result = [p.get_x() for p in ax.patches]
            self.assertEqual(result, [1] * 5)

    @slow
    def test_bar_nan(self):
        df = DataFrame({'A': [10, np.nan, 20], 'B': [5, 10, 20],
                        'C': [1, 2, 3]})
        ax = df.plot.bar()
        expected = [10, 0, 20, 5, 10, 20, 1, 2, 3]
        result = [p.get_height() for p in ax.patches]
        self.assertEqual(result, expected)

        ax = df.plot.bar(stacked=True)
        result = [p.get_height() for p in ax.patches]
        self.assertEqual(result, expected)

        result = [p.get_y() for p in ax.patches]
        expected = [0.0, 0.0, 0.0, 10.0, 0.0, 20.0, 15.0, 10.0, 40.0]
        self.assertEqual(result, expected)

    @slow
    def test_plot_scatter(self):
        df = DataFrame(randn(6, 4),
                       index=list(string.ascii_letters[:6]),
                       columns=['x', 'y', 'z', 'four'])

        _check_plot_works(df.plot.scatter, x='x', y='y')
        _check_plot_works(df.plot.scatter, x=1, y=2)

        with tm.assertRaises(TypeError):
            df.plot.scatter(x='x')
        with tm.assertRaises(TypeError):
            df.plot.scatter(y='y')

        # GH 6951
        axes = df.plot(x='x', y='y', kind='scatter', subplots=True)
        self._check_axes_shape(axes, axes_num=1, layout=(1, 1))

    @slow
    def test_plot_scatter_with_c(self):
        df = DataFrame(randn(6, 4),
                       index=list(string.ascii_letters[:6]),
                       columns=['x', 'y', 'z', 'four'])

        axes = [df.plot.scatter(x='x', y='y', c='z'),
                df.plot.scatter(x=0, y=1, c=2)]
        for ax in axes:
            # default to Greys
            self.assertEqual(ax.collections[0].cmap.name, 'Greys')

            if self.mpl_ge_1_3_1:

                # n.b. there appears to be no public method to get the colorbar
                # label
                self.assertEqual(ax.collections[0].colorbar._label, 'z')

        cm = 'cubehelix'
        ax = df.plot.scatter(x='x', y='y', c='z', colormap=cm)
        self.assertEqual(ax.collections[0].cmap.name, cm)

        # verify turning off colorbar works
        ax = df.plot.scatter(x='x', y='y', c='z', colorbar=False)
        self.assertIs(ax.collections[0].colorbar, None)

        # verify that we can still plot a solid color
        ax = df.plot.scatter(x=0, y=1, c='red')
        self.assertIs(ax.collections[0].colorbar, None)
        self._check_colors(ax.collections, facecolors=['r'])

        # Ensure that we can pass an np.array straight through to matplotlib,
        # this functionality was accidentally removed previously.
        # See https://github.com/pydata/pandas/issues/8852 for bug report
        #
        # Exercise colormap path and non-colormap path as they are independent
        #
        df = DataFrame({'A': [1, 2], 'B': [3, 4]})
        red_rgba = [1.0, 0.0, 0.0, 1.0]
        green_rgba = [0.0, 1.0, 0.0, 1.0]
        rgba_array = np.array([red_rgba, green_rgba])
        ax = df.plot.scatter(x='A', y='B', c=rgba_array)
        # expect the face colors of the points in the non-colormap path to be
        # identical to the values we supplied, normally we'd be on shaky ground
        # comparing floats for equality but here we expect them to be
        # identical.
        self.assertTrue(
            np.array_equal(
                ax.collections[0].get_facecolor(),
                rgba_array))
        # we don't test the colors of the faces in this next plot because they
        # are dependent on the spring colormap, which may change its colors
        # later.
        float_array = np.array([0.0, 1.0])
        df.plot.scatter(x='A', y='B', c=float_array, cmap='spring')

    def test_scatter_colors(self):
        df = DataFrame({'a': [1, 2, 3], 'b': [1, 2, 3], 'c': [1, 2, 3]})
        with tm.assertRaises(TypeError):
            df.plot.scatter(x='a', y='b', c='c', color='green')

        ax = df.plot.scatter(x='a', y='b', c='c')
        tm.assert_numpy_array_equal(ax.collections[0].get_facecolor()[0],
                                    (0, 0, 1, 1))

        ax = df.plot.scatter(x='a', y='b', color='white')
        tm.assert_numpy_array_equal(ax.collections[0].get_facecolor()[0],
                                    (1, 1, 1, 1))

    @slow
    def test_plot_bar(self):
        df = DataFrame(randn(6, 4),
                       index=list(string.ascii_letters[:6]),
                       columns=['one', 'two', 'three', 'four'])

        _check_plot_works(df.plot.bar)
        _check_plot_works(df.plot.bar, legend=False)
        _check_plot_works(df.plot.bar, filterwarnings='ignore', subplots=True)
        _check_plot_works(df.plot.bar, stacked=True)

        df = DataFrame(randn(10, 15),
                       index=list(string.ascii_letters[:10]),
                       columns=lrange(15))
        _check_plot_works(df.plot.bar)

        df = DataFrame({'a': [0, 1], 'b': [1, 0]})
        ax = _check_plot_works(df.plot.bar)
        self._check_ticks_props(ax, xrot=90)

        ax = df.plot.bar(rot=35, fontsize=10)
        self._check_ticks_props(ax, xrot=35, xlabelsize=10, ylabelsize=10)

        ax = _check_plot_works(df.plot.barh)
        self._check_ticks_props(ax, yrot=0)

        ax = df.plot.barh(rot=55, fontsize=11)
        self._check_ticks_props(ax, yrot=55, ylabelsize=11, xlabelsize=11)

    def _check_bar_alignment(self, df, kind='bar', stacked=False,
                             subplots=False, align='center',
                             width=0.5, position=0.5):

        axes = df.plot(kind=kind, stacked=stacked, subplots=subplots,
                       align=align, width=width, position=position,
                       grid=True)

        axes = self._flatten_visible(axes)

        for ax in axes:
            if kind == 'bar':
                axis = ax.xaxis
                ax_min, ax_max = ax.get_xlim()
                min_edge = min([p.get_x() for p in ax.patches])
                max_edge = max([p.get_x() + p.get_width() for p in ax.patches])
            elif kind == 'barh':
                axis = ax.yaxis
                ax_min, ax_max = ax.get_ylim()
                min_edge = min([p.get_y() for p in ax.patches])
                max_edge = max([p.get_y() + p.get_height() for p in ax.patches])
            else:
                raise ValueError

            # GH 7498
            # compare margins between lim and bar edges
            self.assertAlmostEqual(ax_min, min_edge - 0.25)
            self.assertAlmostEqual(ax_max, max_edge + 0.25)

            p = ax.patches[0]
            if kind == 'bar' and (stacked is True or subplots is True):
                edge = p.get_x()
                center = edge + p.get_width() * position
            elif kind == 'bar' and stacked is False:
                center = p.get_x() + p.get_width() * len(df.columns) * position
                edge = p.get_x()
            elif kind == 'barh' and (stacked is True or subplots is True):
                center = p.get_y() + p.get_height() * position
                edge = p.get_y()
            elif kind == 'barh' and stacked is False:
                center = p.get_y() + p.get_height() * len(df.columns) * position
                edge = p.get_y()
            else:
                raise ValueError

            # Check the ticks locates on integer
            self.assertTrue((axis.get_ticklocs() == np.arange(len(df))).all())

            if align == 'center':
                # Check whether the bar locates on center
                self.assertAlmostEqual(axis.get_ticklocs()[0], center)
            elif align == 'edge':
                # Check whether the bar's edge starts from the tick
                self.assertAlmostEqual(axis.get_ticklocs()[0], edge)
            else:
                raise ValueError

        return axes

    @slow
    def test_bar_stacked_center(self):
        # GH2157
        df = DataFrame({'A': [3] * 5, 'B': lrange(5)}, index=lrange(5))
        self._check_bar_alignment(df, kind='bar', stacked=True)
        self._check_bar_alignment(df, kind='bar', stacked=True, width=0.9)
        self._check_bar_alignment(df, kind='barh', stacked=True)
        self._check_bar_alignment(df, kind='barh', stacked=True, width=0.9)

    @slow
    def test_bar_center(self):
        df = DataFrame({'A': [3] * 5, 'B': lrange(5)}, index=lrange(5))
        self._check_bar_alignment(df, kind='bar', stacked=False)
        self._check_bar_alignment(df, kind='bar', stacked=False, width=0.9)
        self._check_bar_alignment(df, kind='barh', stacked=False)
        self._check_bar_alignment(df, kind='barh', stacked=False, width=0.9)

    @slow
    def test_bar_subplots_center(self):
        df = DataFrame({'A': [3] * 5, 'B': lrange(5)}, index=lrange(5))
        self._check_bar_alignment(df, kind='bar', subplots=True)
        self._check_bar_alignment(df, kind='bar', subplots=True, width=0.9)
        self._check_bar_alignment(df, kind='barh', subplots=True)
        self._check_bar_alignment(df, kind='barh', subplots=True, width=0.9)

    @slow
    def test_bar_align_single_column(self):
        df = DataFrame(randn(5))
        self._check_bar_alignment(df, kind='bar', stacked=False)
        self._check_bar_alignment(df, kind='bar', stacked=True)
        self._check_bar_alignment(df, kind='barh', stacked=False)
        self._check_bar_alignment(df, kind='barh', stacked=True)
        self._check_bar_alignment(df, kind='bar', subplots=True)
        self._check_bar_alignment(df, kind='barh', subplots=True)

    @slow
    def test_bar_edge(self):
        df = DataFrame({'A': [3] * 5, 'B': lrange(5)}, index=lrange(5))

        self._check_bar_alignment(df, kind='bar', stacked=True, align='edge')
        self._check_bar_alignment(df, kind='bar', stacked=True,
                                  width=0.9, align='edge')
        self._check_bar_alignment(df, kind='barh', stacked=True, align='edge')
        self._check_bar_alignment(df, kind='barh', stacked=True,
                                  width=0.9, align='edge')

        self._check_bar_alignment(df, kind='bar', stacked=False, align='edge')
        self._check_bar_alignment(df, kind='bar', stacked=False,
                                  width=0.9, align='edge')
        self._check_bar_alignment(df, kind='barh', stacked=False, align='edge')
        self._check_bar_alignment(df, kind='barh', stacked=False,
                                  width=0.9, align='edge')

        self._check_bar_alignment(df, kind='bar', subplots=True, align='edge')
        self._check_bar_alignment(df, kind='bar', subplots=True,
                                  width=0.9, align='edge')
        self._check_bar_alignment(df, kind='barh', subplots=True, align='edge')
        self._check_bar_alignment(df, kind='barh', subplots=True,
                                  width=0.9, align='edge')

    @slow
    def test_bar_log_no_subplots(self):
        # GH3254, GH3298 matplotlib/matplotlib#1882, #1892
        # regressions in 1.2.1
        expected = np.array([1., 10.])

        if not self.mpl_le_1_2_1:
            expected = np.hstack((.1, expected, 100))

        # no subplots
        df = DataFrame({'A': [3] * 5, 'B': lrange(1, 6)}, index=lrange(5))
        ax = df.plot.bar(grid=True, log=True)
        tm.assert_numpy_array_equal(ax.yaxis.get_ticklocs(), expected)

    @slow
    def test_bar_log_subplots(self):
        expected = np.array([1., 10., 100., 1000.])
        if not self.mpl_le_1_2_1:
            expected = np.hstack((.1, expected, 1e4))

        ax = DataFrame([Series([200, 300]),
                        Series([300, 500])]).plot.bar(log=True, subplots=True)

        tm.assert_numpy_array_equal(ax[0].yaxis.get_ticklocs(), expected)
        tm.assert_numpy_array_equal(ax[1].yaxis.get_ticklocs(), expected)

    @slow
    def test_boxplot(self):
        df = self.hist_df
        series = df['height']
        numeric_cols = df._get_numeric_data().columns
        labels = [com.pprint_thing(c) for c in numeric_cols]

        ax = _check_plot_works(df.plot.box)
        self._check_text_labels(ax.get_xticklabels(), labels)
        tm.assert_numpy_array_equal(ax.xaxis.get_ticklocs(),
                           np.arange(1, len(numeric_cols) + 1))
        self.assertEqual(len(ax.lines),
                         self.bp_n_objects * len(numeric_cols))

        # different warning on py3
        if not PY3:
            axes = _check_plot_works(df.plot.box,
                                     subplots=True, logy=True)

            self._check_axes_shape(axes, axes_num=3, layout=(1, 3))
            self._check_ax_scales(axes, yaxis='log')
            for ax, label in zip(axes, labels):
                self._check_text_labels(ax.get_xticklabels(), [label])
                self.assertEqual(len(ax.lines), self.bp_n_objects)

        axes = series.plot.box(rot=40)
        self._check_ticks_props(axes, xrot=40, yrot=0)
        tm.close()

        ax = _check_plot_works(series.plot.box)

        positions = np.array([1, 6, 7])
        ax = df.plot.box(positions=positions)
        numeric_cols = df._get_numeric_data().columns
        labels = [com.pprint_thing(c) for c in numeric_cols]
        self._check_text_labels(ax.get_xticklabels(), labels)
        tm.assert_numpy_array_equal(ax.xaxis.get_ticklocs(), positions)
        self.assertEqual(len(ax.lines), self.bp_n_objects * len(numeric_cols))

    @slow
    def test_boxplot_vertical(self):
        df = self.hist_df
        numeric_cols = df._get_numeric_data().columns
        labels = [com.pprint_thing(c) for c in numeric_cols]

        # if horizontal, yticklabels are rotated
        ax = df.plot.box(rot=50, fontsize=8, vert=False)
        self._check_ticks_props(ax, xrot=0, yrot=50, ylabelsize=8)
        self._check_text_labels(ax.get_yticklabels(), labels)
        self.assertEqual(len(ax.lines), self.bp_n_objects * len(numeric_cols))

        axes = _check_plot_works(df.plot.box, filterwarnings='ignore', subplots=True,
                                 vert=False, logx=True)
        self._check_axes_shape(axes, axes_num=3, layout=(1, 3))
        self._check_ax_scales(axes, xaxis='log')
        for ax, label in zip(axes, labels):
            self._check_text_labels(ax.get_yticklabels(), [label])
            self.assertEqual(len(ax.lines), self.bp_n_objects)

        positions = np.array([3, 2, 8])
        ax = df.plot.box(positions=positions, vert=False)
        self._check_text_labels(ax.get_yticklabels(), labels)
        tm.assert_numpy_array_equal(ax.yaxis.get_ticklocs(), positions)
        self.assertEqual(len(ax.lines), self.bp_n_objects * len(numeric_cols))

    @slow
    def test_boxplot_return_type(self):
        df = DataFrame(randn(6, 4),
                       index=list(string.ascii_letters[:6]),
                       columns=['one', 'two', 'three', 'four'])
        with tm.assertRaises(ValueError):
            df.plot.box(return_type='NOTATYPE')

        result = df.plot.box(return_type='dict')
        self._check_box_return_type(result, 'dict')

        result = df.plot.box(return_type='axes')
        self._check_box_return_type(result, 'axes')

        result = df.plot.box(return_type='both')
        self._check_box_return_type(result, 'both')

    @slow
    def test_boxplot_subplots_return_type(self):
        df = self.hist_df

        # normal style: return_type=None
        result = df.plot.box(subplots=True)
        self.assertIsInstance(result, np.ndarray)
        self._check_box_return_type(result, None,
                                    expected_keys=['height', 'weight', 'category'])

        for t in ['dict', 'axes', 'both']:
            returned = df.plot.box(return_type=t, subplots=True)
            self._check_box_return_type(returned, t,
                                        expected_keys=['height', 'weight', 'category'],
                                        check_ax_title=False)

    @slow
    def test_kde_df(self):
        tm._skip_if_no_scipy()
        _skip_if_no_scipy_gaussian_kde()
        df = DataFrame(randn(100, 4))
        ax = _check_plot_works(df.plot, kind='kde')
        expected = [com.pprint_thing(c) for c in df.columns]
        self._check_legend_labels(ax, labels=expected)
        self._check_ticks_props(ax, xrot=0)

        ax = df.plot(kind='kde', rot=20, fontsize=5)
        self._check_ticks_props(ax, xrot=20, xlabelsize=5, ylabelsize=5)

        axes = _check_plot_works(df.plot, filterwarnings='ignore', kind='kde', subplots=True)
        self._check_axes_shape(axes, axes_num=4, layout=(4, 1))

        axes = df.plot(kind='kde', logy=True, subplots=True)
        self._check_ax_scales(axes, yaxis='log')

    @slow
    def test_kde_missing_vals(self):
        tm._skip_if_no_scipy()
        _skip_if_no_scipy_gaussian_kde()
        df = DataFrame(np.random.uniform(size=(100, 4)))
        df.loc[0, 0] = np.nan
        ax = _check_plot_works(df.plot, kind='kde')

    @slow
    def test_hist_df(self):
        from matplotlib.patches import Rectangle
        if self.mpl_le_1_2_1:
            raise nose.SkipTest("not supported in matplotlib <= 1.2.x")

        df = DataFrame(randn(100, 4))
        series = df[0]

        ax = _check_plot_works(df.plot.hist)
        expected = [com.pprint_thing(c) for c in df.columns]
        self._check_legend_labels(ax, labels=expected)

        axes = _check_plot_works(df.plot.hist, filterwarnings='ignore', subplots=True, logy=True)
        self._check_axes_shape(axes, axes_num=4, layout=(4, 1))
        self._check_ax_scales(axes, yaxis='log')

        axes = series.plot.hist(rot=40)
        self._check_ticks_props(axes, xrot=40, yrot=0)
        tm.close()

        ax = series.plot.hist(normed=True, cumulative=True, bins=4)
        # height of last bin (index 5) must be 1.0
        rects = [x for x in ax.get_children() if isinstance(x, Rectangle)]
        self.assertAlmostEqual(rects[-1].get_height(), 1.0)
        tm.close()

        ax = series.plot.hist(cumulative=True, bins=4)
        rects = [x for x in ax.get_children() if isinstance(x, Rectangle)]

        self.assertAlmostEqual(rects[-2].get_height(), 100.0)
        tm.close()

        # if horizontal, yticklabels are rotated
        axes = df.plot.hist(rot=50, fontsize=8, orientation='horizontal')
        self._check_ticks_props(axes, xrot=0, yrot=50, ylabelsize=8)

    def _check_box_coord(self, patches, expected_y=None, expected_h=None,
                         expected_x=None, expected_w=None):
        result_y = np.array([p.get_y() for p in patches])
        result_height = np.array([p.get_height() for p in patches])
        result_x = np.array([p.get_x() for p in patches])
        result_width = np.array([p.get_width() for p in patches])

        if expected_y is not None:
            self.assert_numpy_array_equal(result_y, expected_y)
        if expected_h is not None:
            self.assert_numpy_array_equal(result_height, expected_h)
        if expected_x is not None:
            self.assert_numpy_array_equal(result_x, expected_x)
        if expected_w is not None:
            self.assert_numpy_array_equal(result_width, expected_w)

    @slow
    def test_hist_df_coord(self):
        normal_df = DataFrame({'A': np.repeat(np.array([1, 2, 3, 4, 5]),
                                              np.array([10, 9, 8, 7, 6])),
                               'B': np.repeat(np.array([1, 2, 3, 4, 5]),
                                              np.array([8, 8, 8, 8, 8])),
                               'C': np.repeat(np.array([1, 2, 3, 4, 5]),
                                              np.array([6, 7, 8, 9, 10]))},
                               columns=['A', 'B', 'C'])

        nan_df = DataFrame({'A': np.repeat(np.array([np.nan, 1, 2, 3, 4, 5]),
                                           np.array([3, 10, 9, 8, 7, 6])),
                            'B': np.repeat(np.array([1, np.nan, 2, 3, 4, 5]),
                                           np.array([8, 3, 8, 8, 8, 8])),
                            'C': np.repeat(np.array([1, 2, 3, np.nan, 4, 5]),
                                           np.array([6, 7, 8, 3, 9, 10]))},
                           columns=['A', 'B', 'C'])

        for df in [normal_df, nan_df]:
            ax = df.plot.hist(bins=5)
            self._check_box_coord(ax.patches[:5], expected_y=np.array([0, 0, 0, 0, 0]),
                                  expected_h=np.array([10, 9, 8, 7, 6]))
            self._check_box_coord(ax.patches[5:10], expected_y=np.array([0, 0, 0, 0, 0]),
                                  expected_h=np.array([8, 8, 8, 8, 8]))
            self._check_box_coord(ax.patches[10:], expected_y=np.array([0, 0, 0, 0, 0]),
                                  expected_h=np.array([6, 7, 8, 9, 10]))

            ax = df.plot.hist(bins=5, stacked=True)
            self._check_box_coord(ax.patches[:5], expected_y=np.array([0, 0, 0, 0, 0]),
                                  expected_h=np.array([10, 9, 8, 7, 6]))
            self._check_box_coord(ax.patches[5:10], expected_y=np.array([10, 9, 8, 7, 6]),
                                  expected_h=np.array([8, 8, 8, 8, 8]))
            self._check_box_coord(ax.patches[10:], expected_y=np.array([18, 17, 16, 15, 14]),
                                  expected_h=np.array([6, 7, 8, 9, 10]))

            axes = df.plot.hist(bins=5, stacked=True, subplots=True)
            self._check_box_coord(axes[0].patches, expected_y=np.array([0, 0, 0, 0, 0]),
                                  expected_h=np.array([10, 9, 8, 7, 6]))
            self._check_box_coord(axes[1].patches, expected_y=np.array([0, 0, 0, 0, 0]),
                                  expected_h=np.array([8, 8, 8, 8, 8]))
            self._check_box_coord(axes[2].patches, expected_y=np.array([0, 0, 0, 0, 0]),
                                  expected_h=np.array([6, 7, 8, 9, 10]))

            if self.mpl_ge_1_3_1:

                # horizontal
                ax = df.plot.hist(bins=5, orientation='horizontal')
                self._check_box_coord(ax.patches[:5], expected_x=np.array([0, 0, 0, 0, 0]),
                                      expected_w=np.array([10, 9, 8, 7, 6]))
                self._check_box_coord(ax.patches[5:10], expected_x=np.array([0, 0, 0, 0, 0]),
                                      expected_w=np.array([8, 8, 8, 8, 8]))
                self._check_box_coord(ax.patches[10:], expected_x=np.array([0, 0, 0, 0, 0]),
                                      expected_w=np.array([6, 7, 8, 9, 10]))

                ax = df.plot.hist(bins=5, stacked=True, orientation='horizontal')
                self._check_box_coord(ax.patches[:5], expected_x=np.array([0, 0, 0, 0, 0]),
                                      expected_w=np.array([10, 9, 8, 7, 6]))
                self._check_box_coord(ax.patches[5:10], expected_x=np.array([10, 9, 8, 7, 6]),
                                      expected_w=np.array([8, 8, 8, 8, 8]))
                self._check_box_coord(ax.patches[10:], expected_x=np.array([18, 17, 16, 15, 14]),
                                      expected_w=np.array([6, 7, 8, 9, 10]))

                axes = df.plot.hist(bins=5, stacked=True,
                               subplots=True, orientation='horizontal')
                self._check_box_coord(axes[0].patches, expected_x=np.array([0, 0, 0, 0, 0]),
                                      expected_w=np.array([10, 9, 8, 7, 6]))
                self._check_box_coord(axes[1].patches, expected_x=np.array([0, 0, 0, 0, 0]),
                                      expected_w=np.array([8, 8, 8, 8, 8]))
                self._check_box_coord(axes[2].patches, expected_x=np.array([0, 0, 0, 0, 0]),
                                      expected_w=np.array([6, 7, 8, 9, 10]))

    @slow
    def test_plot_int_columns(self):
        df = DataFrame(randn(100, 4)).cumsum()
        _check_plot_works(df.plot, legend=True)

    @slow
    def test_df_legend_labels(self):
        kinds = ['line', 'bar', 'barh', 'kde', 'area', 'hist']
        df = DataFrame(rand(3, 3), columns=['a', 'b', 'c'])
        df2 = DataFrame(rand(3, 3), columns=['d', 'e', 'f'])
        df3 = DataFrame(rand(3, 3), columns=['g', 'h', 'i'])
        df4 = DataFrame(rand(3, 3), columns=['j', 'k', 'l'])

        for kind in kinds:
            if not _ok_for_gaussian_kde(kind):
                continue

            ax = df.plot(kind=kind, legend=True)
            self._check_legend_labels(ax, labels=df.columns)

            ax = df2.plot(kind=kind, legend=False, ax=ax)
            self._check_legend_labels(ax, labels=df.columns)

            ax = df3.plot(kind=kind, legend=True, ax=ax)
            self._check_legend_labels(ax, labels=df.columns.union(df3.columns))

            ax = df4.plot(kind=kind, legend='reverse', ax=ax)
            expected = list(df.columns.union(df3.columns)) + list(reversed(df4.columns))
            self._check_legend_labels(ax, labels=expected)

        # Secondary Y
        ax = df.plot(legend=True, secondary_y='b')
        self._check_legend_labels(ax, labels=['a', 'b (right)', 'c'])
        ax = df2.plot(legend=False, ax=ax)
        self._check_legend_labels(ax, labels=['a', 'b (right)', 'c'])
        ax = df3.plot(kind='bar', legend=True, secondary_y='h', ax=ax)
        self._check_legend_labels(ax, labels=['a', 'b (right)', 'c', 'g', 'h (right)', 'i'])

        # Time Series
        ind = date_range('1/1/2014', periods=3)
        df = DataFrame(randn(3, 3), columns=['a', 'b', 'c'], index=ind)
        df2 = DataFrame(randn(3, 3), columns=['d', 'e', 'f'], index=ind)
        df3 = DataFrame(randn(3, 3), columns=['g', 'h', 'i'], index=ind)
        ax = df.plot(legend=True, secondary_y='b')
        self._check_legend_labels(ax, labels=['a', 'b (right)', 'c'])
        ax = df2.plot(legend=False, ax=ax)
        self._check_legend_labels(ax, labels=['a', 'b (right)', 'c'])
        ax = df3.plot(legend=True, ax=ax)
        self._check_legend_labels(ax, labels=['a', 'b (right)', 'c', 'g', 'h', 'i'])

        # scatter
        ax = df.plot.scatter(x='a', y='b', label='data1')
        self._check_legend_labels(ax, labels=['data1'])
        ax = df2.plot.scatter(x='d', y='e', legend=False,
                              label='data2', ax=ax)
        self._check_legend_labels(ax, labels=['data1'])
        ax = df3.plot.scatter(x='g', y='h', label='data3', ax=ax)
        self._check_legend_labels(ax, labels=['data1', 'data3'])

        # ensure label args pass through and
        # index name does not mutate
        # column names don't mutate
        df5 = df.set_index('a')
        ax = df5.plot(y='b')
        self._check_legend_labels(ax, labels=['b'])
        ax = df5.plot(y='b', label='LABEL_b')
        self._check_legend_labels(ax, labels=['LABEL_b'])
        self._check_text_labels(ax.xaxis.get_label(), 'a')
        ax = df5.plot(y='c', label='LABEL_c', ax=ax)
        self._check_legend_labels(ax, labels=['LABEL_b','LABEL_c'])
        self.assertTrue(df5.columns.tolist() == ['b','c'])


    def test_legend_name(self):
        multi = DataFrame(randn(4, 4),
                          columns=[np.array(['a', 'a', 'b', 'b']),
                                   np.array(['x', 'y', 'x', 'y'])])
        multi.columns.names = ['group', 'individual']

        ax = multi.plot()
        leg_title = ax.legend_.get_title()
        self._check_text_labels(leg_title, 'group,individual')

        df = DataFrame(randn(5, 5))
        ax = df.plot(legend=True, ax=ax)
        leg_title = ax.legend_.get_title()
        self._check_text_labels(leg_title, 'group,individual')

        df.columns.name = 'new'
        ax = df.plot(legend=False, ax=ax)
        leg_title = ax.legend_.get_title()
        self._check_text_labels(leg_title, 'group,individual')

        ax = df.plot(legend=True, ax=ax)
        leg_title = ax.legend_.get_title()
        self._check_text_labels(leg_title, 'new')

    @slow
    def test_no_legend(self):
        kinds = ['line', 'bar', 'barh', 'kde', 'area', 'hist']
        df = DataFrame(rand(3, 3), columns=['a', 'b', 'c'])

        for kind in kinds:
            if not _ok_for_gaussian_kde(kind):
                continue

            ax = df.plot(kind=kind, legend=False)
            self._check_legend_labels(ax, visible=False)

    @slow
    def test_style_by_column(self):
        import matplotlib.pyplot as plt
        fig = plt.gcf()

        df = DataFrame(randn(100, 3))
        for markers in [{0: '^', 1: '+', 2: 'o'},
                        {0: '^', 1: '+'},
                        ['^', '+', 'o'],
                        ['^', '+']]:
            fig.clf()
            fig.add_subplot(111)
            ax = df.plot(style=markers)
            for i, l in enumerate(ax.get_lines()[:len(markers)]):
                self.assertEqual(l.get_marker(), markers[i])

    @slow
    def test_line_label_none(self):
        s = Series([1, 2])
        ax = s.plot()
        self.assertEqual(ax.get_legend(), None)

        ax = s.plot(legend=True)
        self.assertEqual(ax.get_legend().get_texts()[0].get_text(),
                         'None')

    @slow
    def test_line_colors(self):
        import sys
        from matplotlib import cm

        custom_colors = 'rgcby'
        df = DataFrame(randn(5, 5))

        ax = df.plot(color=custom_colors)
        self._check_colors(ax.get_lines(), linecolors=custom_colors)

        tmp = sys.stderr
        sys.stderr = StringIO()
        try:
            tm.close()
            ax2 = df.plot(colors=custom_colors)
            lines2 = ax2.get_lines()
            for l1, l2 in zip(ax.get_lines(), lines2):
                self.assertEqual(l1.get_color(), l2.get_color())
        finally:
            sys.stderr = tmp

        tm.close()

        ax = df.plot(colormap='jet')
        rgba_colors = lmap(cm.jet, np.linspace(0, 1, len(df)))
        self._check_colors(ax.get_lines(), linecolors=rgba_colors)
        tm.close()

        ax = df.plot(colormap=cm.jet)
        rgba_colors = lmap(cm.jet, np.linspace(0, 1, len(df)))
        self._check_colors(ax.get_lines(), linecolors=rgba_colors)
        tm.close()

        # make color a list if plotting one column frame
        # handles cases like df.plot(color='DodgerBlue')
        ax = df.ix[:, [0]].plot(color='DodgerBlue')
        self._check_colors(ax.lines, linecolors=['DodgerBlue'])

        ax = df.plot(color='red')
        self._check_colors(ax.get_lines(), linecolors=['red'] * 5)
        tm.close()

        # GH 10299
        custom_colors = ['#FF0000', '#0000FF', '#FFFF00', '#000000', '#FFFFFF']
        ax = df.plot(color=custom_colors)
        self._check_colors(ax.get_lines(), linecolors=custom_colors)
        tm.close()

        with tm.assertRaises(ValueError):
            # Color contains shorthand hex value results in ValueError
            custom_colors = ['#F00', '#00F', '#FF0', '#000', '#FFF']
            # Forced show plot
            _check_plot_works(df.plot, color=custom_colors)

    @slow
    def test_line_colors_and_styles_subplots(self):
        # GH 9894
        from matplotlib import cm
        default_colors = self._maybe_unpack_cycler(self.plt.rcParams)

        df = DataFrame(randn(5, 5))

        axes = df.plot(subplots=True)
        for ax, c in zip(axes, list(default_colors)):
            self._check_colors(ax.get_lines(), linecolors=c)
        tm.close()

        # single color char
        axes = df.plot(subplots=True, color='k')
        for ax in axes:
            self._check_colors(ax.get_lines(), linecolors=['k'])
        tm.close()

        # single color str
        axes = df.plot(subplots=True, color='green')
        for ax in axes:
            self._check_colors(ax.get_lines(), linecolors=['green'])
        tm.close()

        custom_colors = 'rgcby'
        axes = df.plot(color=custom_colors, subplots=True)
        for ax, c in zip(axes, list(custom_colors)):
            self._check_colors(ax.get_lines(), linecolors=[c])
        tm.close()

        axes = df.plot(color=list(custom_colors), subplots=True)
        for ax, c in zip(axes, list(custom_colors)):
            self._check_colors(ax.get_lines(), linecolors=[c])
        tm.close()

        # GH 10299
        custom_colors = ['#FF0000', '#0000FF', '#FFFF00', '#000000', '#FFFFFF']
        axes = df.plot(color=custom_colors, subplots=True)
        for ax, c in zip(axes, list(custom_colors)):
            self._check_colors(ax.get_lines(), linecolors=[c])
        tm.close()

        with tm.assertRaises(ValueError):
            # Color contains shorthand hex value results in ValueError
            custom_colors = ['#F00', '#00F', '#FF0', '#000', '#FFF']
            # Forced show plot
            _check_plot_works(df.plot, color=custom_colors, subplots=True,
                              filterwarnings='ignore')

        rgba_colors = lmap(cm.jet, np.linspace(0, 1, len(df)))
        for cmap in ['jet', cm.jet]:
            axes = df.plot(colormap=cmap, subplots=True)
            for ax, c in zip(axes, rgba_colors):
                self._check_colors(ax.get_lines(), linecolors=[c])
            tm.close()

        # make color a list if plotting one column frame
        # handles cases like df.plot(color='DodgerBlue')
        axes = df.ix[:, [0]].plot(color='DodgerBlue', subplots=True)
        self._check_colors(axes[0].lines, linecolors=['DodgerBlue'])

        # single character style
        axes = df.plot(style='r', subplots=True)
        for ax in axes:
            self._check_colors(ax.get_lines(), linecolors=['r'])
        tm.close()

        # list of styles
        styles = list('rgcby')
        axes = df.plot(style=styles, subplots=True)
        for ax, c in zip(axes, styles):
            self._check_colors(ax.get_lines(), linecolors=[c])
        tm.close()

    @slow
    def test_area_colors(self):
        from matplotlib import cm
        from matplotlib.collections import PolyCollection

        custom_colors = 'rgcby'
        df = DataFrame(rand(5, 5))

        ax = df.plot.area(color=custom_colors)
        self._check_colors(ax.get_lines(), linecolors=custom_colors)
        poly = [o for o in ax.get_children() if isinstance(o, PolyCollection)]
        self._check_colors(poly, facecolors=custom_colors)

        handles, labels = ax.get_legend_handles_labels()
        # legend is stored as Line2D, thus check linecolors
        linehandles = [x for x in handles if not isinstance(x, PolyCollection)]
        self._check_colors(linehandles, linecolors=custom_colors)
        for h in handles:
            self.assertTrue(h.get_alpha() is None)
        tm.close()

        ax = df.plot.area(colormap='jet')
        jet_colors = lmap(cm.jet, np.linspace(0, 1, len(df)))
        self._check_colors(ax.get_lines(), linecolors=jet_colors)
        poly = [o for o in ax.get_children() if isinstance(o, PolyCollection)]
        self._check_colors(poly, facecolors=jet_colors)

        handles, labels = ax.get_legend_handles_labels()
        linehandles = [x for x in handles if not isinstance(x, PolyCollection)]
        self._check_colors(linehandles, linecolors=jet_colors)
        for h in handles:
            self.assertTrue(h.get_alpha() is None)
        tm.close()

        # When stacked=False, alpha is set to 0.5
        ax = df.plot.area(colormap=cm.jet, stacked=False)
        self._check_colors(ax.get_lines(), linecolors=jet_colors)
        poly = [o for o in ax.get_children() if isinstance(o, PolyCollection)]
        jet_with_alpha = [(c[0], c[1], c[2], 0.5) for c in jet_colors]
        self._check_colors(poly, facecolors=jet_with_alpha)

        handles, labels = ax.get_legend_handles_labels()
        # Line2D can't have alpha in its linecolor
        self._check_colors(handles[:len(jet_colors)], linecolors=jet_colors)
        for h in handles:
            self.assertEqual(h.get_alpha(), 0.5)

    @slow
    def test_hist_colors(self):
        default_colors = self._maybe_unpack_cycler(self.plt.rcParams)

        df = DataFrame(randn(5, 5))
        ax = df.plot.hist()
        self._check_colors(ax.patches[::10], facecolors=default_colors[:5])
        tm.close()

        custom_colors = 'rgcby'
        ax = df.plot.hist( color=custom_colors)
        self._check_colors(ax.patches[::10], facecolors=custom_colors)
        tm.close()

        from matplotlib import cm
        # Test str -> colormap functionality
        ax = df.plot.hist( colormap='jet')
        rgba_colors = lmap(cm.jet, np.linspace(0, 1, 5))
        self._check_colors(ax.patches[::10], facecolors=rgba_colors)
        tm.close()

        # Test colormap functionality
        ax = df.plot.hist( colormap=cm.jet)
        rgba_colors = lmap(cm.jet, np.linspace(0, 1, 5))
        self._check_colors(ax.patches[::10], facecolors=rgba_colors)
        tm.close()

        ax = df.ix[:, [0]].plot.hist(color='DodgerBlue')
        self._check_colors([ax.patches[0]], facecolors=['DodgerBlue'])

        ax = df.plot(kind='hist', color='green')
        self._check_colors(ax.patches[::10], facecolors=['green'] * 5)
        tm.close()

    @slow
    def test_kde_colors(self):
        tm._skip_if_no_scipy()
        _skip_if_no_scipy_gaussian_kde()

        from matplotlib import cm

        custom_colors = 'rgcby'
        df = DataFrame(rand(5, 5))

        ax = df.plot.kde(color=custom_colors)
        self._check_colors(ax.get_lines(), linecolors=custom_colors)
        tm.close()

        ax = df.plot.kde(colormap='jet')
        rgba_colors = lmap(cm.jet, np.linspace(0, 1, len(df)))
        self._check_colors(ax.get_lines(), linecolors=rgba_colors)
        tm.close()

        ax = df.plot.kde(colormap=cm.jet)
        rgba_colors = lmap(cm.jet, np.linspace(0, 1, len(df)))
        self._check_colors(ax.get_lines(), linecolors=rgba_colors)

    @slow
    def test_kde_colors_and_styles_subplots(self):
        tm._skip_if_no_scipy()
        _skip_if_no_scipy_gaussian_kde()

        from matplotlib import cm
        default_colors = self._maybe_unpack_cycler(self.plt.rcParams)

        df = DataFrame(randn(5, 5))

        axes = df.plot(kind='kde', subplots=True)
        for ax, c in zip(axes, list(default_colors)):
            self._check_colors(ax.get_lines(), linecolors=[c])
        tm.close()

        # single color char
        axes = df.plot(kind='kde', color='k', subplots=True)
        for ax in axes:
            self._check_colors(ax.get_lines(), linecolors=['k'])
        tm.close()

        # single color str
        axes = df.plot(kind='kde', color='red', subplots=True)
        for ax in axes:
            self._check_colors(ax.get_lines(), linecolors=['red'])
        tm.close()

        custom_colors = 'rgcby'
        axes = df.plot(kind='kde', color=custom_colors, subplots=True)
        for ax, c in zip(axes, list(custom_colors)):
            self._check_colors(ax.get_lines(), linecolors=[c])
        tm.close()

        rgba_colors = lmap(cm.jet, np.linspace(0, 1, len(df)))
        for cmap in ['jet', cm.jet]:
            axes = df.plot(kind='kde', colormap=cmap, subplots=True)
            for ax, c in zip(axes, rgba_colors):
                self._check_colors(ax.get_lines(), linecolors=[c])
            tm.close()

        # make color a list if plotting one column frame
        # handles cases like df.plot(color='DodgerBlue')
        axes = df.ix[:, [0]].plot(kind='kde', color='DodgerBlue', subplots=True)
        self._check_colors(axes[0].lines, linecolors=['DodgerBlue'])

        # single character style
        axes = df.plot(kind='kde', style='r', subplots=True)
        for ax in axes:
            self._check_colors(ax.get_lines(), linecolors=['r'])
        tm.close()

        # list of styles
        styles = list('rgcby')
        axes = df.plot(kind='kde', style=styles, subplots=True)
        for ax, c in zip(axes, styles):
            self._check_colors(ax.get_lines(), linecolors=[c])
        tm.close()

    @slow
    def test_boxplot_colors(self):

        def _check_colors(bp, box_c, whiskers_c, medians_c, caps_c='k', fliers_c='b'):
            self._check_colors(bp['boxes'], linecolors=[box_c] * len(bp['boxes']))
            self._check_colors(bp['whiskers'], linecolors=[whiskers_c] * len(bp['whiskers']))
            self._check_colors(bp['medians'], linecolors=[medians_c] * len(bp['medians']))
            self._check_colors(bp['fliers'], linecolors=[fliers_c] * len(bp['fliers']))
            self._check_colors(bp['caps'], linecolors=[caps_c] * len(bp['caps']))

        default_colors = self._maybe_unpack_cycler(self.plt.rcParams)

        df = DataFrame(randn(5, 5))
        bp = df.plot.box(return_type='dict')
        _check_colors(bp, default_colors[0], default_colors[0], default_colors[2])
        tm.close()

        dict_colors = dict(boxes='#572923', whiskers='#982042',
                           medians='#804823', caps='#123456')
        bp = df.plot.box(color=dict_colors, sym='r+', return_type='dict')
        _check_colors(bp, dict_colors['boxes'], dict_colors['whiskers'],
                      dict_colors['medians'], dict_colors['caps'], 'r')
        tm.close()

        # partial colors
        dict_colors = dict(whiskers='c', medians='m')
        bp = df.plot.box(color=dict_colors, return_type='dict')
        _check_colors(bp, default_colors[0], 'c', 'm')
        tm.close()

        from matplotlib import cm
        # Test str -> colormap functionality
        bp = df.plot.box(colormap='jet', return_type='dict')
        jet_colors = lmap(cm.jet, np.linspace(0, 1, 3))
        _check_colors(bp, jet_colors[0], jet_colors[0], jet_colors[2])
        tm.close()

        # Test colormap functionality
        bp = df.plot.box(colormap=cm.jet, return_type='dict')
        _check_colors(bp, jet_colors[0], jet_colors[0], jet_colors[2])
        tm.close()

        # string color is applied to all artists except fliers
        bp = df.plot.box(color='DodgerBlue', return_type='dict')
        _check_colors(bp, 'DodgerBlue', 'DodgerBlue', 'DodgerBlue',
                      'DodgerBlue')

        # tuple is also applied to all artists except fliers
        bp = df.plot.box(color=(0, 1, 0), sym='#123456', return_type='dict')
        _check_colors(bp, (0, 1, 0), (0, 1, 0), (0, 1, 0), (0, 1, 0), '#123456')

        with tm.assertRaises(ValueError):
            # Color contains invalid key results in ValueError
            df.plot.box(color=dict(boxes='red', xxxx='blue'))

    def test_default_color_cycle(self):
        import matplotlib.pyplot as plt
        colors = list('rgbk')
        if self.mpl_ge_1_5_0:
            import cycler
            plt.rcParams['axes.prop_cycle'] = cycler.cycler('color', colors)
        else:
            plt.rcParams['axes.color_cycle'] = colors

        df = DataFrame(randn(5, 3))
        ax = df.plot()

        expected = self._maybe_unpack_cycler(plt.rcParams)[:3]
        self._check_colors(ax.get_lines(), linecolors=expected)

    def test_unordered_ts(self):
        df = DataFrame(np.array([3.0, 2.0, 1.0]),
                       index=[date(2012, 10, 1),
                              date(2012, 9, 1),
                              date(2012, 8, 1)],
                       columns=['test'])
        ax = df.plot()
        xticks = ax.lines[0].get_xdata()
        self.assertTrue(xticks[0] < xticks[1])
        ydata = ax.lines[0].get_ydata()
        tm.assert_numpy_array_equal(ydata, np.array([1.0, 2.0, 3.0]))

    def test_kind_both_ways(self):
        df = DataFrame({'x': [1, 2, 3]})
        for kind in plotting._common_kinds:
            if not _ok_for_gaussian_kde(kind):
                continue
            df.plot(kind=kind)
            getattr(df.plot, kind)()
        for kind in ['scatter', 'hexbin']:
            df.plot('x', 'x', kind=kind)
            getattr(df.plot, kind)('x', 'x')

    def test_all_invalid_plot_data(self):
        df = DataFrame(list('abcd'))
        for kind in plotting._common_kinds:
            if not _ok_for_gaussian_kde(kind):
                continue
            with tm.assertRaises(TypeError):
                df.plot(kind=kind)

    @slow
    def test_partially_invalid_plot_data(self):
        with tm.RNGContext(42):
            df = DataFrame(randn(10, 2), dtype=object)
            df[np.random.rand(df.shape[0]) > 0.5] = 'a'
            for kind in plotting._common_kinds:
                if not _ok_for_gaussian_kde(kind):
                    continue
                with tm.assertRaises(TypeError):
                    df.plot(kind=kind)

        with tm.RNGContext(42):
            # area plot doesn't support positive/negative mixed data
            kinds = ['area']
            df = DataFrame(rand(10, 2), dtype=object)
            df[np.random.rand(df.shape[0]) > 0.5] = 'a'
            for kind in kinds:
                with tm.assertRaises(TypeError):
                    df.plot(kind=kind)

    def test_invalid_kind(self):
        df = DataFrame(randn(10, 2))
        with tm.assertRaises(ValueError):
            df.plot(kind='aasdf')

    @slow
    def test_hexbin_basic(self):
        df = self.hexbin_df

        ax = df.plot.hexbin(x='A', y='B', gridsize=10)
        # TODO: need better way to test. This just does existence.
        self.assertEqual(len(ax.collections), 1)

        # GH 6951
        axes = df.plot.hexbin(x='A', y='B', subplots=True)
        # hexbin should have 2 axes in the figure, 1 for plotting and another is colorbar
        self.assertEqual(len(axes[0].figure.axes), 2)
        # return value is single axes
        self._check_axes_shape(axes, axes_num=1, layout=(1, 1))

    @slow
    def test_hexbin_with_c(self):
        df = self.hexbin_df

        ax = df.plot.hexbin(x='A', y='B', C='C')
        self.assertEqual(len(ax.collections), 1)

        ax = df.plot.hexbin(x='A', y='B', C='C', reduce_C_function=np.std)
        self.assertEqual(len(ax.collections), 1)

    @slow
    def test_hexbin_cmap(self):
        df = self.hexbin_df

        # Default to BuGn
        ax = df.plot.hexbin(x='A', y='B')
        self.assertEqual(ax.collections[0].cmap.name, 'BuGn')

        cm = 'cubehelix'
        ax = df.plot.hexbin(x='A', y='B', colormap=cm)
        self.assertEqual(ax.collections[0].cmap.name, cm)

    @slow
    def test_no_color_bar(self):
        df = self.hexbin_df

        ax = df.plot.hexbin(x='A', y='B', colorbar=None)
        self.assertIs(ax.collections[0].colorbar, None)

    @slow
    def test_allow_cmap(self):
        df = self.hexbin_df

        ax = df.plot.hexbin(x='A', y='B', cmap='YlGn')
        self.assertEqual(ax.collections[0].cmap.name, 'YlGn')

        with tm.assertRaises(TypeError):
            df.plot.hexbin(x='A', y='B', cmap='YlGn',
                           colormap='BuGn')

    @slow
    def test_pie_df(self):
        df = DataFrame(np.random.rand(5, 3), columns=['X', 'Y', 'Z'],
                       index=['a', 'b', 'c', 'd', 'e'])
        with tm.assertRaises(ValueError):
            df.plot.pie()

        ax = _check_plot_works(df.plot.pie, y='Y')
        self._check_text_labels(ax.texts, df.index)

        ax = _check_plot_works(df.plot.pie, y=2)
        self._check_text_labels(ax.texts, df.index)

        axes = _check_plot_works(df.plot.pie, filterwarnings='ignore', subplots=True)
        self.assertEqual(len(axes), len(df.columns))
        for ax in axes:
            self._check_text_labels(ax.texts, df.index)
        for ax, ylabel in zip(axes, df.columns):
            self.assertEqual(ax.get_ylabel(), ylabel)

        labels = ['A', 'B', 'C', 'D', 'E']
        color_args = ['r', 'g', 'b', 'c', 'm']
        axes = _check_plot_works(df.plot.pie, filterwarnings='ignore', subplots=True,
                                 labels=labels, colors=color_args)
        self.assertEqual(len(axes), len(df.columns))

        for ax in axes:
            self._check_text_labels(ax.texts, labels)
            self._check_colors(ax.patches, facecolors=color_args)

    def test_pie_df_nan(self):
        df = DataFrame(np.random.rand(4, 4))
        for i in range(4):
            df.iloc[i, i] = np.nan
        fig, axes = self.plt.subplots(ncols=4)
        df.plot.pie(subplots=True, ax=axes, legend=True)

        base_expected = ['0', '1', '2', '3']
        for i, ax in enumerate(axes):
            expected = list(base_expected)  # force copy
            expected[i] = ''
            result = [x.get_text() for x in ax.texts]
            self.assertEqual(result, expected)
            # legend labels
            # NaN's not included in legend with subplots
            # see https://github.com/pydata/pandas/issues/8390
            self.assertEqual([x.get_text() for x in
                              ax.get_legend().get_texts()],
                             base_expected[:i] + base_expected[i+1:])

    @slow
    def test_errorbar_plot(self):
        d = {'x': np.arange(12), 'y': np.arange(12, 0, -1)}
        df = DataFrame(d)
        d_err = {'x': np.ones(12)*0.2, 'y': np.ones(12)*0.4}
        df_err = DataFrame(d_err)

        # check line plots
        ax = _check_plot_works(df.plot, yerr=df_err, logy=True)
        self._check_has_errorbars(ax, xerr=0, yerr=2)
        ax = _check_plot_works(df.plot, yerr=df_err, logx=True, logy=True)
        self._check_has_errorbars(ax, xerr=0, yerr=2)
        ax = _check_plot_works(df.plot, yerr=df_err, loglog=True)
        self._check_has_errorbars(ax, xerr=0, yerr=2)

        kinds = ['line', 'bar', 'barh']
        for kind in kinds:
            ax = _check_plot_works(df.plot, yerr=df_err['x'], kind=kind)
            self._check_has_errorbars(ax, xerr=0, yerr=2)
            ax = _check_plot_works(df.plot, yerr=d_err, kind=kind)
            self._check_has_errorbars(ax, xerr=0, yerr=2)
            ax = _check_plot_works(df.plot, yerr=df_err, xerr=df_err, kind=kind)
            self._check_has_errorbars(ax, xerr=2, yerr=2)
            ax = _check_plot_works(df.plot, yerr=df_err['x'], xerr=df_err['x'], kind=kind)
            self._check_has_errorbars(ax, xerr=2, yerr=2)
            ax = _check_plot_works(df.plot, xerr=0.2, yerr=0.2, kind=kind)
            self._check_has_errorbars(ax, xerr=2, yerr=2)
            axes = _check_plot_works(df.plot, filterwarnings='ignore', yerr=df_err,
                                     xerr=df_err, subplots=True, kind=kind)
            self._check_has_errorbars(axes, xerr=1, yerr=1)

        ax = _check_plot_works((df+1).plot, yerr=df_err, xerr=df_err, kind='bar', log=True)
        self._check_has_errorbars(ax, xerr=2, yerr=2)

        # yerr is raw error values
        ax = _check_plot_works(df['y'].plot, yerr=np.ones(12)*0.4)
        self._check_has_errorbars(ax, xerr=0, yerr=1)
        ax = _check_plot_works(df.plot, yerr=np.ones((2, 12))*0.4)
        self._check_has_errorbars(ax, xerr=0, yerr=2)

        # yerr is iterator
        import itertools
        ax = _check_plot_works(df.plot, yerr=itertools.repeat(0.1, len(df)))
        self._check_has_errorbars(ax, xerr=0, yerr=2)

        # yerr is column name
        for yerr in ['yerr', u('誤差')]:
            s_df = df.copy()
            s_df[yerr] = np.ones(12)*0.2
            ax = _check_plot_works(s_df.plot, yerr=yerr)
            self._check_has_errorbars(ax, xerr=0, yerr=2)
            ax = _check_plot_works(s_df.plot, y='y', x='x', yerr=yerr)
            self._check_has_errorbars(ax, xerr=0, yerr=1)

        with tm.assertRaises(ValueError):
            df.plot(yerr=np.random.randn(11))

        df_err = DataFrame({'x': ['zzz']*12, 'y': ['zzz']*12})
        with tm.assertRaises((ValueError, TypeError)):
            df.plot(yerr=df_err)

    @slow
    def test_errorbar_with_integer_column_names(self):
        # test with integer column names
        df = DataFrame(np.random.randn(10, 2))
        df_err = DataFrame(np.random.randn(10, 2))
        ax = _check_plot_works(df.plot, yerr=df_err)
        self._check_has_errorbars(ax, xerr=0, yerr=2)
        ax = _check_plot_works(df.plot, y=0, yerr=1)
        self._check_has_errorbars(ax, xerr=0, yerr=1)

    @slow
    def test_errorbar_with_partial_columns(self):
        df = DataFrame(np.random.randn(10, 3))
        df_err = DataFrame(np.random.randn(10, 2), columns=[0, 2])
        kinds = ['line', 'bar']
        for kind in kinds:
            ax = _check_plot_works(df.plot, yerr=df_err, kind=kind)
            self._check_has_errorbars(ax, xerr=0, yerr=2)

        ix = date_range('1/1/2000', periods=10, freq='M')
        df.set_index(ix, inplace=True)
        df_err.set_index(ix, inplace=True)
        ax = _check_plot_works(df.plot, yerr=df_err, kind='line')
        self._check_has_errorbars(ax, xerr=0, yerr=2)

        d = {'x': np.arange(12), 'y': np.arange(12, 0, -1)}
        df = DataFrame(d)
        d_err = {'x': np.ones(12)*0.2, 'z': np.ones(12)*0.4}
        df_err = DataFrame(d_err)
        for err in [d_err, df_err]:
            ax = _check_plot_works(df.plot, yerr=err)
            self._check_has_errorbars(ax, xerr=0, yerr=1)

    @slow
    def test_errorbar_timeseries(self):

        d = {'x': np.arange(12), 'y': np.arange(12, 0, -1)}
        d_err = {'x': np.ones(12)*0.2, 'y': np.ones(12)*0.4}

        # check time-series plots
        ix = date_range('1/1/2000', '1/1/2001', freq='M')
        tdf = DataFrame(d, index=ix)
        tdf_err = DataFrame(d_err, index=ix)

        kinds = ['line', 'bar', 'barh']
        for kind in kinds:
            ax = _check_plot_works(tdf.plot, yerr=tdf_err, kind=kind)
            self._check_has_errorbars(ax, xerr=0, yerr=2)
            ax = _check_plot_works(tdf.plot, yerr=d_err, kind=kind)
            self._check_has_errorbars(ax, xerr=0, yerr=2)
            ax = _check_plot_works(tdf.plot, y='y', yerr=tdf_err['x'], kind=kind)
            self._check_has_errorbars(ax, xerr=0, yerr=1)
            ax = _check_plot_works(tdf.plot, y='y', yerr='x', kind=kind)
            self._check_has_errorbars(ax, xerr=0, yerr=1)
            ax = _check_plot_works(tdf.plot, yerr=tdf_err, kind=kind)
            self._check_has_errorbars(ax, xerr=0, yerr=2)
            axes = _check_plot_works(tdf.plot, filterwarnings='ignore', kind=kind,
                                     yerr=tdf_err, subplots=True)
            self._check_has_errorbars(axes, xerr=0, yerr=1)

    def test_errorbar_asymmetrical(self):

        np.random.seed(0)
        err = np.random.rand(3, 2, 5)

        data = np.random.randn(5, 3)
        df = DataFrame(data)

        ax = df.plot(yerr=err, xerr=err/2)

        self.assertEqual(ax.lines[7].get_ydata()[0], data[0,1]-err[1,0,0])
        self.assertEqual(ax.lines[8].get_ydata()[0], data[0,1]+err[1,1,0])

        self.assertEqual(ax.lines[5].get_xdata()[0], -err[1,0,0]/2)
        self.assertEqual(ax.lines[6].get_xdata()[0], err[1,1,0]/2)

        with tm.assertRaises(ValueError):
            df.plot(yerr=err.T)

        tm.close()

    def test_table(self):
        df = DataFrame(np.random.rand(10, 3),
                       index=list(string.ascii_letters[:10]))
        _check_plot_works(df.plot, table=True)
        _check_plot_works(df.plot, table=df)

        ax = df.plot()
        self.assertTrue(len(ax.tables) == 0)
        plotting.table(ax, df.T)
        self.assertTrue(len(ax.tables) == 1)

    def test_errorbar_scatter(self):
        df = DataFrame(np.random.randn(5, 2), index=range(5), columns=['x', 'y'])
        df_err = DataFrame(np.random.randn(5, 2) / 5,
                           index=range(5), columns=['x', 'y'])

        ax = _check_plot_works(df.plot.scatter, x='x', y='y')
        self._check_has_errorbars(ax, xerr=0, yerr=0)
        ax = _check_plot_works(df.plot.scatter, x='x', y='y', xerr=df_err)
        self._check_has_errorbars(ax, xerr=1, yerr=0)

        ax = _check_plot_works(df.plot.scatter, x='x', y='y', yerr=df_err)
        self._check_has_errorbars(ax, xerr=0, yerr=1)
        ax = _check_plot_works(df.plot.scatter, x='x', y='y',
                               xerr=df_err, yerr=df_err)
        self._check_has_errorbars(ax, xerr=1, yerr=1)

        def _check_errorbar_color(containers, expected, has_err='has_xerr'):
            errs = [c.lines[1][0] for c in ax.containers if getattr(c, has_err, False)]
            self._check_colors(errs, linecolors=[expected] * len(errs))

        # GH 8081
        df = DataFrame(np.random.randn(10, 5), columns=['a', 'b', 'c', 'd', 'e'])
        ax = df.plot.scatter(x='a', y='b', xerr='d', yerr='e', c='red')
        self._check_has_errorbars(ax, xerr=1, yerr=1)
        _check_errorbar_color(ax.containers, 'red', has_err='has_xerr')
        _check_errorbar_color(ax.containers, 'red', has_err='has_yerr')

        ax = df.plot.scatter(x='a', y='b', yerr='e', color='green')
        self._check_has_errorbars(ax, xerr=0, yerr=1)
        _check_errorbar_color(ax.containers, 'green', has_err='has_yerr')

    @slow
    def test_sharex_and_ax(self):
        # https://github.com/pydata/pandas/issues/9737
        # using gridspec, the axis in fig.get_axis() are sorted differently than pandas expected
        # them, so make sure that only the right ones are removed
        import matplotlib.pyplot as plt
        plt.close('all')
        gs, axes = _generate_4_axes_via_gridspec()

        df = DataFrame({"a": [1, 2, 3, 4, 5, 6],
                        "b": [1, 2, 3, 4, 5, 6],
                        "c": [1, 2, 3, 4, 5, 6],
                        "d": [1, 2, 3, 4, 5, 6]})

        def _check(axes):
            for ax in axes:
                self.assertEqual(len(ax.lines), 1)
                self._check_visible(ax.get_yticklabels(), visible=True)
            for ax in [axes[0], axes[2]]:
                self._check_visible(ax.get_xticklabels(), visible=False)
                self._check_visible(ax.get_xticklabels(minor=True), visible=False)
            for ax in [axes[1], axes[3]]:
                self._check_visible(ax.get_xticklabels(), visible=True)
                self._check_visible(ax.get_xticklabels(minor=True), visible=True)

        for ax in axes:
            df.plot(x="a", y="b", title="title", ax=ax, sharex=True)
        gs.tight_layout(plt.gcf())
        _check(axes)
        tm.close()

        gs, axes = _generate_4_axes_via_gridspec()
        with tm.assert_produces_warning(UserWarning):
            axes = df.plot(subplots=True, ax=axes, sharex=True)
        _check(axes)
        tm.close()

        gs, axes = _generate_4_axes_via_gridspec()
        # without sharex, no labels should be touched!
        for ax in axes:
            df.plot(x="a", y="b", title="title", ax=ax)

        gs.tight_layout(plt.gcf())
        for ax in axes:
            self.assertEqual(len(ax.lines), 1)
            self._check_visible(ax.get_yticklabels(), visible=True)
            self._check_visible(ax.get_xticklabels(), visible=True)
            self._check_visible(ax.get_xticklabels(minor=True), visible=True)
        tm.close()

    @slow
    def test_sharey_and_ax(self):
        # https://github.com/pydata/pandas/issues/9737
        # using gridspec, the axis in fig.get_axis() are sorted differently than pandas expected
        # them, so make sure that only the right ones are removed
        import matplotlib.pyplot as plt

        gs, axes = _generate_4_axes_via_gridspec()

        df = DataFrame({"a": [1, 2, 3, 4, 5, 6],
                        "b": [1, 2, 3, 4, 5, 6],
                        "c": [1, 2, 3, 4, 5, 6],
                        "d": [1, 2, 3, 4, 5, 6]})

        def _check(axes):
            for ax in axes:
                self.assertEqual(len(ax.lines), 1)
                self._check_visible(ax.get_xticklabels(), visible=True)
                self._check_visible(ax.get_xticklabels(minor=True), visible=True)
            for ax in [axes[0], axes[1]]:
                self._check_visible(ax.get_yticklabels(), visible=True)
            for ax in [axes[2], axes[3]]:
                self._check_visible(ax.get_yticklabels(), visible=False)

        for ax in axes:
            df.plot(x="a", y="b", title="title", ax=ax, sharey=True)
        gs.tight_layout(plt.gcf())
        _check(axes)
        tm.close()

        gs, axes = _generate_4_axes_via_gridspec()
        with tm.assert_produces_warning(UserWarning):
            axes = df.plot(subplots=True, ax=axes, sharey=True)

        gs.tight_layout(plt.gcf())
        _check(axes)
        tm.close()

        gs, axes = _generate_4_axes_via_gridspec()
        # without sharex, no labels should be touched!
        for ax in axes:
            df.plot(x="a", y="b", title="title", ax=ax)

        gs.tight_layout(plt.gcf())
        for ax in axes:
            self.assertEqual(len(ax.lines), 1)
            self._check_visible(ax.get_yticklabels(), visible=True)
            self._check_visible(ax.get_xticklabels(), visible=True)
            self._check_visible(ax.get_xticklabels(minor=True), visible=True)

    def test_memory_leak(self):
        """ Check that every plot type gets properly collected. """
        import weakref
        import gc

        results = {}
        for kind in plotting._plot_klass.keys():
            if not _ok_for_gaussian_kde(kind):
                continue
            args = {}
            if kind in ['hexbin', 'scatter', 'pie']:
                df = self.hexbin_df
                args = {'x': 'A', 'y': 'B'}
            elif kind == 'area':
                df = self.tdf.abs()
            else:
                df = self.tdf

            # Use a weakref so we can see if the object gets collected without
            # also preventing it from being collected
            results[kind] = weakref.proxy(df.plot(kind=kind, **args))

        # have matplotlib delete all the figures
        tm.close()
        # force a garbage collection
        gc.collect()
        for key in results:
            # check that every plot was collected
            with tm.assertRaises(ReferenceError):
                # need to actually access something to get an error
                results[key].lines

    @slow
    def test_df_subplots_patterns_minorticks(self):
        # GH 10657
        import matplotlib.pyplot as plt

        df = DataFrame(np.random.randn(10, 2),
                       index=date_range('1/1/2000', periods=10),
                       columns=list('AB'))

        # shared subplots
        fig, axes = plt.subplots(2, 1, sharex=True)
        axes = df.plot(subplots=True, ax=axes)
        for ax in axes:
            self.assertEqual(len(ax.lines), 1)
            self._check_visible(ax.get_yticklabels(), visible=True)
        # xaxis of 1st ax must be hidden
        self._check_visible(axes[0].get_xticklabels(), visible=False)
        self._check_visible(axes[0].get_xticklabels(minor=True), visible=False)
        self._check_visible(axes[1].get_xticklabels(), visible=True)
        self._check_visible(axes[1].get_xticklabels(minor=True), visible=True)
        tm.close()

        fig, axes = plt.subplots(2, 1)
        with tm.assert_produces_warning(UserWarning):
            axes = df.plot(subplots=True, ax=axes, sharex=True)
        for ax in axes:
            self.assertEqual(len(ax.lines), 1)
            self._check_visible(ax.get_yticklabels(), visible=True)
        # xaxis of 1st ax must be hidden
        self._check_visible(axes[0].get_xticklabels(), visible=False)
        self._check_visible(axes[0].get_xticklabels(minor=True), visible=False)
        self._check_visible(axes[1].get_xticklabels(), visible=True)
        self._check_visible(axes[1].get_xticklabels(minor=True), visible=True)
        tm.close()

        # not shared
        fig, axes = plt.subplots(2, 1)
        axes = df.plot(subplots=True, ax=axes)
        for ax in axes:
            self.assertEqual(len(ax.lines), 1)
            self._check_visible(ax.get_yticklabels(), visible=True)
            self._check_visible(ax.get_xticklabels(), visible=True)
            self._check_visible(ax.get_xticklabels(minor=True), visible=True)
        tm.close()

    @slow
    def test_df_gridspec_patterns(self):
        # GH 10819
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        ts = Series(np.random.randn(10),
                    index=date_range('1/1/2000', periods=10))

        df = DataFrame(np.random.randn(10, 2), index=ts.index,
                       columns=list('AB'))

        def _get_vertical_grid():
            gs = gridspec.GridSpec(3, 1)
            fig = plt.figure()
            ax1 = fig.add_subplot(gs[:2, :])
            ax2 = fig.add_subplot(gs[2, :])
            return ax1, ax2

        def _get_horizontal_grid():
            gs = gridspec.GridSpec(1, 3)
            fig = plt.figure()
            ax1 = fig.add_subplot(gs[:, :2])
            ax2 = fig.add_subplot(gs[:, 2])
            return ax1, ax2

        for ax1, ax2 in [_get_vertical_grid(), _get_horizontal_grid()]:
            ax1 = ts.plot(ax=ax1)
            self.assertEqual(len(ax1.lines), 1)
            ax2 = df.plot(ax=ax2)
            self.assertEqual(len(ax2.lines), 2)
            for ax in [ax1, ax2]:
                self._check_visible(ax.get_yticklabels(), visible=True)
                self._check_visible(ax.get_xticklabels(), visible=True)
                self._check_visible(ax.get_xticklabels(minor=True), visible=True)
            tm.close()

        # subplots=True
        for ax1, ax2 in [_get_vertical_grid(), _get_horizontal_grid()]:
            axes = df.plot(subplots=True, ax=[ax1, ax2])
            self.assertEqual(len(ax1.lines), 1)
            self.assertEqual(len(ax2.lines), 1)
            for ax in axes:
                self._check_visible(ax.get_yticklabels(), visible=True)
                self._check_visible(ax.get_xticklabels(), visible=True)
                self._check_visible(ax.get_xticklabels(minor=True), visible=True)
            tm.close()

        # vertical / subplots / sharex=True / sharey=True
        ax1, ax2 = _get_vertical_grid()
        with tm.assert_produces_warning(UserWarning):
            axes = df.plot(subplots=True, ax=[ax1, ax2],
                           sharex=True, sharey=True)
        self.assertEqual(len(axes[0].lines), 1)
        self.assertEqual(len(axes[1].lines), 1)
        for ax in [ax1, ax2]:
            # yaxis are visible because there is only one column
            self._check_visible(ax.get_yticklabels(), visible=True)
        # xaxis of axes0 (top) are hidden
        self._check_visible(axes[0].get_xticklabels(), visible=False)
        self._check_visible(axes[0].get_xticklabels(minor=True), visible=False)
        self._check_visible(axes[1].get_xticklabels(), visible=True)
        self._check_visible(axes[1].get_xticklabels(minor=True), visible=True)
        tm.close()

        # horizontal / subplots / sharex=True / sharey=True
        ax1, ax2 = _get_horizontal_grid()
        with tm.assert_produces_warning(UserWarning):
            axes = df.plot(subplots=True, ax=[ax1, ax2],
                           sharex=True, sharey=True)
        self.assertEqual(len(axes[0].lines), 1)
        self.assertEqual(len(axes[1].lines), 1)
        self._check_visible(axes[0].get_yticklabels(), visible=True)
        # yaxis of axes1 (right) are hidden
        self._check_visible(axes[1].get_yticklabels(), visible=False)
        for ax in [ax1, ax2]:
            # xaxis are visible because there is only one column
            self._check_visible(ax.get_xticklabels(), visible=True)
            self._check_visible(ax.get_xticklabels(minor=True), visible=True)
        tm.close()

        # boxed
        def _get_boxed_grid():
            gs = gridspec.GridSpec(3,3)
            fig = plt.figure()
            ax1 = fig.add_subplot(gs[:2, :2])
            ax2 = fig.add_subplot(gs[:2, 2])
            ax3 = fig.add_subplot(gs[2, :2])
            ax4 = fig.add_subplot(gs[2, 2])
            return ax1, ax2, ax3, ax4

        axes = _get_boxed_grid()
        df = DataFrame(np.random.randn(10, 4),
                      index=ts.index, columns=list('ABCD'))
        axes = df.plot(subplots=True, ax=axes)
        for ax in axes:
            self.assertEqual(len(ax.lines), 1)
            # axis are visible because these are not shared
            self._check_visible(ax.get_yticklabels(), visible=True)
            self._check_visible(ax.get_xticklabels(), visible=True)
            self._check_visible(ax.get_xticklabels(minor=True), visible=True)
        tm.close()

        # subplots / sharex=True / sharey=True
        axes = _get_boxed_grid()
        with tm.assert_produces_warning(UserWarning):
            axes = df.plot(subplots=True, ax=axes, sharex=True, sharey=True)
        for ax in axes:
            self.assertEqual(len(ax.lines), 1)
        for ax in [axes[0], axes[2]]: # left column
            self._check_visible(ax.get_yticklabels(), visible=True)
        for ax in [axes[1], axes[3]]: # right column
            self._check_visible(ax.get_yticklabels(), visible=False)
        for ax in [axes[0], axes[1]]: # top row
            self._check_visible(ax.get_xticklabels(), visible=False)
            self._check_visible(ax.get_xticklabels(minor=True), visible=False)
        for ax in [axes[2], axes[3]]: # bottom row
            self._check_visible(ax.get_xticklabels(), visible=True)
            self._check_visible(ax.get_xticklabels(minor=True), visible=True)
        tm.close()

    @slow
    def test_df_grid_settings(self):
        # Make sure plot defaults to rcParams['axes.grid'] setting, GH 9792
        self._check_grid_settings(DataFrame({'a':[1,2,3],'b':[2,3,4]}),
            plotting._dataframe_kinds, kws={'x':'a','y':'b'})

    def test_option_mpl_style(self):
        set_option('display.mpl_style', 'default')
        set_option('display.mpl_style', None)
        set_option('display.mpl_style', False)

        with tm.assertRaises(ValueError):
            set_option('display.mpl_style', 'default2')

    def test_invalid_colormap(self):
        df = DataFrame(randn(3, 2), columns=['A', 'B'])

        with tm.assertRaises(ValueError):
            df.plot(colormap='invalid_colormap')

    def test_plain_axes(self):

        # supplied ax itself is a SubplotAxes, but figure contains also
        # a plain Axes object (GH11556)
        fig, ax = self.plt.subplots()
        fig.add_axes([0.2, 0.2, 0.2, 0.2])
        Series(rand(10)).plot(ax=ax)

        # suppliad ax itself is a plain Axes, but because the cmap keyword
        # a new ax is created for the colorbar -> also multiples axes (GH11520)
        df = DataFrame({'a': randn(8), 'b': randn(8)})
        fig = self.plt.figure()
        ax = fig.add_axes((0,0,1,1))
        df.plot(kind='scatter', ax=ax, x='a', y='b', c='a', cmap='hsv')

        # other examples
        fig, ax = self.plt.subplots()
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        Series(rand(10)).plot(ax=ax)
        Series(rand(10)).plot(ax=cax)

        fig, ax = self.plt.subplots()
        from mpl_toolkits.axes_grid.inset_locator import inset_axes
        iax = inset_axes(ax, width="30%", height=1., loc=3)
        Series(rand(10)).plot(ax=ax)
        Series(rand(10)).plot(ax=iax)


@tm.mplskip
class TestDataFrameGroupByPlots(TestPlotBase):

    def test_series_groupby_plotting_nominally_works(self):
        n = 10
        weight = Series(np.random.normal(166, 20, size=n))
        height = Series(np.random.normal(60, 10, size=n))
        with tm.RNGContext(42):
            gender = tm.choice(['male', 'female'], size=n)

        weight.groupby(gender).plot()
        tm.close()
        height.groupby(gender).hist()
        tm.close()
        #Regression test for GH8733
        height.groupby(gender).plot(alpha=0.5)
        tm.close()

    def test_plotting_with_float_index_works(self):
        # GH 7025
        df = DataFrame({'def': [1,1,1,2,2,2,3,3,3],
                        'val': np.random.randn(9)},
                       index=[1.0,2.0,3.0,1.0,2.0,3.0,1.0,2.0,3.0])

        df.groupby('def')['val'].plot()
        tm.close()
        df.groupby('def')['val'].apply(lambda x: x.plot())
        tm.close()

    def test_hist_single_row(self):
        # GH10214
        bins = np.arange(80, 100 + 2, 1)
        df = DataFrame({"Name": ["AAA", "BBB"], "ByCol": [1, 2], "Mark": [85, 89]})
        df["Mark"].hist(by=df["ByCol"], bins=bins)
        df = DataFrame({"Name": ["AAA"], "ByCol": [1], "Mark": [85]})
        df["Mark"].hist(by=df["ByCol"], bins=bins)

    def test_plot_submethod_works(self):
        df = DataFrame({'x': [1, 2, 3, 4, 5],
                        'y': [1, 2, 3, 2, 1],
                        'z': list('ababa')})
        df.groupby('z').plot.scatter('x', 'y')
        tm.close()
        df.groupby('z')['x'].plot.line()
        tm.close()


def assert_is_valid_plot_return_object(objs):
    import matplotlib.pyplot as plt
    if isinstance(objs, np.ndarray):
        for el in objs.flat:
            assert isinstance(el, plt.Axes), ('one of \'objs\' is not a '
                                              'matplotlib Axes instance, '
                                              'type encountered {0!r}'
                                              ''.format(el.__class__.__name__))
    else:
        assert isinstance(objs, (plt.Artist, tuple, dict)), \
                ('objs is neither an ndarray of Artist instances nor a '
                 'single Artist instance, tuple, or dict, "objs" is a {0!r} '
                 ''.format(objs.__class__.__name__))


def _check_plot_works(f, filterwarnings='always', **kwargs):
    import matplotlib.pyplot as plt
    ret = None
    with warnings.catch_warnings():
        warnings.simplefilter(filterwarnings)
        try:
            try:
                fig = kwargs['figure']
            except KeyError:
                fig = plt.gcf()

            plt.clf()

            ax = kwargs.get('ax', fig.add_subplot(211))
            ret = f(**kwargs)

            assert_is_valid_plot_return_object(ret)

            try:
                kwargs['ax'] = fig.add_subplot(212)
                ret = f(**kwargs)
            except Exception:
                pass
            else:
                assert_is_valid_plot_return_object(ret)

            with ensure_clean(return_filelike=True) as path:
                plt.savefig(path)
        finally:
            tm.close(fig)

        return ret

def _generate_4_axes_via_gridspec():
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import matplotlib.gridspec

    gs = mpl.gridspec.GridSpec(2, 2)
    ax_tl = plt.subplot(gs[0,0])
    ax_ll = plt.subplot(gs[1,0])
    ax_tr = plt.subplot(gs[0,1])
    ax_lr = plt.subplot(gs[1,1])

    return gs, [ax_tl, ax_ll, ax_tr, ax_lr]


def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
