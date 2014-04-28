#!/usr/bin/env python
# coding: utf-8

import nose
import itertools
import os
import string
from distutils.version import LooseVersion

from datetime import datetime, date, timedelta

from pandas import Series, DataFrame, MultiIndex, PeriodIndex, date_range
from pandas.compat import range, lrange, StringIO, lmap, lzip, u, zip
import pandas.util.testing as tm
from pandas.util.testing import ensure_clean
from pandas.core.config import set_option


import numpy as np
from numpy import random
from numpy.random import rand, randn

from numpy.testing import assert_array_equal
from numpy.testing.decorators import slow
import pandas.tools.plotting as plotting


def _skip_if_no_scipy():
    try:
        import scipy
    except ImportError:
        raise nose.SkipTest("no scipy")

@tm.mplskip
class TestSeriesPlots(tm.TestCase):
    def setUp(self):
        import matplotlib as mpl
        self.mpl_le_1_2_1 = str(mpl.__version__) <= LooseVersion('1.2.1')
        self.ts = tm.makeTimeSeries()
        self.ts.name = 'ts'

        self.series = tm.makeStringSeries()
        self.series.name = 'series'

        self.iseries = tm.makePeriodSeries()
        self.iseries.name = 'iseries'

    def tearDown(self):
        tm.close()

    @slow
    def test_plot(self):
        _check_plot_works(self.ts.plot, label='foo')
        _check_plot_works(self.ts.plot, use_index=False)
        _check_plot_works(self.ts.plot, rot=0)
        _check_plot_works(self.ts.plot, style='.', logy=True)
        _check_plot_works(self.ts.plot, style='.', logx=True)
        _check_plot_works(self.ts.plot, style='.', loglog=True)
        _check_plot_works(self.ts[:10].plot, kind='bar')
        _check_plot_works(self.ts.plot, kind='area', stacked=False)
        _check_plot_works(self.iseries.plot)

        for kind in ['line', 'bar', 'barh', 'kde']:
            _check_plot_works(self.series[:5].plot, kind=kind)

        _check_plot_works(self.series[:10].plot, kind='barh')
        _check_plot_works(Series(randn(10)).plot, kind='bar', color='black')

        # GH 6951
        _check_plot_works(self.ts.plot, subplots=True)

    @slow
    def test_plot_figsize_and_title(self):
        # figsize and title
        import matplotlib.pyplot as plt
        ax = self.series.plot(title='Test', figsize=(16, 8))

        self.assertEqual(ax.title.get_text(), 'Test')
        assert_array_equal(np.round(ax.figure.get_size_inches()),
                           np.array((16., 8.)))

    def test_ts_area_lim(self):
        ax = self.ts.plot(kind='area', stacked=False)
        xmin, xmax = ax.get_xlim()
        lines = ax.get_lines()
        self.assertEqual(xmin, lines[0].get_data(orig=False)[0][0])       
        self.assertEqual(xmax, lines[0].get_data(orig=False)[0][-1])

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
            ax = _check_plot_works(d.plot, kind='area')
            self.assert_numpy_array_equal(ax.lines[0].get_ydata(), expected)
            ax = _check_plot_works(d.plot, kind='area', stacked=False)
            self.assert_numpy_array_equal(ax.lines[0].get_ydata(), expected)

    @slow
    def test_bar_log(self):
        expected = np.array([1., 10., 100., 1000.])

        if not self.mpl_le_1_2_1:
            expected = np.hstack((.1, expected, 1e4))

        ax = Series([200, 500]).plot(log=True, kind='bar')
        assert_array_equal(ax.yaxis.get_ticklocs(), expected)

    @slow
    def test_bar_ignore_index(self):
        df = Series([1, 2, 3, 4], index=['a', 'b', 'c', 'd'])
        ax = df.plot(kind='bar', use_index=False)
        expected = ['0', '1', '2', '3']
        result = [x.get_text() for x in ax.get_xticklabels()]
        self.assertEqual(result, expected)

    def test_rotation(self):
        df = DataFrame(randn(5, 5))
        ax = df.plot(rot=30)
        for l in ax.get_xticklabels():
            self.assertEqual(l.get_rotation(), 30)

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
        ax = _check_plot_works(series.plot, kind='pie')
        for t, expected in zip(ax.texts, series.index):
            self.assertEqual(t.get_text(), expected)
        self.assertEqual(ax.get_ylabel(), 'YLABEL')

        # without wedge labels
        ax = _check_plot_works(series.plot, kind='pie', labels=None)
        for t, expected in zip(ax.texts, [''] * 5):
            self.assertEqual(t.get_text(), expected)

        # with less colors than elements
        color_args = ['r', 'g', 'b']
        ax = _check_plot_works(series.plot, kind='pie', colors=color_args)

        import matplotlib.colors as colors
        conv = colors.colorConverter
        color_expected = ['r', 'g', 'b', 'r', 'g']
        for p, expected in zip(ax.patches, color_expected):
            self.assertEqual(p.get_facecolor(), conv.to_rgba(expected))

        # with labels and colors
        labels = ['A', 'B', 'C', 'D', 'E']
        color_args = ['r', 'g', 'b', 'c', 'm']
        ax = _check_plot_works(series.plot, kind='pie', labels=labels, colors=color_args)

        for t, expected in zip(ax.texts, labels):
            self.assertEqual(t.get_text(), expected)
        for p, expected in zip(ax.patches, color_args):
            self.assertEqual(p.get_facecolor(), conv.to_rgba(expected))

        # with autopct and fontsize
        ax = _check_plot_works(series.plot, kind='pie', colors=color_args,
                               autopct='%.2f', fontsize=7)
        pcts = ['{0:.2f}'.format(s * 100) for s in series.values / float(series.sum())]
        iters = [iter(series.index), iter(pcts)]
        expected_texts = list(it.next() for it in itertools.cycle(iters))
        for t, expected in zip(ax.texts, expected_texts):
            self.assertEqual(t.get_text(), expected)
            self.assertEqual(t.get_fontsize(), 7)

        # includes negative value
        with tm.assertRaises(ValueError):
            series = Series([1, 2, 0, 4, -1], index=['a', 'b', 'c', 'd', 'e'])
            series.plot(kind='pie')

        # includes nan
        series = Series([1, 2, np.nan, 4],
                        index=['a', 'b', 'c', 'd'], name='YLABEL')
        ax = _check_plot_works(series.plot, kind='pie')
        for t, expected in zip(ax.texts, series.index):
            self.assertEqual(t.get_text(), expected)

    @slow
    def test_hist(self):
        _check_plot_works(self.ts.hist)
        _check_plot_works(self.ts.hist, grid=False)
        _check_plot_works(self.ts.hist, figsize=(8, 10))
        _check_plot_works(self.ts.hist, by=self.ts.index.month)
        _check_plot_works(self.ts.hist, by=self.ts.index.month, bins=5)

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 1)
        _check_plot_works(self.ts.hist, ax=ax)
        _check_plot_works(self.ts.hist, ax=ax, figure=fig)
        _check_plot_works(self.ts.hist, figure=fig)
        tm.close()

        fig, (ax1, ax2) = plt.subplots(1, 2)
        _check_plot_works(self.ts.hist, figure=fig, ax=ax1)
        _check_plot_works(self.ts.hist, figure=fig, ax=ax2)

        with tm.assertRaises(ValueError):
            self.ts.hist(by=self.ts.index, figure=fig)

    @slow
    def test_hist_bins(self):
        df = DataFrame(np.random.randn(10, 2))
        ax = df.hist(bins=2)[0][0]
        self.assertEqual(len(ax.patches), 2)

    @slow
    def test_hist_layout(self):
        n = 10
        gender = tm.choice(['Male', 'Female'], size=n)
        df = DataFrame({'gender': gender,
                        'height': random.normal(66, 4, size=n), 'weight':
                        random.normal(161, 32, size=n)})
        with tm.assertRaises(ValueError):
            df.height.hist(layout=(1, 1))

        with tm.assertRaises(ValueError):
            df.height.hist(layout=[1, 1])

    @slow
    def test_hist_layout_with_by(self):
        import matplotlib.pyplot as plt
        n = 10
        gender = tm.choice(['Male', 'Female'], size=n)
        df = DataFrame({'gender': gender,
                        'height': random.normal(66, 4, size=n), 'weight':
                        random.normal(161, 32, size=n),
                        'category': random.randint(4, size=n)})
        _check_plot_works(df.height.hist, by=df.gender, layout=(2, 1))
        tm.close()

        _check_plot_works(df.height.hist, by=df.gender, layout=(1, 2))
        tm.close()

        _check_plot_works(df.weight.hist, by=df.category, layout=(1, 4))
        tm.close()

        _check_plot_works(df.weight.hist, by=df.category, layout=(4, 1))
        tm.close()

    @slow
    def test_hist_no_overlap(self):
        from matplotlib.pyplot import subplot, gcf, close
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
    def test_plot_fails_with_dupe_color_and_style(self):
        x = Series(randn(2))
        with tm.assertRaises(ValueError):
            x.plot(style='k--', color='k')

    @slow
    def test_hist_by_no_extra_plots(self):
        import matplotlib.pyplot as plt
        n = 10
        df = DataFrame({'gender': ['Male'] * 5 + ['Female'] * 5,
                        'height': random.normal(66, 4, size=n)})
        axes = df.height.hist(by=df.gender)
        self.assertEqual(len(plt.get_fignums()), 1)

    def test_plot_fails_when_ax_differs_from_figure(self):
        from pylab import figure, close
        fig1 = figure()
        fig2 = figure()
        ax1 = fig1.add_subplot(111)
        with tm.assertRaises(AssertionError):
            self.ts.hist(ax=ax1, figure=fig2)

    @slow
    def test_kde(self):
        _skip_if_no_scipy()
        _check_plot_works(self.ts.plot, kind='kde')
        _check_plot_works(self.ts.plot, kind='density')
        ax = self.ts.plot(kind='kde', logy=True)
        self.assertEqual(ax.get_yscale(), 'log')

    @slow
    def test_kde_kwargs(self):
        _skip_if_no_scipy()
        from numpy import linspace
        _check_plot_works(self.ts.plot, kind='kde', bw_method=.5, ind=linspace(-100,100,20))
        _check_plot_works(self.ts.plot, kind='density', bw_method=.5, ind=linspace(-100,100,20))
        ax = self.ts.plot(kind='kde', logy=True, bw_method=.5, ind=linspace(-100,100,20))
        self.assertEqual(ax.get_yscale(), 'log')

    @slow
    def test_kde_color(self):
        _skip_if_no_scipy()
        ax = self.ts.plot(kind='kde', logy=True, color='r')
        lines = ax.get_lines()
        self.assertEqual(len(lines), 1)
        self.assertEqual(lines[0].get_color(), 'r')

    @slow
    def test_autocorrelation_plot(self):
        from pandas.tools.plotting import autocorrelation_plot
        _check_plot_works(autocorrelation_plot, self.ts)
        _check_plot_works(autocorrelation_plot, self.ts.values)

        ax = autocorrelation_plot(self.ts, label='Test')
        t = ax.get_legend().get_texts()[0].get_text()
        self.assertEqual(t, 'Test')

    @slow
    def test_lag_plot(self):
        from pandas.tools.plotting import lag_plot
        _check_plot_works(lag_plot, self.ts)
        _check_plot_works(lag_plot, self.ts, lag=5)

    @slow
    def test_bootstrap_plot(self):
        from pandas.tools.plotting import bootstrap_plot
        _check_plot_works(bootstrap_plot, self.ts, size=10)

    def test_invalid_plot_data(self):
        s = Series(list('abcd'))
        for kind in plotting._common_kinds:
            with tm.assertRaises(TypeError):
                s.plot(kind=kind)

    @slow
    def test_valid_object_plot(self):
        s = Series(lrange(10), dtype=object)
        for kind in plotting._common_kinds:
            _check_plot_works(s.plot, kind=kind)

    def test_partially_invalid_plot_data(self):
        s = Series(['a', 'b', 1.0, 2])
        for kind in plotting._common_kinds:
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
            _check_has_errorbars(self, ax, xerr=0, yerr=1)
            ax = _check_plot_works(s.plot, yerr=s_err, kind=kind)
            _check_has_errorbars(self, ax, xerr=0, yerr=1)
            ax = _check_plot_works(s.plot, yerr=s_err.tolist(), kind=kind)
            _check_has_errorbars(self, ax, xerr=0, yerr=1)
            ax = _check_plot_works(s.plot, yerr=d_err, kind=kind)
            _check_has_errorbars(self, ax, xerr=0, yerr=1)
            ax = _check_plot_works(s.plot, xerr=0.2, yerr=0.2, kind=kind)
            _check_has_errorbars(self, ax, xerr=1, yerr=1)

        ax = _check_plot_works(s.plot, xerr=s_err)
        _check_has_errorbars(self, ax, xerr=1, yerr=0)

        # test time series plotting
        ix = date_range('1/1/2000', '1/1/2001', freq='M')
        ts = Series(np.arange(12), index=ix, name='x')
        ts_err = Series(np.random.randn(12), index=ix)
        td_err =  DataFrame(randn(12, 2), index=ix, columns=['x', 'y'])

        ax = _check_plot_works(ts.plot, yerr=ts_err)
        _check_has_errorbars(self, ax, xerr=0, yerr=1)
        ax = _check_plot_works(ts.plot, yerr=td_err)
        _check_has_errorbars(self, ax, xerr=0, yerr=1)

        # check incorrect lengths and types
        with tm.assertRaises(ValueError):
            s.plot(yerr=np.arange(11))

        s_err = ['zzz']*10
        with tm.assertRaises(TypeError):
            s.plot(yerr=s_err)

    def test_table(self):
        _check_plot_works(self.series.plot, table=True)
        _check_plot_works(self.series.plot, table=self.series)


@tm.mplskip
class TestDataFramePlots(tm.TestCase):
    def setUp(self):
        import matplotlib as mpl
        self.mpl_le_1_2_1 = str(mpl.__version__) <= LooseVersion('1.2.1')

    def tearDown(self):
        tm.close()

    @slow
    def test_plot(self):
        df = tm.makeTimeDataFrame()
        _check_plot_works(df.plot, grid=False)
        _check_plot_works(df.plot, subplots=True)
        _check_plot_works(df.plot, subplots=True, use_index=False)

        df = DataFrame({'x': [1, 2], 'y': [3, 4]})
        self._check_plot_fails(df.plot, kind='line', blarg=True)

        df = DataFrame(np.random.rand(10, 3),
                       index=list(string.ascii_letters[:10]))
        _check_plot_works(df.plot, use_index=True)
        _check_plot_works(df.plot, sort_columns=False)
        _check_plot_works(df.plot, yticks=[1, 5, 10])
        _check_plot_works(df.plot, xticks=[1, 5, 10])
        _check_plot_works(df.plot, ylim=(-100, 100), xlim=(-100, 100))
        _check_plot_works(df.plot, subplots=True, title='blah')
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
        _check_plot_works(df.plot, kind='bar', subplots=True)

    def test_nonnumeric_exclude(self):
        import matplotlib.pyplot as plt
        df = DataFrame({'A': ["x", "y", "z"], 'B': [1, 2, 3]})
        ax = df.plot()
        self.assertEqual(len(ax.get_lines()), 1)  # B was plotted

    @slow
    def test_implicit_label(self):
        df = DataFrame(randn(10, 3), columns=['a', 'b', 'c'])
        ax = df.plot(x='a', y='b')
        self.assertEqual(ax.xaxis.get_label().get_text(), 'a')

    @slow
    def test_explicit_label(self):
        df = DataFrame(randn(10, 3), columns=['a', 'b', 'c'])
        ax = df.plot(x='a', y='b', label='LABEL')
        self.assertEqual(ax.xaxis.get_label().get_text(), 'LABEL')

    @slow
    def test_plot_xy(self):
        import matplotlib.pyplot as plt
        # columns.inferred_type == 'string'
        df = tm.makeTimeDataFrame()
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

        self.assertEqual(ax.title.get_text(), 'Test')
        assert_array_equal(np.round(ax.figure.get_size_inches()),
                           np.array((16., 8.)))

        # columns.inferred_type == 'mixed'
        # TODO add MultiIndex test

    @slow
    def test_logscales(self):
        df = DataFrame({'a': np.arange(100)},
                       index=np.arange(100))
        ax = df.plot(logy=True)
        self.assertEqual(ax.xaxis.get_scale(), 'linear')
        self.assertEqual(ax.yaxis.get_scale(), 'log')

        ax = df.plot(logx=True)
        self.assertEqual(ax.xaxis.get_scale(), 'log')
        self.assertEqual(ax.yaxis.get_scale(), 'linear')

        ax = df.plot(loglog=True)
        self.assertEqual(ax.xaxis.get_scale(), 'log')
        self.assertEqual(ax.yaxis.get_scale(), 'log')

    @slow
    def test_xcompat(self):
        import pandas as pd
        import matplotlib.pyplot as plt

        df = tm.makeTimeDataFrame()
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
        tm.assert_isinstance(lines[0].get_xdata(), PeriodIndex)

        tm.close()
        # useful if you're plotting a bunch together
        with pd.plot_params.use('x_compat', True):
            ax = df.plot()
            lines = ax.get_lines()
            self.assertNotIsInstance(lines[0].get_xdata(), PeriodIndex)

        tm.close()
        ax = df.plot()
        lines = ax.get_lines()
        tm.assert_isinstance(lines[0].get_xdata(), PeriodIndex)

    def test_unsorted_index(self):
        df = DataFrame({'y': np.arange(100)},
                       index=np.arange(99, -1, -1), dtype=np.int64)
        ax = df.plot()
        l = ax.get_lines()[0]
        rs = l.get_xydata()
        rs = Series(rs[:, 1], rs[:, 0], dtype=np.int64)
        tm.assert_series_equal(rs, df.y)

    def _check_data(self, xp, rs):
        xp_lines = xp.get_lines()
        rs_lines = rs.get_lines()

        def check_line(xpl, rsl):
            xpdata = xpl.get_xydata()
            rsdata = rsl.get_xydata()
            assert_array_equal(xpdata, rsdata)

        [check_line(xpl, rsl) for xpl, rsl in zip(xp_lines, rs_lines)]
        tm.close()

    @slow
    def test_subplots(self):
        df = DataFrame(np.random.rand(10, 3),
                       index=list(string.ascii_letters[:10]))

        for kind in ['bar', 'barh', 'line', 'area']:
            axes = df.plot(kind=kind, subplots=True, sharex=True, legend=True)

            for ax, column in zip(axes, df.columns):
                self._check_legend_labels(ax, [column])

            axes = df.plot(kind=kind, subplots=True, sharex=True)
            for ax in axes[:-2]:
                [self.assert_(not label.get_visible())
                 for label in ax.get_xticklabels()]
                [self.assert_(label.get_visible())
                 for label in ax.get_yticklabels()]

            [self.assert_(label.get_visible())
             for label in axes[-1].get_xticklabels()]
            [self.assert_(label.get_visible())
             for label in axes[-1].get_yticklabels()]

            axes = df.plot(kind=kind, subplots=True, sharex=False)
            for ax in axes:
                [self.assert_(label.get_visible())
                 for label in ax.get_xticklabels()]
                [self.assert_(label.get_visible())
                 for label in ax.get_yticklabels()]

            axes = df.plot(kind=kind, subplots=True, legend=False)
            for ax in axes:
                self.assertTrue(ax.get_legend() is None)

    def test_negative_log(self):
        df = - DataFrame(rand(6, 4),
                       index=list(string.ascii_letters[:6]),
                       columns=['x', 'y', 'z', 'four'])

        with tm.assertRaises(ValueError):
            df.plot(kind='area', logy=True)
        with tm.assertRaises(ValueError):
            df.plot(kind='area', loglog=True)

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

            ax = _check_plot_works(d.plot, kind='area')
            self.assert_numpy_array_equal(ax.lines[0].get_ydata(), expected1)
            self.assert_numpy_array_equal(ax.lines[1].get_ydata(), expected1 + expected2)

            ax = _check_plot_works(d.plot, kind='area', stacked=False)
            self.assert_numpy_array_equal(ax.lines[0].get_ydata(), expected1)
            self.assert_numpy_array_equal(ax.lines[1].get_ydata(), expected2)

    def test_area_lim(self):
        df = DataFrame(rand(6, 4),
                       columns=['x', 'y', 'z', 'four'])

        neg_df = - df
        for stacked in [True, False]:
            ax = _check_plot_works(df.plot, kind='area', stacked=stacked)
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            lines = ax.get_lines()
            self.assertEqual(xmin, lines[0].get_data()[0][0])        
            self.assertEqual(xmax, lines[0].get_data()[0][-1]) 
            self.assertEqual(ymin, 0) 

            ax = _check_plot_works(neg_df.plot, kind='area', stacked=stacked)
            ymin, ymax = ax.get_ylim()
            self.assertEqual(ymax, 0) 

    @slow
    def test_bar_colors(self):
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors

        default_colors = plt.rcParams.get('axes.color_cycle')
        custom_colors = 'rgcby'

        df = DataFrame(randn(5, 5))
        ax = df.plot(kind='bar')

        rects = ax.patches

        conv = colors.colorConverter
        for i, rect in enumerate(rects[::5]):
            xp = conv.to_rgba(default_colors[i % len(default_colors)])
            rs = rect.get_facecolor()
            self.assertEqual(xp, rs)

        tm.close()

        ax = df.plot(kind='bar', color=custom_colors)

        rects = ax.patches

        conv = colors.colorConverter
        for i, rect in enumerate(rects[::5]):
            xp = conv.to_rgba(custom_colors[i])
            rs = rect.get_facecolor()
            self.assertEqual(xp, rs)

        tm.close()
        from matplotlib import cm

        # Test str -> colormap functionality
        ax = df.plot(kind='bar', colormap='jet')

        rects = ax.patches

        rgba_colors = lmap(cm.jet, np.linspace(0, 1, 5))
        for i, rect in enumerate(rects[::5]):
            xp = rgba_colors[i]
            rs = rect.get_facecolor()
            self.assertEqual(xp, rs)

        tm.close()

        # Test colormap functionality
        ax = df.plot(kind='bar', colormap=cm.jet)

        rects = ax.patches

        rgba_colors = lmap(cm.jet, np.linspace(0, 1, 5))
        for i, rect in enumerate(rects[::5]):
            xp = rgba_colors[i]
            rs = rect.get_facecolor()
            self.assertEqual(xp, rs)

        tm.close()
        df.ix[:, [0]].plot(kind='bar', color='DodgerBlue')

    @slow
    def test_bar_linewidth(self):
        df = DataFrame(randn(5, 5))

        # regular
        ax = df.plot(kind='bar', linewidth=2)
        for r in ax.patches:
            self.assertEqual(r.get_linewidth(), 2)

        # stacked
        ax = df.plot(kind='bar', stacked=True, linewidth=2)
        for r in ax.patches:
            self.assertEqual(r.get_linewidth(), 2)

        # subplots
        axes = df.plot(kind='bar', linewidth=2, subplots=True)
        for ax in axes:
            for r in ax.patches:
                self.assertEqual(r.get_linewidth(), 2)

    @slow
    def test_bar_barwidth(self):
        df = DataFrame(randn(5, 5))

        width = 0.9

        # regular
        ax = df.plot(kind='bar', width=width)
        for r in ax.patches:
            self.assertEqual(r.get_width(), width / len(df.columns))

        # stacked
        ax = df.plot(kind='bar', stacked=True, width=width)
        for r in ax.patches:
            self.assertEqual(r.get_width(), width)

        # horizontal regular
        ax = df.plot(kind='barh', width=width)
        for r in ax.patches:
            self.assertEqual(r.get_height(), width / len(df.columns))

        # horizontal stacked
        ax = df.plot(kind='barh', stacked=True, width=width)
        for r in ax.patches:
            self.assertEqual(r.get_height(), width)

        # subplots
        axes = df.plot(kind='bar', width=width, subplots=True)
        for ax in axes:
            for r in ax.patches:
                self.assertEqual(r.get_width(), width)

        # horizontal subplots
        axes = df.plot(kind='barh', width=width, subplots=True)
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
    def test_plot_scatter(self):
        from matplotlib.pylab import close
        df = DataFrame(randn(6, 4),
                       index=list(string.ascii_letters[:6]),
                       columns=['x', 'y', 'z', 'four'])

        _check_plot_works(df.plot, x='x', y='y', kind='scatter')
        _check_plot_works(df.plot, x=1, y=2, kind='scatter')

        with tm.assertRaises(ValueError):
            df.plot(x='x', kind='scatter')
        with tm.assertRaises(ValueError):
            df.plot(y='y', kind='scatter')

        # GH 6951
        axes = df.plot(x='x', y='y', kind='scatter', subplots=True)
        self.assertEqual(len(axes[0].figure.axes), 1)

    @slow
    def test_plot_bar(self):
        from matplotlib.pylab import close
        df = DataFrame(randn(6, 4),
                       index=list(string.ascii_letters[:6]),
                       columns=['one', 'two', 'three', 'four'])

        _check_plot_works(df.plot, kind='bar')
        _check_plot_works(df.plot, kind='bar', legend=False)
        _check_plot_works(df.plot, kind='bar', subplots=True)
        _check_plot_works(df.plot, kind='bar', stacked=True)

        df = DataFrame(randn(10, 15),
                       index=list(string.ascii_letters[:10]),
                       columns=lrange(15))
        _check_plot_works(df.plot, kind='bar')

        df = DataFrame({'a': [0, 1], 'b': [1, 0]})
        _check_plot_works(df.plot, kind='bar')

    def _check_bar_alignment(self, df, kind='bar', stacked=False,
                             subplots=False, align='center',
                             width=0.5, position=0.5):

        axes = df.plot(kind=kind, stacked=stacked, subplots=subplots,
                       align=align, width=width, position=position,
                       grid=True)

        tick_pos = np.arange(len(df))

        if not isinstance(axes, np.ndarray):
            axes = [axes]

        for ax in axes:
            if kind == 'bar':
                axis = ax.xaxis
                ax_min, ax_max = ax.get_xlim()
            elif kind == 'barh':
                axis = ax.yaxis
                ax_min, ax_max = ax.get_ylim()
            else:
                raise ValueError

            p = ax.patches[0]
            if kind == 'bar' and (stacked is True or subplots is True):
                edge = p.get_x()
                center = edge + p.get_width() * position
                tickoffset = width * position
            elif kind == 'bar' and stacked is False:
                center = p.get_x() + p.get_width() * len(df.columns) * position
                edge = p.get_x()
                if align == 'edge':
                    tickoffset = width * (position - 0.5) + p.get_width() * 1.5
                else:
                    tickoffset = width * position + p.get_width()
            elif kind == 'barh' and (stacked is True or subplots is True):
                center = p.get_y() + p.get_height() * position
                edge = p.get_y()
                tickoffset = width * position
            elif kind == 'barh' and stacked is False:
                center = p.get_y() + p.get_height() * len(df.columns) * position
                edge = p.get_y()
                if align == 'edge':
                    tickoffset = width * (position - 0.5) + p.get_height() * 1.5
                else:
                    tickoffset = width * position + p.get_height()
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

            # Check starting point and axes limit margin
            self.assertEqual(ax_min, tick_pos[0] - tickoffset - 0.25)
            self.assertEqual(ax_max, tick_pos[-1] - tickoffset + 1)
            # Check tick locations and axes limit margin
            t_min = axis.get_ticklocs()[0] - tickoffset
            t_max = axis.get_ticklocs()[-1] - tickoffset
            self.assertAlmostEqual(ax_min, t_min - 0.25)
            self.assertAlmostEqual(ax_max, t_max + 1.0)
        return axes

    @slow
    def test_bar_stacked_center(self):
        # GH2157
        df = DataFrame({'A': [3] * 5, 'B': lrange(5)}, index=lrange(5))
        axes = self._check_bar_alignment(df, kind='bar', stacked=True)
        # Check the axes has the same drawing range before fixing # GH4525
        self.assertEqual(axes[0].get_xlim(), (-0.5, 4.75))

        self._check_bar_alignment(df, kind='bar', stacked=True, width=0.9)

        axes = self._check_bar_alignment(df, kind='barh', stacked=True)
        self.assertEqual(axes[0].get_ylim(), (-0.5, 4.75))

        self._check_bar_alignment(df, kind='barh', stacked=True, width=0.9)

    @slow
    def test_bar_center(self):
        df = DataFrame({'A': [3] * 5, 'B': lrange(5)}, index=lrange(5))
        axes = self._check_bar_alignment(df, kind='bar', stacked=False)
        self.assertEqual(axes[0].get_xlim(), (-0.75, 4.5))

        self._check_bar_alignment(df, kind='bar', stacked=False, width=0.9)

        axes = self._check_bar_alignment(df, kind='barh', stacked=False)
        self.assertEqual(axes[0].get_ylim(), (-0.75, 4.5))

        self._check_bar_alignment(df, kind='barh', stacked=False, width=0.9)

    @slow
    def test_bar_subplots_center(self):
        df = DataFrame({'A': [3] * 5, 'B': lrange(5)}, index=lrange(5))
        axes = self._check_bar_alignment(df, kind='bar', subplots=True)
        for ax in axes:
            self.assertEqual(ax.get_xlim(), (-0.5, 4.75))

        self._check_bar_alignment(df, kind='bar', subplots=True, width=0.9)

        axes = self._check_bar_alignment(df, kind='barh', subplots=True)
        for ax in axes:
            self.assertEqual(ax.get_ylim(), (-0.5, 4.75))

        self._check_bar_alignment(df, kind='barh', subplots=True, width=0.9)

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
        ax = df.plot(kind='bar', grid=True, log=True)
        assert_array_equal(ax.yaxis.get_ticklocs(), expected)

    @slow
    def test_bar_log_subplots(self):
        expected = np.array([1., 10., 100., 1000.])
        if not self.mpl_le_1_2_1:
            expected = np.hstack((.1, expected, 1e4))

        ax = DataFrame([Series([200, 300]),
                        Series([300, 500])]).plot(log=True, kind='bar',
                                                  subplots=True)

        assert_array_equal(ax[0].yaxis.get_ticklocs(), expected)
        assert_array_equal(ax[1].yaxis.get_ticklocs(), expected)

    @slow
    def test_boxplot(self):
        df = DataFrame(randn(6, 4),
                       index=list(string.ascii_letters[:6]),
                       columns=['one', 'two', 'three', 'four'])
        df['indic'] = ['foo', 'bar'] * 3
        df['indic2'] = ['foo', 'bar', 'foo'] * 2

        _check_plot_works(df.boxplot)
        _check_plot_works(df.boxplot, column=['one', 'two'])
        _check_plot_works(df.boxplot, column=['one', 'two'], by='indic')
        _check_plot_works(df.boxplot, column='one', by=['indic', 'indic2'])
        _check_plot_works(df.boxplot, by='indic')
        _check_plot_works(df.boxplot, by=['indic', 'indic2'])
        _check_plot_works(plotting.boxplot, df['one'])
        _check_plot_works(df.boxplot, notch=1)
        _check_plot_works(df.boxplot, by='indic', notch=1)

        df = DataFrame(np.random.rand(10, 2), columns=['Col1', 'Col2'])
        df['X'] = Series(['A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B'])
        _check_plot_works(df.boxplot, by='X')

        # When ax is supplied, existing axes should be used:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        axes = df.boxplot('Col1', by='X', ax=ax)
        self.assertIs(ax.get_axes(), axes)

        # Multiple columns with an ax argument is not supported
        fig, ax = plt.subplots()
        self.assertRaisesRegexp(
            ValueError, 'existing axis', df.boxplot,
            column=['Col1', 'Col2'], by='X', ax=ax
        )

        # When by is None, check that all relevant lines are present in the dict
        fig, ax = plt.subplots()
        d = df.boxplot(ax=ax)
        lines = list(itertools.chain.from_iterable(d.values()))
        self.assertEqual(len(ax.get_lines()), len(lines))

    @slow
    def test_kde(self):
        _skip_if_no_scipy()
        df = DataFrame(randn(100, 4))
        _check_plot_works(df.plot, kind='kde')
        _check_plot_works(df.plot, kind='kde', subplots=True)
        ax = df.plot(kind='kde')
        self._check_legend_labels(ax, df.columns)
        axes = df.plot(kind='kde', logy=True, subplots=True)
        for ax in axes:
            self.assertEqual(ax.get_yscale(), 'log')

    @slow
    def test_hist(self):
        import matplotlib.pyplot as plt
        df = DataFrame(randn(100, 4))
        _check_plot_works(df.hist)
        _check_plot_works(df.hist, grid=False)

        # make sure layout is handled
        df = DataFrame(randn(100, 3))
        _check_plot_works(df.hist)
        axes = df.hist(grid=False)
        self.assert_(not axes[1, 1].get_visible())

        df = DataFrame(randn(100, 1))
        _check_plot_works(df.hist)

        # make sure layout is handled
        df = DataFrame(randn(100, 6))
        _check_plot_works(df.hist)

        # make sure sharex, sharey is handled
        _check_plot_works(df.hist, sharex=True, sharey=True)

        # handle figsize arg
        _check_plot_works(df.hist, figsize=(8, 10))

        # check bins argument
        _check_plot_works(df.hist, bins=5)

        # make sure xlabelsize and xrot are handled
        ser = df[0]
        xf, yf = 20, 20
        xrot, yrot = 30, 30
        ax = ser.hist(xlabelsize=xf, xrot=30, ylabelsize=yf, yrot=30)
        ytick = ax.get_yticklabels()[0]
        xtick = ax.get_xticklabels()[0]
        self.assertAlmostEqual(ytick.get_fontsize(), yf)
        self.assertAlmostEqual(ytick.get_rotation(), yrot)
        self.assertAlmostEqual(xtick.get_fontsize(), xf)
        self.assertAlmostEqual(xtick.get_rotation(), xrot)

        xf, yf = 20, 20
        xrot, yrot = 30, 30
        axes = df.hist(xlabelsize=xf, xrot=30, ylabelsize=yf, yrot=30)
        for i, ax in enumerate(axes.ravel()):
            if i < len(df.columns):
                ytick = ax.get_yticklabels()[0]
                xtick = ax.get_xticklabels()[0]
                self.assertAlmostEqual(ytick.get_fontsize(), yf)
                self.assertAlmostEqual(ytick.get_rotation(), yrot)
                self.assertAlmostEqual(xtick.get_fontsize(), xf)
                self.assertAlmostEqual(xtick.get_rotation(), xrot)

        tm.close()
        # make sure kwargs to hist are handled
        ax = ser.hist(normed=True, cumulative=True, bins=4)
        # height of last bin (index 5) must be 1.0
        self.assertAlmostEqual(ax.get_children()[5].get_height(), 1.0)

        tm.close()
        ax = ser.hist(log=True)
        # scale of y must be 'log'
        self.assertEqual(ax.get_yscale(), 'log')

        tm.close()

        # propagate attr exception from matplotlib.Axes.hist
        with tm.assertRaises(AttributeError):
            ser.hist(foo='bar')

    @slow
    def test_hist_layout(self):
        import matplotlib.pyplot as plt
        df = DataFrame(randn(100, 4))

        layout_to_expected_size = (
            {'layout': None, 'expected_size': (2, 2)},  # default is 2x2
            {'layout': (2, 2), 'expected_size': (2, 2)},
            {'layout': (4, 1), 'expected_size': (4, 1)},
            {'layout': (1, 4), 'expected_size': (1, 4)},
            {'layout': (3, 3), 'expected_size': (3, 3)},
        )

        for layout_test in layout_to_expected_size:
            ax = df.hist(layout=layout_test['layout'])
            self.assertEqual(len(ax), layout_test['expected_size'][0])
            self.assertEqual(len(ax[0]), layout_test['expected_size'][1])

        # layout too small for all 4 plots
        with tm.assertRaises(ValueError):
            df.hist(layout=(1, 1))

        # invalid format for layout
        with tm.assertRaises(ValueError):
            df.hist(layout=(1,))

    @slow
    def test_scatter(self):
        _skip_if_no_scipy()

        df = DataFrame(randn(100, 2))
        import pandas.tools.plotting as plt

        def scat(**kwds):
            return plt.scatter_matrix(df, **kwds)

        _check_plot_works(scat)
        _check_plot_works(scat, marker='+')
        _check_plot_works(scat, vmin=0)
        _check_plot_works(scat, diagonal='kde')
        _check_plot_works(scat, diagonal='density')
        _check_plot_works(scat, diagonal='hist')
        _check_plot_works(scat, range_padding=.1)

        def scat2(x, y, by=None, ax=None, figsize=None):
            return plt.scatter_plot(df, x, y, by, ax, figsize=None)

        _check_plot_works(scat2, 0, 1)
        grouper = Series(np.repeat([1, 2, 3, 4, 5], 20), df.index)
        _check_plot_works(scat2, 0, 1, by=grouper)

    @slow
    def test_andrews_curves(self):
        from pandas import read_csv
        from pandas.tools.plotting import andrews_curves
        from matplotlib import cm
        
        path = os.path.join(curpath(), 'data', 'iris.csv')
        df = read_csv(path)

        _check_plot_works(andrews_curves, df, 'Name')
        _check_plot_works(andrews_curves, df, 'Name',
                          color=('#556270', '#4ECDC4', '#C7F464'))
        _check_plot_works(andrews_curves, df, 'Name',
                          color=['dodgerblue', 'aquamarine', 'seagreen'])
        _check_plot_works(andrews_curves, df, 'Name', colormap=cm.jet)

        colors = ['b', 'g', 'r']
        df = DataFrame({"A": [1, 2, 3],
                        "B": [1, 2, 3],
                        "C": [1, 2, 3],
                        "Name": colors})
        ax = andrews_curves(df, 'Name', color=colors)
        legend_colors = [l.get_color() for l in ax.legend().get_lines()]
        self.assertEqual(colors, legend_colors)
        
        with tm.assert_produces_warning(FutureWarning):
            andrews_curves(data=df, class_column='Name')

    @slow
    def test_parallel_coordinates(self):
        from pandas import read_csv
        from pandas.tools.plotting import parallel_coordinates
        from matplotlib import cm
        path = os.path.join(curpath(), 'data', 'iris.csv')
        df = read_csv(path)
        _check_plot_works(parallel_coordinates, df, 'Name')
        _check_plot_works(parallel_coordinates, df, 'Name',
                          color=('#556270', '#4ECDC4', '#C7F464'))
        _check_plot_works(parallel_coordinates, df, 'Name',
                          color=['dodgerblue', 'aquamarine', 'seagreen'])
        _check_plot_works(parallel_coordinates, df, 'Name', colormap=cm.jet)

        df = read_csv(path, header=None, skiprows=1, names=[1, 2, 4, 8,
                                                            'Name'])
        _check_plot_works(parallel_coordinates, df, 'Name', use_columns=True)
        _check_plot_works(parallel_coordinates, df, 'Name',
                          xticks=[1, 5, 25, 125])
        
        colors = ['b', 'g', 'r']
        df = DataFrame({"A": [1, 2, 3],
                        "B": [1, 2, 3],
                        "C": [1, 2, 3],
                        "Name": colors})
        ax = parallel_coordinates(df, 'Name', color=colors)
        legend_colors = [l.get_color() for l in ax.legend().get_lines()]
        self.assertEqual(colors, legend_colors)
        
        with tm.assert_produces_warning(FutureWarning):
            parallel_coordinates(df, 'Name', colors=colors)
        
        with tm.assert_produces_warning(FutureWarning):
            parallel_coordinates(data=df, class_column='Name')

    @slow
    def test_radviz(self):
        from pandas import read_csv
        from pandas.tools.plotting import radviz
        from matplotlib import cm

        path = os.path.join(curpath(), 'data', 'iris.csv')
        df = read_csv(path)
        _check_plot_works(radviz, df, 'Name')
        _check_plot_works(radviz, df, 'Name',
                          color=('#556270', '#4ECDC4', '#C7F464'))
        _check_plot_works(radviz, df, 'Name',
                          color=['dodgerblue', 'aquamarine', 'seagreen'])
        _check_plot_works(radviz, df, 'Name', colormap=cm.jet)

        colors = [[0., 0., 1., 1.],
                  [0., 0.5, 1., 1.],
                  [1., 0., 0., 1.]]
        df = DataFrame({"A": [1, 2, 3],
                        "B": [2, 1, 3],
                        "C": [3, 2, 1],
                        "Name": ['b', 'g', 'r']})
        ax = radviz(df, 'Name', color=colors)
        legend_colors = [c.get_facecolor().squeeze().tolist()
                         for c in ax.collections]
        self.assertEqual(colors, legend_colors)

    @slow
    def test_plot_int_columns(self):
        df = DataFrame(randn(100, 4)).cumsum()
        _check_plot_works(df.plot, legend=True)

    def _check_legend_labels(self, ax, labels):
        import pandas.core.common as com
        labels = [com.pprint_thing(l) for l in labels]
        self.assertTrue(ax.get_legend() is not None)
        legend_labels = [t.get_text() for t in ax.get_legend().get_texts()]
        self.assertEqual(labels, legend_labels)

    @slow
    def test_df_legend_labels(self):
        kinds = 'line', 'bar', 'barh', 'kde', 'density', 'area'
        df = DataFrame(rand(3, 3), columns=['a', 'b', 'c'])
        df2 = DataFrame(rand(3, 3), columns=['d', 'e', 'f'])
        df3 = DataFrame(rand(3, 3), columns=['g', 'h', 'i'])
        df4 = DataFrame(rand(3, 3), columns=['j', 'k', 'l'])

        for kind in kinds:
            ax = df.plot(kind=kind, legend=True)
            self._check_legend_labels(ax, df.columns)

            ax = df2.plot(kind=kind, legend=False, ax=ax)
            self._check_legend_labels(ax, df.columns)

            ax = df3.plot(kind=kind, legend=True, ax=ax)
            self._check_legend_labels(ax, df.columns + df3.columns)

            ax = df4.plot(kind=kind, legend='reverse', ax=ax)
            expected = list(df.columns + df3.columns) + list(reversed(df4.columns))
            self._check_legend_labels(ax, expected)

        # Secondary Y
        ax = df.plot(legend=True, secondary_y='b')
        self._check_legend_labels(ax, ['a', 'b (right)', 'c'])
        ax = df2.plot(legend=False, ax=ax)
        self._check_legend_labels(ax, ['a', 'b (right)', 'c'])
        ax = df3.plot(kind='bar', legend=True, secondary_y='h', ax=ax)
        self._check_legend_labels(ax, ['a', 'b (right)', 'c', 'g', 'h (right)', 'i'])

        # Time Series
        ind = date_range('1/1/2014', periods=3)
        df = DataFrame(randn(3, 3), columns=['a', 'b', 'c'], index=ind)
        df2 = DataFrame(randn(3, 3), columns=['d', 'e', 'f'], index=ind)
        df3 = DataFrame(randn(3, 3), columns=['g', 'h', 'i'], index=ind)
        ax = df.plot(legend=True, secondary_y='b')
        self._check_legend_labels(ax, ['a', 'b (right)', 'c'])
        ax = df2.plot(legend=False, ax=ax)
        self._check_legend_labels(ax, ['a', 'b (right)', 'c'])
        ax = df3.plot(legend=True, ax=ax)
        self._check_legend_labels(ax, ['a', 'b (right)', 'c', 'g', 'h', 'i'])

        # scatter
        ax = df.plot(kind='scatter', x='a', y='b', label='data1')
        self._check_legend_labels(ax, ['data1'])
        ax = df2.plot(kind='scatter', x='d', y='e', legend=False,
                        label='data2', ax=ax)
        self._check_legend_labels(ax, ['data1'])
        ax = df3.plot(kind='scatter', x='g', y='h', label='data3', ax=ax)
        self._check_legend_labels(ax, ['data1', 'data3'])

    def test_legend_name(self):
        multi = DataFrame(randn(4, 4),
                          columns=[np.array(['a', 'a', 'b', 'b']),
                                   np.array(['x', 'y', 'x', 'y'])])
        multi.columns.names = ['group', 'individual']

        ax = multi.plot()
        leg_title = ax.legend_.get_title()
        self.assertEqual(leg_title.get_text(), 'group,individual')

        df = DataFrame(randn(5, 5))
        ax = df.plot(legend=True, ax=ax)
        leg_title = ax.legend_.get_title()
        self.assertEqual(leg_title.get_text(), 'group,individual')

        df.columns.name = 'new'
        ax = df.plot(legend=False, ax=ax)
        leg_title = ax.legend_.get_title()
        self.assertEqual(leg_title.get_text(), 'group,individual')

        ax = df.plot(legend=True, ax=ax)
        leg_title = ax.legend_.get_title()
        self.assertEqual(leg_title.get_text(), 'new')

    def _check_plot_fails(self, f, *args, **kwargs):
        with tm.assertRaises(Exception):
            f(*args, **kwargs)

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

    def check_line_colors(self, colors, lines):
        for i, l in enumerate(lines):
            xp = colors[i]
            rs = l.get_color()
            self.assertEqual(xp, rs)
            
    def check_collection_colors(self, colors, cols):
        from matplotlib.colors import ColorConverter
        conv = ColorConverter()
        for i, c in enumerate(cols):
            xp = colors[i]
            xp = conv.to_rgba(xp)
            rs = c.get_facecolor()[0]
            for x, y in zip(xp, rs):
                self.assertEqual(x, y)

    @slow
    def test_line_colors(self):
        import matplotlib.pyplot as plt
        import sys
        from matplotlib import cm

        custom_colors = 'rgcby'
        df = DataFrame(randn(5, 5))

        ax = df.plot(color=custom_colors)
        self.check_line_colors(custom_colors, ax.get_lines())

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
        self.check_line_colors(rgba_colors, ax.get_lines())

        tm.close()

        ax = df.plot(colormap=cm.jet)

        rgba_colors = lmap(cm.jet, np.linspace(0, 1, len(df)))
        self.check_line_colors(rgba_colors, ax.get_lines())

        # make color a list if plotting one column frame
        # handles cases like df.plot(color='DodgerBlue')
        tm.close()
        df.ix[:, [0]].plot(color='DodgerBlue')

    @slow
    def test_area_colors(self):
        from matplotlib import cm
        from matplotlib.collections import PolyCollection

        custom_colors = 'rgcby'
        df = DataFrame(rand(5, 5))

        ax = df.plot(kind='area', color=custom_colors)
        self.check_line_colors(custom_colors, ax.get_lines())
        poly = [o for o in ax.get_children() if isinstance(o, PolyCollection)]
        self.check_collection_colors(custom_colors, poly)

        ax = df.plot(kind='area', colormap='jet')
        rgba_colors = lmap(cm.jet, np.linspace(0, 1, len(df)))
        self.check_line_colors(rgba_colors, ax.get_lines())
        poly = [o for o in ax.get_children() if isinstance(o, PolyCollection)]
        self.check_collection_colors(rgba_colors, poly)

        ax = df.plot(kind='area', colormap=cm.jet)
        rgba_colors = lmap(cm.jet, np.linspace(0, 1, len(df)))
        self.check_line_colors(rgba_colors, ax.get_lines())
        poly = [o for o in ax.get_children() if isinstance(o, PolyCollection)]
        self.check_collection_colors(rgba_colors, poly)

    def test_default_color_cycle(self):
        import matplotlib.pyplot as plt
        plt.rcParams['axes.color_cycle'] = list('rgbk')

        df = DataFrame(randn(5, 3))
        ax = df.plot()

        lines = ax.get_lines()
        for i, l in enumerate(lines):
            xp = plt.rcParams['axes.color_cycle'][i]
            rs = l.get_color()
            self.assertEqual(xp, rs)

    def test_unordered_ts(self):
        df = DataFrame(np.array([3.0, 2.0, 1.0]),
                       index=[date(2012, 10, 1),
                              date(2012, 9, 1),
                              date(2012, 8, 1)],
                       columns=['test'])
        ax = df.plot()
        xticks = ax.lines[0].get_xdata()
        self.assert_(xticks[0] < xticks[1])
        ydata = ax.lines[0].get_ydata()
        assert_array_equal(ydata, np.array([1.0, 2.0, 3.0]))

    def test_all_invalid_plot_data(self):
        df = DataFrame(list('abcd'))
        for kind in plotting._common_kinds:
            with tm.assertRaises(TypeError):
                df.plot(kind=kind)

    @slow
    def test_partially_invalid_plot_data(self):
        with tm.RNGContext(42):
            df = DataFrame(randn(10, 2), dtype=object)
            df[np.random.rand(df.shape[0]) > 0.5] = 'a'
            for kind in plotting._common_kinds:
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
        df = DataFrame({"A": np.random.uniform(size=20),
                        "B": np.random.uniform(size=20),
                        "C": np.arange(20) + np.random.uniform(size=20)})

        ax = df.plot(kind='hexbin', x='A', y='B', gridsize=10)
        # TODO: need better way to test. This just does existence.
        self.assertEqual(len(ax.collections), 1)

        # GH 6951
        axes = df.plot(x='A', y='B', kind='hexbin', subplots=True)
        # hexbin should have 2 axes, 1 for plotting and another is colorbar
        self.assertEqual(len(axes[0].figure.axes), 2)

    @slow
    def test_hexbin_with_c(self):
        df = DataFrame({"A": np.random.uniform(size=20),
                        "B": np.random.uniform(size=20),
                        "C": np.arange(20) + np.random.uniform(size=20)})

        ax = df.plot(kind='hexbin', x='A', y='B', C='C')
        self.assertEqual(len(ax.collections), 1)

        ax = df.plot(kind='hexbin', x='A', y='B', C='C',
                          reduce_C_function=np.std)
        self.assertEqual(len(ax.collections), 1)

    @slow
    def test_hexbin_cmap(self):
        df = DataFrame({"A": np.random.uniform(size=20),
                        "B": np.random.uniform(size=20),
                        "C": np.arange(20) + np.random.uniform(size=20)})

        # Default to BuGn
        ax = df.plot(kind='hexbin', x='A', y='B')
        self.assertEquals(ax.collections[0].cmap.name, 'BuGn')

        cm = 'cubehelix'
        ax = df.plot(kind='hexbin', x='A', y='B', colormap=cm)
        self.assertEquals(ax.collections[0].cmap.name, cm)

    @slow
    def test_no_color_bar(self):
        df = DataFrame({"A": np.random.uniform(size=20),
                        "B": np.random.uniform(size=20),
                        "C": np.arange(20) + np.random.uniform(size=20)})

        ax = df.plot(kind='hexbin', x='A', y='B', colorbar=None)
        self.assertIs(ax.collections[0].colorbar, None)

    @slow
    def test_allow_cmap(self):
        df = DataFrame({"A": np.random.uniform(size=20),
                        "B": np.random.uniform(size=20),
                        "C": np.arange(20) + np.random.uniform(size=20)})

        ax = df.plot(kind='hexbin', x='A', y='B', cmap='YlGn')
        self.assertEquals(ax.collections[0].cmap.name, 'YlGn')

        with tm.assertRaises(TypeError):
            df.plot(kind='hexbin', x='A', y='B', cmap='YlGn',
                         colormap='BuGn')

    @slow
    def test_pie_df(self):
        df = DataFrame(np.random.rand(5, 3), columns=['X', 'Y', 'Z'],
                       index=['a', 'b', 'c', 'd', 'e'])
        with tm.assertRaises(ValueError):
            df.plot(kind='pie')

        ax = _check_plot_works(df.plot, kind='pie', y='Y')
        for t, expected in zip(ax.texts, df.index):
            self.assertEqual(t.get_text(), expected)

        axes = _check_plot_works(df.plot, kind='pie', subplots=True)
        self.assertEqual(len(axes), len(df.columns))
        for ax in axes:
            for t, expected in zip(ax.texts, df.index):
                self.assertEqual(t.get_text(), expected)
        for ax, ylabel in zip(axes, df.columns):
            self.assertEqual(ax.get_ylabel(), ylabel)

        labels = ['A', 'B', 'C', 'D', 'E']
        color_args = ['r', 'g', 'b', 'c', 'm']
        axes = _check_plot_works(df.plot, kind='pie', subplots=True,
                                 labels=labels, colors=color_args)
        self.assertEqual(len(axes), len(df.columns))

        import matplotlib.colors as colors
        conv = colors.colorConverter
        for ax in axes:
            for t, expected in zip(ax.texts, labels):
                self.assertEqual(t.get_text(), expected)
            for p, expected in zip(ax.patches, color_args):
                self.assertEqual(p.get_facecolor(), conv.to_rgba(expected))

    def test_errorbar_plot(self):

        d = {'x': np.arange(12), 'y': np.arange(12, 0, -1)}
        df = DataFrame(d)
        d_err = {'x': np.ones(12)*0.2, 'y': np.ones(12)*0.4}
        df_err = DataFrame(d_err)

        # check line plots
        ax = _check_plot_works(df.plot, yerr=df_err, logy=True)
        _check_has_errorbars(self, ax, xerr=0, yerr=2)
        ax = _check_plot_works(df.plot, yerr=df_err, logx=True, logy=True)
        _check_has_errorbars(self, ax, xerr=0, yerr=2)
        ax = _check_plot_works(df.plot, yerr=df_err, loglog=True)
        _check_has_errorbars(self, ax, xerr=0, yerr=2)

        kinds = ['line', 'bar', 'barh']
        for kind in kinds:
            ax = _check_plot_works(df.plot, yerr=df_err['x'], kind=kind)
            _check_has_errorbars(self, ax, xerr=0, yerr=2)
            ax = _check_plot_works(df.plot, yerr=d_err, kind=kind)
            _check_has_errorbars(self, ax, xerr=0, yerr=2)
            ax = _check_plot_works(df.plot, yerr=df_err, xerr=df_err, kind=kind)
            _check_has_errorbars(self, ax, xerr=2, yerr=2)
            ax = _check_plot_works(df.plot, yerr=df_err['x'], xerr=df_err['x'], kind=kind)
            _check_has_errorbars(self, ax, xerr=2, yerr=2)
            ax = _check_plot_works(df.plot, xerr=0.2, yerr=0.2, kind=kind)
            _check_has_errorbars(self, ax, xerr=2, yerr=2)
            axes = _check_plot_works(df.plot, yerr=df_err, xerr=df_err, subplots=True, kind=kind)
            for ax in axes:
                _check_has_errorbars(self, ax, xerr=1, yerr=1)

        ax = _check_plot_works((df+1).plot, yerr=df_err, xerr=df_err, kind='bar', log=True)
        _check_has_errorbars(self, ax, xerr=2, yerr=2)

        # yerr is raw error values
        ax = _check_plot_works(df['y'].plot, yerr=np.ones(12)*0.4)
        _check_has_errorbars(self, ax, xerr=0, yerr=1)
        ax = _check_plot_works(df.plot, yerr=np.ones((2, 12))*0.4)
        _check_has_errorbars(self, ax, xerr=0, yerr=2)

        # yerr is iterator
        import itertools
        ax = _check_plot_works(df.plot, yerr=itertools.repeat(0.1, len(df)))
        _check_has_errorbars(self, ax, xerr=0, yerr=2)

        # yerr is column name
        for yerr in ['yerr', u('誤差')]:
            s_df = df.copy()
            s_df[yerr] = np.ones(12)*0.2
            ax = _check_plot_works(s_df.plot, yerr=yerr)
            _check_has_errorbars(self, ax, xerr=0, yerr=2)
            ax = _check_plot_works(s_df.plot, y='y', x='x', yerr=yerr)
            _check_has_errorbars(self, ax, xerr=0, yerr=1)

        with tm.assertRaises(ValueError):
            df.plot(yerr=np.random.randn(11))

        df_err = DataFrame({'x': ['zzz']*12, 'y': ['zzz']*12})
        with tm.assertRaises(TypeError):
            df.plot(yerr=df_err)

    @slow
    def test_errorbar_with_integer_column_names(self):
        # test with integer column names
        df = DataFrame(np.random.randn(10, 2))
        df_err = DataFrame(np.random.randn(10, 2))
        ax = _check_plot_works(df.plot, yerr=df_err)
        _check_has_errorbars(self, ax, xerr=0, yerr=2)
        ax = _check_plot_works(df.plot, y=0, yerr=1)
        _check_has_errorbars(self, ax, xerr=0, yerr=1)

    @slow
    def test_errorbar_with_partial_columns(self):
        df = DataFrame(np.random.randn(10, 3))
        df_err = DataFrame(np.random.randn(10, 2), columns=[0, 2])
        kinds = ['line', 'bar']
        for kind in kinds:
            ax = _check_plot_works(df.plot, yerr=df_err, kind=kind)
            _check_has_errorbars(self, ax, xerr=0, yerr=2)

        ix = date_range('1/1/2000', periods=10, freq='M')
        df.set_index(ix, inplace=True)
        df_err.set_index(ix, inplace=True)
        ax = _check_plot_works(df.plot, yerr=df_err, kind='line')
        _check_has_errorbars(self, ax, xerr=0, yerr=2)

        d = {'x': np.arange(12), 'y': np.arange(12, 0, -1)}
        df = DataFrame(d)
        d_err = {'x': np.ones(12)*0.2, 'z': np.ones(12)*0.4}
        df_err = DataFrame(d_err)
        for err in [d_err, df_err]:
            ax = _check_plot_works(df.plot, yerr=err)
            _check_has_errorbars(self, ax, xerr=0, yerr=1)            

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
            _check_has_errorbars(self, ax, xerr=0, yerr=2)
            ax = _check_plot_works(tdf.plot, yerr=d_err, kind=kind)
            _check_has_errorbars(self, ax, xerr=0, yerr=2)
            ax = _check_plot_works(tdf.plot, y='y', yerr=tdf_err['x'], kind=kind)
            _check_has_errorbars(self, ax, xerr=0, yerr=1)
            ax = _check_plot_works(tdf.plot, y='y', yerr='x', kind=kind)
            _check_has_errorbars(self, ax, xerr=0, yerr=1)
            ax = _check_plot_works(tdf.plot, yerr=tdf_err, kind=kind)
            _check_has_errorbars(self, ax, xerr=0, yerr=2)
            axes = _check_plot_works(tdf.plot, kind=kind, yerr=tdf_err, subplots=True)
            for ax in axes:
                _check_has_errorbars(self, ax, xerr=0, yerr=1)

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
        self.assert_(len(ax.tables) == 0)
        plotting.table(ax, df.T)
        self.assert_(len(ax.tables) == 1)

    def test_errorbar_scatter(self):
        df = DataFrame(np.random.randn(5, 2), index=range(5), columns=['x', 'y'])
        df_err = DataFrame(np.random.randn(5, 2) / 5,
                           index=range(5), columns=['x', 'y'])
            
        ax = _check_plot_works(df.plot, kind='scatter', x='x', y='y')
        _check_has_errorbars(self, ax, xerr=0, yerr=0)
        ax = _check_plot_works(df.plot, kind='scatter', x='x', y='y', xerr=df_err)
        _check_has_errorbars(self, ax, xerr=1, yerr=0)
        ax = _check_plot_works(df.plot, kind='scatter', x='x', y='y', yerr=df_err)
        _check_has_errorbars(self, ax, xerr=0, yerr=1)
        ax = _check_plot_works(df.plot, kind='scatter', x='x', y='y',
                               xerr=df_err, yerr=df_err)
        _check_has_errorbars(self, ax, xerr=1, yerr=1)


@tm.mplskip
class TestDataFrameGroupByPlots(tm.TestCase):
    def tearDown(self):
        tm.close()

    @slow
    def test_boxplot(self):
        df = DataFrame(np.random.rand(10, 2), columns=['Col1', 'Col2'])
        df['X'] = Series(['A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B'])
        grouped = df.groupby(by='X')
        _check_plot_works(grouped.boxplot)
        _check_plot_works(grouped.boxplot, subplots=False)

        tuples = lzip(string.ascii_letters[:10], range(10))
        df = DataFrame(np.random.rand(10, 3),
                       index=MultiIndex.from_tuples(tuples))
        grouped = df.groupby(level=1)
        _check_plot_works(grouped.boxplot)
        _check_plot_works(grouped.boxplot, subplots=False)

        grouped = df.unstack(level=1).groupby(level=0, axis=1)
        _check_plot_works(grouped.boxplot)
        _check_plot_works(grouped.boxplot, subplots=False)

    def test_series_plot_color_kwargs(self):
        # GH1890
        ax = Series(np.arange(12) + 1).plot(color='green')
        line = ax.get_lines()[0]
        self.assertEqual(line.get_color(), 'green')

    def test_time_series_plot_color_kwargs(self):
        # #1890
        ax = Series(np.arange(12) + 1, index=date_range(
            '1/1/2000', periods=12)).plot(color='green')
        line = ax.get_lines()[0]
        self.assertEqual(line.get_color(), 'green')

    def test_time_series_plot_color_with_empty_kwargs(self):
        import matplotlib as mpl

        def_colors = mpl.rcParams['axes.color_cycle']
        index = date_range('1/1/2000', periods=12)
        s = Series(np.arange(1, 13), index=index)

        ncolors = 3

        for i in range(ncolors):
            ax = s.plot()

        line_colors = [l.get_color() for l in ax.get_lines()]
        self.assertEqual(line_colors, def_colors[:ncolors])

    @slow
    def test_grouped_hist(self):
        import matplotlib.pyplot as plt
        df = DataFrame(randn(500, 2), columns=['A', 'B'])
        df['C'] = np.random.randint(0, 4, 500)
        axes = plotting.grouped_hist(df.A, by=df.C)
        self.assertEqual(len(axes.ravel()), 4)

        tm.close()
        axes = df.hist(by=df.C)
        self.assertEqual(axes.ndim, 2)
        self.assertEqual(len(axes.ravel()), 4)

        for ax in axes.ravel():
            self.assert_(len(ax.patches) > 0)

        tm.close()
        # make sure kwargs to hist are handled
        axes = plotting.grouped_hist(df.A, by=df.C, normed=True,
                                     cumulative=True, bins=4)

        # height of last bin (index 5) must be 1.0
        for ax in axes.ravel():
            height = ax.get_children()[5].get_height()
            self.assertAlmostEqual(height, 1.0)

        tm.close()
        axes = plotting.grouped_hist(df.A, by=df.C, log=True)
        # scale of y must be 'log'
        for ax in axes.ravel():
            self.assertEqual(ax.get_yscale(), 'log')

        tm.close()
        # propagate attr exception from matplotlib.Axes.hist
        with tm.assertRaises(AttributeError):
            plotting.grouped_hist(df.A, by=df.C, foo='bar')

    @slow
    def test_grouped_hist_layout(self):
        import matplotlib.pyplot as plt
        n = 100
        gender = tm.choice(['Male', 'Female'], size=n)
        df = DataFrame({'gender': gender,
                        'height': random.normal(66, 4, size=n),
                        'weight': random.normal(161, 32, size=n),
                        'category': random.randint(4, size=n)})
        self.assertRaises(ValueError, df.hist, column='weight', by=df.gender,
                          layout=(1, 1))
        self.assertRaises(ValueError, df.hist, column='weight', by=df.gender,
                          layout=(1,))
        self.assertRaises(ValueError, df.hist, column='height', by=df.category,
                          layout=(1, 3))
        self.assertRaises(ValueError, df.hist, column='height', by=df.category,
                          layout=(2, 1))
        self.assertEqual(df.hist(column='height', by=df.gender,
                                 layout=(2, 1)).shape, (2,))
        tm.close()
        self.assertEqual(df.hist(column='height', by=df.category,
                                 layout=(4, 1)).shape, (4,))
        tm.close()
        self.assertEqual(df.hist(column='height', by=df.category,
                                 layout=(4, 2)).shape, (4, 2))

    @slow
    def test_axis_share_x(self):
        # GH4089
        n = 100
        df = DataFrame({'gender': tm.choice(['Male', 'Female'], size=n),
                        'height': random.normal(66, 4, size=n),
                        'weight': random.normal(161, 32, size=n)})
        ax1, ax2 = df.hist(column='height', by=df.gender, sharex=True)

        # share x
        self.assertTrue(ax1._shared_x_axes.joined(ax1, ax2))
        self.assertTrue(ax2._shared_x_axes.joined(ax1, ax2))

        # don't share y
        self.assertFalse(ax1._shared_y_axes.joined(ax1, ax2))
        self.assertFalse(ax2._shared_y_axes.joined(ax1, ax2))

    @slow
    def test_axis_share_y(self):
        n = 100
        df = DataFrame({'gender': tm.choice(['Male', 'Female'], size=n),
                        'height': random.normal(66, 4, size=n),
                        'weight': random.normal(161, 32, size=n)})
        ax1, ax2 = df.hist(column='height', by=df.gender, sharey=True)

        # share y
        self.assertTrue(ax1._shared_y_axes.joined(ax1, ax2))
        self.assertTrue(ax2._shared_y_axes.joined(ax1, ax2))

        # don't share x
        self.assertFalse(ax1._shared_x_axes.joined(ax1, ax2))
        self.assertFalse(ax2._shared_x_axes.joined(ax1, ax2))

    @slow
    def test_axis_share_xy(self):
        n = 100
        df = DataFrame({'gender': tm.choice(['Male', 'Female'], size=n),
                        'height': random.normal(66, 4, size=n),
                        'weight': random.normal(161, 32, size=n)})
        ax1, ax2 = df.hist(column='height', by=df.gender, sharex=True,
                           sharey=True)

        # share both x and y
        self.assertTrue(ax1._shared_x_axes.joined(ax1, ax2))
        self.assertTrue(ax2._shared_x_axes.joined(ax1, ax2))

        self.assertTrue(ax1._shared_y_axes.joined(ax1, ax2))
        self.assertTrue(ax2._shared_y_axes.joined(ax1, ax2))

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


def _check_has_errorbars(t, ax, xerr=0, yerr=0):
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
    t.assertEqual(xerr, xerr_count)
    t.assertEqual(yerr, yerr_count)


def _check_plot_works(f, *args, **kwargs):
    import matplotlib.pyplot as plt
    ret = None

    try:
        try:
            fig = kwargs['figure']
        except KeyError:
            fig = plt.gcf()

        plt.clf()

        ax = kwargs.get('ax', fig.add_subplot(211))
        ret = f(*args, **kwargs)

        assert_is_valid_plot_return_object(ret)

        try:
            kwargs['ax'] = fig.add_subplot(212)
            ret = f(*args, **kwargs)
        except Exception:
            pass
        else:
            assert_is_valid_plot_return_object(ret)

        with ensure_clean(return_filelike=True) as path:
            plt.savefig(path)
    finally:
        tm.close(fig)

    return ret


def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
