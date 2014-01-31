import nose
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
from numpy.random import randn

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
        _check_plot_works(self.iseries.plot)
        _check_plot_works(self.series[:5].plot, kind='bar')
        _check_plot_works(self.series[:5].plot, kind='line')
        _check_plot_works(self.series[:5].plot, kind='barh')
        _check_plot_works(self.series[:10].plot, kind='barh')
        _check_plot_works(Series(randn(10)).plot, kind='bar', color='black')

    @slow
    def test_plot_figsize_and_title(self):
        # figsize and title
        import matplotlib.pyplot as plt
        ax = self.series.plot(title='Test', figsize=(16, 8))

        self.assertEqual(ax.title.get_text(), 'Test')
        assert_array_equal(np.round(ax.figure.get_size_inches()),
                           np.array((16., 8.)))

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
    def test_hist(self):
        _check_plot_works(self.ts.hist)
        _check_plot_works(self.ts.hist, grid=False)
        _check_plot_works(self.ts.hist, figsize=(8, 10))
        _check_plot_works(self.ts.hist, by=self.ts.index.month)

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
        df = DataFrame({'gender': tm.choice(['Male', 'Female'], size=n),
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
        kinds = 'line', 'bar', 'barh', 'kde', 'density'

        for kind in kinds:
            with tm.assertRaises(TypeError):
                s.plot(kind=kind)

    @slow
    def test_valid_object_plot(self):
        s = Series(lrange(10), dtype=object)
        kinds = 'line', 'bar', 'barh', 'kde', 'density'

        for kind in kinds:
            _check_plot_works(s.plot, kind=kind)

    def test_partially_invalid_plot_data(self):
        s = Series(['a', 'b', 1.0, 2])
        kinds = 'line', 'bar', 'barh', 'kde', 'density'

        for kind in kinds:
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
    def test_xcompat(self):
        import pandas as pd
        import matplotlib.pyplot as plt

        df = tm.makeTimeDataFrame()
        ax = df.plot(x_compat=True)
        lines = ax.get_lines()
        self.assert_(not isinstance(lines[0].get_xdata(), PeriodIndex))

        tm.close()
        pd.plot_params['xaxis.compat'] = True
        ax = df.plot()
        lines = ax.get_lines()
        self.assert_(not isinstance(lines[0].get_xdata(), PeriodIndex))

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
            self.assert_(not isinstance(lines[0].get_xdata(), PeriodIndex))

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

        axes = df.plot(subplots=True, sharex=True, legend=True)

        for ax in axes:
            self.assert_(ax.get_legend() is not None)

        axes = df.plot(subplots=True, sharex=True)
        for ax in axes[:-2]:
            [self.assert_(not label.get_visible())
             for label in ax.get_xticklabels()]
            [self.assert_(label.get_visible())
             for label in ax.get_yticklabels()]

        [self.assert_(label.get_visible())
         for label in axes[-1].get_xticklabels()]
        [self.assert_(label.get_visible())
         for label in axes[-1].get_yticklabels()]

        axes = df.plot(subplots=True, sharex=False)
        for ax in axes:
            [self.assert_(label.get_visible())
             for label in ax.get_xticklabels()]
            [self.assert_(label.get_visible())
             for label in ax.get_yticklabels()]

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

    def test_bar_stacked_center(self):
        # GH2157
        df = DataFrame({'A': [3] * 5, 'B': lrange(5)}, index=lrange(5))
        ax = df.plot(kind='bar', stacked='True', grid=True)
        self.assertEqual(ax.xaxis.get_ticklocs()[0],
                         ax.patches[0].get_x() + ax.patches[0].get_width() / 2)

    def test_bar_center(self):
        df = DataFrame({'A': [3] * 5, 'B': lrange(5)}, index=lrange(5))
        ax = df.plot(kind='bar', grid=True)
        self.assertEqual(ax.xaxis.get_ticklocs()[0],
                         ax.patches[0].get_x() + ax.patches[0].get_width())

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

    @slow
    def test_kde(self):
        _skip_if_no_scipy()
        df = DataFrame(randn(100, 4))
        _check_plot_works(df.plot, kind='kde')
        _check_plot_works(df.plot, kind='kde', subplots=True)
        ax = df.plot(kind='kde')
        self.assert_(ax.get_legend() is not None)
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

        path = os.path.join(curpath(), 'data', 'iris.csv')
        df = read_csv(path)

        _check_plot_works(andrews_curves, df, 'Name')

    @slow
    def test_parallel_coordinates(self):
        from pandas import read_csv
        from pandas.tools.plotting import parallel_coordinates
        from matplotlib import cm
        path = os.path.join(curpath(), 'data', 'iris.csv')
        df = read_csv(path)
        _check_plot_works(parallel_coordinates, df, 'Name')
        _check_plot_works(parallel_coordinates, df, 'Name',
                          colors=('#556270', '#4ECDC4', '#C7F464'))
        _check_plot_works(parallel_coordinates, df, 'Name',
                          colors=['dodgerblue', 'aquamarine', 'seagreen'])
        _check_plot_works(parallel_coordinates, df, 'Name',
                          colors=('#556270', '#4ECDC4', '#C7F464'))
        _check_plot_works(parallel_coordinates, df, 'Name',
                          colors=['dodgerblue', 'aquamarine', 'seagreen'])
        _check_plot_works(parallel_coordinates, df, 'Name', colormap=cm.jet)

        df = read_csv(path, header=None, skiprows=1, names=[1, 2, 4, 8,
                                                            'Name'])
        _check_plot_works(parallel_coordinates, df, 'Name', use_columns=True)
        _check_plot_works(parallel_coordinates, df, 'Name',
                          xticks=[1, 5, 25, 125])

    @slow
    def test_radviz(self):
        from pandas import read_csv
        from pandas.tools.plotting import radviz
        from matplotlib import cm

        path = os.path.join(curpath(), 'data', 'iris.csv')
        df = read_csv(path)
        _check_plot_works(radviz, df, 'Name')
        _check_plot_works(radviz, df, 'Name', colormap=cm.jet)

    @slow
    def test_plot_int_columns(self):
        df = DataFrame(randn(100, 4)).cumsum()
        _check_plot_works(df.plot, legend=True)

    def test_legend_name(self):
        multi = DataFrame(randn(4, 4),
                          columns=[np.array(['a', 'a', 'b', 'b']),
                                   np.array(['x', 'y', 'x', 'y'])])
        multi.columns.names = ['group', 'individual']

        ax = multi.plot()
        leg_title = ax.legend_.get_title()
        self.assertEqual(leg_title.get_text(), 'group,individual')

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

    @slow
    def test_line_colors(self):
        import matplotlib.pyplot as plt
        import sys
        from matplotlib import cm

        custom_colors = 'rgcby'

        df = DataFrame(randn(5, 5))

        ax = df.plot(color=custom_colors)

        lines = ax.get_lines()
        for i, l in enumerate(lines):
            xp = custom_colors[i]
            rs = l.get_color()
            self.assertEqual(xp, rs)

        tmp = sys.stderr
        sys.stderr = StringIO()
        try:
            tm.close()
            ax2 = df.plot(colors=custom_colors)
            lines2 = ax2.get_lines()
            for l1, l2 in zip(lines, lines2):
                self.assertEqual(l1.get_color(), l2.get_color())
        finally:
            sys.stderr = tmp

        tm.close()

        ax = df.plot(colormap='jet')

        rgba_colors = lmap(cm.jet, np.linspace(0, 1, len(df)))

        lines = ax.get_lines()
        for i, l in enumerate(lines):
            xp = rgba_colors[i]
            rs = l.get_color()
            self.assertEqual(xp, rs)

        tm.close()

        ax = df.plot(colormap=cm.jet)

        rgba_colors = lmap(cm.jet, np.linspace(0, 1, len(df)))

        lines = ax.get_lines()
        for i, l in enumerate(lines):
            xp = rgba_colors[i]
            rs = l.get_color()
            self.assertEqual(xp, rs)

        # make color a list if plotting one column frame
        # handles cases like df.plot(color='DodgerBlue')
        tm.close()
        df.ix[:, [0]].plot(color='DodgerBlue')

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
        kinds = 'line', 'bar', 'barh', 'kde', 'density'
        df = DataFrame(list('abcd'))
        for kind in kinds:
            with tm.assertRaises(TypeError):
                df.plot(kind=kind)

    @slow
    def test_partially_invalid_plot_data(self):
        kinds = 'line', 'bar', 'barh', 'kde', 'density'
        df = DataFrame(randn(10, 2), dtype=object)
        df[np.random.rand(df.shape[0]) > 0.5] = 'a'
        for kind in kinds:
            with tm.assertRaises(TypeError):
                df.plot(kind=kind)

    def test_invalid_kind(self):
        df = DataFrame(randn(10, 2))
        with tm.assertRaises(ValueError):
            df.plot(kind='aasdf')


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


def _check_plot_works(f, *args, **kwargs):
    import matplotlib.pyplot as plt

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


def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
