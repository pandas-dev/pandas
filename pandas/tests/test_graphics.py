import nose
import os
import string
import unittest

from datetime import datetime

from pandas import Series, DataFrame, MultiIndex, PeriodIndex, date_range
import pandas.util.testing as tm

import numpy as np

from numpy.testing import assert_array_equal
from numpy.testing.decorators import slow
import pandas.tools.plotting as plotting

class TestSeriesPlots(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        import sys
        if 'IPython' in sys.modules:
            raise nose.SkipTest

        try:
            import matplotlib as mpl
            mpl.use('Agg', warn=False)
        except ImportError:
            raise nose.SkipTest

    def setUp(self):
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
        _check_plot_works(self.ts.plot, rot=0)
        _check_plot_works(self.ts.plot, style='.', logy=True)
        _check_plot_works(self.ts.plot, style='.', logx=True)
        _check_plot_works(self.ts.plot, style='.', loglog=True)
        _check_plot_works(self.ts[:10].plot, kind='bar')
        _check_plot_works(self.series[:5].plot, kind='bar')
        _check_plot_works(self.series[:5].plot, kind='line')
        _check_plot_works(self.series[:5].plot, kind='barh')
        _check_plot_works(self.series[:10].plot, kind='barh')

        Series(np.random.randn(10)).plot(kind='bar',color='black')

    @slow
    def test_bar_colors(self):
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors

        default_colors = 'brgyk'
        custom_colors = 'rgcby'

        plt.close('all')
        df = DataFrame(np.random.randn(5, 5))
        ax = df.plot(kind='bar')

        rects = ax.patches

        conv = colors.colorConverter
        for i, rect in enumerate(rects[::5]):
            xp = conv.to_rgba(default_colors[i])
            rs = rect.get_facecolor()
            self.assert_(xp == rs)

        plt.close('all')

        ax = df.plot(kind='bar', color=custom_colors)

        rects = ax.patches

        conv = colors.colorConverter
        for i, rect in enumerate(rects[::5]):
            xp = conv.to_rgba(custom_colors[i])
            rs = rect.get_facecolor()
            self.assert_(xp == rs)

    @slow
    def test_bar_linewidth(self):
        df = DataFrame(np.random.randn(5, 5))

        # regular
        ax = df.plot(kind='bar', linewidth=2)
        for r in ax.patches:
            self.assert_(r.get_linewidth() == 2)

        # stacked
        ax = df.plot(kind='bar', stacked=True, linewidth=2)
        for r in ax.patches:
            self.assert_(r.get_linewidth() == 2)

        # subplots
        axes = df.plot(kind='bar', linewidth=2, subplots=True)
        for ax in axes:
            for r in ax.patches:
                self.assert_(r.get_linewidth() == 2)

    @slow
    def test_1rotation(self):
        df = DataFrame(np.random.randn(5, 5))
        ax = df.plot(rot=30)
        for l in ax.get_xticklabels():
            self.assert_(l.get_rotation() == 30)

    @slow
    def test_irregular_datetime(self):
        rng = date_range('1/1/2000', '3/1/2000')
        rng = rng[[0,1,2,3,5,9,10,11,12]]
        ser = Series(np.random.randn(len(rng)), rng)
        ax = ser.plot()
        xp = datetime(1999, 1, 1).toordinal()
        ax.set_xlim('1/1/1999', '1/1/2001')
        self.assert_(xp == ax.get_xlim()[0])

    @slow
    def test_hist(self):
        _check_plot_works(self.ts.hist)
        _check_plot_works(self.ts.hist, grid=False)

    @slow
    def test_kde(self):
        _check_plot_works(self.ts.plot, kind='kde')
        _check_plot_works(self.ts.plot, kind='density')
        ax = self.ts.plot(kind='kde', logy=True)
        self.assert_(ax.get_yscale() == 'log')

    @slow
    def test_autocorrelation_plot(self):
        from pandas.tools.plotting import autocorrelation_plot
        _check_plot_works(autocorrelation_plot, self.ts)
        _check_plot_works(autocorrelation_plot, self.ts.values)

    @slow
    def test_lag_plot(self):
        from pandas.tools.plotting import lag_plot
        _check_plot_works(lag_plot, self.ts)

    @slow
    def test_bootstrap_plot(self):
        from pandas.tools.plotting import bootstrap_plot
        _check_plot_works(bootstrap_plot, self.ts, size=10)

class TestDataFramePlots(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        import sys
        if 'IPython' in sys.modules:
            raise nose.SkipTest

        try:
            import matplotlib as mpl
            mpl.use('Agg', warn=False)
        except ImportError:
            raise nose.SkipTest

    @slow
    def test_plot(self):
        df = tm.makeTimeDataFrame()
        _check_plot_works(df.plot, grid=False)
        _check_plot_works(df.plot, subplots=True)
        _check_plot_works(df.plot, subplots=True, use_index=False)

        df = DataFrame({'x':[1,2], 'y':[3,4]})
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

        tuples = zip(list(string.ascii_letters[:10]), range(10))
        df = DataFrame(np.random.rand(10, 3),
                       index=MultiIndex.from_tuples(tuples))
        _check_plot_works(df.plot, use_index=True)

        # unicode
        index = MultiIndex.from_tuples([(u'\u03b1', 0),
                                        (u'\u03b1', 1),
                                        (u'\u03b2', 2),
                                        (u'\u03b2', 3),
                                        (u'\u03b3', 4),
                                        (u'\u03b3', 5),
                                        (u'\u03b4', 6),
                                        (u'\u03b4', 7)], names=['i0', 'i1'])
        columns = MultiIndex.from_tuples([('bar', u'\u0394'),
            ('bar', u'\u0395')], names=['c0', 'c1'])
        df = DataFrame(np.random.randint(0, 10, (8, 2)),
                       columns=columns,
                       index=index)
        _check_plot_works(df.plot, title=u'\u03A3')

    @slow
    def test_plot_xy(self):
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
        df.columns = range(1, len(df.columns) + 1)
        self._check_data(df.plot(x=1, y=2),
                         df.set_index(1)[2].plot())
        self._check_data(df.plot(x=1), df.set_index(1).plot())
        self._check_data(df.plot(y=1), df[1].plot())

        # columns.inferred_type == 'mixed'
        # TODO add MultiIndex test

    def _check_data(self, xp, rs):
        xp_lines = xp.get_lines()
        rs_lines = rs.get_lines()

        def check_line(xpl, rsl):
            xpdata = xpl.get_xydata()
            rsdata = rsl.get_xydata()
            assert_array_equal(xpdata, rsdata)

        [check_line(xpl, rsl) for xpl, rsl in zip(xp_lines, rs_lines)]

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
    def test_plot_bar(self):
        df = DataFrame(np.random.randn(6, 4),
                       index=list(string.ascii_letters[:6]),
                       columns=['one', 'two', 'three', 'four'])

        _check_plot_works(df.plot, kind='bar')
        _check_plot_works(df.plot, kind='bar', legend=False)
        _check_plot_works(df.plot, kind='bar', subplots=True)
        _check_plot_works(df.plot, kind='bar', stacked=True)

        df = DataFrame(np.random.randn(10, 15),
                       index=list(string.ascii_letters[:10]),
                       columns=range(15))
        _check_plot_works(df.plot, kind='bar')

        df = DataFrame({'a': [0, 1], 'b': [1, 0]})
        _check_plot_works(df.plot, kind='bar')

    @slow
    def test_boxplot(self):
        df = DataFrame(np.random.randn(6, 4),
                       index=list(string.ascii_letters[:6]),
                       columns=['one', 'two', 'three', 'four'])
        df['indic'] = ['foo', 'bar'] * 3
        df['indic2'] = ['foo', 'bar', 'foo'] * 2

        _check_plot_works(df.boxplot)
        _check_plot_works(df.boxplot, column=['one', 'two'])
        _check_plot_works(df.boxplot, column=['one', 'two'],
                          by='indic')
        _check_plot_works(df.boxplot, column='one', by=['indic', 'indic2'])
        _check_plot_works(df.boxplot, by='indic')
        _check_plot_works(df.boxplot, by=['indic', 'indic2'])

        _check_plot_works(lambda x: plotting.boxplot(x), df['one'])

        _check_plot_works(df.boxplot, notch=1)
        _check_plot_works(df.boxplot, by='indic', notch=1)

        df = DataFrame(np.random.rand(10,2), columns=['Col1', 'Col2'] )
        df['X'] = Series(['A','A','A','A','A','B','B','B','B','B'])
        _check_plot_works(df.boxplot, by='X')

    @slow
    def test_kde(self):
        df = DataFrame(np.random.randn(100, 4))
        _check_plot_works(df.plot, kind='kde')
        _check_plot_works(df.plot, kind='kde', subplots=True)
        axes = df.plot(kind='kde', logy=True, subplots=True)
        for ax in axes:
            self.assert_(ax.get_yscale() == 'log')

    @slow
    def test_hist(self):
        df = DataFrame(np.random.randn(100, 4))
        _check_plot_works(df.hist)
        _check_plot_works(df.hist, grid=False)

        #make sure layout is handled
        df = DataFrame(np.random.randn(100, 3))
        _check_plot_works(df.hist)
        axes = df.hist(grid=False)
        self.assert_(not axes[1, 1].get_visible())

        df = DataFrame(np.random.randn(100, 1))
        _check_plot_works(df.hist)

        #make sure layout is handled
        df = DataFrame(np.random.randn(100, 6))
        _check_plot_works(df.hist)

        #make sure sharex, sharey is handled
        _check_plot_works(df.hist, sharex=True, sharey=True)

        #make sure kwargs are handled
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

    @slow
    def test_scatter(self):
        df = DataFrame(np.random.randn(100, 4))
        import pandas.tools.plotting as plt
        def scat(**kwds):
            return plt.scatter_matrix(df, **kwds)
        _check_plot_works(scat)
        _check_plot_works(scat, marker='+')
        _check_plot_works(scat, vmin=0)
        _check_plot_works(scat, diagonal='kde')
        _check_plot_works(scat, diagonal='density')
        _check_plot_works(scat, diagonal='hist')

        def scat2(x, y, by=None, ax=None, figsize=None):
            return plt.scatter_plot(df, x, y, by, ax, figsize=None)

        _check_plot_works(scat2, 0, 1)
        grouper = Series(np.repeat([1, 2, 3, 4, 5], 20), df.index)
        _check_plot_works(scat2, 0, 1, by=grouper)

    @slow
    def test_andrews_curves(self):
        from pandas import read_csv
        from pandas.tools.plotting import andrews_curves
        path = os.path.join(curpath(), 'data/iris.csv')
        df = read_csv(path)
        _check_plot_works(andrews_curves, df, 'Name')

    @slow
    def test_parallel_coordinates(self):
        from pandas import read_csv
        from pandas.tools.plotting import parallel_coordinates
        path = os.path.join(curpath(), 'data/iris.csv')
        df = read_csv(path)
        _check_plot_works(parallel_coordinates, df, 'Name')

    @slow
    def test_radviz(self):
        from pandas import read_csv
        from pandas.tools.plotting import radviz
        path = os.path.join(curpath(), 'data/iris.csv')
        df = read_csv(path)
        _check_plot_works(radviz, df, 'Name')

    @slow
    def test_plot_int_columns(self):
        df = DataFrame(np.random.randn(100, 4)).cumsum()
        _check_plot_works(df.plot, legend=True)

    @slow
    def test_legend_name(self):
        multi = DataFrame(np.random.randn(4, 4),
                          columns=[np.array(['a', 'a', 'b', 'b']),
                                   np.array(['x', 'y', 'x', 'y'])])
        multi.columns.names = ['group', 'individual']

        ax = multi.plot()
        leg_title = ax.legend_.get_title()
        self.assert_(leg_title.get_text(), 'group,individual')

    def _check_plot_fails(self, f, *args, **kwargs):
        self.assertRaises(Exception, f, *args, **kwargs)

    @slow
    def test_style_by_column(self):
        import matplotlib.pyplot as plt
        fig = plt.gcf()

        df = DataFrame(np.random.randn(100, 3))
        for markers in [{0: '^', 1: '+', 2: 'o'},
                        {0: '^', 1: '+'},
                        ['^', '+', 'o'],
                        ['^', '+']]:
            fig.clf()
            fig.add_subplot(111)
            ax = df.plot(style=markers)
            for i, l in enumerate(ax.get_lines()[:len(markers)]):
                self.assertEqual(l.get_marker(), markers[i])

class TestDataFrameGroupByPlots(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        import sys
        if 'IPython' in sys.modules:
            raise nose.SkipTest

        try:
            import matplotlib as mpl
            mpl.use('Agg', warn=False)
        except ImportError:
            raise nose.SkipTest

    @slow
    def test_boxplot(self):
        df = DataFrame(np.random.rand(10,2), columns=['Col1', 'Col2'] )
        df['X'] = Series(['A','A','A','A','A','B','B','B','B','B'])
        grouped = df.groupby(by='X')
        _check_plot_works(grouped.boxplot)
        _check_plot_works(grouped.boxplot, subplots=False)

        tuples = zip(list(string.ascii_letters[:10]), range(10))
        df = DataFrame(np.random.rand(10, 3),
                       index=MultiIndex.from_tuples(tuples))
        grouped = df.groupby(level=1)
        _check_plot_works(grouped.boxplot)
        _check_plot_works(grouped.boxplot, subplots=False)
        grouped = df.unstack(level=1).groupby(level=0, axis=1)
        _check_plot_works(grouped.boxplot)
        _check_plot_works(grouped.boxplot, subplots=False)

    @slow
    def test_series_plot_color_kwargs(self):
        # #1890
        import matplotlib.pyplot as plt

        plt.close('all')
        ax = Series(np.arange(12) + 1).plot(color='green')
        line = ax.get_lines()[0]
        self.assert_(line.get_color() == 'green')

PNG_PATH = 'tmp.png'

def _check_plot_works(f, *args, **kwargs):
    import matplotlib.pyplot as plt

    fig = plt.gcf()
    plt.clf()
    ax = fig.add_subplot(211)
    ret = f(*args, **kwargs)
    assert(ret is not None)  # do something more intelligent

    ax = fig.add_subplot(212)
    try:
        kwargs['ax'] = ax
        ret = f(*args, **kwargs)
        assert(ret is not None)  # do something more intelligent
    except Exception:
        pass
    plt.savefig(PNG_PATH)
    os.remove(PNG_PATH)

def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth

if __name__ == '__main__':
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
