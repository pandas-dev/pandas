import nose
import os
import string
import unittest

from pandas import Series, DataFrame, MultiIndex, PeriodIndex
import pandas.util.testing as tm

import numpy as np

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
    def test_hist(self):
        _check_plot_works(self.ts.hist)
        _check_plot_works(self.ts.hist, grid=False)


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
        _check_plot_works(scat, diagonal='hist')

        def scat2(x, y, by=None, ax=None, figsize=None):
            return plt.scatter_plot(df, x, y, by, ax, figsize=None)

        _check_plot_works(scat2, 0, 1)
        grouper = Series(np.repeat([1, 2, 3, 4, 5], 20), df.index)
        _check_plot_works(scat2, 0, 1, by=grouper)

    @slow
    def test_plot_int_columns(self):
        df = DataFrame(np.random.randn(100, 4)).cumsum()
        _check_plot_works(df.plot, legend=True)

    def _check_plot_fails(self, f, *args, **kwargs):
        self.assertRaises(Exception, f, *args, **kwargs)

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

if __name__ == '__main__':
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
