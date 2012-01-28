import nose
import os
import string
import unittest

from pandas import Series, DataFrame
import pandas.util.testing as tm

import numpy as np

from numpy.testing.decorators import slow

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

    @slow
    def test_plot(self):
        _check_plot_works(self.ts.plot, label='foo')
        _check_plot_works(self.ts.plot, use_index=False)
        _check_plot_works(self.ts.plot, rot=0)
        _check_plot_works(self.ts.plot, style='.', logy=True)
        _check_plot_works(self.ts[:10].plot, kind='bar')
        _check_plot_works(self.series[:5].plot, kind='bar')

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

    @slow
    def test_plot_bar(self):
        df = DataFrame(np.random.randn(6, 4),
                       index=list(string.ascii_letters[:6]),
                       columns=['one', 'two', 'three', 'four'])

        _check_plot_works(df.plot, kind='bar')
        _check_plot_works(df.plot, kind='bar', legend=False)
        _check_plot_works(df.plot, kind='bar', subplots=True)

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

    @slow
    def test_hist(self):
        df = DataFrame(np.random.randn(100, 4))
        _check_plot_works(df.hist)
        _check_plot_works(df.hist, grid=False)

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
