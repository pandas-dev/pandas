import nose
import os
import unittest

from pandas import Series, DataFrame
import pandas.util.testing as tm

import numpy as np

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

    def test_plot(self):
        _check_plot_works(self.ts.plot, label='foo')
        _check_plot_works(self.ts.plot, use_index=False)
        _check_plot_works(self.ts.plot, rot=0)
        _check_plot_works(self.ts.plot, style='.')
        _check_plot_works(self.ts[:10].plot, kind='bar')
        _check_plot_works(self.series[:5].plot, kind='bar')

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

    def test_plot(self):
        df = tm.makeTimeDataFrame()
        _check_plot_works(df.plot, grid=False)
        _check_plot_works(df.plot, subplots=True)
        _check_plot_works(df.plot, subplots=True, use_index=False)

    def test_hist(self):
        df = DataFrame(np.random.randn(100, 4))
        _check_plot_works(df.hist)
        _check_plot_works(df.hist, grid=False)

    def test_plot_int_columns(self):
        df = DataFrame(np.random.randn(100, 4)).cumsum()
        _check_plot_works(df.plot, legend=True)


PNG_PATH = 'tmp.png'

def _check_plot_works(f, *args, **kwargs):
    import matplotlib.pyplot as plt

    fig = plt.gcf()
    plt.clf()
    ax = fig.add_subplot(211)
    f(*args, **kwargs)

    ax = fig.add_subplot(212)
    try:
        kwargs['ax'] = ax
        f(*args, **kwargs)
    except Exception:
        pass
    plt.savefig(PNG_PATH)
    os.remove(PNG_PATH)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
