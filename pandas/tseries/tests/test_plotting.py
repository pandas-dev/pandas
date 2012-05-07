import os
from datetime import datetime, timedelta

import unittest
import nose

import numpy as np
from numpy.testing.decorators import slow

from pandas import Index, Series, DataFrame, isnull, notnull

from pandas.tseries.index import date_range
from pandas.tseries.offsets import Minute, bday
from pandas.tseries.period import period_range
from pandas.tseries.resample import DatetimeIndex, TimeGrouper
import pandas.tseries.offsets as offsets

from pandas.util.testing import assert_series_equal, assert_almost_equal
import pandas.util.testing as tm

class TestTSPlot(unittest.TestCase):

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
        freq = ['S', 'T', 'H', 'D', 'W', 'M', 'Q', 'Y']
        idx = [period_range('12/31/1999', freq=x, periods=100) for x in freq]
        self.period_ser = [Series(np.random.randn(len(x)), x) for x in idx]
        self.period_df = [DataFrame(np.random.randn(len(x), 3), index=x,
                                    columns=['A', 'B', 'C'])
                                    for x in idx]

        freq = ['S', 'T', 'H', 'D', 'W', 'M', 'Q-DEC', 'A']
        idx = [date_range('12/31/1999', freq=x, periods=100) for x in freq]
        self.datetime_ser = [Series(np.random.randn(len(x)), x) for x in idx]
        self.datetime_df = [DataFrame(np.random.randn(len(x), 3), index=x,
                                    columns=['A', 'B', 'C'])
                                    for x in idx]

    @slow
    def test_tsplot(self):
        from pandas.tseries.plotting import tsplot
        import matplotlib.pyplot as plt
        ax = plt.gca()

        f = lambda *args, **kwds: tsplot(s, plt.Axes.plot, *args, **kwds)

        for s in self.period_ser:
            _check_plot_works(f, s.index.freq, ax=ax, series=s)
        for s in self.datetime_ser:
            _check_plot_works(f, s.index.freq.rule_code, ax=ax, series=s)

    @slow
    def test_line_plot_period_series(self):
        for s in self.period_ser:
            _check_plot_works(s.plot, s.index.freq)

    @slow
    def test_line_plot_datetime_series(self):
        for s in self.datetime_ser:
            _check_plot_works(s.plot, s.index.freq.rule_code)

    @slow
    def test_line_plot_period_frame(self):
        for df in self.period_df:
            _check_plot_works(df.plot, df.index.freq)

    @slow
    def test_line_plot_datetime_frame(self):
        for df in self.datetime_df:
            freq = df.index.to_period(df.index.freq.rule_code).freq
            _check_plot_works(df.plot, freq)

    @slow
    def test_line_plot_inferred_freq(self):
        for ser in self.datetime_ser:
            ser = Series(ser.values, Index(np.asarray(ser.index)))
            _check_plot_works(ser.plot, ser.index.inferred_freq)

            ser = ser[[0, 3, 5, 6]]
            _check_plot_works(ser.plot)

    @slow
    def test_aplot_offset_freq(self):
        ser = tm.makeTimeSeries()
        _check_plot_works(ser.plot)

        dr = date_range(ser.index[0], freq='BQS', periods=10)
        ser = Series(np.random.randn(len(dr)), dr)
        _check_plot_works(ser.plot)


PNG_PATH = 'tmp.png'
def _check_plot_works(f, freq=None, series=None, *args, **kwargs):
    import matplotlib.pyplot as plt

    fig = plt.gcf()
    plt.clf()
    ax = fig.add_subplot(211)
    ret = f(*args, **kwargs)
    assert(ret is not None)  # do something more intelligent

    orig_ax = kwargs.pop('ax', plt.gca())
    if series is not None:
        assert(orig_ax.freq == series.index.freq)

    if freq is not None:
        assert(orig_ax.freq == freq)

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

