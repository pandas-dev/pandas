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
import pandas.tseries.plotting as plt

from pandas.util.testing import assert_series_equal, assert_almost_equal
import pandas.util.testing as tm

class TestTSPlot(unittest.TestCase):

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
        import matplotlib.pyplot as pyplot
        ax = pyplot.gca()
        for s in self.period_ser:
            _check_plot_works(plt.tsplot, s.index.freq, axes=ax, series=s)
        for s in self.datetime_ser:
            _check_plot_works(plt.tsplot, s.index.freq.rule_code,
                              axes=ax, series=s)

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
            ser.inferred_freq = None
            _check_plot_works(ser.plot, ser.index.inferred_freq)


PNG_PATH = 'tmp.png'
def _check_plot_works(f, freq, *args, **kwargs):
    import matplotlib.pyplot as plt

    fig = plt.gcf()
    plt.clf()
    ax = fig.add_subplot(211)
    ret = f(*args, **kwargs)
    assert(ret is not None)  # do something more intelligent

    orig_ax = kwargs.pop('axes', plt.gca())
    series = kwargs.pop('series', None)
    if series is not None:
        assert(orig_ax.freq == series.index.freq)

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

