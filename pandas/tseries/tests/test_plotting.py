import os
from datetime import datetime, timedelta

import unittest
import nose

import numpy as np
from numpy.testing.decorators import slow

from pandas import Index, Series, DataFrame, isnull, notnull

from pandas.tseries.index import date_range
from pandas.tseries.offsets import Minute, DateOffset
from pandas.tseries.period import period_range, Period
from pandas.tseries.resample import DatetimeIndex, TimeGrouper
import pandas.tseries.offsets as offsets
import pandas.tseries.frequencies as frequencies

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
    def test_plot_offset_freq(self):
        ser = tm.makeTimeSeries()
        _check_plot_works(ser.plot)

        dr = date_range(ser.index[0], freq='BQS', periods=10)
        ser = Series(np.random.randn(len(dr)), dr)
        _check_plot_works(ser.plot)

    @slow
    def test_plot_multiple_inferred_freq(self):
        dr = Index([datetime(2000, 1, 1),
                    datetime(2000, 1, 6),
                    datetime(2000, 1, 11)])
        ser = Series(np.random.randn(len(dr)), dr)
        _check_plot_works(ser.plot)

    @slow
    def test_irregular_datetime64_repr_bug(self):
        ser = tm.makeTimeSeries()
        ser = ser[[0,1,2,7]]
        import matplotlib.pyplot as plt

        fig = plt.gcf()
        plt.clf()
        ax = fig.add_subplot(211)
        ret = ser.plot()
        assert(ret is not None)

        for rs, xp in zip(ax.get_lines()[0].get_xdata(), ser.index):
            assert(rs == xp)

    @slow
    def test_business_freq(self):
        bts = tm.makePeriodSeries()
        ts = bts.asfreq('D')
        ax = bts.plot()
        self.assert_(ax.get_lines()[0].get_xydata()[0, 0], ts.index[0].ordinal)
        idx = ax.get_lines()[0].get_xdata()
        self.assert_(idx.freqstr == 'D')

    @slow
    def test_dataframe(self):
        bts = DataFrame({'a': tm.makePeriodSeries()})
        ts = bts.asfreq('D')
        ax = bts.plot()
        self.assert_(ax.get_lines()[0].get_xydata()[0, 0], ts.index[0].ordinal)
        idx = ax.get_lines()[0].get_xdata()
        self.assert_(idx.freqstr == 'D')

    @slow
    def test_set_xlim(self):
        ser = tm.makeTimeSeries()
        ax = ser.plot()
        xlim = ax.get_xlim()
        ax.set_xlim(xlim[0] - 5, xlim[1] + 10)
        ax.get_figure().canvas.draw()
        result = ax.get_xlim()
        self.assertEqual(result[0], xlim[0] - 5)
        self.assertEqual(result[1], xlim[1] + 10)

        # string
        expected = (Period('1/1/2000', ax.freq), Period('4/1/2000', ax.freq))
        ax.set_xlim('1/1/2000', '4/1/2000')
        ax.get_figure().canvas.draw()
        result = ax.get_xlim()
        self.assertEqual(int(result[0]), expected[0].ordinal)
        self.assertEqual(int(result[1]), expected[1].ordinal)

        # datetim
        expected = (Period('1/1/2000', ax.freq), Period('4/1/2000', ax.freq))
        ax.set_xlim(datetime(2000, 1, 1), datetime(2000, 4, 1))
        ax.get_figure().canvas.draw()
        result = ax.get_xlim()
        self.assertEqual(int(result[0]), expected[0].ordinal)
        self.assertEqual(int(result[1]), expected[1].ordinal)

    @slow
    def test_finder_quarterly(self):
        xp = Period('1988Q1').ordinal
        yrs = [3.5, 11]
        for n in yrs:
            rng = period_range('1987Q2', periods=int(n * 4), freq='Q')
            ser = Series(np.random.randn(len(rng)), rng)
            ax = ser.plot()
            xaxis = ax.get_xaxis()
            rs = xaxis.get_majorticklocs()[0]
            self.assert_(rs == xp)

    @slow
    def test_finder_monthly(self):
        xp = Period('1988-1').ordinal
        yrs = [1.15, 2.5, 4, 11]
        for n in yrs:
            rng = period_range('1987Q2', periods=int(n * 12), freq='M')
            ser = Series(np.random.randn(len(rng)), rng)
            ax = ser.plot()
            xaxis = ax.get_xaxis()
            rs = xaxis.get_majorticklocs()[0]
            self.assert_(rs == xp)

    @slow
    def test_gaps(self):
        import matplotlib.pyplot as plt
        plt.close('all')
        ts = tm.makeTimeSeries()
        ts[5:25] = np.nan
        ax = ts.plot()
        lines = ax.get_lines()
        self.assert_(len(lines) == 1)
        l = lines[0]
        data = l.get_xydata()
        self.assert_(isinstance(data, np.ma.core.MaskedArray))
        mask = data.mask
        self.assert_(mask[5:25, 1].all())

    @slow
    def test_secondary_y(self):
        import matplotlib.pyplot as plt
        plt.close('all')
        ser = Series(np.random.randn(10))
        ser2 = Series(np.random.randn(10))
        ax = ser.plot(secondary_y=True)
        fig = ax.get_figure()
        axes = fig.get_axes()
        l = ax.get_lines()[0]
        xp = Series(l.get_ydata(), l.get_xdata())
        assert_series_equal(ser, xp)
        self.assert_(ax.get_yaxis().get_ticks_position() == 'right')
        self.assert_(not axes[0].get_yaxis().get_visible())

        ax2 = ser2.plot()
        self.assert_(ax2.get_yaxis().get_ticks_position() == 'left')

        plt.close('all')
        ax = ser2.plot()
        ax2 = ser.plot(secondary_y=True)
        self.assert_(ax.get_yaxis().get_visible())

    @slow
    def test_secondary_y_ts(self):
        import matplotlib.pyplot as plt
        plt.close('all')
        idx = date_range('1/1/2000', periods=10)
        ser = Series(np.random.randn(10), idx)
        ser2 = Series(np.random.randn(10), idx)
        ax = ser.plot(secondary_y=True)
        fig = ax.get_figure()
        axes = fig.get_axes()
        l = ax.get_lines()[0]
        xp = Series(l.get_ydata(), l.get_xdata()).to_timestamp()
        assert_series_equal(ser, xp)
        self.assert_(ax.get_yaxis().get_ticks_position() == 'right')
        self.assert_(not axes[0].get_yaxis().get_visible())

        ax2 = ser2.plot()
        self.assert_(ax2.get_yaxis().get_ticks_position() == 'left')

        plt.close('all')
        ax = ser2.plot()
        ax2 = ser.plot(secondary_y=True)
        self.assert_(ax.get_yaxis().get_visible())

    @slow
    def test_secondary_kde(self):
        import matplotlib.pyplot as plt
        plt.close('all')
        ser = Series(np.random.randn(10))
        ax = ser.plot(secondary_y=True, kind='density')
        fig = ax.get_figure()
        axes = fig.get_axes()
        self.assert_(axes[1].get_yaxis().get_ticks_position() == 'right')

    @slow
    def test_secondary_bar(self):
        import matplotlib.pyplot as plt
        plt.close('all')
        ser = Series(np.random.randn(10))
        ax = ser.plot(secondary_y=True, kind='bar')
        fig = ax.get_figure()
        axes = fig.get_axes()
        self.assert_(axes[1].get_yaxis().get_ticks_position() == 'right')

    @slow
    def test_secondary_frame(self):
        import matplotlib.pyplot as plt
        plt.close('all')
        df = DataFrame(np.random.randn(5, 3), columns=['a', 'b', 'c'])
        axes = df.plot(secondary_y=['a', 'c'], subplots=True)
        self.assert_(axes[0].get_yaxis().get_ticks_position() == 'right')
        self.assert_(axes[1].get_yaxis().get_ticks_position() == 'default')
        self.assert_(axes[2].get_yaxis().get_ticks_position() == 'right')

PNG_PATH = 'tmp.png'
def _check_plot_works(f, freq=None, series=None, *args, **kwargs):
    import matplotlib.pyplot as plt

    fig = plt.gcf()
    plt.clf()
    ax = fig.add_subplot(211)
    ret = f(*args, **kwargs)
    assert(ret is not None)  # do something more intelligent

    orig_ax = kwargs.pop('ax', plt.gca())
    if series is not None: # non-business
        dfreq = series.index.freq
        if isinstance(dfreq, DateOffset):
            dfreq = dfreq.rule_code
        #dfreq = frequencies.offset_to_period_alias(dfreq)
        assert(orig_ax.freq == dfreq)

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

