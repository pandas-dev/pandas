"""Tests suite for Interval handling.

Parts derived from scikits.timeseries code, original authors:
- Pierre Gerard-Marchant & Matt Knox
- pierregm_at_uga_dot_edu - mattknow_ca_at_hotmail_dot_com

"""

from unittest import TestCase
from datetime import datetime

from numpy.ma.testutils import assert_equal
from pandas.core.datetools import Interval
from pandas.core.index import IntervalIndex, DatetimeIndex
import pandas.core.datetools as datetools
import numpy as np

from pandas import Series
from pandas.util.testing import assert_series_equal

class TestIntervalProperties(TestCase):
    "Test properties such as year, month, weekday, etc...."
    #
    def __init__(self, *args, **kwds):
        TestCase.__init__(self, *args, **kwds)

    def test_interval_constructor(self):
        i1 = Interval('1/1/2005', freq='M')
        i2 = Interval('Jan 2005')

        self.assertEquals(i1, i2)

        i1 = Interval('2005', freq='A')
        i2 = Interval('2005')
        i3 = Interval('2005', freq='a')

        self.assertEquals(i1, i2)
        self.assertEquals(i1, i3)

        i4 = Interval('2005', freq='M')
        i5 = Interval('2005', freq='m')

        self.assert_(i1 != i4)
        self.assertEquals(i4, i5)

        i1 = Interval.now('Q')
        i2 = Interval(datetime.now(), freq='Q')
        i3 = Interval.now('q')

        self.assertEquals(i1, i2)
        self.assertEquals(i1, i3)

        # Biz day construction, roll forward if non-weekday
        i1 = Interval('3/10/12', freq='B')
        i2 = Interval('3/12/12', freq='D')
        self.assertEquals(i1, i2.resample('B'))

        i3 = Interval('3/10/12', freq='b')
        self.assertEquals(i1, i3)

        i1 = Interval(year=2005, quarter=1, freq='Q')
        i2 = Interval('1/1/2005', freq='Q')
        self.assertEquals(i1, i2)

        i1 = Interval(year=2005, quarter=3, freq='Q')
        i2 = Interval('9/1/2005', freq='Q')
        self.assertEquals(i1, i2)

        i1 = Interval(year=2005, month=3, day=1, freq='D')
        i2 = Interval('3/1/2005', freq='D')
        self.assertEquals(i1, i2)

        i3 = Interval(year=2005, month=3, day=1, freq='d')
        self.assertEquals(i1, i3)

        i1 = Interval(year=2012, month=3, day=10, freq='B')
        i2 = Interval('3/12/12', freq='B')
        self.assertEquals(i1, i2)

        i1 = Interval('2005Q1')
        i2 = Interval(year=2005, quarter=1, freq='Q')
        i3 = Interval('2005q1')
        self.assertEquals(i1, i2)
        self.assertEquals(i1, i3)

        i1 = Interval('05Q1')
        self.assertEquals(i1, i2)
        lower = Interval('05q1')
        self.assertEquals(i1, lower)

        i1 = Interval('1Q2005')
        self.assertEquals(i1, i2)
        lower = Interval('1q2005')
        self.assertEquals(i1, lower)

        i1 = Interval('1Q05')
        self.assertEquals(i1, i2)
        lower = Interval('1q05')
        self.assertEquals(i1, lower)

        i1 = Interval('4Q1984')
        self.assertEquals(i1.year, 1984)
        lower = Interval('4q1984')
        self.assertEquals(i1, lower)

        i1 = Interval('1982', freq='min')
        i2 = Interval('1982', freq='MIN')
        self.assertEquals(i1, i2)
        i2 = Interval('1982', freq=('Min', 1))
        self.assertEquals(i1, i2)

    def test_freq_str(self):
        i1 = Interval('1982', freq='Min')
        self.assert_(i1.freq[0] != '1')

        i2 = Interval('11/30/2005', freq='2Q')
        self.assertEquals(i2.freq[0], '2')

    def test_to_timestamp(self):
        intv = Interval('1982', freq='A')
        start_ts = intv.to_timestamp(which_end='S')
        aliases = ['s', 'StarT', 'BEGIn']
        for a in aliases:
            self.assertEquals(start_ts, intv.to_timestamp(which_end=a))

        end_ts = intv.to_timestamp(which_end='E')
        aliases = ['e', 'end', 'FINIsH']
        for a in aliases:
            self.assertEquals(end_ts, intv.to_timestamp(which_end=a))

        from_lst = ['A', 'Q', 'M', 'W', 'B',
                    'D', 'H', 'Min', 'S']
        for i, fcode in enumerate(from_lst):
            intv = Interval('1982', freq=fcode)
            result = intv.to_timestamp().to_interval(fcode)
            self.assertEquals(result, intv)

            self.assertEquals(intv.start_time(), intv.to_timestamp('S'))

            self.assertEquals(intv.end_time(), intv.to_timestamp('E'))


    def test_properties_annually(self):
        # Test properties on Intervals with annually frequency.
        a_date = Interval(freq='A', year=2007)
        assert_equal(a_date.year, 2007)

    def test_properties_quarterly(self):
        # Test properties on Intervals with daily frequency.
        qedec_date = Interval(freq="Q@DEC", year=2007, quarter=1)
        qejan_date = Interval(freq="Q@JAN", year=2007, quarter=1)
        qejun_date = Interval(freq="Q@JUN", year=2007, quarter=1)
        #
        for x in range(3):
            for qd in (qedec_date, qejan_date, qejun_date):
                assert_equal((qd + x).qyear, 2007)
                assert_equal((qd + x).quarter, x + 1)


    def test_properties_monthly(self):
        # Test properties on Intervals with daily frequency.
        m_date = Interval(freq='M', year=2007, month=1)
        for x in range(11):
            m_ival_x = m_date + x
            assert_equal(m_ival_x.year, 2007)
            if 1 <= x + 1 <= 3:
                assert_equal(m_ival_x.quarter, 1)
            elif 4 <= x + 1 <= 6:
                assert_equal(m_ival_x.quarter, 2)
            elif 7 <= x + 1 <= 9:
                assert_equal(m_ival_x.quarter, 3)
            elif 10 <= x + 1 <= 12:
                assert_equal(m_ival_x.quarter, 4)
            assert_equal(m_ival_x.month, x + 1)


    def test_properties_weekly(self):
        # Test properties on Intervals with daily frequency.
        w_date = Interval(freq='WK', year=2007, month=1, day=7)
        #
        assert_equal(w_date.year, 2007)
        assert_equal(w_date.quarter, 1)
        assert_equal(w_date.month, 1)
        assert_equal(w_date.week, 1)
        assert_equal((w_date - 1).week, 52)


    def test_properties_daily(self):
        # Test properties on Intervals with daily frequency.
        b_date = Interval(freq='B', year=2007, month=1, day=1)
        #
        assert_equal(b_date.year, 2007)
        assert_equal(b_date.quarter, 1)
        assert_equal(b_date.month, 1)
        assert_equal(b_date.day, 1)
        assert_equal(b_date.weekday, 0)
        assert_equal(b_date.day_of_year, 1)
        #
        d_date = Interval(freq='D', year=2007, month=1, day=1)
        #
        assert_equal(d_date.year, 2007)
        assert_equal(d_date.quarter, 1)
        assert_equal(d_date.month, 1)
        assert_equal(d_date.day, 1)
        assert_equal(d_date.weekday, 0)
        assert_equal(d_date.day_of_year, 1)


    def test_properties_hourly(self):
        # Test properties on Intervals with hourly frequency.
        h_date = Interval(freq='H', year=2007, month=1, day=1, hour=0)
        #
        assert_equal(h_date.year, 2007)
        assert_equal(h_date.quarter, 1)
        assert_equal(h_date.month, 1)
        assert_equal(h_date.day, 1)
        assert_equal(h_date.weekday, 0)
        assert_equal(h_date.day_of_year, 1)
        assert_equal(h_date.hour, 0)
        #


    def test_properties_minutely(self):
        # Test properties on Intervals with minutely frequency.
        t_date = Interval(freq='Min', year=2007, month=1, day=1, hour=0,
                          minute=0)
        #
        assert_equal(t_date.quarter, 1)
        assert_equal(t_date.month, 1)
        assert_equal(t_date.day, 1)
        assert_equal(t_date.weekday, 0)
        assert_equal(t_date.day_of_year, 1)
        assert_equal(t_date.hour, 0)
        assert_equal(t_date.minute, 0)


    def test_properties_secondly(self):
        # Test properties on Intervals with secondly frequency.
        s_date = Interval(freq='Min', year=2007, month=1, day=1,
                                       hour=0, minute=0, second=0)
        #
        assert_equal(s_date.year, 2007)
        assert_equal(s_date.quarter, 1)
        assert_equal(s_date.month, 1)
        assert_equal(s_date.day, 1)
        assert_equal(s_date.weekday, 0)
        assert_equal(s_date.day_of_year, 1)
        assert_equal(s_date.hour, 0)
        assert_equal(s_date.minute, 0)
        assert_equal(s_date.second, 0)

def noWrap(item):
    return item

class TestFreqConversion(TestCase):
    "Test frequency conversion of date objects"

    def __init__(self, *args, **kwds):
        TestCase.__init__(self, *args, **kwds)

    def test_conv_annual(self):
        # frequency conversion tests: from Annual Frequency

        ival_A = Interval(freq='A', year=2007)

        ival_AJAN = Interval(freq="A@JAN", year=2007)
        ival_AJUN = Interval(freq="A@JUN", year=2007)
        ival_ANOV = Interval(freq="A@NOV", year=2007)

        ival_A_to_Q_start = Interval(freq='Q', year=2007, quarter=1)
        ival_A_to_Q_end = Interval(freq='Q', year=2007, quarter=4)
        ival_A_to_M_start = Interval(freq='M', year=2007, month=1)
        ival_A_to_M_end = Interval(freq='M', year=2007, month=12)
        ival_A_to_W_start = Interval(freq='WK', year=2007, month=1, day=1)
        ival_A_to_W_end = Interval(freq='WK', year=2007, month=12, day=31)
        ival_A_to_B_start = Interval(freq='B', year=2007, month=1, day=1)
        ival_A_to_B_end = Interval(freq='B', year=2007, month=12, day=31)
        ival_A_to_D_start = Interval(freq='D', year=2007, month=1, day=1)
        ival_A_to_D_end = Interval(freq='D', year=2007, month=12, day=31)
        ival_A_to_H_start = Interval(freq='H', year=2007, month=1, day=1,
                                    hour=0)
        ival_A_to_H_end = Interval(freq='H', year=2007, month=12, day=31,
                                    hour=23)
        ival_A_to_T_start = Interval(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=0)
        ival_A_to_T_end = Interval(freq='Min', year=2007, month=12, day=31,
                                    hour=23, minute=59)
        ival_A_to_S_start = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_A_to_S_end = Interval(freq='S', year=2007, month=12, day=31,
                                    hour=23, minute=59, second=59)

        ival_AJAN_to_D_end = Interval(freq='D', year=2007, month=1, day=31)
        ival_AJAN_to_D_start = Interval(freq='D', year=2006, month=2, day=1)
        ival_AJUN_to_D_end = Interval(freq='D', year=2007, month=6, day=30)
        ival_AJUN_to_D_start = Interval(freq='D', year=2006, month=7, day=1)
        ival_ANOV_to_D_end = Interval(freq='D', year=2007, month=11, day=30)
        ival_ANOV_to_D_start = Interval(freq='D', year=2006, month=12, day=1)

        assert_equal(ival_A.resample('Q', 'S'), ival_A_to_Q_start)
        assert_equal(ival_A.resample('Q', 'E'), ival_A_to_Q_end)
        assert_equal(ival_A.resample('M', 'S'), ival_A_to_M_start)
        assert_equal(ival_A.resample('M', 'E'), ival_A_to_M_end)
        assert_equal(ival_A.resample('WK', 'S'), ival_A_to_W_start)
        assert_equal(ival_A.resample('WK', 'E'), ival_A_to_W_end)
        assert_equal(ival_A.resample('B', 'S'), ival_A_to_B_start)
        assert_equal(ival_A.resample('B', 'E'), ival_A_to_B_end)
        assert_equal(ival_A.resample('D', 'S'), ival_A_to_D_start)
        assert_equal(ival_A.resample('D', 'E'), ival_A_to_D_end)
        assert_equal(ival_A.resample('H', 'S'), ival_A_to_H_start)
        assert_equal(ival_A.resample('H', 'E'), ival_A_to_H_end)
        assert_equal(ival_A.resample('Min', 'S'), ival_A_to_T_start)
        assert_equal(ival_A.resample('Min', 'E'), ival_A_to_T_end)
        assert_equal(ival_A.resample('S', 'S'), ival_A_to_S_start)
        assert_equal(ival_A.resample('S', 'E'), ival_A_to_S_end)

        assert_equal(ival_AJAN.resample('D', 'S'), ival_AJAN_to_D_start)
        assert_equal(ival_AJAN.resample('D', 'E'), ival_AJAN_to_D_end)

        assert_equal(ival_AJUN.resample('D', 'S'), ival_AJUN_to_D_start)
        assert_equal(ival_AJUN.resample('D', 'E'), ival_AJUN_to_D_end)

        assert_equal(ival_ANOV.resample('D', 'S'), ival_ANOV_to_D_start)
        assert_equal(ival_ANOV.resample('D', 'E'), ival_ANOV_to_D_end)

        assert_equal(ival_A.resample('A'), ival_A)


    def test_conv_quarterly(self):
        # frequency conversion tests: from Quarterly Frequency

        ival_Q = Interval(freq='Q', year=2007, quarter=1)
        ival_Q_end_of_year = Interval(freq='Q', year=2007, quarter=4)

        ival_QEJAN = Interval(freq="Q@JAN", year=2007, quarter=1)
        ival_QEJUN = Interval(freq="Q@JUN", year=2007, quarter=1)

        ival_Q_to_A = Interval(freq='A', year=2007)
        ival_Q_to_M_start = Interval(freq='M', year=2007, month=1)
        ival_Q_to_M_end = Interval(freq='M', year=2007, month=3)
        ival_Q_to_W_start = Interval(freq='WK', year=2007, month=1, day=1)
        ival_Q_to_W_end = Interval(freq='WK', year=2007, month=3, day=31)
        ival_Q_to_B_start = Interval(freq='B', year=2007, month=1, day=1)
        ival_Q_to_B_end = Interval(freq='B', year=2007, month=3, day=30)
        ival_Q_to_D_start = Interval(freq='D', year=2007, month=1, day=1)
        ival_Q_to_D_end = Interval(freq='D', year=2007, month=3, day=31)
        ival_Q_to_H_start = Interval(freq='H', year=2007, month=1, day=1,
                                    hour=0)
        ival_Q_to_H_end = Interval(freq='H', year=2007, month=3, day=31,
                                    hour=23)
        ival_Q_to_T_start = Interval(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=0)
        ival_Q_to_T_end = Interval(freq='Min', year=2007, month=3, day=31,
                                    hour=23, minute=59)
        ival_Q_to_S_start = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_Q_to_S_end = Interval(freq='S', year=2007, month=3, day=31,
                                    hour=23, minute=59, second=59)

        ival_QEJAN_to_D_start = Interval(freq='D', year=2006, month=2, day=1)
        ival_QEJAN_to_D_end = Interval(freq='D', year=2006, month=4, day=30)

        ival_QEJUN_to_D_start = Interval(freq='D', year=2006, month=7, day=1)
        ival_QEJUN_to_D_end = Interval(freq='D', year=2006, month=9, day=30)

        assert_equal(ival_Q.resample('A'), ival_Q_to_A)
        assert_equal(ival_Q_end_of_year.resample('A'), ival_Q_to_A)

        assert_equal(ival_Q.resample('M', 'S'), ival_Q_to_M_start)
        assert_equal(ival_Q.resample('M', 'E'), ival_Q_to_M_end)
        assert_equal(ival_Q.resample('WK', 'S'), ival_Q_to_W_start)
        assert_equal(ival_Q.resample('WK', 'E'), ival_Q_to_W_end)
        assert_equal(ival_Q.resample('B', 'S'), ival_Q_to_B_start)
        assert_equal(ival_Q.resample('B', 'E'), ival_Q_to_B_end)
        assert_equal(ival_Q.resample('D', 'S'), ival_Q_to_D_start)
        assert_equal(ival_Q.resample('D', 'E'), ival_Q_to_D_end)
        assert_equal(ival_Q.resample('H', 'S'), ival_Q_to_H_start)
        assert_equal(ival_Q.resample('H', 'E'), ival_Q_to_H_end)
        assert_equal(ival_Q.resample('Min', 'S'), ival_Q_to_T_start)
        assert_equal(ival_Q.resample('Min', 'E'), ival_Q_to_T_end)
        assert_equal(ival_Q.resample('S', 'S'), ival_Q_to_S_start)
        assert_equal(ival_Q.resample('S', 'E'), ival_Q_to_S_end)

        assert_equal(ival_QEJAN.resample('D', 'S'), ival_QEJAN_to_D_start)
        assert_equal(ival_QEJAN.resample('D', 'E'), ival_QEJAN_to_D_end)
        assert_equal(ival_QEJUN.resample('D', 'S'), ival_QEJUN_to_D_start)
        assert_equal(ival_QEJUN.resample('D', 'E'), ival_QEJUN_to_D_end)

        assert_equal(ival_Q.resample('Q'), ival_Q)


    def test_conv_monthly(self):
        # frequency conversion tests: from Monthly Frequency

        ival_M = Interval(freq='M', year=2007, month=1)
        ival_M_end_of_year = Interval(freq='M', year=2007, month=12)
        ival_M_end_of_quarter = Interval(freq='M', year=2007, month=3)
        ival_M_to_A = Interval(freq='A', year=2007)
        ival_M_to_Q = Interval(freq='Q', year=2007, quarter=1)
        ival_M_to_W_start = Interval(freq='WK', year=2007, month=1, day=1)
        ival_M_to_W_end = Interval(freq='WK', year=2007, month=1, day=31)
        ival_M_to_B_start = Interval(freq='B', year=2007, month=1, day=1)
        ival_M_to_B_end = Interval(freq='B', year=2007, month=1, day=31)
        ival_M_to_D_start = Interval(freq='D', year=2007, month=1, day=1)
        ival_M_to_D_end = Interval(freq='D', year=2007, month=1, day=31)
        ival_M_to_H_start = Interval(freq='H', year=2007, month=1, day=1,
                                    hour=0)
        ival_M_to_H_end = Interval(freq='H', year=2007, month=1, day=31,
                                    hour=23)
        ival_M_to_T_start = Interval(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=0)
        ival_M_to_T_end = Interval(freq='Min', year=2007, month=1, day=31,
                                    hour=23, minute=59)
        ival_M_to_S_start = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_M_to_S_end = Interval(freq='S', year=2007, month=1, day=31,
                                    hour=23, minute=59, second=59)

        assert_equal(ival_M.resample('A'), ival_M_to_A)
        assert_equal(ival_M_end_of_year.resample('A'), ival_M_to_A)
        assert_equal(ival_M.resample('Q'), ival_M_to_Q)
        assert_equal(ival_M_end_of_quarter.resample('Q'), ival_M_to_Q)

        assert_equal(ival_M.resample('WK', 'S'), ival_M_to_W_start)
        assert_equal(ival_M.resample('WK', 'E'), ival_M_to_W_end)
        assert_equal(ival_M.resample('B', 'S'), ival_M_to_B_start)
        assert_equal(ival_M.resample('B', 'E'), ival_M_to_B_end)
        assert_equal(ival_M.resample('D', 'S'), ival_M_to_D_start)
        assert_equal(ival_M.resample('D', 'E'), ival_M_to_D_end)
        assert_equal(ival_M.resample('H', 'S'), ival_M_to_H_start)
        assert_equal(ival_M.resample('H', 'E'), ival_M_to_H_end)
        assert_equal(ival_M.resample('Min', 'S'), ival_M_to_T_start)
        assert_equal(ival_M.resample('Min', 'E'), ival_M_to_T_end)
        assert_equal(ival_M.resample('S', 'S'), ival_M_to_S_start)
        assert_equal(ival_M.resample('S', 'E'), ival_M_to_S_end)

        assert_equal(ival_M.resample('M'), ival_M)


    def test_conv_weekly(self):
        # frequency conversion tests: from Weekly Frequency

        ival_W = Interval(freq='WK', year=2007, month=1, day=1)

        ival_WSUN = Interval(freq='WK', year=2007, month=1, day=7)
        ival_WSAT = Interval(freq='WK@SAT', year=2007, month=1, day=6)
        ival_WFRI = Interval(freq='WK@FRI', year=2007, month=1, day=5)
        ival_WTHU = Interval(freq='WK@THU', year=2007, month=1, day=4)
        ival_WWED = Interval(freq='WK@WED', year=2007, month=1, day=3)
        ival_WTUE = Interval(freq='WK@TUE', year=2007, month=1, day=2)
        ival_WMON = Interval(freq='WK@MON', year=2007, month=1, day=1)

        ival_WSUN_to_D_start = Interval(freq='D', year=2007, month=1, day=1)
        ival_WSUN_to_D_end = Interval(freq='D', year=2007, month=1, day=7)
        ival_WSAT_to_D_start = Interval(freq='D', year=2006, month=12, day=31)
        ival_WSAT_to_D_end = Interval(freq='D', year=2007, month=1, day=6)
        ival_WFRI_to_D_start = Interval(freq='D', year=2006, month=12, day=30)
        ival_WFRI_to_D_end = Interval(freq='D', year=2007, month=1, day=5)
        ival_WTHU_to_D_start = Interval(freq='D', year=2006, month=12, day=29)
        ival_WTHU_to_D_end = Interval(freq='D', year=2007, month=1, day=4)
        ival_WWED_to_D_start = Interval(freq='D', year=2006, month=12, day=28)
        ival_WWED_to_D_end = Interval(freq='D', year=2007, month=1, day=3)
        ival_WTUE_to_D_start = Interval(freq='D', year=2006, month=12, day=27)
        ival_WTUE_to_D_end = Interval(freq='D', year=2007, month=1, day=2)
        ival_WMON_to_D_start = Interval(freq='D', year=2006, month=12, day=26)
        ival_WMON_to_D_end = Interval(freq='D', year=2007, month=1, day=1)

        ival_W_end_of_year = Interval(freq='WK', year=2007, month=12, day=31)
        ival_W_end_of_quarter = Interval(freq='WK', year=2007, month=3, day=31)
        ival_W_end_of_month = Interval(freq='WK', year=2007, month=1, day=31)
        ival_W_to_A = Interval(freq='A', year=2007)
        ival_W_to_Q = Interval(freq='Q', year=2007, quarter=1)
        ival_W_to_M = Interval(freq='M', year=2007, month=1)

        if Interval(freq='D', year=2007, month=12, day=31).weekday == 6:
            ival_W_to_A_end_of_year = Interval(freq='A', year=2007)
        else:
            ival_W_to_A_end_of_year = Interval(freq='A', year=2008)

        if Interval(freq='D', year=2007, month=3, day=31).weekday == 6:
            ival_W_to_Q_end_of_quarter = Interval(freq='Q', year=2007,
                                                  quarter=1)
        else:
            ival_W_to_Q_end_of_quarter = Interval(freq='Q', year=2007,
                                                  quarter=2)

        if Interval(freq='D', year=2007, month=1, day=31).weekday == 6:
            ival_W_to_M_end_of_month = Interval(freq='M', year=2007, month=1)
        else:
            ival_W_to_M_end_of_month = Interval(freq='M', year=2007, month=2)

        ival_W_to_B_start = Interval(freq='B', year=2007, month=1, day=1)
        ival_W_to_B_end = Interval(freq='B', year=2007, month=1, day=5)
        ival_W_to_D_start = Interval(freq='D', year=2007, month=1, day=1)
        ival_W_to_D_end = Interval(freq='D', year=2007, month=1, day=7)
        ival_W_to_H_start = Interval(freq='H', year=2007, month=1, day=1,
                                    hour=0)
        ival_W_to_H_end = Interval(freq='H', year=2007, month=1, day=7,
                                    hour=23)
        ival_W_to_T_start = Interval(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=0)
        ival_W_to_T_end = Interval(freq='Min', year=2007, month=1, day=7,
                                    hour=23, minute=59)
        ival_W_to_S_start = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_W_to_S_end = Interval(freq='S', year=2007, month=1, day=7,
                                    hour=23, minute=59, second=59)

        assert_equal(ival_W.resample('A'), ival_W_to_A)
        assert_equal(ival_W_end_of_year.resample('A'),
                     ival_W_to_A_end_of_year)
        assert_equal(ival_W.resample('Q'), ival_W_to_Q)
        assert_equal(ival_W_end_of_quarter.resample('Q'),
                     ival_W_to_Q_end_of_quarter)
        assert_equal(ival_W.resample('M'), ival_W_to_M)
        assert_equal(ival_W_end_of_month.resample('M'),
                     ival_W_to_M_end_of_month)

        assert_equal(ival_W.resample('B', 'S'), ival_W_to_B_start)
        assert_equal(ival_W.resample('B', 'E'), ival_W_to_B_end)

        assert_equal(ival_W.resample('D', 'S'), ival_W_to_D_start)
        assert_equal(ival_W.resample('D', 'E'), ival_W_to_D_end)

        assert_equal(ival_WSUN.resample('D', 'S'), ival_WSUN_to_D_start)
        assert_equal(ival_WSUN.resample('D', 'E'), ival_WSUN_to_D_end)
        assert_equal(ival_WSAT.resample('D', 'S'), ival_WSAT_to_D_start)
        assert_equal(ival_WSAT.resample('D', 'E'), ival_WSAT_to_D_end)
        assert_equal(ival_WFRI.resample('D', 'S'), ival_WFRI_to_D_start)
        assert_equal(ival_WFRI.resample('D', 'E'), ival_WFRI_to_D_end)
        assert_equal(ival_WTHU.resample('D', 'S'), ival_WTHU_to_D_start)
        assert_equal(ival_WTHU.resample('D', 'E'), ival_WTHU_to_D_end)
        assert_equal(ival_WWED.resample('D', 'S'), ival_WWED_to_D_start)
        assert_equal(ival_WWED.resample('D', 'E'), ival_WWED_to_D_end)
        assert_equal(ival_WTUE.resample('D', 'S'), ival_WTUE_to_D_start)
        assert_equal(ival_WTUE.resample('D', 'E'), ival_WTUE_to_D_end)
        assert_equal(ival_WMON.resample('D', 'S'), ival_WMON_to_D_start)
        assert_equal(ival_WMON.resample('D', 'E'), ival_WMON_to_D_end)

        assert_equal(ival_W.resample('H', 'S'), ival_W_to_H_start)
        assert_equal(ival_W.resample('H', 'E'), ival_W_to_H_end)
        assert_equal(ival_W.resample('Min', 'S'), ival_W_to_T_start)
        assert_equal(ival_W.resample('Min', 'E'), ival_W_to_T_end)
        assert_equal(ival_W.resample('S', 'S'), ival_W_to_S_start)
        assert_equal(ival_W.resample('S', 'E'), ival_W_to_S_end)

        assert_equal(ival_W.resample('WK'), ival_W)


    def test_conv_business(self):
        # frequency conversion tests: from Business Frequency"

        ival_B = Interval(freq='B', year=2007, month=1, day=1)
        ival_B_end_of_year = Interval(freq='B', year=2007, month=12, day=31)
        ival_B_end_of_quarter = Interval(freq='B', year=2007, month=3, day=30)
        ival_B_end_of_month = Interval(freq='B', year=2007, month=1, day=31)
        ival_B_end_of_week = Interval(freq='B', year=2007, month=1, day=5)

        ival_B_to_A = Interval(freq='A', year=2007)
        ival_B_to_Q = Interval(freq='Q', year=2007, quarter=1)
        ival_B_to_M = Interval(freq='M', year=2007, month=1)
        ival_B_to_W = Interval(freq='WK', year=2007, month=1, day=7)
        ival_B_to_D = Interval(freq='D', year=2007, month=1, day=1)
        ival_B_to_H_start = Interval(freq='H', year=2007, month=1, day=1,
                                    hour=0)
        ival_B_to_H_end = Interval(freq='H', year=2007, month=1, day=1,
                                    hour=23)
        ival_B_to_T_start = Interval(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=0)
        ival_B_to_T_end = Interval(freq='Min', year=2007, month=1, day=1,
                                    hour=23, minute=59)
        ival_B_to_S_start = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_B_to_S_end = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=23, minute=59, second=59)

        assert_equal(ival_B.resample('A'), ival_B_to_A)
        assert_equal(ival_B_end_of_year.resample('A'), ival_B_to_A)
        assert_equal(ival_B.resample('Q'), ival_B_to_Q)
        assert_equal(ival_B_end_of_quarter.resample('Q'), ival_B_to_Q)
        assert_equal(ival_B.resample('M'), ival_B_to_M)
        assert_equal(ival_B_end_of_month.resample('M'), ival_B_to_M)
        assert_equal(ival_B.resample('WK'), ival_B_to_W)
        assert_equal(ival_B_end_of_week.resample('WK'), ival_B_to_W)

        assert_equal(ival_B.resample('D'), ival_B_to_D)

        assert_equal(ival_B.resample('H', 'S'), ival_B_to_H_start)
        assert_equal(ival_B.resample('H', 'E'), ival_B_to_H_end)
        assert_equal(ival_B.resample('Min', 'S'), ival_B_to_T_start)
        assert_equal(ival_B.resample('Min', 'E'), ival_B_to_T_end)
        assert_equal(ival_B.resample('S', 'S'), ival_B_to_S_start)
        assert_equal(ival_B.resample('S', 'E'), ival_B_to_S_end)

        assert_equal(ival_B.resample('B'), ival_B)


    def test_conv_daily(self):
        # frequency conversion tests: from Business Frequency"

        ival_D = Interval(freq='D', year=2007, month=1, day=1)
        ival_D_end_of_year = Interval(freq='D', year=2007, month=12, day=31)
        ival_D_end_of_quarter = Interval(freq='D', year=2007, month=3, day=31)
        ival_D_end_of_month = Interval(freq='D', year=2007, month=1, day=31)
        ival_D_end_of_week = Interval(freq='D', year=2007, month=1, day=7)

        ival_D_friday = Interval(freq='D', year=2007, month=1, day=5)
        ival_D_saturday = Interval(freq='D', year=2007, month=1, day=6)
        ival_D_sunday = Interval(freq='D', year=2007, month=1, day=7)
        ival_D_monday = Interval(freq='D', year=2007, month=1, day=8)

        ival_B_friday = Interval(freq='B', year=2007, month=1, day=5)
        ival_B_monday = Interval(freq='B', year=2007, month=1, day=8)

        ival_D_to_A = Interval(freq='A', year=2007)

        ival_Deoq_to_AJAN = Interval(freq='A@JAN', year=2008)
        ival_Deoq_to_AJUN = Interval(freq='A@JUN', year=2007)
        ival_Deoq_to_ADEC = Interval(freq='A@DEC', year=2007)

        ival_D_to_QEJAN = Interval(freq="Q@JAN", year=2007, quarter=4)
        ival_D_to_QEJUN = Interval(freq="Q@JUN", year=2007, quarter=3)
        ival_D_to_QEDEC = Interval(freq="Q@DEC", year=2007, quarter=1)

        ival_D_to_M = Interval(freq='M', year=2007, month=1)
        ival_D_to_W = Interval(freq='WK', year=2007, month=1, day=7)

        ival_D_to_H_start = Interval(freq='H', year=2007, month=1, day=1,
                                    hour=0)
        ival_D_to_H_end = Interval(freq='H', year=2007, month=1, day=1,
                                    hour=23)
        ival_D_to_T_start = Interval(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=0)
        ival_D_to_T_end = Interval(freq='Min', year=2007, month=1, day=1,
                                    hour=23, minute=59)
        ival_D_to_S_start = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_D_to_S_end = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=23, minute=59, second=59)

        assert_equal(ival_D.resample('A'), ival_D_to_A)

        assert_equal(ival_D_end_of_quarter.resample('A@JAN'),
                     ival_Deoq_to_AJAN)
        assert_equal(ival_D_end_of_quarter.resample('A@JUN'),
                     ival_Deoq_to_AJUN)
        assert_equal(ival_D_end_of_quarter.resample('A@DEC'),
                     ival_Deoq_to_ADEC)

        assert_equal(ival_D_end_of_year.resample('A'), ival_D_to_A)
        assert_equal(ival_D_end_of_quarter.resample('Q'), ival_D_to_QEDEC)
        assert_equal(ival_D.resample("Q@JAN"), ival_D_to_QEJAN)
        assert_equal(ival_D.resample("Q@JUN"), ival_D_to_QEJUN)
        assert_equal(ival_D.resample("Q@DEC"), ival_D_to_QEDEC)
        assert_equal(ival_D.resample('M'), ival_D_to_M)
        assert_equal(ival_D_end_of_month.resample('M'), ival_D_to_M)
        assert_equal(ival_D.resample('WK'), ival_D_to_W)
        assert_equal(ival_D_end_of_week.resample('WK'), ival_D_to_W)

        assert_equal(ival_D_friday.resample('B'), ival_B_friday)
        assert_equal(ival_D_saturday.resample('B', 'S'), ival_B_friday)
        assert_equal(ival_D_saturday.resample('B', 'E'), ival_B_monday)
        assert_equal(ival_D_sunday.resample('B', 'S'), ival_B_friday)
        assert_equal(ival_D_sunday.resample('B', 'E'), ival_B_monday)

        assert_equal(ival_D.resample('H', 'S'), ival_D_to_H_start)
        assert_equal(ival_D.resample('H', 'E'), ival_D_to_H_end)
        assert_equal(ival_D.resample('Min', 'S'), ival_D_to_T_start)
        assert_equal(ival_D.resample('Min', 'E'), ival_D_to_T_end)
        assert_equal(ival_D.resample('S', 'S'), ival_D_to_S_start)
        assert_equal(ival_D.resample('S', 'E'), ival_D_to_S_end)

        assert_equal(ival_D.resample('D'), ival_D)

    def test_conv_hourly(self):
        # frequency conversion tests: from Hourly Frequency"

        ival_H = Interval(freq='H', year=2007, month=1, day=1, hour=0)
        ival_H_end_of_year = Interval(freq='H', year=2007, month=12, day=31,
                                    hour=23)
        ival_H_end_of_quarter = Interval(freq='H', year=2007, month=3, day=31,
                                        hour=23)
        ival_H_end_of_month = Interval(freq='H', year=2007, month=1, day=31,
                                    hour=23)
        ival_H_end_of_week = Interval(freq='H', year=2007, month=1, day=7,
                                    hour=23)
        ival_H_end_of_day = Interval(freq='H', year=2007, month=1, day=1,
                                    hour=23)
        ival_H_end_of_bus = Interval(freq='H', year=2007, month=1, day=1,
                                    hour=23)

        ival_H_to_A = Interval(freq='A', year=2007)
        ival_H_to_Q = Interval(freq='Q', year=2007, quarter=1)
        ival_H_to_M = Interval(freq='M', year=2007, month=1)
        ival_H_to_W = Interval(freq='WK', year=2007, month=1, day=7)
        ival_H_to_D = Interval(freq='D', year=2007, month=1, day=1)
        ival_H_to_B = Interval(freq='B', year=2007, month=1, day=1)

        ival_H_to_T_start = Interval(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=0)
        ival_H_to_T_end = Interval(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=59)
        ival_H_to_S_start = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_H_to_S_end = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=59, second=59)

        assert_equal(ival_H.resample('A'), ival_H_to_A)
        assert_equal(ival_H_end_of_year.resample('A'), ival_H_to_A)
        assert_equal(ival_H.resample('Q'), ival_H_to_Q)
        assert_equal(ival_H_end_of_quarter.resample('Q'), ival_H_to_Q)
        assert_equal(ival_H.resample('M'), ival_H_to_M)
        assert_equal(ival_H_end_of_month.resample('M'), ival_H_to_M)
        assert_equal(ival_H.resample('WK'), ival_H_to_W)
        assert_equal(ival_H_end_of_week.resample('WK'), ival_H_to_W)
        assert_equal(ival_H.resample('D'), ival_H_to_D)
        assert_equal(ival_H_end_of_day.resample('D'), ival_H_to_D)
        assert_equal(ival_H.resample('B'), ival_H_to_B)
        assert_equal(ival_H_end_of_bus.resample('B'), ival_H_to_B)

        assert_equal(ival_H.resample('Min', 'S'), ival_H_to_T_start)
        assert_equal(ival_H.resample('Min', 'E'), ival_H_to_T_end)
        assert_equal(ival_H.resample('S', 'S'), ival_H_to_S_start)
        assert_equal(ival_H.resample('S', 'E'), ival_H_to_S_end)

        assert_equal(ival_H.resample('H'), ival_H)

    def test_conv_minutely(self):
        # frequency conversion tests: from Minutely Frequency"

        ival_T = Interval(freq='Min', year=2007, month=1, day=1,
                        hour=0, minute=0)
        ival_T_end_of_year = Interval(freq='Min', year=2007, month=12, day=31,
                                    hour=23, minute=59)
        ival_T_end_of_quarter = Interval(freq='Min', year=2007, month=3, day=31,
                                        hour=23, minute=59)
        ival_T_end_of_month = Interval(freq='Min', year=2007, month=1, day=31,
                                    hour=23, minute=59)
        ival_T_end_of_week = Interval(freq='Min', year=2007, month=1, day=7,
                                    hour=23, minute=59)
        ival_T_end_of_day = Interval(freq='Min', year=2007, month=1, day=1,
                                    hour=23, minute=59)
        ival_T_end_of_bus = Interval(freq='Min', year=2007, month=1, day=1,
                                    hour=23, minute=59)
        ival_T_end_of_hour = Interval(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=59)

        ival_T_to_A = Interval(freq='A', year=2007)
        ival_T_to_Q = Interval(freq='Q', year=2007, quarter=1)
        ival_T_to_M = Interval(freq='M', year=2007, month=1)
        ival_T_to_W = Interval(freq='WK', year=2007, month=1, day=7)
        ival_T_to_D = Interval(freq='D', year=2007, month=1, day=1)
        ival_T_to_B = Interval(freq='B', year=2007, month=1, day=1)
        ival_T_to_H = Interval(freq='H', year=2007, month=1, day=1, hour=0)

        ival_T_to_S_start = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_T_to_S_end = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=59)

        assert_equal(ival_T.resample('A'), ival_T_to_A)
        assert_equal(ival_T_end_of_year.resample('A'), ival_T_to_A)
        assert_equal(ival_T.resample('Q'), ival_T_to_Q)
        assert_equal(ival_T_end_of_quarter.resample('Q'), ival_T_to_Q)
        assert_equal(ival_T.resample('M'), ival_T_to_M)
        assert_equal(ival_T_end_of_month.resample('M'), ival_T_to_M)
        assert_equal(ival_T.resample('WK'), ival_T_to_W)
        assert_equal(ival_T_end_of_week.resample('WK'), ival_T_to_W)
        assert_equal(ival_T.resample('D'), ival_T_to_D)
        assert_equal(ival_T_end_of_day.resample('D'), ival_T_to_D)
        assert_equal(ival_T.resample('B'), ival_T_to_B)
        assert_equal(ival_T_end_of_bus.resample('B'), ival_T_to_B)
        assert_equal(ival_T.resample('H'), ival_T_to_H)
        assert_equal(ival_T_end_of_hour.resample('H'), ival_T_to_H)

        assert_equal(ival_T.resample('S', 'S'), ival_T_to_S_start)
        assert_equal(ival_T.resample('S', 'E'), ival_T_to_S_end)

        assert_equal(ival_T.resample('Min'), ival_T)

    def test_conv_secondly(self):
        # frequency conversion tests: from Secondly Frequency"

        ival_S = Interval(freq='S', year=2007, month=1, day=1,
                        hour=0, minute=0, second=0)
        ival_S_end_of_year = Interval(freq='S', year=2007, month=12, day=31,
                                    hour=23, minute=59, second=59)
        ival_S_end_of_quarter = Interval(freq='S', year=2007, month=3, day=31,
                                        hour=23, minute=59, second=59)
        ival_S_end_of_month = Interval(freq='S', year=2007, month=1, day=31,
                                    hour=23, minute=59, second=59)
        ival_S_end_of_week = Interval(freq='S', year=2007, month=1, day=7,
                                    hour=23, minute=59, second=59)
        ival_S_end_of_day = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=23, minute=59, second=59)
        ival_S_end_of_bus = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=23, minute=59, second=59)
        ival_S_end_of_hour = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=59, second=59)
        ival_S_end_of_minute = Interval(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=59)

        ival_S_to_A = Interval(freq='A', year=2007)
        ival_S_to_Q = Interval(freq='Q', year=2007, quarter=1)
        ival_S_to_M = Interval(freq='M', year=2007, month=1)
        ival_S_to_W = Interval(freq='WK', year=2007, month=1, day=7)
        ival_S_to_D = Interval(freq='D', year=2007, month=1, day=1)
        ival_S_to_B = Interval(freq='B', year=2007, month=1, day=1)
        ival_S_to_H = Interval(freq='H', year=2007, month=1, day=1,
                            hour=0)
        ival_S_to_T = Interval(freq='Min', year=2007, month=1, day=1,
                            hour=0, minute=0)

        assert_equal(ival_S.resample('A'), ival_S_to_A)
        assert_equal(ival_S_end_of_year.resample('A'), ival_S_to_A)
        assert_equal(ival_S.resample('Q'), ival_S_to_Q)
        assert_equal(ival_S_end_of_quarter.resample('Q'), ival_S_to_Q)
        assert_equal(ival_S.resample('M'), ival_S_to_M)
        assert_equal(ival_S_end_of_month.resample('M'), ival_S_to_M)
        assert_equal(ival_S.resample('WK'), ival_S_to_W)
        assert_equal(ival_S_end_of_week.resample('WK'), ival_S_to_W)
        assert_equal(ival_S.resample('D'), ival_S_to_D)
        assert_equal(ival_S_end_of_day.resample('D'), ival_S_to_D)
        assert_equal(ival_S.resample('B'), ival_S_to_B)
        assert_equal(ival_S_end_of_bus.resample('B'), ival_S_to_B)
        assert_equal(ival_S.resample('H'), ival_S_to_H)
        assert_equal(ival_S_end_of_hour.resample('H'), ival_S_to_H)
        assert_equal(ival_S.resample('Min'), ival_S_to_T)
        assert_equal(ival_S_end_of_minute.resample('Min'), ival_S_to_T)

        assert_equal(ival_S.resample('S'), ival_S)

class TestIntervalIndex(TestCase):
    def __init__(self, *args, **kwds):
        TestCase.__init__(self, *args, **kwds)

    def test_constructor(self):
        ii = IntervalIndex(freq='A', start='1/1/2001', end='12/1/2009')
        assert_equal(len(ii), 9)

        ii = IntervalIndex(freq='Q', start='1/1/2001', end='12/1/2009')
        assert_equal(len(ii), 4 * 9)

        ii = IntervalIndex(freq='M', start='1/1/2001', end='12/1/2009')
        assert_equal(len(ii), 12 * 9)

        ii = IntervalIndex(freq='D', start='1/1/2001', end='12/31/2009')
        assert_equal(len(ii), 365 * 9 + 2)

        ii = IntervalIndex(freq='B', start='1/1/2001', end='12/31/2009')
        assert_equal(len(ii), 261 * 9)

        ii = IntervalIndex(freq='H', start='1/1/2001', end='12/31/2001 23:00')
        assert_equal(len(ii), 365 * 24)

        ii = IntervalIndex(freq='Min', start='1/1/2001', end='1/1/2001 23:59')
        assert_equal(len(ii), 24 * 60)

        ii = IntervalIndex(freq='S', start='1/1/2001', end='1/1/2001 23:59:59')
        assert_equal(len(ii), 24 * 60 * 60)

        start = Interval('02-Apr-2005', 'B')
        i1 = IntervalIndex(start=start, periods=20)
        assert_equal(len(i1), 20)
        assert_equal(i1.freq, start.freq)
        assert_equal(i1[0], start)

        end_intv = Interval('2006-12-31', 'W')
        i1 = IntervalIndex(end=end_intv, periods=10)
        assert_equal(len(i1), 10)
        assert_equal(i1.freq, end_intv.freq)
        assert_equal(i1[-1], end_intv)

        end_intv = Interval('2006-12-31', '1w')
        i2 = IntervalIndex(end=end_intv, periods=10)
        assert_equal(len(i1), len(i2))
        self.assert_((i1 == i2).all())
        assert_equal(i1.freq, i2.freq)

        end_intv = Interval('2006-12-31', ('w', 1))
        i2 = IntervalIndex(end=end_intv, periods=10)
        assert_equal(len(i1), len(i2))
        self.assert_((i1 == i2).all())
        assert_equal(i1.freq, i2.freq)

        try:
            IntervalIndex(start=start, end=end_intv)
            raise AssertionError('Cannot allow mixed freq for start and end')
        except ValueError:
            pass

        end_intv = Interval('2005-05-01', 'B')
        i1 = IntervalIndex(start=start, end=end_intv)

        try:
            IntervalIndex(start=start)
            raise AssertionError('Must specify periods if missing start or end')
        except ValueError:
            pass

    def test_shift(self):
        ii1 = IntervalIndex(freq='A', start='1/1/2001', end='12/1/2009')
        ii2 = IntervalIndex(freq='A', start='1/1/2002', end='12/1/2010')
        assert_equal(len(ii1), len(ii2))
        assert_equal(ii1.shift(1).values, ii2.values)

        ii1 = IntervalIndex(freq='A', start='1/1/2001', end='12/1/2009')
        ii2 = IntervalIndex(freq='A', start='1/1/2000', end='12/1/2008')
        assert_equal(len(ii1), len(ii2))
        assert_equal(ii1.shift(-1).values, ii2.values)

        ii1 = IntervalIndex(freq='M', start='1/1/2001', end='12/1/2009')
        ii2 = IntervalIndex(freq='M', start='2/1/2001', end='1/1/2010')
        assert_equal(len(ii1), len(ii2))
        assert_equal(ii1.shift(1).values, ii2.values)

        ii1 = IntervalIndex(freq='M', start='1/1/2001', end='12/1/2009')
        ii2 = IntervalIndex(freq='M', start='12/1/2000', end='11/1/2009')
        assert_equal(len(ii1), len(ii2))
        assert_equal(ii1.shift(-1).values, ii2.values)

        ii1 = IntervalIndex(freq='D', start='1/1/2001', end='12/1/2009')
        ii2 = IntervalIndex(freq='D', start='1/2/2001', end='12/2/2009')
        assert_equal(len(ii1), len(ii2))
        assert_equal(ii1.shift(1).values, ii2.values)

        ii1 = IntervalIndex(freq='D', start='1/1/2001', end='12/1/2009')
        ii2 = IntervalIndex(freq='D', start='12/31/2000', end='11/30/2009')
        assert_equal(len(ii1), len(ii2))
        assert_equal(ii1.shift(-1).values, ii2.values)

    def test_resample(self):
        ii1 = IntervalIndex(freq='A', start='1/1/2001', end='1/1/2001')
        ii2 = IntervalIndex(freq='Q', start='1/1/2001', end='1/1/2001')
        ii3 = IntervalIndex(freq='M', start='1/1/2001', end='1/1/2001')
        ii4 = IntervalIndex(freq='D', start='1/1/2001', end='1/1/2001')
        ii5 = IntervalIndex(freq='H', start='1/1/2001', end='1/1/2001 00:00')
        ii6 = IntervalIndex(freq='Min', start='1/1/2001', end='1/1/2001 00:00')
        ii7 = IntervalIndex(freq='S', start='1/1/2001', end='1/1/2001 00:00:00')

        self.assertEquals(ii1.resample('Q', 'S'), ii2)
        self.assertEquals(ii1.resample('Q', 's'), ii2)
        self.assertEquals(ii1.resample('M', 'start'), ii3)
        self.assertEquals(ii1.resample('D', 'StarT'), ii4)
        self.assertEquals(ii1.resample('H', 'beGIN'), ii5)
        self.assertEquals(ii1.resample('Min', 'S'), ii6)
        self.assertEquals(ii1.resample('S', 'S'), ii7)

        self.assertEquals(ii2.resample('A', 'S'), ii1)
        self.assertEquals(ii2.resample('M', 'S'), ii3)
        self.assertEquals(ii2.resample('D', 'S'), ii4)
        self.assertEquals(ii2.resample('H', 'S'), ii5)
        self.assertEquals(ii2.resample('Min', 'S'), ii6)
        self.assertEquals(ii2.resample('S', 'S'), ii7)

        self.assertEquals(ii3.resample('A', 'S'), ii1)
        self.assertEquals(ii3.resample('Q', 'S'), ii2)
        self.assertEquals(ii3.resample('D', 'S'), ii4)
        self.assertEquals(ii3.resample('H', 'S'), ii5)
        self.assertEquals(ii3.resample('Min', 'S'), ii6)
        self.assertEquals(ii3.resample('S', 'S'), ii7)

        self.assertEquals(ii4.resample('A', 'S'), ii1)
        self.assertEquals(ii4.resample('Q', 'S'), ii2)
        self.assertEquals(ii4.resample('M', 'S'), ii3)
        self.assertEquals(ii4.resample('H', 'S'), ii5)
        self.assertEquals(ii4.resample('Min', 'S'), ii6)
        self.assertEquals(ii4.resample('S', 'S'), ii7)

        self.assertEquals(ii5.resample('A', 'S'), ii1)
        self.assertEquals(ii5.resample('Q', 'S'), ii2)
        self.assertEquals(ii5.resample('M', 'S'), ii3)
        self.assertEquals(ii5.resample('D', 'S'), ii4)
        self.assertEquals(ii5.resample('Min', 'S'), ii6)
        self.assertEquals(ii5.resample('S', 'S'), ii7)

        self.assertEquals(ii6.resample('A', 'S'), ii1)
        self.assertEquals(ii6.resample('Q', 'S'), ii2)
        self.assertEquals(ii6.resample('M', 'S'), ii3)
        self.assertEquals(ii6.resample('D', 'S'), ii4)
        self.assertEquals(ii6.resample('H', 'S'), ii5)
        self.assertEquals(ii6.resample('S', 'S'), ii7)

        self.assertEquals(ii7.resample('A', 'S'), ii1)
        self.assertEquals(ii7.resample('Q', 'S'), ii2)
        self.assertEquals(ii7.resample('M', 'S'), ii3)
        self.assertEquals(ii7.resample('D', 'S'), ii4)
        self.assertEquals(ii7.resample('H', 'S'), ii5)
        self.assertEquals(ii7.resample('Min', 'S'), ii6)

        #self.assertEquals(ii7.resample('A', 'E'), i_end)

    def test_badinput(self):
        self.assertRaises(datetools.DateParseError, Interval, '1/1/-2000', 'A')
        self.assertRaises(ValueError, Interval, -2000, 'A')
        self.assertRaises(ValueError, Interval, 0, 'A')
        self.assertRaises(ValueError, IntervalIndex, [-1, 0, 1], 'A')
        self.assertRaises(ValueError, IntervalIndex, np.array([-1, 0, 1]), 'A')

    def test_dti_to_interval(self):
        dti = DatetimeIndex(start='1/1/2005', end='12/1/2005', freq='M')
        ii1 = dti.to_interval()
        ii2 = dti.to_interval(freq='D')

        self.assertEquals(ii1[0], Interval('Jan 2005', freq='M'))
        self.assertEquals(ii2[0], Interval('1/31/2005', freq='D'))

        self.assertEquals(ii1[-1], Interval('Nov 2005', freq='M'))
        self.assertEquals(ii2[-1], Interval('11/30/2005', freq='D'))

    def test_iindex_slice_index(self):
        ii = IntervalIndex(start='1/1/10', end='12/31/12', freq='M')
        s = Series(np.random.rand(len(ii)), index=ii)
        res = s['2010']
        exp = s[0:12]
        assert_series_equal(res, exp)
        res = s['2011']
        exp = s[12:24]
        assert_series_equal(res, exp)

    def test_iindex_qaccess(self):
        ii = IntervalIndex(['2Q05', '3Q05', '4Q05', '1Q06', '2Q06'], freq='Q')
        s = Series(np.random.rand(len(ii)), index=ii).cumsum()
        # Todo: fix these accessors!
        self.assert_(s['05Q4'] == s[2])

    def test_interval_dt64_round_trip(self):
        dti = DatetimeIndex(['1/1/2002', '1/2/2002', '1/3/2002', '1/4/2002',
                             '1/5/2002', '1/6/2002', '1/7/2002'], freq='B')
        ii = dti.to_interval()
        self.assert_(ii.to_timestamp().equals(dti))

        dti = DatetimeIndex(['1/1/2002', '1/2/2002', '1/3/2002', '1/4/2002',
                             '1/5/2002', '1/6/2002', '1/7/2002'], freq='B')
        ii = dti.to_interval(freq='3H')
        self.assert_(ii.to_timestamp().equals(dti))

    def test_iindex_multiples(self):
        ii = IntervalIndex(start='1/1/10', end='12/31/12', freq='2M')
        self.assertEquals(ii[0], Interval('1/1/10', '2M'))
        self.assertEquals(ii[1], Interval('3/1/10', '2M'))

        self.assertEquals(ii[0].resample('6M'), ii[2].resample('6M'))
        self.assertEquals(ii[0].resample('A'), ii[2].resample('A'))

        self.assertEquals(ii[0].resample('M', how='S'),
                          Interval('Jan 2010', '1M'))
        self.assertEquals(ii[0].resample('M', how='E'),
                          Interval('Feb 2010', '1M'))
        self.assertEquals(ii[1].resample('M', how='S'),
                          Interval('Mar 2010', '1M'))

        i = Interval('1/1/2010 12:05:18', '5S')
        self.assertEquals(i, Interval('1/1/2010 12:05:15', '5S'))

        i = Interval('1/1/2010 12:05:18', '5S')
        self.assertEquals(i.resample('1S', how='E'),
                          Interval('1/1/2010 12:05:19', '1S'))

class TestMethods(TestCase):
    "Base test class for MaskedArrays."

    def __init__(self, *args, **kwds):
        TestCase.__init__(self, *args, **kwds)

    def test_add(self):
        dt1 = Interval(freq='D', year=2008, month=1, day=1)
        dt2 = Interval(freq='D', year=2008, month=1, day=2)
        assert_equal(dt1 + 1, dt2)
        #
        self.assertRaises(ValueError, dt1.__add__, "str")
        self.assertRaises(ValueError, dt1.__add__, dt2)

###############################################################################
#------------------------------------------------------------------------------

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
