"""Tests suite for Period handling.

Parts derived from scikits.timeseries code, original authors:
- Pierre Gerard-Marchant & Matt Knox
- pierregm_at_uga_dot_edu - mattknow_ca_at_hotmail_dot_com

"""

from unittest import TestCase
from datetime import datetime, timedelta

from numpy.ma.testutils import assert_equal

from pandas.tseries.frequencies import MONTHS, DAYS
from pandas.tseries.period import Period, PeriodIndex, period_range
from pandas.tseries.index import DatetimeIndex, date_range
from pandas.tseries.tools import to_datetime

import pandas.core.datetools as datetools
import numpy as np

from pandas import Series, TimeSeries, DataFrame
from pandas.util.testing import assert_series_equal
import pandas.util.testing as tm

class TestPeriodProperties(TestCase):
    "Test properties such as year, month, weekday, etc...."
    #
    def __init__(self, *args, **kwds):
        TestCase.__init__(self, *args, **kwds)

    def test_period_cons_quarterly(self):
        # bugs in scikits.timeseries
        for month in MONTHS:
            freq = 'Q-%s' % month
            exp = Period('1989Q3', freq=freq)
            stamp = exp.to_timestamp('D', how='end')
            p = Period(stamp, freq=freq)
            self.assertEquals(p, exp)

    def test_period_cons_annual(self):
        # bugs in scikits.timeseries
        for month in MONTHS:
            freq = 'A-%s' % month
            exp = Period('1989', freq=freq)
            stamp = exp.to_timestamp('D', how='end') + 30
            p = Period(stamp, freq=freq)
            self.assertEquals(p, exp + 1)

    def test_period_constructor(self):
        i1 = Period('1/1/2005', freq='M')
        i2 = Period('Jan 2005')

        self.assertEquals(i1, i2)

        i1 = Period('2005', freq='A')
        i2 = Period('2005')
        i3 = Period('2005', freq='a')

        self.assertEquals(i1, i2)
        self.assertEquals(i1, i3)

        i4 = Period('2005', freq='M')
        i5 = Period('2005', freq='m')

        self.assert_(i1 != i4)
        self.assertEquals(i4, i5)

        i1 = Period.now('Q')
        i2 = Period(datetime.now(), freq='Q')
        i3 = Period.now('q')

        self.assertEquals(i1, i2)
        self.assertEquals(i1, i3)

        # Biz day construction, roll forward if non-weekday
        i1 = Period('3/10/12', freq='B')
        i2 = Period('3/12/12', freq='D')
        self.assertEquals(i1, i2.asfreq('B'))

        i3 = Period('3/10/12', freq='b')
        self.assertEquals(i1, i3)

        i1 = Period(year=2005, quarter=1, freq='Q')
        i2 = Period('1/1/2005', freq='Q')
        self.assertEquals(i1, i2)

        i1 = Period(year=2005, quarter=3, freq='Q')
        i2 = Period('9/1/2005', freq='Q')
        self.assertEquals(i1, i2)

        i1 = Period(year=2005, month=3, day=1, freq='D')
        i2 = Period('3/1/2005', freq='D')
        self.assertEquals(i1, i2)

        i3 = Period(year=2005, month=3, day=1, freq='d')
        self.assertEquals(i1, i3)

        i1 = Period(year=2012, month=3, day=10, freq='B')
        i2 = Period('3/12/12', freq='B')
        self.assertEquals(i1, i2)

        i1 = Period('2005Q1')
        i2 = Period(year=2005, quarter=1, freq='Q')
        i3 = Period('2005q1')
        self.assertEquals(i1, i2)
        self.assertEquals(i1, i3)

        i1 = Period('05Q1')
        self.assertEquals(i1, i2)
        lower = Period('05q1')
        self.assertEquals(i1, lower)

        i1 = Period('1Q2005')
        self.assertEquals(i1, i2)
        lower = Period('1q2005')
        self.assertEquals(i1, lower)

        i1 = Period('1Q05')
        self.assertEquals(i1, i2)
        lower = Period('1q05')
        self.assertEquals(i1, lower)

        i1 = Period('4Q1984')
        self.assertEquals(i1.year, 1984)
        lower = Period('4q1984')
        self.assertEquals(i1, lower)

        i1 = Period('1982', freq='min')
        i2 = Period('1982', freq='MIN')
        self.assertEquals(i1, i2)
        i2 = Period('1982', freq=('Min', 1))
        self.assertEquals(i1, i2)

        expected = Period('2007-01', freq='M')
        i1 = Period('200701', freq='M')
        self.assertEqual(i1, expected)

        i1 = Period('200701', freq='M')
        self.assertEqual(i1, expected)

        i1 = Period(200701, freq='M')
        self.assertEqual(i1, expected)

        i1 = Period(ordinal=200701, freq='M')
        self.assertEqual(i1.year, 16726)

        self.assertRaises(ValueError, Period, ordinal=200701)

    def test_freq_str(self):
        i1 = Period('1982', freq='Min')
        self.assert_(i1.freq[0] != '1')

        i2 = Period('11/30/2005', freq='2Q')
        self.assertEquals(i2.freq[0], '2')

    def test_to_timestamp(self):
        p = Period('1982', freq='A')
        start_ts = p.to_timestamp(how='S')
        aliases = ['s', 'StarT', 'BEGIn']
        for a in aliases:
            self.assertEquals(start_ts, p.to_timestamp(how=a))

        end_ts = p.to_timestamp(how='E')
        aliases = ['e', 'end', 'FINIsH']
        for a in aliases:
            self.assertEquals(end_ts, p.to_timestamp(how=a))

        from_lst = ['A', 'Q', 'M', 'W', 'B',
                    'D', 'H', 'Min', 'S']
        for i, fcode in enumerate(from_lst):
            p = Period('1982', freq=fcode)
            result = p.to_timestamp().to_period(fcode)
            self.assertEquals(result, p)

            self.assertEquals(p.start_time, p.to_timestamp(how='S'))

            self.assertEquals(p.end_time, p.to_timestamp(how='E'))

        # Frequency other than daily

        p = Period('1985', freq='A')

        result = p.to_timestamp('H', how='end')
        expected = datetime(1985, 12, 31, 23)
        self.assertEquals(result, expected)

        result = p.to_timestamp('T', how='end')
        expected = datetime(1985, 12, 31, 23, 59)
        self.assertEquals(result, expected)

        result = p.to_timestamp('S', how='end')
        expected = datetime(1985, 12, 31, 23, 59, 59)
        self.assertEquals(result, expected)

        expected = datetime(1985, 1, 1)
        result = p.to_timestamp('H', how='start')
        self.assertEquals(result, expected)
        result = p.to_timestamp('T', how='start')
        self.assertEquals(result, expected)
        result = p.to_timestamp('S', how='start')
        self.assertEquals(result, expected)

    def test_properties_annually(self):
        # Test properties on Periods with annually frequency.
        a_date = Period(freq='A', year=2007)
        assert_equal(a_date.year, 2007)

    def test_properties_quarterly(self):
        # Test properties on Periods with daily frequency.
        qedec_date = Period(freq="Q-DEC", year=2007, quarter=1)
        qejan_date = Period(freq="Q-JAN", year=2007, quarter=1)
        qejun_date = Period(freq="Q-JUN", year=2007, quarter=1)
        #
        for x in range(3):
            for qd in (qedec_date, qejan_date, qejun_date):
                assert_equal((qd + x).qyear, 2007)
                assert_equal((qd + x).quarter, x + 1)


    def test_properties_monthly(self):
        # Test properties on Periods with daily frequency.
        m_date = Period(freq='M', year=2007, month=1)
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
        # Test properties on Periods with daily frequency.
        w_date = Period(freq='WK', year=2007, month=1, day=7)
        #
        assert_equal(w_date.year, 2007)
        assert_equal(w_date.quarter, 1)
        assert_equal(w_date.month, 1)
        assert_equal(w_date.week, 1)
        assert_equal((w_date - 1).week, 52)


    def test_properties_daily(self):
        # Test properties on Periods with daily frequency.
        b_date = Period(freq='B', year=2007, month=1, day=1)
        #
        assert_equal(b_date.year, 2007)
        assert_equal(b_date.quarter, 1)
        assert_equal(b_date.month, 1)
        assert_equal(b_date.day, 1)
        assert_equal(b_date.weekday, 0)
        assert_equal(b_date.day_of_year, 1)
        #
        d_date = Period(freq='D', year=2007, month=1, day=1)
        #
        assert_equal(d_date.year, 2007)
        assert_equal(d_date.quarter, 1)
        assert_equal(d_date.month, 1)
        assert_equal(d_date.day, 1)
        assert_equal(d_date.weekday, 0)
        assert_equal(d_date.day_of_year, 1)


    def test_properties_hourly(self):
        # Test properties on Periods with hourly frequency.
        h_date = Period(freq='H', year=2007, month=1, day=1, hour=0)
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
        # Test properties on Periods with minutely frequency.
        t_date = Period(freq='Min', year=2007, month=1, day=1, hour=0,
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
        # Test properties on Periods with secondly frequency.
        s_date = Period(freq='Min', year=2007, month=1, day=1,
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

        ival_A = Period(freq='A', year=2007)

        ival_AJAN = Period(freq="A-JAN", year=2007)
        ival_AJUN = Period(freq="A-JUN", year=2007)
        ival_ANOV = Period(freq="A-NOV", year=2007)

        ival_A_to_Q_start = Period(freq='Q', year=2007, quarter=1)
        ival_A_to_Q_end = Period(freq='Q', year=2007, quarter=4)
        ival_A_to_M_start = Period(freq='M', year=2007, month=1)
        ival_A_to_M_end = Period(freq='M', year=2007, month=12)
        ival_A_to_W_start = Period(freq='WK', year=2007, month=1, day=1)
        ival_A_to_W_end = Period(freq='WK', year=2007, month=12, day=31)
        ival_A_to_B_start = Period(freq='B', year=2007, month=1, day=1)
        ival_A_to_B_end = Period(freq='B', year=2007, month=12, day=31)
        ival_A_to_D_start = Period(freq='D', year=2007, month=1, day=1)
        ival_A_to_D_end = Period(freq='D', year=2007, month=12, day=31)
        ival_A_to_H_start = Period(freq='H', year=2007, month=1, day=1,
                                    hour=0)
        ival_A_to_H_end = Period(freq='H', year=2007, month=12, day=31,
                                    hour=23)
        ival_A_to_T_start = Period(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=0)
        ival_A_to_T_end = Period(freq='Min', year=2007, month=12, day=31,
                                    hour=23, minute=59)
        ival_A_to_S_start = Period(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_A_to_S_end = Period(freq='S', year=2007, month=12, day=31,
                                    hour=23, minute=59, second=59)

        ival_AJAN_to_D_end = Period(freq='D', year=2007, month=1, day=31)
        ival_AJAN_to_D_start = Period(freq='D', year=2006, month=2, day=1)
        ival_AJUN_to_D_end = Period(freq='D', year=2007, month=6, day=30)
        ival_AJUN_to_D_start = Period(freq='D', year=2006, month=7, day=1)
        ival_ANOV_to_D_end = Period(freq='D', year=2007, month=11, day=30)
        ival_ANOV_to_D_start = Period(freq='D', year=2006, month=12, day=1)

        assert_equal(ival_A.asfreq('Q', 'S'), ival_A_to_Q_start)
        assert_equal(ival_A.asfreq('Q', 'e'), ival_A_to_Q_end)
        assert_equal(ival_A.asfreq('M', 's'), ival_A_to_M_start)
        assert_equal(ival_A.asfreq('M', 'E'), ival_A_to_M_end)
        assert_equal(ival_A.asfreq('WK', 'S'), ival_A_to_W_start)
        assert_equal(ival_A.asfreq('WK', 'E'), ival_A_to_W_end)
        assert_equal(ival_A.asfreq('B', 'S'), ival_A_to_B_start)
        assert_equal(ival_A.asfreq('B', 'E'), ival_A_to_B_end)
        assert_equal(ival_A.asfreq('D', 'S'), ival_A_to_D_start)
        assert_equal(ival_A.asfreq('D', 'E'), ival_A_to_D_end)
        assert_equal(ival_A.asfreq('H', 'S'), ival_A_to_H_start)
        assert_equal(ival_A.asfreq('H', 'E'), ival_A_to_H_end)
        assert_equal(ival_A.asfreq('min', 'S'), ival_A_to_T_start)
        assert_equal(ival_A.asfreq('min', 'E'), ival_A_to_T_end)
        assert_equal(ival_A.asfreq('T', 'S'), ival_A_to_T_start)
        assert_equal(ival_A.asfreq('T', 'E'), ival_A_to_T_end)
        assert_equal(ival_A.asfreq('S', 'S'), ival_A_to_S_start)
        assert_equal(ival_A.asfreq('S', 'E'), ival_A_to_S_end)

        assert_equal(ival_AJAN.asfreq('D', 'S'), ival_AJAN_to_D_start)
        assert_equal(ival_AJAN.asfreq('D', 'E'), ival_AJAN_to_D_end)

        assert_equal(ival_AJUN.asfreq('D', 'S'), ival_AJUN_to_D_start)
        assert_equal(ival_AJUN.asfreq('D', 'E'), ival_AJUN_to_D_end)

        assert_equal(ival_ANOV.asfreq('D', 'S'), ival_ANOV_to_D_start)
        assert_equal(ival_ANOV.asfreq('D', 'E'), ival_ANOV_to_D_end)

        assert_equal(ival_A.asfreq('A'), ival_A)


    def test_conv_quarterly(self):
        # frequency conversion tests: from Quarterly Frequency

        ival_Q = Period(freq='Q', year=2007, quarter=1)
        ival_Q_end_of_year = Period(freq='Q', year=2007, quarter=4)

        ival_QEJAN = Period(freq="Q-JAN", year=2007, quarter=1)
        ival_QEJUN = Period(freq="Q-JUN", year=2007, quarter=1)

        ival_Q_to_A = Period(freq='A', year=2007)
        ival_Q_to_M_start = Period(freq='M', year=2007, month=1)
        ival_Q_to_M_end = Period(freq='M', year=2007, month=3)
        ival_Q_to_W_start = Period(freq='WK', year=2007, month=1, day=1)
        ival_Q_to_W_end = Period(freq='WK', year=2007, month=3, day=31)
        ival_Q_to_B_start = Period(freq='B', year=2007, month=1, day=1)
        ival_Q_to_B_end = Period(freq='B', year=2007, month=3, day=30)
        ival_Q_to_D_start = Period(freq='D', year=2007, month=1, day=1)
        ival_Q_to_D_end = Period(freq='D', year=2007, month=3, day=31)
        ival_Q_to_H_start = Period(freq='H', year=2007, month=1, day=1,
                                    hour=0)
        ival_Q_to_H_end = Period(freq='H', year=2007, month=3, day=31,
                                    hour=23)
        ival_Q_to_T_start = Period(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=0)
        ival_Q_to_T_end = Period(freq='Min', year=2007, month=3, day=31,
                                    hour=23, minute=59)
        ival_Q_to_S_start = Period(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_Q_to_S_end = Period(freq='S', year=2007, month=3, day=31,
                                    hour=23, minute=59, second=59)

        ival_QEJAN_to_D_start = Period(freq='D', year=2006, month=2, day=1)
        ival_QEJAN_to_D_end = Period(freq='D', year=2006, month=4, day=30)

        ival_QEJUN_to_D_start = Period(freq='D', year=2006, month=7, day=1)
        ival_QEJUN_to_D_end = Period(freq='D', year=2006, month=9, day=30)

        assert_equal(ival_Q.asfreq('A'), ival_Q_to_A)
        assert_equal(ival_Q_end_of_year.asfreq('A'), ival_Q_to_A)

        assert_equal(ival_Q.asfreq('M', 'S'), ival_Q_to_M_start)
        assert_equal(ival_Q.asfreq('M', 'E'), ival_Q_to_M_end)
        assert_equal(ival_Q.asfreq('WK', 'S'), ival_Q_to_W_start)
        assert_equal(ival_Q.asfreq('WK', 'E'), ival_Q_to_W_end)
        assert_equal(ival_Q.asfreq('B', 'S'), ival_Q_to_B_start)
        assert_equal(ival_Q.asfreq('B', 'E'), ival_Q_to_B_end)
        assert_equal(ival_Q.asfreq('D', 'S'), ival_Q_to_D_start)
        assert_equal(ival_Q.asfreq('D', 'E'), ival_Q_to_D_end)
        assert_equal(ival_Q.asfreq('H', 'S'), ival_Q_to_H_start)
        assert_equal(ival_Q.asfreq('H', 'E'), ival_Q_to_H_end)
        assert_equal(ival_Q.asfreq('Min', 'S'), ival_Q_to_T_start)
        assert_equal(ival_Q.asfreq('Min', 'E'), ival_Q_to_T_end)
        assert_equal(ival_Q.asfreq('S', 'S'), ival_Q_to_S_start)
        assert_equal(ival_Q.asfreq('S', 'E'), ival_Q_to_S_end)

        assert_equal(ival_QEJAN.asfreq('D', 'S'), ival_QEJAN_to_D_start)
        assert_equal(ival_QEJAN.asfreq('D', 'E'), ival_QEJAN_to_D_end)
        assert_equal(ival_QEJUN.asfreq('D', 'S'), ival_QEJUN_to_D_start)
        assert_equal(ival_QEJUN.asfreq('D', 'E'), ival_QEJUN_to_D_end)

        assert_equal(ival_Q.asfreq('Q'), ival_Q)


    def test_conv_monthly(self):
        # frequency conversion tests: from Monthly Frequency

        ival_M = Period(freq='M', year=2007, month=1)
        ival_M_end_of_year = Period(freq='M', year=2007, month=12)
        ival_M_end_of_quarter = Period(freq='M', year=2007, month=3)
        ival_M_to_A = Period(freq='A', year=2007)
        ival_M_to_Q = Period(freq='Q', year=2007, quarter=1)
        ival_M_to_W_start = Period(freq='WK', year=2007, month=1, day=1)
        ival_M_to_W_end = Period(freq='WK', year=2007, month=1, day=31)
        ival_M_to_B_start = Period(freq='B', year=2007, month=1, day=1)
        ival_M_to_B_end = Period(freq='B', year=2007, month=1, day=31)
        ival_M_to_D_start = Period(freq='D', year=2007, month=1, day=1)
        ival_M_to_D_end = Period(freq='D', year=2007, month=1, day=31)
        ival_M_to_H_start = Period(freq='H', year=2007, month=1, day=1,
                                    hour=0)
        ival_M_to_H_end = Period(freq='H', year=2007, month=1, day=31,
                                    hour=23)
        ival_M_to_T_start = Period(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=0)
        ival_M_to_T_end = Period(freq='Min', year=2007, month=1, day=31,
                                    hour=23, minute=59)
        ival_M_to_S_start = Period(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_M_to_S_end = Period(freq='S', year=2007, month=1, day=31,
                                    hour=23, minute=59, second=59)

        assert_equal(ival_M.asfreq('A'), ival_M_to_A)
        assert_equal(ival_M_end_of_year.asfreq('A'), ival_M_to_A)
        assert_equal(ival_M.asfreq('Q'), ival_M_to_Q)
        assert_equal(ival_M_end_of_quarter.asfreq('Q'), ival_M_to_Q)

        assert_equal(ival_M.asfreq('WK', 'S'), ival_M_to_W_start)
        assert_equal(ival_M.asfreq('WK', 'E'), ival_M_to_W_end)
        assert_equal(ival_M.asfreq('B', 'S'), ival_M_to_B_start)
        assert_equal(ival_M.asfreq('B', 'E'), ival_M_to_B_end)
        assert_equal(ival_M.asfreq('D', 'S'), ival_M_to_D_start)
        assert_equal(ival_M.asfreq('D', 'E'), ival_M_to_D_end)
        assert_equal(ival_M.asfreq('H', 'S'), ival_M_to_H_start)
        assert_equal(ival_M.asfreq('H', 'E'), ival_M_to_H_end)
        assert_equal(ival_M.asfreq('Min', 'S'), ival_M_to_T_start)
        assert_equal(ival_M.asfreq('Min', 'E'), ival_M_to_T_end)
        assert_equal(ival_M.asfreq('S', 'S'), ival_M_to_S_start)
        assert_equal(ival_M.asfreq('S', 'E'), ival_M_to_S_end)

        assert_equal(ival_M.asfreq('M'), ival_M)


    def test_conv_weekly(self):
        # frequency conversion tests: from Weekly Frequency

        ival_W = Period(freq='WK', year=2007, month=1, day=1)

        ival_WSUN = Period(freq='WK', year=2007, month=1, day=7)
        ival_WSAT = Period(freq='WK-SAT', year=2007, month=1, day=6)
        ival_WFRI = Period(freq='WK-FRI', year=2007, month=1, day=5)
        ival_WTHU = Period(freq='WK-THU', year=2007, month=1, day=4)
        ival_WWED = Period(freq='WK-WED', year=2007, month=1, day=3)
        ival_WTUE = Period(freq='WK-TUE', year=2007, month=1, day=2)
        ival_WMON = Period(freq='WK-MON', year=2007, month=1, day=1)

        ival_WSUN_to_D_start = Period(freq='D', year=2007, month=1, day=1)
        ival_WSUN_to_D_end = Period(freq='D', year=2007, month=1, day=7)
        ival_WSAT_to_D_start = Period(freq='D', year=2006, month=12, day=31)
        ival_WSAT_to_D_end = Period(freq='D', year=2007, month=1, day=6)
        ival_WFRI_to_D_start = Period(freq='D', year=2006, month=12, day=30)
        ival_WFRI_to_D_end = Period(freq='D', year=2007, month=1, day=5)
        ival_WTHU_to_D_start = Period(freq='D', year=2006, month=12, day=29)
        ival_WTHU_to_D_end = Period(freq='D', year=2007, month=1, day=4)
        ival_WWED_to_D_start = Period(freq='D', year=2006, month=12, day=28)
        ival_WWED_to_D_end = Period(freq='D', year=2007, month=1, day=3)
        ival_WTUE_to_D_start = Period(freq='D', year=2006, month=12, day=27)
        ival_WTUE_to_D_end = Period(freq='D', year=2007, month=1, day=2)
        ival_WMON_to_D_start = Period(freq='D', year=2006, month=12, day=26)
        ival_WMON_to_D_end = Period(freq='D', year=2007, month=1, day=1)

        ival_W_end_of_year = Period(freq='WK', year=2007, month=12, day=31)
        ival_W_end_of_quarter = Period(freq='WK', year=2007, month=3, day=31)
        ival_W_end_of_month = Period(freq='WK', year=2007, month=1, day=31)
        ival_W_to_A = Period(freq='A', year=2007)
        ival_W_to_Q = Period(freq='Q', year=2007, quarter=1)
        ival_W_to_M = Period(freq='M', year=2007, month=1)

        if Period(freq='D', year=2007, month=12, day=31).weekday == 6:
            ival_W_to_A_end_of_year = Period(freq='A', year=2007)
        else:
            ival_W_to_A_end_of_year = Period(freq='A', year=2008)

        if Period(freq='D', year=2007, month=3, day=31).weekday == 6:
            ival_W_to_Q_end_of_quarter = Period(freq='Q', year=2007,
                                                  quarter=1)
        else:
            ival_W_to_Q_end_of_quarter = Period(freq='Q', year=2007,
                                                  quarter=2)

        if Period(freq='D', year=2007, month=1, day=31).weekday == 6:
            ival_W_to_M_end_of_month = Period(freq='M', year=2007, month=1)
        else:
            ival_W_to_M_end_of_month = Period(freq='M', year=2007, month=2)

        ival_W_to_B_start = Period(freq='B', year=2007, month=1, day=1)
        ival_W_to_B_end = Period(freq='B', year=2007, month=1, day=5)
        ival_W_to_D_start = Period(freq='D', year=2007, month=1, day=1)
        ival_W_to_D_end = Period(freq='D', year=2007, month=1, day=7)
        ival_W_to_H_start = Period(freq='H', year=2007, month=1, day=1,
                                    hour=0)
        ival_W_to_H_end = Period(freq='H', year=2007, month=1, day=7,
                                    hour=23)
        ival_W_to_T_start = Period(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=0)
        ival_W_to_T_end = Period(freq='Min', year=2007, month=1, day=7,
                                    hour=23, minute=59)
        ival_W_to_S_start = Period(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_W_to_S_end = Period(freq='S', year=2007, month=1, day=7,
                                    hour=23, minute=59, second=59)

        assert_equal(ival_W.asfreq('A'), ival_W_to_A)
        assert_equal(ival_W_end_of_year.asfreq('A'),
                     ival_W_to_A_end_of_year)
        assert_equal(ival_W.asfreq('Q'), ival_W_to_Q)
        assert_equal(ival_W_end_of_quarter.asfreq('Q'),
                     ival_W_to_Q_end_of_quarter)
        assert_equal(ival_W.asfreq('M'), ival_W_to_M)
        assert_equal(ival_W_end_of_month.asfreq('M'),
                     ival_W_to_M_end_of_month)

        assert_equal(ival_W.asfreq('B', 'S'), ival_W_to_B_start)
        assert_equal(ival_W.asfreq('B', 'E'), ival_W_to_B_end)

        assert_equal(ival_W.asfreq('D', 'S'), ival_W_to_D_start)
        assert_equal(ival_W.asfreq('D', 'E'), ival_W_to_D_end)

        assert_equal(ival_WSUN.asfreq('D', 'S'), ival_WSUN_to_D_start)
        assert_equal(ival_WSUN.asfreq('D', 'E'), ival_WSUN_to_D_end)
        assert_equal(ival_WSAT.asfreq('D', 'S'), ival_WSAT_to_D_start)
        assert_equal(ival_WSAT.asfreq('D', 'E'), ival_WSAT_to_D_end)
        assert_equal(ival_WFRI.asfreq('D', 'S'), ival_WFRI_to_D_start)
        assert_equal(ival_WFRI.asfreq('D', 'E'), ival_WFRI_to_D_end)
        assert_equal(ival_WTHU.asfreq('D', 'S'), ival_WTHU_to_D_start)
        assert_equal(ival_WTHU.asfreq('D', 'E'), ival_WTHU_to_D_end)
        assert_equal(ival_WWED.asfreq('D', 'S'), ival_WWED_to_D_start)
        assert_equal(ival_WWED.asfreq('D', 'E'), ival_WWED_to_D_end)
        assert_equal(ival_WTUE.asfreq('D', 'S'), ival_WTUE_to_D_start)
        assert_equal(ival_WTUE.asfreq('D', 'E'), ival_WTUE_to_D_end)
        assert_equal(ival_WMON.asfreq('D', 'S'), ival_WMON_to_D_start)
        assert_equal(ival_WMON.asfreq('D', 'E'), ival_WMON_to_D_end)

        assert_equal(ival_W.asfreq('H', 'S'), ival_W_to_H_start)
        assert_equal(ival_W.asfreq('H', 'E'), ival_W_to_H_end)
        assert_equal(ival_W.asfreq('Min', 'S'), ival_W_to_T_start)
        assert_equal(ival_W.asfreq('Min', 'E'), ival_W_to_T_end)
        assert_equal(ival_W.asfreq('S', 'S'), ival_W_to_S_start)
        assert_equal(ival_W.asfreq('S', 'E'), ival_W_to_S_end)

        assert_equal(ival_W.asfreq('WK'), ival_W)


    def test_conv_business(self):
        # frequency conversion tests: from Business Frequency"

        ival_B = Period(freq='B', year=2007, month=1, day=1)
        ival_B_end_of_year = Period(freq='B', year=2007, month=12, day=31)
        ival_B_end_of_quarter = Period(freq='B', year=2007, month=3, day=30)
        ival_B_end_of_month = Period(freq='B', year=2007, month=1, day=31)
        ival_B_end_of_week = Period(freq='B', year=2007, month=1, day=5)

        ival_B_to_A = Period(freq='A', year=2007)
        ival_B_to_Q = Period(freq='Q', year=2007, quarter=1)
        ival_B_to_M = Period(freq='M', year=2007, month=1)
        ival_B_to_W = Period(freq='WK', year=2007, month=1, day=7)
        ival_B_to_D = Period(freq='D', year=2007, month=1, day=1)
        ival_B_to_H_start = Period(freq='H', year=2007, month=1, day=1,
                                    hour=0)
        ival_B_to_H_end = Period(freq='H', year=2007, month=1, day=1,
                                    hour=23)
        ival_B_to_T_start = Period(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=0)
        ival_B_to_T_end = Period(freq='Min', year=2007, month=1, day=1,
                                    hour=23, minute=59)
        ival_B_to_S_start = Period(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_B_to_S_end = Period(freq='S', year=2007, month=1, day=1,
                                    hour=23, minute=59, second=59)

        assert_equal(ival_B.asfreq('A'), ival_B_to_A)
        assert_equal(ival_B_end_of_year.asfreq('A'), ival_B_to_A)
        assert_equal(ival_B.asfreq('Q'), ival_B_to_Q)
        assert_equal(ival_B_end_of_quarter.asfreq('Q'), ival_B_to_Q)
        assert_equal(ival_B.asfreq('M'), ival_B_to_M)
        assert_equal(ival_B_end_of_month.asfreq('M'), ival_B_to_M)
        assert_equal(ival_B.asfreq('WK'), ival_B_to_W)
        assert_equal(ival_B_end_of_week.asfreq('WK'), ival_B_to_W)

        assert_equal(ival_B.asfreq('D'), ival_B_to_D)

        assert_equal(ival_B.asfreq('H', 'S'), ival_B_to_H_start)
        assert_equal(ival_B.asfreq('H', 'E'), ival_B_to_H_end)
        assert_equal(ival_B.asfreq('Min', 'S'), ival_B_to_T_start)
        assert_equal(ival_B.asfreq('Min', 'E'), ival_B_to_T_end)
        assert_equal(ival_B.asfreq('S', 'S'), ival_B_to_S_start)
        assert_equal(ival_B.asfreq('S', 'E'), ival_B_to_S_end)

        assert_equal(ival_B.asfreq('B'), ival_B)


    def test_conv_daily(self):
        # frequency conversion tests: from Business Frequency"

        ival_D = Period(freq='D', year=2007, month=1, day=1)
        ival_D_end_of_year = Period(freq='D', year=2007, month=12, day=31)
        ival_D_end_of_quarter = Period(freq='D', year=2007, month=3, day=31)
        ival_D_end_of_month = Period(freq='D', year=2007, month=1, day=31)
        ival_D_end_of_week = Period(freq='D', year=2007, month=1, day=7)

        ival_D_friday = Period(freq='D', year=2007, month=1, day=5)
        ival_D_saturday = Period(freq='D', year=2007, month=1, day=6)
        ival_D_sunday = Period(freq='D', year=2007, month=1, day=7)
        ival_D_monday = Period(freq='D', year=2007, month=1, day=8)

        ival_B_friday = Period(freq='B', year=2007, month=1, day=5)
        ival_B_monday = Period(freq='B', year=2007, month=1, day=8)

        ival_D_to_A = Period(freq='A', year=2007)

        ival_Deoq_to_AJAN = Period(freq='A-JAN', year=2008)
        ival_Deoq_to_AJUN = Period(freq='A-JUN', year=2007)
        ival_Deoq_to_ADEC = Period(freq='A-DEC', year=2007)

        ival_D_to_QEJAN = Period(freq="Q-JAN", year=2007, quarter=4)
        ival_D_to_QEJUN = Period(freq="Q-JUN", year=2007, quarter=3)
        ival_D_to_QEDEC = Period(freq="Q-DEC", year=2007, quarter=1)

        ival_D_to_M = Period(freq='M', year=2007, month=1)
        ival_D_to_W = Period(freq='WK', year=2007, month=1, day=7)

        ival_D_to_H_start = Period(freq='H', year=2007, month=1, day=1,
                                    hour=0)
        ival_D_to_H_end = Period(freq='H', year=2007, month=1, day=1,
                                    hour=23)
        ival_D_to_T_start = Period(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=0)
        ival_D_to_T_end = Period(freq='Min', year=2007, month=1, day=1,
                                    hour=23, minute=59)
        ival_D_to_S_start = Period(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_D_to_S_end = Period(freq='S', year=2007, month=1, day=1,
                                    hour=23, minute=59, second=59)

        assert_equal(ival_D.asfreq('A'), ival_D_to_A)

        assert_equal(ival_D_end_of_quarter.asfreq('A-JAN'),
                     ival_Deoq_to_AJAN)
        assert_equal(ival_D_end_of_quarter.asfreq('A-JUN'),
                     ival_Deoq_to_AJUN)
        assert_equal(ival_D_end_of_quarter.asfreq('A-DEC'),
                     ival_Deoq_to_ADEC)

        assert_equal(ival_D_end_of_year.asfreq('A'), ival_D_to_A)
        assert_equal(ival_D_end_of_quarter.asfreq('Q'), ival_D_to_QEDEC)
        assert_equal(ival_D.asfreq("Q-JAN"), ival_D_to_QEJAN)
        assert_equal(ival_D.asfreq("Q-JUN"), ival_D_to_QEJUN)
        assert_equal(ival_D.asfreq("Q-DEC"), ival_D_to_QEDEC)
        assert_equal(ival_D.asfreq('M'), ival_D_to_M)
        assert_equal(ival_D_end_of_month.asfreq('M'), ival_D_to_M)
        assert_equal(ival_D.asfreq('WK'), ival_D_to_W)
        assert_equal(ival_D_end_of_week.asfreq('WK'), ival_D_to_W)

        assert_equal(ival_D_friday.asfreq('B'), ival_B_friday)
        assert_equal(ival_D_saturday.asfreq('B', 'S'), ival_B_friday)
        assert_equal(ival_D_saturday.asfreq('B', 'E'), ival_B_monday)
        assert_equal(ival_D_sunday.asfreq('B', 'S'), ival_B_friday)
        assert_equal(ival_D_sunday.asfreq('B', 'E'), ival_B_monday)

        assert_equal(ival_D.asfreq('H', 'S'), ival_D_to_H_start)
        assert_equal(ival_D.asfreq('H', 'E'), ival_D_to_H_end)
        assert_equal(ival_D.asfreq('Min', 'S'), ival_D_to_T_start)
        assert_equal(ival_D.asfreq('Min', 'E'), ival_D_to_T_end)
        assert_equal(ival_D.asfreq('S', 'S'), ival_D_to_S_start)
        assert_equal(ival_D.asfreq('S', 'E'), ival_D_to_S_end)

        assert_equal(ival_D.asfreq('D'), ival_D)

    def test_conv_hourly(self):
        # frequency conversion tests: from Hourly Frequency"

        ival_H = Period(freq='H', year=2007, month=1, day=1, hour=0)
        ival_H_end_of_year = Period(freq='H', year=2007, month=12, day=31,
                                    hour=23)
        ival_H_end_of_quarter = Period(freq='H', year=2007, month=3, day=31,
                                        hour=23)
        ival_H_end_of_month = Period(freq='H', year=2007, month=1, day=31,
                                    hour=23)
        ival_H_end_of_week = Period(freq='H', year=2007, month=1, day=7,
                                    hour=23)
        ival_H_end_of_day = Period(freq='H', year=2007, month=1, day=1,
                                    hour=23)
        ival_H_end_of_bus = Period(freq='H', year=2007, month=1, day=1,
                                    hour=23)

        ival_H_to_A = Period(freq='A', year=2007)
        ival_H_to_Q = Period(freq='Q', year=2007, quarter=1)
        ival_H_to_M = Period(freq='M', year=2007, month=1)
        ival_H_to_W = Period(freq='WK', year=2007, month=1, day=7)
        ival_H_to_D = Period(freq='D', year=2007, month=1, day=1)
        ival_H_to_B = Period(freq='B', year=2007, month=1, day=1)

        ival_H_to_T_start = Period(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=0)
        ival_H_to_T_end = Period(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=59)
        ival_H_to_S_start = Period(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_H_to_S_end = Period(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=59, second=59)

        assert_equal(ival_H.asfreq('A'), ival_H_to_A)
        assert_equal(ival_H_end_of_year.asfreq('A'), ival_H_to_A)
        assert_equal(ival_H.asfreq('Q'), ival_H_to_Q)
        assert_equal(ival_H_end_of_quarter.asfreq('Q'), ival_H_to_Q)
        assert_equal(ival_H.asfreq('M'), ival_H_to_M)
        assert_equal(ival_H_end_of_month.asfreq('M'), ival_H_to_M)
        assert_equal(ival_H.asfreq('WK'), ival_H_to_W)
        assert_equal(ival_H_end_of_week.asfreq('WK'), ival_H_to_W)
        assert_equal(ival_H.asfreq('D'), ival_H_to_D)
        assert_equal(ival_H_end_of_day.asfreq('D'), ival_H_to_D)
        assert_equal(ival_H.asfreq('B'), ival_H_to_B)
        assert_equal(ival_H_end_of_bus.asfreq('B'), ival_H_to_B)

        assert_equal(ival_H.asfreq('Min', 'S'), ival_H_to_T_start)
        assert_equal(ival_H.asfreq('Min', 'E'), ival_H_to_T_end)
        assert_equal(ival_H.asfreq('S', 'S'), ival_H_to_S_start)
        assert_equal(ival_H.asfreq('S', 'E'), ival_H_to_S_end)

        assert_equal(ival_H.asfreq('H'), ival_H)

    def test_conv_minutely(self):
        # frequency conversion tests: from Minutely Frequency"

        ival_T = Period(freq='Min', year=2007, month=1, day=1,
                        hour=0, minute=0)
        ival_T_end_of_year = Period(freq='Min', year=2007, month=12, day=31,
                                    hour=23, minute=59)
        ival_T_end_of_quarter = Period(freq='Min', year=2007, month=3, day=31,
                                        hour=23, minute=59)
        ival_T_end_of_month = Period(freq='Min', year=2007, month=1, day=31,
                                    hour=23, minute=59)
        ival_T_end_of_week = Period(freq='Min', year=2007, month=1, day=7,
                                    hour=23, minute=59)
        ival_T_end_of_day = Period(freq='Min', year=2007, month=1, day=1,
                                    hour=23, minute=59)
        ival_T_end_of_bus = Period(freq='Min', year=2007, month=1, day=1,
                                    hour=23, minute=59)
        ival_T_end_of_hour = Period(freq='Min', year=2007, month=1, day=1,
                                    hour=0, minute=59)

        ival_T_to_A = Period(freq='A', year=2007)
        ival_T_to_Q = Period(freq='Q', year=2007, quarter=1)
        ival_T_to_M = Period(freq='M', year=2007, month=1)
        ival_T_to_W = Period(freq='WK', year=2007, month=1, day=7)
        ival_T_to_D = Period(freq='D', year=2007, month=1, day=1)
        ival_T_to_B = Period(freq='B', year=2007, month=1, day=1)
        ival_T_to_H = Period(freq='H', year=2007, month=1, day=1, hour=0)

        ival_T_to_S_start = Period(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=0)
        ival_T_to_S_end = Period(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=59)

        assert_equal(ival_T.asfreq('A'), ival_T_to_A)
        assert_equal(ival_T_end_of_year.asfreq('A'), ival_T_to_A)
        assert_equal(ival_T.asfreq('Q'), ival_T_to_Q)
        assert_equal(ival_T_end_of_quarter.asfreq('Q'), ival_T_to_Q)
        assert_equal(ival_T.asfreq('M'), ival_T_to_M)
        assert_equal(ival_T_end_of_month.asfreq('M'), ival_T_to_M)
        assert_equal(ival_T.asfreq('WK'), ival_T_to_W)
        assert_equal(ival_T_end_of_week.asfreq('WK'), ival_T_to_W)
        assert_equal(ival_T.asfreq('D'), ival_T_to_D)
        assert_equal(ival_T_end_of_day.asfreq('D'), ival_T_to_D)
        assert_equal(ival_T.asfreq('B'), ival_T_to_B)
        assert_equal(ival_T_end_of_bus.asfreq('B'), ival_T_to_B)
        assert_equal(ival_T.asfreq('H'), ival_T_to_H)
        assert_equal(ival_T_end_of_hour.asfreq('H'), ival_T_to_H)

        assert_equal(ival_T.asfreq('S', 'S'), ival_T_to_S_start)
        assert_equal(ival_T.asfreq('S', 'E'), ival_T_to_S_end)

        assert_equal(ival_T.asfreq('Min'), ival_T)

    def test_conv_secondly(self):
        # frequency conversion tests: from Secondly Frequency"

        ival_S = Period(freq='S', year=2007, month=1, day=1,
                        hour=0, minute=0, second=0)
        ival_S_end_of_year = Period(freq='S', year=2007, month=12, day=31,
                                    hour=23, minute=59, second=59)
        ival_S_end_of_quarter = Period(freq='S', year=2007, month=3, day=31,
                                        hour=23, minute=59, second=59)
        ival_S_end_of_month = Period(freq='S', year=2007, month=1, day=31,
                                    hour=23, minute=59, second=59)
        ival_S_end_of_week = Period(freq='S', year=2007, month=1, day=7,
                                    hour=23, minute=59, second=59)
        ival_S_end_of_day = Period(freq='S', year=2007, month=1, day=1,
                                    hour=23, minute=59, second=59)
        ival_S_end_of_bus = Period(freq='S', year=2007, month=1, day=1,
                                    hour=23, minute=59, second=59)
        ival_S_end_of_hour = Period(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=59, second=59)
        ival_S_end_of_minute = Period(freq='S', year=2007, month=1, day=1,
                                    hour=0, minute=0, second=59)

        ival_S_to_A = Period(freq='A', year=2007)
        ival_S_to_Q = Period(freq='Q', year=2007, quarter=1)
        ival_S_to_M = Period(freq='M', year=2007, month=1)
        ival_S_to_W = Period(freq='WK', year=2007, month=1, day=7)
        ival_S_to_D = Period(freq='D', year=2007, month=1, day=1)
        ival_S_to_B = Period(freq='B', year=2007, month=1, day=1)
        ival_S_to_H = Period(freq='H', year=2007, month=1, day=1,
                            hour=0)
        ival_S_to_T = Period(freq='Min', year=2007, month=1, day=1,
                            hour=0, minute=0)

        assert_equal(ival_S.asfreq('A'), ival_S_to_A)
        assert_equal(ival_S_end_of_year.asfreq('A'), ival_S_to_A)
        assert_equal(ival_S.asfreq('Q'), ival_S_to_Q)
        assert_equal(ival_S_end_of_quarter.asfreq('Q'), ival_S_to_Q)
        assert_equal(ival_S.asfreq('M'), ival_S_to_M)
        assert_equal(ival_S_end_of_month.asfreq('M'), ival_S_to_M)
        assert_equal(ival_S.asfreq('WK'), ival_S_to_W)
        assert_equal(ival_S_end_of_week.asfreq('WK'), ival_S_to_W)
        assert_equal(ival_S.asfreq('D'), ival_S_to_D)
        assert_equal(ival_S_end_of_day.asfreq('D'), ival_S_to_D)
        assert_equal(ival_S.asfreq('B'), ival_S_to_B)
        assert_equal(ival_S_end_of_bus.asfreq('B'), ival_S_to_B)
        assert_equal(ival_S.asfreq('H'), ival_S_to_H)
        assert_equal(ival_S_end_of_hour.asfreq('H'), ival_S_to_H)
        assert_equal(ival_S.asfreq('Min'), ival_S_to_T)
        assert_equal(ival_S_end_of_minute.asfreq('Min'), ival_S_to_T)

        assert_equal(ival_S.asfreq('S'), ival_S)

class TestPeriodIndex(TestCase):
    def __init__(self, *args, **kwds):
        TestCase.__init__(self, *args, **kwds)

    def setUp(self):
        pass

    def test_make_time_series(self):
        index = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        series = Series(1, index=index)
        self.assert_(isinstance(series, TimeSeries))

    def test_constructor_use_start_freq(self):
        # GH #1118
        p = Period('4/2/2012', freq='B')
        index = PeriodIndex(start=p, periods=10)
        expected = PeriodIndex(start='4/2/2012', periods=10, freq='B')
        self.assert_(index.equals(expected))

    def test_to_timestamp(self):
        index = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        series = Series(1, index=index, name='foo')

        exp_index = date_range('1/1/2001', end='12/31/2009', freq='A-DEC')
        result = series.to_timestamp('D', 'end')
        self.assert_(result.index.equals(exp_index))
        self.assertEquals(result.name, 'foo')

        exp_index = date_range('1/1/2001', end='1/1/2009', freq='AS-DEC')
        result = series.to_timestamp('D', 'start')
        self.assert_(result.index.equals(exp_index))


        def _get_with_delta(delta, freq='A-DEC'):
            return date_range(to_datetime('1/1/2001') + delta,
                              to_datetime('12/31/2009') + delta, freq=freq)

        delta = timedelta(hours=23)
        result = series.to_timestamp('H', 'end')
        exp_index = _get_with_delta(delta)
        self.assert_(result.index.equals(exp_index))

        delta = timedelta(hours=23, minutes=59)
        result = series.to_timestamp('T', 'end')
        exp_index = _get_with_delta(delta)
        self.assert_(result.index.equals(exp_index))

        result = series.to_timestamp('S', 'end')
        delta = timedelta(hours=23, minutes=59, seconds=59)
        exp_index = _get_with_delta(delta)
        self.assert_(result.index.equals(exp_index))

    def test_index_duplicate_periods(self):
        # monotonic
        idx = PeriodIndex([2000, 2007, 2007, 2009, 2009], freq='A-JUN')
        ts = Series(np.random.randn(len(idx)), index=idx)

        result = ts[2007]
        expected = ts[1:3]
        assert_series_equal(result, expected)
        result[:] = 1
        self.assert_((ts[1:3] == 1).all())

        # not monotonic
        idx = PeriodIndex([2000, 2007, 2007, 2009, 2007], freq='A-JUN')
        ts = Series(np.random.randn(len(idx)), index=idx)

        result = ts[2007]
        expected = ts[idx == 2007]
        assert_series_equal(result, expected)

    def test_constructor(self):
        ii = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        assert_equal(len(ii), 9)

        ii = PeriodIndex(freq='Q', start='1/1/2001', end='12/1/2009')
        assert_equal(len(ii), 4 * 9)

        ii = PeriodIndex(freq='M', start='1/1/2001', end='12/1/2009')
        assert_equal(len(ii), 12 * 9)

        ii = PeriodIndex(freq='D', start='1/1/2001', end='12/31/2009')
        assert_equal(len(ii), 365 * 9 + 2)

        ii = PeriodIndex(freq='B', start='1/1/2001', end='12/31/2009')
        assert_equal(len(ii), 261 * 9)

        ii = PeriodIndex(freq='H', start='1/1/2001', end='12/31/2001 23:00')
        assert_equal(len(ii), 365 * 24)

        ii = PeriodIndex(freq='Min', start='1/1/2001', end='1/1/2001 23:59')
        assert_equal(len(ii), 24 * 60)

        ii = PeriodIndex(freq='S', start='1/1/2001', end='1/1/2001 23:59:59')
        assert_equal(len(ii), 24 * 60 * 60)

        start = Period('02-Apr-2005', 'B')
        i1 = PeriodIndex(start=start, periods=20)
        assert_equal(len(i1), 20)
        assert_equal(i1.freq, start.freq)
        assert_equal(i1[0], start)

        end_intv = Period('2006-12-31', 'W')
        i1 = PeriodIndex(end=end_intv, periods=10)
        assert_equal(len(i1), 10)
        assert_equal(i1.freq, end_intv.freq)
        assert_equal(i1[-1], end_intv)

        end_intv = Period('2006-12-31', '1w')
        i2 = PeriodIndex(end=end_intv, periods=10)
        assert_equal(len(i1), len(i2))
        self.assert_((i1 == i2).all())
        assert_equal(i1.freq, i2.freq)

        end_intv = Period('2006-12-31', ('w', 1))
        i2 = PeriodIndex(end=end_intv, periods=10)
        assert_equal(len(i1), len(i2))
        self.assert_((i1 == i2).all())
        assert_equal(i1.freq, i2.freq)

        try:
            PeriodIndex(start=start, end=end_intv)
            raise AssertionError('Cannot allow mixed freq for start and end')
        except ValueError:
            pass

        end_intv = Period('2005-05-01', 'B')
        i1 = PeriodIndex(start=start, end=end_intv)

        try:
            PeriodIndex(start=start)
            raise AssertionError('Must specify periods if missing start or end')
        except ValueError:
            pass

        # infer freq from first element
        i2 = PeriodIndex([end_intv, Period('2005-05-05', 'B')])
        assert_equal(len(i2), 2)
        assert_equal(i2[0], end_intv)

        i2 = PeriodIndex(np.array([end_intv, Period('2005-05-05', 'B')]))
        assert_equal(len(i2), 2)
        assert_equal(i2[0], end_intv)

        # Mixed freq should fail
        vals = [end_intv, Period('2006-12-31', 'w')]
        self.assertRaises(ValueError, PeriodIndex, vals)
        vals = np.array(vals)
        self.assertRaises(ValueError, PeriodIndex, vals)

    def test_shift(self):
        ii1 = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        ii2 = PeriodIndex(freq='A', start='1/1/2002', end='12/1/2010')
        assert_equal(len(ii1), len(ii2))
        assert_equal(ii1.shift(1).values, ii2.values)

        ii1 = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        ii2 = PeriodIndex(freq='A', start='1/1/2000', end='12/1/2008')
        assert_equal(len(ii1), len(ii2))
        assert_equal(ii1.shift(-1).values, ii2.values)

        ii1 = PeriodIndex(freq='M', start='1/1/2001', end='12/1/2009')
        ii2 = PeriodIndex(freq='M', start='2/1/2001', end='1/1/2010')
        assert_equal(len(ii1), len(ii2))
        assert_equal(ii1.shift(1).values, ii2.values)

        ii1 = PeriodIndex(freq='M', start='1/1/2001', end='12/1/2009')
        ii2 = PeriodIndex(freq='M', start='12/1/2000', end='11/1/2009')
        assert_equal(len(ii1), len(ii2))
        assert_equal(ii1.shift(-1).values, ii2.values)

        ii1 = PeriodIndex(freq='D', start='1/1/2001', end='12/1/2009')
        ii2 = PeriodIndex(freq='D', start='1/2/2001', end='12/2/2009')
        assert_equal(len(ii1), len(ii2))
        assert_equal(ii1.shift(1).values, ii2.values)

        ii1 = PeriodIndex(freq='D', start='1/1/2001', end='12/1/2009')
        ii2 = PeriodIndex(freq='D', start='12/31/2000', end='11/30/2009')
        assert_equal(len(ii1), len(ii2))
        assert_equal(ii1.shift(-1).values, ii2.values)

    def test_asfreq(self):
        ii1 = PeriodIndex(freq='A', start='1/1/2001', end='1/1/2001')
        ii2 = PeriodIndex(freq='Q', start='1/1/2001', end='1/1/2001')
        ii3 = PeriodIndex(freq='M', start='1/1/2001', end='1/1/2001')
        ii4 = PeriodIndex(freq='D', start='1/1/2001', end='1/1/2001')
        ii5 = PeriodIndex(freq='H', start='1/1/2001', end='1/1/2001 00:00')
        ii6 = PeriodIndex(freq='Min', start='1/1/2001', end='1/1/2001 00:00')
        ii7 = PeriodIndex(freq='S', start='1/1/2001', end='1/1/2001 00:00:00')

        self.assertEquals(ii1.asfreq('Q', 'S'), ii2)
        self.assertEquals(ii1.asfreq('Q', 's'), ii2)
        self.assertEquals(ii1.asfreq('M', 'start'), ii3)
        self.assertEquals(ii1.asfreq('D', 'StarT'), ii4)
        self.assertEquals(ii1.asfreq('H', 'beGIN'), ii5)
        self.assertEquals(ii1.asfreq('Min', 'S'), ii6)
        self.assertEquals(ii1.asfreq('S', 'S'), ii7)

        self.assertEquals(ii2.asfreq('A', 'S'), ii1)
        self.assertEquals(ii2.asfreq('M', 'S'), ii3)
        self.assertEquals(ii2.asfreq('D', 'S'), ii4)
        self.assertEquals(ii2.asfreq('H', 'S'), ii5)
        self.assertEquals(ii2.asfreq('Min', 'S'), ii6)
        self.assertEquals(ii2.asfreq('S', 'S'), ii7)

        self.assertEquals(ii3.asfreq('A', 'S'), ii1)
        self.assertEquals(ii3.asfreq('Q', 'S'), ii2)
        self.assertEquals(ii3.asfreq('D', 'S'), ii4)
        self.assertEquals(ii3.asfreq('H', 'S'), ii5)
        self.assertEquals(ii3.asfreq('Min', 'S'), ii6)
        self.assertEquals(ii3.asfreq('S', 'S'), ii7)

        self.assertEquals(ii4.asfreq('A', 'S'), ii1)
        self.assertEquals(ii4.asfreq('Q', 'S'), ii2)
        self.assertEquals(ii4.asfreq('M', 'S'), ii3)
        self.assertEquals(ii4.asfreq('H', 'S'), ii5)
        self.assertEquals(ii4.asfreq('Min', 'S'), ii6)
        self.assertEquals(ii4.asfreq('S', 'S'), ii7)

        self.assertEquals(ii5.asfreq('A', 'S'), ii1)
        self.assertEquals(ii5.asfreq('Q', 'S'), ii2)
        self.assertEquals(ii5.asfreq('M', 'S'), ii3)
        self.assertEquals(ii5.asfreq('D', 'S'), ii4)
        self.assertEquals(ii5.asfreq('Min', 'S'), ii6)
        self.assertEquals(ii5.asfreq('S', 'S'), ii7)

        self.assertEquals(ii6.asfreq('A', 'S'), ii1)
        self.assertEquals(ii6.asfreq('Q', 'S'), ii2)
        self.assertEquals(ii6.asfreq('M', 'S'), ii3)
        self.assertEquals(ii6.asfreq('D', 'S'), ii4)
        self.assertEquals(ii6.asfreq('H', 'S'), ii5)
        self.assertEquals(ii6.asfreq('S', 'S'), ii7)

        self.assertEquals(ii7.asfreq('A', 'S'), ii1)
        self.assertEquals(ii7.asfreq('Q', 'S'), ii2)
        self.assertEquals(ii7.asfreq('M', 'S'), ii3)
        self.assertEquals(ii7.asfreq('D', 'S'), ii4)
        self.assertEquals(ii7.asfreq('H', 'S'), ii5)
        self.assertEquals(ii7.asfreq('Min', 'S'), ii6)

        #self.assertEquals(ii7.asfreq('A', 'E'), i_end)

    def test_ts_repr(self):
        index = PeriodIndex(freq='A', start='1/1/2001', end='12/31/2010')
        ts = Series(np.random.randn(len(index)), index=index)
        repr(ts)

    def test_asfreq_ts(self):
        index = PeriodIndex(freq='A', start='1/1/2001', end='12/31/2010')
        ts = Series(np.random.randn(len(index)), index=index)
        df = DataFrame(np.random.randn(len(index), 3), index=index)

        result = ts.asfreq('D', how='end')
        df_result = df.asfreq('D', how='end')
        exp_index = index.asfreq('D', how='end')
        self.assert_(len(result) == len(ts))
        self.assert_(result.index.equals(exp_index))
        self.assert_(df_result.index.equals(exp_index))

        result = ts.asfreq('D', how='start')
        self.assert_(len(result) == len(ts))
        self.assert_(result.index.equals(index.asfreq('D', how='start')))

    def test_badinput(self):
        self.assertRaises(datetools.DateParseError, Period, '1/1/-2000', 'A')
        self.assertRaises(ValueError, Period, -2000, 'A')
        self.assertRaises(ValueError, Period, 0, 'A')
        self.assertRaises(ValueError, PeriodIndex, [-1, 0, 1], 'A')
        self.assertRaises(ValueError, PeriodIndex, np.array([-1, 0, 1]), 'A')

    def test_dti_to_period(self):
        dti = DatetimeIndex(start='1/1/2005', end='12/1/2005', freq='M')
        ii1 = dti.to_period()
        ii2 = dti.to_period(freq='D')

        self.assertEquals(ii1[0], Period('Jan 2005', freq='M'))
        self.assertEquals(ii2[0], Period('1/31/2005', freq='D'))

        self.assertEquals(ii1[-1], Period('Nov 2005', freq='M'))
        self.assertEquals(ii2[-1], Period('11/30/2005', freq='D'))

    def test_iindex_slice_index(self):
        ii = PeriodIndex(start='1/1/10', end='12/31/12', freq='M')
        s = Series(np.random.rand(len(ii)), index=ii)
        res = s['2010']
        exp = s[0:12]
        assert_series_equal(res, exp)
        res = s['2011']
        exp = s[12:24]
        assert_series_equal(res, exp)

    def test_iindex_qaccess(self):
        ii = PeriodIndex(['2Q05', '3Q05', '4Q05', '1Q06', '2Q06'], freq='Q')
        s = Series(np.random.rand(len(ii)), index=ii).cumsum()
        # Todo: fix these accessors!
        self.assert_(s['05Q4'] == s[2])

    def test_period_dt64_round_trip(self):
        dti = date_range('1/1/2000', '1/7/2002', freq='B')
        ii = dti.to_period()
        self.assert_(ii.to_timestamp().equals(dti))

        dti = date_range('1/1/2000', '1/7/2002', freq='B')
        ii = dti.to_period(freq='H')
        self.assert_(ii.to_timestamp().equals(dti))

    def test_to_period_quarterly(self):
        # make sure we can make the round trip
        for month in MONTHS:
            freq = 'Q-%s' % month
            rng = period_range('1989Q3', '1991Q3', freq=freq)
            stamps = rng.to_timestamp()
            result = stamps.to_period(freq)
            self.assert_(rng.equals(result))

    def test_iindex_multiples(self):
        ii = PeriodIndex(start='1/1/10', end='12/31/12', freq='2M')
        self.assertEquals(ii[0], Period('1/1/10', '2M'))
        self.assertEquals(ii[1], Period('3/1/10', '2M'))

        self.assertEquals(ii[0].asfreq('6M'), ii[2].asfreq('6M'))
        self.assertEquals(ii[0].asfreq('A'), ii[2].asfreq('A'))

        self.assertEquals(ii[0].asfreq('M', how='S'),
                          Period('Jan 2010', '1M'))
        self.assertEquals(ii[0].asfreq('M', how='E'),
                          Period('Feb 2010', '1M'))
        self.assertEquals(ii[1].asfreq('M', how='S'),
                          Period('Mar 2010', '1M'))

        i = Period('1/1/2010 12:05:18', '5S')
        self.assertEquals(i, Period('1/1/2010 12:05:15', '5S'))

        i = Period('1/1/2010 12:05:18', '5S')
        self.assertEquals(i.asfreq('1S', how='E'),
                          Period('1/1/2010 12:05:19', '1S'))

    def test_iteration(self):
        index = PeriodIndex(start='1/1/10', periods=4, freq='B')

        result = list(index)
        self.assert_(isinstance(result[0], Period))
        self.assert_(result[0].freq == index.freq)

    def test_take(self):
        index = PeriodIndex(start='1/1/10', end='12/31/12', freq='D')

        taken = index.take([5, 6, 8, 12])
        taken2 = index[[5, 6, 8, 12]]
        self.assert_(isinstance(taken, PeriodIndex))
        self.assert_(taken.freq == index.freq)
        self.assert_(isinstance(taken2, PeriodIndex))
        self.assert_(taken2.freq == index.freq)

    def test_joins(self):
        index = period_range('1/1/2000', '1/20/2000', freq='D')

        for kind in ['inner', 'outer', 'left', 'right']:
            joined = index.join(index[:-5], how=kind)

            self.assert_(isinstance(joined, PeriodIndex))
            self.assert_(joined.freq == index.freq)

    def test_align_series(self):
        rng = period_range('1/1/2000', '1/1/2010', freq='A')
        ts = Series(np.random.randn(len(rng)), index=rng)

        result = ts + ts[::2]
        expected = ts + ts
        expected[1::2] = np.nan
        assert_series_equal(result, expected)

        result = ts + _permute(ts[::2])
        assert_series_equal(result, expected)

        # it works!
        for kind in ['inner', 'outer', 'left', 'right']:
            ts.align(ts[::2], join=kind)

        self.assertRaises(Exception, ts.__add__,
                          ts.asfreq('D', how='end'))

    def test_align_frame(self):
        rng = period_range('1/1/2000', '1/1/2010', freq='A')
        ts = DataFrame(np.random.randn(len(rng), 3), index=rng)

        result = ts + ts[::2]
        expected = ts + ts
        expected.values[1::2] = np.nan
        tm.assert_frame_equal(result, expected)

        result = ts + _permute(ts[::2])
        tm.assert_frame_equal(result, expected)

    def test_union(self):
        index = period_range('1/1/2000', '1/20/2000', freq='D')

        result = index[:-5].union(index[10:])
        self.assert_(result.equals(index))

        # not in order
        result = _permute(index[:-5]).union(_permute(index[10:]))
        self.assert_(result.equals(index))

        # raise if different frequencies
        index = period_range('1/1/2000', '1/20/2000', freq='D')
        index2 = period_range('1/1/2000', '1/20/2000', freq='W-WED')
        self.assertRaises(Exception, index.union, index2)

    def test_intersection(self):
        index = period_range('1/1/2000', '1/20/2000', freq='D')

        result = index[:-5].intersection(index[10:])
        self.assert_(result.equals(index[10:-5]))

        # not in order
        left = _permute(index[:-5])
        right = _permute(index[10:])
        result = left.intersection(right).order()
        self.assert_(result.equals(index[10:-5]))

        # raise if different frequencies
        index = period_range('1/1/2000', '1/20/2000', freq='D')
        index2 = period_range('1/1/2000', '1/20/2000', freq='W-WED')
        self.assertRaises(Exception, index.intersection, index2)

    def test_fields(self):
        # year, month, day, hour, minute
        # second, weekofyear, week, dayofweek, weekday, dayofyear, quarter
        # qyear
        ii = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        self._check_all_fields(ii)

        ii = PeriodIndex(freq='Q', start='1/1/2001', end='12/1/2003')
        self._check_all_fields(ii)

        ii = PeriodIndex(freq='M', start='1/1/2001', end='1/1/2002')
        self._check_all_fields(ii)

        ii = PeriodIndex(freq='D', start='12/1/2001', end='1/1/2002')
        self._check_all_fields(ii)

        ii = PeriodIndex(freq='B', start='12/1/2001', end='1/1/2002')
        self._check_all_fields(ii)

        ii = PeriodIndex(freq='H', start='12/31/2001', end='1/1/2002 23:00')
        self._check_all_fields(ii)

        ii = PeriodIndex(freq='Min', start='12/31/2001', end='1/1/2002 00:59')
        self._check_all_fields(ii)

        ii = PeriodIndex(freq='S', start='12/31/2001', end='1/1/2001 00:00:01')
        self._check_all_fields(ii)

        end_intv = Period('2006-12-31', 'W')
        i1 = PeriodIndex(end=end_intv, periods=10)
        self._check_all_fields(ii)

    def _check_all_fields(self, periodindex):
        fields = ['year', 'month', 'day', 'hour', 'minute',
                  'second', 'weekofyear', 'week', 'dayofweek',
                  'weekday', 'dayofyear', 'quarter', 'qyear']
        [self._check_field(periodindex, x) for x in fields]

    def _check_field(self, periodindex, fieldname):
        field_idx = getattr(periodindex, fieldname)
        assert_equal(len(periodindex), len(field_idx))
        for x, val in zip(periodindex, field_idx):
            assert_equal(getattr(x, fieldname), val)

def _permute(obj):
    return obj.take(np.random.permutation(len(obj)))

class TestMethods(TestCase):
    "Base test class for MaskedArrays."

    def __init__(self, *args, **kwds):
        TestCase.__init__(self, *args, **kwds)

    def test_add(self):
        dt1 = Period(freq='D', year=2008, month=1, day=1)
        dt2 = Period(freq='D', year=2008, month=1, day=2)
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
