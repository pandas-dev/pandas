"""Tests suite for Interval handling.

Parts derived from scikits.timeseries code, original authors:
- Pierre Gerard-Marchant & Matt Knox
- pierregm_at_uga_dot_edu - mattknow_ca_at_hotmail_dot_com

"""

from unittest import TestCase
from datetime import datetime

from numpy.ma.testutils import assert_equal
from pandas.core.datetools import Interval

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

        self.assertEquals(i1, i2)

        i3 = Interval('2005', freq='M')
        self.assert_(i1 != i3)

        i1 = Interval.now('Q')
        i2 = Interval(datetime.now(), freq='Q')

        self.assertEquals(i1, i2)

        # Biz day construction, roll forward if non-weekday
        i1 = Interval('3/10/12', freq='B')
        i2 = Interval('3/12/12', freq='D')
        self.assertEquals(i1, i2.asfreq('B'))

        i1 = Interval(year=2005, quarter=1, freq='Q')
        i2 = Interval('1/1/2005', freq='Q')
        self.assertEquals(i1, i2)

        i1 = Interval(year=2005, quarter=3, freq='Q')
        i2 = Interval('9/1/2005', freq='Q')
        self.assertEquals(i1, i2)

        i1 = Interval(year=2005, month=3, day=1, freq='D')
        i2 = Interval('3/1/2005', freq='D')
        self.assertEquals(i1, i2)

        i1 = Interval(year=2012, month=3, day=10, freq='B')
        i2 = Interval('3/12/12', freq='B')
        self.assertEquals(i1, i2)

    def test_properties_annually(self):
        "Test properties on Intervals with annually frequency."
        a_date = Interval(freq='A', year=2007)
        assert_equal(a_date.year, 2007)

    def test_properties_quarterly(self):
        "Test properties on Intervals with daily frequency."
        qedec_date = Interval(freq="Q@DEC", year=2007, quarter=1)
        qejan_date = Interval(freq="Q@JAN", year=2007, quarter=1)
        qejun_date = Interval(freq="Q@JUN", year=2007, quarter=1)
        qsdec_date = Interval(freq="QS@DEC", year=2007, quarter=1)
        qsjan_date = Interval(freq="QS@JAN", year=2007, quarter=1)
        qsjun_date = Interval(freq="QS@JUN", year=2007, quarter=1)
        #
        for x in range(3):
            for qd in (qedec_date, qejan_date, qejun_date,
                       qsdec_date, qsjan_date, qsjun_date):
                assert_equal((qd + x).qyear, 2007)
                assert_equal((qd + x).quarter, x + 1)


    def test_properties_monthly(self):
        "Test properties on Intervals with daily frequency."
        m_date = Interval(freq='M', year=2007, month=1)
        for x in range(11):
            m_date_x = m_date + x
            assert_equal(m_date_x.year, 2007)
            if 1 <= x + 1 <= 3:
                assert_equal(m_date_x.quarter, 1)
            elif 4 <= x + 1 <= 6:
                assert_equal(m_date_x.quarter, 2)
            elif 7 <= x + 1 <= 9:
                assert_equal(m_date_x.quarter, 3)
            elif 10 <= x + 1 <= 12:
                assert_equal(m_date_x.quarter, 4)
            assert_equal(m_date_x.month, x + 1)


    def test_properties_weekly(self):
        "Test properties on Intervals with daily frequency."
        w_date = Interval(freq='WK', year=2007, month=1, day=7)
        #
        assert_equal(w_date.year, 2007)
        assert_equal(w_date.quarter, 1)
        assert_equal(w_date.month, 1)
        assert_equal(w_date.week, 1)
        assert_equal((w_date - 1).week, 52)


    def test_properties_daily(self):
        "Test properties on Intervals with daily frequency."
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
        "Test properties on Intervals with hourly frequency."
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
        "Test properties on Intervals with minutely frequency."
        t_date = Interval(freq='Min', year=2007, month=1, day=1, hour=0, minute=0)
        #
        assert_equal(t_date.quarter, 1)
        assert_equal(t_date.month, 1)
        assert_equal(t_date.day, 1)
        assert_equal(t_date.weekday, 0)
        assert_equal(t_date.day_of_year, 1)
        assert_equal(t_date.hour, 0)
        assert_equal(t_date.minute, 0)


    def test_properties_secondly(self):
        "Test properties on Intervals with secondly frequency."
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
        self.dateWrap = [(noWrap, assert_equal)]

    def test_conv_annual(self):
        "frequency conversion tests: from Annual Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_A = dWrap(Interval(freq='A', year=2007))

            date_AJAN = dWrap(Interval(freq="A@JAN", year=2007))
            date_AJUN = dWrap(Interval(freq="A@JUN", year=2007))
            date_ANOV = dWrap(Interval(freq="A@NOV", year=2007))

            date_A_to_Q_start = dWrap(Interval(freq='Q', year=2007, quarter=1))
            date_A_to_Q_end = dWrap(Interval(freq='Q', year=2007, quarter=4))
            date_A_to_M_start = dWrap(Interval(freq='M', year=2007, month=1))
            date_A_to_M_end = dWrap(Interval(freq='M', year=2007, month=12))
            date_A_to_W_start = dWrap(Interval(freq='WK', year=2007, month=1, day=1))
            date_A_to_W_end = dWrap(Interval(freq='WK', year=2007, month=12, day=31))
            date_A_to_B_start = dWrap(Interval(freq='B', year=2007, month=1, day=1))
            date_A_to_B_end = dWrap(Interval(freq='B', year=2007, month=12, day=31))
            date_A_to_D_start = dWrap(Interval(freq='D', year=2007, month=1, day=1))
            date_A_to_D_end = dWrap(Interval(freq='D', year=2007, month=12, day=31))
            date_A_to_H_start = dWrap(Interval(freq='H', year=2007, month=1, day=1,
                                      hour=0))
            date_A_to_H_end = dWrap(Interval(freq='H', year=2007, month=12, day=31,
                                     hour=23))
            date_A_to_T_start = dWrap(Interval(freq='Min', year=2007, month=1, day=1,
                                      hour=0, minute=0))
            date_A_to_T_end = dWrap(Interval(freq='Min', year=2007, month=12, day=31,
                                     hour=23, minute=59))
            date_A_to_S_start = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_A_to_S_end = dWrap(Interval(freq='S', year=2007, month=12, day=31,
                                     hour=23, minute=59, second=59))

            date_AJAN_to_D_end = dWrap(Interval(freq='D', year=2007, month=1, day=31))
            date_AJAN_to_D_start = dWrap(Interval(freq='D', year=2006, month=2, day=1))
            date_AJUN_to_D_end = dWrap(Interval(freq='D', year=2007, month=6, day=30))
            date_AJUN_to_D_start = dWrap(Interval(freq='D', year=2006, month=7, day=1))
            date_ANOV_to_D_end = dWrap(Interval(freq='D', year=2007, month=11, day=30))
            date_ANOV_to_D_start = dWrap(Interval(freq='D', year=2006, month=12, day=1))

            assert_func(date_A.asfreq('Q', 'S'), date_A_to_Q_start)
            assert_func(date_A.asfreq('Q', 'E'), date_A_to_Q_end)
            assert_func(date_A.asfreq('M', 'S'), date_A_to_M_start)
            assert_func(date_A.asfreq('M', 'E'), date_A_to_M_end)
            assert_func(date_A.asfreq('WK', 'S'), date_A_to_W_start)
            assert_func(date_A.asfreq('WK', 'E'), date_A_to_W_end)
            assert_func(date_A.asfreq('B', 'S'), date_A_to_B_start)
            assert_func(date_A.asfreq('B', 'E'), date_A_to_B_end)
            assert_func(date_A.asfreq('D', 'S'), date_A_to_D_start)
            assert_func(date_A.asfreq('D', 'E'), date_A_to_D_end)
            assert_func(date_A.asfreq('H', 'S'), date_A_to_H_start)
            assert_func(date_A.asfreq('H', 'E'), date_A_to_H_end)
            assert_func(date_A.asfreq('Min', 'S'), date_A_to_T_start)
            assert_func(date_A.asfreq('Min', 'E'), date_A_to_T_end)
            assert_func(date_A.asfreq('S', 'S'), date_A_to_S_start)
            assert_func(date_A.asfreq('S', 'E'), date_A_to_S_end)

            assert_func(date_AJAN.asfreq('D', 'S'), date_AJAN_to_D_start)
            assert_func(date_AJAN.asfreq('D', 'E'), date_AJAN_to_D_end)

            assert_func(date_AJUN.asfreq('D', 'S'), date_AJUN_to_D_start)
            assert_func(date_AJUN.asfreq('D', 'E'), date_AJUN_to_D_end)

            assert_func(date_ANOV.asfreq('D', 'S'), date_ANOV_to_D_start)
            assert_func(date_ANOV.asfreq('D', 'E'), date_ANOV_to_D_end)

            assert_func(date_A.asfreq('A'), date_A)


    def test_conv_quarterly(self):
        "frequency conversion tests: from Quarterly Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_Q = dWrap(Interval(freq='Q', year=2007, quarter=1))
            date_Q_end_of_year = dWrap(Interval(freq='Q', year=2007, quarter=4))

            date_QEJAN = dWrap(Interval(freq="Q@JAN", year=2007, quarter=1))
            date_QEJUN = dWrap(Interval(freq="Q@JUN", year=2007, quarter=1))

            date_QSJAN = dWrap(Interval(freq="QS@JAN", year=2007, quarter=1))
            date_QSJUN = dWrap(Interval(freq="QS@JUN", year=2007, quarter=1))
            date_QSDEC = dWrap(Interval(freq="QS@DEC", year=2007, quarter=1))

            date_Q_to_A = dWrap(Interval(freq='A', year=2007))
            date_Q_to_M_start = dWrap(Interval(freq='M', year=2007, month=1))
            date_Q_to_M_end = dWrap(Interval(freq='M', year=2007, month=3))
            date_Q_to_W_start = dWrap(Interval(freq='WK', year=2007, month=1, day=1))
            date_Q_to_W_end = dWrap(Interval(freq='WK', year=2007, month=3, day=31))
            date_Q_to_B_start = dWrap(Interval(freq='B', year=2007, month=1, day=1))
            date_Q_to_B_end = dWrap(Interval(freq='B', year=2007, month=3, day=30))
            date_Q_to_D_start = dWrap(Interval(freq='D', year=2007, month=1, day=1))
            date_Q_to_D_end = dWrap(Interval(freq='D', year=2007, month=3, day=31))
            date_Q_to_H_start = dWrap(Interval(freq='H', year=2007, month=1, day=1,
                                      hour=0))
            date_Q_to_H_end = dWrap(Interval(freq='H', year=2007, month=3, day=31,
                                     hour=23))
            date_Q_to_T_start = dWrap(Interval(freq='Min', year=2007, month=1, day=1,
                                      hour=0, minute=0))
            date_Q_to_T_end = dWrap(Interval(freq='Min', year=2007, month=3, day=31,
                                     hour=23, minute=59))
            date_Q_to_S_start = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_Q_to_S_end = dWrap(Interval(freq='S', year=2007, month=3, day=31,
                                     hour=23, minute=59, second=59))

            date_QEJAN_to_D_start = dWrap(Interval(freq='D', year=2006, month=2, day=1))
            date_QEJAN_to_D_end = dWrap(Interval(freq='D', year=2006, month=4, day=30))

            date_QEJUN_to_D_start = dWrap(Interval(freq='D', year=2006, month=7, day=1))
            date_QEJUN_to_D_end = dWrap(Interval(freq='D', year=2006, month=9, day=30))

            date_QSJAN_to_D_start = dWrap(Interval(freq='D', year=2007, month=2, day=1))
            date_QSJAN_to_D_end = dWrap(Interval(freq='D', year=2007, month=4, day=30))

            date_QSJUN_to_D_start = dWrap(Interval(freq='D', year=2007, month=7, day=1))
            date_QSJUN_to_D_end = dWrap(Interval(freq='D', year=2007, month=9, day=30))

            date_QSDEC_to_D_start = dWrap(Interval(freq='D', year=2007, month=1, day=1))
            date_QSDEC_to_D_end = dWrap(Interval(freq='D', year=2007, month=3, day=31))

            assert_func(date_Q.asfreq('A'), date_Q_to_A)
            assert_func(date_Q_end_of_year.asfreq('A'), date_Q_to_A)

            assert_func(date_Q.asfreq('M', 'S'), date_Q_to_M_start)
            assert_func(date_Q.asfreq('M', 'E'), date_Q_to_M_end)
            assert_func(date_Q.asfreq('WK', 'S'), date_Q_to_W_start)
            assert_func(date_Q.asfreq('WK', 'E'), date_Q_to_W_end)
            assert_func(date_Q.asfreq('B', 'S'), date_Q_to_B_start)
            assert_func(date_Q.asfreq('B', 'E'), date_Q_to_B_end)
            assert_func(date_Q.asfreq('D', 'S'), date_Q_to_D_start)
            assert_func(date_Q.asfreq('D', 'E'), date_Q_to_D_end)
            assert_func(date_Q.asfreq('H', 'S'), date_Q_to_H_start)
            assert_func(date_Q.asfreq('H', 'E'), date_Q_to_H_end)
            assert_func(date_Q.asfreq('Min', 'S'), date_Q_to_T_start)
            assert_func(date_Q.asfreq('Min', 'E'), date_Q_to_T_end)
            assert_func(date_Q.asfreq('S', 'S'), date_Q_to_S_start)
            assert_func(date_Q.asfreq('S', 'E'), date_Q_to_S_end)

            assert_func(date_QEJAN.asfreq('D', 'S'), date_QEJAN_to_D_start)
            assert_func(date_QEJAN.asfreq('D', 'E'), date_QEJAN_to_D_end)
            assert_func(date_QEJUN.asfreq('D', 'S'), date_QEJUN_to_D_start)
            assert_func(date_QEJUN.asfreq('D', 'E'), date_QEJUN_to_D_end)

            assert_func(date_QSJAN.asfreq('D', 'S'), date_QSJAN_to_D_start)
            assert_func(date_QSJAN.asfreq('D', 'E'), date_QSJAN_to_D_end)
            assert_func(date_QSJUN.asfreq('D', 'S'), date_QSJUN_to_D_start)
            assert_func(date_QSJUN.asfreq('D', 'E'), date_QSJUN_to_D_end)
            assert_func(date_QSDEC.asfreq('D', 'S'), date_QSDEC_to_D_start)
            assert_func(date_QSDEC.asfreq('D', 'E'), date_QSDEC_to_D_end)

            assert_func(date_Q.asfreq('Q'), date_Q)


    def test_conv_monthly(self):
        "frequency conversion tests: from Monthly Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_M = dWrap(Interval(freq='M', year=2007, month=1))
            date_M_end_of_year = dWrap(Interval(freq='M', year=2007, month=12))
            date_M_end_of_quarter = dWrap(Interval(freq='M', year=2007, month=3))
            date_M_to_A = dWrap(Interval(freq='A', year=2007))
            date_M_to_Q = dWrap(Interval(freq='Q', year=2007, quarter=1))
            date_M_to_W_start = dWrap(Interval(freq='WK', year=2007, month=1, day=1))
            date_M_to_W_end = dWrap(Interval(freq='WK', year=2007, month=1, day=31))
            date_M_to_B_start = dWrap(Interval(freq='B', year=2007, month=1, day=1))
            date_M_to_B_end = dWrap(Interval(freq='B', year=2007, month=1, day=31))
            date_M_to_D_start = dWrap(Interval(freq='D', year=2007, month=1, day=1))
            date_M_to_D_end = dWrap(Interval(freq='D', year=2007, month=1, day=31))
            date_M_to_H_start = dWrap(Interval(freq='H', year=2007, month=1, day=1,
                                      hour=0))
            date_M_to_H_end = dWrap(Interval(freq='H', year=2007, month=1, day=31,
                                     hour=23))
            date_M_to_T_start = dWrap(Interval(freq='Min', year=2007, month=1, day=1,
                                      hour=0, minute=0))
            date_M_to_T_end = dWrap(Interval(freq='Min', year=2007, month=1, day=31,
                                     hour=23, minute=59))
            date_M_to_S_start = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_M_to_S_end = dWrap(Interval(freq='S', year=2007, month=1, day=31,
                                     hour=23, minute=59, second=59))

            assert_func(date_M.asfreq('A'), date_M_to_A)
            assert_func(date_M_end_of_year.asfreq('A'), date_M_to_A)
            assert_func(date_M.asfreq('Q'), date_M_to_Q)
            assert_func(date_M_end_of_quarter.asfreq('Q'), date_M_to_Q)

            assert_func(date_M.asfreq('WK', 'S'), date_M_to_W_start)
            assert_func(date_M.asfreq('WK', 'E'), date_M_to_W_end)
            assert_func(date_M.asfreq('B', 'S'), date_M_to_B_start)
            assert_func(date_M.asfreq('B', 'E'), date_M_to_B_end)
            assert_func(date_M.asfreq('D', 'S'), date_M_to_D_start)
            assert_func(date_M.asfreq('D', 'E'), date_M_to_D_end)
            assert_func(date_M.asfreq('H', 'S'), date_M_to_H_start)
            assert_func(date_M.asfreq('H', 'E'), date_M_to_H_end)
            assert_func(date_M.asfreq('Min', 'S'), date_M_to_T_start)
            assert_func(date_M.asfreq('Min', 'E'), date_M_to_T_end)
            assert_func(date_M.asfreq('S', 'S'), date_M_to_S_start)
            assert_func(date_M.asfreq('S', 'E'), date_M_to_S_end)

            assert_func(date_M.asfreq('M'), date_M)


    def test_conv_weekly(self):
        "frequency conversion tests: from Weekly Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_W = dWrap(Interval(freq='WK', year=2007, month=1, day=1))

            date_WSUN = dWrap(Interval(freq='WK', year=2007, month=1, day=7))
            date_WSAT = dWrap(Interval(freq='WK@SAT', year=2007, month=1, day=6))
            date_WFRI = dWrap(Interval(freq='WK@FRI', year=2007, month=1, day=5))
            date_WTHU = dWrap(Interval(freq='WK@THU', year=2007, month=1, day=4))
            date_WWED = dWrap(Interval(freq='WK@WED', year=2007, month=1, day=3))
            date_WTUE = dWrap(Interval(freq='WK@TUE', year=2007, month=1, day=2))
            date_WMON = dWrap(Interval(freq='WK@MON', year=2007, month=1, day=1))

            date_WSUN_to_D_start = dWrap(Interval(freq='D', year=2007, month=1, day=1))
            date_WSUN_to_D_end = dWrap(Interval(freq='D', year=2007, month=1, day=7))
            date_WSAT_to_D_start = dWrap(Interval(freq='D', year=2006, month=12, day=31))
            date_WSAT_to_D_end = dWrap(Interval(freq='D', year=2007, month=1, day=6))
            date_WFRI_to_D_start = dWrap(Interval(freq='D', year=2006, month=12, day=30))
            date_WFRI_to_D_end = dWrap(Interval(freq='D', year=2007, month=1, day=5))
            date_WTHU_to_D_start = dWrap(Interval(freq='D', year=2006, month=12, day=29))
            date_WTHU_to_D_end = dWrap(Interval(freq='D', year=2007, month=1, day=4))
            date_WWED_to_D_start = dWrap(Interval(freq='D', year=2006, month=12, day=28))
            date_WWED_to_D_end = dWrap(Interval(freq='D', year=2007, month=1, day=3))
            date_WTUE_to_D_start = dWrap(Interval(freq='D', year=2006, month=12, day=27))
            date_WTUE_to_D_end = dWrap(Interval(freq='D', year=2007, month=1, day=2))
            date_WMON_to_D_start = dWrap(Interval(freq='D', year=2006, month=12, day=26))
            date_WMON_to_D_end = dWrap(Interval(freq='D', year=2007, month=1, day=1))

            date_W_end_of_year = dWrap(Interval(freq='WK', year=2007, month=12, day=31))
            date_W_end_of_quarter = dWrap(Interval(freq='WK', year=2007, month=3, day=31))
            date_W_end_of_month = dWrap(Interval(freq='WK', year=2007, month=1, day=31))
            date_W_to_A = dWrap(Interval(freq='A', year=2007))
            date_W_to_Q = dWrap(Interval(freq='Q', year=2007, quarter=1))
            date_W_to_M = dWrap(Interval(freq='M', year=2007, month=1))

            if Interval(freq='D', year=2007, month=12, day=31).weekday == 6:
                date_W_to_A_end_of_year = dWrap(Interval(freq='A', year=2007))
            else:
                date_W_to_A_end_of_year = dWrap(Interval(freq='A', year=2008))

            if Interval(freq='D', year=2007, month=3, day=31).weekday == 6:
                date_W_to_Q_end_of_quarter = dWrap(Interval(freq='Q', year=2007, quarter=1))
            else:
                date_W_to_Q_end_of_quarter = dWrap(Interval(freq='Q', year=2007, quarter=2))

            if Interval(freq='D', year=2007, month=1, day=31).weekday == 6:
                date_W_to_M_end_of_month = dWrap(Interval(freq='M', year=2007, month=1))
            else:
                date_W_to_M_end_of_month = dWrap(Interval(freq='M', year=2007, month=2))

            date_W_to_B_start = dWrap(Interval(freq='B', year=2007, month=1, day=1))
            date_W_to_B_end = dWrap(Interval(freq='B', year=2007, month=1, day=5))
            date_W_to_D_start = dWrap(Interval(freq='D', year=2007, month=1, day=1))
            date_W_to_D_end = dWrap(Interval(freq='D', year=2007, month=1, day=7))
            date_W_to_H_start = dWrap(Interval(freq='H', year=2007, month=1, day=1,
                                      hour=0))
            date_W_to_H_end = dWrap(Interval(freq='H', year=2007, month=1, day=7,
                                     hour=23))
            date_W_to_T_start = dWrap(Interval(freq='Min', year=2007, month=1, day=1,
                                      hour=0, minute=0))
            date_W_to_T_end = dWrap(Interval(freq='Min', year=2007, month=1, day=7,
                                     hour=23, minute=59))
            date_W_to_S_start = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_W_to_S_end = dWrap(Interval(freq='S', year=2007, month=1, day=7,
                                     hour=23, minute=59, second=59))

            assert_func(date_W.asfreq('A'), date_W_to_A)
            assert_func(date_W_end_of_year.asfreq('A'), date_W_to_A_end_of_year)
            assert_func(date_W.asfreq('Q'), date_W_to_Q)
            assert_func(date_W_end_of_quarter.asfreq('Q'), date_W_to_Q_end_of_quarter)
            assert_func(date_W.asfreq('M'), date_W_to_M)
            assert_func(date_W_end_of_month.asfreq('M'), date_W_to_M_end_of_month)

            assert_func(date_W.asfreq('B', 'S'), date_W_to_B_start)
            assert_func(date_W.asfreq('B', 'E'), date_W_to_B_end)

            assert_func(date_W.asfreq('D', 'S'), date_W_to_D_start)
            assert_func(date_W.asfreq('D', 'E'), date_W_to_D_end)

            assert_func(date_WSUN.asfreq('D', 'S'), date_WSUN_to_D_start)
            assert_func(date_WSUN.asfreq('D', 'E'), date_WSUN_to_D_end)
            assert_func(date_WSAT.asfreq('D', 'S'), date_WSAT_to_D_start)
            assert_func(date_WSAT.asfreq('D', 'E'), date_WSAT_to_D_end)
            assert_func(date_WFRI.asfreq('D', 'S'), date_WFRI_to_D_start)
            assert_func(date_WFRI.asfreq('D', 'E'), date_WFRI_to_D_end)
            assert_func(date_WTHU.asfreq('D', 'S'), date_WTHU_to_D_start)
            assert_func(date_WTHU.asfreq('D', 'E'), date_WTHU_to_D_end)
            assert_func(date_WWED.asfreq('D', 'S'), date_WWED_to_D_start)
            assert_func(date_WWED.asfreq('D', 'E'), date_WWED_to_D_end)
            assert_func(date_WTUE.asfreq('D', 'S'), date_WTUE_to_D_start)
            assert_func(date_WTUE.asfreq('D', 'E'), date_WTUE_to_D_end)
            assert_func(date_WMON.asfreq('D', 'S'), date_WMON_to_D_start)
            assert_func(date_WMON.asfreq('D', 'E'), date_WMON_to_D_end)

            assert_func(date_W.asfreq('H', 'S'), date_W_to_H_start)
            assert_func(date_W.asfreq('H', 'E'), date_W_to_H_end)
            assert_func(date_W.asfreq('Min', 'S'), date_W_to_T_start)
            assert_func(date_W.asfreq('Min', 'E'), date_W_to_T_end)
            assert_func(date_W.asfreq('S', 'S'), date_W_to_S_start)
            assert_func(date_W.asfreq('S', 'E'), date_W_to_S_end)

            assert_func(date_W.asfreq('WK'), date_W)


    def test_conv_business(self):
        "frequency conversion tests: from Business Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_B = dWrap(Interval(freq='B', year=2007, month=1, day=1))
            date_B_end_of_year = dWrap(Interval(freq='B', year=2007, month=12, day=31))
            date_B_end_of_quarter = dWrap(Interval(freq='B', year=2007, month=3, day=30))
            date_B_end_of_month = dWrap(Interval(freq='B', year=2007, month=1, day=31))
            date_B_end_of_week = dWrap(Interval(freq='B', year=2007, month=1, day=5))

            date_B_to_A = dWrap(Interval(freq='A', year=2007))
            date_B_to_Q = dWrap(Interval(freq='Q', year=2007, quarter=1))
            date_B_to_M = dWrap(Interval(freq='M', year=2007, month=1))
            date_B_to_W = dWrap(Interval(freq='WK', year=2007, month=1, day=7))
            date_B_to_D = dWrap(Interval(freq='D', year=2007, month=1, day=1))
            date_B_to_H_start = dWrap(Interval(freq='H', year=2007, month=1, day=1,
                                      hour=0))
            date_B_to_H_end = dWrap(Interval(freq='H', year=2007, month=1, day=1,
                                     hour=23))
            date_B_to_T_start = dWrap(Interval(freq='Min', year=2007, month=1, day=1,
                                      hour=0, minute=0))
            date_B_to_T_end = dWrap(Interval(freq='Min', year=2007, month=1, day=1,
                                     hour=23, minute=59))
            date_B_to_S_start = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_B_to_S_end = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                     hour=23, minute=59, second=59))

            assert_func(date_B.asfreq('A'), date_B_to_A)
            assert_func(date_B_end_of_year.asfreq('A'), date_B_to_A)
            assert_func(date_B.asfreq('Q'), date_B_to_Q)
            assert_func(date_B_end_of_quarter.asfreq('Q'), date_B_to_Q)
            assert_func(date_B.asfreq('M'), date_B_to_M)
            assert_func(date_B_end_of_month.asfreq('M'), date_B_to_M)
            assert_func(date_B.asfreq('WK'), date_B_to_W)
            assert_func(date_B_end_of_week.asfreq('WK'), date_B_to_W)

            assert_func(date_B.asfreq('D'), date_B_to_D)

            assert_func(date_B.asfreq('H', 'S'), date_B_to_H_start)
            assert_func(date_B.asfreq('H', 'E'), date_B_to_H_end)
            assert_func(date_B.asfreq('Min', 'S'), date_B_to_T_start)
            assert_func(date_B.asfreq('Min', 'E'), date_B_to_T_end)
            assert_func(date_B.asfreq('S', 'S'), date_B_to_S_start)
            assert_func(date_B.asfreq('S', 'E'), date_B_to_S_end)

            assert_func(date_B.asfreq('B'), date_B)


    def test_conv_daily(self):
        "frequency conversion tests: from Business Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_D = dWrap(Interval(freq='D', year=2007, month=1, day=1))
            date_D_end_of_year = dWrap(Interval(freq='D', year=2007, month=12, day=31))
            date_D_end_of_quarter = dWrap(Interval(freq='D', year=2007, month=3, day=31))
            date_D_end_of_month = dWrap(Interval(freq='D', year=2007, month=1, day=31))
            date_D_end_of_week = dWrap(Interval(freq='D', year=2007, month=1, day=7))

            date_D_friday = dWrap(Interval(freq='D', year=2007, month=1, day=5))
            date_D_saturday = dWrap(Interval(freq='D', year=2007, month=1, day=6))
            date_D_sunday = dWrap(Interval(freq='D', year=2007, month=1, day=7))
            date_D_monday = dWrap(Interval(freq='D', year=2007, month=1, day=8))

            date_B_friday = dWrap(Interval(freq='B', year=2007, month=1, day=5))
            date_B_monday = dWrap(Interval(freq='B', year=2007, month=1, day=8))

            date_D_to_A = dWrap(Interval(freq='A', year=2007))

            date_Deoq_to_AJAN = dWrap(Interval(freq='A@JAN', year=2008))
            date_Deoq_to_AJUN = dWrap(Interval(freq='A@JUN', year=2007))
            date_Deoq_to_ADEC = dWrap(Interval(freq='A@DEC', year=2007))

            date_D_to_QEJAN = dWrap(Interval(freq="Q@JAN", year=2007, quarter=4))
            date_D_to_QEJUN = dWrap(Interval(freq="Q@JUN", year=2007, quarter=3))
            date_D_to_QEDEC = dWrap(Interval(freq="Q@DEC", year=2007, quarter=1))

            date_D_to_QSJAN = dWrap(Interval(freq="QS@JAN", year=2006, quarter=4))
            date_D_to_QSJUN = dWrap(Interval(freq="QS@JUN", year=2006, quarter=3))
            date_D_to_QSDEC = dWrap(Interval(freq="QS@DEC", year=2007, quarter=1))

            date_D_to_M = dWrap(Interval(freq='M', year=2007, month=1))
            date_D_to_W = dWrap(Interval(freq='WK', year=2007, month=1, day=7))

            date_D_to_H_start = dWrap(Interval(freq='H', year=2007, month=1, day=1,
                                      hour=0))
            date_D_to_H_end = dWrap(Interval(freq='H', year=2007, month=1, day=1,
                                     hour=23))
            date_D_to_T_start = dWrap(Interval(freq='Min', year=2007, month=1, day=1,
                                      hour=0, minute=0))
            date_D_to_T_end = dWrap(Interval(freq='Min', year=2007, month=1, day=1,
                                     hour=23, minute=59))
            date_D_to_S_start = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_D_to_S_end = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                     hour=23, minute=59, second=59))

            assert_func(date_D.asfreq('A'), date_D_to_A)

            assert_func(date_D_end_of_quarter.asfreq('A@JAN'), date_Deoq_to_AJAN)
            assert_func(date_D_end_of_quarter.asfreq('A@JUN'), date_Deoq_to_AJUN)
            assert_func(date_D_end_of_quarter.asfreq('A@DEC'), date_Deoq_to_ADEC)

            assert_func(date_D_end_of_year.asfreq('A'), date_D_to_A)
            assert_func(date_D_end_of_quarter.asfreq('Q'), date_D_to_QEDEC)
            assert_func(date_D.asfreq("Q@JAN"), date_D_to_QEJAN)
            assert_func(date_D.asfreq("Q@JUN"), date_D_to_QEJUN)
            assert_func(date_D.asfreq("Q@DEC"), date_D_to_QEDEC)
            assert_func(date_D.asfreq("QS@JAN"), date_D_to_QSJAN)
            assert_func(date_D.asfreq("QS@JUN"), date_D_to_QSJUN)
            assert_func(date_D.asfreq("QS@DEC"), date_D_to_QSDEC)
            assert_func(date_D.asfreq('M'), date_D_to_M)
            assert_func(date_D_end_of_month.asfreq('M'), date_D_to_M)
            assert_func(date_D.asfreq('WK'), date_D_to_W)
            assert_func(date_D_end_of_week.asfreq('WK'), date_D_to_W)

            assert_func(date_D_friday.asfreq('B'), date_B_friday)
            assert_func(date_D_saturday.asfreq('B', 'S'), date_B_friday)
            assert_func(date_D_saturday.asfreq('B', 'E'), date_B_monday)
            assert_func(date_D_sunday.asfreq('B', 'S'), date_B_friday)
            assert_func(date_D_sunday.asfreq('B', 'E'), date_B_monday)

            assert_func(date_D.asfreq('H', 'S'), date_D_to_H_start)
            assert_func(date_D.asfreq('H', 'E'), date_D_to_H_end)
            assert_func(date_D.asfreq('Min', 'S'), date_D_to_T_start)
            assert_func(date_D.asfreq('Min', 'E'), date_D_to_T_end)
            assert_func(date_D.asfreq('S', 'S'), date_D_to_S_start)
            assert_func(date_D.asfreq('S', 'E'), date_D_to_S_end)

            assert_func(date_D.asfreq('D'), date_D)

    def test_conv_hourly(self):
        "frequency conversion tests: from Hourly Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_H = dWrap(Interval(freq='H', year=2007, month=1, day=1, hour=0))
            date_H_end_of_year = dWrap(Interval(freq='H', year=2007, month=12, day=31,
                                      hour=23))
            date_H_end_of_quarter = dWrap(Interval(freq='H', year=2007, month=3, day=31,
                                         hour=23))
            date_H_end_of_month = dWrap(Interval(freq='H', year=2007, month=1, day=31,
                                       hour=23))
            date_H_end_of_week = dWrap(Interval(freq='H', year=2007, month=1, day=7,
                                      hour=23))
            date_H_end_of_day = dWrap(Interval(freq='H', year=2007, month=1, day=1,
                                     hour=23))
            date_H_end_of_bus = dWrap(Interval(freq='H', year=2007, month=1, day=1,
                                     hour=23))

            date_H_to_A = dWrap(Interval(freq='A', year=2007))
            date_H_to_Q = dWrap(Interval(freq='Q', year=2007, quarter=1))
            date_H_to_M = dWrap(Interval(freq='M', year=2007, month=1))
            date_H_to_W = dWrap(Interval(freq='WK', year=2007, month=1, day=7))
            date_H_to_D = dWrap(Interval(freq='D', year=2007, month=1, day=1))
            date_H_to_B = dWrap(Interval(freq='B', year=2007, month=1, day=1))

            date_H_to_T_start = dWrap(Interval(freq='Min', year=2007, month=1, day=1,
                                      hour=0, minute=0))
            date_H_to_T_end = dWrap(Interval(freq='Min', year=2007, month=1, day=1,
                                     hour=0, minute=59))
            date_H_to_S_start = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_H_to_S_end = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                     hour=0, minute=59, second=59))

            assert_func(date_H.asfreq('A'), date_H_to_A)
            assert_func(date_H_end_of_year.asfreq('A'), date_H_to_A)
            assert_func(date_H.asfreq('Q'), date_H_to_Q)
            assert_func(date_H_end_of_quarter.asfreq('Q'), date_H_to_Q)
            assert_func(date_H.asfreq('M'), date_H_to_M)
            assert_func(date_H_end_of_month.asfreq('M'), date_H_to_M)
            assert_func(date_H.asfreq('WK'), date_H_to_W)
            assert_func(date_H_end_of_week.asfreq('WK'), date_H_to_W)
            assert_func(date_H.asfreq('D'), date_H_to_D)
            assert_func(date_H_end_of_day.asfreq('D'), date_H_to_D)
            assert_func(date_H.asfreq('B'), date_H_to_B)
            assert_func(date_H_end_of_bus.asfreq('B'), date_H_to_B)

            assert_func(date_H.asfreq('Min', 'S'), date_H_to_T_start)
            assert_func(date_H.asfreq('Min', 'E'), date_H_to_T_end)
            assert_func(date_H.asfreq('S', 'S'), date_H_to_S_start)
            assert_func(date_H.asfreq('S', 'E'), date_H_to_S_end)

            assert_func(date_H.asfreq('H'), date_H)

    def test_conv_minutely(self):
        "frequency conversion tests: from Minutely Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_T = dWrap(Interval(freq='Min', year=2007, month=1, day=1,
                          hour=0, minute=0))
            date_T_end_of_year = dWrap(Interval(freq='Min', year=2007, month=12, day=31,
                                      hour=23, minute=59))
            date_T_end_of_quarter = dWrap(Interval(freq='Min', year=2007, month=3, day=31,
                                         hour=23, minute=59))
            date_T_end_of_month = dWrap(Interval(freq='Min', year=2007, month=1, day=31,
                                       hour=23, minute=59))
            date_T_end_of_week = dWrap(Interval(freq='Min', year=2007, month=1, day=7,
                                      hour=23, minute=59))
            date_T_end_of_day = dWrap(Interval(freq='Min', year=2007, month=1, day=1,
                                     hour=23, minute=59))
            date_T_end_of_bus = dWrap(Interval(freq='Min', year=2007, month=1, day=1,
                                     hour=23, minute=59))
            date_T_end_of_hour = dWrap(Interval(freq='Min', year=2007, month=1, day=1,
                                      hour=0, minute=59))

            date_T_to_A = dWrap(Interval(freq='A', year=2007))
            date_T_to_Q = dWrap(Interval(freq='Q', year=2007, quarter=1))
            date_T_to_M = dWrap(Interval(freq='M', year=2007, month=1))
            date_T_to_W = dWrap(Interval(freq='WK', year=2007, month=1, day=7))
            date_T_to_D = dWrap(Interval(freq='D', year=2007, month=1, day=1))
            date_T_to_B = dWrap(Interval(freq='B', year=2007, month=1, day=1))
            date_T_to_H = dWrap(Interval(freq='H', year=2007, month=1, day=1, hour=0))

            date_T_to_S_start = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=0, second=0))
            date_T_to_S_end = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                     hour=0, minute=0, second=59))

            assert_func(date_T.asfreq('A'), date_T_to_A)
            assert_func(date_T_end_of_year.asfreq('A'), date_T_to_A)
            assert_func(date_T.asfreq('Q'), date_T_to_Q)
            assert_func(date_T_end_of_quarter.asfreq('Q'), date_T_to_Q)
            assert_func(date_T.asfreq('M'), date_T_to_M)
            assert_func(date_T_end_of_month.asfreq('M'), date_T_to_M)
            assert_func(date_T.asfreq('WK'), date_T_to_W)
            assert_func(date_T_end_of_week.asfreq('WK'), date_T_to_W)
            assert_func(date_T.asfreq('D'), date_T_to_D)
            assert_func(date_T_end_of_day.asfreq('D'), date_T_to_D)
            assert_func(date_T.asfreq('B'), date_T_to_B)
            assert_func(date_T_end_of_bus.asfreq('B'), date_T_to_B)
            assert_func(date_T.asfreq('H'), date_T_to_H)
            assert_func(date_T_end_of_hour.asfreq('H'), date_T_to_H)

            assert_func(date_T.asfreq('S', 'S'), date_T_to_S_start)
            assert_func(date_T.asfreq('S', 'E'), date_T_to_S_end)

            assert_func(date_T.asfreq('Min'), date_T)

    def test_conv_secondly(self):
        "frequency conversion tests: from Secondly Frequency"

        for dWrap, assert_func in self.dateWrap:
            date_S = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                          hour=0, minute=0, second=0))
            date_S_end_of_year = dWrap(Interval(freq='S', year=2007, month=12, day=31,
                                      hour=23, minute=59, second=59))
            date_S_end_of_quarter = dWrap(Interval(freq='S', year=2007, month=3, day=31,
                                         hour=23, minute=59, second=59))
            date_S_end_of_month = dWrap(Interval(freq='S', year=2007, month=1, day=31,
                                       hour=23, minute=59, second=59))
            date_S_end_of_week = dWrap(Interval(freq='S', year=2007, month=1, day=7,
                                      hour=23, minute=59, second=59))
            date_S_end_of_day = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                     hour=23, minute=59, second=59))
            date_S_end_of_bus = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                     hour=23, minute=59, second=59))
            date_S_end_of_hour = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                      hour=0, minute=59, second=59))
            date_S_end_of_minute = dWrap(Interval(freq='S', year=2007, month=1, day=1,
                                        hour=0, minute=0, second=59))

            date_S_to_A = dWrap(Interval(freq='A', year=2007))
            date_S_to_Q = dWrap(Interval(freq='Q', year=2007, quarter=1))
            date_S_to_M = dWrap(Interval(freq='M', year=2007, month=1))
            date_S_to_W = dWrap(Interval(freq='WK', year=2007, month=1, day=7))
            date_S_to_D = dWrap(Interval(freq='D', year=2007, month=1, day=1))
            date_S_to_B = dWrap(Interval(freq='B', year=2007, month=1, day=1))
            date_S_to_H = dWrap(Interval(freq='H', year=2007, month=1, day=1,
                               hour=0))
            date_S_to_T = dWrap(Interval(freq='Min', year=2007, month=1, day=1,
                               hour=0, minute=0))

            assert_func(date_S.asfreq('A'), date_S_to_A)
            assert_func(date_S_end_of_year.asfreq('A'), date_S_to_A)
            assert_func(date_S.asfreq('Q'), date_S_to_Q)
            assert_func(date_S_end_of_quarter.asfreq('Q'), date_S_to_Q)
            assert_func(date_S.asfreq('M'), date_S_to_M)
            assert_func(date_S_end_of_month.asfreq('M'), date_S_to_M)
            assert_func(date_S.asfreq('WK'), date_S_to_W)
            assert_func(date_S_end_of_week.asfreq('WK'), date_S_to_W)
            assert_func(date_S.asfreq('D'), date_S_to_D)
            assert_func(date_S_end_of_day.asfreq('D'), date_S_to_D)
            assert_func(date_S.asfreq('B'), date_S_to_B)
            assert_func(date_S_end_of_bus.asfreq('B'), date_S_to_B)
            assert_func(date_S.asfreq('H'), date_S_to_H)
            assert_func(date_S_end_of_hour.asfreq('H'), date_S_to_H)
            assert_func(date_S.asfreq('Min'), date_S_to_T)
            assert_func(date_S_end_of_minute.asfreq('Min'), date_S_to_T)

            assert_func(date_S.asfreq('S'), date_S)


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
