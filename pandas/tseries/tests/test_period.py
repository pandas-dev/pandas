"""Tests suite for Period handling.

Parts derived from scikits.timeseries code, original authors:
- Pierre Gerard-Marchant & Matt Knox
- pierregm_at_uga_dot_edu - mattknow_ca_at_hotmail_dot_com

"""

from datetime import datetime, date, timedelta

from numpy.ma.testutils import assert_equal

from pandas import Timestamp
from pandas.tseries.frequencies import MONTHS, DAYS, _period_code_map
from pandas.tseries.period import Period, PeriodIndex, period_range
from pandas.tseries.index import DatetimeIndex, date_range, Index
from pandas.tseries.tools import to_datetime
import pandas.tseries.period as period
import pandas.tseries.offsets as offsets

import pandas.core.datetools as datetools
import pandas as pd
import numpy as np
from numpy.random import randn
from pandas.compat import range, lrange, lmap, zip

from pandas import Series, DataFrame, _np_version_under1p9
from pandas import tslib
from pandas.util.testing import(assert_series_equal, assert_almost_equal,
                                assertRaisesRegexp)
import pandas.util.testing as tm
from pandas import compat


class TestPeriodProperties(tm.TestCase):
    "Test properties such as year, month, weekday, etc...."
    #

    def test_quarterly_negative_ordinals(self):
        p = Period(ordinal=-1, freq='Q-DEC')
        self.assertEqual(p.year, 1969)
        self.assertEqual(p.quarter, 4)

        p = Period(ordinal=-2, freq='Q-DEC')
        self.assertEqual(p.year, 1969)
        self.assertEqual(p.quarter, 3)

        p = Period(ordinal=-2, freq='M')
        self.assertEqual(p.year, 1969)
        self.assertEqual(p.month, 11)

    def test_period_cons_quarterly(self):
        # bugs in scikits.timeseries
        for month in MONTHS:
            freq = 'Q-%s' % month
            exp = Period('1989Q3', freq=freq)
            self.assertIn('1989Q3', str(exp))
            stamp = exp.to_timestamp('D', how='end')
            p = Period(stamp, freq=freq)
            self.assertEqual(p, exp)

            stamp = exp.to_timestamp('3D', how='end')
            p = Period(stamp, freq=freq)
            self.assertEqual(p, exp)

    def test_period_cons_annual(self):
        # bugs in scikits.timeseries
        for month in MONTHS:
            freq = 'A-%s' % month
            exp = Period('1989', freq=freq)
            stamp = exp.to_timestamp('D', how='end') + timedelta(days=30)
            p = Period(stamp, freq=freq)
            self.assertEqual(p, exp + 1)

    def test_period_cons_weekly(self):
        for num in range(10, 17):
            daystr = '2011-02-%d' % num
            for day in DAYS:
                freq = 'W-%s' % day

                result = Period(daystr, freq=freq)
                expected = Period(daystr, freq='D').asfreq(freq)
                self.assertEqual(result, expected)

    def test_period_cons_nat(self):
        p = Period('NaT', freq='M')
        self.assertEqual(p.ordinal, tslib.iNaT)
        self.assertEqual(p.freq, 'M')
        self.assertEqual((p + 1).ordinal, tslib.iNaT)

        p = Period('nat', freq='W-SUN')
        self.assertEqual(p.ordinal, tslib.iNaT)
        self.assertEqual(p.freq, 'W-SUN')
        self.assertEqual((p + 1).ordinal, tslib.iNaT)

        p = Period(tslib.iNaT, freq='D')
        self.assertEqual(p.ordinal, tslib.iNaT)
        self.assertEqual(p.freq, 'D')
        self.assertEqual((p + 1).ordinal, tslib.iNaT)

        p = Period(tslib.iNaT, freq='3D')
        self.assertEqual(p.ordinal, tslib.iNaT)
        self.assertEqual(p.freq, offsets.Day(3))
        self.assertEqual(p.freqstr, '3D')
        self.assertEqual((p + 1).ordinal, tslib.iNaT)

        self.assertRaises(ValueError, Period, 'NaT')

    def test_period_cons_mult(self):
        p1 = Period('2011-01', freq='3M')
        p2 = Period('2011-01', freq='M')
        self.assertEqual(p1.ordinal, p2.ordinal)

        self.assertEqual(p1.freq, offsets.MonthEnd(3))
        self.assertEqual(p1.freqstr, '3M')

        self.assertEqual(p2.freq, offsets.MonthEnd())
        self.assertEqual(p2.freqstr, 'M')

        result = p1 + 1
        self.assertEqual(result.ordinal, (p2 + 3).ordinal)
        self.assertEqual(result.freq, p1.freq)
        self.assertEqual(result.freqstr, '3M')

        result = p1 - 1
        self.assertEqual(result.ordinal, (p2 - 3).ordinal)
        self.assertEqual(result.freq, p1.freq)
        self.assertEqual(result.freqstr, '3M')

        msg = ('Frequency must be positive, because it'
               ' represents span: -3M')
        with tm.assertRaisesRegexp(ValueError, msg):
            Period('2011-01', freq='-3M')

        msg = ('Frequency must be positive, because it'
               ' represents span: 0M')
        with tm.assertRaisesRegexp(ValueError, msg):
            Period('2011-01', freq='0M')

    def test_timestamp_tz_arg(self):
        tm._skip_if_no_pytz()
        import pytz
        for case in ['Europe/Brussels', 'Asia/Tokyo', 'US/Pacific']:
            p = Period('1/1/2005', freq='M').to_timestamp(tz=case)
            exp = Timestamp('1/1/2005', tz='UTC').tz_convert(case)
            exp_zone = pytz.timezone(case).normalize(p)

            self.assertEqual(p, exp)
            self.assertEqual(p.tz, exp_zone.tzinfo)
            self.assertEqual(p.tz, exp.tz)

            p = Period('1/1/2005', freq='3H').to_timestamp(tz=case)
            exp = Timestamp('1/1/2005', tz='UTC').tz_convert(case)
            exp_zone = pytz.timezone(case).normalize(p)

            self.assertEqual(p, exp)
            self.assertEqual(p.tz, exp_zone.tzinfo)
            self.assertEqual(p.tz, exp.tz)

            p = Period('1/1/2005', freq='A').to_timestamp(freq='A', tz=case)
            exp = Timestamp('31/12/2005', tz='UTC').tz_convert(case)
            exp_zone = pytz.timezone(case).normalize(p)

            self.assertEqual(p, exp)
            self.assertEqual(p.tz, exp_zone.tzinfo)
            self.assertEqual(p.tz, exp.tz)

            p = Period('1/1/2005', freq='A').to_timestamp(freq='3H', tz=case)
            exp = Timestamp('1/1/2005', tz='UTC').tz_convert(case)
            exp_zone = pytz.timezone(case).normalize(p)

            self.assertEqual(p, exp)
            self.assertEqual(p.tz, exp_zone.tzinfo)
            self.assertEqual(p.tz, exp.tz)

    def test_timestamp_tz_arg_dateutil(self):
        from pandas.tslib import _dateutil_gettz as gettz
        from pandas.tslib import maybe_get_tz
        for case in ['dateutil/Europe/Brussels', 'dateutil/Asia/Tokyo',
                     'dateutil/US/Pacific']:
            p = Period('1/1/2005', freq='M').to_timestamp(tz=maybe_get_tz(case))
            exp = Timestamp('1/1/2005', tz='UTC').tz_convert(case)
            self.assertEqual(p, exp)
            self.assertEqual(p.tz, gettz(case.split('/', 1)[1]))
            self.assertEqual(p.tz, exp.tz)

            p = Period('1/1/2005', freq='M').to_timestamp(freq='3H', tz=maybe_get_tz(case))
            exp = Timestamp('1/1/2005', tz='UTC').tz_convert(case)
            self.assertEqual(p, exp)
            self.assertEqual(p.tz, gettz(case.split('/', 1)[1]))
            self.assertEqual(p.tz, exp.tz)

    def test_timestamp_tz_arg_dateutil_from_string(self):
        from pandas.tslib import _dateutil_gettz as gettz
        p = Period('1/1/2005', freq='M').to_timestamp(tz='dateutil/Europe/Brussels')
        self.assertEqual(p.tz, gettz('Europe/Brussels'))

    def test_timestamp_nat_tz(self):
        t = Period('NaT', freq='M').to_timestamp()
        self.assertTrue(t is tslib.NaT)

        t = Period('NaT', freq='M').to_timestamp(tz='Asia/Tokyo')
        self.assertTrue(t is tslib.NaT)

    def test_timestamp_mult(self):
        p = pd.Period('2011-01', freq='M')
        self.assertEqual(p.to_timestamp(how='S'), pd.Timestamp('2011-01-01'))
        self.assertEqual(p.to_timestamp(how='E'), pd.Timestamp('2011-01-31'))

        p = pd.Period('2011-01', freq='3M')
        self.assertEqual(p.to_timestamp(how='S'), pd.Timestamp('2011-01-01'))
        self.assertEqual(p.to_timestamp(how='E'), pd.Timestamp('2011-03-31'))

    def test_timestamp_nat_mult(self):
        for freq in ['M', '3M']:
            p = pd.Period('NaT', freq=freq)
            self.assertTrue(p.to_timestamp(how='S') is pd.NaT)
            self.assertTrue(p.to_timestamp(how='E') is pd.NaT)

    def test_period_constructor(self):
        i1 = Period('1/1/2005', freq='M')
        i2 = Period('Jan 2005')

        self.assertEqual(i1, i2)

        i1 = Period('2005', freq='A')
        i2 = Period('2005')
        i3 = Period('2005', freq='a')

        self.assertEqual(i1, i2)
        self.assertEqual(i1, i3)

        i4 = Period('2005', freq='M')
        i5 = Period('2005', freq='m')

        self.assertRaises(ValueError, i1.__ne__, i4)
        self.assertEqual(i4, i5)

        i1 = Period.now('Q')
        i2 = Period(datetime.now(), freq='Q')
        i3 = Period.now('q')

        self.assertEqual(i1, i2)
        self.assertEqual(i1, i3)

        # Biz day construction, roll forward if non-weekday
        i1 = Period('3/10/12', freq='B')
        i2 = Period('3/10/12', freq='D')
        self.assertEqual(i1, i2.asfreq('B'))
        i2 = Period('3/11/12', freq='D')
        self.assertEqual(i1, i2.asfreq('B'))
        i2 = Period('3/12/12', freq='D')
        self.assertEqual(i1, i2.asfreq('B'))

        i3 = Period('3/10/12', freq='b')
        self.assertEqual(i1, i3)

        i1 = Period(year=2005, quarter=1, freq='Q')
        i2 = Period('1/1/2005', freq='Q')
        self.assertEqual(i1, i2)

        i1 = Period(year=2005, quarter=3, freq='Q')
        i2 = Period('9/1/2005', freq='Q')
        self.assertEqual(i1, i2)

        i1 = Period(year=2005, month=3, day=1, freq='D')
        i2 = Period('3/1/2005', freq='D')
        self.assertEqual(i1, i2)

        i3 = Period(year=2005, month=3, day=1, freq='d')
        self.assertEqual(i1, i3)

        i1 = Period(year=2012, month=3, day=10, freq='B')
        i2 = Period('3/12/12', freq='B')
        self.assertEqual(i1, i2)

        i1 = Period('2005Q1')
        i2 = Period(year=2005, quarter=1, freq='Q')
        i3 = Period('2005q1')
        self.assertEqual(i1, i2)
        self.assertEqual(i1, i3)

        i1 = Period('05Q1')
        self.assertEqual(i1, i2)
        lower = Period('05q1')
        self.assertEqual(i1, lower)

        i1 = Period('1Q2005')
        self.assertEqual(i1, i2)
        lower = Period('1q2005')
        self.assertEqual(i1, lower)

        i1 = Period('1Q05')
        self.assertEqual(i1, i2)
        lower = Period('1q05')
        self.assertEqual(i1, lower)

        i1 = Period('4Q1984')
        self.assertEqual(i1.year, 1984)
        lower = Period('4q1984')
        self.assertEqual(i1, lower)

        i1 = Period('1982', freq='min')
        i2 = Period('1982', freq='MIN')
        self.assertEqual(i1, i2)
        i2 = Period('1982', freq=('Min', 1))
        self.assertEqual(i1, i2)

        expected = Period('2007-01', freq='M')
        i1 = Period('200701', freq='M')
        self.assertEqual(i1, expected)

        i1 = Period('200701', freq='M')
        self.assertEqual(i1, expected)

        i1 = Period(200701, freq='M')
        self.assertEqual(i1, expected)

        i1 = Period(ordinal=200701, freq='M')
        self.assertEqual(i1.year, 18695)

        i1 = Period(datetime(2007, 1, 1), freq='M')
        i2 = Period('200701', freq='M')
        self.assertEqual(i1, i2)

        i1 = Period(date(2007, 1, 1), freq='M')
        i2 = Period(datetime(2007, 1, 1), freq='M')
        i3 = Period(np.datetime64('2007-01-01'), freq='M')
        i4 = Period(np.datetime64('2007-01-01 00:00:00Z'), freq='M')
        i5 = Period(np.datetime64('2007-01-01 00:00:00.000Z'), freq='M')
        self.assertEqual(i1, i2)
        self.assertEqual(i1, i3)
        self.assertEqual(i1, i4)
        self.assertEqual(i1, i5)

        i1 = Period('2007-01-01 09:00:00.001')
        expected = Period(datetime(2007, 1, 1, 9, 0, 0, 1000), freq='L')
        self.assertEqual(i1, expected)

        expected = Period(np.datetime64('2007-01-01 09:00:00.001Z'), freq='L')
        self.assertEqual(i1, expected)

        i1 = Period('2007-01-01 09:00:00.00101')
        expected = Period(datetime(2007, 1, 1, 9, 0, 0, 1010), freq='U')
        self.assertEqual(i1, expected)

        expected = Period(np.datetime64('2007-01-01 09:00:00.00101Z'),
                          freq='U')
        self.assertEqual(i1, expected)

        self.assertRaises(ValueError, Period, ordinal=200701)

        self.assertRaises(ValueError, Period, '2007-1-1', freq='X')


    def test_period_constructor_offsets(self):
        self.assertEqual(Period('1/1/2005', freq=offsets.MonthEnd()),
                         Period('1/1/2005', freq='M'))
        self.assertEqual(Period('2005', freq=offsets.YearEnd()),
                         Period('2005', freq='A'))
        self.assertEqual(Period('2005', freq=offsets.MonthEnd()),
                         Period('2005', freq='M'))
        self.assertEqual(Period('3/10/12', freq=offsets.BusinessDay()),
                         Period('3/10/12', freq='B'))
        self.assertEqual(Period('3/10/12', freq=offsets.Day()),
                         Period('3/10/12', freq='D'))

        self.assertEqual(Period(year=2005, quarter=1,
                                freq=offsets.QuarterEnd(startingMonth=12)),
                         Period(year=2005, quarter=1, freq='Q'))
        self.assertEqual(Period(year=2005, quarter=2,
                                freq=offsets.QuarterEnd(startingMonth=12)),
                         Period(year=2005, quarter=2, freq='Q'))

        self.assertEqual(Period(year=2005, month=3, day=1, freq=offsets.Day()),
                         Period(year=2005, month=3, day=1, freq='D'))
        self.assertEqual(Period(year=2012, month=3, day=10, freq=offsets.BDay()),
                         Period(year=2012, month=3, day=10, freq='B'))

        expected = Period('2005-03-01', freq='3D')
        self.assertEqual(Period(year=2005, month=3, day=1, freq=offsets.Day(3)),
                         expected)
        self.assertEqual(Period(year=2005, month=3, day=1, freq='3D'),
                         expected)

        self.assertEqual(Period(year=2012, month=3, day=10, freq=offsets.BDay(3)),
                         Period(year=2012, month=3, day=10, freq='3B'))

        self.assertEqual(Period(200701, freq=offsets.MonthEnd()),
                         Period(200701, freq='M'))

        i1 = Period(ordinal=200701, freq=offsets.MonthEnd())
        i2 = Period(ordinal=200701, freq='M')
        self.assertEqual(i1, i2)
        self.assertEqual(i1.year, 18695)
        self.assertEqual(i2.year, 18695)

        i1 = Period(datetime(2007, 1, 1), freq='M')
        i2 = Period('200701', freq='M')
        self.assertEqual(i1, i2)

        i1 = Period(date(2007, 1, 1), freq='M')
        i2 = Period(datetime(2007, 1, 1), freq='M')
        i3 = Period(np.datetime64('2007-01-01'), freq='M')
        i4 = Period(np.datetime64('2007-01-01 00:00:00Z'), freq='M')
        i5 = Period(np.datetime64('2007-01-01 00:00:00.000Z'), freq='M')
        self.assertEqual(i1, i2)
        self.assertEqual(i1, i3)
        self.assertEqual(i1, i4)
        self.assertEqual(i1, i5)

        i1 = Period('2007-01-01 09:00:00.001')
        expected = Period(datetime(2007, 1, 1, 9, 0, 0, 1000), freq='L')
        self.assertEqual(i1, expected)

        expected = Period(np.datetime64('2007-01-01 09:00:00.001Z'), freq='L')
        self.assertEqual(i1, expected)

        i1 = Period('2007-01-01 09:00:00.00101')
        expected = Period(datetime(2007, 1, 1, 9, 0, 0, 1010), freq='U')
        self.assertEqual(i1, expected)

        expected = Period(np.datetime64('2007-01-01 09:00:00.00101Z'),
                          freq='U')
        self.assertEqual(i1, expected)

        self.assertRaises(ValueError, Period, ordinal=200701)

        self.assertRaises(ValueError, Period, '2007-1-1', freq='X')


    def test_freq_str(self):
        i1 = Period('1982', freq='Min')
        self.assertEqual(i1.freq, offsets.Minute())
        self.assertEqual(i1.freqstr, 'T')

    def test_repr(self):
        p = Period('Jan-2000')
        self.assertIn('2000-01', repr(p))

        p = Period('2000-12-15')
        self.assertIn('2000-12-15', repr(p))

    def test_repr_nat(self):
        p = Period('nat', freq='M')
        self.assertIn(repr(tslib.NaT), repr(p))

    def test_millisecond_repr(self):
        p = Period('2000-01-01 12:15:02.123')

        self.assertEqual("Period('2000-01-01 12:15:02.123', 'L')", repr(p))

    def test_microsecond_repr(self):
        p = Period('2000-01-01 12:15:02.123567')

        self.assertEqual("Period('2000-01-01 12:15:02.123567', 'U')", repr(p))

    def test_strftime(self):
        p = Period('2000-1-1 12:34:12', freq='S')
        res = p.strftime('%Y-%m-%d %H:%M:%S')
        self.assertEqual(res,  '2000-01-01 12:34:12')
        tm.assertIsInstance(res, compat.text_type) # GH3363

    def test_sub_delta(self):
        left, right = Period('2011', freq='A'), Period('2007', freq='A')
        result = left - right
        self.assertEqual(result, 4)

        self.assertRaises(ValueError, left.__sub__,
                          Period('2007-01', freq='M'))

    def test_to_timestamp(self):
        p = Period('1982', freq='A')
        start_ts = p.to_timestamp(how='S')
        aliases = ['s', 'StarT', 'BEGIn']
        for a in aliases:
            self.assertEqual(start_ts, p.to_timestamp('D', how=a))
            # freq with mult should not affect to the result
            self.assertEqual(start_ts, p.to_timestamp('3D', how=a))

        end_ts = p.to_timestamp(how='E')
        aliases = ['e', 'end', 'FINIsH']
        for a in aliases:
            self.assertEqual(end_ts, p.to_timestamp('D', how=a))
            self.assertEqual(end_ts, p.to_timestamp('3D', how=a))

        from_lst = ['A', 'Q', 'M', 'W', 'B',
                    'D', 'H', 'Min', 'S']

        def _ex(p):
            return Timestamp((p + 1).start_time.value - 1)

        for i, fcode in enumerate(from_lst):
            p = Period('1982', freq=fcode)
            result = p.to_timestamp().to_period(fcode)
            self.assertEqual(result, p)

            self.assertEqual(p.start_time, p.to_timestamp(how='S'))

            self.assertEqual(p.end_time, _ex(p))

        # Frequency other than daily

        p = Period('1985', freq='A')

        result = p.to_timestamp('H', how='end')
        expected = datetime(1985, 12, 31, 23)
        self.assertEqual(result, expected)
        result = p.to_timestamp('3H', how='end')
        self.assertEqual(result, expected)

        result = p.to_timestamp('T', how='end')
        expected = datetime(1985, 12, 31, 23, 59)
        self.assertEqual(result, expected)
        result = p.to_timestamp('2T', how='end')
        self.assertEqual(result, expected)


        result = p.to_timestamp(how='end')
        expected = datetime(1985, 12, 31)
        self.assertEqual(result, expected)

        expected = datetime(1985, 1, 1)
        result = p.to_timestamp('H', how='start')
        self.assertEqual(result, expected)
        result = p.to_timestamp('T', how='start')
        self.assertEqual(result, expected)
        result = p.to_timestamp('S', how='start')
        self.assertEqual(result, expected)
        result = p.to_timestamp('3H', how='start')
        self.assertEqual(result, expected)
        result = p.to_timestamp('5S', how='start')
        self.assertEqual(result, expected)

        p = Period('NaT', freq='W')
        self.assertTrue(p.to_timestamp() is tslib.NaT)

    def test_start_time(self):
        freq_lst = ['A', 'Q', 'M', 'D', 'H', 'T', 'S']
        xp = datetime(2012, 1, 1)
        for f in freq_lst:
            p = Period('2012', freq=f)
            self.assertEqual(p.start_time, xp)
        self.assertEqual(Period('2012', freq='B').start_time,
                         datetime(2012, 1, 2))
        self.assertEqual(Period('2012', freq='W').start_time,
                         datetime(2011, 12, 26))

        p = Period('NaT', freq='W')
        self.assertTrue(p.start_time is tslib.NaT)

    def test_end_time(self):
        p = Period('2012', freq='A')

        def _ex(*args):
            return Timestamp(Timestamp(datetime(*args)).value - 1)

        xp = _ex(2013, 1, 1)
        self.assertEqual(xp, p.end_time)

        p = Period('2012', freq='Q')
        xp = _ex(2012, 4, 1)
        self.assertEqual(xp, p.end_time)

        p = Period('2012', freq='M')
        xp = _ex(2012, 2, 1)
        self.assertEqual(xp, p.end_time)

        xp = _ex(2012, 1, 2)
        p = Period('2012', freq='D')
        self.assertEqual(p.end_time, xp)

        xp = _ex(2012, 1, 1, 1)
        p = Period('2012', freq='H')
        self.assertEqual(p.end_time, xp)

        xp = _ex(2012, 1, 3)
        self.assertEqual(Period('2012', freq='B').end_time, xp)

        xp = _ex(2012, 1, 2)
        self.assertEqual(Period('2012', freq='W').end_time, xp)

        p = Period('NaT', freq='W')
        self.assertTrue(p.end_time is tslib.NaT)

    def test_anchor_week_end_time(self):
        def _ex(*args):
            return Timestamp(Timestamp(datetime(*args)).value - 1)

        p = Period('2013-1-1', 'W-SAT')
        xp = _ex(2013, 1, 6)
        self.assertEqual(p.end_time, xp)

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
        w_date = Period(freq='W', year=2007, month=1, day=7)
        #
        assert_equal(w_date.year, 2007)
        assert_equal(w_date.quarter, 1)
        assert_equal(w_date.month, 1)
        assert_equal(w_date.week, 1)
        assert_equal((w_date - 1).week, 52)
        assert_equal(w_date.days_in_month, 31)
        assert_equal(Period(freq='W', year=2012, month=2, day=1).days_in_month, 29)

    def test_properties_weekly_legacy(self):
        # Test properties on Periods with daily frequency.
        with tm.assert_produces_warning(FutureWarning):
            w_date = Period(freq='WK', year=2007, month=1, day=7)
        #
        assert_equal(w_date.year, 2007)
        assert_equal(w_date.quarter, 1)
        assert_equal(w_date.month, 1)
        assert_equal(w_date.week, 1)
        assert_equal((w_date - 1).week, 52)
        assert_equal(w_date.days_in_month, 31)
        with tm.assert_produces_warning(FutureWarning):
            exp = Period(freq='WK', year=2012, month=2, day=1)
        assert_equal(exp.days_in_month, 29)

    def test_properties_daily(self):
        # Test properties on Periods with daily frequency.
        b_date = Period(freq='B', year=2007, month=1, day=1)
        #
        assert_equal(b_date.year, 2007)
        assert_equal(b_date.quarter, 1)
        assert_equal(b_date.month, 1)
        assert_equal(b_date.day, 1)
        assert_equal(b_date.weekday, 0)
        assert_equal(b_date.dayofyear, 1)
        assert_equal(b_date.days_in_month, 31)
        assert_equal(Period(freq='B', year=2012, month=2, day=1).days_in_month, 29)
        #
        d_date = Period(freq='D', year=2007, month=1, day=1)
        #
        assert_equal(d_date.year, 2007)
        assert_equal(d_date.quarter, 1)
        assert_equal(d_date.month, 1)
        assert_equal(d_date.day, 1)
        assert_equal(d_date.weekday, 0)
        assert_equal(d_date.dayofyear, 1)
        assert_equal(d_date.days_in_month, 31)
        assert_equal(Period(freq='D', year=2012, month=2,
                            day=1).days_in_month, 29)

    def test_properties_hourly(self):
        # Test properties on Periods with hourly frequency.
        h_date1 = Period(freq='H', year=2007, month=1, day=1, hour=0)
        h_date2 = Period(freq='2H', year=2007, month=1, day=1, hour=0)

        for h_date in [h_date1, h_date2]:
            assert_equal(h_date.year, 2007)
            assert_equal(h_date.quarter, 1)
            assert_equal(h_date.month, 1)
            assert_equal(h_date.day, 1)
            assert_equal(h_date.weekday, 0)
            assert_equal(h_date.dayofyear, 1)
            assert_equal(h_date.hour, 0)
            assert_equal(h_date.days_in_month, 31)
            assert_equal(Period(freq='H', year=2012, month=2, day=1,
                                hour=0).days_in_month, 29)

    def test_properties_minutely(self):
        # Test properties on Periods with minutely frequency.
        t_date = Period(freq='Min', year=2007, month=1, day=1, hour=0,
                        minute=0)
        #
        assert_equal(t_date.quarter, 1)
        assert_equal(t_date.month, 1)
        assert_equal(t_date.day, 1)
        assert_equal(t_date.weekday, 0)
        assert_equal(t_date.dayofyear, 1)
        assert_equal(t_date.hour, 0)
        assert_equal(t_date.minute, 0)
        assert_equal(t_date.days_in_month, 31)
        assert_equal(Period(freq='D', year=2012, month=2, day=1, hour=0,
                            minute=0).days_in_month, 29)

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
        assert_equal(s_date.dayofyear, 1)
        assert_equal(s_date.hour, 0)
        assert_equal(s_date.minute, 0)
        assert_equal(s_date.second, 0)
        assert_equal(s_date.days_in_month, 31)
        assert_equal(Period(freq='Min', year=2012, month=2, day=1, hour=0,
                            minute=0, second=0).days_in_month, 29)

    def test_properties_nat(self):
        p_nat = Period('NaT', freq='M')
        t_nat = pd.Timestamp('NaT')
        # confirm Period('NaT') work identical with Timestamp('NaT')
        for f in ['year', 'month', 'day', 'hour', 'minute', 'second',
                  'week', 'dayofyear', 'quarter', 'days_in_month']:
            self.assertTrue(np.isnan(getattr(p_nat, f)))
            self.assertTrue(np.isnan(getattr(t_nat, f)))

        for f in ['weekofyear', 'dayofweek', 'weekday', 'qyear']:
            self.assertTrue(np.isnan(getattr(p_nat, f)))

    def test_pnow(self):
        dt = datetime.now()

        val = period.pnow('D')
        exp = Period(dt, freq='D')
        self.assertEqual(val, exp)

        val2 = period.pnow('2D')
        exp2 = Period(dt, freq='2D')
        self.assertEqual(val2, exp2)
        self.assertEqual(val.ordinal, val2.ordinal)
        self.assertEqual(val.ordinal, exp2.ordinal)

    def test_constructor_corner(self):
        expected = Period('2007-01', freq='2M')
        self.assertEqual(Period(year=2007, month=1, freq='2M'), expected)

        self.assertRaises(ValueError, Period, datetime.now())
        self.assertRaises(ValueError, Period, datetime.now().date())
        self.assertRaises(ValueError, Period, 1.6, freq='D')
        self.assertRaises(ValueError, Period, ordinal=1.6, freq='D')
        self.assertRaises(ValueError, Period, ordinal=2, value=1, freq='D')
        self.assertRaises(ValueError, Period)
        self.assertRaises(ValueError, Period, month=1)

        p = Period('2007-01-01', freq='D')

        result = Period(p, freq='A')
        exp = Period('2007', freq='A')
        self.assertEqual(result, exp)

    def test_constructor_infer_freq(self):
        p = Period('2007-01-01')
        self.assertEqual(p.freq, 'D')

        p = Period('2007-01-01 07')
        self.assertEqual(p.freq, 'H')

        p = Period('2007-01-01 07:10')
        self.assertEqual(p.freq, 'T')

        p = Period('2007-01-01 07:10:15')
        self.assertEqual(p.freq, 'S')

        p = Period('2007-01-01 07:10:15.123')
        self.assertEqual(p.freq, 'L')

        p = Period('2007-01-01 07:10:15.123000')
        self.assertEqual(p.freq, 'L')

        p = Period('2007-01-01 07:10:15.123400')
        self.assertEqual(p.freq, 'U')

    def test_asfreq_MS(self):
        initial = Period("2013")

        self.assertEqual(initial.asfreq(freq="M", how="S"), Period('2013-01', 'M'))
        self.assertRaises(ValueError, initial.asfreq, freq="MS", how="S")
        tm.assertRaisesRegexp(ValueError, "Unknown freqstr: MS", pd.Period, '2013-01', 'MS')
        self.assertTrue(_period_code_map.get("MS") is None)

def noWrap(item):
    return item


class TestFreqConversion(tm.TestCase):
    "Test frequency conversion of date objects"

    def test_asfreq_corner(self):
        val = Period(freq='A', year=2007)
        result1 = val.asfreq('5t')
        result2 = val.asfreq('t')
        expected = Period('2007-12-31 23:59', freq='t')
        self.assertEqual(result1.ordinal, expected.ordinal)
        self.assertEqual(result1.freqstr, '5T')
        self.assertEqual(result2.ordinal, expected.ordinal)
        self.assertEqual(result2.freqstr, 'T')

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
        ival_A_to_W_start = Period(freq='W', year=2007, month=1, day=1)
        ival_A_to_W_end = Period(freq='W', year=2007, month=12, day=31)
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
        assert_equal(ival_A.asfreq('W', 'S'), ival_A_to_W_start)
        assert_equal(ival_A.asfreq('W', 'E'), ival_A_to_W_end)
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
        ival_Q_to_W_start = Period(freq='W', year=2007, month=1, day=1)
        ival_Q_to_W_end = Period(freq='W', year=2007, month=3, day=31)
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
        assert_equal(ival_Q.asfreq('W', 'S'), ival_Q_to_W_start)
        assert_equal(ival_Q.asfreq('W', 'E'), ival_Q_to_W_end)
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
        ival_M_to_W_start = Period(freq='W', year=2007, month=1, day=1)
        ival_M_to_W_end = Period(freq='W', year=2007, month=1, day=31)
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

        assert_equal(ival_M.asfreq('W', 'S'), ival_M_to_W_start)
        assert_equal(ival_M.asfreq('W', 'E'), ival_M_to_W_end)
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
        ival_W = Period(freq='W', year=2007, month=1, day=1)

        ival_WSUN = Period(freq='W', year=2007, month=1, day=7)
        ival_WSAT = Period(freq='W-SAT', year=2007, month=1, day=6)
        ival_WFRI = Period(freq='W-FRI', year=2007, month=1, day=5)
        ival_WTHU = Period(freq='W-THU', year=2007, month=1, day=4)
        ival_WWED = Period(freq='W-WED', year=2007, month=1, day=3)
        ival_WTUE = Period(freq='W-TUE', year=2007, month=1, day=2)
        ival_WMON = Period(freq='W-MON', year=2007, month=1, day=1)

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

        ival_W_end_of_year = Period(freq='W', year=2007, month=12, day=31)
        ival_W_end_of_quarter = Period(freq='W', year=2007, month=3, day=31)
        ival_W_end_of_month = Period(freq='W', year=2007, month=1, day=31)
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

        assert_equal(ival_W.asfreq('W'), ival_W)

    def test_conv_weekly_legacy(self):
        # frequency conversion tests: from Weekly Frequency

        with tm.assert_produces_warning(FutureWarning):
            ival_W = Period(freq='WK', year=2007, month=1, day=1)

        with tm.assert_produces_warning(FutureWarning):
            ival_WSUN = Period(freq='WK', year=2007, month=1, day=7)
        with tm.assert_produces_warning(FutureWarning):
            ival_WSAT = Period(freq='WK-SAT', year=2007, month=1, day=6)
        with tm.assert_produces_warning(FutureWarning):
            ival_WFRI = Period(freq='WK-FRI', year=2007, month=1, day=5)
        with tm.assert_produces_warning(FutureWarning):
            ival_WTHU = Period(freq='WK-THU', year=2007, month=1, day=4)
        with tm.assert_produces_warning(FutureWarning):
            ival_WWED = Period(freq='WK-WED', year=2007, month=1, day=3)
        with tm.assert_produces_warning(FutureWarning):
            ival_WTUE = Period(freq='WK-TUE', year=2007, month=1, day=2)
        with tm.assert_produces_warning(FutureWarning):
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

        with tm.assert_produces_warning(FutureWarning):
            ival_W_end_of_year = Period(freq='WK', year=2007, month=12, day=31)
        with tm.assert_produces_warning(FutureWarning):
            ival_W_end_of_quarter = Period(freq='WK', year=2007, month=3, day=31)
        with tm.assert_produces_warning(FutureWarning):
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

        with tm.assert_produces_warning(FutureWarning):
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
        ival_B_to_W = Period(freq='W', year=2007, month=1, day=7)
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
        assert_equal(ival_B.asfreq('W'), ival_B_to_W)
        assert_equal(ival_B_end_of_week.asfreq('W'), ival_B_to_W)

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
        ival_D_to_W = Period(freq='W', year=2007, month=1, day=7)

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
        assert_equal(ival_D.asfreq('W'), ival_D_to_W)
        assert_equal(ival_D_end_of_week.asfreq('W'), ival_D_to_W)

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
        ival_H_to_W = Period(freq='W', year=2007, month=1, day=7)
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
        assert_equal(ival_H.asfreq('W'), ival_H_to_W)
        assert_equal(ival_H_end_of_week.asfreq('W'), ival_H_to_W)
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
        ival_T_to_W = Period(freq='W', year=2007, month=1, day=7)
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
        assert_equal(ival_T.asfreq('W'), ival_T_to_W)
        assert_equal(ival_T_end_of_week.asfreq('W'), ival_T_to_W)
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
        ival_S_to_W = Period(freq='W', year=2007, month=1, day=7)
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
        assert_equal(ival_S.asfreq('W'), ival_S_to_W)
        assert_equal(ival_S_end_of_week.asfreq('W'), ival_S_to_W)
        assert_equal(ival_S.asfreq('D'), ival_S_to_D)
        assert_equal(ival_S_end_of_day.asfreq('D'), ival_S_to_D)
        assert_equal(ival_S.asfreq('B'), ival_S_to_B)
        assert_equal(ival_S_end_of_bus.asfreq('B'), ival_S_to_B)
        assert_equal(ival_S.asfreq('H'), ival_S_to_H)
        assert_equal(ival_S_end_of_hour.asfreq('H'), ival_S_to_H)
        assert_equal(ival_S.asfreq('Min'), ival_S_to_T)
        assert_equal(ival_S_end_of_minute.asfreq('Min'), ival_S_to_T)

        assert_equal(ival_S.asfreq('S'), ival_S)

    def test_asfreq_nat(self):
        p = Period('NaT', freq='A')
        result = p.asfreq('M')
        self.assertEqual(result.ordinal, tslib.iNaT)
        self.assertEqual(result.freq, 'M')

    def test_asfreq_mult(self):
        # normal freq to mult freq
        p = Period(freq='A', year=2007)
        # ordinal will not change
        for freq in ['3A', offsets.YearEnd(3)]:
            result = p.asfreq(freq)
            expected = Period('2007', freq='3A')

            self.assertEqual(result, expected)
            self.assertEqual(result.ordinal, expected.ordinal)
            self.assertEqual(result.freq, expected.freq)
        # ordinal will not change
        for freq in ['3A', offsets.YearEnd(3)]:
            result = p.asfreq(freq, how='S')
            expected = Period('2007', freq='3A')

            self.assertEqual(result, expected)
            self.assertEqual(result.ordinal, expected.ordinal)
            self.assertEqual(result.freq, expected.freq)

        # mult freq to normal freq
        p = Period(freq='3A', year=2007)
        # ordinal will change because how=E is the default
        for freq in ['A', offsets.YearEnd()]:
            result = p.asfreq(freq)
            expected = Period('2009', freq='A')

            self.assertEqual(result, expected)
            self.assertEqual(result.ordinal, expected.ordinal)
            self.assertEqual(result.freq, expected.freq)
        # ordinal will not change
        for freq in ['A', offsets.YearEnd()]:
            result = p.asfreq(freq, how='S')
            expected = Period('2007', freq='A')

            self.assertEqual(result, expected)
            self.assertEqual(result.ordinal, expected.ordinal)
            self.assertEqual(result.freq, expected.freq)

        p = Period(freq='A', year=2007)
        for freq in ['2M', offsets.MonthEnd(2)]:
            result = p.asfreq(freq)
            expected = Period('2007-12', freq='2M')

            self.assertEqual(result, expected)
            self.assertEqual(result.ordinal, expected.ordinal)
            self.assertEqual(result.freq, expected.freq)
        for freq in ['2M', offsets.MonthEnd(2)]:
            result = p.asfreq(freq, how='S')
            expected = Period('2007-01', freq='2M')

            self.assertEqual(result, expected)
            self.assertEqual(result.ordinal, expected.ordinal)
            self.assertEqual(result.freq, expected.freq)

        p = Period(freq='3A', year=2007)
        for freq in ['2M', offsets.MonthEnd(2)]:
            result = p.asfreq(freq)
            expected = Period('2009-12', freq='2M')

            self.assertEqual(result, expected)
            self.assertEqual(result.ordinal, expected.ordinal)
            self.assertEqual(result.freq, expected.freq)
        for freq in ['2M', offsets.MonthEnd(2)]:
            result = p.asfreq(freq, how='S')
            expected = Period('2007-01', freq='2M')

            self.assertEqual(result, expected)
            self.assertEqual(result.ordinal, expected.ordinal)
            self.assertEqual(result.freq, expected.freq)

    def test_asfreq_mult_nat(self):
        # normal freq to mult freq
        for p in [Period('NaT', freq='A'), Period('NaT', freq='3A'),
                  Period('NaT', freq='2M'), Period('NaT', freq='3D')]:
            for freq in ['3A', offsets.YearEnd(3)]:
                result = p.asfreq(freq)
                expected = Period('NaT', freq='3A')
                self.assertEqual(result.ordinal, pd.tslib.iNaT)
                self.assertEqual(result.freq, expected.freq)

                result = p.asfreq(freq, how='S')
                expected = Period('NaT', freq='3A')
                self.assertEqual(result.ordinal, pd.tslib.iNaT)
                self.assertEqual(result.freq, expected.freq)


class TestPeriodIndex(tm.TestCase):

    def setUp(self):
        pass

    def test_hash_error(self):
        index = period_range('20010101', periods=10)
        with tm.assertRaisesRegexp(TypeError,
                                   "unhashable type: %r" %
                                   type(index).__name__):
            hash(index)

    def test_make_time_series(self):
        index = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        series = Series(1, index=index)
        tm.assertIsInstance(series, Series)

    def test_astype(self):
        idx = period_range('1990', '2009', freq='A')

        result = idx.astype('i8')
        self.assert_numpy_array_equal(result, idx.values)

    def test_constructor_use_start_freq(self):
        # GH #1118
        p = Period('4/2/2012', freq='B')
        index = PeriodIndex(start=p, periods=10)
        expected = PeriodIndex(start='4/2/2012', periods=10, freq='B')
        self.assertTrue(index.equals(expected))

    def test_constructor_field_arrays(self):
        # GH #1264

        years = np.arange(1990, 2010).repeat(4)[2:-2]
        quarters = np.tile(np.arange(1, 5), 20)[2:-2]

        index = PeriodIndex(year=years, quarter=quarters, freq='Q-DEC')
        expected = period_range('1990Q3', '2009Q2', freq='Q-DEC')
        self.assertTrue(index.equals(expected))

        index2 = PeriodIndex(year=years, quarter=quarters, freq='2Q-DEC')
        tm.assert_numpy_array_equal(index.asi8, index2.asi8)

        index = PeriodIndex(year=years, quarter=quarters)
        self.assertTrue(index.equals(expected))

        years = [2007, 2007, 2007]
        months = [1, 2]
        self.assertRaises(ValueError, PeriodIndex, year=years, month=months,
                          freq='M')
        self.assertRaises(ValueError, PeriodIndex, year=years, month=months,
                          freq='2M')
        self.assertRaises(ValueError, PeriodIndex, year=years, month=months,
                          freq='M', start=Period('2007-01', freq='M'))

        years = [2007, 2007, 2007]
        months = [1, 2, 3]
        idx = PeriodIndex(year=years, month=months, freq='M')
        exp = period_range('2007-01', periods=3, freq='M')
        self.assertTrue(idx.equals(exp))

    def test_constructor_U(self):
        # U was used as undefined period
        self.assertRaises(ValueError, period_range, '2007-1-1', periods=500,
                          freq='X')

    def test_constructor_arrays_negative_year(self):
        years = np.arange(1960, 2000).repeat(4)
        quarters = np.tile(lrange(1, 5), 40)

        pindex = PeriodIndex(year=years, quarter=quarters)

        self.assert_numpy_array_equal(pindex.year, years)
        self.assert_numpy_array_equal(pindex.quarter, quarters)

    def test_constructor_invalid_quarters(self):
        self.assertRaises(ValueError, PeriodIndex, year=lrange(2000, 2004),
                          quarter=lrange(4), freq='Q-DEC')

    def test_constructor_corner(self):
        self.assertRaises(ValueError, PeriodIndex, periods=10, freq='A')

        start = Period('2007', freq='A-JUN')
        end = Period('2010', freq='A-DEC')
        self.assertRaises(ValueError, PeriodIndex, start=start, end=end)
        self.assertRaises(ValueError, PeriodIndex, start=start)
        self.assertRaises(ValueError, PeriodIndex, end=end)

        result = period_range('2007-01', periods=10.5, freq='M')
        exp = period_range('2007-01', periods=10, freq='M')
        self.assertTrue(result.equals(exp))

    def test_constructor_fromarraylike(self):
        idx = period_range('2007-01', periods=20, freq='M')

        self.assertRaises(ValueError, PeriodIndex, idx.values)
        self.assertRaises(ValueError, PeriodIndex, list(idx.values))
        self.assertRaises(ValueError, PeriodIndex,
                          data=Period('2007', freq='A'))

        result = PeriodIndex(iter(idx))
        self.assertTrue(result.equals(idx))

        result = PeriodIndex(idx)
        self.assertTrue(result.equals(idx))

        result = PeriodIndex(idx, freq='M')
        self.assertTrue(result.equals(idx))

        result = PeriodIndex(idx, freq=offsets.MonthEnd())
        self.assertTrue(result.equals(idx))
        self.assertTrue(result.freq, 'M')

        result = PeriodIndex(idx, freq='2M')
        self.assertTrue(result.equals(idx))
        self.assertTrue(result.freq, '2M')

        result = PeriodIndex(idx, freq=offsets.MonthEnd(2))
        self.assertTrue(result.equals(idx))
        self.assertTrue(result.freq, '2M')

        result = PeriodIndex(idx, freq='D')
        exp = idx.asfreq('D', 'e')
        self.assertTrue(result.equals(exp))

    def test_constructor_datetime64arr(self):
        vals = np.arange(100000, 100000 + 10000, 100, dtype=np.int64)
        vals = vals.view(np.dtype('M8[us]'))

        self.assertRaises(ValueError, PeriodIndex, vals, freq='D')

    def test_constructor_simple_new(self):
        idx = period_range('2007-01', name='p', periods=20, freq='M')
        result = idx._simple_new(idx, 'p', freq=idx.freq)
        self.assertTrue(result.equals(idx))

        result = idx._simple_new(idx.astype('i8'), 'p', freq=idx.freq)
        self.assertTrue(result.equals(idx))

    def test_constructor_nat(self):
        self.assertRaises(
            ValueError, period_range, start='NaT', end='2011-01-01', freq='M')
        self.assertRaises(
            ValueError, period_range, start='2011-01-01', end='NaT', freq='M')

    def test_constructor_year_and_quarter(self):
        year = pd.Series([2001, 2002, 2003])
        quarter = year - 2000
        idx = PeriodIndex(year=year, quarter=quarter)
        strs = ['%dQ%d' % t for t in zip(quarter, year)]
        lops = list(map(Period, strs))
        p = PeriodIndex(lops)
        tm.assert_index_equal(p, idx)

    def test_constructor_freq_mult(self):
        # GH #7811
        for func in [PeriodIndex, period_range]:
            # must be the same, but for sure...
            pidx = func(start='2014-01', freq='2M', periods=4)
            expected = PeriodIndex(['2014-01', '2014-03', '2014-05', '2014-07'], freq='M')
            tm.assert_index_equal(pidx, expected)

            pidx = func(start='2014-01-02', end='2014-01-15', freq='3D')
            expected = PeriodIndex(['2014-01-02', '2014-01-05', '2014-01-08', '2014-01-11',
                                    '2014-01-14'], freq='D')
            tm.assert_index_equal(pidx, expected)

            pidx = func(end='2014-01-01 17:00', freq='4H', periods=3)
            expected = PeriodIndex(['2014-01-01 09:00', '2014-01-01 13:00',
                                    '2014-01-01 17:00'], freq='4H')
            tm.assert_index_equal(pidx, expected)

        msg = ('Frequency must be positive, because it'
               ' represents span: -1M')
        with tm.assertRaisesRegexp(ValueError, msg):
            PeriodIndex(['2011-01'], freq='-1M')

        msg = ('Frequency must be positive, because it'
               ' represents span: 0M')
        with tm.assertRaisesRegexp(ValueError, msg):
            PeriodIndex(['2011-01'], freq='0M')

        msg = ('Frequency must be positive, because it'
               ' represents span: 0M')
        with tm.assertRaisesRegexp(ValueError, msg):
            period_range('2011-01', periods=3, freq='0M')

    def test_constructor_freq_mult_dti_compat(self):
        import itertools
        mults = [1, 2, 3, 4, 5]
        freqs = ['A', 'M', 'D', 'T', 'S']
        for mult, freq in itertools.product(mults, freqs):
            freqstr = str(mult) + freq
            pidx = PeriodIndex(start='2014-04-01', freq=freqstr, periods=10)
            expected = date_range(start='2014-04-01', freq=freqstr, periods=10).to_period(freq)
            tm.assert_index_equal(pidx, expected)

    def test_is_(self):
        create_index = lambda: PeriodIndex(freq='A', start='1/1/2001',
                                           end='12/1/2009')
        index = create_index()
        self.assertEqual(index.is_(index), True)
        self.assertEqual(index.is_(create_index()), False)
        self.assertEqual(index.is_(index.view()), True)
        self.assertEqual(index.is_(index.view().view().view().view().view()), True)
        self.assertEqual(index.view().is_(index), True)
        ind2 = index.view()
        index.name = "Apple"
        self.assertEqual(ind2.is_(index), True)
        self.assertEqual(index.is_(index[:]), False)
        self.assertEqual(index.is_(index.asfreq('M')), False)
        self.assertEqual(index.is_(index.asfreq('A')), False)
        self.assertEqual(index.is_(index - 2), False)
        self.assertEqual(index.is_(index - 0), False)

    def test_comp_period(self):
        idx = period_range('2007-01', periods=20, freq='M')

        result = idx < idx[10]
        exp = idx.values < idx.values[10]
        self.assert_numpy_array_equal(result, exp)

    def test_getitem_ndim2(self):
        idx = period_range('2007-01', periods=3, freq='M')

        result = idx[:, None]
        # MPL kludge
        tm.assertIsInstance(result, PeriodIndex)

    def test_getitem_partial(self):
        rng = period_range('2007-01', periods=50, freq='M')
        ts = Series(np.random.randn(len(rng)), rng)

        self.assertRaises(KeyError, ts.__getitem__, '2006')

        result = ts['2008']
        self.assertTrue((result.index.year == 2008).all())

        result = ts['2008':'2009']
        self.assertEqual(len(result), 24)

        result = ts['2008-1':'2009-12']
        self.assertEqual(len(result), 24)

        result = ts['2008Q1':'2009Q4']
        self.assertEqual(len(result), 24)

        result = ts[:'2009']
        self.assertEqual(len(result), 36)

        result = ts['2009':]
        self.assertEqual(len(result), 50 - 24)

        exp = result
        result = ts[24:]
        assert_series_equal(exp, result)

        ts = ts[10:].append(ts[10:])
        self.assertRaisesRegexp(
            KeyError, "left slice bound for non-unique label: '2008'",
            ts.__getitem__, slice('2008', '2009'))

    def test_getitem_datetime(self):
        rng = period_range(start='2012-01-01', periods=10, freq='W-MON')
        ts = Series(lrange(len(rng)), index=rng)

        dt1 = datetime(2011, 10, 2)
        dt4 = datetime(2012, 4, 20)

        rs = ts[dt1:dt4]
        assert_series_equal(rs, ts)

    def test_slice_with_negative_step(self):
        ts = Series(np.arange(20),
                    period_range('2014-01', periods=20, freq='M'))
        SLC = pd.IndexSlice

        def assert_slices_equivalent(l_slc, i_slc):
            assert_series_equal(ts[l_slc], ts.iloc[i_slc])
            assert_series_equal(ts.loc[l_slc], ts.iloc[i_slc])
            assert_series_equal(ts.ix[l_slc], ts.iloc[i_slc])

        assert_slices_equivalent(SLC[Period('2014-10')::-1], SLC[9::-1])
        assert_slices_equivalent(SLC['2014-10'::-1], SLC[9::-1])

        assert_slices_equivalent(SLC[:Period('2014-10'):-1], SLC[:8:-1])
        assert_slices_equivalent(SLC[:'2014-10':-1], SLC[:8:-1])

        assert_slices_equivalent(SLC['2015-02':'2014-10':-1], SLC[13:8:-1])
        assert_slices_equivalent(SLC[Period('2015-02'):Period('2014-10'):-1], SLC[13:8:-1])
        assert_slices_equivalent(SLC['2015-02':Period('2014-10'):-1], SLC[13:8:-1])
        assert_slices_equivalent(SLC[Period('2015-02'):'2014-10':-1], SLC[13:8:-1])

        assert_slices_equivalent(SLC['2014-10':'2015-02':-1], SLC[:0])

    def test_slice_with_zero_step_raises(self):
        ts = Series(np.arange(20),
                    period_range('2014-01', periods=20, freq='M'))
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: ts[::0])
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: ts.loc[::0])
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: ts.ix[::0])

    def test_contains(self):
        rng = period_range('2007-01', freq='M', periods=10)

        self.assertTrue(Period('2007-01', freq='M') in rng)
        self.assertFalse(Period('2007-01', freq='D') in rng)
        self.assertFalse(Period('2007-01', freq='2M') in rng)

    def test_sub(self):
        rng = period_range('2007-01', periods=50)

        result = rng - 5
        exp = rng + (-5)
        self.assertTrue(result.equals(exp))

    def test_periods_number_check(self):
        self.assertRaises(
            ValueError, period_range, '2011-1-1', '2012-1-1', 'B')

    def test_tolist(self):
        index = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        rs = index.tolist()
        [tm.assertIsInstance(x, Period) for x in rs]

        recon = PeriodIndex(rs)
        self.assertTrue(index.equals(recon))

    def test_to_timestamp(self):
        index = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        series = Series(1, index=index, name='foo')

        exp_index = date_range('1/1/2001', end='12/31/2009', freq='A-DEC')
        result = series.to_timestamp(how='end')
        self.assertTrue(result.index.equals(exp_index))
        self.assertEqual(result.name, 'foo')

        exp_index = date_range('1/1/2001', end='1/1/2009', freq='AS-JAN')
        result = series.to_timestamp(how='start')
        self.assertTrue(result.index.equals(exp_index))

        def _get_with_delta(delta, freq='A-DEC'):
            return date_range(to_datetime('1/1/2001') + delta,
                              to_datetime('12/31/2009') + delta, freq=freq)

        delta = timedelta(hours=23)
        result = series.to_timestamp('H', 'end')
        exp_index = _get_with_delta(delta)
        self.assertTrue(result.index.equals(exp_index))

        delta = timedelta(hours=23, minutes=59)
        result = series.to_timestamp('T', 'end')
        exp_index = _get_with_delta(delta)
        self.assertTrue(result.index.equals(exp_index))

        result = series.to_timestamp('S', 'end')
        delta = timedelta(hours=23, minutes=59, seconds=59)
        exp_index = _get_with_delta(delta)
        self.assertTrue(result.index.equals(exp_index))

        index = PeriodIndex(freq='H', start='1/1/2001', end='1/2/2001')
        series = Series(1, index=index, name='foo')

        exp_index = date_range('1/1/2001 00:59:59', end='1/2/2001 00:59:59',
                               freq='H')
        result = series.to_timestamp(how='end')
        self.assertTrue(result.index.equals(exp_index))
        self.assertEqual(result.name, 'foo')

    def test_to_timestamp_quarterly_bug(self):
        years = np.arange(1960, 2000).repeat(4)
        quarters = np.tile(lrange(1, 5), 40)

        pindex = PeriodIndex(year=years, quarter=quarters)

        stamps = pindex.to_timestamp('D', 'end')
        expected = DatetimeIndex([x.to_timestamp('D', 'end') for x in pindex])
        self.assertTrue(stamps.equals(expected))

    def test_to_timestamp_preserve_name(self):
        index = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009',
                            name='foo')
        self.assertEqual(index.name, 'foo')

        conv = index.to_timestamp('D')
        self.assertEqual(conv.name, 'foo')

    def test_to_timestamp_repr_is_code(self):
        zs=[Timestamp('99-04-17 00:00:00',tz='UTC'),
        Timestamp('2001-04-17 00:00:00',tz='UTC'),
        Timestamp('2001-04-17 00:00:00',tz='America/Los_Angeles'),
        Timestamp('2001-04-17 00:00:00',tz=None)]
        for z in zs:
            self.assertEqual( eval(repr(z)), z)

    def test_to_timestamp_pi_nat(self):
        # GH 7228
        index = PeriodIndex(['NaT', '2011-01', '2011-02'], freq='M', name='idx')

        result = index.to_timestamp('D')
        expected = DatetimeIndex([pd.NaT, datetime(2011, 1, 1),
                                  datetime(2011, 2, 1)], name='idx')
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.name, 'idx')

        result2 = result.to_period(freq='M')
        self.assertTrue(result2.equals(index))
        self.assertEqual(result2.name, 'idx')

        result3 = result.to_period(freq='3M')
        exp = PeriodIndex(['NaT', '2011-01', '2011-02'], freq='3M', name='idx')
        self.assert_index_equal(result3, exp)
        self.assertEqual(result3.freqstr, '3M')

        msg = ('Frequency must be positive, because it'
               ' represents span: -2A')
        with tm.assertRaisesRegexp(ValueError, msg):
            result.to_period(freq='-2A')

    def test_to_timestamp_pi_mult(self):
        idx = PeriodIndex(['2011-01', 'NaT', '2011-02'], freq='2M', name='idx')
        result = idx.to_timestamp()
        expected = DatetimeIndex(['2011-01-01', 'NaT', '2011-02-01'], name='idx')
        self.assert_index_equal(result, expected)
        result = idx.to_timestamp(how='E')
        expected = DatetimeIndex(['2011-02-28', 'NaT', '2011-03-31'], name='idx')
        self.assert_index_equal(result, expected)

    def test_as_frame_columns(self):
        rng = period_range('1/1/2000', periods=5)
        df = DataFrame(randn(10, 5), columns=rng)

        ts = df[rng[0]]
        assert_series_equal(ts, df.ix[:, 0])

        # GH # 1211
        repr(df)

        ts = df['1/1/2000']
        assert_series_equal(ts, df.ix[:, 0])

    def test_indexing(self):

        # GH 4390, iat incorrectly indexing
        index = period_range('1/1/2001', periods=10)
        s = Series(randn(10), index=index)
        expected = s[index[0]]
        result = s.iat[0]
        self.assertEqual(expected, result)

    def test_frame_setitem(self):
        rng = period_range('1/1/2000', periods=5)
        rng.name = 'index'
        df = DataFrame(randn(5, 3), index=rng)

        df['Index'] = rng
        rs = Index(df['Index'])
        self.assertTrue(rs.equals(rng))

        rs = df.reset_index().set_index('index')
        tm.assertIsInstance(rs.index, PeriodIndex)
        self.assertTrue(rs.index.equals(rng))

    def test_period_set_index_reindex(self):
        # GH 6631
        df = DataFrame(np.random.random(6))
        idx1 = period_range('2011/01/01', periods=6, freq='M')
        idx2 = period_range('2013', periods=6, freq='A')

        df = df.set_index(idx1)
        self.assertTrue(df.index.equals(idx1))
        df = df.set_index(idx2)
        self.assertTrue(df.index.equals(idx2))

    def test_frame_to_time_stamp(self):
        K = 5
        index = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        df = DataFrame(randn(len(index), K), index=index)
        df['mix'] = 'a'

        exp_index = date_range('1/1/2001', end='12/31/2009', freq='A-DEC')
        result = df.to_timestamp('D', 'end')
        self.assertTrue(result.index.equals(exp_index))
        assert_almost_equal(result.values, df.values)

        exp_index = date_range('1/1/2001', end='1/1/2009', freq='AS-JAN')
        result = df.to_timestamp('D', 'start')
        self.assertTrue(result.index.equals(exp_index))

        def _get_with_delta(delta, freq='A-DEC'):
            return date_range(to_datetime('1/1/2001') + delta,
                              to_datetime('12/31/2009') + delta, freq=freq)

        delta = timedelta(hours=23)
        result = df.to_timestamp('H', 'end')
        exp_index = _get_with_delta(delta)
        self.assertTrue(result.index.equals(exp_index))

        delta = timedelta(hours=23, minutes=59)
        result = df.to_timestamp('T', 'end')
        exp_index = _get_with_delta(delta)
        self.assertTrue(result.index.equals(exp_index))

        result = df.to_timestamp('S', 'end')
        delta = timedelta(hours=23, minutes=59, seconds=59)
        exp_index = _get_with_delta(delta)
        self.assertTrue(result.index.equals(exp_index))

        # columns
        df = df.T

        exp_index = date_range('1/1/2001', end='12/31/2009', freq='A-DEC')
        result = df.to_timestamp('D', 'end', axis=1)
        self.assertTrue(result.columns.equals(exp_index))
        assert_almost_equal(result.values, df.values)

        exp_index = date_range('1/1/2001', end='1/1/2009', freq='AS-JAN')
        result = df.to_timestamp('D', 'start', axis=1)
        self.assertTrue(result.columns.equals(exp_index))

        delta = timedelta(hours=23)
        result = df.to_timestamp('H', 'end', axis=1)
        exp_index = _get_with_delta(delta)
        self.assertTrue(result.columns.equals(exp_index))

        delta = timedelta(hours=23, minutes=59)
        result = df.to_timestamp('T', 'end', axis=1)
        exp_index = _get_with_delta(delta)
        self.assertTrue(result.columns.equals(exp_index))

        result = df.to_timestamp('S', 'end', axis=1)
        delta = timedelta(hours=23, minutes=59, seconds=59)
        exp_index = _get_with_delta(delta)
        self.assertTrue(result.columns.equals(exp_index))

        # invalid axis
        assertRaisesRegexp(ValueError, 'axis', df.to_timestamp, axis=2)

        result1 = df.to_timestamp('5t', axis=1)
        result2 = df.to_timestamp('t', axis=1)
        expected = pd.date_range('2001-01-01', '2009-01-01', freq='AS')
        self.assertTrue(isinstance(result1.columns, DatetimeIndex))
        self.assertTrue(isinstance(result2.columns, DatetimeIndex))
        self.assert_numpy_array_equal(result1.columns.asi8, expected.asi8)
        self.assert_numpy_array_equal(result2.columns.asi8, expected.asi8)
        # PeriodIndex.to_timestamp always use 'infer'
        self.assertEqual(result1.columns.freqstr, 'AS-JAN')
        self.assertEqual(result2.columns.freqstr, 'AS-JAN')

    def test_index_duplicate_periods(self):
        # monotonic
        idx = PeriodIndex([2000, 2007, 2007, 2009, 2009], freq='A-JUN')
        ts = Series(np.random.randn(len(idx)), index=idx)

        result = ts[2007]
        expected = ts[1:3]
        assert_series_equal(result, expected)
        result[:] = 1
        self.assertTrue((ts[1:3] == 1).all())

        # not monotonic
        idx = PeriodIndex([2000, 2007, 2007, 2009, 2007], freq='A-JUN')
        ts = Series(np.random.randn(len(idx)), index=idx)

        result = ts[2007]
        expected = ts[idx == 2007]
        assert_series_equal(result, expected)

    def test_index_unique(self):
        idx = PeriodIndex([2000, 2007, 2007, 2009, 2009], freq='A-JUN')
        expected = PeriodIndex([2000, 2007, 2009], freq='A-JUN')
        self.assert_numpy_array_equal(idx.unique(), expected.values)
        self.assertEqual(idx.nunique(), 3)

        idx = PeriodIndex([2000, 2007, 2007, 2009, 2007], freq='A-JUN', tz='US/Eastern')
        expected = PeriodIndex([2000, 2007, 2009], freq='A-JUN', tz='US/Eastern')
        self.assert_numpy_array_equal(idx.unique(), expected.values)
        self.assertEqual(idx.nunique(), 3)

    def test_constructor(self):
        pi = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        assert_equal(len(pi), 9)

        pi = PeriodIndex(freq='Q', start='1/1/2001', end='12/1/2009')
        assert_equal(len(pi), 4 * 9)

        pi = PeriodIndex(freq='M', start='1/1/2001', end='12/1/2009')
        assert_equal(len(pi), 12 * 9)

        pi = PeriodIndex(freq='D', start='1/1/2001', end='12/31/2009')
        assert_equal(len(pi), 365 * 9 + 2)

        pi = PeriodIndex(freq='B', start='1/1/2001', end='12/31/2009')
        assert_equal(len(pi), 261 * 9)

        pi = PeriodIndex(freq='H', start='1/1/2001', end='12/31/2001 23:00')
        assert_equal(len(pi), 365 * 24)

        pi = PeriodIndex(freq='Min', start='1/1/2001', end='1/1/2001 23:59')
        assert_equal(len(pi), 24 * 60)

        pi = PeriodIndex(freq='S', start='1/1/2001', end='1/1/2001 23:59:59')
        assert_equal(len(pi), 24 * 60 * 60)

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
        self.assertTrue((i1 == i2).all())
        assert_equal(i1.freq, i2.freq)

        end_intv = Period('2006-12-31', ('w', 1))
        i2 = PeriodIndex(end=end_intv, periods=10)
        assert_equal(len(i1), len(i2))
        self.assertTrue((i1 == i2).all())
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
            raise AssertionError(
                'Must specify periods if missing start or end')
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
        pi1 = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='A', start='1/1/2002', end='12/1/2010')

        self.assertTrue(pi1.shift(0).equals(pi1))

        assert_equal(len(pi1), len(pi2))
        assert_equal(pi1.shift(1).values, pi2.values)

        pi1 = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='A', start='1/1/2000', end='12/1/2008')
        assert_equal(len(pi1), len(pi2))
        assert_equal(pi1.shift(-1).values, pi2.values)

        pi1 = PeriodIndex(freq='M', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='M', start='2/1/2001', end='1/1/2010')
        assert_equal(len(pi1), len(pi2))
        assert_equal(pi1.shift(1).values, pi2.values)

        pi1 = PeriodIndex(freq='M', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='M', start='12/1/2000', end='11/1/2009')
        assert_equal(len(pi1), len(pi2))
        assert_equal(pi1.shift(-1).values, pi2.values)

        pi1 = PeriodIndex(freq='D', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='D', start='1/2/2001', end='12/2/2009')
        assert_equal(len(pi1), len(pi2))
        assert_equal(pi1.shift(1).values, pi2.values)

        pi1 = PeriodIndex(freq='D', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='D', start='12/31/2000', end='11/30/2009')
        assert_equal(len(pi1), len(pi2))
        assert_equal(pi1.shift(-1).values, pi2.values)

    def test_shift_nat(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'], freq='M', name='idx')
        result = idx.shift(1)
        expected = PeriodIndex(['2011-02', '2011-03', 'NaT', '2011-05'], freq='M', name='idx')
        self.assertTrue(result.equals(expected))
        self.assertEqual(result.name, expected.name)

    def test_shift_ndarray(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'], freq='M', name='idx')
        result = idx.shift(np.array([1, 2, 3, 4]))
        expected = PeriodIndex(['2011-02', '2011-04', 'NaT', '2011-08'], freq='M', name='idx')
        self.assertTrue(result.equals(expected))

        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'], freq='M', name='idx')
        result = idx.shift(np.array([1, -2, 3, -4]))
        expected = PeriodIndex(['2011-02', '2010-12', 'NaT', '2010-12'], freq='M', name='idx')
        self.assertTrue(result.equals(expected))

    def test_asfreq(self):
        pi1 = PeriodIndex(freq='A', start='1/1/2001', end='1/1/2001')
        pi2 = PeriodIndex(freq='Q', start='1/1/2001', end='1/1/2001')
        pi3 = PeriodIndex(freq='M', start='1/1/2001', end='1/1/2001')
        pi4 = PeriodIndex(freq='D', start='1/1/2001', end='1/1/2001')
        pi5 = PeriodIndex(freq='H', start='1/1/2001', end='1/1/2001 00:00')
        pi6 = PeriodIndex(freq='Min', start='1/1/2001', end='1/1/2001 00:00')
        pi7 = PeriodIndex(freq='S', start='1/1/2001', end='1/1/2001 00:00:00')

        self.assertEqual(pi1.asfreq('Q', 'S'), pi2)
        self.assertEqual(pi1.asfreq('Q', 's'), pi2)
        self.assertEqual(pi1.asfreq('M', 'start'), pi3)
        self.assertEqual(pi1.asfreq('D', 'StarT'), pi4)
        self.assertEqual(pi1.asfreq('H', 'beGIN'), pi5)
        self.assertEqual(pi1.asfreq('Min', 'S'), pi6)
        self.assertEqual(pi1.asfreq('S', 'S'), pi7)

        self.assertEqual(pi2.asfreq('A', 'S'), pi1)
        self.assertEqual(pi2.asfreq('M', 'S'), pi3)
        self.assertEqual(pi2.asfreq('D', 'S'), pi4)
        self.assertEqual(pi2.asfreq('H', 'S'), pi5)
        self.assertEqual(pi2.asfreq('Min', 'S'), pi6)
        self.assertEqual(pi2.asfreq('S', 'S'), pi7)

        self.assertEqual(pi3.asfreq('A', 'S'), pi1)
        self.assertEqual(pi3.asfreq('Q', 'S'), pi2)
        self.assertEqual(pi3.asfreq('D', 'S'), pi4)
        self.assertEqual(pi3.asfreq('H', 'S'), pi5)
        self.assertEqual(pi3.asfreq('Min', 'S'), pi6)
        self.assertEqual(pi3.asfreq('S', 'S'), pi7)

        self.assertEqual(pi4.asfreq('A', 'S'), pi1)
        self.assertEqual(pi4.asfreq('Q', 'S'), pi2)
        self.assertEqual(pi4.asfreq('M', 'S'), pi3)
        self.assertEqual(pi4.asfreq('H', 'S'), pi5)
        self.assertEqual(pi4.asfreq('Min', 'S'), pi6)
        self.assertEqual(pi4.asfreq('S', 'S'), pi7)

        self.assertEqual(pi5.asfreq('A', 'S'), pi1)
        self.assertEqual(pi5.asfreq('Q', 'S'), pi2)
        self.assertEqual(pi5.asfreq('M', 'S'), pi3)
        self.assertEqual(pi5.asfreq('D', 'S'), pi4)
        self.assertEqual(pi5.asfreq('Min', 'S'), pi6)
        self.assertEqual(pi5.asfreq('S', 'S'), pi7)

        self.assertEqual(pi6.asfreq('A', 'S'), pi1)
        self.assertEqual(pi6.asfreq('Q', 'S'), pi2)
        self.assertEqual(pi6.asfreq('M', 'S'), pi3)
        self.assertEqual(pi6.asfreq('D', 'S'), pi4)
        self.assertEqual(pi6.asfreq('H', 'S'), pi5)
        self.assertEqual(pi6.asfreq('S', 'S'), pi7)

        self.assertEqual(pi7.asfreq('A', 'S'), pi1)
        self.assertEqual(pi7.asfreq('Q', 'S'), pi2)
        self.assertEqual(pi7.asfreq('M', 'S'), pi3)
        self.assertEqual(pi7.asfreq('D', 'S'), pi4)
        self.assertEqual(pi7.asfreq('H', 'S'), pi5)
        self.assertEqual(pi7.asfreq('Min', 'S'), pi6)

        self.assertRaises(ValueError, pi7.asfreq, 'T', 'foo')
        result1 = pi1.asfreq('3M')
        result2 = pi1.asfreq('M')
        expected = PeriodIndex(freq='M', start='2001-12', end='2001-12')
        self.assert_numpy_array_equal(result1.asi8, expected.asi8)
        self.assertEqual(result1.freqstr, '3M')
        self.assert_numpy_array_equal(result2.asi8, expected.asi8)
        self.assertEqual(result2.freqstr, 'M')

    def test_asfreq_nat(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'], freq='M')
        result = idx.asfreq(freq='Q')
        expected = PeriodIndex(['2011Q1', '2011Q1', 'NaT', '2011Q2'], freq='Q')
        self.assertTrue(result.equals(expected))

    def test_asfreq_mult_pi(self):
        pi = PeriodIndex(['2001-01', '2001-02', 'NaT', '2001-03'], freq='2M')

        for freq in ['D', '3D']:
            result = pi.asfreq(freq)
            exp = PeriodIndex(['2001-02-28', '2001-03-31', 'NaT',
                               '2001-04-30'], freq=freq)
            self.assert_index_equal(result, exp)
            self.assertEqual(result.freq, exp.freq)

            result = pi.asfreq(freq, how='S')
            exp = PeriodIndex(['2001-01-01', '2001-02-01', 'NaT',
                               '2001-03-01'], freq=freq)
            self.assert_index_equal(result, exp)
            self.assertEqual(result.freq, exp.freq)

    def test_period_index_length(self):
        pi = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        assert_equal(len(pi), 9)

        pi = PeriodIndex(freq='Q', start='1/1/2001', end='12/1/2009')
        assert_equal(len(pi), 4 * 9)

        pi = PeriodIndex(freq='M', start='1/1/2001', end='12/1/2009')
        assert_equal(len(pi), 12 * 9)

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
        self.assertTrue((i1 == i2).all())
        assert_equal(i1.freq, i2.freq)

        end_intv = Period('2006-12-31', ('w', 1))
        i2 = PeriodIndex(end=end_intv, periods=10)
        assert_equal(len(i1), len(i2))
        self.assertTrue((i1 == i2).all())
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
            raise AssertionError(
                'Must specify periods if missing start or end')
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

    def test_frame_index_to_string(self):
        index = PeriodIndex(['2011-1', '2011-2', '2011-3'], freq='M')
        frame = DataFrame(np.random.randn(3, 4), index=index)

        # it works!
        frame.to_string()

    def test_asfreq_ts(self):
        index = PeriodIndex(freq='A', start='1/1/2001', end='12/31/2010')
        ts = Series(np.random.randn(len(index)), index=index)
        df = DataFrame(np.random.randn(len(index), 3), index=index)

        result = ts.asfreq('D', how='end')
        df_result = df.asfreq('D', how='end')
        exp_index = index.asfreq('D', how='end')
        self.assertEqual(len(result), len(ts))
        self.assertTrue(result.index.equals(exp_index))
        self.assertTrue(df_result.index.equals(exp_index))

        result = ts.asfreq('D', how='start')
        self.assertEqual(len(result), len(ts))
        self.assertTrue(result.index.equals(index.asfreq('D', how='start')))

    def test_badinput(self):
        self.assertRaises(datetools.DateParseError, Period, '1/1/-2000', 'A')
        # self.assertRaises(datetools.DateParseError, Period, '-2000', 'A')
        # self.assertRaises(datetools.DateParseError, Period, '0', 'A')

    def test_negative_ordinals(self):
        p = Period(ordinal=-1000, freq='A')
        p = Period(ordinal=0, freq='A')

        idx1 = PeriodIndex(ordinal=[-1, 0, 1], freq='A')
        idx2 = PeriodIndex(ordinal=np.array([-1, 0, 1]), freq='A')
        tm.assert_numpy_array_equal(idx1,idx2)

    def test_dti_to_period(self):
        dti = DatetimeIndex(start='1/1/2005', end='12/1/2005', freq='M')
        pi1 = dti.to_period()
        pi2 = dti.to_period(freq='D')
        pi3 = dti.to_period(freq='3D')

        self.assertEqual(pi1[0], Period('Jan 2005', freq='M'))
        self.assertEqual(pi2[0], Period('1/31/2005', freq='D'))
        self.assertEqual(pi3[0], Period('1/31/2005', freq='3D'))

        self.assertEqual(pi1[-1], Period('Nov 2005', freq='M'))
        self.assertEqual(pi2[-1], Period('11/30/2005', freq='D'))
        self.assertEqual(pi3[-1], Period('11/30/2005', freq='3D'))

        tm.assert_index_equal(pi1, period_range('1/1/2005', '11/1/2005', freq='M'))
        tm.assert_index_equal(pi2, period_range('1/1/2005', '11/1/2005', freq='M').asfreq('D'))
        tm.assert_index_equal(pi3, period_range('1/1/2005', '11/1/2005', freq='M').asfreq('3D'))

    def test_pindex_slice_index(self):
        pi = PeriodIndex(start='1/1/10', end='12/31/12', freq='M')
        s = Series(np.random.rand(len(pi)), index=pi)
        res = s['2010']
        exp = s[0:12]
        assert_series_equal(res, exp)
        res = s['2011']
        exp = s[12:24]
        assert_series_equal(res, exp)

    def test_getitem_day(self):
        # GH 6716
        # Confirm DatetimeIndex and PeriodIndex works identically
        didx = DatetimeIndex(start='2013/01/01', freq='D', periods=400)
        pidx = PeriodIndex(start='2013/01/01', freq='D', periods=400)

        for idx in [didx, pidx]:
            # getitem against index should raise ValueError
            values = ['2014', '2013/02', '2013/01/02',
                      '2013/02/01 9H', '2013/02/01 09:00']
            for v in values:

                if _np_version_under1p9:
                    with tm.assertRaises(ValueError):
                        idx[v]
                else:
                    # GH7116
                    # these show deprecations as we are trying
                    # to slice with non-integer indexers
                    #with tm.assertRaises(IndexError):
                    #    idx[v]
                    continue

            s = Series(np.random.rand(len(idx)), index=idx)
            assert_series_equal(s['2013/01'], s[0:31])
            assert_series_equal(s['2013/02'], s[31:59])
            assert_series_equal(s['2014'], s[365:])

            invalid = ['2013/02/01 9H', '2013/02/01 09:00']
            for v in invalid:
                with tm.assertRaises(KeyError):
                    s[v]

    def test_range_slice_day(self):
        # GH 6716
        didx = DatetimeIndex(start='2013/01/01', freq='D', periods=400)
        pidx = PeriodIndex(start='2013/01/01', freq='D', periods=400)

        for idx in [didx, pidx]:
            # slices against index should raise IndexError
            values = ['2014', '2013/02', '2013/01/02',
                      '2013/02/01 9H', '2013/02/01 09:00']
            for v in values:
                with tm.assertRaises(IndexError):
                    idx[v:]

            s = Series(np.random.rand(len(idx)), index=idx)

            assert_series_equal(s['2013/01/02':], s[1:])
            assert_series_equal(s['2013/01/02':'2013/01/05'], s[1:5])
            assert_series_equal(s['2013/02':], s[31:])
            assert_series_equal(s['2014':], s[365:])

            invalid = ['2013/02/01 9H', '2013/02/01 09:00']
            for v in invalid:
                with tm.assertRaises(IndexError):
                    idx[v:]

    def test_getitem_seconds(self):
        # GH 6716
        didx = DatetimeIndex(start='2013/01/01 09:00:00', freq='S', periods=4000)
        pidx = PeriodIndex(start='2013/01/01 09:00:00', freq='S', periods=4000)

        for idx in [didx, pidx]:
            # getitem against index should raise ValueError
            values = ['2014', '2013/02', '2013/01/02',
                      '2013/02/01 9H', '2013/02/01 09:00']
            for v in values:
                if _np_version_under1p9:
                    with tm.assertRaises(ValueError):
                        idx[v]
                else:
                    # GH7116
                    # these show deprecations as we are trying
                    # to slice with non-integer indexers
                    #with tm.assertRaises(IndexError):
                    #    idx[v]
                    continue

            s = Series(np.random.rand(len(idx)), index=idx)
            assert_series_equal(s['2013/01/01 10:00'], s[3600:3660])
            assert_series_equal(s['2013/01/01 9H'], s[:3600])
            for d in ['2013/01/01', '2013/01', '2013']:
                assert_series_equal(s[d], s)

    def test_range_slice_seconds(self):
        # GH 6716
        didx = DatetimeIndex(start='2013/01/01 09:00:00', freq='S', periods=4000)
        pidx = PeriodIndex(start='2013/01/01 09:00:00', freq='S', periods=4000)

        for idx in [didx, pidx]:
            # slices against index should raise IndexError
            values = ['2014', '2013/02', '2013/01/02',
                      '2013/02/01 9H', '2013/02/01 09:00']
            for v in values:
                with tm.assertRaises(IndexError):
                    idx[v:]

            s = Series(np.random.rand(len(idx)), index=idx)

            assert_series_equal(s['2013/01/01 09:05':'2013/01/01 09:10'], s[300:660])
            assert_series_equal(s['2013/01/01 10:00':'2013/01/01 10:05'], s[3600:3960])
            assert_series_equal(s['2013/01/01 10H':], s[3600:])
            assert_series_equal(s[:'2013/01/01 09:30'], s[:1860])
            for d in ['2013/01/01', '2013/01', '2013']:
                assert_series_equal(s[d:], s)

    def test_range_slice_outofbounds(self):
        # GH 5407
        didx = DatetimeIndex(start='2013/10/01', freq='D', periods=10)
        pidx = PeriodIndex(start='2013/10/01', freq='D', periods=10)

        for idx in [didx, pidx]:
            df = DataFrame(dict(units=[100 + i for i in range(10)]), index=idx)
            empty = DataFrame(index=idx.__class__([], freq='D'), columns=['units'])
            empty['units'] = empty['units'].astype('int64')

            tm.assert_frame_equal(df['2013/09/01':'2013/09/30'], empty)
            tm.assert_frame_equal(df['2013/09/30':'2013/10/02'], df.iloc[:2])
            tm.assert_frame_equal(df['2013/10/01':'2013/10/02'], df.iloc[:2])
            tm.assert_frame_equal(df['2013/10/02':'2013/09/30'], empty)
            tm.assert_frame_equal(df['2013/10/15':'2013/10/17'], empty)
            tm.assert_frame_equal(df['2013-06':'2013-09'], empty)
            tm.assert_frame_equal(df['2013-11':'2013-12'], empty)

    def test_pindex_fieldaccessor_nat(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2012-03', '2012-04'], freq='D')
        self.assert_numpy_array_equal(idx.year, np.array([2011, 2011, -1, 2012, 2012]))
        self.assert_numpy_array_equal(idx.month, np.array([1, 2, -1, 3, 4]))

    def test_pindex_qaccess(self):
        pi = PeriodIndex(['2Q05', '3Q05', '4Q05', '1Q06', '2Q06'], freq='Q')
        s = Series(np.random.rand(len(pi)), index=pi).cumsum()
        # Todo: fix these accessors!
        self.assertEqual(s['05Q4'], s[2])

    def test_period_dt64_round_trip(self):
        dti = date_range('1/1/2000', '1/7/2002', freq='B')
        pi = dti.to_period()
        self.assertTrue(pi.to_timestamp().equals(dti))

        dti = date_range('1/1/2000', '1/7/2002', freq='B')
        pi = dti.to_period(freq='H')
        self.assertTrue(pi.to_timestamp().equals(dti))

    def test_to_period_quarterly(self):
        # make sure we can make the round trip
        for month in MONTHS:
            freq = 'Q-%s' % month
            rng = period_range('1989Q3', '1991Q3', freq=freq)
            stamps = rng.to_timestamp()
            result = stamps.to_period(freq)
            self.assertTrue(rng.equals(result))

    def test_to_period_quarterlyish(self):
        offsets = ['BQ', 'QS', 'BQS']
        for off in offsets:
            rng = date_range('01-Jan-2012', periods=8, freq=off)
            prng = rng.to_period()
            self.assertEqual(prng.freq, 'Q-DEC')

    def test_to_period_annualish(self):
        offsets = ['BA', 'AS', 'BAS']
        for off in offsets:
            rng = date_range('01-Jan-2012', periods=8, freq=off)
            prng = rng.to_period()
            self.assertEqual(prng.freq, 'A-DEC')

    def test_to_period_monthish(self):
        offsets = ['MS', 'BM']
        for off in offsets:
            rng = date_range('01-Jan-2012', periods=8, freq=off)
            prng = rng.to_period()
            self.assertEqual(prng.freq, 'M')

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            rng = date_range('01-Jan-2012', periods=8, freq='EOM')
        prng = rng.to_period()
        self.assertEqual(prng.freq, 'M')

    def test_multiples(self):
        result1 = Period('1989', freq='2A')
        result2 = Period('1989', freq='A')
        self.assertEqual(result1.ordinal, result2.ordinal)
        self.assertEqual(result1.freqstr, '2A-DEC')
        self.assertEqual(result2.freqstr, 'A-DEC')
        self.assertEqual(result1.freq, offsets.YearEnd(2))
        self.assertEqual(result2.freq, offsets.YearEnd())

        self.assertEqual((result1 + 1).ordinal, result1.ordinal + 2)
        self.assertEqual((result1 - 1).ordinal, result2.ordinal - 2)

    def test_pindex_multiples(self):
        pi = PeriodIndex(start='1/1/11', end='12/31/11', freq='2M')
        expected = PeriodIndex(['2011-01', '2011-03', '2011-05', '2011-07',
                                '2011-09', '2011-11'], freq='M')
        tm.assert_index_equal(pi, expected)
        self.assertEqual(pi.freq, offsets.MonthEnd(2))
        self.assertEqual(pi.freqstr, '2M')

        pi = period_range(start='1/1/11', end='12/31/11', freq='2M')
        tm.assert_index_equal(pi, expected)
        self.assertEqual(pi.freq, offsets.MonthEnd(2))
        self.assertEqual(pi.freqstr, '2M')

        pi = period_range(start='1/1/11', periods=6, freq='2M')
        tm.assert_index_equal(pi, expected)
        self.assertEqual(pi.freq, offsets.MonthEnd(2))
        self.assertEqual(pi.freqstr, '2M')

    def test_iteration(self):
        index = PeriodIndex(start='1/1/10', periods=4, freq='B')

        result = list(index)
        tm.assertIsInstance(result[0], Period)
        self.assertEqual(result[0].freq, index.freq)

    def test_take(self):
        index = PeriodIndex(start='1/1/10', end='12/31/12', freq='D', name='idx')
        expected = PeriodIndex([datetime(2010, 1, 6), datetime(2010, 1, 7),
                                datetime(2010, 1, 9), datetime(2010, 1, 13)],
                                freq='D', name='idx')

        taken1 = index.take([5, 6, 8, 12])
        taken2 = index[[5, 6, 8, 12]]

        for taken in [taken1, taken2]:
            self.assertTrue(taken.equals(expected))
            tm.assertIsInstance(taken, PeriodIndex)
            self.assertEqual(taken.freq, index.freq)
            self.assertEqual(taken.name, expected.name)

    def test_joins(self):
        index = period_range('1/1/2000', '1/20/2000', freq='D')

        for kind in ['inner', 'outer', 'left', 'right']:
            joined = index.join(index[:-5], how=kind)

            tm.assertIsInstance(joined, PeriodIndex)
            self.assertEqual(joined.freq, index.freq)

    def test_join_self(self):
        index = period_range('1/1/2000', '1/20/2000', freq='D')

        for kind in ['inner', 'outer', 'left', 'right']:
            res = index.join(index, how=kind)
            self.assertIs(index, res)

    def test_join_does_not_recur(self):
        df = tm.makeCustomDataframe(3, 2, data_gen_f=lambda *args:
                                    np.random.randint(2), c_idx_type='p',
                                    r_idx_type='dt')
        s = df.iloc[:2, 0]

        res = s.index.join(df.columns, how='outer')
        expected = Index([s.index[0], s.index[1],
                          df.columns[0], df.columns[1]], object)
        tm.assert_index_equal(res, expected)

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
        msg = "Input has different freq=D from PeriodIndex\\(freq=A-DEC\\)"
        with assertRaisesRegexp(ValueError, msg):
            ts + ts.asfreq('D', how="end")

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
        self.assertTrue(result.equals(index))

        # not in order
        result = _permute(index[:-5]).union(_permute(index[10:]))
        self.assertTrue(result.equals(index))

        # raise if different frequencies
        index = period_range('1/1/2000', '1/20/2000', freq='D')
        index2 = period_range('1/1/2000', '1/20/2000', freq='W-WED')
        self.assertRaises(ValueError, index.union, index2)

        self.assertRaises(ValueError, index.join, index.to_timestamp())

        index3 = period_range('1/1/2000', '1/20/2000', freq='2D')
        self.assertRaises(ValueError, index.join, index3)

    def test_intersection(self):
        index = period_range('1/1/2000', '1/20/2000', freq='D')

        result = index[:-5].intersection(index[10:])
        self.assertTrue(result.equals(index[10:-5]))

        # not in order
        left = _permute(index[:-5])
        right = _permute(index[10:])
        result = left.intersection(right).sort_values()
        self.assertTrue(result.equals(index[10:-5]))

        # raise if different frequencies
        index = period_range('1/1/2000', '1/20/2000', freq='D')
        index2 = period_range('1/1/2000', '1/20/2000', freq='W-WED')
        self.assertRaises(ValueError, index.intersection, index2)

        index3 = period_range('1/1/2000', '1/20/2000', freq='2D')
        self.assertRaises(ValueError, index.intersection, index3)

    def test_fields(self):
        # year, month, day, hour, minute
        # second, weekofyear, week, dayofweek, weekday, dayofyear, quarter
        # qyear
        pi = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2005')
        self._check_all_fields(pi)

        pi = PeriodIndex(freq='Q', start='1/1/2001', end='12/1/2002')
        self._check_all_fields(pi)

        pi = PeriodIndex(freq='M', start='1/1/2001', end='1/1/2002')
        self._check_all_fields(pi)

        pi = PeriodIndex(freq='D', start='12/1/2001', end='6/1/2001')
        self._check_all_fields(pi)

        pi = PeriodIndex(freq='B', start='12/1/2001', end='6/1/2001')
        self._check_all_fields(pi)

        pi = PeriodIndex(freq='H', start='12/31/2001', end='1/1/2002 23:00')
        self._check_all_fields(pi)

        pi = PeriodIndex(freq='Min', start='12/31/2001', end='1/1/2002 00:20')
        self._check_all_fields(pi)

        pi = PeriodIndex(freq='S', start='12/31/2001 00:00:00',
                         end='12/31/2001 00:05:00')
        self._check_all_fields(pi)

        end_intv = Period('2006-12-31', 'W')
        i1 = PeriodIndex(end=end_intv, periods=10)
        self._check_all_fields(i1)

    def _check_all_fields(self, periodindex):
        fields = ['year', 'month', 'day', 'hour', 'minute',
                  'second', 'weekofyear', 'week', 'dayofweek',
                  'weekday', 'dayofyear', 'quarter', 'qyear', 'days_in_month']

        periods = list(periodindex)

        for field in fields:
            field_idx = getattr(periodindex, field)
            assert_equal(len(periodindex), len(field_idx))
            for x, val in zip(periods, field_idx):
                assert_equal(getattr(x, field), val)

    def test_is_full(self):
        index = PeriodIndex([2005, 2007, 2009], freq='A')
        self.assertFalse(index.is_full)

        index = PeriodIndex([2005, 2006, 2007], freq='A')
        self.assertTrue(index.is_full)

        index = PeriodIndex([2005, 2005, 2007], freq='A')
        self.assertFalse(index.is_full)

        index = PeriodIndex([2005, 2005, 2006], freq='A')
        self.assertTrue(index.is_full)

        index = PeriodIndex([2006, 2005, 2005], freq='A')
        self.assertRaises(ValueError, getattr, index, 'is_full')

        self.assertTrue(index[:0].is_full)

    def test_map(self):
        index = PeriodIndex([2005, 2007, 2009], freq='A')
        result = index.map(lambda x: x + 1)
        expected = index + 1
        self.assertTrue(result.equals(expected))

        result = index.map(lambda x: x.ordinal)
        exp = [x.ordinal for x in index]
        tm.assert_numpy_array_equal(result, exp)

    def test_map_with_string_constructor(self):
        raw = [2005, 2007, 2009]
        index = PeriodIndex(raw, freq='A')
        types = str,

        if compat.PY3:
            # unicode
            types += compat.text_type,

        for t in types:
            expected = np.array(lmap(t, raw), dtype=object)
            res = index.map(t)

            # should return an array
            tm.assertIsInstance(res, np.ndarray)

            # preserve element types
            self.assertTrue(all(isinstance(resi, t) for resi in res))

            # dtype should be object
            self.assertEqual(res.dtype, np.dtype('object').type)

            # lastly, values should compare equal
            tm.assert_numpy_array_equal(res, expected)

    def test_convert_array_of_periods(self):
        rng = period_range('1/1/2000', periods=20, freq='D')
        periods = list(rng)

        result = pd.Index(periods)
        tm.assertIsInstance(result, PeriodIndex)

    def test_with_multi_index(self):
        # #1705
        index = date_range('1/1/2012', periods=4, freq='12H')
        index_as_arrays = [index.to_period(freq='D'), index.hour]

        s = Series([0, 1, 2, 3], index_as_arrays)

        tm.assertIsInstance(s.index.levels[0], PeriodIndex)

        tm.assertIsInstance(s.index.values[0][0], Period)

    def test_to_datetime_1703(self):
        index = period_range('1/1/2012', periods=4, freq='D')

        result = index.to_datetime()
        self.assertEqual(result[0], Timestamp('1/1/2012'))

    def test_get_loc_msg(self):
        idx = period_range('2000-1-1', freq='A', periods=10)
        bad_period = Period('2012', 'A')
        self.assertRaises(KeyError, idx.get_loc, bad_period)

        try:
            idx.get_loc(bad_period)
        except KeyError as inst:
            self.assertEqual(inst.args[0], bad_period)

    def test_append_concat(self):
        # #1815
        d1 = date_range('12/31/1990', '12/31/1999', freq='A-DEC')
        d2 = date_range('12/31/2000', '12/31/2009', freq='A-DEC')

        s1 = Series(np.random.randn(10), d1)
        s2 = Series(np.random.randn(10), d2)

        s1 = s1.to_period()
        s2 = s2.to_period()

        # drops index
        result = pd.concat([s1, s2])
        tm.assertIsInstance(result.index, PeriodIndex)
        self.assertEqual(result.index[0], s1.index[0])

    def test_pickle_freq(self):
        # GH2891
        prng = period_range('1/1/2011', '1/1/2012', freq='M')
        new_prng = self.round_trip_pickle(prng)
        self.assertEqual(new_prng.freq, offsets.MonthEnd())
        self.assertEqual(new_prng.freqstr, 'M')

    def test_slice_keep_name(self):
        idx = period_range('20010101', periods=10, freq='D', name='bob')
        self.assertEqual(idx.name, idx[1:].name)

    def test_factorize(self):
        idx1 = PeriodIndex(['2014-01', '2014-01', '2014-02', '2014-02',
                            '2014-03', '2014-03'], freq='M')

        exp_arr = np.array([0, 0, 1, 1, 2, 2])
        exp_idx = PeriodIndex(['2014-01', '2014-02', '2014-03'], freq='M')

        arr, idx = idx1.factorize()
        self.assert_numpy_array_equal(arr, exp_arr)
        self.assertTrue(idx.equals(exp_idx))

        arr, idx = idx1.factorize(sort=True)
        self.assert_numpy_array_equal(arr, exp_arr)
        self.assertTrue(idx.equals(exp_idx))

        idx2 = pd.PeriodIndex(['2014-03', '2014-03', '2014-02', '2014-01',
                               '2014-03', '2014-01'], freq='M')

        exp_arr = np.array([2, 2, 1, 0, 2, 0])
        arr, idx = idx2.factorize(sort=True)
        self.assert_numpy_array_equal(arr, exp_arr)
        self.assertTrue(idx.equals(exp_idx))

        exp_arr = np.array([0, 0, 1, 2, 0, 2])
        exp_idx = PeriodIndex(['2014-03', '2014-02', '2014-01'], freq='M')
        arr, idx = idx2.factorize()
        self.assert_numpy_array_equal(arr, exp_arr)
        self.assertTrue(idx.equals(exp_idx))

    def test_recreate_from_data(self):
        for o in ['M', 'Q', 'A', 'D', 'B', 'T', 'S', 'L', 'U', 'N', 'H']:
            org = PeriodIndex(start='2001/04/01', freq=o, periods=1)
            idx = PeriodIndex(org.values, freq=o)
            self.assertTrue(idx.equals(org))

    def test_combine_first(self):
        # GH 3367
        didx = pd.DatetimeIndex(start='1950-01-31', end='1950-07-31', freq='M')
        pidx = pd.PeriodIndex(start=pd.Period('1950-1'), end=pd.Period('1950-7'), freq='M')
        # check to be consistent with DatetimeIndex
        for idx in [didx, pidx]:
            a = pd.Series([1, np.nan, np.nan, 4, 5, np.nan, 7], index=idx)
            b = pd.Series([9, 9, 9, 9, 9, 9, 9], index=idx)

            result = a.combine_first(b)
            expected = pd.Series([1, 9, 9, 4, 5, 9, 7], index=idx, dtype=np.float64)
            tm.assert_series_equal(result, expected)

    def test_searchsorted(self):
        for freq in ['D', '2D']:
            pidx = pd.PeriodIndex(['2014-01-01', '2014-01-02', '2014-01-03',
                                   '2014-01-04', '2014-01-05'], freq=freq)

            p1 = pd.Period('2014-01-01', freq=freq)
            self.assertEqual(pidx.searchsorted(p1), 0)

            p2 = pd.Period('2014-01-04', freq=freq)
            self.assertEqual(pidx.searchsorted(p2), 3)

            msg = "Input has different freq=H from PeriodIndex"
            with self.assertRaisesRegexp(ValueError, msg):
                pidx.searchsorted(pd.Period('2014-01-01', freq='H'))

            msg = "Input has different freq=5D from PeriodIndex"
            with self.assertRaisesRegexp(ValueError, msg):
                pidx.searchsorted(pd.Period('2014-01-01', freq='5D'))


    def test_round_trip(self):

        p = Period('2000Q1')
        new_p = self.round_trip_pickle(p)
        self.assertEqual(new_p, p)

def _permute(obj):
    return obj.take(np.random.permutation(len(obj)))


class TestMethods(tm.TestCase):
    "Base test class for MaskedArrays."

    def test_add(self):
        dt1 = Period(freq='D', year=2008, month=1, day=1)
        dt2 = Period(freq='D', year=2008, month=1, day=2)
        assert_equal(dt1 + 1, dt2)
        #
        # GH 4731
        msg = "unsupported operand type\(s\)"
        with tm.assertRaisesRegexp(TypeError, msg):
            dt1 + "str"

        with tm.assertRaisesRegexp(TypeError, msg):
            dt1 + dt2

    def test_add_offset(self):
        # freq is DateOffset
        for freq in ['A', '2A', '3A']:
            p = Period('2011', freq=freq)
            self.assertEqual(p + offsets.YearEnd(2), Period('2013', freq=freq))

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(365, 'D'), timedelta(365)]:
                with tm.assertRaises(ValueError):
                    p + o

        for freq in ['M', '2M', '3M']:
            p = Period('2011-03', freq=freq)
            self.assertEqual(p + offsets.MonthEnd(2), Period('2011-05', freq=freq))
            self.assertEqual(p + offsets.MonthEnd(12), Period('2012-03', freq=freq))

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(365, 'D'), timedelta(365)]:
                with tm.assertRaises(ValueError):
                    p + o

        # freq is Tick
        for freq in ['D', '2D', '3D']:
            p = Period('2011-04-01', freq=freq)
            self.assertEqual(p + offsets.Day(5), Period('2011-04-06', freq=freq))
            self.assertEqual(p + offsets.Hour(24), Period('2011-04-02', freq=freq))
            self.assertEqual(p + np.timedelta64(2, 'D'), Period('2011-04-03', freq=freq))
            self.assertEqual(p + np.timedelta64(3600 * 24, 's'), Period('2011-04-02', freq=freq))
            self.assertEqual(p + timedelta(-2), Period('2011-03-30', freq=freq))
            self.assertEqual(p + timedelta(hours=48), Period('2011-04-03', freq=freq))

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(4, 'h'), timedelta(hours=23)]:
                with tm.assertRaises(ValueError):
                    p + o

        for freq in ['H', '2H', '3H']:
            p = Period('2011-04-01 09:00', freq=freq)
            self.assertEqual(p + offsets.Day(2), Period('2011-04-03 09:00', freq=freq))
            self.assertEqual(p + offsets.Hour(3), Period('2011-04-01 12:00', freq=freq))
            self.assertEqual(p + np.timedelta64(3, 'h'), Period('2011-04-01 12:00', freq=freq))
            self.assertEqual(p + np.timedelta64(3600, 's'), Period('2011-04-01 10:00', freq=freq))
            self.assertEqual(p + timedelta(minutes=120), Period('2011-04-01 11:00', freq=freq))
            self.assertEqual(p + timedelta(days=4, minutes=180), Period('2011-04-05 12:00', freq=freq))

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(3200, 's'), timedelta(hours=23, minutes=30)]:
                with tm.assertRaises(ValueError):
                    p + o

    def test_add_offset_nat(self):
        # freq is DateOffset
        for freq in ['A', '2A', '3A']:
            p = Period('NaT', freq=freq)
            for o in [offsets.YearEnd(2)]:
                self.assertEqual((p + o).ordinal, tslib.iNaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(365, 'D'), timedelta(365)]:
                with tm.assertRaises(ValueError):
                    p + o

        for freq in ['M', '2M', '3M']:
            p = Period('NaT', freq=freq)
            for o in [offsets.MonthEnd(2), offsets.MonthEnd(12)]:
                self.assertEqual((p + o).ordinal, tslib.iNaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(365, 'D'), timedelta(365)]:
                with tm.assertRaises(ValueError):
                    p + o

        # freq is Tick
        for freq in ['D', '2D', '3D']:
            p = Period('NaT', freq=freq)
            for o in [offsets.Day(5), offsets.Hour(24), np.timedelta64(2, 'D'),
                      np.timedelta64(3600 * 24, 's'), timedelta(-2), timedelta(hours=48)]:
                self.assertEqual((p + o).ordinal, tslib.iNaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(4, 'h'), timedelta(hours=23)]:
                with tm.assertRaises(ValueError):
                    p + o

        for freq in ['H', '2H', '3H']:
            p = Period('NaT', freq=freq)
            for o in [offsets.Day(2), offsets.Hour(3), np.timedelta64(3, 'h'),
                      np.timedelta64(3600, 's'), timedelta(minutes=120),
                      timedelta(days=4, minutes=180)]:
                self.assertEqual((p + o).ordinal, tslib.iNaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(3200, 's'), timedelta(hours=23, minutes=30)]:
                with tm.assertRaises(ValueError):
                    p + o

    def test_sub_offset(self):
        # freq is DateOffset
        for freq in ['A', '2A', '3A']:
            p = Period('2011', freq=freq)
            self.assertEqual(p - offsets.YearEnd(2), Period('2009', freq=freq))

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(365, 'D'), timedelta(365)]:
                with tm.assertRaises(ValueError):
                    p - o

        for freq in ['M', '2M', '3M']:
            p = Period('2011-03', freq=freq)
            self.assertEqual(p - offsets.MonthEnd(2), Period('2011-01', freq=freq))
            self.assertEqual(p - offsets.MonthEnd(12), Period('2010-03', freq=freq))

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(365, 'D'), timedelta(365)]:
                with tm.assertRaises(ValueError):
                    p - o

        # freq is Tick
        for freq in ['D', '2D', '3D']:
            p = Period('2011-04-01', freq=freq)
            self.assertEqual(p - offsets.Day(5), Period('2011-03-27', freq=freq))
            self.assertEqual(p - offsets.Hour(24), Period('2011-03-31', freq=freq))
            self.assertEqual(p - np.timedelta64(2, 'D'), Period('2011-03-30', freq=freq))
            self.assertEqual(p - np.timedelta64(3600 * 24, 's'), Period('2011-03-31', freq=freq))
            self.assertEqual(p - timedelta(-2), Period('2011-04-03', freq=freq))
            self.assertEqual(p - timedelta(hours=48), Period('2011-03-30', freq=freq))

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(4, 'h'), timedelta(hours=23)]:
                with tm.assertRaises(ValueError):
                    p - o

        for freq in ['H', '2H', '3H']:
            p = Period('2011-04-01 09:00', freq=freq)
            self.assertEqual(p - offsets.Day(2), Period('2011-03-30 09:00', freq=freq))
            self.assertEqual(p - offsets.Hour(3), Period('2011-04-01 06:00', freq=freq))
            self.assertEqual(p - np.timedelta64(3, 'h'), Period('2011-04-01 06:00', freq=freq))
            self.assertEqual(p - np.timedelta64(3600, 's'), Period('2011-04-01 08:00', freq=freq))
            self.assertEqual(p - timedelta(minutes=120), Period('2011-04-01 07:00', freq=freq))
            self.assertEqual(p - timedelta(days=4, minutes=180), Period('2011-03-28 06:00', freq=freq))

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(3200, 's'), timedelta(hours=23, minutes=30)]:
                with tm.assertRaises(ValueError):
                    p - o

    def test_sub_offset_nat(self):
        # freq is DateOffset
        for freq in ['A', '2A', '3A']:
            p = Period('NaT', freq=freq)
            for o in [offsets.YearEnd(2)]:
                self.assertEqual((p - o).ordinal, tslib.iNaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(365, 'D'), timedelta(365)]:
                with tm.assertRaises(ValueError):
                    p - o

        for freq in ['M', '2M', '3M']:
            p = Period('NaT', freq=freq)
            for o in [offsets.MonthEnd(2), offsets.MonthEnd(12)]:
                self.assertEqual((p - o).ordinal, tslib.iNaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(365, 'D'), timedelta(365)]:
                with tm.assertRaises(ValueError):
                    p - o

        # freq is Tick
        for freq in ['D', '2D', '3D']:
            p = Period('NaT', freq=freq)
            for o in [offsets.Day(5), offsets.Hour(24), np.timedelta64(2, 'D'),
                      np.timedelta64(3600 * 24, 's'), timedelta(-2), timedelta(hours=48)]:
                self.assertEqual((p - o).ordinal, tslib.iNaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(4, 'h'), timedelta(hours=23)]:
                with tm.assertRaises(ValueError):
                    p - o

        for freq in ['H', '2H', '3H']:
            p = Period('NaT', freq=freq)
            for o in [offsets.Day(2), offsets.Hour(3), np.timedelta64(3, 'h'),
                      np.timedelta64(3600, 's'), timedelta(minutes=120),
                      timedelta(days=4, minutes=180)]:
                self.assertEqual((p - o).ordinal, tslib.iNaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1), offsets.Minute(),
                      np.timedelta64(3200, 's'), timedelta(hours=23, minutes=30)]:
                with tm.assertRaises(ValueError):
                    p - o

    def test_nat_ops(self):
        for freq in ['M', '2M', '3M']:
            p = Period('NaT', freq=freq)
            self.assertEqual((p + 1).ordinal, tslib.iNaT)
            self.assertEqual((p - 1).ordinal, tslib.iNaT)
            self.assertEqual((p - Period('2011-01', freq=freq)).ordinal, tslib.iNaT)
            self.assertEqual((Period('2011-01', freq=freq) - p).ordinal, tslib.iNaT)

    def test_pi_ops_nat(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'], freq='M', name='idx')
        result = idx + 2
        expected = PeriodIndex(['2011-03', '2011-04', 'NaT', '2011-06'], freq='M', name='idx')
        self.assertTrue(result.equals(expected))

        result2 = result - 2
        self.assertTrue(result2.equals(idx))

        msg = "unsupported operand type\(s\)"
        with tm.assertRaisesRegexp(TypeError, msg):
            idx + "str"

    def test_pi_ops_array(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-04'], freq='M', name='idx')
        result = idx + np.array([1, 2, 3, 4])
        exp = PeriodIndex(['2011-02', '2011-04', 'NaT', '2011-08'], freq='M', name='idx')
        self.assert_index_equal(result, exp)

        result = np.add(idx, np.array([4, -1, 1, 2]))
        exp = PeriodIndex(['2011-05', '2011-01', 'NaT', '2011-06'], freq='M', name='idx')
        self.assert_index_equal(result, exp)

        result = idx - np.array([1, 2, 3, 4])
        exp = PeriodIndex(['2010-12', '2010-12', 'NaT', '2010-12'], freq='M', name='idx')
        self.assert_index_equal(result, exp)

        result = np.subtract(idx, np.array([3, 2, 3, -2]))
        exp = PeriodIndex(['2010-10', '2010-12', 'NaT', '2011-06'], freq='M', name='idx')
        self.assert_index_equal(result, exp)

        # incompatible freq
        msg = "Input has different freq from PeriodIndex\(freq=M\)"
        with tm.assertRaisesRegexp(ValueError, msg):
            idx + np.array([np.timedelta64(1, 'D')] * 4)

        idx = PeriodIndex(['2011-01-01 09:00', '2011-01-01 10:00', 'NaT',
                           '2011-01-01 12:00'], freq='H', name='idx')
        result = idx + np.array([np.timedelta64(1, 'D')] * 4)
        exp = PeriodIndex(['2011-01-02 09:00', '2011-01-02 10:00', 'NaT',
                           '2011-01-02 12:00'], freq='H', name='idx')
        self.assert_index_equal(result, exp)

        result = idx - np.array([np.timedelta64(1, 'h')] * 4)
        exp = PeriodIndex(['2011-01-01 08:00', '2011-01-01 09:00', 'NaT',
                           '2011-01-01 11:00'], freq='H', name='idx')
        self.assert_index_equal(result, exp)

        msg = "Input has different freq from PeriodIndex\(freq=H\)"
        with tm.assertRaisesRegexp(ValueError, msg):
            idx + np.array([np.timedelta64(1, 's')] * 4)

        idx = PeriodIndex(['2011-01-01 09:00:00', '2011-01-01 10:00:00', 'NaT',
                           '2011-01-01 12:00:00'], freq='S', name='idx')
        result = idx + np.array([np.timedelta64(1, 'h'), np.timedelta64(30, 's'),
                                 np.timedelta64(2, 'h'), np.timedelta64(15, 'm')])
        exp = PeriodIndex(['2011-01-01 10:00:00', '2011-01-01 10:00:30', 'NaT',
                           '2011-01-01 12:15:00'], freq='S', name='idx')
        self.assert_index_equal(result, exp)


class TestPeriodRepresentation(tm.TestCase):
    """
    Wish to match NumPy units
    """

    def test_annual(self):
        self._check_freq('A', 1970)

    def test_monthly(self):
        self._check_freq('M', '1970-01')

    def test_weekly(self):
        self._check_freq('W-THU', '1970-01-01')

    def test_daily(self):
        self._check_freq('D', '1970-01-01')

    def test_business_daily(self):
        self._check_freq('B', '1970-01-01')

    def test_hourly(self):
        self._check_freq('H', '1970-01-01')

    def test_minutely(self):
        self._check_freq('T', '1970-01-01')

    def test_secondly(self):
        self._check_freq('S', '1970-01-01')

    def test_millisecondly(self):
        self._check_freq('L', '1970-01-01')

    def test_microsecondly(self):
        self._check_freq('U', '1970-01-01')

    def test_nanosecondly(self):
        self._check_freq('N', '1970-01-01')

    def _check_freq(self, freq, base_date):
        rng = PeriodIndex(start=base_date, periods=10, freq=freq)
        exp = np.arange(10, dtype=np.int64)
        self.assert_numpy_array_equal(rng.values, exp)

    def test_negone_ordinals(self):
        freqs = ['A', 'M', 'Q', 'D', 'H', 'T', 'S']

        period = Period(ordinal=-1, freq='D')
        for freq in freqs:
            repr(period.asfreq(freq))

        for freq in freqs:
            period = Period(ordinal=-1, freq=freq)
            repr(period)
            self.assertEqual(period.year, 1969)

        period = Period(ordinal=-1, freq='B')
        repr(period)
        period = Period(ordinal=-1, freq='W')
        repr(period)


class TestComparisons(tm.TestCase):
    def setUp(self):
        self.january1 = Period('2000-01', 'M')
        self.january2 = Period('2000-01', 'M')
        self.february = Period('2000-02', 'M')
        self.march = Period('2000-03', 'M')
        self.day = Period('2012-01-01', 'D')

    def test_equal(self):
        self.assertEqual(self.january1, self.january2)

    def test_equal_Raises_Value(self):
        with tm.assertRaises(ValueError):
            self.january1 == self.day

    def test_notEqual(self):
        self.assertNotEqual(self.january1, 1)
        self.assertNotEqual(self.january1, self.february)

    def test_greater(self):
        self.assertTrue(self.february > self.january1)

    def test_greater_Raises_Value(self):
        with tm.assertRaises(ValueError):
            self.january1 > self.day

    def test_greater_Raises_Type(self):
        with tm.assertRaises(TypeError):
            self.january1 > 1

    def test_greaterEqual(self):
        self.assertTrue(self.january1 >= self.january2)

    def test_greaterEqual_Raises_Value(self):
        with tm.assertRaises(ValueError):
            self.january1 >= self.day
        with tm.assertRaises(TypeError):
            print(self.january1 >= 1)

    def test_smallerEqual(self):
        self.assertTrue(self.january1 <= self.january2)

    def test_smallerEqual_Raises_Value(self):
        with tm.assertRaises(ValueError):
            self.january1 <= self.day

    def test_smallerEqual_Raises_Type(self):
        with tm.assertRaises(TypeError):
            self.january1 <= 1

    def test_smaller(self):
        self.assertTrue(self.january1 < self.february)

    def test_smaller_Raises_Value(self):
        with tm.assertRaises(ValueError):
            self.january1 < self.day

    def test_smaller_Raises_Type(self):
        with tm.assertRaises(TypeError):
            self.january1 < 1

    def test_sort(self):
        periods = [self.march, self.january1, self.february]
        correctPeriods = [self.january1, self.february, self.march]
        self.assertEqual(sorted(periods), correctPeriods)

    def test_period_nat_comp(self):
        p_nat = Period('NaT', freq='D')
        p = Period('2011-01-01', freq='D')

        nat = pd.Timestamp('NaT')
        t = pd.Timestamp('2011-01-01')
        # confirm Period('NaT') work identical with Timestamp('NaT')
        for left, right in [(p_nat, p), (p, p_nat), (p_nat, p_nat),
                            (nat, t), (t, nat), (nat, nat)]:
            self.assertEqual(left < right, False)
            self.assertEqual(left > right, False)
            self.assertEqual(left == right, False)
            self.assertEqual(left != right, True)
            self.assertEqual(left <= right, False)
            self.assertEqual(left >= right, False)

    def test_pi_pi_comp(self):

        for freq in ['M', '2M', '3M']:
            base = PeriodIndex(['2011-01', '2011-02',
                                '2011-03', '2011-04'], freq=freq)
            p = Period('2011-02', freq=freq)

            exp = np.array([False, True, False, False])
            self.assert_numpy_array_equal(base == p, exp)

            exp = np.array([True, False, True, True])
            self.assert_numpy_array_equal(base != p, exp)

            exp = np.array([False, False, True, True])
            self.assert_numpy_array_equal(base > p, exp)

            exp = np.array([True, False, False, False])
            self.assert_numpy_array_equal(base < p, exp)

            exp = np.array([False, True, True, True])
            self.assert_numpy_array_equal(base >= p, exp)

            exp = np.array([True, True, False, False])
            self.assert_numpy_array_equal(base <= p, exp)

            idx = PeriodIndex(['2011-02', '2011-01', '2011-03', '2011-05'], freq=freq)

            exp = np.array([False, False, True, False])
            self.assert_numpy_array_equal(base == idx, exp)

            exp = np.array([True, True, False, True])
            self.assert_numpy_array_equal(base != idx, exp)

            exp = np.array([False, True, False, False])
            self.assert_numpy_array_equal(base > idx, exp)

            exp = np.array([True, False, False, True])
            self.assert_numpy_array_equal(base < idx, exp)

            exp = np.array([False, True, True, False])
            self.assert_numpy_array_equal(base >= idx, exp)

            exp = np.array([True, False, True, True])
            self.assert_numpy_array_equal(base <= idx, exp)

            # different base freq
            msg = "Input has different freq=A-DEC from PeriodIndex"
            with tm.assertRaisesRegexp(ValueError, msg):
                base <= Period('2011', freq='A')

            with tm.assertRaisesRegexp(ValueError, msg):
                idx = PeriodIndex(['2011', '2012', '2013', '2014'], freq='A')
                base <= idx

            # different mult
            msg = "Input has different freq=4M from PeriodIndex"
            with tm.assertRaisesRegexp(ValueError, msg):
                base <= Period('2011', freq='4M')

            with tm.assertRaisesRegexp(ValueError, msg):
                idx = PeriodIndex(['2011', '2012', '2013', '2014'], freq='4M')
                base <= idx

    def test_pi_nat_comp(self):
        for freq in ['M', '2M', '3M']:
            idx1 = PeriodIndex(['2011-01', '2011-02', 'NaT', '2011-05'], freq=freq)

            result = idx1 > Period('2011-02', freq=freq)
            exp = np.array([False, False, False, True])
            self.assert_numpy_array_equal(result, exp)

            result = idx1 == Period('NaT', freq=freq)
            exp = np.array([False, False, False, False])
            self.assert_numpy_array_equal(result, exp)

            result = idx1 != Period('NaT', freq=freq)
            exp = np.array([True, True, True, True])
            self.assert_numpy_array_equal(result, exp)

            idx2 = PeriodIndex(['2011-02', '2011-01', '2011-04', 'NaT'], freq=freq)
            result = idx1 < idx2
            exp = np.array([True, False, False, False])
            self.assert_numpy_array_equal(result, exp)

            result = idx1 == idx2
            exp = np.array([False, False, False, False])
            self.assert_numpy_array_equal(result, exp)

            result = idx1 != idx2
            exp = np.array([True, True, True, True])
            self.assert_numpy_array_equal(result, exp)

            result = idx1 == idx1
            exp = np.array([True, True, False, True])
            self.assert_numpy_array_equal(result, exp)

            result = idx1 != idx1
            exp = np.array([False, False, True, False])
            self.assert_numpy_array_equal(result, exp)

            diff = PeriodIndex(['2011-02', '2011-01', '2011-04', 'NaT'], freq='4M')
            msg = "Input has different freq=4M from PeriodIndex"
            with tm.assertRaisesRegexp(ValueError, msg):
                idx1 > diff
            with tm.assertRaisesRegexp(ValueError, msg):
                idx1 == diff


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
