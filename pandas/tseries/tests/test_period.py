"""Tests suite for Period handling.

Parts derived from scikits.timeseries code, original authors:
- Pierre Gerard-Marchant & Matt Knox
- pierregm_at_uga_dot_edu - mattknow_ca_at_hotmail_dot_com

"""

from datetime import datetime, date, timedelta

from pandas import Timestamp, _period
from pandas.tseries.frequencies import MONTHS, DAYS, _period_code_map
from pandas.tseries.period import Period, PeriodIndex, period_range
from pandas.tseries.index import DatetimeIndex, date_range, Index
from pandas.tseries.tools import to_datetime
import pandas.tseries.period as period
import pandas.tseries.offsets as offsets

import pandas as pd
import numpy as np
from numpy.random import randn
from pandas.compat import range, lrange, lmap, zip, text_type, PY3, iteritems
from pandas.compat.numpy import np_datetime64_compat

from pandas import (Series, DataFrame,
                    _np_version_under1p9, _np_version_under1p10,
                    _np_version_under1p12)
from pandas import tslib
import pandas.util.testing as tm


class TestPeriodProperties(tm.TestCase):
    "Test properties such as year, month, weekday, etc...."

    def test_quarterly_negative_ordinals(self):
        p = Period(ordinal=-1, freq='Q-DEC')
        self.assertEqual(p.year, 1969)
        self.assertEqual(p.quarter, 4)
        self.assertIsInstance(p, Period)

        p = Period(ordinal=-2, freq='Q-DEC')
        self.assertEqual(p.year, 1969)
        self.assertEqual(p.quarter, 3)
        self.assertIsInstance(p, Period)

        p = Period(ordinal=-2, freq='M')
        self.assertEqual(p.year, 1969)
        self.assertEqual(p.month, 11)
        self.assertIsInstance(p, Period)

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
            self.assertIsInstance(p, Period)

    def test_period_cons_weekly(self):
        for num in range(10, 17):
            daystr = '2011-02-%d' % num
            for day in DAYS:
                freq = 'W-%s' % day

                result = Period(daystr, freq=freq)
                expected = Period(daystr, freq='D').asfreq(freq)
                self.assertEqual(result, expected)
                self.assertIsInstance(result, Period)

    def test_period_from_ordinal(self):
        p = pd.Period('2011-01', freq='M')
        res = pd.Period._from_ordinal(p.ordinal, freq='M')
        self.assertEqual(p, res)
        self.assertIsInstance(res, Period)

    def test_period_cons_nat(self):
        p = Period('NaT', freq='M')
        self.assertIs(p, pd.NaT)

        p = Period('nat', freq='W-SUN')
        self.assertIs(p, pd.NaT)

        p = Period(tslib.iNaT, freq='D')
        self.assertIs(p, pd.NaT)

        p = Period(tslib.iNaT, freq='3D')
        self.assertIs(p, pd.NaT)

        p = Period(tslib.iNaT, freq='1D1H')
        self.assertIs(p, pd.NaT)

        p = Period('NaT')
        self.assertIs(p, pd.NaT)

        p = Period(tslib.iNaT)
        self.assertIs(p, pd.NaT)

    def test_cons_null_like(self):
        # check Timestamp compat
        self.assertIs(Timestamp('NaT'), pd.NaT)
        self.assertIs(Period('NaT'), pd.NaT)

        self.assertIs(Timestamp(None), pd.NaT)
        self.assertIs(Period(None), pd.NaT)

        self.assertIs(Timestamp(float('nan')), pd.NaT)
        self.assertIs(Period(float('nan')), pd.NaT)

        self.assertIs(Timestamp(np.nan), pd.NaT)
        self.assertIs(Period(np.nan), pd.NaT)

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

        msg = ('Frequency must be positive, because it' ' represents span: 0M')
        with tm.assertRaisesRegexp(ValueError, msg):
            Period('2011-01', freq='0M')

    def test_period_cons_combined(self):
        p = [(Period('2011-01', freq='1D1H'),
              Period('2011-01', freq='1H1D'),
              Period('2011-01', freq='H')),
             (Period(ordinal=1, freq='1D1H'),
              Period(ordinal=1, freq='1H1D'),
              Period(ordinal=1, freq='H'))]

        for p1, p2, p3 in p:
            self.assertEqual(p1.ordinal, p3.ordinal)
            self.assertEqual(p2.ordinal, p3.ordinal)

            self.assertEqual(p1.freq, offsets.Hour(25))
            self.assertEqual(p1.freqstr, '25H')

            self.assertEqual(p2.freq, offsets.Hour(25))
            self.assertEqual(p2.freqstr, '25H')

            self.assertEqual(p3.freq, offsets.Hour())
            self.assertEqual(p3.freqstr, 'H')

            result = p1 + 1
            self.assertEqual(result.ordinal, (p3 + 25).ordinal)
            self.assertEqual(result.freq, p1.freq)
            self.assertEqual(result.freqstr, '25H')

            result = p2 + 1
            self.assertEqual(result.ordinal, (p3 + 25).ordinal)
            self.assertEqual(result.freq, p2.freq)
            self.assertEqual(result.freqstr, '25H')

            result = p1 - 1
            self.assertEqual(result.ordinal, (p3 - 25).ordinal)
            self.assertEqual(result.freq, p1.freq)
            self.assertEqual(result.freqstr, '25H')

            result = p2 - 1
            self.assertEqual(result.ordinal, (p3 - 25).ordinal)
            self.assertEqual(result.freq, p2.freq)
            self.assertEqual(result.freqstr, '25H')

        msg = ('Frequency must be positive, because it'
               ' represents span: -25H')
        with tm.assertRaisesRegexp(ValueError, msg):
            Period('2011-01', freq='-1D1H')
        with tm.assertRaisesRegexp(ValueError, msg):
            Period('2011-01', freq='-1H1D')
        with tm.assertRaisesRegexp(ValueError, msg):
            Period(ordinal=1, freq='-1D1H')
        with tm.assertRaisesRegexp(ValueError, msg):
            Period(ordinal=1, freq='-1H1D')

        msg = ('Frequency must be positive, because it'
               ' represents span: 0D')
        with tm.assertRaisesRegexp(ValueError, msg):
            Period('2011-01', freq='0D0H')
        with tm.assertRaisesRegexp(ValueError, msg):
            Period(ordinal=1, freq='0D0H')

        # You can only combine together day and intraday offsets
        msg = ('Invalid frequency: 1W1D')
        with tm.assertRaisesRegexp(ValueError, msg):
            Period('2011-01', freq='1W1D')
        msg = ('Invalid frequency: 1D1W')
        with tm.assertRaisesRegexp(ValueError, msg):
            Period('2011-01', freq='1D1W')

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
            p = Period('1/1/2005', freq='M').to_timestamp(
                tz=maybe_get_tz(case))
            exp = Timestamp('1/1/2005', tz='UTC').tz_convert(case)
            self.assertEqual(p, exp)
            self.assertEqual(p.tz, gettz(case.split('/', 1)[1]))
            self.assertEqual(p.tz, exp.tz)

            p = Period('1/1/2005',
                       freq='M').to_timestamp(freq='3H', tz=maybe_get_tz(case))
            exp = Timestamp('1/1/2005', tz='UTC').tz_convert(case)
            self.assertEqual(p, exp)
            self.assertEqual(p.tz, gettz(case.split('/', 1)[1]))
            self.assertEqual(p.tz, exp.tz)

    def test_timestamp_tz_arg_dateutil_from_string(self):
        from pandas.tslib import _dateutil_gettz as gettz
        p = Period('1/1/2005',
                   freq='M').to_timestamp(tz='dateutil/Europe/Brussels')
        self.assertEqual(p.tz, gettz('Europe/Brussels'))

    def test_timestamp_mult(self):
        p = pd.Period('2011-01', freq='M')
        self.assertEqual(p.to_timestamp(how='S'), pd.Timestamp('2011-01-01'))
        self.assertEqual(p.to_timestamp(how='E'), pd.Timestamp('2011-01-31'))

        p = pd.Period('2011-01', freq='3M')
        self.assertEqual(p.to_timestamp(how='S'), pd.Timestamp('2011-01-01'))
        self.assertEqual(p.to_timestamp(how='E'), pd.Timestamp('2011-03-31'))

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
        i4 = Period(np_datetime64_compat('2007-01-01 00:00:00Z'), freq='M')
        i5 = Period(np_datetime64_compat('2007-01-01 00:00:00.000Z'), freq='M')
        self.assertEqual(i1, i2)
        self.assertEqual(i1, i3)
        self.assertEqual(i1, i4)
        self.assertEqual(i1, i5)

        i1 = Period('2007-01-01 09:00:00.001')
        expected = Period(datetime(2007, 1, 1, 9, 0, 0, 1000), freq='L')
        self.assertEqual(i1, expected)

        expected = Period(np_datetime64_compat(
            '2007-01-01 09:00:00.001Z'), freq='L')
        self.assertEqual(i1, expected)

        i1 = Period('2007-01-01 09:00:00.00101')
        expected = Period(datetime(2007, 1, 1, 9, 0, 0, 1010), freq='U')
        self.assertEqual(i1, expected)

        expected = Period(np_datetime64_compat('2007-01-01 09:00:00.00101Z'),
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
        self.assertEqual(Period(year=2012, month=3, day=10,
                                freq=offsets.BDay()),
                         Period(year=2012, month=3, day=10, freq='B'))

        expected = Period('2005-03-01', freq='3D')
        self.assertEqual(Period(year=2005, month=3, day=1,
                                freq=offsets.Day(3)), expected)
        self.assertEqual(Period(year=2005, month=3, day=1, freq='3D'),
                         expected)

        self.assertEqual(Period(year=2012, month=3, day=10,
                                freq=offsets.BDay(3)),
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
        i4 = Period(np_datetime64_compat('2007-01-01 00:00:00Z'), freq='M')
        i5 = Period(np_datetime64_compat('2007-01-01 00:00:00.000Z'), freq='M')
        self.assertEqual(i1, i2)
        self.assertEqual(i1, i3)
        self.assertEqual(i1, i4)
        self.assertEqual(i1, i5)

        i1 = Period('2007-01-01 09:00:00.001')
        expected = Period(datetime(2007, 1, 1, 9, 0, 0, 1000), freq='L')
        self.assertEqual(i1, expected)

        expected = Period(np_datetime64_compat(
            '2007-01-01 09:00:00.001Z'), freq='L')
        self.assertEqual(i1, expected)

        i1 = Period('2007-01-01 09:00:00.00101')
        expected = Period(datetime(2007, 1, 1, 9, 0, 0, 1010), freq='U')
        self.assertEqual(i1, expected)

        expected = Period(np_datetime64_compat('2007-01-01 09:00:00.00101Z'),
                          freq='U')
        self.assertEqual(i1, expected)

        self.assertRaises(ValueError, Period, ordinal=200701)

        self.assertRaises(ValueError, Period, '2007-1-1', freq='X')

    def test_freq_str(self):
        i1 = Period('1982', freq='Min')
        self.assertEqual(i1.freq, offsets.Minute())
        self.assertEqual(i1.freqstr, 'T')

    def test_period_deprecated_freq(self):
        cases = {"M": ["MTH", "MONTH", "MONTHLY", "Mth", "month", "monthly"],
                 "B": ["BUS", "BUSINESS", "BUSINESSLY", "WEEKDAY", "bus"],
                 "D": ["DAY", "DLY", "DAILY", "Day", "Dly", "Daily"],
                 "H": ["HR", "HOUR", "HRLY", "HOURLY", "hr", "Hour", "HRly"],
                 "T": ["minute", "MINUTE", "MINUTELY", "minutely"],
                 "S": ["sec", "SEC", "SECOND", "SECONDLY", "second"],
                 "L": ["MILLISECOND", "MILLISECONDLY", "millisecond"],
                 "U": ["MICROSECOND", "MICROSECONDLY", "microsecond"],
                 "N": ["NANOSECOND", "NANOSECONDLY", "nanosecond"]}

        msg = pd.tseries.frequencies._INVALID_FREQ_ERROR
        for exp, freqs in iteritems(cases):
            for freq in freqs:
                with self.assertRaisesRegexp(ValueError, msg):
                    Period('2016-03-01 09:00', freq=freq)
                with self.assertRaisesRegexp(ValueError, msg):
                    Period(ordinal=1, freq=freq)

            # check supported freq-aliases still works
            p1 = Period('2016-03-01 09:00', freq=exp)
            p2 = Period(ordinal=1, freq=exp)
            tm.assertIsInstance(p1, Period)
            tm.assertIsInstance(p2, Period)

    def test_hash(self):
        self.assertEqual(hash(Period('2011-01', freq='M')),
                         hash(Period('2011-01', freq='M')))

        self.assertNotEqual(hash(Period('2011-01-01', freq='D')),
                            hash(Period('2011-01', freq='M')))

        self.assertNotEqual(hash(Period('2011-01', freq='3M')),
                            hash(Period('2011-01', freq='2M')))

        self.assertNotEqual(hash(Period('2011-01', freq='M')),
                            hash(Period('2011-02', freq='M')))

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
        self.assertEqual(res, '2000-01-01 12:34:12')
        tm.assertIsInstance(res, text_type)  # GH3363

    def test_sub_delta(self):
        left, right = Period('2011', freq='A'), Period('2007', freq='A')
        result = left - right
        self.assertEqual(result, 4)

        with self.assertRaises(period.IncompatibleFrequency):
            left - Period('2007-01', freq='M')

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

        from_lst = ['A', 'Q', 'M', 'W', 'B', 'D', 'H', 'Min', 'S']

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

        p = Period('2012', freq='D')
        xp = _ex(2012, 1, 2)
        self.assertEqual(xp, p.end_time)

        p = Period('2012', freq='H')
        xp = _ex(2012, 1, 1, 1)
        self.assertEqual(xp, p.end_time)

        p = Period('2012', freq='B')
        xp = _ex(2012, 1, 3)
        self.assertEqual(xp, p.end_time)

        p = Period('2012', freq='W')
        xp = _ex(2012, 1, 2)
        self.assertEqual(xp, p.end_time)

        # Test for GH 11738
        p = Period('2012', freq='15D')
        xp = _ex(2012, 1, 16)
        self.assertEqual(xp, p.end_time)

        p = Period('2012', freq='1D1H')
        xp = _ex(2012, 1, 2, 1)
        self.assertEqual(xp, p.end_time)

        p = Period('2012', freq='1H1D')
        xp = _ex(2012, 1, 2, 1)
        self.assertEqual(xp, p.end_time)

    def test_anchor_week_end_time(self):
        def _ex(*args):
            return Timestamp(Timestamp(datetime(*args)).value - 1)

        p = Period('2013-1-1', 'W-SAT')
        xp = _ex(2013, 1, 6)
        self.assertEqual(p.end_time, xp)

    def test_properties_annually(self):
        # Test properties on Periods with annually frequency.
        a_date = Period(freq='A', year=2007)
        self.assertEqual(a_date.year, 2007)

    def test_properties_quarterly(self):
        # Test properties on Periods with daily frequency.
        qedec_date = Period(freq="Q-DEC", year=2007, quarter=1)
        qejan_date = Period(freq="Q-JAN", year=2007, quarter=1)
        qejun_date = Period(freq="Q-JUN", year=2007, quarter=1)
        #
        for x in range(3):
            for qd in (qedec_date, qejan_date, qejun_date):
                self.assertEqual((qd + x).qyear, 2007)
                self.assertEqual((qd + x).quarter, x + 1)

    def test_properties_monthly(self):
        # Test properties on Periods with daily frequency.
        m_date = Period(freq='M', year=2007, month=1)
        for x in range(11):
            m_ival_x = m_date + x
            self.assertEqual(m_ival_x.year, 2007)
            if 1 <= x + 1 <= 3:
                self.assertEqual(m_ival_x.quarter, 1)
            elif 4 <= x + 1 <= 6:
                self.assertEqual(m_ival_x.quarter, 2)
            elif 7 <= x + 1 <= 9:
                self.assertEqual(m_ival_x.quarter, 3)
            elif 10 <= x + 1 <= 12:
                self.assertEqual(m_ival_x.quarter, 4)
            self.assertEqual(m_ival_x.month, x + 1)

    def test_properties_weekly(self):
        # Test properties on Periods with daily frequency.
        w_date = Period(freq='W', year=2007, month=1, day=7)
        #
        self.assertEqual(w_date.year, 2007)
        self.assertEqual(w_date.quarter, 1)
        self.assertEqual(w_date.month, 1)
        self.assertEqual(w_date.week, 1)
        self.assertEqual((w_date - 1).week, 52)
        self.assertEqual(w_date.days_in_month, 31)
        self.assertEqual(Period(freq='W', year=2012,
                                month=2, day=1).days_in_month, 29)

    def test_properties_weekly_legacy(self):
        # Test properties on Periods with daily frequency.
        w_date = Period(freq='W', year=2007, month=1, day=7)
        self.assertEqual(w_date.year, 2007)
        self.assertEqual(w_date.quarter, 1)
        self.assertEqual(w_date.month, 1)
        self.assertEqual(w_date.week, 1)
        self.assertEqual((w_date - 1).week, 52)
        self.assertEqual(w_date.days_in_month, 31)

        exp = Period(freq='W', year=2012, month=2, day=1)
        self.assertEqual(exp.days_in_month, 29)

        msg = pd.tseries.frequencies._INVALID_FREQ_ERROR
        with self.assertRaisesRegexp(ValueError, msg):
            Period(freq='WK', year=2007, month=1, day=7)

    def test_properties_daily(self):
        # Test properties on Periods with daily frequency.
        b_date = Period(freq='B', year=2007, month=1, day=1)
        #
        self.assertEqual(b_date.year, 2007)
        self.assertEqual(b_date.quarter, 1)
        self.assertEqual(b_date.month, 1)
        self.assertEqual(b_date.day, 1)
        self.assertEqual(b_date.weekday, 0)
        self.assertEqual(b_date.dayofyear, 1)
        self.assertEqual(b_date.days_in_month, 31)
        self.assertEqual(Period(freq='B', year=2012,
                                month=2, day=1).days_in_month, 29)
        #
        d_date = Period(freq='D', year=2007, month=1, day=1)
        #
        self.assertEqual(d_date.year, 2007)
        self.assertEqual(d_date.quarter, 1)
        self.assertEqual(d_date.month, 1)
        self.assertEqual(d_date.day, 1)
        self.assertEqual(d_date.weekday, 0)
        self.assertEqual(d_date.dayofyear, 1)
        self.assertEqual(d_date.days_in_month, 31)
        self.assertEqual(Period(freq='D', year=2012, month=2,
                                day=1).days_in_month, 29)

    def test_properties_hourly(self):
        # Test properties on Periods with hourly frequency.
        h_date1 = Period(freq='H', year=2007, month=1, day=1, hour=0)
        h_date2 = Period(freq='2H', year=2007, month=1, day=1, hour=0)

        for h_date in [h_date1, h_date2]:
            self.assertEqual(h_date.year, 2007)
            self.assertEqual(h_date.quarter, 1)
            self.assertEqual(h_date.month, 1)
            self.assertEqual(h_date.day, 1)
            self.assertEqual(h_date.weekday, 0)
            self.assertEqual(h_date.dayofyear, 1)
            self.assertEqual(h_date.hour, 0)
            self.assertEqual(h_date.days_in_month, 31)
            self.assertEqual(Period(freq='H', year=2012, month=2, day=1,
                                    hour=0).days_in_month, 29)

    def test_properties_minutely(self):
        # Test properties on Periods with minutely frequency.
        t_date = Period(freq='Min', year=2007, month=1, day=1, hour=0,
                        minute=0)
        #
        self.assertEqual(t_date.quarter, 1)
        self.assertEqual(t_date.month, 1)
        self.assertEqual(t_date.day, 1)
        self.assertEqual(t_date.weekday, 0)
        self.assertEqual(t_date.dayofyear, 1)
        self.assertEqual(t_date.hour, 0)
        self.assertEqual(t_date.minute, 0)
        self.assertEqual(t_date.days_in_month, 31)
        self.assertEqual(Period(freq='D', year=2012, month=2, day=1, hour=0,
                                minute=0).days_in_month, 29)

    def test_properties_secondly(self):
        # Test properties on Periods with secondly frequency.
        s_date = Period(freq='Min', year=2007, month=1, day=1, hour=0,
                        minute=0, second=0)
        #
        self.assertEqual(s_date.year, 2007)
        self.assertEqual(s_date.quarter, 1)
        self.assertEqual(s_date.month, 1)
        self.assertEqual(s_date.day, 1)
        self.assertEqual(s_date.weekday, 0)
        self.assertEqual(s_date.dayofyear, 1)
        self.assertEqual(s_date.hour, 0)
        self.assertEqual(s_date.minute, 0)
        self.assertEqual(s_date.second, 0)
        self.assertEqual(s_date.days_in_month, 31)
        self.assertEqual(Period(freq='Min', year=2012, month=2, day=1, hour=0,
                                minute=0, second=0).days_in_month, 29)

    def test_properties_nat(self):
        p_nat = Period('NaT', freq='M')
        t_nat = pd.Timestamp('NaT')
        self.assertIs(p_nat, t_nat)

        # confirm Period('NaT') work identical with Timestamp('NaT')
        for f in ['year', 'month', 'day', 'hour', 'minute', 'second', 'week',
                  'dayofyear', 'quarter', 'days_in_month']:
            self.assertTrue(np.isnan(getattr(p_nat, f)))
            self.assertTrue(np.isnan(getattr(t_nat, f)))

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
        self.assertIs(Period(None), pd.NaT)
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

        self.assertEqual(initial.asfreq(freq="M", how="S"),
                         Period('2013-01', 'M'))

        msg = pd.tseries.frequencies._INVALID_FREQ_ERROR
        with self.assertRaisesRegexp(ValueError, msg):
            initial.asfreq(freq="MS", how="S")

        with tm.assertRaisesRegexp(ValueError, msg):
            pd.Period('2013-01', 'MS')

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
        ival_A_to_H_start = Period(freq='H', year=2007, month=1, day=1, hour=0)
        ival_A_to_H_end = Period(freq='H', year=2007, month=12, day=31,
                                 hour=23)
        ival_A_to_T_start = Period(freq='Min', year=2007, month=1, day=1,
                                   hour=0, minute=0)
        ival_A_to_T_end = Period(freq='Min', year=2007, month=12, day=31,
                                 hour=23, minute=59)
        ival_A_to_S_start = Period(freq='S', year=2007, month=1, day=1, hour=0,
                                   minute=0, second=0)
        ival_A_to_S_end = Period(freq='S', year=2007, month=12, day=31,
                                 hour=23, minute=59, second=59)

        ival_AJAN_to_D_end = Period(freq='D', year=2007, month=1, day=31)
        ival_AJAN_to_D_start = Period(freq='D', year=2006, month=2, day=1)
        ival_AJUN_to_D_end = Period(freq='D', year=2007, month=6, day=30)
        ival_AJUN_to_D_start = Period(freq='D', year=2006, month=7, day=1)
        ival_ANOV_to_D_end = Period(freq='D', year=2007, month=11, day=30)
        ival_ANOV_to_D_start = Period(freq='D', year=2006, month=12, day=1)

        self.assertEqual(ival_A.asfreq('Q', 'S'), ival_A_to_Q_start)
        self.assertEqual(ival_A.asfreq('Q', 'e'), ival_A_to_Q_end)
        self.assertEqual(ival_A.asfreq('M', 's'), ival_A_to_M_start)
        self.assertEqual(ival_A.asfreq('M', 'E'), ival_A_to_M_end)
        self.assertEqual(ival_A.asfreq('W', 'S'), ival_A_to_W_start)
        self.assertEqual(ival_A.asfreq('W', 'E'), ival_A_to_W_end)
        self.assertEqual(ival_A.asfreq('B', 'S'), ival_A_to_B_start)
        self.assertEqual(ival_A.asfreq('B', 'E'), ival_A_to_B_end)
        self.assertEqual(ival_A.asfreq('D', 'S'), ival_A_to_D_start)
        self.assertEqual(ival_A.asfreq('D', 'E'), ival_A_to_D_end)
        self.assertEqual(ival_A.asfreq('H', 'S'), ival_A_to_H_start)
        self.assertEqual(ival_A.asfreq('H', 'E'), ival_A_to_H_end)
        self.assertEqual(ival_A.asfreq('min', 'S'), ival_A_to_T_start)
        self.assertEqual(ival_A.asfreq('min', 'E'), ival_A_to_T_end)
        self.assertEqual(ival_A.asfreq('T', 'S'), ival_A_to_T_start)
        self.assertEqual(ival_A.asfreq('T', 'E'), ival_A_to_T_end)
        self.assertEqual(ival_A.asfreq('S', 'S'), ival_A_to_S_start)
        self.assertEqual(ival_A.asfreq('S', 'E'), ival_A_to_S_end)

        self.assertEqual(ival_AJAN.asfreq('D', 'S'), ival_AJAN_to_D_start)
        self.assertEqual(ival_AJAN.asfreq('D', 'E'), ival_AJAN_to_D_end)

        self.assertEqual(ival_AJUN.asfreq('D', 'S'), ival_AJUN_to_D_start)
        self.assertEqual(ival_AJUN.asfreq('D', 'E'), ival_AJUN_to_D_end)

        self.assertEqual(ival_ANOV.asfreq('D', 'S'), ival_ANOV_to_D_start)
        self.assertEqual(ival_ANOV.asfreq('D', 'E'), ival_ANOV_to_D_end)

        self.assertEqual(ival_A.asfreq('A'), ival_A)

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
        ival_Q_to_H_start = Period(freq='H', year=2007, month=1, day=1, hour=0)
        ival_Q_to_H_end = Period(freq='H', year=2007, month=3, day=31, hour=23)
        ival_Q_to_T_start = Period(freq='Min', year=2007, month=1, day=1,
                                   hour=0, minute=0)
        ival_Q_to_T_end = Period(freq='Min', year=2007, month=3, day=31,
                                 hour=23, minute=59)
        ival_Q_to_S_start = Period(freq='S', year=2007, month=1, day=1, hour=0,
                                   minute=0, second=0)
        ival_Q_to_S_end = Period(freq='S', year=2007, month=3, day=31, hour=23,
                                 minute=59, second=59)

        ival_QEJAN_to_D_start = Period(freq='D', year=2006, month=2, day=1)
        ival_QEJAN_to_D_end = Period(freq='D', year=2006, month=4, day=30)

        ival_QEJUN_to_D_start = Period(freq='D', year=2006, month=7, day=1)
        ival_QEJUN_to_D_end = Period(freq='D', year=2006, month=9, day=30)

        self.assertEqual(ival_Q.asfreq('A'), ival_Q_to_A)
        self.assertEqual(ival_Q_end_of_year.asfreq('A'), ival_Q_to_A)

        self.assertEqual(ival_Q.asfreq('M', 'S'), ival_Q_to_M_start)
        self.assertEqual(ival_Q.asfreq('M', 'E'), ival_Q_to_M_end)
        self.assertEqual(ival_Q.asfreq('W', 'S'), ival_Q_to_W_start)
        self.assertEqual(ival_Q.asfreq('W', 'E'), ival_Q_to_W_end)
        self.assertEqual(ival_Q.asfreq('B', 'S'), ival_Q_to_B_start)
        self.assertEqual(ival_Q.asfreq('B', 'E'), ival_Q_to_B_end)
        self.assertEqual(ival_Q.asfreq('D', 'S'), ival_Q_to_D_start)
        self.assertEqual(ival_Q.asfreq('D', 'E'), ival_Q_to_D_end)
        self.assertEqual(ival_Q.asfreq('H', 'S'), ival_Q_to_H_start)
        self.assertEqual(ival_Q.asfreq('H', 'E'), ival_Q_to_H_end)
        self.assertEqual(ival_Q.asfreq('Min', 'S'), ival_Q_to_T_start)
        self.assertEqual(ival_Q.asfreq('Min', 'E'), ival_Q_to_T_end)
        self.assertEqual(ival_Q.asfreq('S', 'S'), ival_Q_to_S_start)
        self.assertEqual(ival_Q.asfreq('S', 'E'), ival_Q_to_S_end)

        self.assertEqual(ival_QEJAN.asfreq('D', 'S'), ival_QEJAN_to_D_start)
        self.assertEqual(ival_QEJAN.asfreq('D', 'E'), ival_QEJAN_to_D_end)
        self.assertEqual(ival_QEJUN.asfreq('D', 'S'), ival_QEJUN_to_D_start)
        self.assertEqual(ival_QEJUN.asfreq('D', 'E'), ival_QEJUN_to_D_end)

        self.assertEqual(ival_Q.asfreq('Q'), ival_Q)

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
        ival_M_to_H_start = Period(freq='H', year=2007, month=1, day=1, hour=0)
        ival_M_to_H_end = Period(freq='H', year=2007, month=1, day=31, hour=23)
        ival_M_to_T_start = Period(freq='Min', year=2007, month=1, day=1,
                                   hour=0, minute=0)
        ival_M_to_T_end = Period(freq='Min', year=2007, month=1, day=31,
                                 hour=23, minute=59)
        ival_M_to_S_start = Period(freq='S', year=2007, month=1, day=1, hour=0,
                                   minute=0, second=0)
        ival_M_to_S_end = Period(freq='S', year=2007, month=1, day=31, hour=23,
                                 minute=59, second=59)

        self.assertEqual(ival_M.asfreq('A'), ival_M_to_A)
        self.assertEqual(ival_M_end_of_year.asfreq('A'), ival_M_to_A)
        self.assertEqual(ival_M.asfreq('Q'), ival_M_to_Q)
        self.assertEqual(ival_M_end_of_quarter.asfreq('Q'), ival_M_to_Q)

        self.assertEqual(ival_M.asfreq('W', 'S'), ival_M_to_W_start)
        self.assertEqual(ival_M.asfreq('W', 'E'), ival_M_to_W_end)
        self.assertEqual(ival_M.asfreq('B', 'S'), ival_M_to_B_start)
        self.assertEqual(ival_M.asfreq('B', 'E'), ival_M_to_B_end)
        self.assertEqual(ival_M.asfreq('D', 'S'), ival_M_to_D_start)
        self.assertEqual(ival_M.asfreq('D', 'E'), ival_M_to_D_end)
        self.assertEqual(ival_M.asfreq('H', 'S'), ival_M_to_H_start)
        self.assertEqual(ival_M.asfreq('H', 'E'), ival_M_to_H_end)
        self.assertEqual(ival_M.asfreq('Min', 'S'), ival_M_to_T_start)
        self.assertEqual(ival_M.asfreq('Min', 'E'), ival_M_to_T_end)
        self.assertEqual(ival_M.asfreq('S', 'S'), ival_M_to_S_start)
        self.assertEqual(ival_M.asfreq('S', 'E'), ival_M_to_S_end)

        self.assertEqual(ival_M.asfreq('M'), ival_M)

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
            ival_W_to_Q_end_of_quarter = Period(freq='Q', year=2007, quarter=1)
        else:
            ival_W_to_Q_end_of_quarter = Period(freq='Q', year=2007, quarter=2)

        if Period(freq='D', year=2007, month=1, day=31).weekday == 6:
            ival_W_to_M_end_of_month = Period(freq='M', year=2007, month=1)
        else:
            ival_W_to_M_end_of_month = Period(freq='M', year=2007, month=2)

        ival_W_to_B_start = Period(freq='B', year=2007, month=1, day=1)
        ival_W_to_B_end = Period(freq='B', year=2007, month=1, day=5)
        ival_W_to_D_start = Period(freq='D', year=2007, month=1, day=1)
        ival_W_to_D_end = Period(freq='D', year=2007, month=1, day=7)
        ival_W_to_H_start = Period(freq='H', year=2007, month=1, day=1, hour=0)
        ival_W_to_H_end = Period(freq='H', year=2007, month=1, day=7, hour=23)
        ival_W_to_T_start = Period(freq='Min', year=2007, month=1, day=1,
                                   hour=0, minute=0)
        ival_W_to_T_end = Period(freq='Min', year=2007, month=1, day=7,
                                 hour=23, minute=59)
        ival_W_to_S_start = Period(freq='S', year=2007, month=1, day=1, hour=0,
                                   minute=0, second=0)
        ival_W_to_S_end = Period(freq='S', year=2007, month=1, day=7, hour=23,
                                 minute=59, second=59)

        self.assertEqual(ival_W.asfreq('A'), ival_W_to_A)
        self.assertEqual(ival_W_end_of_year.asfreq('A'),
                         ival_W_to_A_end_of_year)
        self.assertEqual(ival_W.asfreq('Q'), ival_W_to_Q)
        self.assertEqual(ival_W_end_of_quarter.asfreq('Q'),
                         ival_W_to_Q_end_of_quarter)
        self.assertEqual(ival_W.asfreq('M'), ival_W_to_M)
        self.assertEqual(ival_W_end_of_month.asfreq('M'),
                         ival_W_to_M_end_of_month)

        self.assertEqual(ival_W.asfreq('B', 'S'), ival_W_to_B_start)
        self.assertEqual(ival_W.asfreq('B', 'E'), ival_W_to_B_end)

        self.assertEqual(ival_W.asfreq('D', 'S'), ival_W_to_D_start)
        self.assertEqual(ival_W.asfreq('D', 'E'), ival_W_to_D_end)

        self.assertEqual(ival_WSUN.asfreq('D', 'S'), ival_WSUN_to_D_start)
        self.assertEqual(ival_WSUN.asfreq('D', 'E'), ival_WSUN_to_D_end)
        self.assertEqual(ival_WSAT.asfreq('D', 'S'), ival_WSAT_to_D_start)
        self.assertEqual(ival_WSAT.asfreq('D', 'E'), ival_WSAT_to_D_end)
        self.assertEqual(ival_WFRI.asfreq('D', 'S'), ival_WFRI_to_D_start)
        self.assertEqual(ival_WFRI.asfreq('D', 'E'), ival_WFRI_to_D_end)
        self.assertEqual(ival_WTHU.asfreq('D', 'S'), ival_WTHU_to_D_start)
        self.assertEqual(ival_WTHU.asfreq('D', 'E'), ival_WTHU_to_D_end)
        self.assertEqual(ival_WWED.asfreq('D', 'S'), ival_WWED_to_D_start)
        self.assertEqual(ival_WWED.asfreq('D', 'E'), ival_WWED_to_D_end)
        self.assertEqual(ival_WTUE.asfreq('D', 'S'), ival_WTUE_to_D_start)
        self.assertEqual(ival_WTUE.asfreq('D', 'E'), ival_WTUE_to_D_end)
        self.assertEqual(ival_WMON.asfreq('D', 'S'), ival_WMON_to_D_start)
        self.assertEqual(ival_WMON.asfreq('D', 'E'), ival_WMON_to_D_end)

        self.assertEqual(ival_W.asfreq('H', 'S'), ival_W_to_H_start)
        self.assertEqual(ival_W.asfreq('H', 'E'), ival_W_to_H_end)
        self.assertEqual(ival_W.asfreq('Min', 'S'), ival_W_to_T_start)
        self.assertEqual(ival_W.asfreq('Min', 'E'), ival_W_to_T_end)
        self.assertEqual(ival_W.asfreq('S', 'S'), ival_W_to_S_start)
        self.assertEqual(ival_W.asfreq('S', 'E'), ival_W_to_S_end)

        self.assertEqual(ival_W.asfreq('W'), ival_W)

        msg = pd.tseries.frequencies._INVALID_FREQ_ERROR
        with self.assertRaisesRegexp(ValueError, msg):
            ival_W.asfreq('WK')

    def test_conv_weekly_legacy(self):
        # frequency conversion tests: from Weekly Frequency
        msg = pd.tseries.frequencies._INVALID_FREQ_ERROR
        with self.assertRaisesRegexp(ValueError, msg):
            Period(freq='WK', year=2007, month=1, day=1)

        with self.assertRaisesRegexp(ValueError, msg):
            Period(freq='WK-SAT', year=2007, month=1, day=6)
        with self.assertRaisesRegexp(ValueError, msg):
            Period(freq='WK-FRI', year=2007, month=1, day=5)
        with self.assertRaisesRegexp(ValueError, msg):
            Period(freq='WK-THU', year=2007, month=1, day=4)
        with self.assertRaisesRegexp(ValueError, msg):
            Period(freq='WK-WED', year=2007, month=1, day=3)
        with self.assertRaisesRegexp(ValueError, msg):
            Period(freq='WK-TUE', year=2007, month=1, day=2)
        with self.assertRaisesRegexp(ValueError, msg):
            Period(freq='WK-MON', year=2007, month=1, day=1)

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
        ival_B_to_H_start = Period(freq='H', year=2007, month=1, day=1, hour=0)
        ival_B_to_H_end = Period(freq='H', year=2007, month=1, day=1, hour=23)
        ival_B_to_T_start = Period(freq='Min', year=2007, month=1, day=1,
                                   hour=0, minute=0)
        ival_B_to_T_end = Period(freq='Min', year=2007, month=1, day=1,
                                 hour=23, minute=59)
        ival_B_to_S_start = Period(freq='S', year=2007, month=1, day=1, hour=0,
                                   minute=0, second=0)
        ival_B_to_S_end = Period(freq='S', year=2007, month=1, day=1, hour=23,
                                 minute=59, second=59)

        self.assertEqual(ival_B.asfreq('A'), ival_B_to_A)
        self.assertEqual(ival_B_end_of_year.asfreq('A'), ival_B_to_A)
        self.assertEqual(ival_B.asfreq('Q'), ival_B_to_Q)
        self.assertEqual(ival_B_end_of_quarter.asfreq('Q'), ival_B_to_Q)
        self.assertEqual(ival_B.asfreq('M'), ival_B_to_M)
        self.assertEqual(ival_B_end_of_month.asfreq('M'), ival_B_to_M)
        self.assertEqual(ival_B.asfreq('W'), ival_B_to_W)
        self.assertEqual(ival_B_end_of_week.asfreq('W'), ival_B_to_W)

        self.assertEqual(ival_B.asfreq('D'), ival_B_to_D)

        self.assertEqual(ival_B.asfreq('H', 'S'), ival_B_to_H_start)
        self.assertEqual(ival_B.asfreq('H', 'E'), ival_B_to_H_end)
        self.assertEqual(ival_B.asfreq('Min', 'S'), ival_B_to_T_start)
        self.assertEqual(ival_B.asfreq('Min', 'E'), ival_B_to_T_end)
        self.assertEqual(ival_B.asfreq('S', 'S'), ival_B_to_S_start)
        self.assertEqual(ival_B.asfreq('S', 'E'), ival_B_to_S_end)

        self.assertEqual(ival_B.asfreq('B'), ival_B)

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

        # TODO: unused?
        # ival_D_monday = Period(freq='D', year=2007, month=1, day=8)

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

        ival_D_to_H_start = Period(freq='H', year=2007, month=1, day=1, hour=0)
        ival_D_to_H_end = Period(freq='H', year=2007, month=1, day=1, hour=23)
        ival_D_to_T_start = Period(freq='Min', year=2007, month=1, day=1,
                                   hour=0, minute=0)
        ival_D_to_T_end = Period(freq='Min', year=2007, month=1, day=1,
                                 hour=23, minute=59)
        ival_D_to_S_start = Period(freq='S', year=2007, month=1, day=1, hour=0,
                                   minute=0, second=0)
        ival_D_to_S_end = Period(freq='S', year=2007, month=1, day=1, hour=23,
                                 minute=59, second=59)

        self.assertEqual(ival_D.asfreq('A'), ival_D_to_A)

        self.assertEqual(ival_D_end_of_quarter.asfreq('A-JAN'),
                         ival_Deoq_to_AJAN)
        self.assertEqual(ival_D_end_of_quarter.asfreq('A-JUN'),
                         ival_Deoq_to_AJUN)
        self.assertEqual(ival_D_end_of_quarter.asfreq('A-DEC'),
                         ival_Deoq_to_ADEC)

        self.assertEqual(ival_D_end_of_year.asfreq('A'), ival_D_to_A)
        self.assertEqual(ival_D_end_of_quarter.asfreq('Q'), ival_D_to_QEDEC)
        self.assertEqual(ival_D.asfreq("Q-JAN"), ival_D_to_QEJAN)
        self.assertEqual(ival_D.asfreq("Q-JUN"), ival_D_to_QEJUN)
        self.assertEqual(ival_D.asfreq("Q-DEC"), ival_D_to_QEDEC)
        self.assertEqual(ival_D.asfreq('M'), ival_D_to_M)
        self.assertEqual(ival_D_end_of_month.asfreq('M'), ival_D_to_M)
        self.assertEqual(ival_D.asfreq('W'), ival_D_to_W)
        self.assertEqual(ival_D_end_of_week.asfreq('W'), ival_D_to_W)

        self.assertEqual(ival_D_friday.asfreq('B'), ival_B_friday)
        self.assertEqual(ival_D_saturday.asfreq('B', 'S'), ival_B_friday)
        self.assertEqual(ival_D_saturday.asfreq('B', 'E'), ival_B_monday)
        self.assertEqual(ival_D_sunday.asfreq('B', 'S'), ival_B_friday)
        self.assertEqual(ival_D_sunday.asfreq('B', 'E'), ival_B_monday)

        self.assertEqual(ival_D.asfreq('H', 'S'), ival_D_to_H_start)
        self.assertEqual(ival_D.asfreq('H', 'E'), ival_D_to_H_end)
        self.assertEqual(ival_D.asfreq('Min', 'S'), ival_D_to_T_start)
        self.assertEqual(ival_D.asfreq('Min', 'E'), ival_D_to_T_end)
        self.assertEqual(ival_D.asfreq('S', 'S'), ival_D_to_S_start)
        self.assertEqual(ival_D.asfreq('S', 'E'), ival_D_to_S_end)

        self.assertEqual(ival_D.asfreq('D'), ival_D)

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
        ival_H_to_T_end = Period(freq='Min', year=2007, month=1, day=1, hour=0,
                                 minute=59)
        ival_H_to_S_start = Period(freq='S', year=2007, month=1, day=1, hour=0,
                                   minute=0, second=0)
        ival_H_to_S_end = Period(freq='S', year=2007, month=1, day=1, hour=0,
                                 minute=59, second=59)

        self.assertEqual(ival_H.asfreq('A'), ival_H_to_A)
        self.assertEqual(ival_H_end_of_year.asfreq('A'), ival_H_to_A)
        self.assertEqual(ival_H.asfreq('Q'), ival_H_to_Q)
        self.assertEqual(ival_H_end_of_quarter.asfreq('Q'), ival_H_to_Q)
        self.assertEqual(ival_H.asfreq('M'), ival_H_to_M)
        self.assertEqual(ival_H_end_of_month.asfreq('M'), ival_H_to_M)
        self.assertEqual(ival_H.asfreq('W'), ival_H_to_W)
        self.assertEqual(ival_H_end_of_week.asfreq('W'), ival_H_to_W)
        self.assertEqual(ival_H.asfreq('D'), ival_H_to_D)
        self.assertEqual(ival_H_end_of_day.asfreq('D'), ival_H_to_D)
        self.assertEqual(ival_H.asfreq('B'), ival_H_to_B)
        self.assertEqual(ival_H_end_of_bus.asfreq('B'), ival_H_to_B)

        self.assertEqual(ival_H.asfreq('Min', 'S'), ival_H_to_T_start)
        self.assertEqual(ival_H.asfreq('Min', 'E'), ival_H_to_T_end)
        self.assertEqual(ival_H.asfreq('S', 'S'), ival_H_to_S_start)
        self.assertEqual(ival_H.asfreq('S', 'E'), ival_H_to_S_end)

        self.assertEqual(ival_H.asfreq('H'), ival_H)

    def test_conv_minutely(self):
        # frequency conversion tests: from Minutely Frequency"

        ival_T = Period(freq='Min', year=2007, month=1, day=1, hour=0,
                        minute=0)
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

        ival_T_to_S_start = Period(freq='S', year=2007, month=1, day=1, hour=0,
                                   minute=0, second=0)
        ival_T_to_S_end = Period(freq='S', year=2007, month=1, day=1, hour=0,
                                 minute=0, second=59)

        self.assertEqual(ival_T.asfreq('A'), ival_T_to_A)
        self.assertEqual(ival_T_end_of_year.asfreq('A'), ival_T_to_A)
        self.assertEqual(ival_T.asfreq('Q'), ival_T_to_Q)
        self.assertEqual(ival_T_end_of_quarter.asfreq('Q'), ival_T_to_Q)
        self.assertEqual(ival_T.asfreq('M'), ival_T_to_M)
        self.assertEqual(ival_T_end_of_month.asfreq('M'), ival_T_to_M)
        self.assertEqual(ival_T.asfreq('W'), ival_T_to_W)
        self.assertEqual(ival_T_end_of_week.asfreq('W'), ival_T_to_W)
        self.assertEqual(ival_T.asfreq('D'), ival_T_to_D)
        self.assertEqual(ival_T_end_of_day.asfreq('D'), ival_T_to_D)
        self.assertEqual(ival_T.asfreq('B'), ival_T_to_B)
        self.assertEqual(ival_T_end_of_bus.asfreq('B'), ival_T_to_B)
        self.assertEqual(ival_T.asfreq('H'), ival_T_to_H)
        self.assertEqual(ival_T_end_of_hour.asfreq('H'), ival_T_to_H)

        self.assertEqual(ival_T.asfreq('S', 'S'), ival_T_to_S_start)
        self.assertEqual(ival_T.asfreq('S', 'E'), ival_T_to_S_end)

        self.assertEqual(ival_T.asfreq('Min'), ival_T)

    def test_conv_secondly(self):
        # frequency conversion tests: from Secondly Frequency"

        ival_S = Period(freq='S', year=2007, month=1, day=1, hour=0, minute=0,
                        second=0)
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
        ival_S_to_H = Period(freq='H', year=2007, month=1, day=1, hour=0)
        ival_S_to_T = Period(freq='Min', year=2007, month=1, day=1, hour=0,
                             minute=0)

        self.assertEqual(ival_S.asfreq('A'), ival_S_to_A)
        self.assertEqual(ival_S_end_of_year.asfreq('A'), ival_S_to_A)
        self.assertEqual(ival_S.asfreq('Q'), ival_S_to_Q)
        self.assertEqual(ival_S_end_of_quarter.asfreq('Q'), ival_S_to_Q)
        self.assertEqual(ival_S.asfreq('M'), ival_S_to_M)
        self.assertEqual(ival_S_end_of_month.asfreq('M'), ival_S_to_M)
        self.assertEqual(ival_S.asfreq('W'), ival_S_to_W)
        self.assertEqual(ival_S_end_of_week.asfreq('W'), ival_S_to_W)
        self.assertEqual(ival_S.asfreq('D'), ival_S_to_D)
        self.assertEqual(ival_S_end_of_day.asfreq('D'), ival_S_to_D)
        self.assertEqual(ival_S.asfreq('B'), ival_S_to_B)
        self.assertEqual(ival_S_end_of_bus.asfreq('B'), ival_S_to_B)
        self.assertEqual(ival_S.asfreq('H'), ival_S_to_H)
        self.assertEqual(ival_S_end_of_hour.asfreq('H'), ival_S_to_H)
        self.assertEqual(ival_S.asfreq('Min'), ival_S_to_T)
        self.assertEqual(ival_S_end_of_minute.asfreq('Min'), ival_S_to_T)

        self.assertEqual(ival_S.asfreq('S'), ival_S)

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

    def test_asfreq_combined(self):
        # normal freq to combined freq
        p = Period('2007', freq='H')

        # ordinal will not change
        expected = Period('2007', freq='25H')
        for freq, how in zip(['1D1H', '1H1D'], ['E', 'S']):
            result = p.asfreq(freq, how=how)
            self.assertEqual(result, expected)
            self.assertEqual(result.ordinal, expected.ordinal)
            self.assertEqual(result.freq, expected.freq)

        # combined freq to normal freq
        p1 = Period(freq='1D1H', year=2007)
        p2 = Period(freq='1H1D', year=2007)

        # ordinal will change because how=E is the default
        result1 = p1.asfreq('H')
        result2 = p2.asfreq('H')
        expected = Period('2007-01-02', freq='H')
        self.assertEqual(result1, expected)
        self.assertEqual(result1.ordinal, expected.ordinal)
        self.assertEqual(result1.freq, expected.freq)
        self.assertEqual(result2, expected)
        self.assertEqual(result2.ordinal, expected.ordinal)
        self.assertEqual(result2.freq, expected.freq)

        # ordinal will not change
        result1 = p1.asfreq('H', how='S')
        result2 = p2.asfreq('H', how='S')
        expected = Period('2007-01-01', freq='H')
        self.assertEqual(result1, expected)
        self.assertEqual(result1.ordinal, expected.ordinal)
        self.assertEqual(result1.freq, expected.freq)
        self.assertEqual(result2, expected)
        self.assertEqual(result2.ordinal, expected.ordinal)
        self.assertEqual(result2.freq, expected.freq)

    def test_is_leap_year(self):
        # GH 13727
        for freq in ['A', 'M', 'D', 'H']:
            p = Period('2000-01-01 00:00:00', freq=freq)
            self.assertTrue(p.is_leap_year)
            self.assertIsInstance(p.is_leap_year, bool)

            p = Period('1999-01-01 00:00:00', freq=freq)
            self.assertFalse(p.is_leap_year)

            p = Period('2004-01-01 00:00:00', freq=freq)
            self.assertTrue(p.is_leap_year)

            p = Period('2100-01-01 00:00:00', freq=freq)
            self.assertFalse(p.is_leap_year)


class TestPeriodIndex(tm.TestCase):
    def setUp(self):
        pass

    def test_hash_error(self):
        index = period_range('20010101', periods=10)
        with tm.assertRaisesRegexp(TypeError, "unhashable type: %r" %
                                   type(index).__name__):
            hash(index)

    def test_make_time_series(self):
        index = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        series = Series(1, index=index)
        tm.assertIsInstance(series, Series)

    def test_constructor_use_start_freq(self):
        # GH #1118
        p = Period('4/2/2012', freq='B')
        index = PeriodIndex(start=p, periods=10)
        expected = PeriodIndex(start='4/2/2012', periods=10, freq='B')
        tm.assert_index_equal(index, expected)

    def test_constructor_field_arrays(self):
        # GH #1264

        years = np.arange(1990, 2010).repeat(4)[2:-2]
        quarters = np.tile(np.arange(1, 5), 20)[2:-2]

        index = PeriodIndex(year=years, quarter=quarters, freq='Q-DEC')
        expected = period_range('1990Q3', '2009Q2', freq='Q-DEC')
        tm.assert_index_equal(index, expected)

        index2 = PeriodIndex(year=years, quarter=quarters, freq='2Q-DEC')
        tm.assert_numpy_array_equal(index.asi8, index2.asi8)

        index = PeriodIndex(year=years, quarter=quarters)
        tm.assert_index_equal(index, expected)

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
        tm.assert_index_equal(idx, exp)

    def test_constructor_U(self):
        # U was used as undefined period
        self.assertRaises(ValueError, period_range, '2007-1-1', periods=500,
                          freq='X')

    def test_constructor_nano(self):
        idx = period_range(start=Period(ordinal=1, freq='N'),
                           end=Period(ordinal=4, freq='N'), freq='N')
        exp = PeriodIndex([Period(ordinal=1, freq='N'),
                           Period(ordinal=2, freq='N'),
                           Period(ordinal=3, freq='N'),
                           Period(ordinal=4, freq='N')], freq='N')
        tm.assert_index_equal(idx, exp)

    def test_constructor_arrays_negative_year(self):
        years = np.arange(1960, 2000, dtype=np.int64).repeat(4)
        quarters = np.tile(np.array([1, 2, 3, 4], dtype=np.int64), 40)

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
        tm.assert_index_equal(result, exp)

    def test_constructor_fromarraylike(self):
        idx = period_range('2007-01', periods=20, freq='M')

        # values is an array of Period, thus can retrieve freq
        tm.assert_index_equal(PeriodIndex(idx.values), idx)
        tm.assert_index_equal(PeriodIndex(list(idx.values)), idx)

        self.assertRaises(ValueError, PeriodIndex, idx._values)
        self.assertRaises(ValueError, PeriodIndex, list(idx._values))
        self.assertRaises(ValueError, PeriodIndex,
                          data=Period('2007', freq='A'))

        result = PeriodIndex(iter(idx))
        tm.assert_index_equal(result, idx)

        result = PeriodIndex(idx)
        tm.assert_index_equal(result, idx)

        result = PeriodIndex(idx, freq='M')
        tm.assert_index_equal(result, idx)

        result = PeriodIndex(idx, freq=offsets.MonthEnd())
        tm.assert_index_equal(result, idx)
        self.assertTrue(result.freq, 'M')

        result = PeriodIndex(idx, freq='2M')
        tm.assert_index_equal(result, idx.asfreq('2M'))
        self.assertTrue(result.freq, '2M')

        result = PeriodIndex(idx, freq=offsets.MonthEnd(2))
        tm.assert_index_equal(result, idx.asfreq('2M'))
        self.assertTrue(result.freq, '2M')

        result = PeriodIndex(idx, freq='D')
        exp = idx.asfreq('D', 'e')
        tm.assert_index_equal(result, exp)

    def test_constructor_datetime64arr(self):
        vals = np.arange(100000, 100000 + 10000, 100, dtype=np.int64)
        vals = vals.view(np.dtype('M8[us]'))

        self.assertRaises(ValueError, PeriodIndex, vals, freq='D')

    def test_constructor_dtype(self):
        # passing a dtype with a tz should localize
        idx = PeriodIndex(['2013-01', '2013-03'], dtype='period[M]')
        exp = PeriodIndex(['2013-01', '2013-03'], freq='M')
        tm.assert_index_equal(idx, exp)
        self.assertEqual(idx.dtype, 'period[M]')

        idx = PeriodIndex(['2013-01-05', '2013-03-05'], dtype='period[3D]')
        exp = PeriodIndex(['2013-01-05', '2013-03-05'], freq='3D')
        tm.assert_index_equal(idx, exp)
        self.assertEqual(idx.dtype, 'period[3D]')

        # if we already have a freq and its not the same, then asfreq
        # (not changed)
        idx = PeriodIndex(['2013-01-01', '2013-01-02'], freq='D')

        res = PeriodIndex(idx, dtype='period[M]')
        exp = PeriodIndex(['2013-01', '2013-01'], freq='M')
        tm.assert_index_equal(res, exp)
        self.assertEqual(res.dtype, 'period[M]')

        res = PeriodIndex(idx, freq='M')
        tm.assert_index_equal(res, exp)
        self.assertEqual(res.dtype, 'period[M]')

        msg = 'specified freq and dtype are different'
        with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
            PeriodIndex(['2011-01'], freq='M', dtype='period[D]')

    def test_constructor_empty(self):
        idx = pd.PeriodIndex([], freq='M')
        tm.assertIsInstance(idx, PeriodIndex)
        self.assertEqual(len(idx), 0)
        self.assertEqual(idx.freq, 'M')

        with tm.assertRaisesRegexp(ValueError, 'freq not specified'):
            pd.PeriodIndex([])

    def test_constructor_pi_nat(self):
        idx = PeriodIndex([Period('2011-01', freq='M'), pd.NaT,
                           Period('2011-01', freq='M')])
        exp = PeriodIndex(['2011-01', 'NaT', '2011-01'], freq='M')
        tm.assert_index_equal(idx, exp)

        idx = PeriodIndex(np.array([Period('2011-01', freq='M'), pd.NaT,
                                    Period('2011-01', freq='M')]))
        tm.assert_index_equal(idx, exp)

        idx = PeriodIndex([pd.NaT, pd.NaT, Period('2011-01', freq='M'),
                           Period('2011-01', freq='M')])
        exp = PeriodIndex(['NaT', 'NaT', '2011-01', '2011-01'], freq='M')
        tm.assert_index_equal(idx, exp)

        idx = PeriodIndex(np.array([pd.NaT, pd.NaT,
                                    Period('2011-01', freq='M'),
                                    Period('2011-01', freq='M')]))
        tm.assert_index_equal(idx, exp)

        idx = PeriodIndex([pd.NaT, pd.NaT, '2011-01', '2011-01'], freq='M')
        tm.assert_index_equal(idx, exp)

        with tm.assertRaisesRegexp(ValueError, 'freq not specified'):
            PeriodIndex([pd.NaT, pd.NaT])

        with tm.assertRaisesRegexp(ValueError, 'freq not specified'):
            PeriodIndex(np.array([pd.NaT, pd.NaT]))

        with tm.assertRaisesRegexp(ValueError, 'freq not specified'):
            PeriodIndex(['NaT', 'NaT'])

        with tm.assertRaisesRegexp(ValueError, 'freq not specified'):
            PeriodIndex(np.array(['NaT', 'NaT']))

    def test_constructor_incompat_freq(self):
        msg = "Input has different freq=D from PeriodIndex\\(freq=M\\)"

        with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
            PeriodIndex([Period('2011-01', freq='M'), pd.NaT,
                         Period('2011-01', freq='D')])

        with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
            PeriodIndex(np.array([Period('2011-01', freq='M'), pd.NaT,
                                  Period('2011-01', freq='D')]))

        # first element is pd.NaT
        with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
            PeriodIndex([pd.NaT, Period('2011-01', freq='M'),
                         Period('2011-01', freq='D')])

        with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
            PeriodIndex(np.array([pd.NaT, Period('2011-01', freq='M'),
                                  Period('2011-01', freq='D')]))

    def test_constructor_mixed(self):
        idx = PeriodIndex(['2011-01', pd.NaT, Period('2011-01', freq='M')])
        exp = PeriodIndex(['2011-01', 'NaT', '2011-01'], freq='M')
        tm.assert_index_equal(idx, exp)

        idx = PeriodIndex(['NaT', pd.NaT, Period('2011-01', freq='M')])
        exp = PeriodIndex(['NaT', 'NaT', '2011-01'], freq='M')
        tm.assert_index_equal(idx, exp)

        idx = PeriodIndex([Period('2011-01-01', freq='D'), pd.NaT,
                           '2012-01-01'])
        exp = PeriodIndex(['2011-01-01', 'NaT', '2012-01-01'], freq='D')
        tm.assert_index_equal(idx, exp)

    def test_constructor_simple_new(self):
        idx = period_range('2007-01', name='p', periods=2, freq='M')
        result = idx._simple_new(idx, 'p', freq=idx.freq)
        tm.assert_index_equal(result, idx)

        result = idx._simple_new(idx.astype('i8'), 'p', freq=idx.freq)
        tm.assert_index_equal(result, idx)

        result = idx._simple_new([pd.Period('2007-01', freq='M'),
                                  pd.Period('2007-02', freq='M')],
                                 'p', freq=idx.freq)
        self.assert_index_equal(result, idx)

        result = idx._simple_new(np.array([pd.Period('2007-01', freq='M'),
                                           pd.Period('2007-02', freq='M')]),
                                 'p', freq=idx.freq)
        self.assert_index_equal(result, idx)

    def test_constructor_simple_new_empty(self):
        # GH13079
        idx = PeriodIndex([], freq='M', name='p')
        result = idx._simple_new(idx, name='p', freq='M')
        tm.assert_index_equal(result, idx)

    def test_constructor_simple_new_floats(self):
        # GH13079
        for floats in [[1.1], np.array([1.1])]:
            with self.assertRaises(TypeError):
                pd.PeriodIndex._simple_new(floats, freq='M')

    def test_shallow_copy_empty(self):

        # GH13067
        idx = PeriodIndex([], freq='M')
        result = idx._shallow_copy()
        expected = idx

        tm.assert_index_equal(result, expected)

    def test_constructor_nat(self):
        self.assertRaises(ValueError, period_range, start='NaT',
                          end='2011-01-01', freq='M')
        self.assertRaises(ValueError, period_range, start='2011-01-01',
                          end='NaT', freq='M')

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
            expected = PeriodIndex(['2014-01', '2014-03',
                                    '2014-05', '2014-07'], freq='2M')
            tm.assert_index_equal(pidx, expected)

            pidx = func(start='2014-01-02', end='2014-01-15', freq='3D')
            expected = PeriodIndex(['2014-01-02', '2014-01-05',
                                    '2014-01-08', '2014-01-11',
                                    '2014-01-14'], freq='3D')
            tm.assert_index_equal(pidx, expected)

            pidx = func(end='2014-01-01 17:00', freq='4H', periods=3)
            expected = PeriodIndex(['2014-01-01 09:00', '2014-01-01 13:00',
                                    '2014-01-01 17:00'], freq='4H')
            tm.assert_index_equal(pidx, expected)

        msg = ('Frequency must be positive, because it'
               ' represents span: -1M')
        with tm.assertRaisesRegexp(ValueError, msg):
            PeriodIndex(['2011-01'], freq='-1M')

        msg = ('Frequency must be positive, because it' ' represents span: 0M')
        with tm.assertRaisesRegexp(ValueError, msg):
            PeriodIndex(['2011-01'], freq='0M')

        msg = ('Frequency must be positive, because it' ' represents span: 0M')
        with tm.assertRaisesRegexp(ValueError, msg):
            period_range('2011-01', periods=3, freq='0M')

    def test_constructor_freq_mult_dti_compat(self):
        import itertools
        mults = [1, 2, 3, 4, 5]
        freqs = ['A', 'M', 'D', 'T', 'S']
        for mult, freq in itertools.product(mults, freqs):
            freqstr = str(mult) + freq
            pidx = PeriodIndex(start='2014-04-01', freq=freqstr, periods=10)
            expected = date_range(start='2014-04-01', freq=freqstr,
                                  periods=10).to_period(freqstr)
            tm.assert_index_equal(pidx, expected)

    def test_constructor_freq_combined(self):
        for freq in ['1D1H', '1H1D']:
            pidx = PeriodIndex(['2016-01-01', '2016-01-02'], freq=freq)
            expected = PeriodIndex(['2016-01-01 00:00', '2016-01-02 00:00'],
                                   freq='25H')
        for freq, func in zip(['1D1H', '1H1D'], [PeriodIndex, period_range]):
            pidx = func(start='2016-01-01', periods=2, freq=freq)
            expected = PeriodIndex(['2016-01-01 00:00', '2016-01-02 01:00'],
                                   freq='25H')
            tm.assert_index_equal(pidx, expected)

    def test_dtype_str(self):
        pi = pd.PeriodIndex([], freq='M')
        self.assertEqual(pi.dtype_str, 'period[M]')
        self.assertEqual(pi.dtype_str, str(pi.dtype))

        pi = pd.PeriodIndex([], freq='3M')
        self.assertEqual(pi.dtype_str, 'period[3M]')
        self.assertEqual(pi.dtype_str, str(pi.dtype))

    def test_view_asi8(self):
        idx = pd.PeriodIndex([], freq='M')

        exp = np.array([], dtype=np.int64)
        tm.assert_numpy_array_equal(idx.view('i8'), exp)
        tm.assert_numpy_array_equal(idx.asi8, exp)

        idx = pd.PeriodIndex(['2011-01', pd.NaT], freq='M')

        exp = np.array([492, -9223372036854775808], dtype=np.int64)
        tm.assert_numpy_array_equal(idx.view('i8'), exp)
        tm.assert_numpy_array_equal(idx.asi8, exp)

        exp = np.array([14975, -9223372036854775808], dtype=np.int64)
        idx = pd.PeriodIndex(['2011-01-01', pd.NaT], freq='D')
        tm.assert_numpy_array_equal(idx.view('i8'), exp)
        tm.assert_numpy_array_equal(idx.asi8, exp)

    def test_values(self):
        idx = pd.PeriodIndex([], freq='M')

        exp = np.array([], dtype=np.object)
        tm.assert_numpy_array_equal(idx.values, exp)
        tm.assert_numpy_array_equal(idx.get_values(), exp)
        exp = np.array([], dtype=np.int64)
        tm.assert_numpy_array_equal(idx._values, exp)

        idx = pd.PeriodIndex(['2011-01', pd.NaT], freq='M')

        exp = np.array([pd.Period('2011-01', freq='M'), pd.NaT], dtype=object)
        tm.assert_numpy_array_equal(idx.values, exp)
        tm.assert_numpy_array_equal(idx.get_values(), exp)
        exp = np.array([492, -9223372036854775808], dtype=np.int64)
        tm.assert_numpy_array_equal(idx._values, exp)

        idx = pd.PeriodIndex(['2011-01-01', pd.NaT], freq='D')

        exp = np.array([pd.Period('2011-01-01', freq='D'), pd.NaT],
                       dtype=object)
        tm.assert_numpy_array_equal(idx.values, exp)
        tm.assert_numpy_array_equal(idx.get_values(), exp)
        exp = np.array([14975, -9223372036854775808], dtype=np.int64)
        tm.assert_numpy_array_equal(idx._values, exp)

    def test_asobject_like(self):
        idx = pd.PeriodIndex([], freq='M')

        exp = np.array([], dtype=object)
        tm.assert_numpy_array_equal(idx.asobject.values, exp)
        tm.assert_numpy_array_equal(idx._mpl_repr(), exp)

        idx = pd.PeriodIndex(['2011-01', pd.NaT], freq='M')

        exp = np.array([pd.Period('2011-01', freq='M'), pd.NaT], dtype=object)
        tm.assert_numpy_array_equal(idx.asobject.values, exp)
        tm.assert_numpy_array_equal(idx._mpl_repr(), exp)

        exp = np.array([pd.Period('2011-01-01', freq='D'), pd.NaT],
                       dtype=object)
        idx = pd.PeriodIndex(['2011-01-01', pd.NaT], freq='D')
        tm.assert_numpy_array_equal(idx.asobject.values, exp)
        tm.assert_numpy_array_equal(idx._mpl_repr(), exp)

    def test_is_(self):
        create_index = lambda: PeriodIndex(freq='A', start='1/1/2001',
                                           end='12/1/2009')
        index = create_index()
        self.assertEqual(index.is_(index), True)
        self.assertEqual(index.is_(create_index()), False)
        self.assertEqual(index.is_(index.view()), True)
        self.assertEqual(
            index.is_(index.view().view().view().view().view()), True)
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

    def test_getitem_index(self):
        idx = period_range('2007-01', periods=10, freq='M', name='x')

        result = idx[[1, 3, 5]]
        exp = pd.PeriodIndex(['2007-02', '2007-04', '2007-06'],
                             freq='M', name='x')
        tm.assert_index_equal(result, exp)

        result = idx[[True, True, False, False, False,
                      True, True, False, False, False]]
        exp = pd.PeriodIndex(['2007-01', '2007-02', '2007-06', '2007-07'],
                             freq='M', name='x')
        tm.assert_index_equal(result, exp)

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
        tm.assert_series_equal(exp, result)

        ts = ts[10:].append(ts[10:])
        self.assertRaisesRegexp(KeyError,
                                "left slice bound for non-unique "
                                "label: '2008'",
                                ts.__getitem__, slice('2008', '2009'))

    def test_getitem_datetime(self):
        rng = period_range(start='2012-01-01', periods=10, freq='W-MON')
        ts = Series(lrange(len(rng)), index=rng)

        dt1 = datetime(2011, 10, 2)
        dt4 = datetime(2012, 4, 20)

        rs = ts[dt1:dt4]
        tm.assert_series_equal(rs, ts)

    def test_getitem_nat(self):
        idx = pd.PeriodIndex(['2011-01', 'NaT', '2011-02'], freq='M')
        self.assertEqual(idx[0], pd.Period('2011-01', freq='M'))
        self.assertIs(idx[1], tslib.NaT)

        s = pd.Series([0, 1, 2], index=idx)
        self.assertEqual(s[pd.NaT], 1)

        s = pd.Series(idx, index=idx)
        self.assertEqual(s[pd.Period('2011-01', freq='M')],
                         pd.Period('2011-01', freq='M'))
        self.assertIs(s[pd.NaT], tslib.NaT)

    def test_getitem_list_periods(self):
        # GH 7710
        rng = period_range(start='2012-01-01', periods=10, freq='D')
        ts = Series(lrange(len(rng)), index=rng)
        exp = ts.iloc[[1]]
        tm.assert_series_equal(ts[[Period('2012-01-02', freq='D')]], exp)

    def test_slice_with_negative_step(self):
        ts = Series(np.arange(20),
                    period_range('2014-01', periods=20, freq='M'))
        SLC = pd.IndexSlice

        def assert_slices_equivalent(l_slc, i_slc):
            tm.assert_series_equal(ts[l_slc], ts.iloc[i_slc])
            tm.assert_series_equal(ts.loc[l_slc], ts.iloc[i_slc])
            tm.assert_series_equal(ts.ix[l_slc], ts.iloc[i_slc])

        assert_slices_equivalent(SLC[Period('2014-10')::-1], SLC[9::-1])
        assert_slices_equivalent(SLC['2014-10'::-1], SLC[9::-1])

        assert_slices_equivalent(SLC[:Period('2014-10'):-1], SLC[:8:-1])
        assert_slices_equivalent(SLC[:'2014-10':-1], SLC[:8:-1])

        assert_slices_equivalent(SLC['2015-02':'2014-10':-1], SLC[13:8:-1])
        assert_slices_equivalent(SLC[Period('2015-02'):Period('2014-10'):-1],
                                 SLC[13:8:-1])
        assert_slices_equivalent(SLC['2015-02':Period('2014-10'):-1],
                                 SLC[13:8:-1])
        assert_slices_equivalent(SLC[Period('2015-02'):'2014-10':-1],
                                 SLC[13:8:-1])

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

    def test_contains_nat(self):
        # GH13582
        idx = period_range('2007-01', freq='M', periods=10)
        self.assertFalse(pd.NaT in idx)
        self.assertFalse(None in idx)
        self.assertFalse(float('nan') in idx)
        self.assertFalse(np.nan in idx)

        idx = pd.PeriodIndex(['2011-01', 'NaT', '2011-02'], freq='M')
        self.assertTrue(pd.NaT in idx)
        self.assertTrue(None in idx)
        self.assertTrue(float('nan') in idx)
        self.assertTrue(np.nan in idx)

    def test_sub(self):
        rng = period_range('2007-01', periods=50)

        result = rng - 5
        exp = rng + (-5)
        tm.assert_index_equal(result, exp)

    def test_periods_number_check(self):
        with tm.assertRaises(ValueError):
            period_range('2011-1-1', '2012-1-1', 'B')

    def test_tolist(self):
        index = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        rs = index.tolist()
        [tm.assertIsInstance(x, Period) for x in rs]

        recon = PeriodIndex(rs)
        tm.assert_index_equal(index, recon)

    def test_to_timestamp(self):
        index = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        series = Series(1, index=index, name='foo')

        exp_index = date_range('1/1/2001', end='12/31/2009', freq='A-DEC')
        result = series.to_timestamp(how='end')
        tm.assert_index_equal(result.index, exp_index)
        self.assertEqual(result.name, 'foo')

        exp_index = date_range('1/1/2001', end='1/1/2009', freq='AS-JAN')
        result = series.to_timestamp(how='start')
        tm.assert_index_equal(result.index, exp_index)

        def _get_with_delta(delta, freq='A-DEC'):
            return date_range(to_datetime('1/1/2001') + delta,
                              to_datetime('12/31/2009') + delta, freq=freq)

        delta = timedelta(hours=23)
        result = series.to_timestamp('H', 'end')
        exp_index = _get_with_delta(delta)
        tm.assert_index_equal(result.index, exp_index)

        delta = timedelta(hours=23, minutes=59)
        result = series.to_timestamp('T', 'end')
        exp_index = _get_with_delta(delta)
        tm.assert_index_equal(result.index, exp_index)

        result = series.to_timestamp('S', 'end')
        delta = timedelta(hours=23, minutes=59, seconds=59)
        exp_index = _get_with_delta(delta)
        tm.assert_index_equal(result.index, exp_index)

        index = PeriodIndex(freq='H', start='1/1/2001', end='1/2/2001')
        series = Series(1, index=index, name='foo')

        exp_index = date_range('1/1/2001 00:59:59', end='1/2/2001 00:59:59',
                               freq='H')
        result = series.to_timestamp(how='end')
        tm.assert_index_equal(result.index, exp_index)
        self.assertEqual(result.name, 'foo')

    def test_to_timestamp_quarterly_bug(self):
        years = np.arange(1960, 2000).repeat(4)
        quarters = np.tile(lrange(1, 5), 40)

        pindex = PeriodIndex(year=years, quarter=quarters)

        stamps = pindex.to_timestamp('D', 'end')
        expected = DatetimeIndex([x.to_timestamp('D', 'end') for x in pindex])
        tm.assert_index_equal(stamps, expected)

    def test_to_timestamp_preserve_name(self):
        index = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009',
                            name='foo')
        self.assertEqual(index.name, 'foo')

        conv = index.to_timestamp('D')
        self.assertEqual(conv.name, 'foo')

    def test_to_timestamp_repr_is_code(self):
        zs = [Timestamp('99-04-17 00:00:00', tz='UTC'),
              Timestamp('2001-04-17 00:00:00', tz='UTC'),
              Timestamp('2001-04-17 00:00:00', tz='America/Los_Angeles'),
              Timestamp('2001-04-17 00:00:00', tz=None)]
        for z in zs:
            self.assertEqual(eval(repr(z)), z)

    def test_to_timestamp_pi_nat(self):
        # GH 7228
        index = PeriodIndex(['NaT', '2011-01', '2011-02'], freq='M',
                            name='idx')

        result = index.to_timestamp('D')
        expected = DatetimeIndex([pd.NaT, datetime(2011, 1, 1),
                                  datetime(2011, 2, 1)], name='idx')
        tm.assert_index_equal(result, expected)
        self.assertEqual(result.name, 'idx')

        result2 = result.to_period(freq='M')
        tm.assert_index_equal(result2, index)
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
        expected = DatetimeIndex(
            ['2011-01-01', 'NaT', '2011-02-01'], name='idx')
        self.assert_index_equal(result, expected)
        result = idx.to_timestamp(how='E')
        expected = DatetimeIndex(
            ['2011-02-28', 'NaT', '2011-03-31'], name='idx')
        self.assert_index_equal(result, expected)

    def test_to_timestamp_pi_combined(self):
        idx = PeriodIndex(start='2011', periods=2, freq='1D1H', name='idx')
        result = idx.to_timestamp()
        expected = DatetimeIndex(
            ['2011-01-01 00:00', '2011-01-02 01:00'], name='idx')
        self.assert_index_equal(result, expected)
        result = idx.to_timestamp(how='E')
        expected = DatetimeIndex(
            ['2011-01-02 00:59:59', '2011-01-03 01:59:59'], name='idx')
        self.assert_index_equal(result, expected)
        result = idx.to_timestamp(how='E', freq='H')
        expected = DatetimeIndex(
            ['2011-01-02 00:00', '2011-01-03 01:00'], name='idx')
        self.assert_index_equal(result, expected)

    def test_to_timestamp_to_period_astype(self):
        idx = DatetimeIndex([pd.NaT, '2011-01-01', '2011-02-01'], name='idx')

        res = idx.astype('period[M]')
        exp = PeriodIndex(['NaT', '2011-01', '2011-02'], freq='M', name='idx')
        tm.assert_index_equal(res, exp)

        res = idx.astype('period[3M]')
        exp = PeriodIndex(['NaT', '2011-01', '2011-02'], freq='3M', name='idx')
        self.assert_index_equal(res, exp)

    def test_start_time(self):
        index = PeriodIndex(freq='M', start='2016-01-01', end='2016-05-31')
        expected_index = date_range('2016-01-01', end='2016-05-31', freq='MS')
        tm.assert_index_equal(index.start_time, expected_index)

    def test_end_time(self):
        index = PeriodIndex(freq='M', start='2016-01-01', end='2016-05-31')
        expected_index = date_range('2016-01-01', end='2016-05-31', freq='M')
        tm.assert_index_equal(index.end_time, expected_index)

    def test_as_frame_columns(self):
        rng = period_range('1/1/2000', periods=5)
        df = DataFrame(randn(10, 5), columns=rng)

        ts = df[rng[0]]
        tm.assert_series_equal(ts, df.ix[:, 0])

        # GH # 1211
        repr(df)

        ts = df['1/1/2000']
        tm.assert_series_equal(ts, df.ix[:, 0])

    def test_indexing(self):

        # GH 4390, iat incorrectly indexing
        index = period_range('1/1/2001', periods=10)
        s = Series(randn(10), index=index)
        expected = s[index[0]]
        result = s.iat[0]
        self.assertEqual(expected, result)

    def test_frame_setitem(self):
        rng = period_range('1/1/2000', periods=5, name='index')
        df = DataFrame(randn(5, 3), index=rng)

        df['Index'] = rng
        rs = Index(df['Index'])
        tm.assert_index_equal(rs, rng, check_names=False)
        self.assertEqual(rs.name, 'Index')
        self.assertEqual(rng.name, 'index')

        rs = df.reset_index().set_index('index')
        tm.assertIsInstance(rs.index, PeriodIndex)
        tm.assert_index_equal(rs.index, rng)

    def test_period_set_index_reindex(self):
        # GH 6631
        df = DataFrame(np.random.random(6))
        idx1 = period_range('2011/01/01', periods=6, freq='M')
        idx2 = period_range('2013', periods=6, freq='A')

        df = df.set_index(idx1)
        tm.assert_index_equal(df.index, idx1)
        df = df.set_index(idx2)
        tm.assert_index_equal(df.index, idx2)

    def test_frame_to_time_stamp(self):
        K = 5
        index = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        df = DataFrame(randn(len(index), K), index=index)
        df['mix'] = 'a'

        exp_index = date_range('1/1/2001', end='12/31/2009', freq='A-DEC')
        result = df.to_timestamp('D', 'end')
        tm.assert_index_equal(result.index, exp_index)
        tm.assert_numpy_array_equal(result.values, df.values)

        exp_index = date_range('1/1/2001', end='1/1/2009', freq='AS-JAN')
        result = df.to_timestamp('D', 'start')
        tm.assert_index_equal(result.index, exp_index)

        def _get_with_delta(delta, freq='A-DEC'):
            return date_range(to_datetime('1/1/2001') + delta,
                              to_datetime('12/31/2009') + delta, freq=freq)

        delta = timedelta(hours=23)
        result = df.to_timestamp('H', 'end')
        exp_index = _get_with_delta(delta)
        tm.assert_index_equal(result.index, exp_index)

        delta = timedelta(hours=23, minutes=59)
        result = df.to_timestamp('T', 'end')
        exp_index = _get_with_delta(delta)
        tm.assert_index_equal(result.index, exp_index)

        result = df.to_timestamp('S', 'end')
        delta = timedelta(hours=23, minutes=59, seconds=59)
        exp_index = _get_with_delta(delta)
        tm.assert_index_equal(result.index, exp_index)

        # columns
        df = df.T

        exp_index = date_range('1/1/2001', end='12/31/2009', freq='A-DEC')
        result = df.to_timestamp('D', 'end', axis=1)
        tm.assert_index_equal(result.columns, exp_index)
        tm.assert_numpy_array_equal(result.values, df.values)

        exp_index = date_range('1/1/2001', end='1/1/2009', freq='AS-JAN')
        result = df.to_timestamp('D', 'start', axis=1)
        tm.assert_index_equal(result.columns, exp_index)

        delta = timedelta(hours=23)
        result = df.to_timestamp('H', 'end', axis=1)
        exp_index = _get_with_delta(delta)
        tm.assert_index_equal(result.columns, exp_index)

        delta = timedelta(hours=23, minutes=59)
        result = df.to_timestamp('T', 'end', axis=1)
        exp_index = _get_with_delta(delta)
        tm.assert_index_equal(result.columns, exp_index)

        result = df.to_timestamp('S', 'end', axis=1)
        delta = timedelta(hours=23, minutes=59, seconds=59)
        exp_index = _get_with_delta(delta)
        tm.assert_index_equal(result.columns, exp_index)

        # invalid axis
        tm.assertRaisesRegexp(ValueError, 'axis', df.to_timestamp, axis=2)

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
        tm.assert_series_equal(result, expected)
        result[:] = 1
        self.assertTrue((ts[1:3] == 1).all())

        # not monotonic
        idx = PeriodIndex([2000, 2007, 2007, 2009, 2007], freq='A-JUN')
        ts = Series(np.random.randn(len(idx)), index=idx)

        result = ts[2007]
        expected = ts[idx == 2007]
        tm.assert_series_equal(result, expected)

    def test_index_unique(self):
        idx = PeriodIndex([2000, 2007, 2007, 2009, 2009], freq='A-JUN')
        expected = PeriodIndex([2000, 2007, 2009], freq='A-JUN')
        self.assert_index_equal(idx.unique(), expected)
        self.assertEqual(idx.nunique(), 3)

        idx = PeriodIndex([2000, 2007, 2007, 2009, 2007], freq='A-JUN',
                          tz='US/Eastern')
        expected = PeriodIndex([2000, 2007, 2009], freq='A-JUN',
                               tz='US/Eastern')
        self.assert_index_equal(idx.unique(), expected)
        self.assertEqual(idx.nunique(), 3)

    def test_constructor(self):
        pi = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        self.assertEqual(len(pi), 9)

        pi = PeriodIndex(freq='Q', start='1/1/2001', end='12/1/2009')
        self.assertEqual(len(pi), 4 * 9)

        pi = PeriodIndex(freq='M', start='1/1/2001', end='12/1/2009')
        self.assertEqual(len(pi), 12 * 9)

        pi = PeriodIndex(freq='D', start='1/1/2001', end='12/31/2009')
        self.assertEqual(len(pi), 365 * 9 + 2)

        pi = PeriodIndex(freq='B', start='1/1/2001', end='12/31/2009')
        self.assertEqual(len(pi), 261 * 9)

        pi = PeriodIndex(freq='H', start='1/1/2001', end='12/31/2001 23:00')
        self.assertEqual(len(pi), 365 * 24)

        pi = PeriodIndex(freq='Min', start='1/1/2001', end='1/1/2001 23:59')
        self.assertEqual(len(pi), 24 * 60)

        pi = PeriodIndex(freq='S', start='1/1/2001', end='1/1/2001 23:59:59')
        self.assertEqual(len(pi), 24 * 60 * 60)

        start = Period('02-Apr-2005', 'B')
        i1 = PeriodIndex(start=start, periods=20)
        self.assertEqual(len(i1), 20)
        self.assertEqual(i1.freq, start.freq)
        self.assertEqual(i1[0], start)

        end_intv = Period('2006-12-31', 'W')
        i1 = PeriodIndex(end=end_intv, periods=10)
        self.assertEqual(len(i1), 10)
        self.assertEqual(i1.freq, end_intv.freq)
        self.assertEqual(i1[-1], end_intv)

        end_intv = Period('2006-12-31', '1w')
        i2 = PeriodIndex(end=end_intv, periods=10)
        self.assertEqual(len(i1), len(i2))
        self.assertTrue((i1 == i2).all())
        self.assertEqual(i1.freq, i2.freq)

        end_intv = Period('2006-12-31', ('w', 1))
        i2 = PeriodIndex(end=end_intv, periods=10)
        self.assertEqual(len(i1), len(i2))
        self.assertTrue((i1 == i2).all())
        self.assertEqual(i1.freq, i2.freq)

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
        self.assertEqual(len(i2), 2)
        self.assertEqual(i2[0], end_intv)

        i2 = PeriodIndex(np.array([end_intv, Period('2005-05-05', 'B')]))
        self.assertEqual(len(i2), 2)
        self.assertEqual(i2[0], end_intv)

        # Mixed freq should fail
        vals = [end_intv, Period('2006-12-31', 'w')]
        self.assertRaises(ValueError, PeriodIndex, vals)
        vals = np.array(vals)
        self.assertRaises(ValueError, PeriodIndex, vals)

    def test_numpy_repeat(self):
        index = period_range('20010101', periods=2)
        expected = PeriodIndex([Period('2001-01-01'), Period('2001-01-01'),
                                Period('2001-01-02'), Period('2001-01-02')])

        tm.assert_index_equal(np.repeat(index, 2), expected)

        msg = "the 'axis' parameter is not supported"
        tm.assertRaisesRegexp(ValueError, msg, np.repeat, index, 2, axis=1)

    def test_shift(self):
        pi1 = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='A', start='1/1/2002', end='12/1/2010')

        tm.assert_index_equal(pi1.shift(0), pi1)

        self.assertEqual(len(pi1), len(pi2))
        self.assert_index_equal(pi1.shift(1), pi2)

        pi1 = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='A', start='1/1/2000', end='12/1/2008')
        self.assertEqual(len(pi1), len(pi2))
        self.assert_index_equal(pi1.shift(-1), pi2)

        pi1 = PeriodIndex(freq='M', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='M', start='2/1/2001', end='1/1/2010')
        self.assertEqual(len(pi1), len(pi2))
        self.assert_index_equal(pi1.shift(1), pi2)

        pi1 = PeriodIndex(freq='M', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='M', start='12/1/2000', end='11/1/2009')
        self.assertEqual(len(pi1), len(pi2))
        self.assert_index_equal(pi1.shift(-1), pi2)

        pi1 = PeriodIndex(freq='D', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='D', start='1/2/2001', end='12/2/2009')
        self.assertEqual(len(pi1), len(pi2))
        self.assert_index_equal(pi1.shift(1), pi2)

        pi1 = PeriodIndex(freq='D', start='1/1/2001', end='12/1/2009')
        pi2 = PeriodIndex(freq='D', start='12/31/2000', end='11/30/2009')
        self.assertEqual(len(pi1), len(pi2))
        self.assert_index_equal(pi1.shift(-1), pi2)

    def test_shift_nat(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT',
                           '2011-04'], freq='M', name='idx')
        result = idx.shift(1)
        expected = PeriodIndex(['2011-02', '2011-03', 'NaT',
                                '2011-05'], freq='M', name='idx')
        tm.assert_index_equal(result, expected)
        self.assertEqual(result.name, expected.name)

    def test_shift_ndarray(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT',
                           '2011-04'], freq='M', name='idx')
        result = idx.shift(np.array([1, 2, 3, 4]))
        expected = PeriodIndex(['2011-02', '2011-04', 'NaT',
                                '2011-08'], freq='M', name='idx')
        tm.assert_index_equal(result, expected)

        idx = PeriodIndex(['2011-01', '2011-02', 'NaT',
                           '2011-04'], freq='M', name='idx')
        result = idx.shift(np.array([1, -2, 3, -4]))
        expected = PeriodIndex(['2011-02', '2010-12', 'NaT',
                                '2010-12'], freq='M', name='idx')
        tm.assert_index_equal(result, expected)

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
        tm.assert_index_equal(result, expected)

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

    def test_asfreq_combined_pi(self):
        pi = pd.PeriodIndex(['2001-01-01 00:00', '2001-01-02 02:00', 'NaT'],
                            freq='H')
        exp = PeriodIndex(['2001-01-01 00:00', '2001-01-02 02:00', 'NaT'],
                          freq='25H')
        for freq, how in zip(['1D1H', '1H1D'], ['S', 'E']):
            result = pi.asfreq(freq, how=how)
            self.assert_index_equal(result, exp)
            self.assertEqual(result.freq, exp.freq)

        for freq in ['1D1H', '1H1D']:
            pi = pd.PeriodIndex(['2001-01-01 00:00', '2001-01-02 02:00',
                                 'NaT'], freq=freq)
            result = pi.asfreq('H')
            exp = PeriodIndex(['2001-01-02 00:00', '2001-01-03 02:00', 'NaT'],
                              freq='H')
            self.assert_index_equal(result, exp)
            self.assertEqual(result.freq, exp.freq)

            pi = pd.PeriodIndex(['2001-01-01 00:00', '2001-01-02 02:00',
                                 'NaT'], freq=freq)
            result = pi.asfreq('H', how='S')
            exp = PeriodIndex(['2001-01-01 00:00', '2001-01-02 02:00', 'NaT'],
                              freq='H')
            self.assert_index_equal(result, exp)
            self.assertEqual(result.freq, exp.freq)

    def test_period_index_length(self):
        pi = PeriodIndex(freq='A', start='1/1/2001', end='12/1/2009')
        self.assertEqual(len(pi), 9)

        pi = PeriodIndex(freq='Q', start='1/1/2001', end='12/1/2009')
        self.assertEqual(len(pi), 4 * 9)

        pi = PeriodIndex(freq='M', start='1/1/2001', end='12/1/2009')
        self.assertEqual(len(pi), 12 * 9)

        start = Period('02-Apr-2005', 'B')
        i1 = PeriodIndex(start=start, periods=20)
        self.assertEqual(len(i1), 20)
        self.assertEqual(i1.freq, start.freq)
        self.assertEqual(i1[0], start)

        end_intv = Period('2006-12-31', 'W')
        i1 = PeriodIndex(end=end_intv, periods=10)
        self.assertEqual(len(i1), 10)
        self.assertEqual(i1.freq, end_intv.freq)
        self.assertEqual(i1[-1], end_intv)

        end_intv = Period('2006-12-31', '1w')
        i2 = PeriodIndex(end=end_intv, periods=10)
        self.assertEqual(len(i1), len(i2))
        self.assertTrue((i1 == i2).all())
        self.assertEqual(i1.freq, i2.freq)

        end_intv = Period('2006-12-31', ('w', 1))
        i2 = PeriodIndex(end=end_intv, periods=10)
        self.assertEqual(len(i1), len(i2))
        self.assertTrue((i1 == i2).all())
        self.assertEqual(i1.freq, i2.freq)

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
        self.assertEqual(len(i2), 2)
        self.assertEqual(i2[0], end_intv)

        i2 = PeriodIndex(np.array([end_intv, Period('2005-05-05', 'B')]))
        self.assertEqual(len(i2), 2)
        self.assertEqual(i2[0], end_intv)

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
        tm.assert_index_equal(result.index, exp_index)
        tm.assert_index_equal(df_result.index, exp_index)

        result = ts.asfreq('D', how='start')
        self.assertEqual(len(result), len(ts))
        tm.assert_index_equal(result.index, index.asfreq('D', how='start'))

    def test_badinput(self):
        self.assertRaises(ValueError, Period, '-2000', 'A')
        self.assertRaises(tslib.DateParseError, Period, '0', 'A')
        self.assertRaises(tslib.DateParseError, Period, '1/1/-2000', 'A')

    def test_negative_ordinals(self):
        Period(ordinal=-1000, freq='A')
        Period(ordinal=0, freq='A')

        idx1 = PeriodIndex(ordinal=[-1, 0, 1], freq='A')
        idx2 = PeriodIndex(ordinal=np.array([-1, 0, 1]), freq='A')
        tm.assert_index_equal(idx1, idx2)

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

        tm.assert_index_equal(pi1, period_range('1/1/2005', '11/1/2005',
                                                freq='M'))
        tm.assert_index_equal(pi2, period_range('1/1/2005', '11/1/2005',
                                                freq='M').asfreq('D'))
        tm.assert_index_equal(pi3, period_range('1/1/2005', '11/1/2005',
                                                freq='M').asfreq('3D'))

    def test_pindex_slice_index(self):
        pi = PeriodIndex(start='1/1/10', end='12/31/12', freq='M')
        s = Series(np.random.rand(len(pi)), index=pi)
        res = s['2010']
        exp = s[0:12]
        tm.assert_series_equal(res, exp)
        res = s['2011']
        exp = s[12:24]
        tm.assert_series_equal(res, exp)

    def test_getitem_day(self):
        # GH 6716
        # Confirm DatetimeIndex and PeriodIndex works identically
        didx = DatetimeIndex(start='2013/01/01', freq='D', periods=400)
        pidx = PeriodIndex(start='2013/01/01', freq='D', periods=400)

        for idx in [didx, pidx]:
            # getitem against index should raise ValueError
            values = ['2014', '2013/02', '2013/01/02', '2013/02/01 9H',
                      '2013/02/01 09:00']
            for v in values:

                if _np_version_under1p9:
                    with tm.assertRaises(ValueError):
                        idx[v]
                else:
                    # GH7116
                    # these show deprecations as we are trying
                    # to slice with non-integer indexers
                    # with tm.assertRaises(IndexError):
                    #    idx[v]
                    continue

            s = Series(np.random.rand(len(idx)), index=idx)
            tm.assert_series_equal(s['2013/01'], s[0:31])
            tm.assert_series_equal(s['2013/02'], s[31:59])
            tm.assert_series_equal(s['2014'], s[365:])

            invalid = ['2013/02/01 9H', '2013/02/01 09:00']
            for v in invalid:
                with tm.assertRaises(KeyError):
                    s[v]

    def test_range_slice_day(self):
        # GH 6716
        didx = DatetimeIndex(start='2013/01/01', freq='D', periods=400)
        pidx = PeriodIndex(start='2013/01/01', freq='D', periods=400)

        # changed to TypeError in 1.12
        # https://github.com/numpy/numpy/pull/6271
        exc = IndexError if _np_version_under1p12 else TypeError

        for idx in [didx, pidx]:
            # slices against index should raise IndexError
            values = ['2014', '2013/02', '2013/01/02', '2013/02/01 9H',
                      '2013/02/01 09:00']
            for v in values:
                with tm.assertRaises(exc):
                    idx[v:]

            s = Series(np.random.rand(len(idx)), index=idx)

            tm.assert_series_equal(s['2013/01/02':], s[1:])
            tm.assert_series_equal(s['2013/01/02':'2013/01/05'], s[1:5])
            tm.assert_series_equal(s['2013/02':], s[31:])
            tm.assert_series_equal(s['2014':], s[365:])

            invalid = ['2013/02/01 9H', '2013/02/01 09:00']
            for v in invalid:
                with tm.assertRaises(exc):
                    idx[v:]

    def test_getitem_seconds(self):
        # GH 6716
        didx = DatetimeIndex(start='2013/01/01 09:00:00', freq='S',
                             periods=4000)
        pidx = PeriodIndex(start='2013/01/01 09:00:00', freq='S', periods=4000)

        for idx in [didx, pidx]:
            # getitem against index should raise ValueError
            values = ['2014', '2013/02', '2013/01/02', '2013/02/01 9H',
                      '2013/02/01 09:00']
            for v in values:
                if _np_version_under1p9:
                    with tm.assertRaises(ValueError):
                        idx[v]
                else:
                    # GH7116
                    # these show deprecations as we are trying
                    # to slice with non-integer indexers
                    # with tm.assertRaises(IndexError):
                    #    idx[v]
                    continue

            s = Series(np.random.rand(len(idx)), index=idx)
            tm.assert_series_equal(s['2013/01/01 10:00'], s[3600:3660])
            tm.assert_series_equal(s['2013/01/01 9H'], s[:3600])
            for d in ['2013/01/01', '2013/01', '2013']:
                tm.assert_series_equal(s[d], s)

    def test_range_slice_seconds(self):
        # GH 6716
        didx = DatetimeIndex(start='2013/01/01 09:00:00', freq='S',
                             periods=4000)
        pidx = PeriodIndex(start='2013/01/01 09:00:00', freq='S', periods=4000)

        # changed to TypeError in 1.12
        # https://github.com/numpy/numpy/pull/6271
        exc = IndexError if _np_version_under1p12 else TypeError

        for idx in [didx, pidx]:
            # slices against index should raise IndexError
            values = ['2014', '2013/02', '2013/01/02', '2013/02/01 9H',
                      '2013/02/01 09:00']
            for v in values:
                with tm.assertRaises(exc):
                    idx[v:]

            s = Series(np.random.rand(len(idx)), index=idx)

            tm.assert_series_equal(s['2013/01/01 09:05':'2013/01/01 09:10'],
                                   s[300:660])
            tm.assert_series_equal(s['2013/01/01 10:00':'2013/01/01 10:05'],
                                   s[3600:3960])
            tm.assert_series_equal(s['2013/01/01 10H':], s[3600:])
            tm.assert_series_equal(s[:'2013/01/01 09:30'], s[:1860])
            for d in ['2013/01/01', '2013/01', '2013']:
                tm.assert_series_equal(s[d:], s)

    def test_range_slice_outofbounds(self):
        # GH 5407
        didx = DatetimeIndex(start='2013/10/01', freq='D', periods=10)
        pidx = PeriodIndex(start='2013/10/01', freq='D', periods=10)

        for idx in [didx, pidx]:
            df = DataFrame(dict(units=[100 + i for i in range(10)]), index=idx)
            empty = DataFrame(index=idx.__class__([], freq='D'),
                              columns=['units'])
            empty['units'] = empty['units'].astype('int64')

            tm.assert_frame_equal(df['2013/09/01':'2013/09/30'], empty)
            tm.assert_frame_equal(df['2013/09/30':'2013/10/02'], df.iloc[:2])
            tm.assert_frame_equal(df['2013/10/01':'2013/10/02'], df.iloc[:2])
            tm.assert_frame_equal(df['2013/10/02':'2013/09/30'], empty)
            tm.assert_frame_equal(df['2013/10/15':'2013/10/17'], empty)
            tm.assert_frame_equal(df['2013-06':'2013-09'], empty)
            tm.assert_frame_equal(df['2013-11':'2013-12'], empty)

    def test_astype_asfreq(self):
        pi1 = PeriodIndex(['2011-01-01', '2011-02-01', '2011-03-01'], freq='D')
        exp = PeriodIndex(['2011-01', '2011-02', '2011-03'], freq='M')
        tm.assert_index_equal(pi1.asfreq('M'), exp)
        tm.assert_index_equal(pi1.astype('period[M]'), exp)

        exp = PeriodIndex(['2011-01', '2011-02', '2011-03'], freq='3M')
        tm.assert_index_equal(pi1.asfreq('3M'), exp)
        tm.assert_index_equal(pi1.astype('period[3M]'), exp)

    def test_pindex_fieldaccessor_nat(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT',
                           '2012-03', '2012-04'], freq='D')

        exp = np.array([2011, 2011, -1, 2012, 2012], dtype=np.int64)
        self.assert_numpy_array_equal(idx.year, exp)
        exp = np.array([1, 2, -1, 3, 4], dtype=np.int64)
        self.assert_numpy_array_equal(idx.month, exp)

    def test_pindex_qaccess(self):
        pi = PeriodIndex(['2Q05', '3Q05', '4Q05', '1Q06', '2Q06'], freq='Q')
        s = Series(np.random.rand(len(pi)), index=pi).cumsum()
        # Todo: fix these accessors!
        self.assertEqual(s['05Q4'], s[2])

    def test_period_dt64_round_trip(self):
        dti = date_range('1/1/2000', '1/7/2002', freq='B')
        pi = dti.to_period()
        tm.assert_index_equal(pi.to_timestamp(), dti)

        dti = date_range('1/1/2000', '1/7/2002', freq='B')
        pi = dti.to_period(freq='H')
        tm.assert_index_equal(pi.to_timestamp(), dti)

    def test_period_astype_to_timestamp(self):
        pi = pd.PeriodIndex(['2011-01', '2011-02', '2011-03'], freq='M')

        exp = pd.DatetimeIndex(['2011-01-01', '2011-02-01', '2011-03-01'])
        tm.assert_index_equal(pi.astype('datetime64[ns]'), exp)

        exp = pd.DatetimeIndex(['2011-01-31', '2011-02-28', '2011-03-31'])
        tm.assert_index_equal(pi.astype('datetime64[ns]', how='end'), exp)

        exp = pd.DatetimeIndex(['2011-01-01', '2011-02-01', '2011-03-01'],
                               tz='US/Eastern')
        res = pi.astype('datetime64[ns, US/Eastern]')
        tm.assert_index_equal(pi.astype('datetime64[ns, US/Eastern]'), exp)

        exp = pd.DatetimeIndex(['2011-01-31', '2011-02-28', '2011-03-31'],
                               tz='US/Eastern')
        res = pi.astype('datetime64[ns, US/Eastern]', how='end')
        tm.assert_index_equal(res, exp)

    def test_to_period_quarterly(self):
        # make sure we can make the round trip
        for month in MONTHS:
            freq = 'Q-%s' % month
            rng = period_range('1989Q3', '1991Q3', freq=freq)
            stamps = rng.to_timestamp()
            result = stamps.to_period(freq)
            tm.assert_index_equal(rng, result)

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

        rng = date_range('01-Jan-2012', periods=8, freq='M')
        prng = rng.to_period()
        self.assertEqual(prng.freq, 'M')

        msg = pd.tseries.frequencies._INVALID_FREQ_ERROR
        with self.assertRaisesRegexp(ValueError, msg):
            date_range('01-Jan-2012', periods=8, freq='EOM')

    def test_multiples(self):
        result1 = Period('1989', freq='2A')
        result2 = Period('1989', freq='A')
        self.assertEqual(result1.ordinal, result2.ordinal)
        self.assertEqual(result1.freqstr, '2A-DEC')
        self.assertEqual(result2.freqstr, 'A-DEC')
        self.assertEqual(result1.freq, offsets.YearEnd(2))
        self.assertEqual(result2.freq, offsets.YearEnd())

        self.assertEqual((result1 + 1).ordinal, result1.ordinal + 2)
        self.assertEqual((1 + result1).ordinal, result1.ordinal + 2)
        self.assertEqual((result1 - 1).ordinal, result2.ordinal - 2)
        self.assertEqual((-1 + result1).ordinal, result2.ordinal - 2)

    def test_pindex_multiples(self):
        pi = PeriodIndex(start='1/1/11', end='12/31/11', freq='2M')
        expected = PeriodIndex(['2011-01', '2011-03', '2011-05', '2011-07',
                                '2011-09', '2011-11'], freq='2M')
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
        index = PeriodIndex(start='1/1/10', end='12/31/12', freq='D',
                            name='idx')
        expected = PeriodIndex([datetime(2010, 1, 6), datetime(2010, 1, 7),
                                datetime(2010, 1, 9), datetime(2010, 1, 13)],
                               freq='D', name='idx')

        taken1 = index.take([5, 6, 8, 12])
        taken2 = index[[5, 6, 8, 12]]

        for taken in [taken1, taken2]:
            tm.assert_index_equal(taken, expected)
            tm.assertIsInstance(taken, PeriodIndex)
            self.assertEqual(taken.freq, index.freq)
            self.assertEqual(taken.name, expected.name)

    def test_take_fill_value(self):
        # GH 12631
        idx = pd.PeriodIndex(['2011-01-01', '2011-02-01', '2011-03-01'],
                             name='xxx', freq='D')
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.PeriodIndex(['2011-02-01', '2011-01-01', '2011-03-01'],
                                  name='xxx', freq='D')
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.PeriodIndex(['2011-02-01', '2011-01-01', 'NaT'],
                                  name='xxx', freq='D')
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        expected = pd.PeriodIndex(['2011-02-01', '2011-01-01', '2011-03-01'],
                                  name='xxx', freq='D')
        tm.assert_index_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with tm.assertRaises(IndexError):
            idx.take(np.array([1, -5]))

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
        df = tm.makeCustomDataframe(
            3, 2, data_gen_f=lambda *args: np.random.randint(2),
            c_idx_type='p', r_idx_type='dt')
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
        tm.assert_series_equal(result, expected)

        result = ts + _permute(ts[::2])
        tm.assert_series_equal(result, expected)

        # it works!
        for kind in ['inner', 'outer', 'left', 'right']:
            ts.align(ts[::2], join=kind)
        msg = "Input has different freq=D from PeriodIndex\\(freq=A-DEC\\)"
        with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
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
        tm.assert_index_equal(result, index)

        # not in order
        result = _permute(index[:-5]).union(_permute(index[10:]))
        tm.assert_index_equal(result, index)

        # raise if different frequencies
        index = period_range('1/1/2000', '1/20/2000', freq='D')
        index2 = period_range('1/1/2000', '1/20/2000', freq='W-WED')
        with tm.assertRaises(period.IncompatibleFrequency):
            index.union(index2)

        msg = 'can only call with other PeriodIndex-ed objects'
        with tm.assertRaisesRegexp(ValueError, msg):
            index.join(index.to_timestamp())

        index3 = period_range('1/1/2000', '1/20/2000', freq='2D')
        with tm.assertRaises(period.IncompatibleFrequency):
            index.join(index3)

    def test_union_dataframe_index(self):
        rng1 = pd.period_range('1/1/1999', '1/1/2012', freq='M')
        s1 = pd.Series(np.random.randn(len(rng1)), rng1)

        rng2 = pd.period_range('1/1/1980', '12/1/2001', freq='M')
        s2 = pd.Series(np.random.randn(len(rng2)), rng2)
        df = pd.DataFrame({'s1': s1, 's2': s2})

        exp = pd.period_range('1/1/1980', '1/1/2012', freq='M')
        self.assert_index_equal(df.index, exp)

    def test_intersection(self):
        index = period_range('1/1/2000', '1/20/2000', freq='D')

        result = index[:-5].intersection(index[10:])
        tm.assert_index_equal(result, index[10:-5])

        # not in order
        left = _permute(index[:-5])
        right = _permute(index[10:])
        result = left.intersection(right).sort_values()
        tm.assert_index_equal(result, index[10:-5])

        # raise if different frequencies
        index = period_range('1/1/2000', '1/20/2000', freq='D')
        index2 = period_range('1/1/2000', '1/20/2000', freq='W-WED')
        with tm.assertRaises(period.IncompatibleFrequency):
            index.intersection(index2)

        index3 = period_range('1/1/2000', '1/20/2000', freq='2D')
        with tm.assertRaises(period.IncompatibleFrequency):
            index.intersection(index3)

    def test_intersection_cases(self):
        base = period_range('6/1/2000', '6/30/2000', freq='D', name='idx')

        # if target has the same name, it is preserved
        rng2 = period_range('5/15/2000', '6/20/2000', freq='D', name='idx')
        expected2 = period_range('6/1/2000', '6/20/2000', freq='D',
                                 name='idx')

        # if target name is different, it will be reset
        rng3 = period_range('5/15/2000', '6/20/2000', freq='D', name='other')
        expected3 = period_range('6/1/2000', '6/20/2000', freq='D',
                                 name=None)

        rng4 = period_range('7/1/2000', '7/31/2000', freq='D', name='idx')
        expected4 = PeriodIndex([], name='idx', freq='D')

        for (rng, expected) in [(rng2, expected2), (rng3, expected3),
                                (rng4, expected4)]:
            result = base.intersection(rng)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)

        # non-monotonic
        base = PeriodIndex(['2011-01-05', '2011-01-04', '2011-01-02',
                            '2011-01-03'], freq='D', name='idx')

        rng2 = PeriodIndex(['2011-01-04', '2011-01-02',
                            '2011-02-02', '2011-02-03'],
                           freq='D', name='idx')
        expected2 = PeriodIndex(['2011-01-04', '2011-01-02'], freq='D',
                                name='idx')

        rng3 = PeriodIndex(['2011-01-04', '2011-01-02', '2011-02-02',
                            '2011-02-03'],
                           freq='D', name='other')
        expected3 = PeriodIndex(['2011-01-04', '2011-01-02'], freq='D',
                                name=None)

        rng4 = period_range('7/1/2000', '7/31/2000', freq='D', name='idx')
        expected4 = PeriodIndex([], freq='D', name='idx')

        for (rng, expected) in [(rng2, expected2), (rng3, expected3),
                                (rng4, expected4)]:
            result = base.intersection(rng)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, 'D')

        # empty same freq
        rng = date_range('6/1/2000', '6/15/2000', freq='T')
        result = rng[0:0].intersection(rng)
        self.assertEqual(len(result), 0)

        result = rng.intersection(rng[0:0])
        self.assertEqual(len(result), 0)

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
        fields = ['year', 'month', 'day', 'hour', 'minute', 'second',
                  'weekofyear', 'week', 'dayofweek', 'weekday', 'dayofyear',
                  'quarter', 'qyear', 'days_in_month', 'is_leap_year']

        periods = list(periodindex)
        s = pd.Series(periodindex)

        for field in fields:
            field_idx = getattr(periodindex, field)
            self.assertEqual(len(periodindex), len(field_idx))
            for x, val in zip(periods, field_idx):
                self.assertEqual(getattr(x, field), val)

            if len(s) == 0:
                continue

            field_s = getattr(s.dt, field)
            self.assertEqual(len(periodindex), len(field_s))
            for x, val in zip(periods, field_s):
                self.assertEqual(getattr(x, field), val)

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
        tm.assert_index_equal(result, expected)

        result = index.map(lambda x: x.ordinal)
        exp = np.array([x.ordinal for x in index], dtype=np.int64)
        tm.assert_numpy_array_equal(result, exp)

    def test_map_with_string_constructor(self):
        raw = [2005, 2007, 2009]
        index = PeriodIndex(raw, freq='A')
        types = str,

        if PY3:
            # unicode
            types += text_type,

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

    def test_to_timestamp_1703(self):
        index = period_range('1/1/2012', periods=4, freq='D')

        result = index.to_timestamp()
        self.assertEqual(result[0], Timestamp('1/1/2012'))

    def test_to_datetime_depr(self):
        index = period_range('1/1/2012', periods=4, freq='D')

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
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

    def test_get_loc_nat(self):
        didx = DatetimeIndex(['2011-01-01', 'NaT', '2011-01-03'])
        pidx = PeriodIndex(['2011-01-01', 'NaT', '2011-01-03'], freq='M')

        # check DatetimeIndex compat
        for idx in [didx, pidx]:
            self.assertEqual(idx.get_loc(pd.NaT), 1)
            self.assertEqual(idx.get_loc(None), 1)
            self.assertEqual(idx.get_loc(float('nan')), 1)
            self.assertEqual(idx.get_loc(np.nan), 1)

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

        exp_arr = np.array([0, 0, 1, 1, 2, 2], dtype=np.intp)
        exp_idx = PeriodIndex(['2014-01', '2014-02', '2014-03'], freq='M')

        arr, idx = idx1.factorize()
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        arr, idx = idx1.factorize(sort=True)
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        idx2 = pd.PeriodIndex(['2014-03', '2014-03', '2014-02', '2014-01',
                               '2014-03', '2014-01'], freq='M')

        exp_arr = np.array([2, 2, 1, 0, 2, 0], dtype=np.intp)
        arr, idx = idx2.factorize(sort=True)
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        exp_arr = np.array([0, 0, 1, 2, 0, 2], dtype=np.intp)
        exp_idx = PeriodIndex(['2014-03', '2014-02', '2014-01'], freq='M')
        arr, idx = idx2.factorize()
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

    def test_recreate_from_data(self):
        for o in ['M', 'Q', 'A', 'D', 'B', 'T', 'S', 'L', 'U', 'N', 'H']:
            org = PeriodIndex(start='2001/04/01', freq=o, periods=1)
            idx = PeriodIndex(org.values, freq=o)
            tm.assert_index_equal(idx, org)

    def test_combine_first(self):
        # GH 3367
        didx = pd.DatetimeIndex(start='1950-01-31', end='1950-07-31', freq='M')
        pidx = pd.PeriodIndex(start=pd.Period('1950-1'),
                              end=pd.Period('1950-7'), freq='M')
        # check to be consistent with DatetimeIndex
        for idx in [didx, pidx]:
            a = pd.Series([1, np.nan, np.nan, 4, 5, np.nan, 7], index=idx)
            b = pd.Series([9, 9, 9, 9, 9, 9, 9], index=idx)

            result = a.combine_first(b)
            expected = pd.Series([1, 9, 9, 4, 5, 9, 7], index=idx,
                                 dtype=np.float64)
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
            with self.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                pidx.searchsorted(pd.Period('2014-01-01', freq='H'))

            msg = "Input has different freq=5D from PeriodIndex"
            with self.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                pidx.searchsorted(pd.Period('2014-01-01', freq='5D'))

    def test_round_trip(self):

        p = Period('2000Q1')
        new_p = self.round_trip_pickle(p)
        self.assertEqual(new_p, p)


def _permute(obj):
    return obj.take(np.random.permutation(len(obj)))


class TestMethods(tm.TestCase):

    def test_add(self):
        dt1 = Period(freq='D', year=2008, month=1, day=1)
        dt2 = Period(freq='D', year=2008, month=1, day=2)
        self.assertEqual(dt1 + 1, dt2)
        self.assertEqual(1 + dt1, dt2)

    def test_add_pdnat(self):
        p = pd.Period('2011-01', freq='M')
        self.assertIs(p + pd.NaT, pd.NaT)
        self.assertIs(pd.NaT + p, pd.NaT)

        p = pd.Period('NaT', freq='M')
        self.assertIs(p + pd.NaT, pd.NaT)
        self.assertIs(pd.NaT + p, pd.NaT)

    def test_add_raises(self):
        # GH 4731
        dt1 = Period(freq='D', year=2008, month=1, day=1)
        dt2 = Period(freq='D', year=2008, month=1, day=2)
        msg = r"unsupported operand type\(s\)"
        with tm.assertRaisesRegexp(TypeError, msg):
            dt1 + "str"

        msg = r"unsupported operand type\(s\)"
        with tm.assertRaisesRegexp(TypeError, msg):
            "str" + dt1

        with tm.assertRaisesRegexp(TypeError, msg):
            dt1 + dt2

    def test_sub(self):
        dt1 = Period('2011-01-01', freq='D')
        dt2 = Period('2011-01-15', freq='D')

        self.assertEqual(dt1 - dt2, -14)
        self.assertEqual(dt2 - dt1, 14)

        msg = r"Input has different freq=M from Period\(freq=D\)"
        with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
            dt1 - pd.Period('2011-02', freq='M')

    def test_add_offset(self):
        # freq is DateOffset
        for freq in ['A', '2A', '3A']:
            p = Period('2011', freq=freq)
            exp = Period('2013', freq=freq)
            self.assertEqual(p + offsets.YearEnd(2), exp)
            self.assertEqual(offsets.YearEnd(2) + p, exp)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(365, 'D'),
                      timedelta(365)]:
                with tm.assertRaises(period.IncompatibleFrequency):
                    p + o

                if isinstance(o, np.timedelta64):
                    with tm.assertRaises(TypeError):
                        o + p
                else:
                    with tm.assertRaises(period.IncompatibleFrequency):
                        o + p

        for freq in ['M', '2M', '3M']:
            p = Period('2011-03', freq=freq)
            exp = Period('2011-05', freq=freq)
            self.assertEqual(p + offsets.MonthEnd(2), exp)
            self.assertEqual(offsets.MonthEnd(2) + p, exp)

            exp = Period('2012-03', freq=freq)
            self.assertEqual(p + offsets.MonthEnd(12), exp)
            self.assertEqual(offsets.MonthEnd(12) + p, exp)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(365, 'D'),
                      timedelta(365)]:
                with tm.assertRaises(period.IncompatibleFrequency):
                    p + o

                if isinstance(o, np.timedelta64):
                    with tm.assertRaises(TypeError):
                        o + p
                else:
                    with tm.assertRaises(period.IncompatibleFrequency):
                        o + p

        # freq is Tick
        for freq in ['D', '2D', '3D']:
            p = Period('2011-04-01', freq=freq)

            exp = Period('2011-04-06', freq=freq)
            self.assertEqual(p + offsets.Day(5), exp)
            self.assertEqual(offsets.Day(5) + p, exp)

            exp = Period('2011-04-02', freq=freq)
            self.assertEqual(p + offsets.Hour(24), exp)
            self.assertEqual(offsets.Hour(24) + p, exp)

            exp = Period('2011-04-03', freq=freq)
            self.assertEqual(p + np.timedelta64(2, 'D'), exp)
            with tm.assertRaises(TypeError):
                np.timedelta64(2, 'D') + p

            exp = Period('2011-04-02', freq=freq)
            self.assertEqual(p + np.timedelta64(3600 * 24, 's'), exp)
            with tm.assertRaises(TypeError):
                np.timedelta64(3600 * 24, 's') + p

            exp = Period('2011-03-30', freq=freq)
            self.assertEqual(p + timedelta(-2), exp)
            self.assertEqual(timedelta(-2) + p, exp)

            exp = Period('2011-04-03', freq=freq)
            self.assertEqual(p + timedelta(hours=48), exp)
            self.assertEqual(timedelta(hours=48) + p, exp)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(4, 'h'),
                      timedelta(hours=23)]:
                with tm.assertRaises(period.IncompatibleFrequency):
                    p + o

                if isinstance(o, np.timedelta64):
                    with tm.assertRaises(TypeError):
                        o + p
                else:
                    with tm.assertRaises(period.IncompatibleFrequency):
                        o + p

        for freq in ['H', '2H', '3H']:
            p = Period('2011-04-01 09:00', freq=freq)

            exp = Period('2011-04-03 09:00', freq=freq)
            self.assertEqual(p + offsets.Day(2), exp)
            self.assertEqual(offsets.Day(2) + p, exp)

            exp = Period('2011-04-01 12:00', freq=freq)
            self.assertEqual(p + offsets.Hour(3), exp)
            self.assertEqual(offsets.Hour(3) + p, exp)

            exp = Period('2011-04-01 12:00', freq=freq)
            self.assertEqual(p + np.timedelta64(3, 'h'), exp)
            with tm.assertRaises(TypeError):
                np.timedelta64(3, 'h') + p

            exp = Period('2011-04-01 10:00', freq=freq)
            self.assertEqual(p + np.timedelta64(3600, 's'), exp)
            with tm.assertRaises(TypeError):
                np.timedelta64(3600, 's') + p

            exp = Period('2011-04-01 11:00', freq=freq)
            self.assertEqual(p + timedelta(minutes=120), exp)
            self.assertEqual(timedelta(minutes=120) + p, exp)

            exp = Period('2011-04-05 12:00', freq=freq)
            self.assertEqual(p + timedelta(days=4, minutes=180), exp)
            self.assertEqual(timedelta(days=4, minutes=180) + p, exp)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(3200, 's'),
                      timedelta(hours=23, minutes=30)]:
                with tm.assertRaises(period.IncompatibleFrequency):
                    p + o

                if isinstance(o, np.timedelta64):
                    with tm.assertRaises(TypeError):
                        o + p
                else:
                    with tm.assertRaises(period.IncompatibleFrequency):
                        o + p

    def test_add_offset_nat(self):
        # freq is DateOffset
        for freq in ['A', '2A', '3A']:
            p = Period('NaT', freq=freq)
            for o in [offsets.YearEnd(2)]:
                self.assertIs(p + o, tslib.NaT)
                self.assertIs(o + p, tslib.NaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(365, 'D'),
                      timedelta(365)]:
                self.assertIs(p + o, tslib.NaT)

                if isinstance(o, np.timedelta64):
                    with tm.assertRaises(TypeError):
                        o + p
                else:
                    self.assertIs(o + p, tslib.NaT)

        for freq in ['M', '2M', '3M']:
            p = Period('NaT', freq=freq)
            for o in [offsets.MonthEnd(2), offsets.MonthEnd(12)]:
                self.assertIs(p + o, tslib.NaT)

                if isinstance(o, np.timedelta64):
                    with tm.assertRaises(TypeError):
                        o + p
                else:
                    self.assertIs(o + p, tslib.NaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(365, 'D'),
                      timedelta(365)]:
                self.assertIs(p + o, tslib.NaT)

                if isinstance(o, np.timedelta64):
                    with tm.assertRaises(TypeError):
                        o + p
                else:
                    self.assertIs(o + p, tslib.NaT)

        # freq is Tick
        for freq in ['D', '2D', '3D']:
            p = Period('NaT', freq=freq)
            for o in [offsets.Day(5), offsets.Hour(24), np.timedelta64(2, 'D'),
                      np.timedelta64(3600 * 24, 's'), timedelta(-2),
                      timedelta(hours=48)]:
                self.assertIs(p + o, tslib.NaT)

                if isinstance(o, np.timedelta64):
                    with tm.assertRaises(TypeError):
                        o + p
                else:
                    self.assertIs(o + p, tslib.NaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(4, 'h'),
                      timedelta(hours=23)]:
                self.assertIs(p + o, tslib.NaT)

                if isinstance(o, np.timedelta64):
                    with tm.assertRaises(TypeError):
                        o + p
                else:
                    self.assertIs(o + p, tslib.NaT)

        for freq in ['H', '2H', '3H']:
            p = Period('NaT', freq=freq)
            for o in [offsets.Day(2), offsets.Hour(3), np.timedelta64(3, 'h'),
                      np.timedelta64(3600, 's'), timedelta(minutes=120),
                      timedelta(days=4, minutes=180)]:
                self.assertIs(p + o, tslib.NaT)

                if not isinstance(o, np.timedelta64):
                    self.assertIs(o + p, tslib.NaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(3200, 's'),
                      timedelta(hours=23, minutes=30)]:
                self.assertIs(p + o, tslib.NaT)

                if isinstance(o, np.timedelta64):
                    with tm.assertRaises(TypeError):
                        o + p
                else:
                    self.assertIs(o + p, tslib.NaT)

    def test_sub_pdnat(self):
        # GH 13071
        p = pd.Period('2011-01', freq='M')
        self.assertIs(p - pd.NaT, pd.NaT)
        self.assertIs(pd.NaT - p, pd.NaT)

        p = pd.Period('NaT', freq='M')
        self.assertIs(p - pd.NaT, pd.NaT)
        self.assertIs(pd.NaT - p, pd.NaT)

    def test_sub_offset(self):
        # freq is DateOffset
        for freq in ['A', '2A', '3A']:
            p = Period('2011', freq=freq)
            self.assertEqual(p - offsets.YearEnd(2), Period('2009', freq=freq))

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(365, 'D'),
                      timedelta(365)]:
                with tm.assertRaises(period.IncompatibleFrequency):
                    p - o

        for freq in ['M', '2M', '3M']:
            p = Period('2011-03', freq=freq)
            self.assertEqual(p - offsets.MonthEnd(2),
                             Period('2011-01', freq=freq))
            self.assertEqual(p - offsets.MonthEnd(12),
                             Period('2010-03', freq=freq))

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(365, 'D'),
                      timedelta(365)]:
                with tm.assertRaises(period.IncompatibleFrequency):
                    p - o

        # freq is Tick
        for freq in ['D', '2D', '3D']:
            p = Period('2011-04-01', freq=freq)
            self.assertEqual(p - offsets.Day(5),
                             Period('2011-03-27', freq=freq))
            self.assertEqual(p - offsets.Hour(24),
                             Period('2011-03-31', freq=freq))
            self.assertEqual(p - np.timedelta64(2, 'D'),
                             Period('2011-03-30', freq=freq))
            self.assertEqual(p - np.timedelta64(3600 * 24, 's'),
                             Period('2011-03-31', freq=freq))
            self.assertEqual(p - timedelta(-2),
                             Period('2011-04-03', freq=freq))
            self.assertEqual(p - timedelta(hours=48),
                             Period('2011-03-30', freq=freq))

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(4, 'h'),
                      timedelta(hours=23)]:
                with tm.assertRaises(period.IncompatibleFrequency):
                    p - o

        for freq in ['H', '2H', '3H']:
            p = Period('2011-04-01 09:00', freq=freq)
            self.assertEqual(p - offsets.Day(2),
                             Period('2011-03-30 09:00', freq=freq))
            self.assertEqual(p - offsets.Hour(3),
                             Period('2011-04-01 06:00', freq=freq))
            self.assertEqual(p - np.timedelta64(3, 'h'),
                             Period('2011-04-01 06:00', freq=freq))
            self.assertEqual(p - np.timedelta64(3600, 's'),
                             Period('2011-04-01 08:00', freq=freq))
            self.assertEqual(p - timedelta(minutes=120),
                             Period('2011-04-01 07:00', freq=freq))
            self.assertEqual(p - timedelta(days=4, minutes=180),
                             Period('2011-03-28 06:00', freq=freq))

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(3200, 's'),
                      timedelta(hours=23, minutes=30)]:
                with tm.assertRaises(period.IncompatibleFrequency):
                    p - o

    def test_sub_offset_nat(self):
        # freq is DateOffset
        for freq in ['A', '2A', '3A']:
            p = Period('NaT', freq=freq)
            for o in [offsets.YearEnd(2)]:
                self.assertIs(p - o, tslib.NaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(365, 'D'),
                      timedelta(365)]:
                self.assertIs(p - o, tslib.NaT)

        for freq in ['M', '2M', '3M']:
            p = Period('NaT', freq=freq)
            for o in [offsets.MonthEnd(2), offsets.MonthEnd(12)]:
                self.assertIs(p - o, tslib.NaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(365, 'D'),
                      timedelta(365)]:
                self.assertIs(p - o, tslib.NaT)

        # freq is Tick
        for freq in ['D', '2D', '3D']:
            p = Period('NaT', freq=freq)
            for o in [offsets.Day(5), offsets.Hour(24), np.timedelta64(2, 'D'),
                      np.timedelta64(3600 * 24, 's'), timedelta(-2),
                      timedelta(hours=48)]:
                self.assertIs(p - o, tslib.NaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(4, 'h'),
                      timedelta(hours=23)]:
                self.assertIs(p - o, tslib.NaT)

        for freq in ['H', '2H', '3H']:
            p = Period('NaT', freq=freq)
            for o in [offsets.Day(2), offsets.Hour(3), np.timedelta64(3, 'h'),
                      np.timedelta64(3600, 's'), timedelta(minutes=120),
                      timedelta(days=4, minutes=180)]:
                self.assertIs(p - o, tslib.NaT)

            for o in [offsets.YearBegin(2), offsets.MonthBegin(1),
                      offsets.Minute(), np.timedelta64(3200, 's'),
                      timedelta(hours=23, minutes=30)]:
                self.assertIs(p - o, tslib.NaT)

    def test_nat_ops(self):
        for freq in ['M', '2M', '3M']:
            p = Period('NaT', freq=freq)
            self.assertIs(p + 1, tslib.NaT)
            self.assertIs(1 + p, tslib.NaT)
            self.assertIs(p - 1, tslib.NaT)
            self.assertIs(p - Period('2011-01', freq=freq), tslib.NaT)
            self.assertIs(Period('2011-01', freq=freq) - p, tslib.NaT)

    def test_period_ops_offset(self):
        p = Period('2011-04-01', freq='D')
        result = p + offsets.Day()
        exp = pd.Period('2011-04-02', freq='D')
        self.assertEqual(result, exp)

        result = p - offsets.Day(2)
        exp = pd.Period('2011-03-30', freq='D')
        self.assertEqual(result, exp)

        msg = r"Input cannot be converted to Period\(freq=D\)"
        with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
            p + offsets.Hour(2)

        with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
            p - offsets.Hour(2)


class TestPeriodIndexSeriesMethods(tm.TestCase):
    """ Test PeriodIndex and Period Series Ops consistency """

    def _check(self, values, func, expected):
        idx = pd.PeriodIndex(values)
        result = func(idx)
        if isinstance(expected, pd.Index):
            tm.assert_index_equal(result, expected)
        else:
            # comp op results in bool
            tm.assert_numpy_array_equal(result, expected)

        s = pd.Series(values)
        result = func(s)

        exp = pd.Series(expected, name=values.name)
        tm.assert_series_equal(result, exp)

    def test_pi_ops(self):
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03',
                           '2011-04'], freq='M', name='idx')

        expected = PeriodIndex(['2011-03', '2011-04',
                                '2011-05', '2011-06'], freq='M', name='idx')
        self._check(idx, lambda x: x + 2, expected)
        self._check(idx, lambda x: 2 + x, expected)

        self._check(idx + 2, lambda x: x - 2, idx)
        result = idx - Period('2011-01', freq='M')
        exp = pd.Index([0, 1, 2, 3], name='idx')
        tm.assert_index_equal(result, exp)

        result = Period('2011-01', freq='M') - idx
        exp = pd.Index([0, -1, -2, -3], name='idx')
        tm.assert_index_equal(result, exp)

    def test_pi_ops_errors(self):
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03',
                           '2011-04'], freq='M', name='idx')
        s = pd.Series(idx)

        msg = r"unsupported operand type\(s\)"

        for obj in [idx, s]:
            for ng in ["str", 1.5]:
                with tm.assertRaisesRegexp(TypeError, msg):
                    obj + ng

                with tm.assertRaises(TypeError):
                    # error message differs between PY2 and 3
                    ng + obj

                with tm.assertRaisesRegexp(TypeError, msg):
                    obj - ng

                with tm.assertRaises(TypeError):
                    np.add(obj, ng)

                if _np_version_under1p10:
                    self.assertIs(np.add(ng, obj), NotImplemented)
                else:
                    with tm.assertRaises(TypeError):
                        np.add(ng, obj)

                with tm.assertRaises(TypeError):
                    np.subtract(obj, ng)

                if _np_version_under1p10:
                    self.assertIs(np.subtract(ng, obj), NotImplemented)
                else:
                    with tm.assertRaises(TypeError):
                        np.subtract(ng, obj)

    def test_pi_ops_nat(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT',
                           '2011-04'], freq='M', name='idx')
        expected = PeriodIndex(['2011-03', '2011-04',
                                'NaT', '2011-06'], freq='M', name='idx')
        self._check(idx, lambda x: x + 2, expected)
        self._check(idx, lambda x: 2 + x, expected)
        self._check(idx, lambda x: np.add(x, 2), expected)

        self._check(idx + 2, lambda x: x - 2, idx)
        self._check(idx + 2, lambda x: np.subtract(x, 2), idx)

        # freq with mult
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT',
                           '2011-04'], freq='2M', name='idx')
        expected = PeriodIndex(['2011-07', '2011-08',
                                'NaT', '2011-10'], freq='2M', name='idx')
        self._check(idx, lambda x: x + 3, expected)
        self._check(idx, lambda x: 3 + x, expected)
        self._check(idx, lambda x: np.add(x, 3), expected)

        self._check(idx + 3, lambda x: x - 3, idx)
        self._check(idx + 3, lambda x: np.subtract(x, 3), idx)

    def test_pi_ops_array_int(self):
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT',
                           '2011-04'], freq='M', name='idx')
        f = lambda x: x + np.array([1, 2, 3, 4])
        exp = PeriodIndex(['2011-02', '2011-04', 'NaT',
                           '2011-08'], freq='M', name='idx')
        self._check(idx, f, exp)

        f = lambda x: np.add(x, np.array([4, -1, 1, 2]))
        exp = PeriodIndex(['2011-05', '2011-01', 'NaT',
                           '2011-06'], freq='M', name='idx')
        self._check(idx, f, exp)

        f = lambda x: x - np.array([1, 2, 3, 4])
        exp = PeriodIndex(['2010-12', '2010-12', 'NaT',
                           '2010-12'], freq='M', name='idx')
        self._check(idx, f, exp)

        f = lambda x: np.subtract(x, np.array([3, 2, 3, -2]))
        exp = PeriodIndex(['2010-10', '2010-12', 'NaT',
                           '2011-06'], freq='M', name='idx')
        self._check(idx, f, exp)

    def test_pi_ops_offset(self):
        idx = PeriodIndex(['2011-01-01', '2011-02-01', '2011-03-01',
                           '2011-04-01'], freq='D', name='idx')
        f = lambda x: x + offsets.Day()
        exp = PeriodIndex(['2011-01-02', '2011-02-02', '2011-03-02',
                           '2011-04-02'], freq='D', name='idx')
        self._check(idx, f, exp)

        f = lambda x: x + offsets.Day(2)
        exp = PeriodIndex(['2011-01-03', '2011-02-03', '2011-03-03',
                           '2011-04-03'], freq='D', name='idx')
        self._check(idx, f, exp)

        f = lambda x: x - offsets.Day(2)
        exp = PeriodIndex(['2010-12-30', '2011-01-30', '2011-02-27',
                           '2011-03-30'], freq='D', name='idx')
        self._check(idx, f, exp)

    def test_pi_offset_errors(self):
        idx = PeriodIndex(['2011-01-01', '2011-02-01', '2011-03-01',
                           '2011-04-01'], freq='D', name='idx')
        s = pd.Series(idx)

        # Series op is applied per Period instance, thus error is raised
        # from Period
        msg_idx = r"Input has different freq from PeriodIndex\(freq=D\)"
        msg_s = r"Input cannot be converted to Period\(freq=D\)"
        for obj, msg in [(idx, msg_idx), (s, msg_s)]:
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                obj + offsets.Hour(2)

            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                offsets.Hour(2) + obj

            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                obj - offsets.Hour(2)

    def test_pi_sub_period(self):
        # GH 13071
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03',
                           '2011-04'], freq='M', name='idx')

        result = idx - pd.Period('2012-01', freq='M')
        exp = pd.Index([-12, -11, -10, -9], name='idx')
        tm.assert_index_equal(result, exp)

        result = np.subtract(idx, pd.Period('2012-01', freq='M'))
        tm.assert_index_equal(result, exp)

        result = pd.Period('2012-01', freq='M') - idx
        exp = pd.Index([12, 11, 10, 9], name='idx')
        tm.assert_index_equal(result, exp)

        result = np.subtract(pd.Period('2012-01', freq='M'), idx)
        if _np_version_under1p10:
            self.assertIs(result, NotImplemented)
        else:
            tm.assert_index_equal(result, exp)

        exp = pd.TimedeltaIndex([np.nan, np.nan, np.nan, np.nan], name='idx')
        tm.assert_index_equal(idx - pd.Period('NaT', freq='M'), exp)
        tm.assert_index_equal(pd.Period('NaT', freq='M') - idx, exp)

    def test_pi_sub_pdnat(self):
        # GH 13071
        idx = PeriodIndex(['2011-01', '2011-02', 'NaT',
                           '2011-04'], freq='M', name='idx')
        exp = pd.TimedeltaIndex([pd.NaT] * 4, name='idx')
        tm.assert_index_equal(pd.NaT - idx, exp)
        tm.assert_index_equal(idx - pd.NaT, exp)

    def test_pi_sub_period_nat(self):
        # GH 13071
        idx = PeriodIndex(['2011-01', 'NaT', '2011-03',
                           '2011-04'], freq='M', name='idx')

        result = idx - pd.Period('2012-01', freq='M')
        exp = pd.Index([-12, np.nan, -10, -9], name='idx')
        tm.assert_index_equal(result, exp)

        result = pd.Period('2012-01', freq='M') - idx
        exp = pd.Index([12, np.nan, 10, 9], name='idx')
        tm.assert_index_equal(result, exp)

        exp = pd.TimedeltaIndex([np.nan, np.nan, np.nan, np.nan], name='idx')
        tm.assert_index_equal(idx - pd.Period('NaT', freq='M'), exp)
        tm.assert_index_equal(pd.Period('NaT', freq='M') - idx, exp)

    def test_pi_comp_period(self):
        idx = PeriodIndex(['2011-01', '2011-02', '2011-03',
                           '2011-04'], freq='M', name='idx')

        f = lambda x: x == pd.Period('2011-03', freq='M')
        exp = np.array([False, False, True, False], dtype=np.bool)
        self._check(idx, f, exp)
        f = lambda x: pd.Period('2011-03', freq='M') == x
        self._check(idx, f, exp)

        f = lambda x: x != pd.Period('2011-03', freq='M')
        exp = np.array([True, True, False, True], dtype=np.bool)
        self._check(idx, f, exp)
        f = lambda x: pd.Period('2011-03', freq='M') != x
        self._check(idx, f, exp)

        f = lambda x: pd.Period('2011-03', freq='M') >= x
        exp = np.array([True, True, True, False], dtype=np.bool)
        self._check(idx, f, exp)

        f = lambda x: x > pd.Period('2011-03', freq='M')
        exp = np.array([False, False, False, True], dtype=np.bool)
        self._check(idx, f, exp)

        f = lambda x: pd.Period('2011-03', freq='M') >= x
        exp = np.array([True, True, True, False], dtype=np.bool)
        self._check(idx, f, exp)

    def test_pi_comp_period_nat(self):
        idx = PeriodIndex(['2011-01', 'NaT', '2011-03',
                           '2011-04'], freq='M', name='idx')

        f = lambda x: x == pd.Period('2011-03', freq='M')
        exp = np.array([False, False, True, False], dtype=np.bool)
        self._check(idx, f, exp)
        f = lambda x: pd.Period('2011-03', freq='M') == x
        self._check(idx, f, exp)

        f = lambda x: x == tslib.NaT
        exp = np.array([False, False, False, False], dtype=np.bool)
        self._check(idx, f, exp)
        f = lambda x: tslib.NaT == x
        self._check(idx, f, exp)

        f = lambda x: x != pd.Period('2011-03', freq='M')
        exp = np.array([True, True, False, True], dtype=np.bool)
        self._check(idx, f, exp)
        f = lambda x: pd.Period('2011-03', freq='M') != x
        self._check(idx, f, exp)

        f = lambda x: x != tslib.NaT
        exp = np.array([True, True, True, True], dtype=np.bool)
        self._check(idx, f, exp)
        f = lambda x: tslib.NaT != x
        self._check(idx, f, exp)

        f = lambda x: pd.Period('2011-03', freq='M') >= x
        exp = np.array([True, False, True, False], dtype=np.bool)
        self._check(idx, f, exp)

        f = lambda x: x < pd.Period('2011-03', freq='M')
        exp = np.array([True, False, False, False], dtype=np.bool)
        self._check(idx, f, exp)

        f = lambda x: x > tslib.NaT
        exp = np.array([False, False, False, False], dtype=np.bool)
        self._check(idx, f, exp)

        f = lambda x: tslib.NaT >= x
        exp = np.array([False, False, False, False], dtype=np.bool)
        self._check(idx, f, exp)


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
        self.assert_numpy_array_equal(rng._values, exp)
        self.assert_numpy_array_equal(rng.asi8, exp)

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
        with tm.assertRaises(period.IncompatibleFrequency):
            self.january1 == self.day

    def test_notEqual(self):
        self.assertNotEqual(self.january1, 1)
        self.assertNotEqual(self.january1, self.february)

    def test_greater(self):
        self.assertTrue(self.february > self.january1)

    def test_greater_Raises_Value(self):
        with tm.assertRaises(period.IncompatibleFrequency):
            self.january1 > self.day

    def test_greater_Raises_Type(self):
        with tm.assertRaises(TypeError):
            self.january1 > 1

    def test_greaterEqual(self):
        self.assertTrue(self.january1 >= self.january2)

    def test_greaterEqual_Raises_Value(self):
        with tm.assertRaises(period.IncompatibleFrequency):
            self.january1 >= self.day

        with tm.assertRaises(TypeError):
            print(self.january1 >= 1)

    def test_smallerEqual(self):
        self.assertTrue(self.january1 <= self.january2)

    def test_smallerEqual_Raises_Value(self):
        with tm.assertRaises(period.IncompatibleFrequency):
            self.january1 <= self.day

    def test_smallerEqual_Raises_Type(self):
        with tm.assertRaises(TypeError):
            self.january1 <= 1

    def test_smaller(self):
        self.assertTrue(self.january1 < self.february)

    def test_smaller_Raises_Value(self):
        with tm.assertRaises(period.IncompatibleFrequency):
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
        for left, right in [(p_nat, p), (p, p_nat), (p_nat, p_nat), (nat, t),
                            (t, nat), (nat, nat)]:
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
            self.assert_numpy_array_equal(p == base, exp)

            exp = np.array([True, False, True, True])
            self.assert_numpy_array_equal(base != p, exp)
            self.assert_numpy_array_equal(p != base, exp)

            exp = np.array([False, False, True, True])
            self.assert_numpy_array_equal(base > p, exp)
            self.assert_numpy_array_equal(p < base, exp)

            exp = np.array([True, False, False, False])
            self.assert_numpy_array_equal(base < p, exp)
            self.assert_numpy_array_equal(p > base, exp)

            exp = np.array([False, True, True, True])
            self.assert_numpy_array_equal(base >= p, exp)
            self.assert_numpy_array_equal(p <= base, exp)

            exp = np.array([True, True, False, False])
            self.assert_numpy_array_equal(base <= p, exp)
            self.assert_numpy_array_equal(p >= base, exp)

            idx = PeriodIndex(['2011-02', '2011-01', '2011-03',
                               '2011-05'], freq=freq)

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
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                base <= Period('2011', freq='A')

            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                Period('2011', freq='A') >= base

            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                idx = PeriodIndex(['2011', '2012', '2013', '2014'], freq='A')
                base <= idx

            # different mult
            msg = "Input has different freq=4M from PeriodIndex"
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                base <= Period('2011', freq='4M')

            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                Period('2011', freq='4M') >= base

            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                idx = PeriodIndex(['2011', '2012', '2013', '2014'], freq='4M')
                base <= idx

    def test_pi_nat_comp(self):
        for freq in ['M', '2M', '3M']:
            idx1 = PeriodIndex(
                ['2011-01', '2011-02', 'NaT', '2011-05'], freq=freq)

            result = idx1 > Period('2011-02', freq=freq)
            exp = np.array([False, False, False, True])
            self.assert_numpy_array_equal(result, exp)
            result = Period('2011-02', freq=freq) < idx1
            self.assert_numpy_array_equal(result, exp)

            result = idx1 == Period('NaT', freq=freq)
            exp = np.array([False, False, False, False])
            self.assert_numpy_array_equal(result, exp)
            result = Period('NaT', freq=freq) == idx1
            self.assert_numpy_array_equal(result, exp)

            result = idx1 != Period('NaT', freq=freq)
            exp = np.array([True, True, True, True])
            self.assert_numpy_array_equal(result, exp)
            result = Period('NaT', freq=freq) != idx1
            self.assert_numpy_array_equal(result, exp)

            idx2 = PeriodIndex(['2011-02', '2011-01', '2011-04',
                                'NaT'], freq=freq)
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

            diff = PeriodIndex(['2011-02', '2011-01', '2011-04',
                                'NaT'], freq='4M')
            msg = "Input has different freq=4M from PeriodIndex"
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                idx1 > diff

            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                idx1 == diff


class TestSeriesPeriod(tm.TestCase):

    def setUp(self):
        self.series = Series(period_range('2000-01-01', periods=10, freq='D'))

    def test_auto_conversion(self):
        series = Series(list(period_range('2000-01-01', periods=10, freq='D')))
        self.assertEqual(series.dtype, 'object')

        series = pd.Series([pd.Period('2011-01-01', freq='D'),
                            pd.Period('2011-02-01', freq='D')])
        self.assertEqual(series.dtype, 'object')

    def test_getitem(self):
        self.assertEqual(self.series[1], pd.Period('2000-01-02', freq='D'))

        result = self.series[[2, 4]]
        exp = pd.Series([pd.Period('2000-01-03', freq='D'),
                         pd.Period('2000-01-05', freq='D')],
                        index=[2, 4])
        self.assert_series_equal(result, exp)
        self.assertEqual(result.dtype, 'object')

    def test_constructor_cant_cast_period(self):
        with tm.assertRaises(TypeError):
            Series(period_range('2000-01-01', periods=10, freq='D'),
                   dtype=float)

    def test_constructor_cast_object(self):
        s = Series(period_range('1/1/2000', periods=10), dtype=object)
        exp = Series(period_range('1/1/2000', periods=10))
        tm.assert_series_equal(s, exp)

    def test_isnull(self):
        # GH 13737
        s = Series([pd.Period('2011-01', freq='M'),
                    pd.Period('NaT', freq='M')])
        tm.assert_series_equal(s.isnull(), Series([False, True]))
        tm.assert_series_equal(s.notnull(), Series([True, False]))

    def test_fillna(self):
        # GH 13737
        s = Series([pd.Period('2011-01', freq='M'),
                    pd.Period('NaT', freq='M')])

        res = s.fillna(pd.Period('2012-01', freq='M'))
        exp = Series([pd.Period('2011-01', freq='M'),
                      pd.Period('2012-01', freq='M')])
        tm.assert_series_equal(res, exp)
        self.assertEqual(res.dtype, 'object')

        res = s.fillna('XXX')
        exp = Series([pd.Period('2011-01', freq='M'), 'XXX'])
        tm.assert_series_equal(res, exp)
        self.assertEqual(res.dtype, 'object')

    def test_dropna(self):
        # GH 13737
        s = Series([pd.Period('2011-01', freq='M'),
                    pd.Period('NaT', freq='M')])
        tm.assert_series_equal(s.dropna(),
                               Series([pd.Period('2011-01', freq='M')]))

    def test_series_comparison_scalars(self):
        val = pd.Period('2000-01-04', freq='D')
        result = self.series > val
        expected = pd.Series([x > val for x in self.series])
        tm.assert_series_equal(result, expected)

        val = self.series[5]
        result = self.series > val
        expected = pd.Series([x > val for x in self.series])
        tm.assert_series_equal(result, expected)

    def test_between(self):
        left, right = self.series[[2, 7]]
        result = self.series.between(left, right)
        expected = (self.series >= left) & (self.series <= right)
        tm.assert_series_equal(result, expected)

    # ---------------------------------------------------------------------
    # NaT support

    """
    # ToDo: Enable when support period dtype
    def test_NaT_scalar(self):
        series = Series([0, 1000, 2000, iNaT], dtype='period[D]')

        val = series[3]
        self.assertTrue(isnull(val))

        series[2] = val
        self.assertTrue(isnull(series[2]))

    def test_NaT_cast(self):
        result = Series([np.nan]).astype('period[D]')
        expected = Series([NaT])
        tm.assert_series_equal(result, expected)
    """

    def test_set_none_nan(self):
        # currently Period is stored as object dtype, not as NaT
        self.series[3] = None
        self.assertIs(self.series[3], None)

        self.series[3:5] = None
        self.assertIs(self.series[4], None)

        self.series[5] = np.nan
        self.assertTrue(np.isnan(self.series[5]))

        self.series[5:7] = np.nan
        self.assertTrue(np.isnan(self.series[6]))

    def test_intercept_astype_object(self):
        expected = self.series.astype('object')

        df = DataFrame({'a': self.series,
                        'b': np.random.randn(len(self.series))})

        result = df.values.squeeze()
        self.assertTrue((result[:, 0] == expected.values).all())

        df = DataFrame({'a': self.series, 'b': ['foo'] * len(self.series)})

        result = df.values.squeeze()
        self.assertTrue((result[:, 0] == expected.values).all())

    def test_ops_series_timedelta(self):
        # GH 13043
        s = pd.Series([pd.Period('2015-01-01', freq='D'),
                       pd.Period('2015-01-02', freq='D')], name='xxx')
        self.assertEqual(s.dtype, object)

        exp = pd.Series([pd.Period('2015-01-02', freq='D'),
                         pd.Period('2015-01-03', freq='D')], name='xxx')
        tm.assert_series_equal(s + pd.Timedelta('1 days'), exp)
        tm.assert_series_equal(pd.Timedelta('1 days') + s, exp)

        tm.assert_series_equal(s + pd.tseries.offsets.Day(), exp)
        tm.assert_series_equal(pd.tseries.offsets.Day() + s, exp)

    def test_ops_series_period(self):
        # GH 13043
        s = pd.Series([pd.Period('2015-01-01', freq='D'),
                       pd.Period('2015-01-02', freq='D')], name='xxx')
        self.assertEqual(s.dtype, object)

        p = pd.Period('2015-01-10', freq='D')
        # dtype will be object because of original dtype
        exp = pd.Series([9, 8], name='xxx', dtype=object)
        tm.assert_series_equal(p - s, exp)
        tm.assert_series_equal(s - p, -exp)

        s2 = pd.Series([pd.Period('2015-01-05', freq='D'),
                        pd.Period('2015-01-04', freq='D')], name='xxx')
        self.assertEqual(s2.dtype, object)

        exp = pd.Series([4, 2], name='xxx', dtype=object)
        tm.assert_series_equal(s2 - s, exp)
        tm.assert_series_equal(s - s2, -exp)

    def test_comp_series_period_scalar(self):
        # GH 13200
        for freq in ['M', '2M', '3M']:
            base = Series([Period(x, freq=freq) for x in
                           ['2011-01', '2011-02', '2011-03', '2011-04']])
            p = Period('2011-02', freq=freq)

            exp = pd.Series([False, True, False, False])
            tm.assert_series_equal(base == p, exp)
            tm.assert_series_equal(p == base, exp)

            exp = pd.Series([True, False, True, True])
            tm.assert_series_equal(base != p, exp)
            tm.assert_series_equal(p != base, exp)

            exp = pd.Series([False, False, True, True])
            tm.assert_series_equal(base > p, exp)
            tm.assert_series_equal(p < base, exp)

            exp = pd.Series([True, False, False, False])
            tm.assert_series_equal(base < p, exp)
            tm.assert_series_equal(p > base, exp)

            exp = pd.Series([False, True, True, True])
            tm.assert_series_equal(base >= p, exp)
            tm.assert_series_equal(p <= base, exp)

            exp = pd.Series([True, True, False, False])
            tm.assert_series_equal(base <= p, exp)
            tm.assert_series_equal(p >= base, exp)

            # different base freq
            msg = "Input has different freq=A-DEC from Period"
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                base <= Period('2011', freq='A')

            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                Period('2011', freq='A') >= base

    def test_comp_series_period_series(self):
        # GH 13200
        for freq in ['M', '2M', '3M']:
            base = Series([Period(x, freq=freq) for x in
                           ['2011-01', '2011-02', '2011-03', '2011-04']])

            s = Series([Period(x, freq=freq) for x in
                        ['2011-02', '2011-01', '2011-03', '2011-05']])

            exp = Series([False, False, True, False])
            tm.assert_series_equal(base == s, exp)

            exp = Series([True, True, False, True])
            tm.assert_series_equal(base != s, exp)

            exp = Series([False, True, False, False])
            tm.assert_series_equal(base > s, exp)

            exp = Series([True, False, False, True])
            tm.assert_series_equal(base < s, exp)

            exp = Series([False, True, True, False])
            tm.assert_series_equal(base >= s, exp)

            exp = Series([True, False, True, True])
            tm.assert_series_equal(base <= s, exp)

            s2 = Series([Period(x, freq='A') for x in
                         ['2011', '2011', '2011', '2011']])

            # different base freq
            msg = "Input has different freq=A-DEC from Period"
            with tm.assertRaisesRegexp(period.IncompatibleFrequency, msg):
                base <= s2

    def test_comp_series_period_object(self):
        # GH 13200
        base = Series([Period('2011', freq='A'), Period('2011-02', freq='M'),
                       Period('2013', freq='A'), Period('2011-04', freq='M')])

        s = Series([Period('2012', freq='A'), Period('2011-01', freq='M'),
                    Period('2013', freq='A'), Period('2011-05', freq='M')])

        exp = Series([False, False, True, False])
        tm.assert_series_equal(base == s, exp)

        exp = Series([True, True, False, True])
        tm.assert_series_equal(base != s, exp)

        exp = Series([False, True, False, False])
        tm.assert_series_equal(base > s, exp)

        exp = Series([True, False, False, True])
        tm.assert_series_equal(base < s, exp)

        exp = Series([False, True, True, False])
        tm.assert_series_equal(base >= s, exp)

        exp = Series([True, False, True, True])
        tm.assert_series_equal(base <= s, exp)

    def test_ops_frame_period(self):
        # GH 13043
        df = pd.DataFrame({'A': [pd.Period('2015-01', freq='M'),
                                 pd.Period('2015-02', freq='M')],
                           'B': [pd.Period('2014-01', freq='M'),
                                 pd.Period('2014-02', freq='M')]})
        self.assertEqual(df['A'].dtype, object)
        self.assertEqual(df['B'].dtype, object)

        p = pd.Period('2015-03', freq='M')
        # dtype will be object because of original dtype
        exp = pd.DataFrame({'A': np.array([2, 1], dtype=object),
                            'B': np.array([14, 13], dtype=object)})
        tm.assert_frame_equal(p - df, exp)
        tm.assert_frame_equal(df - p, -exp)

        df2 = pd.DataFrame({'A': [pd.Period('2015-05', freq='M'),
                                  pd.Period('2015-06', freq='M')],
                            'B': [pd.Period('2015-05', freq='M'),
                                  pd.Period('2015-06', freq='M')]})
        self.assertEqual(df2['A'].dtype, object)
        self.assertEqual(df2['B'].dtype, object)

        exp = pd.DataFrame({'A': np.array([4, 4], dtype=object),
                            'B': np.array([16, 16], dtype=object)})
        tm.assert_frame_equal(df2 - df, exp)
        tm.assert_frame_equal(df - df2, -exp)


class TestPeriodField(tm.TestCase):
    def test_get_period_field_raises_on_out_of_range(self):
        self.assertRaises(ValueError, _period.get_period_field, -1, 0, 0)

    def test_get_period_field_array_raises_on_out_of_range(self):
        self.assertRaises(ValueError, _period.get_period_field_arr, -1,
                          np.empty(1), 0)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
