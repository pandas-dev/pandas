import numpy as np
from datetime import datetime, date, timedelta

import pandas as pd
import pandas.util.testing as tm
import pandas.tseries.period as period
from pandas.compat import text_type, iteritems
from pandas.compat.numpy import np_datetime64_compat
from pandas import Period, Timedelta, Timestamp, tslib, offsets, _period
from pandas.tseries.frequencies import DAYS, MONTHS, _period_code_map


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

    def test_badinput(self):
        self.assertRaises(ValueError, Period, '-2000', 'A')
        self.assertRaises(tslib.DateParseError, Period, '0', 'A')
        self.assertRaises(tslib.DateParseError, Period, '1/1/-2000', 'A')

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

    def test_round_trip(self):

        p = Period('2000Q1')
        new_p = self.round_trip_pickle(p)
        self.assertEqual(new_p, p)


class TestPeriodField(tm.TestCase):

    def test_get_period_field_raises_on_out_of_range(self):
        self.assertRaises(ValueError, _period.get_period_field, -1, 0, 0)

    def test_get_period_field_array_raises_on_out_of_range(self):
        self.assertRaises(ValueError, _period.get_period_field_arr, -1,
                          np.empty(1), 0)


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
