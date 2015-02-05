from datetime import datetime, time, timedelta
from pandas.compat import range
import sys
import os

import nose

import numpy as np

from pandas import Index, DatetimeIndex, Timestamp, Series, date_range, period_range

import pandas.tseries.frequencies as frequencies
from pandas.tseries.tools import to_datetime

import pandas.tseries.offsets as offsets
from pandas.tseries.period import PeriodIndex
import pandas.compat as compat

import pandas.util.testing as tm
from pandas import Timedelta

def test_to_offset_multiple():
    freqstr = '2h30min'
    freqstr2 = '2h 30min'

    result = frequencies.to_offset(freqstr)
    assert(result == frequencies.to_offset(freqstr2))
    expected = offsets.Minute(150)
    assert(result == expected)

    freqstr = '2h30min15s'
    result = frequencies.to_offset(freqstr)
    expected = offsets.Second(150 * 60 + 15)
    assert(result == expected)

    freqstr = '2h 60min'
    result = frequencies.to_offset(freqstr)
    expected = offsets.Hour(3)
    assert(result == expected)

    freqstr = '15l500u'
    result = frequencies.to_offset(freqstr)
    expected = offsets.Micro(15500)
    assert(result == expected)

    freqstr = '10s75L'
    result = frequencies.to_offset(freqstr)
    expected = offsets.Milli(10075)
    assert(result == expected)

    freqstr = '2800N'
    result = frequencies.to_offset(freqstr)
    expected = offsets.Nano(2800)
    assert(result == expected)

    # malformed
    try:
        frequencies.to_offset('2h20m')
    except ValueError:
        pass
    else:
        assert(False)


def test_to_offset_negative():
    freqstr = '-1S'
    result = frequencies.to_offset(freqstr)
    assert(result.n == -1)

    freqstr = '-5min10s'
    result = frequencies.to_offset(freqstr)
    assert(result.n == -310)


def test_to_offset_leading_zero():
    freqstr = '00H 00T 01S'
    result = frequencies.to_offset(freqstr)
    assert(result.n == 1)

    freqstr = '-00H 03T 14S'
    result = frequencies.to_offset(freqstr)
    assert(result.n == -194)


def test_to_offset_pd_timedelta():
    # Tests for #9064
    td = Timedelta(days=1, seconds=1)
    result = frequencies.to_offset(td)
    expected = offsets.Second(86401)
    assert(expected==result)

    td = Timedelta(days=-1, seconds=1)
    result = frequencies.to_offset(td)
    expected = offsets.Second(-86399)
    assert(expected==result)

    td = Timedelta(hours=1, minutes=10)
    result = frequencies.to_offset(td)
    expected = offsets.Minute(70)
    assert(expected==result)

    td = Timedelta(hours=1, minutes=-10)
    result = frequencies.to_offset(td)
    expected = offsets.Minute(50)
    assert(expected==result)

    td = Timedelta(weeks=1)
    result = frequencies.to_offset(td)
    expected = offsets.Day(7)
    assert(expected==result)

    td1 = Timedelta(hours=1)
    result1 = frequencies.to_offset(td1)
    result2 = frequencies.to_offset('60min')
    assert(result1 == result2)

    td = Timedelta(microseconds=1)
    result = frequencies.to_offset(td)
    expected = offsets.Micro(1)
    assert(expected == result)

    td = Timedelta(microseconds=0)
    tm.assertRaises(ValueError, lambda: frequencies.to_offset(td))


def test_anchored_shortcuts():
    result = frequencies.to_offset('W')
    expected = frequencies.to_offset('W-SUN')
    assert(result == expected)

    result = frequencies.to_offset('Q')
    expected = frequencies.to_offset('Q-DEC')
    assert(result == expected)


_dti = DatetimeIndex


class TestFrequencyInference(tm.TestCase):

    def test_raise_if_period_index(self):
        index = PeriodIndex(start="1/1/1990", periods=20, freq="M")
        self.assertRaises(TypeError, frequencies.infer_freq, index)

    def test_raise_if_too_few(self):
        index = _dti(['12/31/1998', '1/3/1999'])
        self.assertRaises(ValueError, frequencies.infer_freq, index)

    def test_business_daily(self):
        index = _dti(['12/31/1998', '1/3/1999', '1/4/1999'])
        self.assertEqual(frequencies.infer_freq(index), 'B')

    def test_day(self):
        self._check_tick(timedelta(1), 'D')

    def test_day_corner(self):
        index = _dti(['1/1/2000', '1/2/2000', '1/3/2000'])
        self.assertEqual(frequencies.infer_freq(index), 'D')

    def test_non_datetimeindex(self):
        dates = to_datetime(['1/1/2000', '1/2/2000', '1/3/2000'])
        self.assertEqual(frequencies.infer_freq(dates), 'D')

    def test_hour(self):
        self._check_tick(timedelta(hours=1), 'H')

    def test_minute(self):
        self._check_tick(timedelta(minutes=1), 'T')

    def test_second(self):
        self._check_tick(timedelta(seconds=1), 'S')

    def test_millisecond(self):
        self._check_tick(timedelta(microseconds=1000), 'L')

    def test_microsecond(self):
        self._check_tick(timedelta(microseconds=1), 'U')

    def test_nanosecond(self):
        self._check_tick(np.timedelta64(1, 'ns'), 'N')

    def _check_tick(self, base_delta, code):
        b = Timestamp(datetime.now())
        for i in range(1, 5):
            inc = base_delta * i
            index = _dti([b + inc * j for j in range(3)])
            if i > 1:
                exp_freq = '%d%s' % (i, code)
            else:
                exp_freq = code
            self.assertEqual(frequencies.infer_freq(index), exp_freq)

        index = _dti([b + base_delta * 7] +
                     [b + base_delta * j for j in range(3)])
        self.assertIsNone(frequencies.infer_freq(index))

        index = _dti([b + base_delta * j for j in range(3)] +
                     [b + base_delta * 7])

        self.assertIsNone(frequencies.infer_freq(index))

    def test_weekly(self):
        days = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']

        for day in days:
            self._check_generated_range('1/1/2000', 'W-%s' % day)

    def test_week_of_month(self):
        days = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']

        for day in days:
            for i in range(1, 5):
                self._check_generated_range('1/1/2000', 'WOM-%d%s' % (i, day))

    def test_fifth_week_of_month(self):
        # Only supports freq up to WOM-4. See #9425
        func = lambda: date_range('2014-01-01', freq='WOM-5MON')
        self.assertRaises(ValueError, func)

    def test_fifth_week_of_month_infer(self):
        # Only attempts to infer up to WOM-4. See #9425
        index = DatetimeIndex(["2014-03-31", "2014-06-30", "2015-03-30"])
        assert frequencies.infer_freq(index) is None

    def test_week_of_month_fake(self):
        #All of these dates are on same day of week and are 4 or 5 weeks apart
        index = DatetimeIndex(["2013-08-27","2013-10-01","2013-10-29","2013-11-26"])
        assert frequencies.infer_freq(index) != 'WOM-4TUE'

    def test_monthly(self):
        self._check_generated_range('1/1/2000', 'M')

    def test_monthly_ambiguous(self):
        rng = _dti(['1/31/2000', '2/29/2000', '3/31/2000'])
        self.assertEqual(rng.inferred_freq, 'M')

    def test_business_monthly(self):
        self._check_generated_range('1/1/2000', 'BM')

    def test_business_start_monthly(self):
        self._check_generated_range('1/1/2000', 'BMS')

    def test_quarterly(self):
        for month in ['JAN', 'FEB', 'MAR']:
            self._check_generated_range('1/1/2000', 'Q-%s' % month)

    def test_annual(self):
        for month in MONTHS:
            self._check_generated_range('1/1/2000', 'A-%s' % month)

    def test_business_annual(self):
        for month in MONTHS:
            self._check_generated_range('1/1/2000', 'BA-%s' % month)

    def test_annual_ambiguous(self):
        rng = _dti(['1/31/2000', '1/31/2001', '1/31/2002'])
        self.assertEqual(rng.inferred_freq, 'A-JAN')

    def _check_generated_range(self, start, freq):
        freq = freq.upper()

        gen = date_range(start, periods=7, freq=freq)
        index = _dti(gen.values)
        if not freq.startswith('Q-'):
            self.assertEqual(frequencies.infer_freq(index), gen.freqstr)
        else:
            inf_freq = frequencies.infer_freq(index)
            self.assertTrue((inf_freq == 'Q-DEC' and
                             gen.freqstr in ('Q', 'Q-DEC', 'Q-SEP', 'Q-JUN',
                                             'Q-MAR'))
                            or
                            (inf_freq == 'Q-NOV' and
                             gen.freqstr in ('Q-NOV', 'Q-AUG', 'Q-MAY', 'Q-FEB'))
                            or
                            (inf_freq == 'Q-OCT' and
                             gen.freqstr in ('Q-OCT', 'Q-JUL', 'Q-APR', 'Q-JAN')))

        gen = date_range(start, periods=5, freq=freq)
        index = _dti(gen.values)
        if not freq.startswith('Q-'):
            self.assertEqual(frequencies.infer_freq(index), gen.freqstr)
        else:
            inf_freq = frequencies.infer_freq(index)
            self.assertTrue((inf_freq == 'Q-DEC' and
                             gen.freqstr in ('Q', 'Q-DEC', 'Q-SEP', 'Q-JUN',
                                             'Q-MAR'))
                            or
                            (inf_freq == 'Q-NOV' and
                             gen.freqstr in ('Q-NOV', 'Q-AUG', 'Q-MAY', 'Q-FEB'))
                            or
                            (inf_freq == 'Q-OCT' and
                             gen.freqstr in ('Q-OCT', 'Q-JUL', 'Q-APR', 'Q-JAN')))

    def test_infer_freq(self):
        rng = period_range('1959Q2', '2009Q3', freq='Q')
        rng = Index(rng.to_timestamp('D', how='e').asobject)
        self.assertEqual(rng.inferred_freq, 'Q-DEC')

        rng = period_range('1959Q2', '2009Q3', freq='Q-NOV')
        rng = Index(rng.to_timestamp('D', how='e').asobject)
        self.assertEqual(rng.inferred_freq, 'Q-NOV')

        rng = period_range('1959Q2', '2009Q3', freq='Q-OCT')
        rng = Index(rng.to_timestamp('D', how='e').asobject)
        self.assertEqual(rng.inferred_freq, 'Q-OCT')

    def test_infer_freq_tz(self):

        freqs = {'AS-JAN': ['2009-01-01', '2010-01-01', '2011-01-01', '2012-01-01'],
                 'Q-OCT': ['2009-01-31', '2009-04-30', '2009-07-31', '2009-10-31'],
                 'M': ['2010-11-30', '2010-12-31', '2011-01-31', '2011-02-28'],
                 'W-SAT': ['2010-12-25', '2011-01-01', '2011-01-08', '2011-01-15'],
                 'D': ['2011-01-01', '2011-01-02', '2011-01-03', '2011-01-04'],
                 'H': ['2011-12-31 22:00', '2011-12-31 23:00', '2012-01-01 00:00', '2012-01-01 01:00']
        }

        # GH 7310
        for tz in [None, 'Australia/Sydney', 'Asia/Tokyo', 'Europe/Paris',
                   'US/Pacific', 'US/Eastern']:
            for expected, dates in compat.iteritems(freqs):
                idx = DatetimeIndex(dates, tz=tz)
                self.assertEqual(idx.inferred_freq, expected)

    def test_infer_freq_tz_transition(self):
        # Tests for #8772
        date_pairs = [['2013-11-02', '2013-11-5'], #Fall DST
                      ['2014-03-08', '2014-03-11'], #Spring DST
                      ['2014-01-01', '2014-01-03']] #Regular Time
        freqs = ['3H', '10T', '3601S', '3600001L', '3600000001U', '3600000000001N']

        for tz in [None, 'Australia/Sydney', 'Asia/Tokyo', 'Europe/Paris',
                   'US/Pacific', 'US/Eastern']:
            for date_pair in date_pairs:
                for freq in freqs:
                    idx = date_range(date_pair[0], date_pair[1], freq=freq, tz=tz)
                    print(idx)
                    self.assertEqual(idx.inferred_freq, freq)

        index = date_range("2013-11-03", periods=5, freq="3H").tz_localize("America/Chicago")
        self.assertIsNone(index.inferred_freq)

    def test_infer_freq_businesshour(self):
        # GH 7905
        idx = DatetimeIndex(['2014-07-01 09:00', '2014-07-01 10:00', '2014-07-01 11:00',
                             '2014-07-01 12:00', '2014-07-01 13:00', '2014-07-01 14:00'])
        # hourly freq in a day must result in 'H'
        self.assertEqual(idx.inferred_freq, 'H')

        idx = DatetimeIndex(['2014-07-01 09:00', '2014-07-01 10:00', '2014-07-01 11:00',
                             '2014-07-01 12:00', '2014-07-01 13:00', '2014-07-01 14:00',
                             '2014-07-01 15:00', '2014-07-01 16:00',
                             '2014-07-02 09:00', '2014-07-02 10:00', '2014-07-02 11:00'])
        self.assertEqual(idx.inferred_freq, 'BH')

        idx = DatetimeIndex(['2014-07-04 09:00', '2014-07-04 10:00', '2014-07-04 11:00',
                             '2014-07-04 12:00', '2014-07-04 13:00', '2014-07-04 14:00',
                             '2014-07-04 15:00', '2014-07-04 16:00',
                             '2014-07-07 09:00', '2014-07-07 10:00', '2014-07-07 11:00'])
        self.assertEqual(idx.inferred_freq, 'BH')

        idx = DatetimeIndex(['2014-07-04 09:00', '2014-07-04 10:00', '2014-07-04 11:00',
                             '2014-07-04 12:00', '2014-07-04 13:00', '2014-07-04 14:00',
                             '2014-07-04 15:00', '2014-07-04 16:00',
                             '2014-07-07 09:00', '2014-07-07 10:00', '2014-07-07 11:00',
                             '2014-07-07 12:00', '2014-07-07 13:00', '2014-07-07 14:00',
                             '2014-07-07 15:00', '2014-07-07 16:00',
                             '2014-07-08 09:00', '2014-07-08 10:00', '2014-07-08 11:00',
                             '2014-07-08 12:00', '2014-07-08 13:00', '2014-07-08 14:00',
                             '2014-07-08 15:00', '2014-07-08 16:00'])
        self.assertEqual(idx.inferred_freq, 'BH')

    def test_not_monotonic(self):
        rng = _dti(['1/31/2000', '1/31/2001', '1/31/2002'])
        rng = rng[::-1]
        self.assertIsNone(rng.inferred_freq)

    def test_non_datetimeindex(self):
        rng = _dti(['1/31/2000', '1/31/2001', '1/31/2002'])

        vals = rng.to_pydatetime()

        result = frequencies.infer_freq(vals)
        self.assertEqual(result, rng.inferred_freq)

    def test_invalid_index_types(self):

        # test all index types
        for i in [ tm.makeIntIndex(10),
                   tm.makeFloatIndex(10),
                   tm.makePeriodIndex(10) ]:
            self.assertRaises(TypeError, lambda : frequencies.infer_freq(i))

        for i in [ tm.makeStringIndex(10),
                   tm.makeUnicodeIndex(10) ]:
            self.assertRaises(ValueError, lambda : frequencies.infer_freq(i))

    def test_string_datetimelike_compat(self):

        # GH 6463
        expected = frequencies.infer_freq(['2004-01', '2004-02', '2004-03', '2004-04'])
        result = frequencies.infer_freq(Index(['2004-01', '2004-02', '2004-03', '2004-04']))
        self.assertEqual(result,expected)

    def test_series(self):

        # GH6407
        # inferring series

        # invalid type of Series
        for s in [ Series(np.arange(10)),
                   Series(np.arange(10.))]:
            self.assertRaises(TypeError, lambda : frequencies.infer_freq(s))

        # a non-convertible string
        self.assertRaises(ValueError, lambda : frequencies.infer_freq(Series(['foo','bar'])))

        # cannot infer on PeriodIndex
        for freq in [None, 'L', 'Y']:
            s = Series(period_range('2013',periods=10,freq=freq))
            self.assertRaises(TypeError, lambda : frequencies.infer_freq(s))

        # DateTimeIndex
        for freq in ['M', 'L', 'S']:
            s = Series(date_range('20130101',periods=10,freq=freq))
            inferred = frequencies.infer_freq(s)
            self.assertEqual(inferred,freq)

        s = Series(date_range('20130101','20130110'))
        inferred = frequencies.infer_freq(s)
        self.assertEqual(inferred,'D')

MONTHS = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP',
          'OCT', 'NOV', 'DEC']


def test_is_superperiod_subperiod():
    assert(frequencies.is_superperiod(offsets.YearEnd(), offsets.MonthEnd()))
    assert(frequencies.is_subperiod(offsets.MonthEnd(), offsets.YearEnd()))

    assert(frequencies.is_superperiod(offsets.Hour(), offsets.Minute()))
    assert(frequencies.is_subperiod(offsets.Minute(), offsets.Hour()))

    assert(frequencies.is_superperiod(offsets.Second(), offsets.Milli()))
    assert(frequencies.is_subperiod(offsets.Milli(), offsets.Second()))

    assert(frequencies.is_superperiod(offsets.Milli(), offsets.Micro()))
    assert(frequencies.is_subperiod(offsets.Micro(), offsets.Milli()))

    assert(frequencies.is_superperiod(offsets.Micro(), offsets.Nano()))
    assert(frequencies.is_subperiod(offsets.Nano(), offsets.Micro()))


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
