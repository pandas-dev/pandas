from datetime import datetime, time, timedelta
from pandas.compat import range
import sys
import os

import nose

import numpy as np

from pandas import Index, DatetimeIndex, Timestamp, Series, date_range, period_range

from pandas.tseries.frequencies import to_offset, infer_freq
from pandas.tseries.tools import to_datetime
import pandas.tseries.frequencies as fmod
import pandas.tseries.offsets as offsets
from pandas.tseries.period import PeriodIndex

import pandas.lib as lib

from pandas import _np_version_under1p7
import pandas.util.testing as tm

def test_to_offset_multiple():
    freqstr = '2h30min'
    freqstr2 = '2h 30min'

    result = to_offset(freqstr)
    assert(result == to_offset(freqstr2))
    expected = offsets.Minute(150)
    assert(result == expected)

    freqstr = '2h30min15s'
    result = to_offset(freqstr)
    expected = offsets.Second(150 * 60 + 15)
    assert(result == expected)

    freqstr = '2h 60min'
    result = to_offset(freqstr)
    expected = offsets.Hour(3)
    assert(result == expected)

    freqstr = '15l500u'
    result = to_offset(freqstr)
    expected = offsets.Micro(15500)
    assert(result == expected)

    freqstr = '10s75L'
    result = to_offset(freqstr)
    expected = offsets.Milli(10075)
    assert(result == expected)

    if not _np_version_under1p7:
        freqstr = '2800N'
        result = to_offset(freqstr)
        expected = offsets.Nano(2800)
        assert(result == expected)

    # malformed
    try:
        to_offset('2h20m')
    except ValueError:
        pass
    else:
        assert(False)


def test_to_offset_negative():
    freqstr = '-1S'
    result = to_offset(freqstr)
    assert(result.n == -1)

    freqstr = '-5min10s'
    result = to_offset(freqstr)
    assert(result.n == -310)


def test_to_offset_leading_zero():
    freqstr = '00H 00T 01S'
    result = to_offset(freqstr)
    assert(result.n == 1)

    freqstr = '-00H 03T 14S'
    result = to_offset(freqstr)
    assert(result.n == -194)


def test_anchored_shortcuts():
    result = to_offset('W')
    expected = to_offset('W-SUN')
    assert(result == expected)

    result = to_offset('Q')
    expected = to_offset('Q-DEC')
    assert(result == expected)


_dti = DatetimeIndex


class TestFrequencyInference(tm.TestCase):

    def test_raise_if_period_index(self):
        index = PeriodIndex(start="1/1/1990", periods=20, freq="M")
        self.assertRaises(TypeError, infer_freq, index)

    def test_raise_if_too_few(self):
        index = _dti(['12/31/1998', '1/3/1999'])
        self.assertRaises(ValueError, infer_freq, index)

    def test_business_daily(self):
        index = _dti(['12/31/1998', '1/3/1999', '1/4/1999'])
        self.assertEqual(infer_freq(index), 'B')

    def test_day(self):
        self._check_tick(timedelta(1), 'D')

    def test_day_corner(self):
        index = _dti(['1/1/2000', '1/2/2000', '1/3/2000'])
        self.assertEqual(infer_freq(index), 'D')

    def test_non_datetimeindex(self):
        dates = to_datetime(['1/1/2000', '1/2/2000', '1/3/2000'])
        self.assertEqual(infer_freq(dates), 'D')

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
        if _np_version_under1p7:
            raise nose.SkipTest("requires numpy >= 1.7 to run")
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
            self.assertEqual(infer_freq(index), exp_freq)

        index = _dti([b + base_delta * 7] +
                     [b + base_delta * j for j in range(3)])
        self.assertIsNone(infer_freq(index))

        index = _dti([b + base_delta * j for j in range(3)] +
                     [b + base_delta * 7])
        self.assertIsNone(infer_freq(index))

    def test_weekly(self):
        days = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']

        for day in days:
            self._check_generated_range('1/1/2000', 'W-%s' % day)

    def test_week_of_month(self):
        days = ['MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN']

        for day in days:
            for i in range(1, 5):
                self._check_generated_range('1/1/2000', 'WOM-%d%s' % (i, day))

    def test_week_of_month_fake(self):
        #All of these dates are on same day of week and are 4 or 5 weeks apart
        index = DatetimeIndex(["2013-08-27","2013-10-01","2013-10-29","2013-11-26"])
        assert infer_freq(index) != 'WOM-4TUE'

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
            self.assertEqual(infer_freq(index), gen.freqstr)
        else:
            inf_freq = infer_freq(index)
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
            self.assertEqual(infer_freq(index), gen.freqstr)
        else:
            inf_freq = infer_freq(index)
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

        # GH 7310
        for tz in [None, 'Asia/Tokyo', 'US/Pacific', 'Europe/Paris']:
            dates = ['2010-11-30', '2010-12-31', '2011-01-31', '2011-02-28']
            idx = DatetimeIndex(dates)
            self.assertEqual(idx.inferred_freq, 'M')

            dates = ['2011-01-01', '2011-01-02', '2011-01-03', '2011-01-04']
            idx = DatetimeIndex(dates)
            self.assertEqual(idx.inferred_freq, 'D')

            dates = ['2011-12-31 22:00', '2011-12-31 23:00', '2012-01-01 00:00', '2012-01-01 01:00']
            idx = DatetimeIndex(dates)
            self.assertEqual(idx.inferred_freq, 'H')

    def test_not_monotonic(self):
        rng = _dti(['1/31/2000', '1/31/2001', '1/31/2002'])
        rng = rng[::-1]
        self.assertIsNone(rng.inferred_freq)

    def test_non_datetimeindex(self):
        rng = _dti(['1/31/2000', '1/31/2001', '1/31/2002'])

        vals = rng.to_pydatetime()

        result = infer_freq(vals)
        self.assertEqual(result, rng.inferred_freq)

    def test_invalid_index_types(self):

        # test all index types
        for i in [ tm.makeIntIndex(10),
                   tm.makeFloatIndex(10),
                   tm.makePeriodIndex(10) ]:
            self.assertRaises(TypeError, lambda : infer_freq(i))

        for i in [ tm.makeStringIndex(10),
                   tm.makeUnicodeIndex(10) ]:
            self.assertRaises(ValueError, lambda : infer_freq(i))

    def test_string_datetimelike_compat(self):

        # GH 6463
        expected = infer_freq(['2004-01', '2004-02', '2004-03', '2004-04'])
        result = infer_freq(Index(['2004-01', '2004-02', '2004-03', '2004-04']))
        self.assertEqual(result,expected)

    def test_series(self):

        # GH6407
        # inferring series

        # invalid type of Series
        for s in [ Series(np.arange(10)),
                   Series(np.arange(10.))]:
            self.assertRaises(TypeError, lambda : infer_freq(s))

        # a non-convertible string
        self.assertRaises(ValueError, lambda : infer_freq(Series(['foo','bar'])))

        # cannot infer on PeriodIndex
        for freq in [None, 'L', 'Y']:
            s = Series(period_range('2013',periods=10,freq=freq))
            self.assertRaises(TypeError, lambda : infer_freq(s))

        # DateTimeIndex
        for freq in ['M', 'L', 'S']:
            s = Series(date_range('20130101',periods=10,freq=freq))
            inferred = infer_freq(s)
            self.assertEqual(inferred,freq)

        s = Series(date_range('20130101','20130110'))
        inferred = infer_freq(s)
        self.assertEqual(inferred,'D')

MONTHS = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP',
          'OCT', 'NOV', 'DEC']


def test_is_superperiod_subperiod():
    assert(fmod.is_superperiod(offsets.YearEnd(), offsets.MonthEnd()))
    assert(fmod.is_subperiod(offsets.MonthEnd(), offsets.YearEnd()))

    assert(fmod.is_superperiod(offsets.Hour(), offsets.Minute()))
    assert(fmod.is_subperiod(offsets.Minute(), offsets.Hour()))

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
