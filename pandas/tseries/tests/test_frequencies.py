from datetime import datetime, time, timedelta
import sys
import os
import unittest

import nose

import numpy as np

from pandas import Index, DatetimeIndex, date_range, period_range

from pandas.tseries.frequencies import to_offset, infer_freq
from pandas.tseries.tools import to_datetime
import pandas.tseries.frequencies as fmod
import pandas.tseries.offsets as offsets

import pandas.lib as lib


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


def test_anchored_shortcuts():
    result = to_offset('W')
    expected = to_offset('W-SUN')
    assert(result == expected)

    result = to_offset('Q')
    expected = to_offset('Q-DEC')
    assert(result == expected)


_dti = DatetimeIndex


class TestFrequencyInference(unittest.TestCase):

    def test_raise_if_too_few(self):
        index = _dti(['12/31/1998', '1/3/1999'])
        self.assertRaises(ValueError, infer_freq, index)

    def test_business_daily(self):
        index = _dti(['12/31/1998', '1/3/1999', '1/4/1999'])
        self.assert_(infer_freq(index) == 'B')

    def test_day(self):
        self._check_tick(timedelta(1), 'D')

    def test_day_corner(self):
        index = _dti(['1/1/2000', '1/2/2000', '1/3/2000'])
        self.assert_(infer_freq(index) == 'D')

    def test_non_datetimeindex(self):
        dates = to_datetime(['1/1/2000', '1/2/2000', '1/3/2000'])
        self.assert_(infer_freq(dates) == 'D')

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
        idx = DatetimeIndex(np.arange(0, 100, 10))
        inferred = idx.inferred_freq

        self.assert_(inferred == '10N')

    def _check_tick(self, base_delta, code):
        b = datetime.now()
        for i in range(1, 5):
            inc = base_delta * i
            index = _dti([b + inc * j for j in range(3)])
            if i > 1:
                exp_freq = '%d%s' % (i, code)
            else:
                exp_freq = code
            self.assert_(infer_freq(index) == exp_freq)

        index = _dti([b + base_delta * 7] +
                     [b + base_delta * j for j in range(3)])
        self.assert_(infer_freq(index) is None)

        index = _dti([b + base_delta * j for j in range(3)] +
                     [b + base_delta * 7])
        self.assert_(infer_freq(index) is None)

    def test_weekly(self):
        days = ['MON', 'TUE', 'WED', 'THU', 'FRI']

        for day in days:
            self._check_generated_range('1/1/2000', 'W-%s' % day)

    def test_week_of_month(self):
        days = ['MON', 'TUE', 'WED', 'THU', 'FRI']

        for day in days:
            for i in range(1, 5):
                self._check_generated_range('1/1/2000', 'WOM-%d%s' % (i, day))

    def test_monthly(self):
        self._check_generated_range('1/1/2000', 'M')

    def test_monthly_ambiguous(self):
        rng = _dti(['1/31/2000', '2/29/2000', '3/31/2000'])
        self.assert_(rng.inferred_freq == 'M')

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
        self.assert_(rng.inferred_freq == 'A-JAN')

    def _check_generated_range(self, start, freq):
        freq = freq.upper()

        gen = date_range(start, periods=7, freq=freq)
        index = _dti(gen.values)
        if not freq.startswith('Q-'):
            self.assert_(infer_freq(index) == gen.freqstr)
        else:
            inf_freq = infer_freq(index)
            self.assert_((inf_freq == 'Q-DEC' and
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
            self.assert_(infer_freq(index) == gen.freqstr)
        else:
            inf_freq = infer_freq(index)
            self.assert_((inf_freq == 'Q-DEC' and
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
        self.assert_(rng.inferred_freq == 'Q-DEC')

        rng = period_range('1959Q2', '2009Q3', freq='Q-NOV')
        rng = Index(rng.to_timestamp('D', how='e').asobject)
        self.assert_(rng.inferred_freq == 'Q-NOV')

        rng = period_range('1959Q2', '2009Q3', freq='Q-OCT')
        rng = Index(rng.to_timestamp('D', how='e').asobject)
        self.assert_(rng.inferred_freq == 'Q-OCT')

    def test_not_monotonic(self):
        rng = _dti(['1/31/2000', '1/31/2001', '1/31/2002'])
        rng = rng[::-1]
        self.assert_(rng.inferred_freq is None)

    def test_non_datetimeindex(self):
        rng = _dti(['1/31/2000', '1/31/2001', '1/31/2002'])

        vals = rng.to_pydatetime()

        result = infer_freq(vals)
        self.assertEqual(result, rng.inferred_freq)

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
