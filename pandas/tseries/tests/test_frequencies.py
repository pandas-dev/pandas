from datetime import datetime, time, timedelta
import sys
import os
import unittest

import nose

import numpy as np

from pandas import Index, DatetimeIndex, date_range

from pandas.tseries.frequencies import to_offset, infer_freq
import pandas.tseries.offsets as offsets

import pandas._tseries as lib

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

    def test_monthly(self):
        self._check_generated_range('1/1/2000', 'M')

    def test_monthly_ambiguous(self):
        rng = _dti(['1/31/2000', '2/29/2000', '3/31/2000'])
        self.assert_(rng.inferred_freq == 'M')

    def test_business_monthly(self):
        self._check_generated_range('1/1/2000', 'BM')

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
        self.assert_(infer_freq(index) == gen.freqstr)

        gen = date_range(start, periods=5, freq=freq)
        index = _dti(gen.values)
        self.assert_(infer_freq(index) == gen.freqstr)

MONTHS = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP',
          'OCT', 'NOV', 'DEC']

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

