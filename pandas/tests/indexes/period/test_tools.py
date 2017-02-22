import numpy as np
from datetime import datetime, timedelta

import pandas as pd
import pandas.util.testing as tm
import pandas.tseries.period as period
from pandas.compat import lrange
from pandas.tseries.frequencies import get_freq, MONTHS
from pandas._period import period_ordinal, period_asfreq
from pandas import (PeriodIndex, Period, DatetimeIndex, Timestamp, Series,
                    date_range, to_datetime, period_range)


class TestPeriodRepresentation(tm.TestCase):
    """
    Wish to match NumPy units
    """

    def _check_freq(self, freq, base_date):
        rng = PeriodIndex(start=base_date, periods=10, freq=freq)
        exp = np.arange(10, dtype=np.int64)
        self.assert_numpy_array_equal(rng._values, exp)
        self.assert_numpy_array_equal(rng.asi8, exp)

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


class TestTslib(tm.TestCase):
    def test_intraday_conversion_factors(self):
        self.assertEqual(period_asfreq(
            1, get_freq('D'), get_freq('H'), False), 24)
        self.assertEqual(period_asfreq(
            1, get_freq('D'), get_freq('T'), False), 1440)
        self.assertEqual(period_asfreq(
            1, get_freq('D'), get_freq('S'), False), 86400)
        self.assertEqual(period_asfreq(1, get_freq(
            'D'), get_freq('L'), False), 86400000)
        self.assertEqual(period_asfreq(1, get_freq(
            'D'), get_freq('U'), False), 86400000000)
        self.assertEqual(period_asfreq(1, get_freq(
            'D'), get_freq('N'), False), 86400000000000)

        self.assertEqual(period_asfreq(
            1, get_freq('H'), get_freq('T'), False), 60)
        self.assertEqual(period_asfreq(
            1, get_freq('H'), get_freq('S'), False), 3600)
        self.assertEqual(period_asfreq(1, get_freq('H'),
                                       get_freq('L'), False), 3600000)
        self.assertEqual(period_asfreq(1, get_freq(
            'H'), get_freq('U'), False), 3600000000)
        self.assertEqual(period_asfreq(1, get_freq(
            'H'), get_freq('N'), False), 3600000000000)

        self.assertEqual(period_asfreq(
            1, get_freq('T'), get_freq('S'), False), 60)
        self.assertEqual(period_asfreq(
            1, get_freq('T'), get_freq('L'), False), 60000)
        self.assertEqual(period_asfreq(1, get_freq(
            'T'), get_freq('U'), False), 60000000)
        self.assertEqual(period_asfreq(1, get_freq(
            'T'), get_freq('N'), False), 60000000000)

        self.assertEqual(period_asfreq(
            1, get_freq('S'), get_freq('L'), False), 1000)
        self.assertEqual(period_asfreq(1, get_freq('S'),
                                       get_freq('U'), False), 1000000)
        self.assertEqual(period_asfreq(1, get_freq(
            'S'), get_freq('N'), False), 1000000000)

        self.assertEqual(period_asfreq(
            1, get_freq('L'), get_freq('U'), False), 1000)
        self.assertEqual(period_asfreq(1, get_freq('L'),
                                       get_freq('N'), False), 1000000)

        self.assertEqual(period_asfreq(
            1, get_freq('U'), get_freq('N'), False), 1000)

    def test_period_ordinal_start_values(self):
        # information for 1.1.1970
        self.assertEqual(0, period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0,
                                           get_freq('A')))
        self.assertEqual(0, period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0,
                                           get_freq('M')))
        self.assertEqual(1, period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0,
                                           get_freq('W')))
        self.assertEqual(0, period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0,
                                           get_freq('D')))
        self.assertEqual(0, period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0,
                                           get_freq('B')))

    def test_period_ordinal_week(self):
        self.assertEqual(1, period_ordinal(1970, 1, 4, 0, 0, 0, 0, 0,
                                           get_freq('W')))
        self.assertEqual(2, period_ordinal(1970, 1, 5, 0, 0, 0, 0, 0,
                                           get_freq('W')))

        self.assertEqual(2284, period_ordinal(2013, 10, 6, 0, 0, 0, 0, 0,
                                              get_freq('W')))
        self.assertEqual(2285, period_ordinal(2013, 10, 7, 0, 0, 0, 0, 0,
                                              get_freq('W')))

    def test_period_ordinal_business_day(self):
        # Thursday
        self.assertEqual(11415, period_ordinal(2013, 10, 3, 0, 0, 0, 0, 0,
                                               get_freq('B')))
        # Friday
        self.assertEqual(11416, period_ordinal(2013, 10, 4, 0, 0, 0, 0, 0,
                                               get_freq('B')))
        # Saturday
        self.assertEqual(11417, period_ordinal(2013, 10, 5, 0, 0, 0, 0, 0,
                                               get_freq('B')))
        # Sunday
        self.assertEqual(11417, period_ordinal(2013, 10, 6, 0, 0, 0, 0, 0,
                                               get_freq('B')))
        # Monday
        self.assertEqual(11417, period_ordinal(2013, 10, 7, 0, 0, 0, 0, 0,
                                               get_freq('B')))
        # Tuesday
        self.assertEqual(11418, period_ordinal(2013, 10, 8, 0, 0, 0, 0, 0,
                                               get_freq('B')))


class TestPeriodIndex(tm.TestCase):

    def setUp(self):
        pass

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

    def test_period_dt64_round_trip(self):
        dti = date_range('1/1/2000', '1/7/2002', freq='B')
        pi = dti.to_period()
        tm.assert_index_equal(pi.to_timestamp(), dti)

        dti = date_range('1/1/2000', '1/7/2002', freq='B')
        pi = dti.to_period(freq='H')
        tm.assert_index_equal(pi.to_timestamp(), dti)

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

            with tm.assert_produces_warning(FutureWarning):
                pidx.searchsorted(key=p2)
