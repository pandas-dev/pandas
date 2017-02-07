""" test the scalar Timedelta """
import numpy as np
from datetime import timedelta

import pandas as pd
import pandas.util.testing as tm
from pandas.tseries.timedeltas import _coerce_scalar_to_timedelta_type as ct
from pandas import (Timedelta, TimedeltaIndex, timedelta_range, Series,
                    to_timedelta, tslib, compat, isnull)

iNaT = tslib.iNaT


class TestTimedeltas(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        pass

    def test_construction(self):

        expected = np.timedelta64(10, 'D').astype('m8[ns]').view('i8')
        self.assertEqual(Timedelta(10, unit='d').value, expected)
        self.assertEqual(Timedelta(10.0, unit='d').value, expected)
        self.assertEqual(Timedelta('10 days').value, expected)
        self.assertEqual(Timedelta(days=10).value, expected)
        self.assertEqual(Timedelta(days=10.0).value, expected)

        expected += np.timedelta64(10, 's').astype('m8[ns]').view('i8')
        self.assertEqual(Timedelta('10 days 00:00:10').value, expected)
        self.assertEqual(Timedelta(days=10, seconds=10).value, expected)
        self.assertEqual(
            Timedelta(days=10, milliseconds=10 * 1000).value, expected)
        self.assertEqual(
            Timedelta(days=10, microseconds=10 * 1000 * 1000).value, expected)

        # test construction with np dtypes
        # GH 8757
        timedelta_kwargs = {'days': 'D',
                            'seconds': 's',
                            'microseconds': 'us',
                            'milliseconds': 'ms',
                            'minutes': 'm',
                            'hours': 'h',
                            'weeks': 'W'}
        npdtypes = [np.int64, np.int32, np.int16, np.float64, np.float32,
                    np.float16]
        for npdtype in npdtypes:
            for pykwarg, npkwarg in timedelta_kwargs.items():
                expected = np.timedelta64(1,
                                          npkwarg).astype('m8[ns]').view('i8')
                self.assertEqual(
                    Timedelta(**{pykwarg: npdtype(1)}).value, expected)

        # rounding cases
        self.assertEqual(Timedelta(82739999850000).value, 82739999850000)
        self.assertTrue('0 days 22:58:59.999850' in str(Timedelta(
            82739999850000)))
        self.assertEqual(Timedelta(123072001000000).value, 123072001000000)
        self.assertTrue('1 days 10:11:12.001' in str(Timedelta(
            123072001000000)))

        # string conversion with/without leading zero
        # GH 9570
        self.assertEqual(Timedelta('0:00:00'), timedelta(hours=0))
        self.assertEqual(Timedelta('00:00:00'), timedelta(hours=0))
        self.assertEqual(Timedelta('-1:00:00'), -timedelta(hours=1))
        self.assertEqual(Timedelta('-01:00:00'), -timedelta(hours=1))

        # more strings & abbrevs
        # GH 8190
        self.assertEqual(Timedelta('1 h'), timedelta(hours=1))
        self.assertEqual(Timedelta('1 hour'), timedelta(hours=1))
        self.assertEqual(Timedelta('1 hr'), timedelta(hours=1))
        self.assertEqual(Timedelta('1 hours'), timedelta(hours=1))
        self.assertEqual(Timedelta('-1 hours'), -timedelta(hours=1))
        self.assertEqual(Timedelta('1 m'), timedelta(minutes=1))
        self.assertEqual(Timedelta('1.5 m'), timedelta(seconds=90))
        self.assertEqual(Timedelta('1 minute'), timedelta(minutes=1))
        self.assertEqual(Timedelta('1 minutes'), timedelta(minutes=1))
        self.assertEqual(Timedelta('1 s'), timedelta(seconds=1))
        self.assertEqual(Timedelta('1 second'), timedelta(seconds=1))
        self.assertEqual(Timedelta('1 seconds'), timedelta(seconds=1))
        self.assertEqual(Timedelta('1 ms'), timedelta(milliseconds=1))
        self.assertEqual(Timedelta('1 milli'), timedelta(milliseconds=1))
        self.assertEqual(Timedelta('1 millisecond'), timedelta(milliseconds=1))
        self.assertEqual(Timedelta('1 us'), timedelta(microseconds=1))
        self.assertEqual(Timedelta('1 micros'), timedelta(microseconds=1))
        self.assertEqual(Timedelta('1 microsecond'), timedelta(microseconds=1))
        self.assertEqual(Timedelta('1.5 microsecond'),
                         Timedelta('00:00:00.000001500'))
        self.assertEqual(Timedelta('1 ns'), Timedelta('00:00:00.000000001'))
        self.assertEqual(Timedelta('1 nano'), Timedelta('00:00:00.000000001'))
        self.assertEqual(Timedelta('1 nanosecond'),
                         Timedelta('00:00:00.000000001'))

        # combos
        self.assertEqual(Timedelta('10 days 1 hour'),
                         timedelta(days=10, hours=1))
        self.assertEqual(Timedelta('10 days 1 h'), timedelta(days=10, hours=1))
        self.assertEqual(Timedelta('10 days 1 h 1m 1s'), timedelta(
            days=10, hours=1, minutes=1, seconds=1))
        self.assertEqual(Timedelta('-10 days 1 h 1m 1s'), -
                         timedelta(days=10, hours=1, minutes=1, seconds=1))
        self.assertEqual(Timedelta('-10 days 1 h 1m 1s'), -
                         timedelta(days=10, hours=1, minutes=1, seconds=1))
        self.assertEqual(Timedelta('-10 days 1 h 1m 1s 3us'), -
                         timedelta(days=10, hours=1, minutes=1,
                                   seconds=1, microseconds=3))
        self.assertEqual(Timedelta('-10 days 1 h 1.5m 1s 3us'), -
                         timedelta(days=10, hours=1, minutes=1,
                                   seconds=31, microseconds=3))

        # currently invalid as it has a - on the hhmmdd part (only allowed on
        # the days)
        self.assertRaises(ValueError,
                          lambda: Timedelta('-10 days -1 h 1.5m 1s 3us'))

        # only leading neg signs are allowed
        self.assertRaises(ValueError,
                          lambda: Timedelta('10 days -1 h 1.5m 1s 3us'))

        # no units specified
        self.assertRaises(ValueError, lambda: Timedelta('3.1415'))

        # invalid construction
        tm.assertRaisesRegexp(ValueError, "cannot construct a Timedelta",
                              lambda: Timedelta())
        tm.assertRaisesRegexp(ValueError, "unit abbreviation w/o a number",
                              lambda: Timedelta('foo'))
        tm.assertRaisesRegexp(ValueError,
                              "cannot construct a Timedelta from the passed "
                              "arguments, allowed keywords are ",
                              lambda: Timedelta(day=10))

        # roundtripping both for string and value
        for v in ['1s', '-1s', '1us', '-1us', '1 day', '-1 day',
                  '-23:59:59.999999', '-1 days +23:59:59.999999', '-1ns',
                  '1ns', '-23:59:59.999999999']:

            td = Timedelta(v)
            self.assertEqual(Timedelta(td.value), td)

            # str does not normally display nanos
            if not td.nanoseconds:
                self.assertEqual(Timedelta(str(td)), td)
            self.assertEqual(Timedelta(td._repr_base(format='all')), td)

        # floats
        expected = np.timedelta64(
            10, 's').astype('m8[ns]').view('i8') + np.timedelta64(
                500, 'ms').astype('m8[ns]').view('i8')
        self.assertEqual(Timedelta(10.5, unit='s').value, expected)

        # nat
        self.assertEqual(Timedelta('').value, iNaT)
        self.assertEqual(Timedelta('nat').value, iNaT)
        self.assertEqual(Timedelta('NAT').value, iNaT)
        self.assertEqual(Timedelta(None).value, iNaT)
        self.assertEqual(Timedelta(np.nan).value, iNaT)
        self.assertTrue(isnull(Timedelta('nat')))

        # offset
        self.assertEqual(to_timedelta(pd.offsets.Hour(2)),
                         Timedelta('0 days, 02:00:00'))
        self.assertEqual(Timedelta(pd.offsets.Hour(2)),
                         Timedelta('0 days, 02:00:00'))
        self.assertEqual(Timedelta(pd.offsets.Second(2)),
                         Timedelta('0 days, 00:00:02'))

        # unicode
        # GH 11995
        expected = Timedelta('1H')
        result = pd.Timedelta(u'1H')
        self.assertEqual(result, expected)
        self.assertEqual(to_timedelta(pd.offsets.Hour(2)),
                         Timedelta(u'0 days, 02:00:00'))

        self.assertRaises(ValueError, lambda: Timedelta(u'foo bar'))

    def test_overflow_on_construction(self):
        # xref https://github.com/statsmodels/statsmodels/issues/3374
        value = pd.Timedelta('1day').value * 20169940
        self.assertRaises(OverflowError, pd.Timedelta, value)

    def test_total_seconds_scalar(self):
        # GH 10939
        rng = Timedelta('1 days, 10:11:12.100123456')
        expt = 1 * 86400 + 10 * 3600 + 11 * 60 + 12 + 100123456. / 1e9
        tm.assert_almost_equal(rng.total_seconds(), expt)

        rng = Timedelta(np.nan)
        self.assertTrue(np.isnan(rng.total_seconds()))

    def test_repr(self):

        self.assertEqual(repr(Timedelta(10, unit='d')),
                         "Timedelta('10 days 00:00:00')")
        self.assertEqual(repr(Timedelta(10, unit='s')),
                         "Timedelta('0 days 00:00:10')")
        self.assertEqual(repr(Timedelta(10, unit='ms')),
                         "Timedelta('0 days 00:00:00.010000')")
        self.assertEqual(repr(Timedelta(-10, unit='ms')),
                         "Timedelta('-1 days +23:59:59.990000')")

    def test_conversion(self):

        for td in [Timedelta(10, unit='d'),
                   Timedelta('1 days, 10:11:12.012345')]:
            pydt = td.to_pytimedelta()
            self.assertTrue(td == Timedelta(pydt))
            self.assertEqual(td, pydt)
            self.assertTrue(isinstance(pydt, timedelta) and not isinstance(
                pydt, Timedelta))

            self.assertEqual(td, np.timedelta64(td.value, 'ns'))
            td64 = td.to_timedelta64()
            self.assertEqual(td64, np.timedelta64(td.value, 'ns'))
            self.assertEqual(td, td64)
            self.assertTrue(isinstance(td64, np.timedelta64))

        # this is NOT equal and cannot be roundtriped (because of the nanos)
        td = Timedelta('1 days, 10:11:12.012345678')
        self.assertTrue(td != td.to_pytimedelta())

    def test_freq_conversion(self):

        td = Timedelta('1 days 2 hours 3 ns')
        result = td / np.timedelta64(1, 'D')
        self.assertEqual(result, td.value / float(86400 * 1e9))
        result = td / np.timedelta64(1, 's')
        self.assertEqual(result, td.value / float(1e9))
        result = td / np.timedelta64(1, 'ns')
        self.assertEqual(result, td.value)

    def test_fields(self):
        def check(value):
            # that we are int/long like
            self.assertTrue(isinstance(value, (int, compat.long)))

        # compat to datetime.timedelta
        rng = to_timedelta('1 days, 10:11:12')
        self.assertEqual(rng.days, 1)
        self.assertEqual(rng.seconds, 10 * 3600 + 11 * 60 + 12)
        self.assertEqual(rng.microseconds, 0)
        self.assertEqual(rng.nanoseconds, 0)

        self.assertRaises(AttributeError, lambda: rng.hours)
        self.assertRaises(AttributeError, lambda: rng.minutes)
        self.assertRaises(AttributeError, lambda: rng.milliseconds)

        # GH 10050
        check(rng.days)
        check(rng.seconds)
        check(rng.microseconds)
        check(rng.nanoseconds)

        td = Timedelta('-1 days, 10:11:12')
        self.assertEqual(abs(td), Timedelta('13:48:48'))
        self.assertTrue(str(td) == "-1 days +10:11:12")
        self.assertEqual(-td, Timedelta('0 days 13:48:48'))
        self.assertEqual(-Timedelta('-1 days, 10:11:12').value, 49728000000000)
        self.assertEqual(Timedelta('-1 days, 10:11:12').value, -49728000000000)

        rng = to_timedelta('-1 days, 10:11:12.100123456')
        self.assertEqual(rng.days, -1)
        self.assertEqual(rng.seconds, 10 * 3600 + 11 * 60 + 12)
        self.assertEqual(rng.microseconds, 100 * 1000 + 123)
        self.assertEqual(rng.nanoseconds, 456)
        self.assertRaises(AttributeError, lambda: rng.hours)
        self.assertRaises(AttributeError, lambda: rng.minutes)
        self.assertRaises(AttributeError, lambda: rng.milliseconds)

        # components
        tup = pd.to_timedelta(-1, 'us').components
        self.assertEqual(tup.days, -1)
        self.assertEqual(tup.hours, 23)
        self.assertEqual(tup.minutes, 59)
        self.assertEqual(tup.seconds, 59)
        self.assertEqual(tup.milliseconds, 999)
        self.assertEqual(tup.microseconds, 999)
        self.assertEqual(tup.nanoseconds, 0)

        # GH 10050
        check(tup.days)
        check(tup.hours)
        check(tup.minutes)
        check(tup.seconds)
        check(tup.milliseconds)
        check(tup.microseconds)
        check(tup.nanoseconds)

        tup = Timedelta('-1 days 1 us').components
        self.assertEqual(tup.days, -2)
        self.assertEqual(tup.hours, 23)
        self.assertEqual(tup.minutes, 59)
        self.assertEqual(tup.seconds, 59)
        self.assertEqual(tup.milliseconds, 999)
        self.assertEqual(tup.microseconds, 999)
        self.assertEqual(tup.nanoseconds, 0)

    def test_nat_converters(self):
        self.assertEqual(to_timedelta(
            'nat', box=False).astype('int64'), tslib.iNaT)
        self.assertEqual(to_timedelta(
            'nan', box=False).astype('int64'), tslib.iNaT)

        def testit(unit, transform):

            # array
            result = to_timedelta(np.arange(5), unit=unit)
            expected = TimedeltaIndex([np.timedelta64(i, transform(unit))
                                       for i in np.arange(5).tolist()])
            tm.assert_index_equal(result, expected)

            # scalar
            result = to_timedelta(2, unit=unit)
            expected = Timedelta(np.timedelta64(2, transform(unit)).astype(
                'timedelta64[ns]'))
            self.assertEqual(result, expected)

        # validate all units
        # GH 6855
        for unit in ['Y', 'M', 'W', 'D', 'y', 'w', 'd']:
            testit(unit, lambda x: x.upper())
        for unit in ['days', 'day', 'Day', 'Days']:
            testit(unit, lambda x: 'D')
        for unit in ['h', 'm', 's', 'ms', 'us', 'ns', 'H', 'S', 'MS', 'US',
                     'NS']:
            testit(unit, lambda x: x.lower())

        # offsets

        # m
        testit('T', lambda x: 'm')

        # ms
        testit('L', lambda x: 'ms')

    def test_numeric_conversions(self):
        self.assertEqual(ct(0), np.timedelta64(0, 'ns'))
        self.assertEqual(ct(10), np.timedelta64(10, 'ns'))
        self.assertEqual(ct(10, unit='ns'), np.timedelta64(
            10, 'ns').astype('m8[ns]'))

        self.assertEqual(ct(10, unit='us'), np.timedelta64(
            10, 'us').astype('m8[ns]'))
        self.assertEqual(ct(10, unit='ms'), np.timedelta64(
            10, 'ms').astype('m8[ns]'))
        self.assertEqual(ct(10, unit='s'), np.timedelta64(
            10, 's').astype('m8[ns]'))
        self.assertEqual(ct(10, unit='d'), np.timedelta64(
            10, 'D').astype('m8[ns]'))

    def test_timedelta_conversions(self):
        self.assertEqual(ct(timedelta(seconds=1)),
                         np.timedelta64(1, 's').astype('m8[ns]'))
        self.assertEqual(ct(timedelta(microseconds=1)),
                         np.timedelta64(1, 'us').astype('m8[ns]'))
        self.assertEqual(ct(timedelta(days=1)),
                         np.timedelta64(1, 'D').astype('m8[ns]'))

    def test_round(self):

        t1 = Timedelta('1 days 02:34:56.789123456')
        t2 = Timedelta('-1 days 02:34:56.789123456')

        for (freq, s1, s2) in [('N', t1, t2),
                               ('U', Timedelta('1 days 02:34:56.789123000'),
                                Timedelta('-1 days 02:34:56.789123000')),
                               ('L', Timedelta('1 days 02:34:56.789000000'),
                                Timedelta('-1 days 02:34:56.789000000')),
                               ('S', Timedelta('1 days 02:34:57'),
                                Timedelta('-1 days 02:34:57')),
                               ('2S', Timedelta('1 days 02:34:56'),
                                Timedelta('-1 days 02:34:56')),
                               ('5S', Timedelta('1 days 02:34:55'),
                                Timedelta('-1 days 02:34:55')),
                               ('T', Timedelta('1 days 02:35:00'),
                                Timedelta('-1 days 02:35:00')),
                               ('12T', Timedelta('1 days 02:36:00'),
                                Timedelta('-1 days 02:36:00')),
                               ('H', Timedelta('1 days 03:00:00'),
                                Timedelta('-1 days 03:00:00')),
                               ('d', Timedelta('1 days'),
                                Timedelta('-1 days'))]:
            r1 = t1.round(freq)
            self.assertEqual(r1, s1)
            r2 = t2.round(freq)
            self.assertEqual(r2, s2)

        # invalid
        for freq in ['Y', 'M', 'foobar']:
            self.assertRaises(ValueError, lambda: t1.round(freq))

        t1 = timedelta_range('1 days', periods=3, freq='1 min 2 s 3 us')
        t2 = -1 * t1
        t1a = timedelta_range('1 days', periods=3, freq='1 min 2 s')
        t1c = pd.TimedeltaIndex([1, 1, 1], unit='D')

        # note that negative times round DOWN! so don't give whole numbers
        for (freq, s1, s2) in [('N', t1, t2),
                               ('U', t1, t2),
                               ('L', t1a,
                                TimedeltaIndex(['-1 days +00:00:00',
                                                '-2 days +23:58:58',
                                                '-2 days +23:57:56'],
                                               dtype='timedelta64[ns]',
                                               freq=None)
                                ),
                               ('S', t1a,
                                TimedeltaIndex(['-1 days +00:00:00',
                                                '-2 days +23:58:58',
                                                '-2 days +23:57:56'],
                                               dtype='timedelta64[ns]',
                                               freq=None)
                                ),
                               ('12T', t1c,
                                TimedeltaIndex(['-1 days',
                                                '-1 days',
                                                '-1 days'],
                                               dtype='timedelta64[ns]',
                                               freq=None)
                                ),
                               ('H', t1c,
                                TimedeltaIndex(['-1 days',
                                                '-1 days',
                                                '-1 days'],
                                               dtype='timedelta64[ns]',
                                               freq=None)
                                ),
                               ('d', t1c,
                                pd.TimedeltaIndex([-1, -1, -1], unit='D')
                                )]:

            r1 = t1.round(freq)
            tm.assert_index_equal(r1, s1)
            r2 = t2.round(freq)
        tm.assert_index_equal(r2, s2)

        # invalid
        for freq in ['Y', 'M', 'foobar']:
            self.assertRaises(ValueError, lambda: t1.round(freq))

    def test_contains(self):
        # Checking for any NaT-like objects
        # GH 13603
        td = to_timedelta(range(5), unit='d') + pd.offsets.Hour(1)
        for v in [pd.NaT, None, float('nan'), np.nan]:
            self.assertFalse((v in td))

        td = to_timedelta([pd.NaT])
        for v in [pd.NaT, None, float('nan'), np.nan]:
            self.assertTrue((v in td))

    def test_identity(self):

        td = Timedelta(10, unit='d')
        self.assertTrue(isinstance(td, Timedelta))
        self.assertTrue(isinstance(td, timedelta))

    def test_short_format_converters(self):
        def conv(v):
            return v.astype('m8[ns]')

        self.assertEqual(ct('10'), np.timedelta64(10, 'ns'))
        self.assertEqual(ct('10ns'), np.timedelta64(10, 'ns'))
        self.assertEqual(ct('100'), np.timedelta64(100, 'ns'))
        self.assertEqual(ct('100ns'), np.timedelta64(100, 'ns'))

        self.assertEqual(ct('1000'), np.timedelta64(1000, 'ns'))
        self.assertEqual(ct('1000ns'), np.timedelta64(1000, 'ns'))
        self.assertEqual(ct('1000NS'), np.timedelta64(1000, 'ns'))

        self.assertEqual(ct('10us'), np.timedelta64(10000, 'ns'))
        self.assertEqual(ct('100us'), np.timedelta64(100000, 'ns'))
        self.assertEqual(ct('1000us'), np.timedelta64(1000000, 'ns'))
        self.assertEqual(ct('1000Us'), np.timedelta64(1000000, 'ns'))
        self.assertEqual(ct('1000uS'), np.timedelta64(1000000, 'ns'))

        self.assertEqual(ct('1ms'), np.timedelta64(1000000, 'ns'))
        self.assertEqual(ct('10ms'), np.timedelta64(10000000, 'ns'))
        self.assertEqual(ct('100ms'), np.timedelta64(100000000, 'ns'))
        self.assertEqual(ct('1000ms'), np.timedelta64(1000000000, 'ns'))

        self.assertEqual(ct('-1s'), -np.timedelta64(1000000000, 'ns'))
        self.assertEqual(ct('1s'), np.timedelta64(1000000000, 'ns'))
        self.assertEqual(ct('10s'), np.timedelta64(10000000000, 'ns'))
        self.assertEqual(ct('100s'), np.timedelta64(100000000000, 'ns'))
        self.assertEqual(ct('1000s'), np.timedelta64(1000000000000, 'ns'))

        self.assertEqual(ct('1d'), conv(np.timedelta64(1, 'D')))
        self.assertEqual(ct('-1d'), -conv(np.timedelta64(1, 'D')))
        self.assertEqual(ct('1D'), conv(np.timedelta64(1, 'D')))
        self.assertEqual(ct('10D'), conv(np.timedelta64(10, 'D')))
        self.assertEqual(ct('100D'), conv(np.timedelta64(100, 'D')))
        self.assertEqual(ct('1000D'), conv(np.timedelta64(1000, 'D')))
        self.assertEqual(ct('10000D'), conv(np.timedelta64(10000, 'D')))

        # space
        self.assertEqual(ct(' 10000D '), conv(np.timedelta64(10000, 'D')))
        self.assertEqual(ct(' - 10000D '), -conv(np.timedelta64(10000, 'D')))

        # invalid
        self.assertRaises(ValueError, ct, '1foo')
        self.assertRaises(ValueError, ct, 'foo')

    def test_full_format_converters(self):
        def conv(v):
            return v.astype('m8[ns]')

        d1 = np.timedelta64(1, 'D')

        self.assertEqual(ct('1days'), conv(d1))
        self.assertEqual(ct('1days,'), conv(d1))
        self.assertEqual(ct('- 1days,'), -conv(d1))

        self.assertEqual(ct('00:00:01'), conv(np.timedelta64(1, 's')))
        self.assertEqual(ct('06:00:01'), conv(
            np.timedelta64(6 * 3600 + 1, 's')))
        self.assertEqual(ct('06:00:01.0'), conv(
            np.timedelta64(6 * 3600 + 1, 's')))
        self.assertEqual(ct('06:00:01.01'), conv(
            np.timedelta64(1000 * (6 * 3600 + 1) + 10, 'ms')))

        self.assertEqual(ct('- 1days, 00:00:01'),
                         conv(-d1 + np.timedelta64(1, 's')))
        self.assertEqual(ct('1days, 06:00:01'), conv(
            d1 + np.timedelta64(6 * 3600 + 1, 's')))
        self.assertEqual(ct('1days, 06:00:01.01'), conv(
            d1 + np.timedelta64(1000 * (6 * 3600 + 1) + 10, 'ms')))

        # invalid
        self.assertRaises(ValueError, ct, '- 1days, 00')

    def test_overflow(self):
        # GH 9442
        s = Series(pd.date_range('20130101', periods=100000, freq='H'))
        s[0] += pd.Timedelta('1s 1ms')

        # mean
        result = (s - s.min()).mean()
        expected = pd.Timedelta((pd.DatetimeIndex((s - s.min())).asi8 / len(s)
                                 ).sum())

        # the computation is converted to float so might be some loss of
        # precision
        self.assertTrue(np.allclose(result.value / 1000, expected.value /
                                    1000))

        # sum
        self.assertRaises(ValueError, lambda: (s - s.min()).sum())
        s1 = s[0:10000]
        self.assertRaises(ValueError, lambda: (s1 - s1.min()).sum())
        s2 = s[0:1000]
        result = (s2 - s2.min()).sum()

    def test_pickle(self):

        v = Timedelta('1 days 10:11:12.0123456')
        v_p = self.round_trip_pickle(v)
        self.assertEqual(v, v_p)

    def test_timedelta_hash_equality(self):
        # GH 11129
        v = Timedelta(1, 'D')
        td = timedelta(days=1)
        self.assertEqual(hash(v), hash(td))

        d = {td: 2}
        self.assertEqual(d[v], 2)

        tds = timedelta_range('1 second', periods=20)
        self.assertTrue(all(hash(td) == hash(td.to_pytimedelta()) for td in
                            tds))

        # python timedeltas drop ns resolution
        ns_td = Timedelta(1, 'ns')
        self.assertNotEqual(hash(ns_td), hash(ns_td.to_pytimedelta()))

    def test_implementation_limits(self):
        min_td = Timedelta(Timedelta.min)
        max_td = Timedelta(Timedelta.max)

        # GH 12727
        # timedelta limits correspond to int64 boundaries
        self.assertTrue(min_td.value == np.iinfo(np.int64).min + 1)
        self.assertTrue(max_td.value == np.iinfo(np.int64).max)

        # Beyond lower limit, a NAT before the Overflow
        self.assertIsInstance(min_td - Timedelta(1, 'ns'),
                              pd.tslib.NaTType)

        with tm.assertRaises(OverflowError):
            min_td - Timedelta(2, 'ns')

        with tm.assertRaises(OverflowError):
            max_td + Timedelta(1, 'ns')

        # Same tests using the internal nanosecond values
        td = Timedelta(min_td.value - 1, 'ns')
        self.assertIsInstance(td, pd.tslib.NaTType)

        with tm.assertRaises(OverflowError):
            Timedelta(min_td.value - 2, 'ns')

        with tm.assertRaises(OverflowError):
            Timedelta(max_td.value + 1, 'ns')

    def test_timedelta_arithmetic(self):
        data = pd.Series(['nat', '32 days'], dtype='timedelta64[ns]')
        deltas = [timedelta(days=1), Timedelta(1, unit='D')]
        for delta in deltas:
            result_method = data.add(delta)
            result_operator = data + delta
            expected = pd.Series(['nat', '33 days'], dtype='timedelta64[ns]')
            tm.assert_series_equal(result_operator, expected)
            tm.assert_series_equal(result_method, expected)

            result_method = data.sub(delta)
            result_operator = data - delta
            expected = pd.Series(['nat', '31 days'], dtype='timedelta64[ns]')
            tm.assert_series_equal(result_operator, expected)
            tm.assert_series_equal(result_method, expected)
            # GH 9396
            result_method = data.div(delta)
            result_operator = data / delta
            expected = pd.Series([np.nan, 32.], dtype='float64')
            tm.assert_series_equal(result_operator, expected)
            tm.assert_series_equal(result_method, expected)

    def test_apply_to_timedelta(self):
        timedelta_NaT = pd.to_timedelta('NaT')

        list_of_valid_strings = ['00:00:01', '00:00:02']
        a = pd.to_timedelta(list_of_valid_strings)
        b = Series(list_of_valid_strings).apply(pd.to_timedelta)
        # Can't compare until apply on a Series gives the correct dtype
        # assert_series_equal(a, b)

        list_of_strings = ['00:00:01', np.nan, pd.NaT, timedelta_NaT]

        # TODO: unused?
        a = pd.to_timedelta(list_of_strings)  # noqa
        b = Series(list_of_strings).apply(pd.to_timedelta)  # noqa
        # Can't compare until apply on a Series gives the correct dtype
        # assert_series_equal(a, b)

    def test_components(self):
        rng = timedelta_range('1 days, 10:11:12', periods=2, freq='s')
        rng.components

        # with nat
        s = Series(rng)
        s[1] = np.nan

        result = s.dt.components
        self.assertFalse(result.iloc[0].isnull().all())
        self.assertTrue(result.iloc[1].isnull().all())

    def test_isoformat(self):
        td = Timedelta(days=6, minutes=50, seconds=3,
                       milliseconds=10, microseconds=10, nanoseconds=12)
        expected = 'P6DT0H50M3.010010012S'
        result = td.isoformat()
        self.assertEqual(result, expected)

        td = Timedelta(days=4, hours=12, minutes=30, seconds=5)
        result = td.isoformat()
        expected = 'P4DT12H30M5S'
        self.assertEqual(result, expected)

        td = Timedelta(nanoseconds=123)
        result = td.isoformat()
        expected = 'P0DT0H0M0.000000123S'
        self.assertEqual(result, expected)

        # trim nano
        td = Timedelta(microseconds=10)
        result = td.isoformat()
        expected = 'P0DT0H0M0.00001S'
        self.assertEqual(result, expected)

        # trim micro
        td = Timedelta(milliseconds=1)
        result = td.isoformat()
        expected = 'P0DT0H0M0.001S'
        self.assertEqual(result, expected)

        # NaT
        result = Timedelta('NaT').isoformat()
        expected = 'NaT'
        self.assertEqual(result, expected)

        # don't strip every 0
        result = Timedelta(minutes=1).isoformat()
        expected = 'P0DT0H1M0S'
        self.assertEqual(result, expected)

    def test_ops_error_str(self):
        # GH 13624
        td = Timedelta('1 day')

        for l, r in [(td, 'a'), ('a', td)]:

            with tm.assertRaises(TypeError):
                l + r

            with tm.assertRaises(TypeError):
                l > r

            self.assertFalse(l == r)
            self.assertTrue(l != r)
