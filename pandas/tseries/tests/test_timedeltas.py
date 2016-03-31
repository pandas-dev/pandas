# pylint: disable-msg=E1101,W0612

from __future__ import division
from datetime import timedelta, time
import nose

from distutils.version import LooseVersion
import numpy as np
import pandas as pd

from pandas import (Index, Series, DataFrame, Timestamp, Timedelta,
                    TimedeltaIndex, isnull, date_range,
                    timedelta_range, Int64Index)
from pandas.compat import range
from pandas import compat, to_timedelta, tslib
from pandas.tseries.timedeltas import _coerce_scalar_to_timedelta_type as ct
from pandas.util.testing import (assert_series_equal, assert_frame_equal,
                                 assert_almost_equal, assert_index_equal)
from numpy.testing import assert_allclose
from pandas.tseries.offsets import Day, Second
import pandas.util.testing as tm
from numpy.random import randn
from pandas import _np_version_under1p8

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
        tm.assertRaisesRegexp(ValueError, "cannot construct a TimeDelta",
                              lambda: Timedelta())
        tm.assertRaisesRegexp(ValueError, "unit abbreviation w/o a number",
                              lambda: Timedelta('foo'))
        tm.assertRaisesRegexp(ValueError,
                              "cannot construct a TimeDelta from the passed "
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
        self.assertTrue(isnull(Timestamp('nat')))
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

    def test_repr(self):

        self.assertEqual(repr(Timedelta(10, unit='d')),
                         "Timedelta('10 days 00:00:00')")
        self.assertEqual(repr(Timedelta(10, unit='s')),
                         "Timedelta('0 days 00:00:10')")
        self.assertEqual(repr(Timedelta(10, unit='ms')),
                         "Timedelta('0 days 00:00:00.010000')")
        self.assertEqual(repr(Timedelta(-10, unit='ms')),
                         "Timedelta('-1 days +23:59:59.990000')")

    def test_identity(self):

        td = Timedelta(10, unit='d')
        self.assertTrue(isinstance(td, Timedelta))
        self.assertTrue(isinstance(td, timedelta))

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

    def test_ops(self):

        td = Timedelta(10, unit='d')
        self.assertEqual(-td, Timedelta(-10, unit='d'))
        self.assertEqual(+td, Timedelta(10, unit='d'))
        self.assertEqual(td - td, Timedelta(0, unit='ns'))
        self.assertTrue((td - pd.NaT) is pd.NaT)
        self.assertEqual(td + td, Timedelta(20, unit='d'))
        self.assertTrue((td + pd.NaT) is pd.NaT)
        self.assertEqual(td * 2, Timedelta(20, unit='d'))
        self.assertTrue((td * pd.NaT) is pd.NaT)
        self.assertEqual(td / 2, Timedelta(5, unit='d'))
        self.assertEqual(abs(td), td)
        self.assertEqual(abs(-td), td)
        self.assertEqual(td / td, 1)
        self.assertTrue((td / pd.NaT) is np.nan)

        # invert
        self.assertEqual(-td, Timedelta('-10d'))
        self.assertEqual(td * -1, Timedelta('-10d'))
        self.assertEqual(-1 * td, Timedelta('-10d'))
        self.assertEqual(abs(-td), Timedelta('10d'))

        # invalid
        self.assertRaises(TypeError, lambda: Timedelta(11, unit='d') // 2)

        # invalid multiply with another timedelta
        self.assertRaises(TypeError, lambda: td * td)

        # can't operate with integers
        self.assertRaises(TypeError, lambda: td + 2)
        self.assertRaises(TypeError, lambda: td - 2)

    def test_ops_offsets(self):
        td = Timedelta(10, unit='d')
        self.assertEqual(Timedelta(241, unit='h'), td + pd.offsets.Hour(1))
        self.assertEqual(Timedelta(241, unit='h'), pd.offsets.Hour(1) + td)
        self.assertEqual(240, td / pd.offsets.Hour(1))
        self.assertEqual(1 / 240.0, pd.offsets.Hour(1) / td)
        self.assertEqual(Timedelta(239, unit='h'), td - pd.offsets.Hour(1))
        self.assertEqual(Timedelta(-239, unit='h'), pd.offsets.Hour(1) - td)

    def test_freq_conversion(self):

        td = Timedelta('1 days 2 hours 3 ns')
        result = td / np.timedelta64(1, 'D')
        self.assertEqual(result, td.value / float(86400 * 1e9))
        result = td / np.timedelta64(1, 's')
        self.assertEqual(result, td.value / float(1e9))
        result = td / np.timedelta64(1, 'ns')
        self.assertEqual(result, td.value)

    def test_ops_ndarray(self):
        td = Timedelta('1 day')

        # timedelta, timedelta
        other = pd.to_timedelta(['1 day']).values
        expected = pd.to_timedelta(['2 days']).values
        self.assert_numpy_array_equal(td + other, expected)
        if LooseVersion(np.__version__) >= '1.8':
            self.assert_numpy_array_equal(other + td, expected)
        self.assertRaises(TypeError, lambda: td + np.array([1]))
        self.assertRaises(TypeError, lambda: np.array([1]) + td)

        expected = pd.to_timedelta(['0 days']).values
        self.assert_numpy_array_equal(td - other, expected)
        if LooseVersion(np.__version__) >= '1.8':
            self.assert_numpy_array_equal(-other + td, expected)
        self.assertRaises(TypeError, lambda: td - np.array([1]))
        self.assertRaises(TypeError, lambda: np.array([1]) - td)

        expected = pd.to_timedelta(['2 days']).values
        self.assert_numpy_array_equal(td * np.array([2]), expected)
        self.assert_numpy_array_equal(np.array([2]) * td, expected)
        self.assertRaises(TypeError, lambda: td * other)
        self.assertRaises(TypeError, lambda: other * td)

        self.assert_numpy_array_equal(td / other, np.array([1]))
        if LooseVersion(np.__version__) >= '1.8':
            self.assert_numpy_array_equal(other / td, np.array([1]))

        # timedelta, datetime
        other = pd.to_datetime(['2000-01-01']).values
        expected = pd.to_datetime(['2000-01-02']).values
        self.assert_numpy_array_equal(td + other, expected)
        if LooseVersion(np.__version__) >= '1.8':
            self.assert_numpy_array_equal(other + td, expected)

        expected = pd.to_datetime(['1999-12-31']).values
        self.assert_numpy_array_equal(-td + other, expected)
        if LooseVersion(np.__version__) >= '1.8':
            self.assert_numpy_array_equal(other - td, expected)

    def test_ops_series(self):
        # regression test for GH8813
        td = Timedelta('1 day')
        other = pd.Series([1, 2])
        expected = pd.Series(pd.to_timedelta(['1 day', '2 days']))
        tm.assert_series_equal(expected, td * other)
        tm.assert_series_equal(expected, other * td)

    def test_compare_timedelta_series(self):
        # regresssion test for GH5963
        s = pd.Series([timedelta(days=1), timedelta(days=2)])
        actual = s > timedelta(days=1)
        expected = pd.Series([False, True])
        tm.assert_series_equal(actual, expected)

    def test_compare_timedelta_ndarray(self):
        # GH11835
        periods = [Timedelta('0 days 01:00:00'), Timedelta('0 days 01:00:00')]
        arr = np.array(periods)
        result = arr[0] > arr
        expected = np.array([False, False])
        self.assert_numpy_array_equal(result, expected)

    def test_ops_notimplemented(self):
        class Other:
            pass

        other = Other()

        td = Timedelta('1 day')
        self.assertTrue(td.__add__(other) is NotImplemented)
        self.assertTrue(td.__sub__(other) is NotImplemented)
        self.assertTrue(td.__truediv__(other) is NotImplemented)
        self.assertTrue(td.__mul__(other) is NotImplemented)
        self.assertTrue(td.__floordiv__(td) is NotImplemented)

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

    def test_timedelta_range(self):

        expected = to_timedelta(np.arange(5), unit='D')
        result = timedelta_range('0 days', periods=5, freq='D')
        tm.assert_index_equal(result, expected)

        expected = to_timedelta(np.arange(11), unit='D')
        result = timedelta_range('0 days', '10 days', freq='D')
        tm.assert_index_equal(result, expected)

        expected = to_timedelta(np.arange(5), unit='D') + Second(2) + Day()
        result = timedelta_range('1 days, 00:00:02', '5 days, 00:00:02',
                                 freq='D')
        tm.assert_index_equal(result, expected)

        expected = to_timedelta([1, 3, 5, 7, 9], unit='D') + Second(2)
        result = timedelta_range('1 days, 00:00:02', periods=5, freq='2D')
        tm.assert_index_equal(result, expected)

        expected = to_timedelta(np.arange(50), unit='T') * 30
        result = timedelta_range('0 days', freq='30T', periods=50)
        tm.assert_index_equal(result, expected)

        # GH 11776
        arr = np.arange(10).reshape(2, 5)
        df = pd.DataFrame(np.arange(10).reshape(2, 5))
        for arg in (arr, df):
            with tm.assertRaisesRegexp(TypeError, "1-d array"):
                to_timedelta(arg)
            for errors in ['ignore', 'raise', 'coerce']:
                with tm.assertRaisesRegexp(TypeError, "1-d array"):
                    to_timedelta(arg, errors=errors)

        # issue10583
        df = pd.DataFrame(np.random.normal(size=(10, 4)))
        df.index = pd.timedelta_range(start='0s', periods=10, freq='s')
        expected = df.loc[pd.Timedelta('0s'):, :]
        result = df.loc['0s':, :]
        assert_frame_equal(expected, result)

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

    def test_nat_converters(self):
        self.assertEqual(to_timedelta(
            'nat', box=False).astype('int64'), tslib.iNaT)
        self.assertEqual(to_timedelta(
            'nan', box=False).astype('int64'), tslib.iNaT)

    def test_to_timedelta(self):
        def conv(v):
            return v.astype('m8[ns]')

        d1 = np.timedelta64(1, 'D')

        self.assertEqual(to_timedelta('1 days 06:05:01.00003', box=False),
                         conv(d1 + np.timedelta64(6 * 3600 +
                                                  5 * 60 + 1, 's') +
                              np.timedelta64(30, 'us')))
        self.assertEqual(to_timedelta('15.5us', box=False),
                         conv(np.timedelta64(15500, 'ns')))

        # empty string
        result = to_timedelta('', box=False)
        self.assertEqual(result.astype('int64'), tslib.iNaT)

        result = to_timedelta(['', ''])
        self.assertTrue(isnull(result).all())

        # pass thru
        result = to_timedelta(np.array([np.timedelta64(1, 's')]))
        expected = pd.Index(np.array([np.timedelta64(1, 's')]))
        tm.assert_index_equal(result, expected)

        # ints
        result = np.timedelta64(0, 'ns')
        expected = to_timedelta(0, box=False)
        self.assertEqual(result, expected)

        # Series
        expected = Series([timedelta(days=1), timedelta(days=1, seconds=1)])
        result = to_timedelta(Series(['1d', '1days 00:00:01']))
        tm.assert_series_equal(result, expected)

        # with units
        result = TimedeltaIndex([np.timedelta64(0, 'ns'), np.timedelta64(
            10, 's').astype('m8[ns]')])
        expected = to_timedelta([0, 10], unit='s')
        tm.assert_index_equal(result, expected)

        # single element conversion
        v = timedelta(seconds=1)
        result = to_timedelta(v, box=False)
        expected = np.timedelta64(timedelta(seconds=1))
        self.assertEqual(result, expected)

        v = np.timedelta64(timedelta(seconds=1))
        result = to_timedelta(v, box=False)
        expected = np.timedelta64(timedelta(seconds=1))
        self.assertEqual(result, expected)

        # arrays of various dtypes
        arr = np.array([1] * 5, dtype='int64')
        result = to_timedelta(arr, unit='s')
        expected = TimedeltaIndex([np.timedelta64(1, 's')] * 5)
        tm.assert_index_equal(result, expected)

        arr = np.array([1] * 5, dtype='int64')
        result = to_timedelta(arr, unit='m')
        expected = TimedeltaIndex([np.timedelta64(1, 'm')] * 5)
        tm.assert_index_equal(result, expected)

        arr = np.array([1] * 5, dtype='int64')
        result = to_timedelta(arr, unit='h')
        expected = TimedeltaIndex([np.timedelta64(1, 'h')] * 5)
        tm.assert_index_equal(result, expected)

        arr = np.array([1] * 5, dtype='timedelta64[s]')
        result = to_timedelta(arr)
        expected = TimedeltaIndex([np.timedelta64(1, 's')] * 5)
        tm.assert_index_equal(result, expected)

        arr = np.array([1] * 5, dtype='timedelta64[D]')
        result = to_timedelta(arr)
        expected = TimedeltaIndex([np.timedelta64(1, 'D')] * 5)
        tm.assert_index_equal(result, expected)

        # Test with lists as input when box=false
        expected = np.array(np.arange(3) * 1000000000, dtype='timedelta64[ns]')
        result = to_timedelta(range(3), unit='s', box=False)
        tm.assert_numpy_array_equal(expected, result)

        result = to_timedelta(np.arange(3), unit='s', box=False)
        tm.assert_numpy_array_equal(expected, result)

        result = to_timedelta([0, 1, 2], unit='s', box=False)
        tm.assert_numpy_array_equal(expected, result)

        # Tests with fractional seconds as input:
        expected = np.array(
            [0, 500000000, 800000000, 1200000000], dtype='timedelta64[ns]')
        result = to_timedelta([0., 0.5, 0.8, 1.2], unit='s', box=False)
        tm.assert_numpy_array_equal(expected, result)

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

    def test_to_timedelta_invalid(self):

        # these will error
        self.assertRaises(ValueError, lambda: to_timedelta([1, 2], unit='foo'))
        self.assertRaises(ValueError, lambda: to_timedelta(1, unit='foo'))

        # time not supported ATM
        self.assertRaises(ValueError, lambda: to_timedelta(time(second=1)))
        self.assertTrue(to_timedelta(
            time(second=1), errors='coerce') is pd.NaT)

        self.assertRaises(ValueError, lambda: to_timedelta(['foo', 'bar']))
        tm.assert_index_equal(TimedeltaIndex([pd.NaT, pd.NaT]),
                              to_timedelta(['foo', 'bar'], errors='coerce'))

        tm.assert_index_equal(TimedeltaIndex(['1 day', pd.NaT, '1 min']),
                              to_timedelta(['1 day', 'bar', '1 min'],
                                           errors='coerce'))

    def test_to_timedelta_via_apply(self):
        # GH 5458
        expected = Series([np.timedelta64(1, 's')])
        result = Series(['00:00:01']).apply(to_timedelta)
        tm.assert_series_equal(result, expected)

        result = Series([to_timedelta('00:00:01')])
        tm.assert_series_equal(result, expected)

    def test_timedelta_ops(self):
        # GH4984
        # make sure ops return Timedelta
        s = Series([Timestamp('20130101') + timedelta(seconds=i * i)
                    for i in range(10)])
        td = s.diff()

        result = td.mean()
        expected = to_timedelta(timedelta(seconds=9))
        self.assertEqual(result, expected)

        result = td.to_frame().mean()
        self.assertEqual(result[0], expected)

        result = td.quantile(.1)
        expected = Timedelta(np.timedelta64(2600, 'ms'))
        self.assertEqual(result, expected)

        result = td.median()
        expected = to_timedelta('00:00:09')
        self.assertEqual(result, expected)

        result = td.to_frame().median()
        self.assertEqual(result[0], expected)

        # GH 6462
        # consistency in returned values for sum
        result = td.sum()
        expected = to_timedelta('00:01:21')
        self.assertEqual(result, expected)

        result = td.to_frame().sum()
        self.assertEqual(result[0], expected)

        # std
        result = td.std()
        expected = to_timedelta(Series(td.dropna().values).std())
        self.assertEqual(result, expected)

        result = td.to_frame().std()
        self.assertEqual(result[0], expected)

        # invalid ops
        for op in ['skew', 'kurt', 'sem', 'prod']:
            self.assertRaises(TypeError, getattr(td, op))

        # GH 10040
        # make sure NaT is properly handled by median()
        s = Series([Timestamp('2015-02-03'), Timestamp('2015-02-07')])
        self.assertEqual(s.diff().median(), timedelta(days=4))

        s = Series([Timestamp('2015-02-03'), Timestamp('2015-02-07'),
                    Timestamp('2015-02-15')])
        self.assertEqual(s.diff().median(), timedelta(days=6))

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

    def test_timedelta_ops_scalar(self):
        # GH 6808
        base = pd.to_datetime('20130101 09:01:12.123456')
        expected_add = pd.to_datetime('20130101 09:01:22.123456')
        expected_sub = pd.to_datetime('20130101 09:01:02.123456')

        for offset in [pd.to_timedelta(10, unit='s'), timedelta(seconds=10),
                       np.timedelta64(10, 's'),
                       np.timedelta64(10000000000, 'ns'),
                       pd.offsets.Second(10)]:
            result = base + offset
            self.assertEqual(result, expected_add)

            result = base - offset
            self.assertEqual(result, expected_sub)

        base = pd.to_datetime('20130102 09:01:12.123456')
        expected_add = pd.to_datetime('20130103 09:01:22.123456')
        expected_sub = pd.to_datetime('20130101 09:01:02.123456')

        for offset in [pd.to_timedelta('1 day, 00:00:10'),
                       pd.to_timedelta('1 days, 00:00:10'),
                       timedelta(days=1, seconds=10),
                       np.timedelta64(1, 'D') + np.timedelta64(10, 's'),
                       pd.offsets.Day() + pd.offsets.Second(10)]:
            result = base + offset
            self.assertEqual(result, expected_add)

            result = base - offset
            self.assertEqual(result, expected_sub)

    def test_to_timedelta_on_missing_values(self):
        # GH5438
        timedelta_NaT = np.timedelta64('NaT')

        actual = pd.to_timedelta(Series(['00:00:01', np.nan]))
        expected = Series([np.timedelta64(1000000000, 'ns'),
                           timedelta_NaT], dtype='<m8[ns]')
        assert_series_equal(actual, expected)

        actual = pd.to_timedelta(Series(['00:00:01', pd.NaT]))
        assert_series_equal(actual, expected)

        actual = pd.to_timedelta(np.nan)
        self.assertEqual(actual.value, timedelta_NaT.astype('int64'))

        actual = pd.to_timedelta(pd.NaT)
        self.assertEqual(actual.value, timedelta_NaT.astype('int64'))

    def test_to_timedelta_on_nanoseconds(self):
        # GH 9273
        result = Timedelta(nanoseconds=100)
        expected = Timedelta('100ns')
        self.assertEqual(result, expected)

        result = Timedelta(days=1, hours=1, minutes=1, weeks=1, seconds=1,
                           milliseconds=1, microseconds=1, nanoseconds=1)
        expected = Timedelta(694861001001001)
        self.assertEqual(result, expected)

        result = Timedelta(microseconds=1) + Timedelta(nanoseconds=1)
        expected = Timedelta('1us1ns')
        self.assertEqual(result, expected)

        result = Timedelta(microseconds=1) - Timedelta(nanoseconds=1)
        expected = Timedelta('999ns')
        self.assertEqual(result, expected)

        result = Timedelta(microseconds=1) + 5 * Timedelta(nanoseconds=-2)
        expected = Timedelta('990ns')
        self.assertEqual(result, expected)

        self.assertRaises(TypeError, lambda: Timedelta(nanoseconds='abc'))

    def test_timedelta_ops_with_missing_values(self):
        # setup
        s1 = pd.to_timedelta(Series(['00:00:01']))
        s2 = pd.to_timedelta(Series(['00:00:02']))
        sn = pd.to_timedelta(Series([pd.NaT]))
        df1 = DataFrame(['00:00:01']).apply(pd.to_timedelta)
        df2 = DataFrame(['00:00:02']).apply(pd.to_timedelta)
        dfn = DataFrame([pd.NaT]).apply(pd.to_timedelta)
        scalar1 = pd.to_timedelta('00:00:01')
        scalar2 = pd.to_timedelta('00:00:02')
        timedelta_NaT = pd.to_timedelta('NaT')
        NA = np.nan

        actual = scalar1 + scalar1
        self.assertEqual(actual, scalar2)
        actual = scalar2 - scalar1
        self.assertEqual(actual, scalar1)

        actual = s1 + s1
        assert_series_equal(actual, s2)
        actual = s2 - s1
        assert_series_equal(actual, s1)

        actual = s1 + scalar1
        assert_series_equal(actual, s2)
        actual = scalar1 + s1
        assert_series_equal(actual, s2)
        actual = s2 - scalar1
        assert_series_equal(actual, s1)
        actual = -scalar1 + s2
        assert_series_equal(actual, s1)

        actual = s1 + timedelta_NaT
        assert_series_equal(actual, sn)
        actual = timedelta_NaT + s1
        assert_series_equal(actual, sn)
        actual = s1 - timedelta_NaT
        assert_series_equal(actual, sn)
        actual = -timedelta_NaT + s1
        assert_series_equal(actual, sn)

        actual = s1 + NA
        assert_series_equal(actual, sn)
        actual = NA + s1
        assert_series_equal(actual, sn)
        actual = s1 - NA
        assert_series_equal(actual, sn)
        actual = -NA + s1
        assert_series_equal(actual, sn)

        actual = s1 + pd.NaT
        assert_series_equal(actual, sn)
        actual = s2 - pd.NaT
        assert_series_equal(actual, sn)

        actual = s1 + df1
        assert_frame_equal(actual, df2)
        actual = s2 - df1
        assert_frame_equal(actual, df1)
        actual = df1 + s1
        assert_frame_equal(actual, df2)
        actual = df2 - s1
        assert_frame_equal(actual, df1)

        actual = df1 + df1
        assert_frame_equal(actual, df2)
        actual = df2 - df1
        assert_frame_equal(actual, df1)

        actual = df1 + scalar1
        assert_frame_equal(actual, df2)
        actual = df2 - scalar1
        assert_frame_equal(actual, df1)

        actual = df1 + timedelta_NaT
        assert_frame_equal(actual, dfn)
        actual = df1 - timedelta_NaT
        assert_frame_equal(actual, dfn)

        actual = df1 + NA
        assert_frame_equal(actual, dfn)
        actual = df1 - NA
        assert_frame_equal(actual, dfn)

        actual = df1 + pd.NaT  # NaT is datetime, not timedelta
        assert_frame_equal(actual, dfn)
        actual = df1 - pd.NaT
        assert_frame_equal(actual, dfn)

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


class TestTimedeltaIndex(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_pass_TimedeltaIndex_to_index(self):

        rng = timedelta_range('1 days', '10 days')
        idx = Index(rng, dtype=object)

        expected = Index(rng.to_pytimedelta(), dtype=object)

        self.assert_numpy_array_equal(idx.values, expected.values)

    def test_pickle(self):

        rng = timedelta_range('1 days', periods=10)
        rng_p = self.round_trip_pickle(rng)
        tm.assert_index_equal(rng, rng_p)

    def test_hash_error(self):
        index = timedelta_range('1 days', periods=10)
        with tm.assertRaisesRegexp(TypeError, "unhashable type: %r" %
                                   type(index).__name__):
            hash(index)

    def test_append_join_nondatetimeindex(self):
        rng = timedelta_range('1 days', periods=10)
        idx = Index(['a', 'b', 'c', 'd'])

        result = rng.append(idx)
        tm.assertIsInstance(result[0], Timedelta)

        # it works
        rng.join(idx, how='outer')

    def test_append_numpy_bug_1681(self):

        td = timedelta_range('1 days', '10 days', freq='2D')
        a = DataFrame()
        c = DataFrame({'A': 'foo', 'B': td}, index=td)
        str(c)

        result = a.append(c)
        self.assertTrue((result['B'] == td).all())

    def test_astype(self):
        rng = timedelta_range('1 days', periods=10)

        result = rng.astype('i8')
        self.assert_numpy_array_equal(result, rng.asi8)

    def test_fields(self):
        rng = timedelta_range('1 days, 10:11:12.100123456', periods=2,
                              freq='s')
        self.assert_numpy_array_equal(rng.days, np.array(
            [1, 1], dtype='int64'))
        self.assert_numpy_array_equal(
            rng.seconds,
            np.array([10 * 3600 + 11 * 60 + 12, 10 * 3600 + 11 * 60 + 13],
                     dtype='int64'))
        self.assert_numpy_array_equal(rng.microseconds, np.array(
            [100 * 1000 + 123, 100 * 1000 + 123], dtype='int64'))
        self.assert_numpy_array_equal(rng.nanoseconds, np.array(
            [456, 456], dtype='int64'))

        self.assertRaises(AttributeError, lambda: rng.hours)
        self.assertRaises(AttributeError, lambda: rng.minutes)
        self.assertRaises(AttributeError, lambda: rng.milliseconds)

        # with nat
        s = Series(rng)
        s[1] = np.nan

        tm.assert_series_equal(s.dt.days, Series([1, np.nan], index=[0, 1]))
        tm.assert_series_equal(s.dt.seconds, Series(
            [10 * 3600 + 11 * 60 + 12, np.nan], index=[0, 1]))

    def test_total_seconds(self):
        # GH 10939
        # test index
        rng = timedelta_range('1 days, 10:11:12.100123456', periods=2,
                              freq='s')
        expt = [1 * 86400 + 10 * 3600 + 11 * 60 + 12 + 100123456. / 1e9,
                1 * 86400 + 10 * 3600 + 11 * 60 + 13 + 100123456. / 1e9]
        assert_allclose(rng.total_seconds(), expt, atol=1e-10, rtol=0)

        # test Series
        s = Series(rng)
        s_expt = Series(expt, index=[0, 1])
        tm.assert_series_equal(s.dt.total_seconds(), s_expt)

        # with nat
        s[1] = np.nan
        s_expt = Series([1 * 86400 + 10 * 3600 + 11 * 60 +
                         12 + 100123456. / 1e9, np.nan], index=[0, 1])
        tm.assert_series_equal(s.dt.total_seconds(), s_expt)

        # with both nat
        s = Series([np.nan, np.nan], dtype='timedelta64[ns]')
        tm.assert_series_equal(s.dt.total_seconds(), Series(
            [np.nan, np.nan], index=[0, 1]))

    def test_total_seconds_scalar(self):
        # GH 10939
        rng = Timedelta('1 days, 10:11:12.100123456')
        expt = 1 * 86400 + 10 * 3600 + 11 * 60 + 12 + 100123456. / 1e9
        assert_allclose(rng.total_seconds(), expt, atol=1e-10, rtol=0)

        rng = Timedelta(np.nan)
        self.assertTrue(np.isnan(rng.total_seconds()))

    def test_components(self):
        rng = timedelta_range('1 days, 10:11:12', periods=2, freq='s')
        rng.components

        # with nat
        s = Series(rng)
        s[1] = np.nan

        result = s.dt.components
        self.assertFalse(result.iloc[0].isnull().all())
        self.assertTrue(result.iloc[1].isnull().all())

    def test_constructor(self):
        expected = TimedeltaIndex(['1 days', '1 days 00:00:05', '2 days',
                                   '2 days 00:00:02', '0 days 00:00:03'])
        result = TimedeltaIndex(['1 days', '1 days, 00:00:05', np.timedelta64(
            2, 'D'), timedelta(days=2, seconds=2), pd.offsets.Second(3)])
        tm.assert_index_equal(result, expected)

        # unicode
        result = TimedeltaIndex([u'1 days', '1 days, 00:00:05', np.timedelta64(
            2, 'D'), timedelta(days=2, seconds=2), pd.offsets.Second(3)])

        expected = TimedeltaIndex(['0 days 00:00:00', '0 days 00:00:01',
                                   '0 days 00:00:02'])
        tm.assert_index_equal(TimedeltaIndex(range(3), unit='s'), expected)
        expected = TimedeltaIndex(['0 days 00:00:00', '0 days 00:00:05',
                                   '0 days 00:00:09'])
        tm.assert_index_equal(TimedeltaIndex([0, 5, 9], unit='s'), expected)
        expected = TimedeltaIndex(
            ['0 days 00:00:00.400', '0 days 00:00:00.450',
             '0 days 00:00:01.200'])
        tm.assert_index_equal(TimedeltaIndex([400, 450, 1200], unit='ms'),
                              expected)

    def test_constructor_coverage(self):
        rng = timedelta_range('1 days', periods=10.5)
        exp = timedelta_range('1 days', periods=10)
        self.assertTrue(rng.equals(exp))

        self.assertRaises(ValueError, TimedeltaIndex, start='1 days',
                          periods='foo', freq='D')

        self.assertRaises(ValueError, TimedeltaIndex, start='1 days',
                          end='10 days')

        self.assertRaises(ValueError, TimedeltaIndex, '1 days')

        # generator expression
        gen = (timedelta(i) for i in range(10))
        result = TimedeltaIndex(gen)
        expected = TimedeltaIndex([timedelta(i) for i in range(10)])
        self.assertTrue(result.equals(expected))

        # NumPy string array
        strings = np.array(['1 days', '2 days', '3 days'])
        result = TimedeltaIndex(strings)
        expected = to_timedelta([1, 2, 3], unit='d')
        self.assertTrue(result.equals(expected))

        from_ints = TimedeltaIndex(expected.asi8)
        self.assertTrue(from_ints.equals(expected))

        # non-conforming freq
        self.assertRaises(ValueError, TimedeltaIndex,
                          ['1 days', '2 days', '4 days'], freq='D')

        self.assertRaises(ValueError, TimedeltaIndex, periods=10, freq='D')

    def test_constructor_name(self):
        idx = TimedeltaIndex(start='1 days', periods=1, freq='D', name='TEST')
        self.assertEqual(idx.name, 'TEST')

        # GH10025
        idx2 = TimedeltaIndex(idx, name='something else')
        self.assertEqual(idx2.name, 'something else')

    def test_freq_conversion(self):

        # doc example

        # series
        td = Series(date_range('20130101', periods=4)) - \
            Series(date_range('20121201', periods=4))
        td[2] += timedelta(minutes=5, seconds=3)
        td[3] = np.nan

        result = td / np.timedelta64(1, 'D')
        expected = Series([31, 31, (31 * 86400 + 5 * 60 + 3) / 86400.0, np.nan
                           ])
        assert_series_equal(result, expected)

        result = td.astype('timedelta64[D]')
        expected = Series([31, 31, 31, np.nan])
        assert_series_equal(result, expected)

        result = td / np.timedelta64(1, 's')
        expected = Series([31 * 86400, 31 * 86400, 31 * 86400 + 5 * 60 + 3,
                           np.nan])
        assert_series_equal(result, expected)

        result = td.astype('timedelta64[s]')
        assert_series_equal(result, expected)

        # tdi
        td = TimedeltaIndex(td)

        result = td / np.timedelta64(1, 'D')
        expected = Index([31, 31, (31 * 86400 + 5 * 60 + 3) / 86400.0, np.nan])
        assert_index_equal(result, expected)

        result = td.astype('timedelta64[D]')
        expected = Index([31, 31, 31, np.nan])
        assert_index_equal(result, expected)

        result = td / np.timedelta64(1, 's')
        expected = Index([31 * 86400, 31 * 86400, 31 * 86400 + 5 * 60 + 3,
                          np.nan])
        assert_index_equal(result, expected)

        result = td.astype('timedelta64[s]')
        assert_index_equal(result, expected)

    def test_comparisons_coverage(self):
        rng = timedelta_range('1 days', periods=10)

        result = rng < rng[3]
        exp = np.array([True, True, True] + [False] * 7)
        self.assert_numpy_array_equal(result, exp)

        # raise TypeError for now
        self.assertRaises(TypeError, rng.__lt__, rng[3].value)

        result = rng == list(rng)
        exp = rng == rng
        self.assert_numpy_array_equal(result, exp)

    def test_comparisons_nat(self):

        tdidx1 = pd.TimedeltaIndex(['1 day', pd.NaT, '1 day 00:00:01', pd.NaT,
                                    '1 day 00:00:01', '5 day 00:00:03'])
        tdidx2 = pd.TimedeltaIndex(['2 day', '2 day', pd.NaT, pd.NaT,
                                    '1 day 00:00:02', '5 days 00:00:03'])
        tdarr = np.array([np.timedelta64(2, 'D'),
                          np.timedelta64(2, 'D'), np.timedelta64('nat'),
                          np.timedelta64('nat'),
                          np.timedelta64(1, 'D') + np.timedelta64(2, 's'),
                          np.timedelta64(5, 'D') + np.timedelta64(3, 's')])

        if _np_version_under1p8:
            # cannot test array because np.datetime('nat') returns today's date
            cases = [(tdidx1, tdidx2)]
        else:
            cases = [(tdidx1, tdidx2), (tdidx1, tdarr)]

        # Check pd.NaT is handles as the same as np.nan
        for idx1, idx2 in cases:

            result = idx1 < idx2
            expected = np.array([True, False, False, False, True, False])
            self.assert_numpy_array_equal(result, expected)

            result = idx2 > idx1
            expected = np.array([True, False, False, False, True, False])
            self.assert_numpy_array_equal(result, expected)

            result = idx1 <= idx2
            expected = np.array([True, False, False, False, True, True])
            self.assert_numpy_array_equal(result, expected)

            result = idx2 >= idx1
            expected = np.array([True, False, False, False, True, True])
            self.assert_numpy_array_equal(result, expected)

            result = idx1 == idx2
            expected = np.array([False, False, False, False, False, True])
            self.assert_numpy_array_equal(result, expected)

            result = idx1 != idx2
            expected = np.array([True, True, True, True, True, False])
            self.assert_numpy_array_equal(result, expected)

    def test_map(self):

        rng = timedelta_range('1 day', periods=10)

        f = lambda x: x.days
        result = rng.map(f)
        exp = [f(x) for x in rng]
        self.assert_numpy_array_equal(result, exp)

    def test_misc_coverage(self):

        rng = timedelta_range('1 day', periods=5)
        result = rng.groupby(rng.days)
        tm.assertIsInstance(list(result.values())[0][0], Timedelta)

        idx = TimedeltaIndex(['3d', '1d', '2d'])
        self.assertTrue(idx.equals(list(idx)))

        non_td = Index(list('abc'))
        self.assertFalse(idx.equals(list(non_td)))

    def test_union(self):

        i1 = timedelta_range('1day', periods=5)
        i2 = timedelta_range('3day', periods=5)
        result = i1.union(i2)
        expected = timedelta_range('1day', periods=7)
        self.assert_numpy_array_equal(result, expected)

        i1 = Int64Index(np.arange(0, 20, 2))
        i2 = TimedeltaIndex(start='1 day', periods=10, freq='D')
        i1.union(i2)  # Works
        i2.union(i1)  # Fails with "AttributeError: can't set attribute"

    def test_union_coverage(self):

        idx = TimedeltaIndex(['3d', '1d', '2d'])
        ordered = TimedeltaIndex(idx.sort_values(), freq='infer')
        result = ordered.union(idx)
        self.assertTrue(result.equals(ordered))

        result = ordered[:0].union(ordered)
        self.assertTrue(result.equals(ordered))
        self.assertEqual(result.freq, ordered.freq)

    def test_union_bug_1730(self):

        rng_a = timedelta_range('1 day', periods=4, freq='3H')
        rng_b = timedelta_range('1 day', periods=4, freq='4H')

        result = rng_a.union(rng_b)
        exp = TimedeltaIndex(sorted(set(list(rng_a)) | set(list(rng_b))))
        self.assertTrue(result.equals(exp))

    def test_union_bug_1745(self):

        left = TimedeltaIndex(['1 day 15:19:49.695000'])
        right = TimedeltaIndex(
            ['2 day 13:04:21.322000', '1 day 15:27:24.873000',
             '1 day 15:31:05.350000'])

        result = left.union(right)
        exp = TimedeltaIndex(sorted(set(list(left)) | set(list(right))))
        self.assertTrue(result.equals(exp))

    def test_union_bug_4564(self):

        left = timedelta_range("1 day", "30d")
        right = left + pd.offsets.Minute(15)

        result = left.union(right)
        exp = TimedeltaIndex(sorted(set(list(left)) | set(list(right))))
        self.assertTrue(result.equals(exp))

    def test_intersection_bug_1708(self):
        index_1 = timedelta_range('1 day', periods=4, freq='h')
        index_2 = index_1 + pd.offsets.Hour(5)

        result = index_1 & index_2
        self.assertEqual(len(result), 0)

        index_1 = timedelta_range('1 day', periods=4, freq='h')
        index_2 = index_1 + pd.offsets.Hour(1)

        result = index_1 & index_2
        expected = timedelta_range('1 day 01:00:00', periods=3, freq='h')
        tm.assert_index_equal(result, expected)

    def test_get_duplicates(self):
        idx = TimedeltaIndex(['1 day', '2 day', '2 day', '3 day', '3day',
                              '4day'])

        result = idx.get_duplicates()
        ex = TimedeltaIndex(['2 day', '3day'])
        self.assertTrue(result.equals(ex))

    def test_argmin_argmax(self):
        idx = TimedeltaIndex(['1 day 00:00:05', '1 day 00:00:01',
                              '1 day 00:00:02'])
        self.assertEqual(idx.argmin(), 1)
        self.assertEqual(idx.argmax(), 0)

    def test_sort_values(self):

        idx = TimedeltaIndex(['4d', '1d', '2d'])

        ordered = idx.sort_values()
        self.assertTrue(ordered.is_monotonic)

        ordered = idx.sort_values(ascending=False)
        self.assertTrue(ordered[::-1].is_monotonic)

        ordered, dexer = idx.sort_values(return_indexer=True)
        self.assertTrue(ordered.is_monotonic)
        self.assert_numpy_array_equal(dexer, [1, 2, 0])

        ordered, dexer = idx.sort_values(return_indexer=True, ascending=False)
        self.assertTrue(ordered[::-1].is_monotonic)
        self.assert_numpy_array_equal(dexer, [0, 2, 1])

    def test_insert(self):

        idx = TimedeltaIndex(['4day', '1day', '2day'], name='idx')

        result = idx.insert(2, timedelta(days=5))
        exp = TimedeltaIndex(['4day', '1day', '5day', '2day'], name='idx')
        self.assertTrue(result.equals(exp))

        # insertion of non-datetime should coerce to object index
        result = idx.insert(1, 'inserted')
        expected = Index([Timedelta('4day'), 'inserted', Timedelta('1day'),
                          Timedelta('2day')], name='idx')
        self.assertNotIsInstance(result, TimedeltaIndex)
        tm.assert_index_equal(result, expected)
        self.assertEqual(result.name, expected.name)

        idx = timedelta_range('1day 00:00:01', periods=3, freq='s', name='idx')

        # preserve freq
        expected_0 = TimedeltaIndex(['1day', '1day 00:00:01', '1day 00:00:02',
                                     '1day 00:00:03'],
                                    name='idx', freq='s')
        expected_3 = TimedeltaIndex(['1day 00:00:01', '1day 00:00:02',
                                     '1day 00:00:03', '1day 00:00:04'],
                                    name='idx', freq='s')

        # reset freq to None
        expected_1_nofreq = TimedeltaIndex(['1day 00:00:01', '1day 00:00:01',
                                            '1day 00:00:02', '1day 00:00:03'],
                                           name='idx', freq=None)
        expected_3_nofreq = TimedeltaIndex(['1day 00:00:01', '1day 00:00:02',
                                            '1day 00:00:03', '1day 00:00:05'],
                                           name='idx', freq=None)

        cases = [(0, Timedelta('1day'), expected_0),
                 (-3, Timedelta('1day'), expected_0),
                 (3, Timedelta('1day 00:00:04'), expected_3),
                 (1, Timedelta('1day 00:00:01'), expected_1_nofreq),
                 (3, Timedelta('1day 00:00:05'), expected_3_nofreq)]

        for n, d, expected in cases:
            result = idx.insert(n, d)
            self.assertTrue(result.equals(expected))
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)

    def test_delete(self):
        idx = timedelta_range(start='1 Days', periods=5, freq='D', name='idx')

        # prserve freq
        expected_0 = timedelta_range(start='2 Days', periods=4, freq='D',
                                     name='idx')
        expected_4 = timedelta_range(start='1 Days', periods=4, freq='D',
                                     name='idx')

        # reset freq to None
        expected_1 = TimedeltaIndex(
            ['1 day', '3 day', '4 day', '5 day'], freq=None, name='idx')

        cases = {0: expected_0,
                 -5: expected_0,
                 -1: expected_4,
                 4: expected_4,
                 1: expected_1}
        for n, expected in compat.iteritems(cases):
            result = idx.delete(n)
            self.assertTrue(result.equals(expected))
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)

        with tm.assertRaises((IndexError, ValueError)):
            # either depeidnig on numpy version
            result = idx.delete(5)

    def test_delete_slice(self):
        idx = timedelta_range(start='1 days', periods=10, freq='D', name='idx')

        # prserve freq
        expected_0_2 = timedelta_range(start='4 days', periods=7, freq='D',
                                       name='idx')
        expected_7_9 = timedelta_range(start='1 days', periods=7, freq='D',
                                       name='idx')

        # reset freq to None
        expected_3_5 = TimedeltaIndex(['1 d', '2 d', '3 d',
                                       '7 d', '8 d', '9 d', '10d'],
                                      freq=None, name='idx')

        cases = {(0, 1, 2): expected_0_2,
                 (7, 8, 9): expected_7_9,
                 (3, 4, 5): expected_3_5}
        for n, expected in compat.iteritems(cases):
            result = idx.delete(n)
            self.assertTrue(result.equals(expected))
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)

            result = idx.delete(slice(n[0], n[-1] + 1))
            self.assertTrue(result.equals(expected))
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)

    def test_take(self):

        tds = ['1day 02:00:00', '1 day 04:00:00', '1 day 10:00:00']
        idx = TimedeltaIndex(start='1d', end='2d', freq='H', name='idx')
        expected = TimedeltaIndex(tds, freq=None, name='idx')

        taken1 = idx.take([2, 4, 10])
        taken2 = idx[[2, 4, 10]]

        for taken in [taken1, taken2]:
            self.assertTrue(taken.equals(expected))
            tm.assertIsInstance(taken, TimedeltaIndex)
            self.assertIsNone(taken.freq)
            self.assertEqual(taken.name, expected.name)

    def test_take_fill_value(self):
        # GH 12631
        idx = pd.TimedeltaIndex(['1 days', '2 days', '3 days'],
                                name='xxx')
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.TimedeltaIndex(['2 days', '1 days', '3 days'],
                                     name='xxx')
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.TimedeltaIndex(['2 days', '1 days', 'NaT'],
                                     name='xxx')
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        expected = pd.TimedeltaIndex(['2 days', '1 days', '3 days'],
                                     name='xxx')
        tm.assert_index_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with tm.assertRaises(IndexError):
            idx.take(np.array([1, -5]))

    def test_isin(self):

        index = tm.makeTimedeltaIndex(4)
        result = index.isin(index)
        self.assertTrue(result.all())

        result = index.isin(list(index))
        self.assertTrue(result.all())

        assert_almost_equal(index.isin([index[2], 5]),
                            [False, False, True, False])

    def test_does_not_convert_mixed_integer(self):
        df = tm.makeCustomDataframe(10, 10,
                                    data_gen_f=lambda *args, **kwargs: randn(),
                                    r_idx_type='i', c_idx_type='td')
        str(df)

        cols = df.columns.join(df.index, how='outer')
        joined = cols.join(df.columns)
        self.assertEqual(cols.dtype, np.dtype('O'))
        self.assertEqual(cols.dtype, joined.dtype)
        tm.assert_index_equal(cols, joined)

    def test_slice_keeps_name(self):

        # GH4226
        dr = pd.timedelta_range('1d', '5d', freq='H', name='timebucket')
        self.assertEqual(dr[1:].name, dr.name)

    def test_join_self(self):

        index = timedelta_range('1 day', periods=10)
        kinds = 'outer', 'inner', 'left', 'right'
        for kind in kinds:
            joined = index.join(index, how=kind)
            self.assertIs(index, joined)

    def test_factorize(self):
        idx1 = TimedeltaIndex(['1 day', '1 day', '2 day', '2 day', '3 day',
                               '3 day'])

        exp_arr = np.array([0, 0, 1, 1, 2, 2])
        exp_idx = TimedeltaIndex(['1 day', '2 day', '3 day'])

        arr, idx = idx1.factorize()
        self.assert_numpy_array_equal(arr, exp_arr)
        self.assertTrue(idx.equals(exp_idx))

        arr, idx = idx1.factorize(sort=True)
        self.assert_numpy_array_equal(arr, exp_arr)
        self.assertTrue(idx.equals(exp_idx))

        # freq must be preserved
        idx3 = timedelta_range('1 day', periods=4, freq='s')
        exp_arr = np.array([0, 1, 2, 3])
        arr, idx = idx3.factorize()
        self.assert_numpy_array_equal(arr, exp_arr)
        self.assertTrue(idx.equals(idx3))


class TestSlicing(tm.TestCase):
    def test_partial_slice(self):
        rng = timedelta_range('1 day 10:11:12', freq='h', periods=500)
        s = Series(np.arange(len(rng)), index=rng)

        result = s['5 day':'6 day']
        expected = s.iloc[86:134]
        assert_series_equal(result, expected)

        result = s['5 day':]
        expected = s.iloc[86:]
        assert_series_equal(result, expected)

        result = s[:'6 day']
        expected = s.iloc[:134]
        assert_series_equal(result, expected)

        result = s['6 days, 23:11:12']
        self.assertEqual(result, s.iloc[133])

        self.assertRaises(KeyError, s.__getitem__, '50 days')

    def test_partial_slice_high_reso(self):

        # higher reso
        rng = timedelta_range('1 day 10:11:12', freq='us', periods=2000)
        s = Series(np.arange(len(rng)), index=rng)

        result = s['1 day 10:11:12':]
        expected = s.iloc[0:]
        assert_series_equal(result, expected)

        result = s['1 day 10:11:12.001':]
        expected = s.iloc[1000:]
        assert_series_equal(result, expected)

        result = s['1 days, 10:11:12.001001']
        self.assertEqual(result, s.iloc[1001])

    def test_slice_with_negative_step(self):
        ts = Series(np.arange(20), timedelta_range('0', periods=20, freq='H'))
        SLC = pd.IndexSlice

        def assert_slices_equivalent(l_slc, i_slc):
            assert_series_equal(ts[l_slc], ts.iloc[i_slc])
            assert_series_equal(ts.loc[l_slc], ts.iloc[i_slc])
            assert_series_equal(ts.ix[l_slc], ts.iloc[i_slc])

        assert_slices_equivalent(SLC[Timedelta(hours=7)::-1], SLC[7::-1])
        assert_slices_equivalent(SLC['7 hours'::-1], SLC[7::-1])

        assert_slices_equivalent(SLC[:Timedelta(hours=7):-1], SLC[:6:-1])
        assert_slices_equivalent(SLC[:'7 hours':-1], SLC[:6:-1])

        assert_slices_equivalent(SLC['15 hours':'7 hours':-1], SLC[15:6:-1])
        assert_slices_equivalent(SLC[Timedelta(hours=15):Timedelta(hours=7):-
                                     1], SLC[15:6:-1])
        assert_slices_equivalent(SLC['15 hours':Timedelta(hours=7):-1],
                                 SLC[15:6:-1])
        assert_slices_equivalent(SLC[Timedelta(hours=15):'7 hours':-1],
                                 SLC[15:6:-1])

        assert_slices_equivalent(SLC['7 hours':'15 hours':-1], SLC[:0])

    def test_slice_with_zero_step_raises(self):
        ts = Series(np.arange(20), timedelta_range('0', periods=20, freq='H'))
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: ts[::0])
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: ts.loc[::0])
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: ts.ix[::0])

    def test_tdi_ops_attributes(self):
        rng = timedelta_range('2 days', periods=5, freq='2D', name='x')

        result = rng + 1
        exp = timedelta_range('4 days', periods=5, freq='2D', name='x')
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.freq, '2D')

        result = rng - 2
        exp = timedelta_range('-2 days', periods=5, freq='2D', name='x')
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.freq, '2D')

        result = rng * 2
        exp = timedelta_range('4 days', periods=5, freq='4D', name='x')
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.freq, '4D')

        result = rng / 2
        exp = timedelta_range('1 days', periods=5, freq='D', name='x')
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.freq, 'D')

        result = -rng
        exp = timedelta_range('-2 days', periods=5, freq='-2D', name='x')
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.freq, '-2D')

        rng = pd.timedelta_range('-2 days', periods=5, freq='D', name='x')

        result = abs(rng)
        exp = TimedeltaIndex(['2 days', '1 days', '0 days', '1 days',
                              '2 days'], name='x')
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.freq, None)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
