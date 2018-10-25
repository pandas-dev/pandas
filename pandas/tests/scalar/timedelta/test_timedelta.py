"""
Test the scalar Timedelta
These tests are intended to isolate pandas/_libs/tslibs; tests that
depend on implementations in Index/Series/etc belong in test_vector_compat.py
"""
from datetime import timedelta

import pytest
import numpy as np

import pandas.util.testing as tm
from pandas.core.tools.timedeltas import _coerce_scalar_to_timedelta_type as ct
from pandas import Timedelta, to_timedelta, compat, Timestamp, NaT


# TODO: belongs in test_arithmetic.py?
class TestTimedeltaArithmetic(object):

    def test_arithmetic_overflow(self):
        with pytest.raises(OverflowError):
            Timestamp('1700-01-01') + Timedelta(13 * 19999, unit='D')

        with pytest.raises(OverflowError):
            Timestamp('1700-01-01') + timedelta(days=13 * 19999)

    def test_ops_error_str(self):
        # GH#13624
        td = Timedelta('1 day')

        for left, right in [(td, 'a'), ('a', td)]:

            with pytest.raises(TypeError):
                left + right

            # GH 20829: python 2 comparison naturally does not raise TypeError
            if compat.PY3:
                with pytest.raises(TypeError):
                    left > right

            assert not left == right
            assert left != right

    def test_ops_notimplemented(self):
        class Other(object):
            pass

        other = Other()

        td = Timedelta('1 day')
        assert td.__add__(other) is NotImplemented
        assert td.__sub__(other) is NotImplemented
        assert td.__truediv__(other) is NotImplemented
        assert td.__mul__(other) is NotImplemented
        assert td.__floordiv__(other) is NotImplemented

    def test_unary_ops(self):
        td = Timedelta(10, unit='d')

        # __neg__, __pos__
        assert -td == Timedelta(-10, unit='d')
        assert -td == Timedelta('-10d')
        assert +td == Timedelta(10, unit='d')

        # __abs__, __abs__(__neg__)
        assert abs(td) == td
        assert abs(-td) == td
        assert abs(-td) == Timedelta('10d')


class TestTimedeltaComparison(object):
    def test_comparison_object_array(self):
        # analogous to GH#15183
        td = Timedelta('2 days')
        other = Timedelta('3 hours')

        arr = np.array([other, td], dtype=object)
        res = arr == td
        expected = np.array([False, True], dtype=bool)
        assert (res == expected).all()

        # 2D case
        arr = np.array([[other, td],
                        [td, other]],
                       dtype=object)
        res = arr != td
        expected = np.array([[True, False], [False, True]], dtype=bool)
        assert res.shape == expected.shape
        assert (res == expected).all()

    def test_compare_timedelta_ndarray(self):
        # GH11835
        periods = [Timedelta('0 days 01:00:00'), Timedelta('0 days 01:00:00')]
        arr = np.array(periods)
        result = arr[0] > arr
        expected = np.array([False, False])
        tm.assert_numpy_array_equal(result, expected)

    def test_compare_custom_object(self):
        """Make sure non supported operations on Timedelta returns NonImplemented
        and yields to other operand (GH20829)."""
        class CustomClass(object):

            def __init__(self, cmp_result=None):
                self.cmp_result = cmp_result

            def generic_result(self):
                if self.cmp_result is None:
                    return NotImplemented
                else:
                    return self.cmp_result

            def __eq__(self, other):
                return self.generic_result()

            def __gt__(self, other):
                return self.generic_result()

        t = Timedelta('1s')

        assert not (t == "string")
        assert not (t == 1)
        assert not (t == CustomClass())
        assert not (t == CustomClass(cmp_result=False))

        assert t < CustomClass(cmp_result=True)
        assert not (t < CustomClass(cmp_result=False))

        assert t == CustomClass(cmp_result=True)

    @pytest.mark.skipif(compat.PY2,
                        reason="python 2 does not raise TypeError for "
                               "comparisons of different types")
    @pytest.mark.parametrize("val", ["string", 1])
    def test_compare_unknown_type(self, val):
        # GH#20829
        t = Timedelta('1s')
        with pytest.raises(TypeError):
            t >= val
        with pytest.raises(TypeError):
            t > val
        with pytest.raises(TypeError):
            t <= val
        with pytest.raises(TypeError):
            t < val


class TestTimedeltas(object):

    @pytest.mark.parametrize("unit, value, expected", [
        ('us', 9.999, 9999), ('ms', 9.999999, 9999999),
        ('s', 9.999999999, 9999999999)])
    def test_rounding_on_int_unit_construction(self, unit, value, expected):
        # GH 12690
        result = Timedelta(value, unit=unit)
        assert result.value == expected
        result = Timedelta(str(value) + unit)
        assert result.value == expected

    def test_total_seconds_scalar(self):
        # see gh-10939
        rng = Timedelta('1 days, 10:11:12.100123456')
        expt = 1 * 86400 + 10 * 3600 + 11 * 60 + 12 + 100123456. / 1e9
        tm.assert_almost_equal(rng.total_seconds(), expt)

        rng = Timedelta(np.nan)
        assert np.isnan(rng.total_seconds())

    def test_conversion(self):

        for td in [Timedelta(10, unit='d'),
                   Timedelta('1 days, 10:11:12.012345')]:
            pydt = td.to_pytimedelta()
            assert td == Timedelta(pydt)
            assert td == pydt
            assert (isinstance(pydt, timedelta) and not isinstance(
                pydt, Timedelta))

            assert td == np.timedelta64(td.value, 'ns')
            td64 = td.to_timedelta64()

            assert td64 == np.timedelta64(td.value, 'ns')
            assert td == td64

            assert isinstance(td64, np.timedelta64)

        # this is NOT equal and cannot be roundtriped (because of the nanos)
        td = Timedelta('1 days, 10:11:12.012345678')
        assert td != td.to_pytimedelta()

    def test_freq_conversion(self):

        # truediv
        td = Timedelta('1 days 2 hours 3 ns')
        result = td / np.timedelta64(1, 'D')
        assert result == td.value / float(86400 * 1e9)
        result = td / np.timedelta64(1, 's')
        assert result == td.value / float(1e9)
        result = td / np.timedelta64(1, 'ns')
        assert result == td.value

        # floordiv
        td = Timedelta('1 days 2 hours 3 ns')
        result = td // np.timedelta64(1, 'D')
        assert result == 1
        result = td // np.timedelta64(1, 's')
        assert result == 93600
        result = td // np.timedelta64(1, 'ns')
        assert result == td.value

    def test_fields(self):
        def check(value):
            # that we are int/long like
            assert isinstance(value, (int, compat.long))

        # compat to datetime.timedelta
        rng = to_timedelta('1 days, 10:11:12')
        assert rng.days == 1
        assert rng.seconds == 10 * 3600 + 11 * 60 + 12
        assert rng.microseconds == 0
        assert rng.nanoseconds == 0

        pytest.raises(AttributeError, lambda: rng.hours)
        pytest.raises(AttributeError, lambda: rng.minutes)
        pytest.raises(AttributeError, lambda: rng.milliseconds)

        # GH 10050
        check(rng.days)
        check(rng.seconds)
        check(rng.microseconds)
        check(rng.nanoseconds)

        td = Timedelta('-1 days, 10:11:12')
        assert abs(td) == Timedelta('13:48:48')
        assert str(td) == "-1 days +10:11:12"
        assert -td == Timedelta('0 days 13:48:48')
        assert -Timedelta('-1 days, 10:11:12').value == 49728000000000
        assert Timedelta('-1 days, 10:11:12').value == -49728000000000

        rng = to_timedelta('-1 days, 10:11:12.100123456')
        assert rng.days == -1
        assert rng.seconds == 10 * 3600 + 11 * 60 + 12
        assert rng.microseconds == 100 * 1000 + 123
        assert rng.nanoseconds == 456
        pytest.raises(AttributeError, lambda: rng.hours)
        pytest.raises(AttributeError, lambda: rng.minutes)
        pytest.raises(AttributeError, lambda: rng.milliseconds)

        # components
        tup = to_timedelta(-1, 'us').components
        assert tup.days == -1
        assert tup.hours == 23
        assert tup.minutes == 59
        assert tup.seconds == 59
        assert tup.milliseconds == 999
        assert tup.microseconds == 999
        assert tup.nanoseconds == 0

        # GH 10050
        check(tup.days)
        check(tup.hours)
        check(tup.minutes)
        check(tup.seconds)
        check(tup.milliseconds)
        check(tup.microseconds)
        check(tup.nanoseconds)

        tup = Timedelta('-1 days 1 us').components
        assert tup.days == -2
        assert tup.hours == 23
        assert tup.minutes == 59
        assert tup.seconds == 59
        assert tup.milliseconds == 999
        assert tup.microseconds == 999
        assert tup.nanoseconds == 0

    def test_iso_conversion(self):
        # GH #21877
        expected = Timedelta(1, unit='s')
        assert to_timedelta('P0DT0H0M1S') == expected

    def test_numeric_conversions(self):
        assert ct(0) == np.timedelta64(0, 'ns')
        assert ct(10) == np.timedelta64(10, 'ns')
        assert ct(10, unit='ns') == np.timedelta64(10, 'ns').astype('m8[ns]')

        assert ct(10, unit='us') == np.timedelta64(10, 'us').astype('m8[ns]')
        assert ct(10, unit='ms') == np.timedelta64(10, 'ms').astype('m8[ns]')
        assert ct(10, unit='s') == np.timedelta64(10, 's').astype('m8[ns]')
        assert ct(10, unit='d') == np.timedelta64(10, 'D').astype('m8[ns]')

    def test_timedelta_conversions(self):
        assert (ct(timedelta(seconds=1)) ==
                np.timedelta64(1, 's').astype('m8[ns]'))
        assert (ct(timedelta(microseconds=1)) ==
                np.timedelta64(1, 'us').astype('m8[ns]'))
        assert (ct(timedelta(days=1)) ==
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
            assert r1 == s1
            r2 = t2.round(freq)
            assert r2 == s2

        # invalid
        for freq in ['Y', 'M', 'foobar']:
            with pytest.raises(ValueError):
                t1.round(freq)

    def test_identity(self):

        td = Timedelta(10, unit='d')
        assert isinstance(td, Timedelta)
        assert isinstance(td, timedelta)

    def test_short_format_converters(self):
        def conv(v):
            return v.astype('m8[ns]')

        assert ct('10') == np.timedelta64(10, 'ns')
        assert ct('10ns') == np.timedelta64(10, 'ns')
        assert ct('100') == np.timedelta64(100, 'ns')
        assert ct('100ns') == np.timedelta64(100, 'ns')

        assert ct('1000') == np.timedelta64(1000, 'ns')
        assert ct('1000ns') == np.timedelta64(1000, 'ns')
        assert ct('1000NS') == np.timedelta64(1000, 'ns')

        assert ct('10us') == np.timedelta64(10000, 'ns')
        assert ct('100us') == np.timedelta64(100000, 'ns')
        assert ct('1000us') == np.timedelta64(1000000, 'ns')
        assert ct('1000Us') == np.timedelta64(1000000, 'ns')
        assert ct('1000uS') == np.timedelta64(1000000, 'ns')

        assert ct('1ms') == np.timedelta64(1000000, 'ns')
        assert ct('10ms') == np.timedelta64(10000000, 'ns')
        assert ct('100ms') == np.timedelta64(100000000, 'ns')
        assert ct('1000ms') == np.timedelta64(1000000000, 'ns')

        assert ct('-1s') == -np.timedelta64(1000000000, 'ns')
        assert ct('1s') == np.timedelta64(1000000000, 'ns')
        assert ct('10s') == np.timedelta64(10000000000, 'ns')
        assert ct('100s') == np.timedelta64(100000000000, 'ns')
        assert ct('1000s') == np.timedelta64(1000000000000, 'ns')

        assert ct('1d') == conv(np.timedelta64(1, 'D'))
        assert ct('-1d') == -conv(np.timedelta64(1, 'D'))
        assert ct('1D') == conv(np.timedelta64(1, 'D'))
        assert ct('10D') == conv(np.timedelta64(10, 'D'))
        assert ct('100D') == conv(np.timedelta64(100, 'D'))
        assert ct('1000D') == conv(np.timedelta64(1000, 'D'))
        assert ct('10000D') == conv(np.timedelta64(10000, 'D'))

        # space
        assert ct(' 10000D ') == conv(np.timedelta64(10000, 'D'))
        assert ct(' - 10000D ') == -conv(np.timedelta64(10000, 'D'))

        # invalid
        pytest.raises(ValueError, ct, '1foo')
        pytest.raises(ValueError, ct, 'foo')

    def test_full_format_converters(self):
        def conv(v):
            return v.astype('m8[ns]')

        d1 = np.timedelta64(1, 'D')

        assert ct('1days') == conv(d1)
        assert ct('1days,') == conv(d1)
        assert ct('- 1days,') == -conv(d1)

        assert ct('00:00:01') == conv(np.timedelta64(1, 's'))
        assert ct('06:00:01') == conv(np.timedelta64(6 * 3600 + 1, 's'))
        assert ct('06:00:01.0') == conv(np.timedelta64(6 * 3600 + 1, 's'))
        assert ct('06:00:01.01') == conv(np.timedelta64(
            1000 * (6 * 3600 + 1) + 10, 'ms'))

        assert (ct('- 1days, 00:00:01') ==
                conv(-d1 + np.timedelta64(1, 's')))
        assert (ct('1days, 06:00:01') ==
                conv(d1 + np.timedelta64(6 * 3600 + 1, 's')))
        assert (ct('1days, 06:00:01.01') ==
                conv(d1 + np.timedelta64(1000 * (6 * 3600 + 1) + 10, 'ms')))

        # invalid
        pytest.raises(ValueError, ct, '- 1days, 00')

    def test_pickle(self):

        v = Timedelta('1 days 10:11:12.0123456')
        v_p = tm.round_trip_pickle(v)
        assert v == v_p

    def test_timedelta_hash_equality(self):
        # GH#11129
        v = Timedelta(1, 'D')
        td = timedelta(days=1)
        assert hash(v) == hash(td)

        d = {td: 2}
        assert d[v] == 2

        # python timedeltas drop ns resolution
        ns_td = Timedelta(1, 'ns')
        assert hash(ns_td) != hash(ns_td.to_pytimedelta())

    def test_implementation_limits(self):
        min_td = Timedelta(Timedelta.min)
        max_td = Timedelta(Timedelta.max)

        # GH 12727
        # timedelta limits correspond to int64 boundaries
        assert min_td.value == np.iinfo(np.int64).min + 1
        assert max_td.value == np.iinfo(np.int64).max

        # Beyond lower limit, a NAT before the Overflow
        assert (min_td - Timedelta(1, 'ns')) is NaT

        with pytest.raises(OverflowError):
            min_td - Timedelta(2, 'ns')

        with pytest.raises(OverflowError):
            max_td + Timedelta(1, 'ns')

        # Same tests using the internal nanosecond values
        td = Timedelta(min_td.value - 1, 'ns')
        assert td is NaT

        with pytest.raises(OverflowError):
            Timedelta(min_td.value - 2, 'ns')

        with pytest.raises(OverflowError):
            Timedelta(max_td.value + 1, 'ns')

    def test_total_seconds_precision(self):
        # GH#19458
        assert Timedelta('30S').total_seconds() == 30.0
        assert Timedelta('0').total_seconds() == 0.0
        assert Timedelta('-2S').total_seconds() == -2.0
        assert Timedelta('5.324S').total_seconds() == 5.324
        assert (Timedelta('30S').total_seconds() - 30.0) < 1e-20
        assert (30.0 - Timedelta('30S').total_seconds()) < 1e-20


@pytest.mark.parametrize('value, expected', [
    (Timedelta('10S'), True),
    (Timedelta('-10S'), True),
    (Timedelta(10, unit='ns'), True),
    (Timedelta(0, unit='ns'), False),
    (Timedelta(-10, unit='ns'), True),
    (Timedelta(None), True),
    (NaT, True),
])
def test_truthiness(value, expected):
    # https://github.com/pandas-dev/pandas/issues/21484
    assert bool(value) is expected
