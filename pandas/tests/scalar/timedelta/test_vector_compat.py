"""
Tests for pd.Timedelta behavior that depend on the implementations
of Index/Series
"""
from datetime import timedelta

import numpy as np
import pytest

import pandas.util.testing as tm

from pandas._libs.tslib import iNaT
from pandas import (
    NaT, Timedelta,
    Series, DatetimeIndex, TimedeltaIndex,
    date_range, to_timedelta, timedelta_range, offsets
)


class TestTimedeltaArithmeticArrayCompat(object):
    def test_td_floordiv_numeric_series(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=4)
        ser = Series([1], dtype=np.int64)
        res = td // ser
        assert res.dtype.kind == 'm'

    def test_array_timedelta_floordiv(self):
        # GH#19761
        ints = date_range('2012-10-08', periods=4, freq='D').view('i8')
        msg = r"Use 'array // timedelta.value'"
        with tm.assert_produces_warning(FutureWarning) as m:
            result = ints // Timedelta(1, unit='s')

        assert msg in str(m[0].message)
        expected = np.array([1349654400, 1349740800, 1349827200, 1349913600],
                            dtype='i8')
        tm.assert_numpy_array_equal(result, expected)

    def test_td_rfloordiv_numeric_series(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=3)
        ser = Series([1], dtype=np.int64)
        res = td.__rfloordiv__(ser)
        assert res is NotImplemented
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            # TODO: GH#19761. Change to TypeError.
            ser // td

    @pytest.mark.parametrize('delta', [timedelta(days=1),
                                       Timedelta(1, unit='D')])
    def test_timedelta_arithmetic(self, delta):
        data = Series(['nat', '32 days'], dtype='timedelta64[ns]')

        result_method = data.add(delta)
        result_operator = data + delta
        expected = Series(['nat', '33 days'], dtype='timedelta64[ns]')
        tm.assert_series_equal(result_operator, expected)
        tm.assert_series_equal(result_method, expected)

        result_method = data.sub(delta)
        result_operator = data - delta
        expected = Series(['nat', '31 days'], dtype='timedelta64[ns]')
        tm.assert_series_equal(result_operator, expected)
        tm.assert_series_equal(result_method, expected)
        # GH#9396
        result_method = data.div(delta)
        result_operator = data / delta
        expected = Series([np.nan, 32.], dtype='float64')
        tm.assert_series_equal(result_operator, expected)
        tm.assert_series_equal(result_method, expected)


class TestTimedeltaArrayCompat(object):
    def test_overflow(self):
        # GH#9442
        s = Series(date_range('20130101', periods=100000, freq='H'))
        s[0] += Timedelta('1s 1ms')

        # mean
        result = (s - s.min()).mean()
        expected = Timedelta((DatetimeIndex((s - s.min())).asi8 / len(s)
                              ).sum())

        # the computation is converted to float so
        # might be some loss of precision
        assert np.allclose(result.value / 1000, expected.value / 1000)

        # sum
        with pytest.raises(ValueError):
            (s - s.min()).sum()

        s1 = s[0:10000]
        with pytest.raises(ValueError):
            (s1 - s1.min()).sum()

        s2 = s[0:1000]
        result = (s2 - s2.min()).sum()

    def test_apply_to_timedelta(self):
        timedelta_NaT = to_timedelta('NaT')

        list_of_valid_strings = ['00:00:01', '00:00:02']
        a = to_timedelta(list_of_valid_strings)
        b = Series(list_of_valid_strings).apply(to_timedelta)
        # Can't compare until apply on a Series gives the correct dtype
        # assert_series_equal(a, b)

        list_of_strings = ['00:00:01', np.nan, NaT, timedelta_NaT]

        # TODO: unused?
        a = to_timedelta(list_of_strings)  # noqa
        b = Series(list_of_strings).apply(to_timedelta)  # noqa
        # Can't compare until apply on a Series gives the correct dtype
        # assert_series_equal(a, b)

    def test_components(self):
        rng = timedelta_range('1 days, 10:11:12', periods=2, freq='s')
        rng.components

        # with nat
        s = Series(rng)
        s[1] = np.nan

        result = s.dt.components
        assert not result.iloc[0].isna().all()
        assert result.iloc[1].isna().all()

    def test_nat_converters(self):
        assert to_timedelta('nat', box=False).astype('int64') == iNaT
        assert to_timedelta('nan', box=False).astype('int64') == iNaT

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
            assert result == expected

        # validate all units
        # GH#6855
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

    def test_round(self):
        t1 = timedelta_range('1 days', periods=3, freq='1 min 2 s 3 us')
        t2 = -1 * t1
        t1a = timedelta_range('1 days', periods=3, freq='1 min 2 s')
        t1c = TimedeltaIndex([1, 1, 1], unit='D')

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
                                TimedeltaIndex([-1, -1, -1], unit='D')
                                )]:

            r1 = t1.round(freq)
            tm.assert_index_equal(r1, s1)
            r2 = t2.round(freq)
        tm.assert_index_equal(r2, s2)  # TODO: Should this be indented?

        # invalid
        for freq in ['Y', 'M', 'foobar']:
            with pytest.raises(ValueError):
                t1.round(freq)

    def test_timedelta_hash_equality(self):
        # GH#11129
        tds = timedelta_range('1 second', periods=20)
        assert all(hash(td) == hash(td.to_pytimedelta()) for td in tds)

    @pytest.mark.parametrize('naval', [NaT, None, float('nan'), np.nan])
    def test_contains(self, naval):
        # Checking for any NaT-like objects
        # GH#13603
        td = to_timedelta(range(5), unit='d') + offsets.Hour(1)
        assert not (naval in td)

        td = to_timedelta([NaT])
        assert (naval in td)
