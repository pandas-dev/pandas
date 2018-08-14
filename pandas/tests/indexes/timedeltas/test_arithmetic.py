# -*- coding: utf-8 -*-

import pytest
import numpy as np
from datetime import timedelta

import pandas as pd
import pandas.util.testing as tm
from pandas import (DatetimeIndex, TimedeltaIndex, Int64Index,
                    timedelta_range, date_range,
                    Series,
                    Timedelta)
from pandas.errors import NullFrequencyError


@pytest.fixture(params=[pd.offsets.Hour(2), timedelta(hours=2),
                        np.timedelta64(2, 'h'), Timedelta(hours=2)],
                ids=str)
def delta(request):
    # Several ways of representing two hours
    return request.param


@pytest.fixture(params=['B', 'D'])
def freq(request):
    return request.param


class TestTimedeltaIndexArithmetic(object):
    # Addition and Subtraction Operations

    # -------------------------------------------------------------
    # TimedeltaIndex.shift is used by __add__/__sub__

    def test_tdi_shift_empty(self):
        # GH#9903
        idx = pd.TimedeltaIndex([], name='xxx')
        tm.assert_index_equal(idx.shift(0, freq='H'), idx)
        tm.assert_index_equal(idx.shift(3, freq='H'), idx)

    def test_tdi_shift_hours(self):
        # GH#9903
        idx = pd.TimedeltaIndex(['5 hours', '6 hours', '9 hours'], name='xxx')
        tm.assert_index_equal(idx.shift(0, freq='H'), idx)
        exp = pd.TimedeltaIndex(['8 hours', '9 hours', '12 hours'], name='xxx')
        tm.assert_index_equal(idx.shift(3, freq='H'), exp)
        exp = pd.TimedeltaIndex(['2 hours', '3 hours', '6 hours'], name='xxx')
        tm.assert_index_equal(idx.shift(-3, freq='H'), exp)

    def test_tdi_shift_minutes(self):
        # GH#9903
        idx = pd.TimedeltaIndex(['5 hours', '6 hours', '9 hours'], name='xxx')
        tm.assert_index_equal(idx.shift(0, freq='T'), idx)
        exp = pd.TimedeltaIndex(['05:03:00', '06:03:00', '9:03:00'],
                                name='xxx')
        tm.assert_index_equal(idx.shift(3, freq='T'), exp)
        exp = pd.TimedeltaIndex(['04:57:00', '05:57:00', '8:57:00'],
                                name='xxx')
        tm.assert_index_equal(idx.shift(-3, freq='T'), exp)

    def test_tdi_shift_int(self):
        # GH#8083
        trange = pd.to_timedelta(range(5), unit='d') + pd.offsets.Hour(1)
        result = trange.shift(1)
        expected = TimedeltaIndex(['1 days 01:00:00', '2 days 01:00:00',
                                   '3 days 01:00:00',
                                   '4 days 01:00:00', '5 days 01:00:00'],
                                  freq='D')
        tm.assert_index_equal(result, expected)

    def test_tdi_shift_nonstandard_freq(self):
        # GH#8083
        trange = pd.to_timedelta(range(5), unit='d') + pd.offsets.Hour(1)
        result = trange.shift(3, freq='2D 1s')
        expected = TimedeltaIndex(['6 days 01:00:03', '7 days 01:00:03',
                                   '8 days 01:00:03', '9 days 01:00:03',
                                   '10 days 01:00:03'], freq='D')
        tm.assert_index_equal(result, expected)

    def test_shift_no_freq(self):
        # GH#19147
        tdi = TimedeltaIndex(['1 days 01:00:00', '2 days 01:00:00'], freq=None)
        with pytest.raises(NullFrequencyError):
            tdi.shift(2)

    # -------------------------------------------------------------
    # Binary operations TimedeltaIndex and integer

    def test_tdi_add_int(self, one):
        # Variants of `one` for #19012
        rng = timedelta_range('1 days 09:00:00', freq='H', periods=10)
        result = rng + one
        expected = timedelta_range('1 days 10:00:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)

    def test_tdi_iadd_int(self, one):
        rng = timedelta_range('1 days 09:00:00', freq='H', periods=10)
        expected = timedelta_range('1 days 10:00:00', freq='H', periods=10)
        rng += one
        tm.assert_index_equal(rng, expected)

    def test_tdi_sub_int(self, one):
        rng = timedelta_range('1 days 09:00:00', freq='H', periods=10)
        result = rng - one
        expected = timedelta_range('1 days 08:00:00', freq='H', periods=10)
        tm.assert_index_equal(result, expected)

    def test_tdi_isub_int(self, one):
        rng = timedelta_range('1 days 09:00:00', freq='H', periods=10)
        expected = timedelta_range('1 days 08:00:00', freq='H', periods=10)
        rng -= one
        tm.assert_index_equal(rng, expected)

    # -------------------------------------------------------------
    # __add__/__sub__ with integer arrays

    @pytest.mark.parametrize('box', [np.array, pd.Index])
    def test_tdi_add_integer_array(self, box):
        # GH#19959
        rng = timedelta_range('1 days 09:00:00', freq='H', periods=3)
        other = box([4, 3, 2])
        expected = TimedeltaIndex(['1 day 13:00:00'] * 3)
        result = rng + other
        tm.assert_index_equal(result, expected)
        result = other + rng
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize('box', [np.array, pd.Index])
    def test_tdi_sub_integer_array(self, box):
        # GH#19959
        rng = timedelta_range('9H', freq='H', periods=3)
        other = box([4, 3, 2])
        expected = TimedeltaIndex(['5H', '7H', '9H'])
        result = rng - other
        tm.assert_index_equal(result, expected)
        result = other - rng
        tm.assert_index_equal(result, -expected)

    @pytest.mark.parametrize('box', [np.array, pd.Index])
    def test_tdi_addsub_integer_array_no_freq(self, box):
        # GH#19959
        tdi = TimedeltaIndex(['1 Day', 'NaT', '3 Hours'])
        other = box([14, -1, 16])
        with pytest.raises(NullFrequencyError):
            tdi + other
        with pytest.raises(NullFrequencyError):
            other + tdi
        with pytest.raises(NullFrequencyError):
            tdi - other
        with pytest.raises(NullFrequencyError):
            other - tdi

    # -------------------------------------------------------------
    # Binary operations TimedeltaIndex and timedelta-like
    # Note: add and sub are tested in tests.test_arithmetic

    def test_tdi_iadd_timedeltalike(self, delta):
        # only test adding/sub offsets as + is now numeric
        rng = timedelta_range('1 days', '10 days')
        expected = timedelta_range('1 days 02:00:00', '10 days 02:00:00',
                                   freq='D')
        rng += delta
        tm.assert_index_equal(rng, expected)

    def test_tdi_isub_timedeltalike(self, delta):
        # only test adding/sub offsets as - is now numeric
        rng = timedelta_range('1 days', '10 days')
        expected = timedelta_range('0 days 22:00:00', '9 days 22:00:00')
        rng -= delta
        tm.assert_index_equal(rng, expected)

    # -------------------------------------------------------------

    def test_addition_ops(self):
        # with datetimes/timedelta and tdi/dti
        tdi = TimedeltaIndex(['1 days', pd.NaT, '2 days'], name='foo')

        # random indexes
        pytest.raises(NullFrequencyError, lambda: tdi + Int64Index([1, 2, 3]))

        # this is a union!
        # pytest.raises(TypeError, lambda : Int64Index([1,2,3]) + tdi)

    def test_tdi_ops_attributes(self):
        rng = timedelta_range('2 days', periods=5, freq='2D', name='x')

        result = rng + 1
        exp = timedelta_range('4 days', periods=5, freq='2D', name='x')
        tm.assert_index_equal(result, exp)
        assert result.freq == '2D'

        result = rng - 2
        exp = timedelta_range('-2 days', periods=5, freq='2D', name='x')
        tm.assert_index_equal(result, exp)
        assert result.freq == '2D'

        result = rng * 2
        exp = timedelta_range('4 days', periods=5, freq='4D', name='x')
        tm.assert_index_equal(result, exp)
        assert result.freq == '4D'

        result = rng / 2
        exp = timedelta_range('1 days', periods=5, freq='D', name='x')
        tm.assert_index_equal(result, exp)
        assert result.freq == 'D'

        result = -rng
        exp = timedelta_range('-2 days', periods=5, freq='-2D', name='x')
        tm.assert_index_equal(result, exp)
        assert result.freq == '-2D'

        rng = pd.timedelta_range('-2 days', periods=5, freq='D', name='x')

        result = abs(rng)
        exp = TimedeltaIndex(['2 days', '1 days', '0 days', '1 days',
                              '2 days'], name='x')
        tm.assert_index_equal(result, exp)
        assert result.freq is None

    # TODO: Needs more informative name, probably split up into
    # more targeted tests
    def test_timedelta(self, freq):
        index = date_range('1/1/2000', periods=50, freq=freq)

        shifted = index + timedelta(1)
        back = shifted + timedelta(-1)
        tm.assert_index_equal(index, back)

        if freq == 'D':
            expected = pd.tseries.offsets.Day(1)
            assert index.freq == expected
            assert shifted.freq == expected
            assert back.freq == expected
        else:  # freq == 'B'
            assert index.freq == pd.tseries.offsets.BusinessDay(1)
            assert shifted.freq is None
            assert back.freq == pd.tseries.offsets.BusinessDay(1)

        result = index - timedelta(1)
        expected = index + timedelta(-1)
        tm.assert_index_equal(result, expected)

        # GH4134, buggy with timedeltas
        rng = date_range('2013', '2014')
        s = Series(rng)
        result1 = rng - pd.offsets.Hour(1)
        result2 = DatetimeIndex(s - np.timedelta64(100000000))
        result3 = rng - np.timedelta64(100000000)
        result4 = DatetimeIndex(s - pd.offsets.Hour(1))
        tm.assert_index_equal(result1, result4)
        tm.assert_index_equal(result2, result3)
