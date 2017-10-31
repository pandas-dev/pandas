# -*- coding: utf-8 -*-
import warnings
from datetime import datetime, timedelta

import pytest

import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas.errors import PerformanceWarning
from pandas import (Timestamp, Timedelta, Series,
                    DatetimeIndex, TimedeltaIndex,
                    date_range)


class TestDatetimeIndexArithmetic(object):
    tz = [None, 'UTC', 'Asia/Tokyo', 'US/Eastern', 'dateutil/Asia/Singapore',
          'dateutil/US/Pacific']

    def test_add_iadd(self):
        for tz in self.tz:

            # offset
            offsets = [pd.offsets.Hour(2), timedelta(hours=2),
                       np.timedelta64(2, 'h'), Timedelta(hours=2)]

            for delta in offsets:
                rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
                result = rng + delta
                expected = pd.date_range('2000-01-01 02:00',
                                         '2000-02-01 02:00', tz=tz)
                tm.assert_index_equal(result, expected)
                rng += delta
                tm.assert_index_equal(rng, expected)

            # int
            rng = pd.date_range('2000-01-01 09:00', freq='H', periods=10,
                                tz=tz)
            result = rng + 1
            expected = pd.date_range('2000-01-01 10:00', freq='H', periods=10,
                                     tz=tz)
            tm.assert_index_equal(result, expected)
            rng += 1
            tm.assert_index_equal(rng, expected)

        idx = DatetimeIndex(['2011-01-01', '2011-01-02'])
        msg = "cannot add DatetimeIndex and Timestamp"
        with tm.assert_raises_regex(TypeError, msg):
            idx + Timestamp('2011-01-01')

        with tm.assert_raises_regex(TypeError, msg):
            Timestamp('2011-01-01') + idx

    def test_sub_isub(self):
        for tz in self.tz:

            # offset
            offsets = [pd.offsets.Hour(2), timedelta(hours=2),
                       np.timedelta64(2, 'h'), Timedelta(hours=2)]

            for delta in offsets:
                rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
                expected = pd.date_range('1999-12-31 22:00',
                                         '2000-01-31 22:00', tz=tz)

                result = rng - delta
                tm.assert_index_equal(result, expected)
                rng -= delta
                tm.assert_index_equal(rng, expected)

            # int
            rng = pd.date_range('2000-01-01 09:00', freq='H', periods=10,
                                tz=tz)
            result = rng - 1
            expected = pd.date_range('2000-01-01 08:00', freq='H', periods=10,
                                     tz=tz)
            tm.assert_index_equal(result, expected)
            rng -= 1
            tm.assert_index_equal(rng, expected)

    @pytest.mark.parametrize('addend', [
        datetime(2011, 1, 1),
        DatetimeIndex(['2011-01-01', '2011-01-02']),
        DatetimeIndex(['2011-01-01', '2011-01-02']).tz_localize('US/Eastern'),
        np.datetime64('2011-01-01'),
        Timestamp('2011-01-01')])
    def test_add_datetimelike_and_dti(self, addend):
        # GH#9631
        dti = DatetimeIndex(['2011-01-01', '2011-01-02'])
        msg = 'cannot add DatetimeIndex and {0}'.format(
            type(addend).__name__)
        with tm.assert_raises_regex(TypeError, msg):
            dti + addend
        with tm.assert_raises_regex(TypeError, msg):
            addend + dti

    @pytest.mark.parametrize('addend', [
        datetime(2011, 1, 1),
        DatetimeIndex(['2011-01-01', '2011-01-02']),
        DatetimeIndex(['2011-01-01', '2011-01-02']).tz_localize('US/Eastern'),
        np.datetime64('2011-01-01'),
        Timestamp('2011-01-01')])
    def test_add_datetimelike_and_dti_tz(self, addend):
        # GH#9631
        dti_tz = DatetimeIndex(['2011-01-01',
                                '2011-01-02']).tz_localize('US/Eastern')
        msg = 'cannot add DatetimeIndex and {0}'.format(
            type(addend).__name__)
        with tm.assert_raises_regex(TypeError, msg):
            dti_tz + addend
        with tm.assert_raises_regex(TypeError, msg):
            addend + dti_tz

    def test_sub_dti_dti(self):
        # previously performed setop (deprecated in 0.16.0), now changed to
        # return subtraction -> TimeDeltaIndex (GH ...)

        dti = date_range('20130101', periods=3)
        dti_tz = date_range('20130101', periods=3).tz_localize('US/Eastern')
        dti_tz2 = date_range('20130101', periods=3).tz_localize('UTC')
        expected = TimedeltaIndex([0, 0, 0])

        result = dti - dti
        tm.assert_index_equal(result, expected)

        result = dti_tz - dti_tz
        tm.assert_index_equal(result, expected)

        with pytest.raises(TypeError):
            dti_tz - dti

        with pytest.raises(TypeError):
            dti - dti_tz

        with pytest.raises(TypeError):
            dti_tz - dti_tz2

        # isub
        dti -= dti
        tm.assert_index_equal(dti, expected)

        # different length raises ValueError
        dti1 = date_range('20130101', periods=3)
        dti2 = date_range('20130101', periods=4)
        with pytest.raises(ValueError):
            dti1 - dti2

        # NaN propagation
        dti1 = DatetimeIndex(['2012-01-01', np.nan, '2012-01-03'])
        dti2 = DatetimeIndex(['2012-01-02', '2012-01-03', np.nan])
        expected = TimedeltaIndex(['1 days', np.nan, np.nan])
        result = dti2 - dti1
        tm.assert_index_equal(result, expected)

    def test_sub_period(self):
        # GH 13078
        # not supported, check TypeError
        p = pd.Period('2011-01-01', freq='D')

        for freq in [None, 'D']:
            idx = pd.DatetimeIndex(['2011-01-01', '2011-01-02'], freq=freq)

            with pytest.raises(TypeError):
                idx - p

            with pytest.raises(TypeError):
                p - idx

    def test_ufunc_coercions(self):
        idx = date_range('2011-01-01', periods=3, freq='2D', name='x')

        delta = np.timedelta64(1, 'D')
        for result in [idx + delta, np.add(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = date_range('2011-01-02', periods=3, freq='2D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == '2D'

        for result in [idx - delta, np.subtract(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = date_range('2010-12-31', periods=3, freq='2D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == '2D'

        delta = np.array([np.timedelta64(1, 'D'), np.timedelta64(2, 'D'),
                          np.timedelta64(3, 'D')])
        for result in [idx + delta, np.add(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = DatetimeIndex(['2011-01-02', '2011-01-05', '2011-01-08'],
                                freq='3D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == '3D'

        for result in [idx - delta, np.subtract(idx, delta)]:
            assert isinstance(result, DatetimeIndex)
            exp = DatetimeIndex(['2010-12-31', '2011-01-01', '2011-01-02'],
                                freq='D', name='x')
            tm.assert_index_equal(result, exp)
            assert result.freq == 'D'

    def test_overflow_offset(self):
        # xref https://github.com/statsmodels/statsmodels/issues/3374
        # ends up multiplying really large numbers which overflow

        t = Timestamp('2017-01-13 00:00:00', freq='D')
        offset = 20169940 * pd.offsets.Day(1)

        def f():
            t + offset
        pytest.raises(OverflowError, f)

        def f():
            offset + t
        pytest.raises(OverflowError, f)

        def f():
            t - offset
        pytest.raises(OverflowError, f)


# GH 10699
@pytest.mark.parametrize('klass,assert_func', zip([Series, DatetimeIndex],
                                                  [tm.assert_series_equal,
                                                   tm.assert_index_equal]))
def test_datetime64_with_DateOffset(klass, assert_func):
    s = klass(date_range('2000-01-01', '2000-01-31'), name='a')
    result = s + pd.DateOffset(years=1)
    result2 = pd.DateOffset(years=1) + s
    exp = klass(date_range('2001-01-01', '2001-01-31'), name='a')
    assert_func(result, exp)
    assert_func(result2, exp)

    result = s - pd.DateOffset(years=1)
    exp = klass(date_range('1999-01-01', '1999-01-31'), name='a')
    assert_func(result, exp)

    s = klass([Timestamp('2000-01-15 00:15:00', tz='US/Central'),
               pd.Timestamp('2000-02-15', tz='US/Central')], name='a')
    result = s + pd.offsets.Day()
    result2 = pd.offsets.Day() + s
    exp = klass([Timestamp('2000-01-16 00:15:00', tz='US/Central'),
                 Timestamp('2000-02-16', tz='US/Central')], name='a')
    assert_func(result, exp)
    assert_func(result2, exp)

    s = klass([Timestamp('2000-01-15 00:15:00', tz='US/Central'),
               pd.Timestamp('2000-02-15', tz='US/Central')], name='a')
    result = s + pd.offsets.MonthEnd()
    result2 = pd.offsets.MonthEnd() + s
    exp = klass([Timestamp('2000-01-31 00:15:00', tz='US/Central'),
                 Timestamp('2000-02-29', tz='US/Central')], name='a')
    assert_func(result, exp)
    assert_func(result2, exp)

    # array of offsets - valid for Series only
    if klass is Series:
        with tm.assert_produces_warning(PerformanceWarning):
            s = klass([Timestamp('2000-1-1'), Timestamp('2000-2-1')])
            result = s + Series([pd.offsets.DateOffset(years=1),
                                 pd.offsets.MonthEnd()])
            exp = klass([Timestamp('2001-1-1'), Timestamp('2000-2-29')
                         ])
            assert_func(result, exp)

            # same offset
            result = s + Series([pd.offsets.DateOffset(years=1),
                                 pd.offsets.DateOffset(years=1)])
            exp = klass([Timestamp('2001-1-1'), Timestamp('2001-2-1')])
            assert_func(result, exp)

    s = klass([Timestamp('2000-01-05 00:15:00'),
               Timestamp('2000-01-31 00:23:00'),
               Timestamp('2000-01-01'),
               Timestamp('2000-03-31'),
               Timestamp('2000-02-29'),
               Timestamp('2000-12-31'),
               Timestamp('2000-05-15'),
               Timestamp('2001-06-15')])

    # DateOffset relativedelta fastpath
    relative_kwargs = [('years', 2), ('months', 5), ('days', 3),
                       ('hours', 5), ('minutes', 10), ('seconds', 2),
                       ('microseconds', 5)]
    for i, kwd in enumerate(relative_kwargs):
        op = pd.DateOffset(**dict([kwd]))
        assert_func(klass([x + op for x in s]), s + op)
        assert_func(klass([x - op for x in s]), s - op)
        op = pd.DateOffset(**dict(relative_kwargs[:i + 1]))
        assert_func(klass([x + op for x in s]), s + op)
        assert_func(klass([x - op for x in s]), s - op)

    # assert these are equal on a piecewise basis
    offsets = ['YearBegin', ('YearBegin', {'month': 5}),
               'YearEnd', ('YearEnd', {'month': 5}),
               'MonthBegin', 'MonthEnd',
               'SemiMonthEnd', 'SemiMonthBegin',
               'Week', ('Week', {'weekday': 3}),
               'BusinessDay', 'BDay', 'QuarterEnd', 'QuarterBegin',
               'CustomBusinessDay', 'CDay', 'CBMonthEnd',
               'CBMonthBegin', 'BMonthBegin', 'BMonthEnd',
               'BusinessHour', 'BYearBegin', 'BYearEnd',
               'BQuarterBegin', ('LastWeekOfMonth', {'weekday': 2}),
               ('FY5253Quarter', {'qtr_with_extra_week': 1,
                                  'startingMonth': 1,
                                  'weekday': 2,
                                  'variation': 'nearest'}),
               ('FY5253', {'weekday': 0,
                           'startingMonth': 2,
                           'variation':
                           'nearest'}),
               ('WeekOfMonth', {'weekday': 2,
                                'week': 2}),
               'Easter', ('DateOffset', {'day': 4}),
               ('DateOffset', {'month': 5})]

    with warnings.catch_warnings(record=True):
        for normalize in (True, False):
            for do in offsets:
                if isinstance(do, tuple):
                    do, kwargs = do
                else:
                    do = do
                    kwargs = {}

                    for n in [0, 5]:
                        if (do in ['WeekOfMonth', 'LastWeekOfMonth',
                                   'FY5253Quarter', 'FY5253'] and n == 0):
                            continue
                    op = getattr(pd.offsets, do)(n,
                                                 normalize=normalize,
                                                 **kwargs)
                    assert_func(klass([x + op for x in s]), s + op)
                    assert_func(klass([x - op for x in s]), s - op)
                    assert_func(klass([op + x for x in s]), op + s)
