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


@pytest.fixture(params=[None, 'UTC', 'Asia/Tokyo',
                        'US/Eastern', 'dateutil/Asia/Singapore',
                        'dateutil/US/Pacific'])
def tz(request):
    return request.param


@pytest.fixture(params=[pd.offsets.Hour(2), timedelta(hours=2),
                        np.timedelta64(2, 'h'), Timedelta(hours=2)],
                ids=str)
def delta(request):
    # Several ways of representing two hours
    return request.param


@pytest.fixture(
    params=[
        datetime(2011, 1, 1),
        DatetimeIndex(['2011-01-01', '2011-01-02']),
        DatetimeIndex(['2011-01-01', '2011-01-02']).tz_localize('US/Eastern'),
        np.datetime64('2011-01-01'),
        Timestamp('2011-01-01')],
    ids=lambda x: type(x).__name__)
def addend(request):
    return request.param


class TestDatetimeIndexArithmetic(object):

    def test_dti_add_timestamp_raises(self):
        idx = DatetimeIndex(['2011-01-01', '2011-01-02'])
        msg = "cannot add DatetimeIndex and Timestamp"
        with tm.assert_raises_regex(TypeError, msg):
            idx + Timestamp('2011-01-01')

    def test_dti_radd_timestamp_raises(self):
        idx = DatetimeIndex(['2011-01-01', '2011-01-02'])
        msg = "cannot add DatetimeIndex and Timestamp"
        with tm.assert_raises_regex(TypeError, msg):
            Timestamp('2011-01-01') + idx

    # -------------------------------------------------------------
    # Binary operations DatetimeIndex and int

    def test_dti_add_int(self, tz, one):
        # Variants of `one` for #19012
        rng = pd.date_range('2000-01-01 09:00', freq='H',
                            periods=10, tz=tz)
        result = rng + one
        expected = pd.date_range('2000-01-01 10:00', freq='H',
                                 periods=10, tz=tz)
        tm.assert_index_equal(result, expected)

    def test_dti_iadd_int(self, tz, one):
        rng = pd.date_range('2000-01-01 09:00', freq='H',
                            periods=10, tz=tz)
        expected = pd.date_range('2000-01-01 10:00', freq='H',
                                 periods=10, tz=tz)
        rng += one
        tm.assert_index_equal(rng, expected)

    def test_dti_sub_int(self, tz, one):
        rng = pd.date_range('2000-01-01 09:00', freq='H',
                            periods=10, tz=tz)
        result = rng - one
        expected = pd.date_range('2000-01-01 08:00', freq='H',
                                 periods=10, tz=tz)
        tm.assert_index_equal(result, expected)

    def test_dti_isub_int(self, tz, one):
        rng = pd.date_range('2000-01-01 09:00', freq='H',
                            periods=10, tz=tz)
        expected = pd.date_range('2000-01-01 08:00', freq='H',
                                 periods=10, tz=tz)
        rng -= one
        tm.assert_index_equal(rng, expected)

    # -------------------------------------------------------------
    # Binary operations DatetimeIndex and timedelta-like

    def test_dti_add_timedeltalike(self, tz, delta):
        rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
        result = rng + delta
        expected = pd.date_range('2000-01-01 02:00',
                                 '2000-02-01 02:00', tz=tz)
        tm.assert_index_equal(result, expected)

    def test_dti_iadd_timedeltalike(self, tz, delta):
        rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
        expected = pd.date_range('2000-01-01 02:00',
                                 '2000-02-01 02:00', tz=tz)
        rng += delta
        tm.assert_index_equal(rng, expected)

    def test_dti_sub_timedeltalike(self, tz, delta):
        rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
        expected = pd.date_range('1999-12-31 22:00',
                                 '2000-01-31 22:00', tz=tz)
        result = rng - delta
        tm.assert_index_equal(result, expected)

    def test_dti_isub_timedeltalike(self, tz, delta):
        rng = pd.date_range('2000-01-01', '2000-02-01', tz=tz)
        expected = pd.date_range('1999-12-31 22:00',
                                 '2000-01-31 22:00', tz=tz)
        rng -= delta
        tm.assert_index_equal(rng, expected)

    # -------------------------------------------------------------
    # Binary operations DatetimeIndex and TimedeltaIndex/array
    def test_dti_add_tdi(self, tz):
        # GH 17558
        dti = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        tdi = pd.timedelta_range('0 days', periods=10)
        expected = pd.date_range('2017-01-01', periods=10, tz=tz)

        # add with TimdeltaIndex
        result = dti + tdi
        tm.assert_index_equal(result, expected)

        result = tdi + dti
        tm.assert_index_equal(result, expected)

        # add with timedelta64 array
        result = dti + tdi.values
        tm.assert_index_equal(result, expected)

        result = tdi.values + dti
        tm.assert_index_equal(result, expected)

    def test_dti_iadd_tdi(self, tz):
        # GH 17558
        dti = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        tdi = pd.timedelta_range('0 days', periods=10)
        expected = pd.date_range('2017-01-01', periods=10, tz=tz)

        # iadd with TimdeltaIndex
        result = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        result += tdi
        tm.assert_index_equal(result, expected)

        result = pd.timedelta_range('0 days', periods=10)
        result += dti
        tm.assert_index_equal(result, expected)

        # iadd with timedelta64 array
        result = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        result += tdi.values
        tm.assert_index_equal(result, expected)

        result = pd.timedelta_range('0 days', periods=10)
        result += dti
        tm.assert_index_equal(result, expected)

    def test_dti_sub_tdi(self, tz):
        # GH 17558
        dti = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        tdi = pd.timedelta_range('0 days', periods=10)
        expected = pd.date_range('2017-01-01', periods=10, tz=tz, freq='-1D')

        # sub with TimedeltaIndex
        result = dti - tdi
        tm.assert_index_equal(result, expected)

        msg = 'cannot subtract TimedeltaIndex and DatetimeIndex'
        with tm.assert_raises_regex(TypeError, msg):
            tdi - dti

        # sub with timedelta64 array
        result = dti - tdi.values
        tm.assert_index_equal(result, expected)

        msg = 'cannot perform __neg__ with this index type:'
        with tm.assert_raises_regex(TypeError, msg):
            tdi.values - dti

    def test_dti_isub_tdi(self, tz):
        # GH 17558
        dti = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        tdi = pd.timedelta_range('0 days', periods=10)
        expected = pd.date_range('2017-01-01', periods=10, tz=tz, freq='-1D')

        # isub with TimedeltaIndex
        result = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        result -= tdi
        tm.assert_index_equal(result, expected)

        msg = 'cannot subtract TimedeltaIndex and DatetimeIndex'
        with tm.assert_raises_regex(TypeError, msg):
            tdi -= dti

        # isub with timedelta64 array
        result = DatetimeIndex([Timestamp('2017-01-01', tz=tz)] * 10)
        result -= tdi.values
        tm.assert_index_equal(result, expected)

        msg = '|'.join(['cannot perform __neg__ with this index type:',
                        'ufunc subtract cannot use operands with types'])
        with tm.assert_raises_regex(TypeError, msg):
            tdi.values -= dti

    # -------------------------------------------------------------
    # Binary Operations DatetimeIndex and datetime-like
    # TODO: A couple other tests belong in this section.  Move them in
    # A PR where there isn't already a giant diff.

    def test_add_datetimelike_and_dti(self, addend):
        # GH#9631
        dti = DatetimeIndex(['2011-01-01', '2011-01-02'])
        msg = 'cannot add DatetimeIndex and {0}'.format(
            type(addend).__name__)
        with tm.assert_raises_regex(TypeError, msg):
            dti + addend
        with tm.assert_raises_regex(TypeError, msg):
            addend + dti

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

    # -------------------------------------------------------------

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

    def test_datetimeindex_sub_timestamp_overflow(self):
        dtimax = pd.to_datetime(['now', pd.Timestamp.max])
        dtimin = pd.to_datetime(['now', pd.Timestamp.min])

        tsneg = Timestamp('1950-01-01')
        ts_neg_variants = [tsneg,
                           tsneg.to_pydatetime(),
                           tsneg.to_datetime64().astype('datetime64[ns]'),
                           tsneg.to_datetime64().astype('datetime64[D]')]

        tspos = Timestamp('1980-01-01')
        ts_pos_variants = [tspos,
                           tspos.to_pydatetime(),
                           tspos.to_datetime64().astype('datetime64[ns]'),
                           tspos.to_datetime64().astype('datetime64[D]')]

        for variant in ts_neg_variants:
            with pytest.raises(OverflowError):
                dtimax - variant

        expected = pd.Timestamp.max.value - tspos.value
        for variant in ts_pos_variants:
            res = dtimax - variant
            assert res[1].value == expected

        expected = pd.Timestamp.min.value - tsneg.value
        for variant in ts_neg_variants:
            res = dtimin - variant
            assert res[1].value == expected

        for variant in ts_pos_variants:
            with pytest.raises(OverflowError):
                dtimin - variant

    @pytest.mark.parametrize('box', [np.array, pd.Index])
    def test_dti_add_offset_array(self, tz, box):
        # GH#18849
        dti = pd.date_range('2017-01-01', periods=2, tz=tz)
        other = box([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti + other
        expected = DatetimeIndex([dti[n] + other[n] for n in range(len(dti))],
                                 name=dti.name, freq='infer')
        tm.assert_index_equal(res, expected)

        with tm.assert_produces_warning(PerformanceWarning):
            res2 = other + dti
        tm.assert_index_equal(res2, expected)

    @pytest.mark.parametrize('box', [np.array, pd.Index])
    def test_dti_sub_offset_array(self, tz, box):
        # GH#18824
        dti = pd.date_range('2017-01-01', periods=2, tz=tz)
        other = box([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti - other
        expected = DatetimeIndex([dti[n] - other[n] for n in range(len(dti))],
                                 name=dti.name, freq='infer')
        tm.assert_index_equal(res, expected)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('foo', 'foo', 'foo')])
    def test_dti_with_offset_series(self, tz, names):
        # GH#18849
        dti = pd.date_range('2017-01-01', periods=2, tz=tz, name=names[0])
        other = Series([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)],
                       name=names[1])

        expected_add = Series([dti[n] + other[n] for n in range(len(dti))],
                              name=names[2])

        with tm.assert_produces_warning(PerformanceWarning):
            res = dti + other
        tm.assert_series_equal(res, expected_add)

        with tm.assert_produces_warning(PerformanceWarning):
            res2 = other + dti
        tm.assert_series_equal(res2, expected_add)

        expected_sub = Series([dti[n] - other[n] for n in range(len(dti))],
                              name=names[2])

        with tm.assert_produces_warning(PerformanceWarning):
            res3 = dti - other
        tm.assert_series_equal(res3, expected_sub)


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
