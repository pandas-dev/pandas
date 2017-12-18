# -*- coding: utf-8 -*-
"""
Behavioral based tests for offsets and date_range.
"""
from datetime import timedelta

import pytest
from hypothesis import given, assume
import hypothesis.strategies as st
import hypothesis.extra.numpy as hen
import hypothesis.extra.pytz as hepytz  # hypothesis[pytz]

import pandas as pd

from pandas.tseries.offsets import (Hour, Minute, Second, Milli, Micro, Nano,
                                    MonthEnd, MonthBegin,
                                    BMonthEnd, BMonthBegin,
                                    QuarterEnd, QuarterBegin,
                                    BQuarterEnd, BQuarterBegin,
                                    YearEnd, YearBegin,
                                    BYearEnd, BYearBegin,
                                    Week, LastWeekOfMonth, WeekOfMonth,
                                    SemiMonthBegin, SemiMonthEnd,
                                    Easter,
                                    FY5253, FY5253Quarter,
                                    DateOffset)
# TODO:
# BusinessDay, BusinessHour, CustomBusinessDay, CustomBusinessHour,
# CustomBusinessMonthEnd, CustomBusinessMonthBegin


tick_classes = [Hour, Minute, Second, Milli, Micro, Nano]
yqm_classes = [MonthBegin, MonthEnd, BMonthBegin, BMonthEnd,
               QuarterBegin, QuarterEnd, BQuarterBegin, BQuarterEnd,
               YearBegin, YearEnd, BYearBegin, BYearEnd]
offset_types = [Week, LastWeekOfMonth, WeekOfMonth, SemiMonthEnd,
                SemiMonthBegin, FY5253Quarter, FY5253,
                Easter, DateOffset] + tick_classes + yqm_classes

# ----------------------------------------------------------------
# Helpers for generating random data

dt_max = pd.Timestamp.max.replace(nanosecond=0).to_pydatetime()
td_max = timedelta(106751, 85636, 854775)
td_min = -td_max - timedelta(microseconds=1)

n_strategy = st.integers(min_value=-999, max_value=999)
# TODO: Choose these bounds systematically.  (-999, 999) is arbitrarily chosen
# to get rid of OverflowErrors in development
month_strategy = st.integers(min_value=1, max_value=12)
weekday_strategy = st.integers(min_value=0, max_value=6)


def gen_dst_crossing():
    # Generate either a pair of Timestamps or a date_range that is known
    # to cross a DST transition
    raise NotImplementedError


def gen_date_range_freq():
    # return a freq str or offset object suitable for passing as
    # `freq` kwarg to date_range
    return st.sampled_from(['Y', 'Q', 'M', 'D', 'H',
                            'T', 's', 'ms', 'us', 'ns'])
    # TODO: Add the rest; business, multiples, ...


@st.composite
def gen_random_date_range(draw):
    # TODO: Choose the min/max values more systematically
    start = st.datetimes(min_value=pd.Timestamp(1900, 1, 1).to_pydatetime(),
                         max_value=pd.Timestamp(2100, 1, 1).to_pydatetime())
    periods = st.integers(min_value=10, max_value=100)
    freq = gen_date_range_freq()
    tz = gen_random_tz()

    dti = pd.date_range(start=draw(start), tz=draw(tz),
                        freq=draw(freq), periods=draw(periods))
    return dti


def gen_random_tz():
    # Allows None
    return st.one_of(st.none(), hepytz.timezones())
    # TODO: Weighting between naive and timezones?
    # TODO: Get datetuil timezones?


gen_random_datetime = st.datetimes(min_value=pd.Timestamp.min.to_pydatetime(),
                                   max_value=pd.Timestamp.max.to_pydatetime(),
                                   timezones=gen_random_tz())


def gen_random_timestamp():
    nano = st.integers(min_value=0, max_value=999)
    dt = st.datetimes(min_value=pd.Timestamp.min.to_pydatetime(),
                      max_value=pd.Timestamp.max.to_pydatetime(),
                      timezones=gen_random_tz())
    ts = pd.Timestamp(dt)

    if dt != dt_max:
        ts.replace(nanosecond=nano)
    else:
        ts = ts.replace(nanosecond=min(nano, pd.Timestamp.max.nanosecond))

    # TODO: worry about timezones near min/max?
    return ts


def gen_random_datelike():
    # py_dates = st.dates()
    py_datetimes = gen_random_datetime

    # dt64_dtypes = hen.datetime64_dtypes()
    # np_dates = hen.arrays(dtype=dt64_dtypes, shape=())
    # TODO: Allow for non-scalar versions?
    # FIXME: dt64.__add__(offset) does not get dispatched to
    # offset.__radd__(dt64), just raises TypeError

    any_dates = st.one_of(py_datetimes)
    return any_dates


def gen_timedeltalike():
    py_timedeltas = st.timedeltas(min_value=td_min, max_value=td_max)
    pd_timedeltas = py_timedeltas.map(pd.Timedelta)
    # TODO: get those last few nanoseconds?

    td64_dtypes = hen.timedelta64_dtypes()
    np_timedeltas = hen.arrays(dtype=td64_dtypes, shape=())
    # TODO: Allow for non-scalar versions?

    # TODO: Week
    # TODO: Tick
    any_tds = st.one_of(py_timedeltas, pd_timedeltas, np_timedeltas)
    return any_tds


@st.composite
def gen_random_relativedelta_DateOffset(draw):
    relativedelta_kwds = set([
        'years', 'months', 'weeks', 'days',
        'year', 'month', 'week', 'day', 'weekday',
        'hour', 'minute', 'second', 'microsecond',
        'nanosecond', 'nanoseconds',
        'hours', 'minutes', 'seconds', 'milliseconds', 'microseconds'])
    kwargs = {kwd: st.integers() for kwd in relativedelta_kwds}
    kwargs['n'] = st.integers()
    kwargs['normalize'] = st.booleans()
    kwargs = {key: draw(kwargs[key]) for key in kwargs}
    return DateOffset(**kwargs)


@st.composite
def gen_random_offset(draw, cls):
    # Note: `draw` is a dummy argument that gets supplied by the composite
    # decorator
    n = n_strategy
    normalize = st.booleans()

    if cls in tick_classes + [MonthBegin, MonthEnd, BMonthBegin, BMonthEnd,
                              Easter]:
        n = n.filter(lambda x: abs(x) < 100)  # TODO: avoid arbitrary cutoff
        tup = st.tuples(n, normalize)

    elif cls in [QuarterBegin, QuarterEnd, BQuarterBegin, BQuarterEnd]:
        n = n.filter(lambda x: abs(x) < 25)  # TODO: avoid arbitrary cutoff
        startingMonth = month_strategy
        tup = st.tuples(n, normalize, startingMonth)

    elif cls in [YearBegin, YearEnd, BYearBegin, BYearEnd]:
        n = n.filter(lambda x: abs(x) < 6)  # TODO: avoid arbitrary cutoff
        month = month_strategy
        tup = st.tuples(n, normalize, month)

    elif cls == Week:
        n = n.filter(lambda x: abs(x) < 400)  # TODO: avoid arbitrary cutoff
        weekday = st.sampled_from([None, 0, 1, 2, 3, 4, 5, 6])
        tup = st.tuples(n, normalize, weekday)

    elif cls == LastWeekOfMonth:
        n = n.filter(lambda x: abs(x) < 400)  # TODO: avoid arbitrary cutoff
        n = n.filter(lambda x: x != 0)
        weekday = weekday_strategy
        tup = st.tuples(n, normalize, weekday)

    elif cls == WeekOfMonth:
        n = n.filter(lambda x: abs(x) < 400)  # TODO: avoid arbitrary cutoff
        n = n.filter(lambda x: x != 0)
        week = st.integers(min_value=0, max_value=3)
        weekday = weekday_strategy
        tup = st.tuples(n, normalize, week, weekday)

    elif cls in [SemiMonthBegin, SemiMonthEnd]:
        n = n.filter(lambda x: abs(x) < 800)  # TODO: avoid arbitrary cutoff
        day_of_month = st.integers(min_value=cls._min_day_of_month,
                                   max_value=27)
        tup = st.tuples(n, normalize, day_of_month)

    elif cls is FY5253:
        n = n.filter(lambda x: abs(x) < 6)  # TODO: avoid arbitrary cutoff
        n = n.filter(lambda x: x != 0)
        weekday = weekday_strategy
        startingMonth = month_strategy
        variation = st.sampled_from(["nearest", "last"])
        tup = st.tuples(n, normalize, weekday, startingMonth, variation)

    elif cls is FY5253Quarter:
        n = n.filter(lambda x: abs(x) < 24)  # TODO: avoid arbitrary cutoff
        n = n.filter(lambda x: x != 0)
        weekday = weekday_strategy
        startingMonth = month_strategy
        qtr_with_extra_week = st.integers(min_value=1, max_value=4)
        variation = st.sampled_from(["nearest", "last"])
        tup = st.tuples(n, normalize, weekday, startingMonth,
                        qtr_with_extra_week, variation)

    elif cls is DateOffset:
        # klass = cls(days=value, normalize=normalize)
        return gen_random_relativedelta_DateOffset()

    else:
        raise NotImplementedError(cls)

    args = draw(tup)
    return cls(*args)

# ----------------------------------------------------------------
# Tick-specific behavior tests


@given(n=n_strategy, m=n_strategy)
@pytest.mark.parametrize('cls', tick_classes)
def test_tick_add_sub(cls, n, m):
    # For all Tick subclasses and all integers n, m, we should have
    # tick(n) + tick(m) == tick(n+m)
    # tick(n) - tick(m) == tick(n-m)
    left = cls(n)
    right = cls(m)
    expected = cls(n + m)

    assert left + right == expected
    assert left.apply(right) == expected

    expected = cls(n - m)
    assert left - right == expected


@given(n=n_strategy, m=n_strategy)
@pytest.mark.parametrize('cls', tick_classes)
def test_tick_equality(cls, n, m):
    # tick == tock iff tick.n == tock.n
    left = cls(n)
    right = cls(m)
    if n == m:
        assert left == right
        assert not (left != right)
    else:
        assert left != right
        assert not (left == right)


# ----------------------------------------------------------------

@given(dt=gen_random_datelike(), data=st.data())
@pytest.mark.parametrize('cls', offset_types)
def test_on_offset_implementations(cls, dt, data):
    # check that the class-specific implementations of onOffset match
    # the general case definition:
    #   (dt + offset) - offset == dt

    offset = data.draw(gen_random_offset(cls), label='offset')
    # TODO: Is there a more performant way to do this?

    assume(not offset.normalize)
    compare = (dt + offset) - offset
    expected = compare == dt

    res = offset.onOffset(dt)
    assert res == expected


@given(data=st.data())
@pytest.mark.parametrize('cls', yqm_classes)
def test_apply_index_implementations(cls, data):
    # offset.apply_index(dti)[i] should match dti[i] + offset

    offset = data.draw(gen_random_offset(cls), label='offset')
    assume(offset.n != 0)  # TODO: test for that case separately

    # rng = pd.date_range(start='1/1/2000', periods=100000, freq='T')
    rng = data.draw(gen_random_date_range(), label='rng')
    ser = pd.Series(rng)

    res = rng + offset
    res_v2 = offset.apply_index(rng)
    assert (res == res_v2).all()

    assert res[0] == rng[0] + offset
    assert res[-1] == rng[-1] + offset
    res2 = ser + offset
    # apply_index is only for indexes, not series, so no res2_v2
    assert res2.iloc[0] == ser.iloc[0] + offset
    assert res2.iloc[-1] == ser.iloc[-1] + offset
    # TODO: Check randomly assorted entries, not just first/last


@given(freq=gen_date_range_freq())
def test_range_matches_addition(freq):

    raise pytest.skip('Need to generate date_range args')
    dr = pd.date_range('2016-10-30 12:00:00', freq=freq,
                       periods=20, tz='US/Eastern')
    assert dr[-1] > pd.Timestamp('2016-11-10')  # DST transition is crossed

    res = dr + freq
    assert res[:-1].equals(dr[1:])


@given(data=st.data())
@pytest.mark.parametrize('cls', yqm_classes)
def test_shift_across_dst(cls, data):
    # GH#18319 check that 1) timezone is correctly normalized and
    # 2) that hour is not incorrectly changed by this normalization

    raise pytest.skip('Need to generate date_range args')
    offset = data.draw(gen_random_offset(cls), label='offset')
    dti = pd.date_range(start='2017-10-30 12:00:00', end='2017-11-06',
                        freq='D', tz='US/Eastern')
    # dti includes a transition across DST boundary
    assert (dti.hour == 12).all()  # we haven't screwed up yet

    res = dti + offset
    assert (res.hour == 12).all()
