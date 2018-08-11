# -*- coding: utf-8 -*-
"""
Behavioral based tests for offsets and date_range.

This file is adapted from https://github.com/pandas-dev/pandas/pull/18761 -
which was more ambitious but less idiomatic in its use of Hypothesis.

You may wish to consult the previous version for inspiration on further
tests, or when trying to pin down the bugs exposed by the tests below.
"""

import pytest
from hypothesis import given, assume, strategies as st
from hypothesis.extra.pytz import timezones as pytz_timezones
from hypothesis.extra.dateutil import timezones as dateutil_timezones

import pandas as pd

from pandas.tseries.offsets import (
    Hour, Minute, Second, Milli, Micro, Nano,
    MonthEnd, MonthBegin, BMonthEnd, BMonthBegin,
    QuarterEnd, QuarterBegin, BQuarterEnd, BQuarterBegin,
    YearEnd, YearBegin, BYearEnd, BYearBegin,
)


tick_classes = [Hour, Minute, Second, Milli, Micro, Nano]
yqm_classes = [MonthBegin, MonthEnd, BMonthBegin, BMonthEnd,
               QuarterBegin, QuarterEnd, BQuarterBegin, BQuarterEnd,
               YearBegin, YearEnd, BYearBegin, BYearEnd]

# ----------------------------------------------------------------
# Helpers for generating random data

gen_date_range = st.builds(
    pd.date_range,
    start=st.datetimes(
        # TODO: Choose the min/max values more systematically
        min_value=pd.Timestamp(1900, 1, 1).to_pydatetime(),
        max_value=pd.Timestamp(2100, 1, 1).to_pydatetime()
    ),
    periods=st.integers(min_value=2, max_value=100),
    freq=st.sampled_from('Y Q M D H T s ms us ns'.split()),
    tz=st.one_of(st.none(), dateutil_timezones(), pytz_timezones()),
)

gen_random_datetime = st.datetimes(
    min_value=pd.Timestamp.min.to_pydatetime(),
    max_value=pd.Timestamp.max.to_pydatetime(),
    timezones=st.one_of(st.none(), dateutil_timezones(), pytz_timezones())
)

# Register the various offset classes so st.from_type can create instances.
# We *could* just append the strategies to a list, but this provides a nice
# demo and enables future tests to use a simple e.g. `from_type(Hour)`.
for cls in tick_classes + [MonthBegin, MonthEnd, BMonthBegin, BMonthEnd]:
    st.register_type_strategy(cls, st.builds(
        cls,
        n=st.integers(-99, 99),
        normalize=st.booleans(),
    ))

for cls in [YearBegin, YearEnd, BYearBegin, BYearEnd]:
    st.register_type_strategy(cls, st.builds(
        cls,
        n=st.integers(-5, 5),
        normalize=st.booleans(),
        month=st.integers(min_value=1, max_value=12),
    ))

for cls in [QuarterBegin, QuarterEnd, BQuarterBegin, BQuarterEnd]:
    st.register_type_strategy(cls, st.builds(
        cls,
        n=st.integers(-24, 24),
        normalize=st.booleans(),
        startingMonth=st.integers(min_value=1, max_value=12)
    ))

# This strategy can generate any kind of Offset in `tick_classes` or
# `yqm_classes`, with arguments as specified directly above in registration.
gen_yqm_offset = st.one_of([st.from_type(cls) for cls in yqm_classes])


# ----------------------------------------------------------------
# Offset-specific behaviour tests


# Based on CI runs: Always passes on OSX, fails on Linux, sometimes on Windows
@pytest.mark.xfail(strict=False, reason='inconsistent between OSs, Pythons')
@given(gen_random_datetime, gen_yqm_offset)
def test_on_offset_implementations(dt, offset):
    assume(not offset.normalize)
    # check that the class-specific implementations of onOffset match
    # the general case definition:
    #   (dt + offset) - offset == dt
    compare = (dt + offset) - offset
    assert offset.onOffset(dt) == (compare == dt)


@pytest.mark.xfail(strict=True)
@given(gen_yqm_offset, gen_date_range)
def test_apply_index_implementations(offset, rng):
    # offset.apply_index(dti)[i] should match dti[i] + offset
    assume(offset.n != 0)  # TODO: test for that case separately

    # rng = pd.date_range(start='1/1/2000', periods=100000, freq='T')
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


@pytest.mark.xfail(strict=True)
@given(gen_yqm_offset)
def test_shift_across_dst(offset):
    # GH#18319 check that 1) timezone is correctly normalized and
    # 2) that hour is not incorrectly changed by this normalization
    # Note that dti includes a transition across DST boundary
    dti = pd.date_range(start='2017-10-30 12:00:00', end='2017-11-06',
                        freq='D', tz='US/Eastern')
    assert (dti.hour == 12).all()  # we haven't screwed up yet

    res = dti + offset
    assert (res.hour == 12).all()
