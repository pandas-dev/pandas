from datetime import (
    UTC,
    datetime,
)
from zoneinfo import ZoneInfo

import pytest

from pandas.compat import IS64

import pandas as pd

from pandas.tseries.offsets import (
    BMonthBegin,
    BMonthEnd,
    BQuarterBegin,
    BQuarterEnd,
    BYearBegin,
    BYearEnd,
    MonthBegin,
    MonthEnd,
    QuarterBegin,
    QuarterEnd,
    YearBegin,
    YearEnd,
)


@pytest.fixture(
    params=[
        MonthBegin(-2),
        MonthBegin(-1),
        MonthBegin(1),
        MonthBegin(2),
        MonthEnd(-2),
        MonthEnd(-1),
        MonthEnd(1),
        MonthEnd(2),
        BMonthBegin(-2),
        BMonthBegin(-1),
        BMonthBegin(1),
        BMonthBegin(2),
        BMonthEnd(-2),
        BMonthEnd(-1),
        BMonthEnd(1),
        BMonthEnd(2),
        QuarterBegin(-2, startingMonth=1),
        QuarterBegin(-2, startingMonth=2),
        QuarterBegin(-2, startingMonth=3),
        QuarterBegin(-1, startingMonth=1),
        QuarterBegin(-1, startingMonth=2),
        QuarterBegin(-1, startingMonth=3),
        QuarterBegin(1, startingMonth=1),
        QuarterBegin(1, startingMonth=2),
        QuarterBegin(1, startingMonth=3),
        QuarterBegin(2, startingMonth=1),
        QuarterBegin(2, startingMonth=2),
        QuarterBegin(2, startingMonth=3),
        QuarterEnd(-2, startingMonth=1),
        QuarterEnd(-2, startingMonth=2),
        QuarterEnd(-2, startingMonth=3),
        QuarterEnd(-1, startingMonth=1),
        QuarterEnd(-1, startingMonth=2),
        QuarterEnd(-1, startingMonth=3),
        QuarterEnd(1, startingMonth=1),
        QuarterEnd(1, startingMonth=2),
        QuarterEnd(1, startingMonth=3),
        QuarterEnd(2, startingMonth=1),
        QuarterEnd(2, startingMonth=2),
        QuarterEnd(2, startingMonth=3),
        BQuarterBegin(-2, startingMonth=1),
        BQuarterBegin(-2, startingMonth=2),
        BQuarterBegin(-2, startingMonth=3),
        BQuarterBegin(-1, startingMonth=1),
        BQuarterBegin(-1, startingMonth=2),
        BQuarterBegin(-1, startingMonth=3),
        BQuarterBegin(1, startingMonth=1),
        BQuarterBegin(1, startingMonth=2),
        BQuarterBegin(1, startingMonth=3),
        BQuarterBegin(2, startingMonth=1),
        BQuarterBegin(2, startingMonth=2),
        BQuarterBegin(2, startingMonth=3),
        BQuarterEnd(-2, startingMonth=1),
        BQuarterEnd(-2, startingMonth=2),
        BQuarterEnd(-2, startingMonth=3),
        BQuarterEnd(-1, startingMonth=1),
        BQuarterEnd(-1, startingMonth=2),
        BQuarterEnd(-1, startingMonth=3),
        BQuarterEnd(1, startingMonth=1),
        BQuarterEnd(1, startingMonth=2),
        BQuarterEnd(1, startingMonth=3),
        BQuarterEnd(2, startingMonth=1),
        BQuarterEnd(2, startingMonth=2),
        BQuarterEnd(2, startingMonth=3),
        YearBegin(-2, month=1),
        YearBegin(-2, month=12),
        YearBegin(-1, month=1),
        YearBegin(-1, month=12),
        YearBegin(1, month=1),
        YearBegin(1, month=12),
        YearBegin(2, month=1),
        YearBegin(2, month=12),
        YearEnd(-2, month=1),
        YearEnd(-2, month=12),
        YearEnd(-1, month=1),
        YearEnd(-1, month=12),
        YearEnd(1, month=1),
        YearEnd(1, month=12),
        YearEnd(2, month=1),
        YearEnd(2, month=12),
        BYearBegin(-2, month=1),
        BYearBegin(-2, month=12),
        BYearBegin(-1, month=1),
        BYearBegin(-1, month=12),
        BYearBegin(1, month=1),
        BYearBegin(1, month=12),
        BYearBegin(2, month=1),
        BYearBegin(2, month=12),
        BYearEnd(-2, month=1),
        BYearEnd(-2, month=12),
        BYearEnd(-1, month=1),
        BYearEnd(-1, month=12),
        BYearEnd(1, month=1),
        BYearEnd(1, month=12),
        BYearEnd(2, month=1),
        BYearEnd(2, month=12),
    ]
)
def offset_cases(request):
    return request.param


@pytest.mark.parametrize(
    "dt",
    [
        datetime(1900, 1, 1),
        datetime(1900, 1, 1, tzinfo=UTC),
        pytest.param(
            datetime(1900, 1, 1, tzinfo=ZoneInfo("US/Eastern")),
            marks=pytest.mark.skipif(
                not IS64,
                reason=(
                    "stdlib datetime.fromtimestamp fails on 32-bit platforms with "
                    "overflow"
                ),
            ),
        ),
        datetime(1900, 1, 1, tzinfo=ZoneInfo("Africa/Kinshasa")),
    ],
)
def test_on_offset_implementations(dt, offset_cases):
    # check that the class-specific implementations of is_on_offset match
    # the general case definition:
    #   (dt + offset) - offset == dt
    compare = (dt + offset_cases) - offset_cases
    assert offset_cases.is_on_offset(dt) == (compare == dt)


def test_shift_across_dst(offset_cases):
    # GH#18319 check that 1) timezone is correctly normalized and
    # 2) that hour is not incorrectly changed by this normalization
    # Note that dti includes a transition across DST boundary
    dti = pd.date_range(
        start="2017-10-30 12:00:00", end="2017-11-06", freq="D", tz="US/Eastern"
    )
    assert (dti.hour == 12).all()  # we haven't screwed up yet

    res = dti + offset_cases
    assert (res.hour == 12).all()
