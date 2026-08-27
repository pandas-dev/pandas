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

_ANCHORS = [
    *((klass, {}) for klass in (MonthBegin, MonthEnd, BMonthBegin, BMonthEnd)),
    *(
        (klass, {"startingMonth": month})
        for klass in (QuarterBegin, QuarterEnd, BQuarterBegin, BQuarterEnd)
        for month in (1, 2, 3)
    ),
    *(
        (klass, {"month": month})
        for klass in (YearBegin, YearEnd, BYearBegin, BYearEnd)
        for month in (1, 12)
    ),
]


@pytest.fixture(
    params=[klass(n, **kwargs) for klass, kwargs in _ANCHORS for n in (-2, -1, 1, 2)],  # type: ignore[arg-type]
    ids=lambda off: off.freqstr,
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
