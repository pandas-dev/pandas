"""
Behavioral based tests for offsets and date_range.

This file is adapted from https://github.com/pandas-dev/pandas/pull/18761 -
which was more ambitious but less idiomatic in its use of Hypothesis.

You may wish to consult the previous version for inspiration on further
tests, or when trying to pin down the bugs exposed by the tests below.
"""

from datetime import (
    UTC,
    datetime,
)
import zoneinfo

import pytest

from pandas.compat import WASM

import pandas as pd

from pandas.tseries.offsets import (
    MonthBegin,
    MonthEnd,
    QuarterBegin,
    QuarterEnd,
    YearBegin,
    YearEnd,
)

YQM_OFFSETS = [
    MonthBegin(1),
    MonthEnd(1),
    QuarterBegin(1),
    QuarterEnd(1),
    YearBegin(1),
    YearEnd(1),
    MonthBegin(-1),
    QuarterEnd(-2),
]
DATETIMES = [
    datetime(1900, 1, 1),
    datetime(1900, 1, 1, tzinfo=UTC),
    datetime(1900, 1, 1, tzinfo=zoneinfo.ZoneInfo("Africa/Kinshasa")),
]

# ----------------------------------------------------------------
# Offset-specific behaviour tests


@pytest.mark.parametrize("dt", DATETIMES)
@pytest.mark.parametrize("offset", YQM_OFFSETS)
def test_on_offset_implementations(dt, offset):
    if offset.normalize:
        return
    # This case is flaky in CI 2024-11-04
    if (
        WASM
        and isinstance(dt.tzinfo, zoneinfo.ZoneInfo)
        and dt.tzinfo.key == "Indian/Cocos"
        and isinstance(offset, pd.offsets.MonthBegin)
    ):
        return
    # check that the class-specific implementations of is_on_offset match
    # the general case definition:
    #   (dt + offset) - offset == dt
    try:
        compare = (dt + offset) - offset
    except ValueError:
        # When dt + offset does not exist or is DST-ambiguous, skip
        # DST-ambiguous example (GH41906):
        # dt = datetime.datetime(1900, 1, 1, tzinfo=ZoneInfo('Africa/Kinshasa'))
        # offset = MonthBegin(66)
        return

    assert offset.is_on_offset(dt) == (compare == dt)


@pytest.mark.parametrize("offset", YQM_OFFSETS)
def test_shift_across_dst(offset):
    # GH#18319 check that 1) timezone is correctly normalized and
    # 2) that hour is not incorrectly changed by this normalization
    if offset.normalize:
        return

    # Note that dti includes a transition across DST boundary
    dti = pd.date_range(
        start="2017-10-30 12:00:00", end="2017-11-06", freq="D", tz="US/Eastern"
    )
    assert (dti.hour == 12).all()  # we haven't screwed up yet

    res = dti + offset
    assert (res.hour == 12).all()
