"""
Hypothesis data generator helpers.
"""

from hypothesis import strategies as st
from hypothesis.extra.dateutil import timezones as dateutil_timezones

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

DATETIME_JAN_1_1900_OPTIONAL_TZ = st.datetimes(  # pyright: ignore[reportCallIssue]
    min_value=pd.Timestamp(1900, 1, 1).to_pydatetime(),  # pyright: ignore[reportArgumentType]
    max_value=pd.Timestamp(1900, 1, 1).to_pydatetime(),  # pyright: ignore[reportArgumentType]
    timezones=st.one_of(st.none(), dateutil_timezones(), st.timezones()),
)

DATETIME_IN_PD_TIMESTAMP_RANGE_NO_TZ = st.datetimes(
    min_value=pd.Timestamp.min.to_pydatetime(warn=False),
    max_value=pd.Timestamp.max.to_pydatetime(warn=False),
)

INT_NEG_999_TO_POS_999 = st.integers(-999, 999)

# The strategy for each type is registered in conftest.py, as they don't carry
# enough runtime information (e.g. type hints) to infer how to build them.
YQM_OFFSET = st.one_of(
    *map(
        st.from_type,
        [
            MonthBegin,
            MonthEnd,
            BMonthBegin,
            BMonthEnd,
            QuarterBegin,
            QuarterEnd,
            BQuarterBegin,
            BQuarterEnd,
            YearBegin,
            YearEnd,
            BYearBegin,
            BYearEnd,
        ],
    )
)
