"""
Hypothesis data generator helpers.
"""
from datetime import datetime

from hypothesis import strategies as st
from hypothesis.extra.dateutil import timezones as dateutil_timezones
from hypothesis.extra.pytz import timezones as pytz_timezones

from pandas.compat import is_platform_windows

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

OPTIONAL_INTS = st.lists(st.one_of(st.integers(), st.none()), max_size=10, min_size=3)

OPTIONAL_FLOATS = st.lists(st.one_of(st.floats(), st.none()), max_size=10, min_size=3)

OPTIONAL_TEXT = st.lists(st.one_of(st.none(), st.text()), max_size=10, min_size=3)

OPTIONAL_DICTS = st.lists(
    st.one_of(st.none(), st.dictionaries(st.text(), st.integers())),
    max_size=10,
    min_size=3,
)

OPTIONAL_LISTS = st.lists(
    st.one_of(st.none(), st.lists(st.text(), max_size=10, min_size=3)),
    max_size=10,
    min_size=3,
)

if is_platform_windows():
    _min_timestamp = datetime(1900, 1, 1)
else:
    _min_timestamp = None

DATETIME_NO_TZ = st.datetimes(min_value=_min_timestamp)

DATETIME_OPTIONAL_TZ = st.datetimes(
    min_value=_min_timestamp,
    timezones=st.one_of(st.none(), dateutil_timezones(), pytz_timezones()),
)

DATE_RANGE = st.builds(
    pd.date_range,
    start=st.datetimes(
        min_value=pd.Timestamp.min.to_pydatetime(warn=False),
        max_value=pd.Timestamp.max.to_pydatetime(warn=False),
    ),
    periods=st.integers(min_value=2, max_value=100),
    freq=st.sampled_from("Y Q M D H T s ms us ns".split()),
    tz=st.one_of(st.none(), dateutil_timezones(), pytz_timezones()),
)

INTEGER = st.integers()

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
