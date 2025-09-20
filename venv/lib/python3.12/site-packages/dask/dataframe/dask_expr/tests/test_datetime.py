from __future__ import annotations

import pytest

from dask.dataframe.dask_expr._collection import from_pandas
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

pd = _backend_library()


@pytest.fixture()
def ser():
    return pd.Series(pd.date_range(start="2020-01-01", end="2020-03-03"))


@pytest.fixture()
def dser(ser):
    return from_pandas(ser, npartitions=3)


@pytest.fixture()
def ser_td():
    return pd.Series(pd.timedelta_range(start="0s", end="1000s", freq="s"))


@pytest.fixture()
def dser_td(ser_td):
    return from_pandas(ser_td, npartitions=3)


@pytest.fixture()
def ser_pr():
    return pd.Series(pd.period_range(start="2020-01-01", end="2020-12-31", freq="D"))


@pytest.fixture()
def dser_pr(ser_pr):
    return from_pandas(ser_pr, npartitions=3)


@pytest.mark.parametrize(
    "func, args",
    [
        ("ceil", ("D",)),
        ("day_name", ()),
        ("floor", ("D",)),
        ("isocalendar", ()),
        ("month_name", ()),
        ("normalize", ()),
        ("round", ("D",)),
        ("strftime", ("%B %d, %Y, %r",)),
        ("to_period", ("D",)),
    ],
)
def test_datetime_accessor_methods(ser, dser, func, args):
    assert_eq(getattr(ser.dt, func)(*args), getattr(dser.dt, func)(*args))


@pytest.mark.parametrize(
    "func",
    [
        "date",
        "day",
        "day_of_week",
        "day_of_year",
        "dayofweek",
        "dayofyear",
        "days_in_month",
        "daysinmonth",
        "hour",
        "is_leap_year",
        "is_month_end",
        "is_month_start",
        "is_quarter_end",
        "is_quarter_start",
        "is_year_end",
        "is_year_start",
        "microsecond",
        "minute",
        "month",
        "nanosecond",
        "quarter",
        "second",
        "time",
        "timetz",
        "weekday",
        "year",
    ],
)
def test_datetime_accessor_properties(ser, dser, func):
    assert_eq(getattr(ser.dt, func), getattr(dser.dt, func))


@pytest.mark.parametrize(
    "func",
    [
        "components",
        "days",
        "microseconds",
        "nanoseconds",
        "seconds",
    ],
)
def test_timedelta_accessor_properties(ser_td, dser_td, func):
    assert_eq(getattr(ser_td.dt, func), getattr(dser_td.dt, func))


@pytest.mark.parametrize(
    "func",
    [
        "end_time",
        "start_time",
        "qyear",
        "week",
        "weekofyear",
    ],
)
def test_period_accessor_properties(ser_pr, dser_pr, func):
    assert_eq(getattr(ser_pr.dt, func), getattr(dser_pr.dt, func))
