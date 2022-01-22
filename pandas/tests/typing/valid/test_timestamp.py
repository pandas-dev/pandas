# flake8: noqa: F841
# TODO: many functions need return types annotations for pyright
# to run with reportGeneralTypeIssues = true

import datetime as dt

import pandas as pd


def test_types_init() -> None:
    ts: pd.Timestamp = pd.Timestamp("2021-03-01T12")
    ts1: pd.Timestamp = pd.Timestamp(dt.date(2021, 3, 15))
    ts2: pd.Timestamp = pd.Timestamp(dt.datetime(2021, 3, 10, 12))
    ts3: pd.Timestamp = pd.Timestamp(pd.Timestamp("2021-03-01T12"))
    ts4: pd.Timestamp = pd.Timestamp(1515590000.1, unit="s")
    ts5: pd.Timestamp = pd.Timestamp(1515590000.1, unit="s", tz="US/Pacific")
    ts6: pd.Timestamp = pd.Timestamp(1515590000100000000)  # plain integer (nanosecond)
    ts7: pd.Timestamp = pd.Timestamp(2021, 3, 10, 12)
    ts8: pd.Timestamp = pd.Timestamp(year=2021, month=3, day=10, hour=12)
    ts9: pd.Timestamp = pd.Timestamp(
        year=2021, month=3, day=10, hour=12, tz="US/Pacific"
    )


def test_types_arithmetic() -> None:
    # error: Incompatible types in assignment (expression has type "datetime", variable
    # has type "Timestamp")
    # error: Argument 1 to "to_datetime" has incompatible type "str"; expected
    # "datetime"
    ts: pd.Timestamp = pd.to_datetime("2021-03-01")  # type:ignore[assignment,arg-type]
    # error: Incompatible types in assignment (expression has type "datetime", variable
    # has type "Timestamp")
    # error: Argument 1 to "to_datetime" has incompatible type "str"; expected
    # "datetime"
    ts2: pd.Timestamp = pd.to_datetime("2021-01-01")  # type:ignore[assignment,arg-type]
    delta: pd.Timedelta = pd.to_timedelta("1 day")

    # error: Incompatible types in assignment (expression has type "timedelta", variable
    # has type "Timedelta")
    tsr: pd.Timedelta = ts - ts2  # type: ignore[assignment]
    tsr2: pd.Timestamp = ts + delta


def test_types_comparison() -> None:
    # Incompatible types in assignment (expression has type "datetime", variable has
    # type "Timestamp")
    # error: Argument 1 to "to_datetime" has incompatible type "str"; expected
    # "datetime"
    ts: pd.Timestamp = pd.to_datetime("2021-03-01")  # type: ignore[assignment,arg-type]
    # Incompatible types in assignment (expression has type "datetime", variable has
    # type "Timestamp")
    # error: Argument 1 to "to_datetime" has incompatible type "str"; expected
    # "datetime"
    ts2: pd.Timestamp = pd.to_datetime(  # type: ignore[assignment]
        "2021-01-01"  # type: ignore[arg-type]
    )

    tsr: bool = ts < ts2
    tsr2: bool = ts > ts2


def test_types_pydatetime() -> None:
    ts: pd.Timestamp = pd.Timestamp("2021-03-01T12")

    datet: dt.datetime = ts.to_pydatetime()
    datet2: dt.datetime = ts.to_pydatetime(False)
    datet3: dt.datetime = ts.to_pydatetime(warn=True)
