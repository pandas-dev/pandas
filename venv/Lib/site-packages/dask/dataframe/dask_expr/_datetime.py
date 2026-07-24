from __future__ import annotations

from dask.dataframe.dask_expr._accessor import Accessor


class DatetimeAccessor(Accessor):
    """Accessor object for datetimelike properties of the Series values.

    Examples
    --------

    >>> s.dt.microsecond  # doctest: +SKIP
    """

    _accessor_name = "dt"

    _accessor_methods = (
        "ceil",
        "day_name",
        "floor",
        "isocalendar",
        "month_name",
        "normalize",
        "round",
        "strftime",
        "to_period",
        "to_pydatetime",
        "to_pytimedelta",
        "to_timestamp",
        "total_seconds",
        "tz_convert",
        "tz_localize",
    )

    _accessor_properties = (
        "components",
        "date",
        "day",
        "day_of_week",
        "day_of_year",
        "dayofweek",
        "dayofyear",
        "days",
        "days_in_month",
        "daysinmonth",
        "end_time",
        "freq",
        "hour",
        "is_leap_year",
        "is_month_end",
        "is_month_start",
        "is_quarter_end",
        "is_quarter_start",
        "is_year_end",
        "is_year_start",
        "microsecond",
        "microseconds",
        "minute",
        "month",
        "nanosecond",
        "nanoseconds",
        "quarter",
        "qyear",
        "second",
        "seconds",
        "start_time",
        "time",
        "timetz",
        "tz",
        "week",
        "weekday",
        "weekofyear",
        "year",
    )
