"""Helpful functions used internally within arrow."""

import datetime
from typing import Any, Optional, cast

from dateutil.rrule import WEEKLY, rrule

from arrow.constants import (
    MAX_ORDINAL,
    MAX_TIMESTAMP,
    MAX_TIMESTAMP_MS,
    MAX_TIMESTAMP_US,
    MIN_ORDINAL,
)


def next_weekday(
    start_date: Optional[datetime.date], weekday: int
) -> datetime.datetime:
    """Get next weekday from the specified start date.

    :param start_date: Datetime object representing the start date.
    :param weekday: Next weekday to obtain. Can be a value between 0 (Monday) and 6 (Sunday).
    :return: Datetime object corresponding to the next weekday after start_date.

    Usage::

        # Get first Monday after epoch
        >>> next_weekday(datetime(1970, 1, 1), 0)
        1970-01-05 00:00:00

        # Get first Thursday after epoch
        >>> next_weekday(datetime(1970, 1, 1), 3)
        1970-01-01 00:00:00

        # Get first Sunday after epoch
        >>> next_weekday(datetime(1970, 1, 1), 6)
        1970-01-04 00:00:00
    """
    if weekday < 0 or weekday > 6:
        raise ValueError("Weekday must be between 0 (Monday) and 6 (Sunday).")
    return cast(
        datetime.datetime,
        rrule(freq=WEEKLY, dtstart=start_date, byweekday=weekday, count=1)[0],
    )


def is_timestamp(value: Any) -> bool:
    """Check if value is a valid timestamp."""
    if isinstance(value, bool):
        return False
    if not isinstance(value, (int, float, str)):
        return False
    try:
        float(value)
        return True
    except ValueError:
        return False


def validate_ordinal(value: Any) -> None:
    """Raise an exception if value is an invalid Gregorian ordinal.

    :param value: the input to be checked

    """
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"Ordinal must be an integer (got type {type(value)}).")
    if not (MIN_ORDINAL <= value <= MAX_ORDINAL):
        raise ValueError(f"Ordinal {value} is out of range.")


def normalize_timestamp(timestamp: float) -> float:
    """Normalize millisecond and microsecond timestamps into normal timestamps."""
    if timestamp > MAX_TIMESTAMP:
        if timestamp < MAX_TIMESTAMP_MS:
            timestamp /= 1000
        elif timestamp < MAX_TIMESTAMP_US:
            timestamp /= 1_000_000
        else:
            raise ValueError(f"The specified timestamp {timestamp!r} is too large.")
    return timestamp


# Credit to https://stackoverflow.com/a/1700069
def iso_to_gregorian(iso_year: int, iso_week: int, iso_day: int) -> datetime.date:
    """Converts an ISO week date into a datetime object.

    :param iso_year: the year
    :param iso_week: the week number, each year has either 52 or 53 weeks
    :param iso_day: the day numbered 1 through 7, beginning with Monday

    """

    if not 1 <= iso_week <= 53:
        raise ValueError("ISO Calendar week value must be between 1-53.")

    if not 1 <= iso_day <= 7:
        raise ValueError("ISO Calendar day value must be between 1-7")

    # The first week of the year always contains 4 Jan.
    fourth_jan = datetime.date(iso_year, 1, 4)
    delta = datetime.timedelta(fourth_jan.isoweekday() - 1)
    year_start = fourth_jan - delta
    gregorian = year_start + datetime.timedelta(days=iso_day - 1, weeks=iso_week - 1)

    return gregorian


def validate_bounds(bounds: str) -> None:
    if bounds != "()" and bounds != "(]" and bounds != "[)" and bounds != "[]":
        raise ValueError(
            "Invalid bounds. Please select between '()', '(]', '[)', or '[]'."
        )


__all__ = ["next_weekday", "is_timestamp", "validate_ordinal", "iso_to_gregorian"]
