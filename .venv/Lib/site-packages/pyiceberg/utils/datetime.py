#  Licensed to the Apache Software Foundation (ASF) under one
#  or more contributor license agreements.  See the NOTICE file
#  distributed with this work for additional information
#  regarding copyright ownership.  The ASF licenses this file
#  to you under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance
#  with the License.  You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing,
#  software distributed under the License is distributed on an
#  "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
#  KIND, either express or implied.  See the License for the
#  specific language governing permissions and limitations
#  under the License.
"""Helper methods for working with date/time representations."""

from __future__ import annotations

import re
from datetime import (
    date,
    datetime,
    time,
    timedelta,
)

EPOCH_DATE = date.fromisoformat("1970-01-01")
EPOCH_TIMESTAMP = datetime.fromisoformat("1970-01-01T00:00:00.000000")
ISO_TIMESTAMP = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(.\d{1,6})?")
EPOCH_TIMESTAMPTZ = datetime.fromisoformat("1970-01-01T00:00:00.000000+00:00")
ISO_TIMESTAMPTZ = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(.\d{1,6})?[-+]\d{2}:\d{2}")


def micros_to_days(timestamp: int) -> int:
    """Convert a timestamp in microseconds to a date in days."""
    return timedelta(microseconds=timestamp).days


def micros_to_time(micros: int) -> time:
    """Convert a timestamp in microseconds to a time."""
    micros, microseconds = divmod(micros, 1000000)
    micros, seconds = divmod(micros, 60)
    micros, minutes = divmod(micros, 60)
    hours = micros
    return time(hour=hours, minute=minutes, second=seconds, microsecond=microseconds)


def date_str_to_days(date_str: str) -> int:
    """Convert an ISO-8601 formatted date to days from 1970-01-01."""
    return (date.fromisoformat(date_str) - EPOCH_DATE).days


def date_to_days(date_val: date) -> int:
    """Convert a Python date object to days from 1970-01-01."""
    return (date_val - EPOCH_DATE).days


def days_to_date(days: int) -> date:
    """Create a date from the number of days from 1970-01-01."""
    return EPOCH_DATE + timedelta(days)


def time_str_to_micros(time_str: str) -> int:
    """Convert an ISO-8601 formatted time to microseconds from midnight."""
    return time_to_micros(time.fromisoformat(time_str))


def time_to_micros(t: time) -> int:
    """Convert a datetime.time object to microseconds from midnight."""
    return (((t.hour * 60 + t.minute) * 60) + t.second) * 1_000_000 + t.microsecond


def datetime_to_micros(dt: datetime) -> int:
    """Convert a datetime to microseconds from 1970-01-01T00:00:00.000000."""
    if dt.tzinfo:
        delta = dt - EPOCH_TIMESTAMPTZ
    else:
        delta = dt - EPOCH_TIMESTAMP
    return (delta.days * 86400 + delta.seconds) * 1_000_000 + delta.microseconds


def timestamp_to_micros(timestamp_str: str) -> int:
    """Convert an ISO-9601 formatted timestamp without zone to microseconds from 1970-01-01T00:00:00.000000."""
    if ISO_TIMESTAMP.fullmatch(timestamp_str):
        return datetime_to_micros(datetime.fromisoformat(timestamp_str))
    if ISO_TIMESTAMPTZ.fullmatch(timestamp_str):
        # When we can match a timestamp without a zone, we can give a more specific error
        raise ValueError(f"Zone offset provided, but not expected: {timestamp_str}")
    raise ValueError(f"Invalid timestamp without zone: {timestamp_str} (must be ISO-8601)")


def datetime_to_millis(dt: datetime) -> int:
    """Convert a datetime to milliseconds from 1970-01-01T00:00:00.000000."""
    if dt.tzinfo:
        delta = dt - EPOCH_TIMESTAMPTZ
    else:
        delta = dt - EPOCH_TIMESTAMP
    return (delta.days * 86400 + delta.seconds) * 1_000 + delta.microseconds // 1_000


def millis_to_datetime(millis: int) -> datetime:
    """Convert milliseconds from epoch to a timestamp."""
    dt = timedelta(milliseconds=millis)
    return EPOCH_TIMESTAMP + dt


def timestamptz_to_micros(timestamptz_str: str) -> int:
    """Convert an ISO-8601 formatted timestamp with zone to microseconds from 1970-01-01T00:00:00.000000+00:00."""
    if ISO_TIMESTAMPTZ.fullmatch(timestamptz_str):
        return datetime_to_micros(datetime.fromisoformat(timestamptz_str))
    if ISO_TIMESTAMP.fullmatch(timestamptz_str):
        # When we can match a timestamp without a zone, we can give a more specific error
        raise ValueError(f"Missing zone offset: {timestamptz_str} (must be ISO-8601)")
    raise ValueError(f"Invalid timestamp with zone: {timestamptz_str} (must be ISO-8601)")


def micros_to_timestamp(micros: int) -> datetime:
    """Convert microseconds from epoch to a timestamp."""
    dt = timedelta(microseconds=micros)
    return EPOCH_TIMESTAMP + dt


def micros_to_timestamptz(micros: int) -> datetime:
    """Convert microseconds from epoch to an utc timestamp."""
    dt = timedelta(microseconds=micros)
    return EPOCH_TIMESTAMPTZ + dt


def to_human_year(year_ordinal: int) -> str:
    """Convert a DateType value to human string."""
    return f"{EPOCH_TIMESTAMP.year + year_ordinal:0=4d}"


def to_human_month(month_ordinal: int) -> str:
    """Convert a DateType value to human string."""
    return f"{EPOCH_TIMESTAMP.year + month_ordinal // 12:0=4d}-{1 + month_ordinal % 12:0=2d}"


def to_human_day(day_ordinal: int) -> str:
    """Convert a DateType value to human string."""
    return (EPOCH_DATE + timedelta(days=day_ordinal)).isoformat()


def to_human_hour(hour_ordinal: int) -> str:
    """Convert a DateType value to human string."""
    return (EPOCH_TIMESTAMP + timedelta(hours=hour_ordinal)).isoformat("-", "hours")


def to_human_time(micros_from_midnight: int) -> str:
    """Convert a TimeType value to human string."""
    return micros_to_time(micros_from_midnight).isoformat()


def to_human_timestamptz(timestamp_micros: int) -> str:
    """Convert a TimestamptzType value to human string."""
    return (EPOCH_TIMESTAMPTZ + timedelta(microseconds=timestamp_micros)).isoformat()


def to_human_timestamp(timestamp_micros: int) -> str:
    """Convert a TimestampType value to human string."""
    return (EPOCH_TIMESTAMP + timedelta(microseconds=timestamp_micros)).isoformat()


def micros_to_hours(micros: int) -> int:
    """Convert a timestamp in microseconds to hours from 1970-01-01T00:00."""
    return micros // 3_600_000_000


def days_to_months(days: int) -> int:
    d = days_to_date(days)
    return (d.year - EPOCH_DATE.year) * 12 + (d.month - EPOCH_DATE.month)


def micros_to_months(micros: int) -> int:
    dt = micros_to_timestamp(micros)
    return (dt.year - EPOCH_TIMESTAMP.year) * 12 + (dt.month - EPOCH_TIMESTAMP.month)


def days_to_years(days: int) -> int:
    return days_to_date(days).year - EPOCH_DATE.year


def micros_to_years(micros: int) -> int:
    return micros_to_timestamp(micros).year - EPOCH_TIMESTAMP.year
