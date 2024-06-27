from __future__ import annotations

from typing import TYPE_CHECKING

from isoduration.constants import PERIOD_PREFIX, TIME_PREFIX
from isoduration.formatter.checking import validate_date_duration

if TYPE_CHECKING:  # pragma: no cover
    from isoduration.types import DateDuration, TimeDuration


def format_date(date_duration: DateDuration, global_sign: int) -> str:
    validate_date_duration(date_duration)

    date_duration_str = PERIOD_PREFIX

    if date_duration.weeks != 0:
        date_duration_str += f"{(date_duration.weeks * global_sign):g}W"

    if date_duration.years != 0:
        date_duration_str += f"{(date_duration.years * global_sign):g}Y"
    if date_duration.months != 0:
        date_duration_str += f"{(date_duration.months * global_sign):g}M"
    if date_duration.days != 0:
        date_duration_str += f"{(date_duration.days * global_sign):g}D"

    return date_duration_str


def format_time(time_duration: TimeDuration, global_sign: int) -> str:
    time_duration_str = TIME_PREFIX

    if time_duration.hours != 0:
        time_duration_str += f"{(time_duration.hours * global_sign):g}H"
    if time_duration.minutes != 0:
        time_duration_str += f"{(time_duration.minutes * global_sign):g}M"
    if time_duration.seconds != 0:
        time_duration_str += f"{(time_duration.seconds * global_sign):g}S"

    if time_duration_str == TIME_PREFIX:
        return ""
    return time_duration_str
