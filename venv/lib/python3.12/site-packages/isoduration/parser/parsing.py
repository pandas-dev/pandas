from collections import OrderedDict
from decimal import Decimal, InvalidOperation

import arrow  # type: ignore

from isoduration.parser.exceptions import (
    IncorrectDesignator,
    NoTime,
    OutOfDesignators,
    UnknownToken,
    UnparseableValue,
)
from isoduration.parser.util import (
    is_letter,
    is_number,
    is_time,
    is_week,
    parse_designator,
)
from isoduration.types import DateDuration, Duration, TimeDuration


def parse_datetime_duration(duration_str: str, sign: int) -> Duration:
    try:
        duration: arrow.Arrow = arrow.get(duration_str)
    except (arrow.ParserError, ValueError):
        raise UnparseableValue(f"Value could not be parsed as datetime: {duration_str}")

    return Duration(
        DateDuration(
            years=sign * duration.year,
            months=sign * duration.month,
            days=sign * duration.day,
        ),
        TimeDuration(
            hours=sign * duration.hour,
            minutes=sign * duration.minute,
            seconds=sign * duration.second,
        ),
    )


def parse_date_duration(date_str: str, sign: int) -> Duration:
    date_designators = OrderedDict(
        (("Y", "years"), ("M", "months"), ("D", "days"), ("W", "weeks"))
    )

    duration = DateDuration()
    tmp_value = ""

    for idx, ch in enumerate(date_str):
        if is_time(ch):
            if tmp_value != "" and tmp_value == date_str[:idx]:
                # PYYYY-MM-DDThh:mm:ss
                # PYYYYMMDDThhmmss
                return parse_datetime_duration(date_str, sign)

            time_idx = idx + 1
            time_str = date_str[time_idx:]

            if time_str == "":
                raise NoTime("Wanted time, no time provided")

            return Duration(duration, parse_time_duration(time_str, sign))

        if is_letter(ch):
            try:
                key = parse_designator(date_designators, ch)
                value = sign * Decimal(tmp_value)
            except OutOfDesignators as exc:
                raise IncorrectDesignator(
                    f"Wrong date designator, or designator in the wrong order: {ch}"
                ) from exc
            except InvalidOperation as exc:
                raise UnparseableValue(
                    f"Value could not be parsed as decimal: {tmp_value}"
                ) from exc

            if is_week(ch) and duration != DateDuration():
                raise IncorrectDesignator(
                    "Week is incompatible with any other date designator"
                )

            setattr(duration, key, value)
            tmp_value = ""

            continue

        if is_number(ch):
            if ch == ",":
                tmp_value += "."
            else:
                tmp_value += ch

            continue

        raise UnknownToken(f"Token not recognizable: {ch}")

    return Duration(duration, TimeDuration())


def parse_time_duration(time_str: str, sign: int) -> TimeDuration:
    time_designators = OrderedDict((("H", "hours"), ("M", "minutes"), ("S", "seconds")))

    duration = TimeDuration()
    tmp_value = ""

    for ch in time_str:
        if is_letter(ch):
            try:
                key = parse_designator(time_designators, ch)
                value = sign * Decimal(tmp_value)
            except OutOfDesignators as exc:
                raise IncorrectDesignator(
                    f"Wrong time designator, or designator in the wrong order: {ch}"
                ) from exc
            except InvalidOperation as exc:
                raise UnparseableValue(
                    f"Value could not be parsed as decimal: {tmp_value}"
                ) from exc

            setattr(duration, key, value)
            tmp_value = ""

            continue

        if is_number(ch):
            if ch == ",":
                tmp_value += "."
            else:
                tmp_value += ch

            continue

        raise UnknownToken(f"Token not recognizable: {ch}")

    return duration
