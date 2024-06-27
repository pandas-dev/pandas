from __future__ import annotations

from datetime import datetime
from decimal import ROUND_HALF_UP, Decimal
from typing import TYPE_CHECKING

from isoduration.operations.util import max_day_in_month, mod2, mod3, quot2, quot3

if TYPE_CHECKING:  # pragma: no cover
    from isoduration.types import Duration


def add(start: datetime, duration: Duration) -> datetime:
    """
    https://www.w3.org/TR/xmlschema-2/#adding-durations-to-dateTimes
    """

    # Months.
    temp = Decimal(start.month) + duration.date.months
    end_month = mod3(temp, Decimal(1), Decimal(13))
    carry = quot3(temp, Decimal(1), Decimal(13))

    # Years.
    end_year = Decimal(start.year) + duration.date.years + carry

    # Zone.
    end_tzinfo = start.tzinfo

    # Seconds.
    temp = Decimal(start.second) + duration.time.seconds
    end_second = mod2(temp, Decimal("60"))
    carry = quot2(temp, Decimal("60"))

    # Minutes.
    temp = Decimal(start.minute) + duration.time.minutes + carry
    end_minute = mod2(temp, Decimal("60"))
    carry = quot2(temp, Decimal("60"))

    # Hours.
    temp = Decimal(start.hour) + duration.time.hours + carry
    end_hour = mod2(temp, Decimal("24"))
    carry = quot2(temp, Decimal("24"))

    # Days.
    end_max_day_in_month = max_day_in_month(end_year, end_month)

    if start.day > end_max_day_in_month:
        temp = end_max_day_in_month
    else:
        temp = Decimal(start.day)

    end_day = temp + duration.date.days + (7 * duration.date.weeks) + carry

    while True:
        if end_day < 1:
            end_day += max_day_in_month(end_year, end_month - 1)
            carry = Decimal(-1)
        elif end_day > max_day_in_month(end_year, end_month):
            end_day -= max_day_in_month(end_year, end_month)
            carry = Decimal(1)
        else:
            break

        temp = end_month + carry
        end_month = mod3(temp, Decimal(1), Decimal(13))
        end_year = end_year + quot3(temp, Decimal(1), Decimal(13))

    return datetime(
        year=int(end_year.to_integral_value(ROUND_HALF_UP)),
        month=int(end_month.to_integral_value(ROUND_HALF_UP)),
        day=int(end_day.to_integral_value(ROUND_HALF_UP)),
        hour=int(end_hour.to_integral_value(ROUND_HALF_UP)),
        minute=int(end_minute.to_integral_value(ROUND_HALF_UP)),
        second=int(end_second.to_integral_value(ROUND_HALF_UP)),
        tzinfo=end_tzinfo,
    )
