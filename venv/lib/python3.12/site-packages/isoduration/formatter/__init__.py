from __future__ import annotations

from typing import TYPE_CHECKING

from isoduration.constants import PERIOD_PREFIX
from isoduration.formatter.checking import check_global_sign
from isoduration.formatter.formatting import format_date, format_time

if TYPE_CHECKING:  # pragma: no cover
    from isoduration.types import Duration


def format_duration(duration: Duration) -> str:
    global_sign = check_global_sign(duration)
    date_duration_str = format_date(duration.date, global_sign)
    time_duration_str = format_time(duration.time, global_sign)

    duration_str = f"{date_duration_str}{time_duration_str}"
    sign_str = "-" if global_sign < 0 else ""

    if duration_str == PERIOD_PREFIX:
        return f"{PERIOD_PREFIX}0D"

    return f"{sign_str}{duration_str}"
