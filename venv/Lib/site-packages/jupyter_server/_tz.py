"""
Timezone utilities

Just UTC-awareness right now
"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

from datetime import datetime, timedelta, timezone, tzinfo

# constant for zero offset
ZERO = timedelta(0)


class tzUTC(tzinfo):  # noqa: N801
    """tzinfo object for UTC (zero offset)"""

    def utcoffset(self, d: datetime | None) -> timedelta:
        """Compute utcoffset."""
        return ZERO

    def dst(self, d: datetime | None) -> timedelta:
        """Compute dst."""
        return ZERO


def utcnow() -> datetime:
    """Return timezone-aware UTC timestamp"""
    return datetime.now(timezone.utc)


def utcfromtimestamp(timestamp: float) -> datetime:
    return datetime.fromtimestamp(timestamp, timezone.utc)


UTC = tzUTC()  # type:ignore[abstract]


def isoformat(dt: datetime) -> str:
    """Return iso-formatted timestamp

    Like .isoformat(), but uses Z for UTC instead of +00:00
    """
    return dt.isoformat().replace("+00:00", "Z")
