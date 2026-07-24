import datetime
from collections.abc import Callable
from typing import Final

__all__ = ["main"]

START: Final[datetime.datetime]
END: Final[datetime.datetime]
DISTANCE_FROM_TIMEZONE_CHANGE: Final[datetime.timedelta]

DTS: Final[list[datetime.datetime]]

def main(create_timezones: list[Callable[[str], datetime.tzinfo]], name: str) -> None: ...
