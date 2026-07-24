import datetime
import optparse
from collections.abc import Sequence
from typing import Final, Literal

from .base import Component

version: Final[str]

def change_tz(
    cal: Component,
    new_timezone: datetime._TzInfo | None,
    default: datetime._TzInfo | None,
    utc_only: bool = False,
    utc_tz: datetime._TzInfo | None = ...,
) -> None: ...
def show_timezones() -> None: ...
def convert_events(utc_only: bool, args: Sequence[str]) -> None: ...
def main() -> None: ...
def get_options() -> tuple[optparse.Values, Literal[False]] | tuple[optparse.Values, list[str]]: ...
