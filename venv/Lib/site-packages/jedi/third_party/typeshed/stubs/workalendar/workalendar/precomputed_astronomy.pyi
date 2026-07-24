import datetime
import pathlib
from collections.abc import Callable
from typing import Final

TZAwareDate = datetime.date
YEAR_INTERVAL: Final = 30
TIME_ZONES: Final = ("America/Santiago", "Asia/Hong_Kong", "Asia/Taipei", "Asia/Tokyo")
pre_computed_equinoxes_path: Final[pathlib.Path]
pre_computed_solar_terms_path: Final[pathlib.Path]

def fromisoformat(iso: str) -> datetime.date: ...
def create_astronomical_data(progress: Callable[[int], int] | None = None) -> None: ...
def calculate_equinoxes(year: int, timezone: str = "UTC") -> tuple[TZAwareDate, TZAwareDate]: ...
def solar_term(year: int, degrees: int, timezone: str = "UTC") -> TZAwareDate: ...
