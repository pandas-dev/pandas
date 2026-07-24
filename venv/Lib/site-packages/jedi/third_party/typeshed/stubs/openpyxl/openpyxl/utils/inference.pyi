from re import Pattern
from typing import Final

PERCENT_REGEX: Final[Pattern[str]]
TIME_REGEX: Final[Pattern[str]]
NUMBER_REGEX: Final[Pattern[str]]

def cast_numeric(value): ...
def cast_percentage(value): ...
def cast_time(value): ...
