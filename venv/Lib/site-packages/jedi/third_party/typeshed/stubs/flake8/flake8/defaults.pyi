from re import Pattern
from typing import Final

EXCLUDE: Final[tuple[str, ...]]
IGNORE: Final[tuple[str, ...]]
MAX_LINE_LENGTH: Final = 79
INDENT_SIZE: Final = 4
WHITESPACE: Final[frozenset[str]]
STATISTIC_NAMES: Final[tuple[str, ...]]
NOQA_INLINE_REGEXP: Final[Pattern[str]]
NOQA_FILE: Final[Pattern[str]]
VALID_CODE_PREFIX: Final[Pattern[str]]
