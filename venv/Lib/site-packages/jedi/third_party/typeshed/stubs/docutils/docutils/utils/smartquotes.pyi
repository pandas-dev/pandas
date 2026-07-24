from collections.abc import Generator, Iterable
from re import Pattern
from typing import ClassVar, Final, Literal

options: Final[str]

class smartchars:
    endash: ClassVar[str]
    emdash: ClassVar[str]
    ellipsis: ClassVar[str]
    apostrophe: ClassVar[str]
    quotes: ClassVar[dict[str, str | tuple[str, str, str, str]]]
    language: str
    def __init__(self, language: str = "en") -> None: ...

class RegularExpressions:
    START_SINGLE: ClassVar[Pattern[str]]
    START_DOUBLE: ClassVar[Pattern[str]]
    ADJACENT_1: ClassVar[Pattern[str]]
    ADJACENT_2: ClassVar[Pattern[str]]
    OPEN_SINGLE: ClassVar[Pattern[str]]
    OPEN_DOUBLE: ClassVar[Pattern[str]]
    DECADE: ClassVar[Pattern[str]]
    APOSTROPHE: ClassVar[Pattern[str]]
    OPENING_SECONDARY: ClassVar[Pattern[str]]
    CLOSING_SECONDARY: ClassVar[Pattern[str]]
    OPENING_PRIMARY: ClassVar[Pattern[str]]
    CLOSING_PRIMARY: ClassVar[Pattern[str]]

regexes: RegularExpressions
default_smartypants_attr: Final = "1"

def smartyPants(text: str, attr="1", language: str = "en") -> str: ...
def educate_tokens(text_tokens: Iterable[tuple[str, str]], attr="1", language: str = "en") -> Generator[str]: ...
def educateQuotes(text: str, language: str = "en") -> str: ...
def educateBackticks(text: str, language: str = "en") -> str: ...
def educateSingleBackticks(text: str, language: str = "en") -> str: ...
def educateDashes(text: str) -> str: ...
def educateDashesOldSchool(text: str) -> str: ...
def educateDashesOldSchoolInverted(text: str) -> str: ...
def educateEllipses(text: str) -> str: ...
def stupefyEntities(text: str, language: str = "en") -> str: ...
def processEscapes(text: str, restore: bool = False) -> str: ...
def tokenize(text: str) -> Generator[tuple[Literal["tag", "text"], str]]: ...
