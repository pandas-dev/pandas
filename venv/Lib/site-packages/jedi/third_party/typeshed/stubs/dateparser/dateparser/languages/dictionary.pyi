import re
from _typeshed import Incomplete
from itertools import chain
from typing import Final, overload

from dateparser.conf import Settings

PARSER_HARDCODED_TOKENS: Final[list[str]]
PARSER_KNOWN_TOKENS: Final[list[str]]
ALWAYS_KEEP_TOKENS: Final[list[str]]
KNOWN_WORD_TOKENS: Final[list[str]]
PARENTHESES_PATTERN: Final[re.Pattern[str]]
NUMERAL_PATTERN: Final[re.Pattern[str]]
KEEP_TOKEN_PATTERN: Final[re.Pattern[str]]

class UnknownTokenError(Exception): ...

class Dictionary:
    info: dict[str, Incomplete]
    def __init__(self, locale_info: dict[str, Incomplete], settings: Settings | None = None) -> None: ...
    def __contains__(self, key: str) -> bool: ...
    def __getitem__(self, key: str): ...
    def __iter__(self) -> chain[str]: ...
    def are_tokens_valid(self, tokens: list[str]) -> bool: ...
    @overload
    def split(self, string: None, keep_formatting: bool = False) -> None: ...
    @overload
    def split(self, string: str, keep_formatting: bool = False) -> list[str]: ...

class NormalizedDictionary(Dictionary):
    def __init__(self, locale_info: dict[str, Incomplete], settings: Settings | None = None) -> None: ...
