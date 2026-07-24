from collections.abc import Mapping
from typing import Any, Final
from typing_extensions import Self

class _TokenType(tuple[str, ...]):
    parent: _TokenType | None
    def split(self) -> list[_TokenType]: ...
    subtypes: set[_TokenType]
    def __contains__(self, val: _TokenType) -> bool: ...  # type: ignore[override]
    def __getattr__(self, name: str) -> _TokenType: ...
    def __copy__(self) -> Self: ...
    def __deepcopy__(self, memo: Any) -> Self: ...

Token: _TokenType
Text: _TokenType
Whitespace: _TokenType
Escape: _TokenType
Error: _TokenType
Other: _TokenType
Keyword: _TokenType
Name: _TokenType
Literal: _TokenType
String: _TokenType
Number: _TokenType
Punctuation: _TokenType
Operator: _TokenType
Comment: _TokenType
Generic: _TokenType

def is_token_subtype(ttype: _TokenType, other: _TokenType) -> bool: ...
def string_to_tokentype(s: str | _TokenType) -> _TokenType: ...

# dict, but shouldn't be mutated
STANDARD_TYPES: Final[Mapping[_TokenType, str]]
