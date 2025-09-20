import abc
import re
from typing import Final, Generic, Literal, Self, TypeAlias, TypeVar, type_check_only

import optype.numpy as onp

__all__ = ["BadFortranFormat", "ExpFormat", "FortranFormatParser", "IntFormat"]

_NumberT = TypeVar("_NumberT", int, float)
_TokenType: TypeAlias = Literal["INT", "INT_ID", "EXP_ID", "DOT", "LPAR", "RPAR"]

TOKENS: Final[dict[_TokenType, str]]

class BadFortranFormat(SyntaxError): ...

@type_check_only
class _NumberFormat(Generic[_NumberT], metaclass=abc.ABCMeta):
    width: Final[int]
    repeat: Final[int | None]
    min: Final[int | None]
    @property
    def fortran_format(self, /) -> str: ...
    @property
    def python_format(self, /) -> str: ...
    @classmethod
    def from_number(cls, n: _NumberT, min: int | None = None) -> Self: ...

class IntFormat(_NumberFormat[int]):
    def __init__(self, /, width: int, min: int | None = None, repeat: int | None = None) -> None: ...

class ExpFormat(_NumberFormat[float]):
    significand: Final[int]
    def __init__(self, /, width: int, significand: int, min: int | None = None, repeat: int | None = None) -> None: ...

class Token:
    type: Final[_TokenType]
    value: Final[str]
    pos: Final[int]
    def __init__(self, /, type: _TokenType, value: str, pos: int) -> None: ...

class Tokenizer:
    tokens: Final[list[_TokenType]]
    res: Final[list[re.Pattern[str]]]
    data: str
    curpos: int
    len: int
    def input(self, /, s: str) -> None: ...
    def next_token(self, /) -> Token: ...

class FortranFormatParser:
    tokenizer: Final[Tokenizer]
    def parse(self, /, s: str) -> IntFormat | ExpFormat: ...

def number_digits(n: onp.ToInt) -> int: ...
