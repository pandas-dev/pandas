from _typeshed import Incomplete
from re import Pattern
from typing import Final, Literal
from typing_extensions import TypeAlias

_TokenTypesNotOperand: TypeAlias = Literal[
    "LITERAL", "FUNC", "ARRAY", "PAREN", "SEP", "OPERATOR-PREFIX", "OPERATOR-INFIX", "OPERATOR-POSTFIX", "WHITE-SPACE"
]
_TokenTypes: TypeAlias = Literal["OPERAND", _TokenTypesNotOperand]
_TokenOperandSubtypes: TypeAlias = Literal["TEXT", "NUMBER", "LOGICAL", "ERROR", "RANGE"]
_TokenSubtypes: TypeAlias = Literal["", _TokenOperandSubtypes, "OPEN", "CLOSE", "ARG", "ROW"]

class TokenizerError(Exception): ...

class Tokenizer:
    SN_RE: Final[Pattern[str]]
    WSPACE_RE: Final[Pattern[str]]
    STRING_REGEXES: Final[dict[str, Pattern[str]]]
    ERROR_CODES: Final[tuple[str, ...]]
    TOKEN_ENDERS: Final = ",;}) +-*/^&=><%"
    formula: Incomplete
    items: Incomplete
    token_stack: Incomplete
    offset: int
    token: Incomplete
    def __init__(self, formula) -> None: ...
    def check_scientific_notation(self): ...
    def assert_empty_token(self, can_follow=()) -> None: ...
    def save_token(self) -> None: ...
    def render(self): ...

class Token:
    __slots__ = ["value", "type", "subtype"]
    LITERAL: Final = "LITERAL"
    OPERAND: Final = "OPERAND"
    FUNC: Final = "FUNC"
    ARRAY: Final = "ARRAY"
    PAREN: Final = "PAREN"
    SEP: Final = "SEP"
    OP_PRE: Final = "OPERATOR-PREFIX"
    OP_IN: Final = "OPERATOR-INFIX"
    OP_POST: Final = "OPERATOR-POSTFIX"
    WSPACE: Final = "WHITE-SPACE"
    value: Incomplete
    type: _TokenTypes
    subtype: _TokenSubtypes
    def __init__(self, value, type_: _TokenTypes, subtype: _TokenSubtypes = "") -> None: ...
    TEXT: Final = "TEXT"
    NUMBER: Final = "NUMBER"
    LOGICAL: Final = "LOGICAL"
    ERROR: Final = "ERROR"
    RANGE: Final = "RANGE"
    @classmethod
    def make_operand(cls, value): ...
    OPEN: Final = "OPEN"
    CLOSE: Final = "CLOSE"
    @classmethod
    def make_subexp(cls, value, func: bool = False): ...
    def get_closer(self): ...
    ARG: Final = "ARG"
    ROW: Final = "ROW"
    @classmethod
    def make_separator(cls, value): ...
