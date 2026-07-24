from collections.abc import Sequence
from typing import Any

class JMESPathError(ValueError): ...

class ParseError(JMESPathError):
    lex_position: int
    token_value: str
    token_type: str
    msg: str
    expression: str | None
    def __init__(
        self, lex_position: int, token_value: str, token_type: str, msg: str = "Invalid jmespath expression"
    ) -> None: ...

class IncompleteExpressionError(ParseError):
    # When ParseError is used directly, the token always have a non-null value and type
    token_value: str | None  # type: ignore[assignment]
    token_type: str | None  # type: ignore[assignment]
    expression: str
    def set_expression(self, expression: str) -> None: ...

class LexerError(ParseError):
    lexer_position: int
    lexer_value: str
    message: str
    def __init__(self, lexer_position: int, lexer_value: str, message: str, expression: str | None = None) -> None: ...

class ArityError(ParseError):
    expected_arity: int
    actual_arity: int
    function_name: str
    def __init__(self, expected: int, actual: int, name: str) -> None: ...

class VariadictArityError(ArityError): ...

class JMESPathTypeError(JMESPathError):
    function_name: str
    current_value: Any
    actual_type: str
    expected_types: Sequence[str]
    def __init__(self, function_name: str, current_value: Any, actual_type: str, expected_types: Sequence[str]) -> None: ...

class EmptyExpressionError(JMESPathError):
    def __init__(self) -> None: ...

class UnknownFunctionError(JMESPathError): ...
