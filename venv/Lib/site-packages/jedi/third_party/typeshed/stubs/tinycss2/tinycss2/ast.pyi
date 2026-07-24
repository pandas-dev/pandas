from typing import Literal

class Node:
    source_line: int
    source_column: int
    type: str

    def __init__(self, source_line: int, source_column: int) -> None: ...
    def serialize(self) -> str: ...

class ParseError(Node):
    type: Literal["error"]
    kind: str
    message: str
    repr_format: str
    def __init__(self, line: int, column: int, kind: str, message: str) -> None: ...

class Comment(Node):
    type: Literal["comment"]
    value: str
    repr_format: str
    def __init__(self, line: int, column: int, value: str) -> None: ...

class WhitespaceToken(Node):
    type: Literal["whitespace"]
    value: str
    repr_format: str
    def __init__(self, line: int, column: int, value: str) -> None: ...

class LiteralToken(Node):
    type: Literal["literal"]
    value: str
    repr_format: str
    def __init__(self, line: int, column: int, value: str) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...

class IdentToken(Node):
    type: Literal["ident"]
    value: str
    lower_value: str
    repr_format: str
    def __init__(self, line: int, column: int, value: str) -> None: ...

class AtKeywordToken(Node):
    type: Literal["at-keyword"]
    value: str
    lower_value: str
    repr_format: str
    def __init__(self, line: int, column: int, value: str) -> None: ...

class HashToken(Node):
    type: Literal["hash"]
    value: str
    is_identifier: bool
    repr_format: str
    def __init__(self, line: int, column: int, value: str, is_identifier: bool) -> None: ...

class StringToken(Node):
    type: Literal["string"]
    value: str
    representation: str
    repr_format: str
    def __init__(self, line: int, column: int, value: str, representation: str) -> None: ...

class URLToken(Node):
    type: Literal["url"]
    value: str
    representation: str
    repr_format: str
    def __init__(self, line: int, column: int, value: str, representation: str) -> None: ...

class UnicodeRangeToken(Node):
    type: Literal["unicode-range"]
    start: int
    end: int
    repr_format: str
    def __init__(self, line: int, column: int, start: int, end: int) -> None: ...

class NumberToken(Node):
    type: Literal["number"]
    value: float
    int_value: int | None
    is_integer: bool
    representation: str
    repr_format: str
    def __init__(self, line: int, column: int, value: float, int_value: int | None, representation: str) -> None: ...

class PercentageToken(Node):
    type: Literal["percentage"]
    value: float
    int_value: int | None
    is_integer: bool
    representation: str
    repr_format: str
    def __init__(self, line: int, column: int, value: float, int_value: int | None, representation: str) -> None: ...

class DimensionToken(Node):
    type: Literal["dimension"]
    value: float
    int_value: int | None
    is_integer: bool
    representation: str
    unit: str
    lower_unit: str
    repr_format: str
    def __init__(self, line: int, column: int, value: float, int_value: int | None, representation: str, unit: str) -> None: ...

class ParenthesesBlock(Node):
    type: Literal["() block"]
    content: list[Node]
    repr_format: str
    def __init__(self, line: int, column: int, content: list[Node]) -> None: ...

class SquareBracketsBlock(Node):
    type: Literal["[] block"]
    content: list[Node]
    repr_format: str
    def __init__(self, line: int, column: int, content: list[Node]) -> None: ...

class CurlyBracketsBlock(Node):
    type: Literal["{} block"]
    content: list[Node]
    repr_format: str
    def __init__(self, line: int, column: int, content: list[Node]) -> None: ...

class FunctionBlock(Node):
    type: Literal["function"]
    name: str
    lower_name: str
    arguments: list[Node]
    repr_format: str
    def __init__(self, line: int, column: int, name: str, arguments: list[Node]) -> None: ...

class Declaration(Node):
    type: Literal["declaration"]
    name: str
    lower_name: str
    value: list[Node]
    important: bool
    repr_format: str
    def __init__(self, line: int, column: int, name: str, lower_name: str, value: list[Node], important: bool) -> None: ...

class QualifiedRule(Node):
    type: Literal["qualified-rule"]
    prelude: list[Node]
    content: list[Node]
    repr_format: str
    def __init__(self, line: int, column: int, prelude: list[Node], content: list[Node]) -> None: ...

class AtRule(Node):
    type: Literal["at-rule"]
    at_keyword: str
    lower_at_keyword: str
    prelude: list[Node]
    content: list[Node] | None
    repr_format: str
    def __init__(
        self, line: int, column: int, at_keyword: str, lower_at_keyword: str, prelude: list[Node], content: list[Node] | None
    ) -> None: ...
