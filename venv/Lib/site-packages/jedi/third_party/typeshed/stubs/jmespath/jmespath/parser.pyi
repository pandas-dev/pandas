from collections.abc import Iterator
from typing import Any, ClassVar

from jmespath.lexer import _LexerTokenizeResult
from jmespath.visitor import Options, _TreeNode

class Parser:
    BINDING_POWER: ClassVar[dict[str, int]]
    tokenizer: Iterator[_LexerTokenizeResult] | None
    def __init__(self, lookahead: int = 2) -> None: ...
    def parse(self, expression: str) -> ParsedResult: ...
    @classmethod
    def purge(cls) -> None: ...

class ParsedResult:
    expression: str
    parsed: _TreeNode
    def __init__(self, expression: str, parsed: _TreeNode) -> None: ...
    def search(self, value: Any, options: Options | None = None) -> Any: ...
