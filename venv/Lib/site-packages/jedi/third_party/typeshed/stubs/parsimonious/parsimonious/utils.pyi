import ast
from typing import Any

class StrAndRepr: ...

def evaluate_string(string: str | ast.AST) -> Any: ...

class Token(StrAndRepr):
    __slots__ = ["type"]
    type: str
    def __init__(self, type: str) -> None: ...
