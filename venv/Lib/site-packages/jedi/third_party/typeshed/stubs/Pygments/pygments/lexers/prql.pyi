from _typeshed import Incomplete
from typing import ClassVar

from ..lexer import RegexLexer

__all__ = ["PrqlLexer"]

class PrqlLexer(RegexLexer):
    builtinTypes: ClassVar[Incomplete]
    def innerstring_rules(ttype) -> list[tuple[str, Incomplete]]: ...
    def fstring_rules(ttype) -> list[tuple[str, Incomplete]]: ...
