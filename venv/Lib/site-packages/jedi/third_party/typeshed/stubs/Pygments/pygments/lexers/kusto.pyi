from typing import Final

from ..lexer import RegexLexer

__all__ = ["KustoLexer"]

KUSTO_KEYWORDS: Final[list[str]]
KUSTO_PUNCTUATION: Final[list[str]]

class KustoLexer(RegexLexer): ...
