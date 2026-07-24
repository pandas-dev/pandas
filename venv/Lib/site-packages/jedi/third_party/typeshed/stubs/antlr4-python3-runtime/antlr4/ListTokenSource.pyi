from _typeshed import Incomplete

from antlr4.CommonTokenFactory import CommonTokenFactory as CommonTokenFactory
from antlr4.Lexer import TokenSource as TokenSource
from antlr4.Token import Token as Token

class ListTokenSource(TokenSource):
    __slots__ = ("tokens", "sourceName", "pos", "eofToken", "_factory")
    tokens: Incomplete
    sourceName: Incomplete
    pos: int
    eofToken: Incomplete
    def __init__(self, tokens: list[Token], sourceName: str | None = None) -> None: ...
    @property
    def column(self): ...
    def nextToken(self): ...
    @property
    def line(self): ...
    def getInputStream(self): ...
    def getSourceName(self): ...
