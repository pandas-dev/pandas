from antlr4.InputStream import InputStream
from antlr4.Lexer import TokenSource
from antlr4.Token import CommonToken as CommonToken

class TokenFactory: ...

class CommonTokenFactory(TokenFactory):
    __slots__ = "copyText"
    DEFAULT: CommonTokenFactory | None
    copyText: bool
    def __init__(self, copyText: bool = False) -> None: ...
    def create(
        self,
        source: tuple[TokenSource, InputStream],
        type: int,
        text: str,
        channel: int,
        start: int,
        stop: int,
        line: int,
        column: int,
    ) -> CommonToken: ...
    def createThin(self, type: int, text: str) -> CommonToken: ...
