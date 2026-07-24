from antlr4.InputStream import InputStream
from antlr4.Lexer import TokenSource

class Token:
    __slots__ = ("source", "type", "channel", "start", "stop", "tokenIndex", "line", "column", "_text")
    INVALID_TYPE: int
    EPSILON: int
    MIN_USER_TOKEN_TYPE: int
    EOF: int
    DEFAULT_CHANNEL: int
    HIDDEN_CHANNEL: int
    source: tuple[TokenSource | None, InputStream | None]
    type: int
    channel: int
    start: int
    stop: int
    tokenIndex: int | None
    line: int
    column: int
    def __init__(self) -> None: ...
    @property
    def text(self) -> str: ...
    @text.setter
    def text(self, text: str) -> None: ...
    def getTokenSource(self) -> TokenSource | None: ...
    def getInputStream(self) -> InputStream | None: ...

class CommonToken(Token):
    EMPTY_SOURCE: tuple[None, None]
    source: tuple[TokenSource | None, InputStream | None]
    type: int
    channel: int
    start: int
    stop: int
    tokenIndex: int
    line: int
    column: int
    def __init__(
        self,
        source: tuple[TokenSource | None, InputStream | None] = (None, None),
        type: int | None = None,
        channel: int = 0,
        start: int = -1,
        stop: int = -1,
    ) -> None: ...
    def clone(self) -> CommonToken: ...
    @property
    def text(self) -> str: ...
    @text.setter
    def text(self, text: str) -> None: ...
