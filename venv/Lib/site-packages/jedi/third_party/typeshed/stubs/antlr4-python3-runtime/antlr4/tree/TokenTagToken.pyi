from _typeshed import Incomplete

from antlr4.Token import CommonToken as CommonToken

class TokenTagToken(CommonToken):
    __slots__ = ("tokenName", "label")
    tokenName: Incomplete
    label: Incomplete
    def __init__(self, tokenName: str, type: int, label: str | None = None) -> None: ...
    def getText(self): ...
