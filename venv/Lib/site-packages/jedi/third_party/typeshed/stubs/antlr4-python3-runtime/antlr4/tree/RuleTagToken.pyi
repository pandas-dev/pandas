from _typeshed import Incomplete

from antlr4.Token import Token as Token

class RuleTagToken(Token):
    __slots__ = ("label", "ruleName")
    source: Incomplete
    type: Incomplete
    channel: Incomplete
    start: int
    stop: int
    tokenIndex: int
    line: int
    column: int
    label: Incomplete
    ruleName: Incomplete
    def __init__(self, ruleName: str, bypassTokenType: int, label: str | None = None) -> None: ...
    def getText(self): ...
