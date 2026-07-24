from _typeshed import Incomplete

from antlr4.InputStream import InputStream as InputStream
from antlr4.ParserRuleContext import ParserRuleContext as ParserRuleContext
from antlr4.Recognizer import Recognizer as Recognizer

class UnsupportedOperationException(Exception):
    def __init__(self, msg: str) -> None: ...

class IllegalStateException(Exception):
    def __init__(self, msg: str) -> None: ...

class CancellationException(IllegalStateException):
    def __init__(self, msg: str) -> None: ...

class RecognitionException(Exception):
    message: Incomplete
    recognizer: Incomplete
    input: Incomplete
    ctx: Incomplete
    offendingToken: Incomplete
    offendingState: int
    def __init__(
        self, message: str | None = None, recognizer: Recognizer | None = None, input: InputStream | None = None, ctx=None
    ) -> None: ...
    def getExpectedTokens(self): ...

class LexerNoViableAltException(RecognitionException):
    startIndex: Incomplete
    deadEndConfigs: Incomplete
    message: str
    def __init__(self, lexer, input: InputStream, startIndex: int, deadEndConfigs) -> None: ...

class NoViableAltException(RecognitionException):
    deadEndConfigs: Incomplete
    startToken: Incomplete
    offendingToken: Incomplete
    def __init__(
        self,
        recognizer,
        input=None,
        startToken=None,
        offendingToken=None,
        deadEndConfigs=None,
        ctx: ParserRuleContext | None = None,
    ) -> None: ...

class InputMismatchException(RecognitionException):
    offendingToken: Incomplete
    def __init__(self, recognizer) -> None: ...

class FailedPredicateException(RecognitionException):
    ruleIndex: Incomplete
    predicateIndex: Incomplete
    predicate: Incomplete
    offendingToken: Incomplete
    def __init__(self, recognizer, predicate: str | None = None, message: str | None = None) -> None: ...
    def formatMessage(self, predicate: str, message: str): ...

class ParseCancellationException(CancellationException): ...
