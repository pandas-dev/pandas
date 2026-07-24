from _typeshed import Incomplete

from antlr4.atn.ATNState import ATNState as ATNState, DecisionState as DecisionState
from antlr4.atn.LexerActionExecutor import LexerActionExecutor as LexerActionExecutor
from antlr4.atn.SemanticContext import SemanticContext as SemanticContext
from antlr4.PredictionContext import PredictionContext as PredictionContext

class ATNConfig:
    __slots__ = ("state", "alt", "context", "semanticContext", "reachesIntoOuterContext", "precedenceFilterSuppressed")
    state: Incomplete
    alt: Incomplete
    context: Incomplete
    semanticContext: Incomplete
    reachesIntoOuterContext: Incomplete
    precedenceFilterSuppressed: Incomplete
    def __init__(
        self,
        state: ATNState | None = None,
        alt: int | None = None,
        context: PredictionContext | None = None,
        semantic: SemanticContext | None = None,
        config: ATNConfig | None = None,
    ) -> None: ...
    def __eq__(self, other): ...
    def __hash__(self): ...
    def hashCodeForConfigSet(self): ...
    def equalsForConfigSet(self, other): ...

class LexerATNConfig(ATNConfig):
    __slots__ = ("lexerActionExecutor", "passedThroughNonGreedyDecision")
    lexerActionExecutor: Incomplete
    passedThroughNonGreedyDecision: Incomplete
    def __init__(
        self,
        state: ATNState,
        alt: int | None = None,
        context: PredictionContext | None = None,
        semantic: SemanticContext = ...,
        lexerActionExecutor: LexerActionExecutor | None = None,
        config: LexerATNConfig | None = None,
    ) -> None: ...
    def __hash__(self): ...
    def __eq__(self, other): ...
    def hashCodeForConfigSet(self): ...
    def equalsForConfigSet(self, other): ...
    def checkNonGreedyDecision(self, source: LexerATNConfig, target: ATNState): ...
