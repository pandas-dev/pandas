from _typeshed import Incomplete

from antlr4.atn.ATN import ATN as ATN
from antlr4.atn.ATNConfigSet import ATNConfigSet as ATNConfigSet
from antlr4.dfa.DFAState import DFAState as DFAState
from antlr4.PredictionContext import (
    PredictionContext as PredictionContext,
    PredictionContextCache as PredictionContextCache,
    getCachedPredictionContext as getCachedPredictionContext,
)

class ATNSimulator:
    __slots__ = ("atn", "sharedContextCache", "__dict__")
    ERROR: Incomplete
    atn: Incomplete
    sharedContextCache: Incomplete
    def __init__(self, atn: ATN, sharedContextCache: PredictionContextCache) -> None: ...
    def getCachedContext(self, context: PredictionContext): ...
