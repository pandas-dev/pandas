from _typeshed import Incomplete

from antlr4.atn.ATN import ATN as ATN
from antlr4.atn.ATNConfig import ATNConfig as ATNConfig
from antlr4.atn.ATNState import ATNState as ATNState, RuleStopState as RuleStopState
from antlr4.atn.Transition import (
    AbstractPredicateTransition as AbstractPredicateTransition,
    NotSetTransition as NotSetTransition,
    RuleTransition as RuleTransition,
    WildcardTransition as WildcardTransition,
)
from antlr4.IntervalSet import IntervalSet as IntervalSet
from antlr4.PredictionContext import (
    PredictionContext as PredictionContext,
    PredictionContextFromRuleContext as PredictionContextFromRuleContext,
    SingletonPredictionContext as SingletonPredictionContext,
)
from antlr4.RuleContext import RuleContext as RuleContext
from antlr4.Token import Token as Token

class LL1Analyzer:
    __slots__ = "atn"
    HIT_PRED: Incomplete
    atn: Incomplete
    def __init__(self, atn: ATN) -> None: ...
    def getDecisionLookahead(self, s: ATNState): ...
    def LOOK(self, s: ATNState, stopState: ATNState | None = None, ctx: RuleContext | None = None): ...
