from _typeshed import Incomplete

from antlr4.atn.ATNState import ATNState as ATNState, DecisionState as DecisionState
from antlr4.atn.ATNType import ATNType as ATNType
from antlr4.IntervalSet import IntervalSet as IntervalSet
from antlr4.RuleContext import RuleContext as RuleContext
from antlr4.Token import Token as Token

class ATN:
    __slots__ = (
        "grammarType",
        "maxTokenType",
        "states",
        "decisionToState",
        "ruleToStartState",
        "ruleToStopState",
        "modeNameToStartState",
        "ruleToTokenType",
        "lexerActions",
        "modeToStartState",
    )
    INVALID_ALT_NUMBER: int
    grammarType: Incomplete
    maxTokenType: Incomplete
    states: Incomplete
    decisionToState: Incomplete
    ruleToStartState: Incomplete
    ruleToStopState: Incomplete
    modeNameToStartState: Incomplete
    ruleToTokenType: Incomplete
    lexerActions: Incomplete
    modeToStartState: Incomplete
    def __init__(self, grammarType: ATNType, maxTokenType: int) -> None: ...
    def nextTokensInContext(self, s: ATNState, ctx: RuleContext): ...
    def nextTokensNoContext(self, s: ATNState): ...
    def nextTokens(self, s: ATNState, ctx: RuleContext | None = None): ...
    def addState(self, state: ATNState): ...
    def removeState(self, state: ATNState): ...
    def defineDecisionState(self, s: DecisionState): ...
    def getDecisionState(self, decision: int): ...
    def getExpectedTokens(self, stateNumber: int, ctx: RuleContext): ...
