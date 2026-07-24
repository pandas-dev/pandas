from _typeshed import Incomplete

from antlr4.atn.ATN import ATN as ATN
from antlr4.atn.ATNConfig import LexerATNConfig as LexerATNConfig
from antlr4.atn.ATNConfigSet import ATNConfigSet as ATNConfigSet, OrderedATNConfigSet as OrderedATNConfigSet
from antlr4.atn.ATNSimulator import ATNSimulator as ATNSimulator
from antlr4.atn.ATNState import ATNState as ATNState, RuleStopState as RuleStopState
from antlr4.atn.LexerActionExecutor import LexerActionExecutor as LexerActionExecutor
from antlr4.atn.Transition import Transition as Transition
from antlr4.dfa.DFA import DFA
from antlr4.dfa.DFAState import DFAState as DFAState
from antlr4.error.Errors import (
    LexerNoViableAltException as LexerNoViableAltException,
    UnsupportedOperationException as UnsupportedOperationException,
)
from antlr4.InputStream import InputStream as InputStream
from antlr4.PredictionContext import (
    PredictionContext as PredictionContext,
    PredictionContextCache as PredictionContextCache,
    SingletonPredictionContext as SingletonPredictionContext,
)
from antlr4.Token import Token as Token

class SimState:
    __slots__ = ("index", "line", "column", "dfaState")
    def __init__(self) -> None: ...
    index: int
    line: int
    column: int
    dfaState: Incomplete
    def reset(self) -> None: ...

class LexerATNSimulator(ATNSimulator):
    __slots__ = ("decisionToDFA", "recog", "startIndex", "line", "column", "mode", "DEFAULT_MODE", "MAX_CHAR_VALUE", "prevAccept")
    debug: bool
    dfa_debug: bool
    MIN_DFA_EDGE: int
    MAX_DFA_EDGE: int
    ERROR: Incomplete
    decisionToDFA: Incomplete
    recog: Incomplete
    startIndex: int
    line: int
    column: int
    mode: Incomplete
    DEFAULT_MODE: Incomplete
    MAX_CHAR_VALUE: Incomplete
    prevAccept: Incomplete
    def __init__(self, recog, atn: ATN, decisionToDFA: list[DFA], sharedContextCache: PredictionContextCache) -> None: ...
    def copyState(self, simulator: LexerATNSimulator): ...
    def match(self, input: InputStream, mode: int): ...
    def reset(self) -> None: ...
    def matchATN(self, input: InputStream): ...
    def execATN(self, input: InputStream, ds0: DFAState): ...
    def getExistingTargetState(self, s: DFAState, t: int): ...
    def computeTargetState(self, input: InputStream, s: DFAState, t: int): ...
    def failOrAccept(self, prevAccept: SimState, input: InputStream, reach: ATNConfigSet, t: int): ...
    def getReachableConfigSet(self, input: InputStream, closure: ATNConfigSet, reach: ATNConfigSet, t: int): ...
    def accept(
        self, input: InputStream, lexerActionExecutor: LexerActionExecutor, startIndex: int, index: int, line: int, charPos: int
    ): ...
    def getReachableTarget(self, trans: Transition, t: int): ...
    def computeStartState(self, input: InputStream, p: ATNState): ...
    def closure(
        self,
        input: InputStream,
        config: LexerATNConfig,
        configs: ATNConfigSet,
        currentAltReachedAcceptState: bool,
        speculative: bool,
        treatEofAsEpsilon: bool,
    ): ...
    def getEpsilonTarget(
        self,
        input: InputStream,
        config: LexerATNConfig,
        t: Transition,
        configs: ATNConfigSet,
        speculative: bool,
        treatEofAsEpsilon: bool,
    ): ...
    def evaluatePredicate(self, input: InputStream, ruleIndex: int, predIndex: int, speculative: bool): ...
    def captureSimState(self, settings: SimState, input: InputStream, dfaState: DFAState): ...
    def addDFAEdge(self, from_: DFAState, tk: int, to: DFAState | None = None, cfgs: ATNConfigSet | None = None) -> DFAState: ...
    def addDFAState(self, configs: ATNConfigSet) -> DFAState: ...
    def getDFA(self, mode: int): ...
    def getText(self, input: InputStream): ...
    def consume(self, input: InputStream): ...
    def getTokenName(self, t: int): ...
