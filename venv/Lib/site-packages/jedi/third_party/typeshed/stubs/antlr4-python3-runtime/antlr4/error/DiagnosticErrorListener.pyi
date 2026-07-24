from _typeshed import Incomplete

from antlr4 import DFA as DFA
from antlr4.atn.ATNConfigSet import ATNConfigSet as ATNConfigSet
from antlr4.error.ErrorListener import ErrorListener as ErrorListener

class DiagnosticErrorListener(ErrorListener):
    exactOnly: Incomplete
    def __init__(self, exactOnly: bool = True) -> None: ...
    def reportAmbiguity(
        self, recognizer, dfa: DFA, startIndex: int, stopIndex: int, exact: bool, ambigAlts: set[int], configs: ATNConfigSet
    ): ...
    def reportAttemptingFullContext(
        self, recognizer, dfa: DFA, startIndex: int, stopIndex: int, conflictingAlts: set[int], configs: ATNConfigSet
    ): ...
    def reportContextSensitivity(
        self, recognizer, dfa: DFA, startIndex: int, stopIndex: int, prediction: int, configs: ATNConfigSet
    ): ...
    def getDecisionDescription(self, recognizer, dfa: DFA): ...
    def getConflictingAlts(self, reportedAlts: set[int], configs: ATNConfigSet): ...
