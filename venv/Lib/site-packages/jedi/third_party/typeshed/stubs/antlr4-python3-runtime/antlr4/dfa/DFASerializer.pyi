from _typeshed import Incomplete

from antlr4 import DFA as DFA
from antlr4.dfa.DFAState import DFAState as DFAState
from antlr4.Utils import str_list as str_list

class DFASerializer:
    __slots__ = ("dfa", "literalNames", "symbolicNames")
    dfa: Incomplete
    literalNames: list[str] | None
    symbolicNames: list[str] | None
    def __init__(self, dfa: DFA, literalNames: list[str] | None = None, symbolicNames: list[str] | None = None) -> None: ...
    def getEdgeLabel(self, i: int): ...
    def getStateString(self, s: DFAState): ...

class LexerDFASerializer(DFASerializer):
    def __init__(self, dfa: DFA) -> None: ...
    def getEdgeLabel(self, i: int): ...
