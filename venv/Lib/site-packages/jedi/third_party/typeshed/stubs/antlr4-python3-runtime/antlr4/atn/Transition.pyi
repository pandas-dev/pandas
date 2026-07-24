from _typeshed import Incomplete

from antlr4.atn.ATNState import RuleStartState
from antlr4.IntervalSet import IntervalSet

class Transition:
    __slots__ = ("target", "isEpsilon", "label")
    EPSILON: int
    RANGE: int
    RULE: int
    PREDICATE: int
    ATOM: int
    ACTION: int
    SET: int
    NOT_SET: int
    WILDCARD: int
    PRECEDENCE: int
    serializationNames: Incomplete
    serializationTypes: Incomplete
    target: Incomplete
    isEpsilon: bool
    label: Incomplete
    def __init__(self, target) -> None: ...

class AtomTransition(Transition):
    __slots__ = ("label_", "serializationType")
    label_: Incomplete
    label: Incomplete
    serializationType: Incomplete
    def __init__(self, target, label: int) -> None: ...
    def makeLabel(self): ...
    def matches(self, symbol: int, minVocabSymbol: int, maxVocabSymbol: int): ...

class RuleTransition(Transition):
    __slots__ = ("ruleIndex", "precedence", "followState", "serializationType")
    ruleIndex: Incomplete
    precedence: Incomplete
    followState: Incomplete
    serializationType: Incomplete
    isEpsilon: bool
    def __init__(self, ruleStart: RuleStartState, ruleIndex: int, precedence: int, followState) -> None: ...
    def matches(self, symbol: int, minVocabSymbol: int, maxVocabSymbol: int): ...

class EpsilonTransition(Transition):
    __slots__ = ("serializationType", "outermostPrecedenceReturn")
    serializationType: Incomplete
    isEpsilon: bool
    outermostPrecedenceReturn: Incomplete
    def __init__(self, target, outermostPrecedenceReturn: int = -1) -> None: ...
    def matches(self, symbol: int, minVocabSymbol: int, maxVocabSymbol: int): ...

class RangeTransition(Transition):
    __slots__ = ("serializationType", "start", "stop")
    serializationType: Incomplete
    start: Incomplete
    stop: Incomplete
    label: Incomplete
    def __init__(self, target, start: int, stop: int) -> None: ...
    def makeLabel(self): ...
    def matches(self, symbol: int, minVocabSymbol: int, maxVocabSymbol: int): ...

class AbstractPredicateTransition(Transition):
    def __init__(self, target) -> None: ...

class PredicateTransition(AbstractPredicateTransition):
    __slots__ = ("serializationType", "ruleIndex", "predIndex", "isCtxDependent")
    serializationType: Incomplete
    ruleIndex: Incomplete
    predIndex: Incomplete
    isCtxDependent: Incomplete
    isEpsilon: bool
    def __init__(self, target, ruleIndex: int, predIndex: int, isCtxDependent: bool) -> None: ...
    def matches(self, symbol: int, minVocabSymbol: int, maxVocabSymbol: int): ...
    def getPredicate(self): ...

class ActionTransition(Transition):
    __slots__ = ("serializationType", "ruleIndex", "actionIndex", "isCtxDependent")
    serializationType: Incomplete
    ruleIndex: Incomplete
    actionIndex: Incomplete
    isCtxDependent: Incomplete
    isEpsilon: bool
    def __init__(self, target, ruleIndex: int, actionIndex: int = -1, isCtxDependent: bool = False) -> None: ...
    def matches(self, symbol: int, minVocabSymbol: int, maxVocabSymbol: int): ...

class SetTransition(Transition):
    __slots__ = "serializationType"
    serializationType: Incomplete
    label: Incomplete
    def __init__(self, target, set: IntervalSet) -> None: ...
    def matches(self, symbol: int, minVocabSymbol: int, maxVocabSymbol: int): ...

class NotSetTransition(SetTransition):
    serializationType: Incomplete
    def __init__(self, target, set: IntervalSet) -> None: ...
    def matches(self, symbol: int, minVocabSymbol: int, maxVocabSymbol: int): ...

class WildcardTransition(Transition):
    __slots__ = "serializationType"
    serializationType: Incomplete
    def __init__(self, target) -> None: ...
    def matches(self, symbol: int, minVocabSymbol: int, maxVocabSymbol: int): ...

class PrecedencePredicateTransition(AbstractPredicateTransition):
    __slots__ = ("serializationType", "precedence")
    serializationType: Incomplete
    precedence: Incomplete
    isEpsilon: bool
    def __init__(self, target, precedence: int) -> None: ...
    def matches(self, symbol: int, minVocabSymbol: int, maxVocabSymbol: int): ...
    def getPredicate(self): ...
