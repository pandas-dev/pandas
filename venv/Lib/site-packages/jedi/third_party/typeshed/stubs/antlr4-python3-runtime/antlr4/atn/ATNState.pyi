from _typeshed import Incomplete

from antlr4.atn.Transition import Transition as Transition

INITIAL_NUM_TRANSITIONS: int

class ATNState:
    __slots__ = ("atn", "stateNumber", "stateType", "ruleIndex", "epsilonOnlyTransitions", "transitions", "nextTokenWithinRule")
    INVALID_TYPE: int
    BASIC: int
    RULE_START: int
    BLOCK_START: int
    PLUS_BLOCK_START: int
    STAR_BLOCK_START: int
    TOKEN_START: int
    RULE_STOP: int
    BLOCK_END: int
    STAR_LOOP_BACK: int
    STAR_LOOP_ENTRY: int
    PLUS_LOOP_BACK: int
    LOOP_END: int
    serializationNames: Incomplete
    INVALID_STATE_NUMBER: int
    atn: Incomplete
    stateNumber: Incomplete
    stateType: Incomplete
    ruleIndex: int
    epsilonOnlyTransitions: bool
    transitions: Incomplete
    nextTokenWithinRule: Incomplete
    def __init__(self) -> None: ...
    def __hash__(self): ...
    def __eq__(self, other): ...
    def onlyHasEpsilonTransitions(self): ...
    def isNonGreedyExitState(self): ...
    def addTransition(self, trans: Transition, index: int = -1): ...

class BasicState(ATNState):
    stateType: Incomplete
    def __init__(self) -> None: ...

class DecisionState(ATNState):
    __slots__ = ("decision", "nonGreedy")
    decision: int
    nonGreedy: bool
    def __init__(self) -> None: ...

class BlockStartState(DecisionState):
    __slots__ = "endState"
    endState: Incomplete
    def __init__(self) -> None: ...

class BasicBlockStartState(BlockStartState):
    stateType: Incomplete
    def __init__(self) -> None: ...

class BlockEndState(ATNState):
    __slots__ = "startState"
    stateType: Incomplete
    startState: Incomplete
    def __init__(self) -> None: ...

class RuleStopState(ATNState):
    stateType: Incomplete
    def __init__(self) -> None: ...

class RuleStartState(ATNState):
    __slots__ = ("stopState", "isPrecedenceRule")
    stateType: Incomplete
    stopState: Incomplete
    isPrecedenceRule: bool
    def __init__(self) -> None: ...

class PlusLoopbackState(DecisionState):
    stateType: Incomplete
    def __init__(self) -> None: ...

class PlusBlockStartState(BlockStartState):
    __slots__ = "loopBackState"
    stateType: Incomplete
    loopBackState: Incomplete
    def __init__(self) -> None: ...

class StarBlockStartState(BlockStartState):
    stateType: Incomplete
    def __init__(self) -> None: ...

class StarLoopbackState(ATNState):
    stateType: Incomplete
    def __init__(self) -> None: ...

class StarLoopEntryState(DecisionState):
    __slots__ = ("loopBackState", "isPrecedenceDecision")
    stateType: Incomplete
    loopBackState: Incomplete
    isPrecedenceDecision: Incomplete
    def __init__(self) -> None: ...

class LoopEndState(ATNState):
    __slots__ = "loopBackState"
    stateType: Incomplete
    loopBackState: Incomplete
    def __init__(self) -> None: ...

class TokensStartState(DecisionState):
    stateType: Incomplete
    def __init__(self) -> None: ...
