from _typeshed import Incomplete

class ATNDeserializationOptions:
    __slots__ = ("readonly", "verifyATN", "generateRuleBypassTransitions")
    defaultOptions: Incomplete
    readonly: bool
    verifyATN: Incomplete
    generateRuleBypassTransitions: Incomplete
    def __init__(self, copyFrom: ATNDeserializationOptions | None = None) -> None: ...
    def __setattr__(self, key, value) -> None: ...
