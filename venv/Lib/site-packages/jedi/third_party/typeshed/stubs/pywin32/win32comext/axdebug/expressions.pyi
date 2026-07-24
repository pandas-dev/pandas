from _typeshed import Incomplete

from win32com.server.util import ListEnumeratorGateway
from win32comext.axdebug import gateways

def MakeNiceString(ob): ...

class ProvideExpressionContexts(gateways.ProvideExpressionContexts): ...

class ExpressionContext(gateways.DebugExpressionContext):
    frame: Incomplete
    def __init__(self, frame) -> None: ...
    def ParseLanguageText(self, code, radix, delim, flags): ...
    def GetLanguageInfo(self): ...

class Expression(gateways.DebugExpression):
    callback: Incomplete
    frame: Incomplete
    code: Incomplete
    radix: Incomplete
    delim: Incomplete
    flags: Incomplete
    isComplete: int
    result: Incomplete
    hresult: Incomplete
    def __init__(self, frame, code, radix, delim, flags) -> None: ...
    def Start(self, callback): ...
    def Abort(self) -> None: ...
    def QueryIsComplete(self): ...
    def GetResultAsString(self): ...
    def GetResultAsDebugProperty(self): ...

def MakeEnumDebugProperty(object, dwFieldSpec, nRadix, iid, stackFrame: Incomplete | None = ...): ...
def GetPropertyInfo(
    obname,
    obvalue,
    dwFieldSpec,
    nRadix,
    hresult: int = ...,
    dictionary: Incomplete | None = ...,
    stackFrame: Incomplete | None = ...,
): ...

class EnumDebugPropertyInfo(ListEnumeratorGateway):
    def GetCount(self): ...

class DebugProperty:
    name: Incomplete
    value: Incomplete
    parent: Incomplete
    hresult: Incomplete
    dictionary: Incomplete
    stackFrame: Incomplete
    def __init__(
        self,
        name,
        value,
        parent: Incomplete | None = ...,
        hresult: int = ...,
        dictionary: Incomplete | None = ...,
        stackFrame: Incomplete | None = ...,
    ) -> None: ...
    def GetPropertyInfo(self, dwFieldSpec, nRadix): ...
    def GetExtendedInfo(self) -> None: ...
    def SetValueAsString(self, value, radix) -> None: ...
    def EnumMembers(self, dwFieldSpec, nRadix, iid): ...
    def GetParent(self) -> None: ...
