from _typeshed import Incomplete
from collections.abc import Iterable
from typing import ClassVar, Final, Literal

class NotSupportedException(Exception): ...

DropIndirection: Final = "DropIndirection"
NoTranslateTypes: Final[list[int]]
NoTranslateMap: Final[set[int]]

class MapEntry:
    dispid: Incomplete
    desc: Incomplete | None
    names: Incomplete
    doc: Incomplete
    resultCLSID: Incomplete
    resultDocumentation: Incomplete | None
    wasProperty: Incomplete
    hidden: bool | Literal[0, 1]
    def __init__(
        self, desc_or_id, names=None, doc=None, resultCLSID=..., resultDoc=None, hidden: bool | Literal[0, 1] = 0
    ) -> None: ...
    def GetResultCLSID(self) -> Incomplete | None: ...
    def GetResultCLSIDStr(self) -> str: ...
    def GetResultName(self) -> Incomplete | None: ...

class OleItem:
    typename: ClassVar[str]
    doc: Incomplete
    python_name: Incomplete | None
    bWritten: bool | Literal[0, 1]
    bIsDispatch: bool | Literal[0, 1]
    bIsSink: bool | Literal[0, 1]
    clsid: Incomplete | None
    co_class: Incomplete | None
    def __init__(self, doc=None) -> None: ...

class DispatchItem(OleItem):
    typename: ClassVar[str]
    propMap: dict[Incomplete, Incomplete]
    propMapGet: dict[Incomplete, Incomplete]
    propMapPut: dict[Incomplete, Incomplete]
    mapFuncs: dict[Incomplete, Incomplete]
    defaultDispatchName: Incomplete | None
    hidden: bool | Literal[0, 1]
    def __init__(self, typeinfo=None, attr=None, doc=None, bForUser: bool | Literal[0, 1] = 1) -> None: ...
    clsid: Incomplete
    bIsDispatch: bool | Literal[0, 1]
    def Build(self, typeinfo, attr, bForUser: bool | Literal[0, 1] = 1) -> None: ...
    def CountInOutOptArgs(self, argTuple: Iterable[Incomplete]) -> tuple[int, int, int]: ...
    def MakeFuncMethod(self, entry, name: str, bMakeClass: bool | Literal[0, 1] = 1) -> list[str]: ...
    def MakeDispatchFuncMethod(self, entry, name: str, bMakeClass: bool | Literal[0, 1] = 1) -> list[str]: ...
    def MakeVarArgsFuncMethod(self, entry, name: str, bMakeClass: bool | Literal[0, 1] = 1) -> list[str]: ...

class VTableItem(DispatchItem):
    vtableFuncs: list[tuple[Incomplete, Incomplete, Incomplete]]
    def Build(self, typeinfo, attr, bForUser: bool | Literal[0, 1] = 1) -> None: ...

class LazyDispatchItem(DispatchItem):
    typename: ClassVar[str]
    clsid: Incomplete
    def __init__(self, attr, doc) -> None: ...

typeSubstMap: Final[dict[int, int]]
valid_identifier_chars: Final[str]

def demunge_leading_underscores(className: str) -> str: ...
def MakePublicAttributeName(className: str, is_global: bool = False) -> str: ...
def MakeDefaultArgRepr(defArgVal) -> str | None: ...
def BuildCallList(fdesc, names, defNamedOptArg, defNamedNotOptArg, defUnnamedArg, defOutArg, is_comment: bool = False) -> str: ...
