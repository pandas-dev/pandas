from _typeshed import ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal, NoReturn, overload

from openpyxl.descriptors import Strict
from openpyxl.descriptors.base import Alias, Bool, Integer, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedInteger, NestedText
from openpyxl.descriptors.serialisable import Serialisable

from ..xml._functions_overloads import _HasTagAndGet

class NumFmt(Serialisable):
    formatCode: String[Literal[False]]
    sourceLinked: Bool[Literal[False]]
    def __init__(self, formatCode: str, sourceLinked: _ConvertibleToBool = False) -> None: ...

class NumberValueDescriptor(NestedText[Incomplete, Incomplete]):
    allow_none: bool
    expected_type: type[Incomplete]
    def __set__(self, instance: Serialisable | Strict, value) -> None: ...  # type: ignore[override]

class NumVal(Serialisable):
    idx: Integer[Literal[False]]
    formatCode: NestedText[str, Literal[True]]
    v: Incomplete
    def __init__(self, idx: ConvertibleToInt, formatCode: object = None, v=None) -> None: ...

class NumData(Serialisable):
    formatCode: NestedText[str, Literal[True]]
    ptCount: NestedInteger[Literal[True]]
    pt: Incomplete
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        formatCode: object = None,
        ptCount: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        pt=(),
        extLst: Unused = None,
    ) -> None: ...

class NumRef(Serialisable):
    f: NestedText[str, Literal[False]]
    ref: Alias
    numCache: Typed[NumData, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, f: object = None, numCache: NumData | None = None, extLst: Unused = None) -> None: ...

class StrVal(Serialisable):
    tagname: ClassVar[str]
    idx: Integer[Literal[False]]
    v: NestedText[str, Literal[False]]
    def __init__(self, idx: ConvertibleToInt = 0, v: object = None) -> None: ...

class StrData(Serialisable):
    tagname: ClassVar[str]
    ptCount: NestedInteger[Literal[True]]
    pt: Incomplete
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self, ptCount: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None, pt=(), extLst: Unused = None
    ) -> None: ...

class StrRef(Serialisable):
    tagname: ClassVar[str]
    f: NestedText[str, Literal[True]]
    strCache: Typed[StrData, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, f: object = None, strCache: StrData | None = None, extLst: Unused = None) -> None: ...

class NumDataSource(Serialisable):
    numRef: Typed[NumRef, Literal[True]]
    numLit: Typed[NumData, Literal[True]]
    def __init__(self, numRef: NumRef | None = None, numLit: NumData | None = None) -> None: ...

class Level(Serialisable):
    tagname: ClassVar[str]
    pt: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, pt=()) -> None: ...

class MultiLevelStrData(Serialisable):
    tagname: ClassVar[str]
    ptCount: Integer[Literal[True]]
    lvl: Incomplete
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, ptCount: ConvertibleToInt | None = None, lvl=(), extLst: Unused = None) -> None: ...

class MultiLevelStrRef(Serialisable):
    tagname: ClassVar[str]
    f: NestedText[str, Literal[False]]
    multiLvlStrCache: Typed[MultiLevelStrData, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, f: object = None, multiLvlStrCache: MultiLevelStrData | None = None, extLst: Unused = None) -> None: ...

class AxDataSource(Serialisable):
    tagname: ClassVar[str]
    numRef: Typed[NumRef, Literal[True]]
    numLit: Typed[NumData, Literal[True]]
    strRef: Typed[StrRef, Literal[True]]
    strLit: Typed[StrData, Literal[True]]
    multiLvlStrRef: Typed[MultiLevelStrRef, Literal[True]]
    @overload
    def __init__(
        self, numRef: None = None, numLit: None = None, strRef: None = None, strLit: None = None, multiLvlStrRef: None = None
    ) -> NoReturn: ...
    @overload
    def __init__(
        self,
        numRef: NumRef | None = None,
        numLit: NumData | None = None,
        strRef: StrRef | None = None,
        strLit: StrData | None = None,
        multiLvlStrRef: MultiLevelStrRef | None = None,
    ) -> None: ...
