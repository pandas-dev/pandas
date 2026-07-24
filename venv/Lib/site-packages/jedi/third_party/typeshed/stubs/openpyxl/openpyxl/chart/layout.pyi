from _typeshed import ConvertibleToFloat, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Alias, Typed
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedMinMax, NestedNoneSet, NestedSet, _NestedNoneSetParam
from openpyxl.descriptors.serialisable import Serialisable

from ..xml._functions_overloads import _HasTagAndGet

_ManualLayoutMode: TypeAlias = Literal["edge", "factor"]
_ManualLayoutLayoutTarget: TypeAlias = Literal["inner", "outer"]

class ManualLayout(Serialisable):
    tagname: ClassVar[str]
    layoutTarget: NestedNoneSet[_ManualLayoutLayoutTarget]
    xMode: NestedNoneSet[_ManualLayoutMode]
    yMode: NestedNoneSet[_ManualLayoutMode]
    wMode: NestedSet[_ManualLayoutMode]
    hMode: NestedSet[_ManualLayoutMode]
    x: NestedMinMax[float, Literal[True]]
    y: NestedMinMax[float, Literal[True]]
    w: NestedMinMax[float, Literal[True]]
    width: Alias
    h: NestedMinMax[float, Literal[True]]
    height: Alias
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        layoutTarget: _NestedNoneSetParam[_ManualLayoutLayoutTarget] = None,
        xMode: _NestedNoneSetParam[_ManualLayoutMode] = None,
        yMode: _NestedNoneSetParam[_ManualLayoutMode] = None,
        wMode: _HasTagAndGet[_ManualLayoutMode] | _ManualLayoutMode = "factor",
        hMode: _HasTagAndGet[_ManualLayoutMode] | _ManualLayoutMode = "factor",
        x: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        y: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        w: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        h: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        extLst: Unused = None,
    ) -> None: ...

class Layout(Serialisable):
    tagname: ClassVar[str]
    manualLayout: Typed[ManualLayout, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, manualLayout: ManualLayout | None = None, extLst: Unused = None) -> None: ...
