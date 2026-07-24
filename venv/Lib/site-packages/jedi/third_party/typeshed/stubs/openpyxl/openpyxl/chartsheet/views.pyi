from _typeshed import ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal

from openpyxl.descriptors.base import Bool, Integer, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.serialisable import Serialisable

class ChartsheetView(Serialisable):
    tagname: ClassVar[str]
    tabSelected: Bool[Literal[True]]
    zoomScale: Integer[Literal[True]]
    workbookViewId: Integer[Literal[False]]
    zoomToFit: Bool[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        tabSelected: _ConvertibleToBool | None = None,
        zoomScale: ConvertibleToInt | None = None,
        workbookViewId: ConvertibleToInt = 0,
        zoomToFit: _ConvertibleToBool | None = True,
        extLst: Unused = None,
    ) -> None: ...

class ChartsheetViewList(Serialisable):
    tagname: ClassVar[str]
    sheetView: Incomplete
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, sheetView=None, extLst: Unused = None) -> None: ...
