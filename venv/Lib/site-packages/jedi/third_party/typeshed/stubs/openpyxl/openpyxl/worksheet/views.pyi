from _typeshed import ConvertibleToFloat, ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Bool, Float, Integer, NoneSet, Set, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.sequence import Sequence
from openpyxl.descriptors.serialisable import Serialisable

_Pane: TypeAlias = Literal["bottomRight", "topRight", "bottomLeft", "topLeft"]
_SheetViewView: TypeAlias = Literal["normal", "pageBreakPreview", "pageLayout"]
_PaneState: TypeAlias = Literal["split", "frozen", "frozenSplit"]

class Pane(Serialisable):
    xSplit: Float[Literal[True]]
    ySplit: Float[Literal[True]]
    topLeftCell: String[Literal[True]]
    activePane: Set[_Pane]
    state: Set[_PaneState]
    def __init__(
        self,
        xSplit: ConvertibleToFloat | None = None,
        ySplit: ConvertibleToFloat | None = None,
        topLeftCell: str | None = None,
        activePane: _Pane = "topLeft",
        state: _PaneState = "split",
    ) -> None: ...

class Selection(Serialisable):
    pane: NoneSet[_Pane]
    activeCell: String[Literal[True]]
    activeCellId: Integer[Literal[True]]
    sqref: String[Literal[True]]
    def __init__(
        self,
        pane: _Pane | Literal["none"] | None = None,
        activeCell: str | None = "A1",
        activeCellId: ConvertibleToInt | None = None,
        sqref: str | None = "A1",
    ) -> None: ...

class SheetView(Serialisable):
    tagname: ClassVar[str]
    windowProtection: Bool[Literal[True]]
    showFormulas: Bool[Literal[True]]
    showGridLines: Bool[Literal[True]]
    showRowColHeaders: Bool[Literal[True]]
    showZeros: Bool[Literal[True]]
    rightToLeft: Bool[Literal[True]]
    tabSelected: Bool[Literal[True]]
    showRuler: Bool[Literal[True]]
    showOutlineSymbols: Bool[Literal[True]]
    defaultGridColor: Bool[Literal[True]]
    showWhiteSpace: Bool[Literal[True]]
    view: NoneSet[_SheetViewView]
    topLeftCell: String[Literal[True]]
    colorId: Integer[Literal[True]]
    zoomScale: Integer[Literal[True]]
    zoomScaleNormal: Integer[Literal[True]]
    zoomScaleSheetLayoutView: Integer[Literal[True]]
    zoomScalePageLayoutView: Integer[Literal[True]]
    zoomToFit: Bool[Literal[True]]
    workbookViewId: Integer[Literal[True]]
    selection: Incomplete
    pane: Typed[Pane, Literal[True]]
    def __init__(
        self,
        windowProtection: _ConvertibleToBool | None = None,
        showFormulas: _ConvertibleToBool | None = None,
        showGridLines: _ConvertibleToBool | None = None,
        showRowColHeaders: _ConvertibleToBool | None = None,
        showZeros: _ConvertibleToBool | None = None,
        rightToLeft: _ConvertibleToBool | None = None,
        tabSelected: _ConvertibleToBool | None = None,
        showRuler: _ConvertibleToBool | None = None,
        showOutlineSymbols: _ConvertibleToBool | None = None,
        defaultGridColor: _ConvertibleToBool | None = None,
        showWhiteSpace: _ConvertibleToBool | None = None,
        view: _SheetViewView | Literal["none"] | None = None,
        topLeftCell: str | None = None,
        colorId: ConvertibleToInt | None = None,
        zoomScale: ConvertibleToInt | None = None,
        zoomScaleNormal: ConvertibleToInt | None = None,
        zoomScaleSheetLayoutView: ConvertibleToInt | None = None,
        zoomScalePageLayoutView: ConvertibleToInt | None = None,
        zoomToFit: _ConvertibleToBool | None = None,
        workbookViewId: ConvertibleToInt | None = 0,
        selection=None,
        pane: Pane | None = None,
    ) -> None: ...

class SheetViewList(Serialisable):
    tagname: ClassVar[str]
    sheetView: Sequence[list[SheetView]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, sheetView: SheetView | None = None, extLst: Unused = None) -> None: ...
    @property
    def active(self) -> SheetView: ...
