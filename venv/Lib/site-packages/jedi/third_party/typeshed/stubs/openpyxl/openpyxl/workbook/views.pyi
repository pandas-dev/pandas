from _typeshed import ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal, overload
from typing_extensions import TypeAlias

from openpyxl import _VisibilityType
from openpyxl.descriptors.base import Bool, Integer, NoneSet, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.serialisable import Serialisable

_CustomWorkbookViewShowComments: TypeAlias = Literal["commNone", "commIndicator", "commIndAndComment"]
_CustomWorkbookViewShowObjects: TypeAlias = Literal["all", "placeholders"]

class BookView(Serialisable):
    tagname: ClassVar[str]
    visibility: NoneSet[_VisibilityType]
    minimized: Bool[Literal[True]]
    showHorizontalScroll: Bool[Literal[True]]
    showVerticalScroll: Bool[Literal[True]]
    showSheetTabs: Bool[Literal[True]]
    xWindow: Integer[Literal[True]]
    yWindow: Integer[Literal[True]]
    windowWidth: Integer[Literal[True]]
    windowHeight: Integer[Literal[True]]
    tabRatio: Integer[Literal[True]]
    firstSheet: Integer[Literal[True]]
    activeTab: Integer[Literal[True]]
    autoFilterDateGrouping: Bool[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        visibility: _VisibilityType | Literal["none"] | None = "visible",
        minimized: _ConvertibleToBool | None = False,
        showHorizontalScroll: _ConvertibleToBool | None = True,
        showVerticalScroll: _ConvertibleToBool | None = True,
        showSheetTabs: _ConvertibleToBool | None = True,
        xWindow: ConvertibleToInt | None = None,
        yWindow: ConvertibleToInt | None = None,
        windowWidth: ConvertibleToInt | None = None,
        windowHeight: ConvertibleToInt | None = None,
        tabRatio: ConvertibleToInt | None = 600,
        firstSheet: ConvertibleToInt | None = 0,
        activeTab: ConvertibleToInt | None = 0,
        autoFilterDateGrouping: _ConvertibleToBool | None = True,
        extLst: Unused = None,
    ) -> None: ...

class CustomWorkbookView(Serialisable):
    tagname: ClassVar[str]
    name: String[Literal[False]]
    guid: Incomplete
    autoUpdate: Bool[Literal[True]]
    mergeInterval: Integer[Literal[True]]
    changesSavedWin: Bool[Literal[True]]
    onlySync: Bool[Literal[True]]
    personalView: Bool[Literal[True]]
    includePrintSettings: Bool[Literal[True]]
    includeHiddenRowCol: Bool[Literal[True]]
    maximized: Bool[Literal[True]]
    minimized: Bool[Literal[True]]
    showHorizontalScroll: Bool[Literal[True]]
    showVerticalScroll: Bool[Literal[True]]
    showSheetTabs: Bool[Literal[True]]
    xWindow: Integer[Literal[False]]
    yWindow: Integer[Literal[False]]
    windowWidth: Integer[Literal[False]]
    windowHeight: Integer[Literal[False]]
    tabRatio: Integer[Literal[True]]
    activeSheetId: Integer[Literal[False]]
    showFormulaBar: Bool[Literal[True]]
    showStatusbar: Bool[Literal[True]]
    showComments: NoneSet[_CustomWorkbookViewShowComments]
    showObjects: NoneSet[_CustomWorkbookViewShowObjects]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        name: str,
        guid=None,
        autoUpdate: _ConvertibleToBool | None = None,
        mergeInterval: ConvertibleToInt | None = None,
        changesSavedWin: _ConvertibleToBool | None = None,
        onlySync: _ConvertibleToBool | None = None,
        personalView: _ConvertibleToBool | None = None,
        includePrintSettings: _ConvertibleToBool | None = None,
        includeHiddenRowCol: _ConvertibleToBool | None = None,
        maximized: _ConvertibleToBool | None = None,
        minimized: _ConvertibleToBool | None = None,
        showHorizontalScroll: _ConvertibleToBool | None = None,
        showVerticalScroll: _ConvertibleToBool | None = None,
        showSheetTabs: _ConvertibleToBool | None = None,
        *,
        xWindow: ConvertibleToInt,
        yWindow: ConvertibleToInt,
        windowWidth: ConvertibleToInt,
        windowHeight: ConvertibleToInt,
        tabRatio: ConvertibleToInt | None = None,
        activeSheetId: ConvertibleToInt,
        showFormulaBar: _ConvertibleToBool | None = None,
        showStatusbar: _ConvertibleToBool | None = None,
        showComments: _CustomWorkbookViewShowComments | Literal["none"] | None = "commIndicator",
        showObjects: _CustomWorkbookViewShowObjects | Literal["none"] | None = "all",
        extLst: Unused = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        name: str,
        guid: Incomplete | None,
        autoUpdate: _ConvertibleToBool | None,
        mergeInterval: ConvertibleToInt | None,
        changesSavedWin: _ConvertibleToBool | None,
        onlySync: _ConvertibleToBool | None,
        personalView: _ConvertibleToBool | None,
        includePrintSettings: _ConvertibleToBool | None,
        includeHiddenRowCol: _ConvertibleToBool | None,
        maximized: _ConvertibleToBool | None,
        minimized: _ConvertibleToBool | None,
        showHorizontalScroll: _ConvertibleToBool | None,
        showVerticalScroll: _ConvertibleToBool | None,
        showSheetTabs: _ConvertibleToBool | None,
        xWindow: ConvertibleToInt,
        yWindow: ConvertibleToInt,
        windowWidth: ConvertibleToInt,
        windowHeight: ConvertibleToInt,
        tabRatio: ConvertibleToInt | None,
        activeSheetId: ConvertibleToInt,
        showFormulaBar: _ConvertibleToBool | None = None,
        showStatusbar: _ConvertibleToBool | None = None,
        showComments: _CustomWorkbookViewShowComments | Literal["none"] | None = "commIndicator",
        showObjects: _CustomWorkbookViewShowObjects | Literal["none"] | None = "all",
        extLst: Unused = None,
    ) -> None: ...
