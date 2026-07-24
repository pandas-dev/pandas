from _typeshed import ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal, overload
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Bool, Integer, NoneSet, Set, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.worksheet.filters import AutoFilter
from openpyxl.xml.functions import Element

_PivotAreaType: TypeAlias = Literal["normal", "data", "all", "origin", "button", "topEnd", "topRight"]
_PivotAxis: TypeAlias = Literal["axisRow", "axisCol", "axisPage", "axisValues"]
_ConditionalFormatType: TypeAlias = Literal["all", "row", "column"]
_FormatAction: TypeAlias = Literal["blank", "formatting", "drill", "formula"]
_PivotFilterType: TypeAlias = Literal[
    "unknown",
    "count",
    "percent",
    "sum",
    "captionEqual",
    "captionNotEqual",
    "captionBeginsWith",
    "captionNotBeginsWith",
    "captionEndsWith",
    "captionNotEndsWith",
    "captionContains",
    "captionNotContains",
    "captionGreaterThan",
    "captionGreaterThanOrEqual",
    "captionLessThan",
    "captionLessThanOrEqual",
    "captionBetween",
    "captionNotBetween",
    "valueEqual",
    "valueNotEqual",
    "valueGreaterThan",
    "valueGreaterThanOrEqual",
    "valueLessThan",
    "valueLessThanOrEqual",
    "valueBetween",
    "valueNotBetween",
    "dateEqual",
    "dateNotEqual",
    "dateOlderThan",
    "dateOlderThanOrEqual",
    "dateNewerThan",
    "dateNewerThanOrEqual",
    "dateBetween",
    "dateNotBetween",
    "tomorrow",
    "today",
    "yesterday",
    "nextWeek",
    "thisWeek",
    "lastWeek",
    "nextMonth",
    "thisMonth",
    "lastMonth",
    "nextQuarter",
    "thisQuarter",
    "lastQuarter",
    "nextYear",
    "thisYear",
    "lastYear",
    "yearToDate",
    "Q1",
    "Q2",
    "Q3",
    "Q4",
    "M1",
    "M2",
    "M3",
    "M4",
    "M5",
    "M6",
    "M7",
    "M8",
    "M9",
    "M10",
    "M11",
    "M12",
]
_ConditionalFormatScope: TypeAlias = Literal["selection", "data", "field"]
_DataFieldSubtotal: TypeAlias = Literal[
    "average", "count", "countNums", "max", "min", "product", "stdDev", "stdDevp", "sum", "var", "varp"
]
_DataFieldShowDataAs: TypeAlias = Literal[
    "normal", "difference", "percent", "percentDiff", "runTotal", "percentOfRow", "percentOfCol", "percentOfTotal", "index"
]
_ItemType: TypeAlias = Literal[
    "data",
    "default",
    "sum",
    "countA",
    "avg",
    "max",
    "min",
    "product",
    "count",
    "stdDev",
    "stdDevP",
    "var",
    "varP",
    "grand",
    "blank",
]
_PivotFieldSortType: TypeAlias = Literal["manual", "ascending", "descending"]

class HierarchyUsage(Serialisable):
    tagname: ClassVar[str]
    hierarchyUsage: Integer[Literal[False]]
    def __init__(self, hierarchyUsage: ConvertibleToInt) -> None: ...

class ColHierarchiesUsage(Serialisable):
    tagname: ClassVar[str]
    colHierarchyUsage: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    __attrs__: ClassVar[tuple[str, ...]]
    def __init__(self, count=None, colHierarchyUsage=()) -> None: ...
    @property
    def count(self) -> int: ...

class RowHierarchiesUsage(Serialisable):
    tagname: ClassVar[str]
    rowHierarchyUsage: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    __attrs__: ClassVar[tuple[str, ...]]
    def __init__(self, count=None, rowHierarchyUsage=()) -> None: ...
    @property
    def count(self) -> int: ...

class PivotFilter(Serialisable):
    tagname: ClassVar[str]
    fld: Integer[Literal[False]]
    mpFld: Integer[Literal[True]]
    type: Set[_PivotFilterType]
    evalOrder: Integer[Literal[True]]
    id: Integer[Literal[False]]
    iMeasureHier: Integer[Literal[True]]
    iMeasureFld: Integer[Literal[True]]
    name: String[Literal[True]]
    description: String[Literal[True]]
    stringValue1: String[Literal[True]]
    stringValue2: String[Literal[True]]
    autoFilter: Typed[AutoFilter, Literal[False]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        fld: ConvertibleToInt,
        mpFld: ConvertibleToInt | None = None,
        *,
        type: _PivotFilterType,
        evalOrder: ConvertibleToInt | None = None,
        id: ConvertibleToInt,
        iMeasureHier: ConvertibleToInt | None = None,
        iMeasureFld: ConvertibleToInt | None = None,
        name: str | None = None,
        description: str | None = None,
        stringValue1: str | None = None,
        stringValue2: str | None = None,
        autoFilter: AutoFilter,
        extLst: ExtensionList | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        fld: ConvertibleToInt,
        mpFld: ConvertibleToInt | None,
        type: _PivotFilterType,
        evalOrder: ConvertibleToInt | None,
        id: ConvertibleToInt,
        iMeasureHier: ConvertibleToInt | None,
        iMeasureFld: ConvertibleToInt | None,
        name: str | None,
        description: str | None,
        stringValue1: str | None,
        stringValue2: str | None,
        autoFilter: AutoFilter,
        extLst: ExtensionList | None = None,
    ) -> None: ...

class PivotFilters(Serialisable):
    count: Integer[Literal[False]]
    filter: Typed[PivotFilter, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, count: ConvertibleToInt, filter: PivotFilter | None = None) -> None: ...

class PivotTableStyle(Serialisable):
    tagname: ClassVar[str]
    name: String[Literal[True]]
    showRowHeaders: Bool[Literal[False]]
    showColHeaders: Bool[Literal[False]]
    showRowStripes: Bool[Literal[False]]
    showColStripes: Bool[Literal[False]]
    showLastColumn: Bool[Literal[False]]
    def __init__(
        self,
        name: str | None = None,
        showRowHeaders: _ConvertibleToBool = None,
        showColHeaders: _ConvertibleToBool = None,
        showRowStripes: _ConvertibleToBool = None,
        showColStripes: _ConvertibleToBool = None,
        showLastColumn: _ConvertibleToBool = None,
    ) -> None: ...

class MemberList(Serialisable):
    tagname: ClassVar[str]
    level: Integer[Literal[True]]
    member: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, count=None, level: ConvertibleToInt | None = None, member=()) -> None: ...
    @property
    def count(self) -> int: ...

class MemberProperty(Serialisable):
    tagname: ClassVar[str]
    name: String[Literal[True]]
    showCell: Bool[Literal[True]]
    showTip: Bool[Literal[True]]
    showAsCaption: Bool[Literal[True]]
    nameLen: Integer[Literal[True]]
    pPos: Integer[Literal[True]]
    pLen: Integer[Literal[True]]
    level: Integer[Literal[True]]
    field: Integer[Literal[False]]
    @overload
    def __init__(
        self,
        name: str | None = None,
        showCell: _ConvertibleToBool | None = None,
        showTip: _ConvertibleToBool | None = None,
        showAsCaption: _ConvertibleToBool | None = None,
        nameLen: ConvertibleToInt | None = None,
        pPos: ConvertibleToInt | None = None,
        pLen: ConvertibleToInt | None = None,
        level: ConvertibleToInt | None = None,
        *,
        field: ConvertibleToInt,
    ) -> None: ...
    @overload
    def __init__(
        self,
        name: str | None,
        showCell: _ConvertibleToBool | None,
        showTip: _ConvertibleToBool | None,
        showAsCaption: _ConvertibleToBool | None,
        nameLen: ConvertibleToInt | None,
        pPos: ConvertibleToInt | None,
        pLen: ConvertibleToInt | None,
        level: ConvertibleToInt | None,
        field: ConvertibleToInt,
    ) -> None: ...

class PivotHierarchy(Serialisable):
    tagname: ClassVar[str]
    outline: Bool[Literal[False]]
    multipleItemSelectionAllowed: Bool[Literal[False]]
    subtotalTop: Bool[Literal[False]]
    showInFieldList: Bool[Literal[False]]
    dragToRow: Bool[Literal[False]]
    dragToCol: Bool[Literal[False]]
    dragToPage: Bool[Literal[False]]
    dragToData: Bool[Literal[False]]
    dragOff: Bool[Literal[False]]
    includeNewItemsInFilter: Bool[Literal[False]]
    caption: String[Literal[True]]
    mps: Incomplete
    members: Typed[MemberList, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        outline: _ConvertibleToBool = None,
        multipleItemSelectionAllowed: _ConvertibleToBool = None,
        subtotalTop: _ConvertibleToBool = None,
        showInFieldList: _ConvertibleToBool = None,
        dragToRow: _ConvertibleToBool = None,
        dragToCol: _ConvertibleToBool = None,
        dragToPage: _ConvertibleToBool = None,
        dragToData: _ConvertibleToBool = None,
        dragOff: _ConvertibleToBool = None,
        includeNewItemsInFilter: _ConvertibleToBool = None,
        caption: str | None = None,
        mps=(),
        members: MemberList | None = None,
        extLst: ExtensionList | None = None,
    ) -> None: ...

class Reference(Serialisable):
    tagname: ClassVar[str]
    field: Integer[Literal[True]]
    selected: Bool[Literal[True]]
    byPosition: Bool[Literal[True]]
    relative: Bool[Literal[True]]
    defaultSubtotal: Bool[Literal[True]]
    sumSubtotal: Bool[Literal[True]]
    countASubtotal: Bool[Literal[True]]
    avgSubtotal: Bool[Literal[True]]
    maxSubtotal: Bool[Literal[True]]
    minSubtotal: Bool[Literal[True]]
    productSubtotal: Bool[Literal[True]]
    countSubtotal: Bool[Literal[True]]
    stdDevSubtotal: Bool[Literal[True]]
    stdDevPSubtotal: Bool[Literal[True]]
    varSubtotal: Bool[Literal[True]]
    varPSubtotal: Bool[Literal[True]]
    x: Incomplete
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        field: ConvertibleToInt | None = None,
        count: Unused = None,
        selected: _ConvertibleToBool | None = None,
        byPosition: _ConvertibleToBool | None = None,
        relative: _ConvertibleToBool | None = None,
        defaultSubtotal: _ConvertibleToBool | None = None,
        sumSubtotal: _ConvertibleToBool | None = None,
        countASubtotal: _ConvertibleToBool | None = None,
        avgSubtotal: _ConvertibleToBool | None = None,
        maxSubtotal: _ConvertibleToBool | None = None,
        minSubtotal: _ConvertibleToBool | None = None,
        productSubtotal: _ConvertibleToBool | None = None,
        countSubtotal: _ConvertibleToBool | None = None,
        stdDevSubtotal: _ConvertibleToBool | None = None,
        stdDevPSubtotal: _ConvertibleToBool | None = None,
        varSubtotal: _ConvertibleToBool | None = None,
        varPSubtotal: _ConvertibleToBool | None = None,
        x: Incomplete | None = (),
        extLst: ExtensionList | None = None,
    ) -> None: ...
    @property
    def count(self) -> int: ...

class PivotArea(Serialisable):
    tagname: ClassVar[str]
    references: Incomplete
    extLst: Typed[ExtensionList, Literal[True]]
    field: Integer[Literal[True]]
    type: NoneSet[_PivotAreaType]
    dataOnly: Bool[Literal[True]]
    labelOnly: Bool[Literal[True]]
    grandRow: Bool[Literal[True]]
    grandCol: Bool[Literal[True]]
    cacheIndex: Bool[Literal[True]]
    outline: Bool[Literal[True]]
    offset: String[Literal[True]]
    collapsedLevelsAreSubtotals: Bool[Literal[True]]
    axis: NoneSet[_PivotAxis]
    fieldPosition: Integer[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        references=(),
        extLst: ExtensionList | None = None,
        field: ConvertibleToInt | None = None,
        type: _PivotAreaType | Literal["none"] | None = "normal",
        dataOnly: _ConvertibleToBool | None = True,
        labelOnly: _ConvertibleToBool | None = None,
        grandRow: _ConvertibleToBool | None = None,
        grandCol: _ConvertibleToBool | None = None,
        cacheIndex: _ConvertibleToBool | None = None,
        outline: _ConvertibleToBool | None = True,
        offset: str | None = None,
        collapsedLevelsAreSubtotals: _ConvertibleToBool | None = None,
        axis: _PivotAxis | Literal["none"] | None = None,
        fieldPosition: ConvertibleToInt | None = None,
    ) -> None: ...

class ChartFormat(Serialisable):
    tagname: ClassVar[str]
    chart: Integer[Literal[False]]
    format: Integer[Literal[False]]
    series: Bool[Literal[False]]
    pivotArea: Typed[PivotArea, Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self, chart: ConvertibleToInt, format: ConvertibleToInt, series: _ConvertibleToBool = None, *, pivotArea: PivotArea
    ) -> None: ...
    @overload
    def __init__(
        self, chart: ConvertibleToInt, format: ConvertibleToInt, series: _ConvertibleToBool, pivotArea: PivotArea
    ) -> None: ...

class ConditionalFormat(Serialisable):
    tagname: ClassVar[str]
    scope: Set[_ConditionalFormatScope]
    type: NoneSet[_ConditionalFormatType]
    priority: Integer[Literal[False]]
    pivotAreas: Incomplete
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        scope: _ConditionalFormatScope = "selection",
        type: _ConditionalFormatType | Literal["none"] | None = None,
        *,
        priority: ConvertibleToInt,
        pivotAreas=(),
        extLst: ExtensionList | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        scope: _ConditionalFormatScope,
        type: _ConditionalFormatType | Literal["none"] | None,
        priority: ConvertibleToInt,
        pivotAreas=(),
        extLst: ExtensionList | None = None,
    ) -> None: ...

class ConditionalFormatList(Serialisable):
    tagname: ClassVar[str]
    conditionalFormat: Incomplete
    __attrs__: ClassVar[tuple[str, ...]]
    def __init__(self, conditionalFormat=(), count=None) -> None: ...
    def by_priority(self): ...
    @property
    def count(self) -> int: ...
    def to_tree(self, tagname: str | None = None) -> Element: ...  # type: ignore[override]

class Format(Serialisable):
    tagname: ClassVar[str]
    action: NoneSet[_FormatAction]
    dxfId: Integer[Literal[True]]
    pivotArea: Typed[PivotArea, Literal[False]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        action: _FormatAction | Literal["none"] | None = "formatting",
        dxfId: ConvertibleToInt | None = None,
        *,
        pivotArea: PivotArea,
        extLst: ExtensionList | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        action: _FormatAction | Literal["none"] | None,
        dxfId: ConvertibleToInt | None,
        pivotArea: PivotArea,
        extLst: ExtensionList | None = None,
    ) -> None: ...

class DataField(Serialisable):
    tagname: ClassVar[str]
    name: String[Literal[True]]
    fld: Integer[Literal[False]]
    subtotal: Set[_DataFieldSubtotal]
    showDataAs: Set[_DataFieldShowDataAs]
    baseField: Integer[Literal[False]]
    baseItem: Integer[Literal[False]]
    numFmtId: Integer[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        name: str | None = None,
        *,
        fld: ConvertibleToInt,
        subtotal: str = "sum",
        showDataAs: str = "normal",
        baseField: ConvertibleToInt = -1,
        baseItem: ConvertibleToInt = 1048832,
        numFmtId: ConvertibleToInt | None = None,
        extLst: ExtensionList | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        name: str | None,
        fld: ConvertibleToInt,
        subtotal: str = "sum",
        showDataAs: str = "normal",
        baseField: ConvertibleToInt = -1,
        baseItem: ConvertibleToInt = 1048832,
        numFmtId: ConvertibleToInt | None = None,
        extLst: ExtensionList | None = None,
    ) -> None: ...

class PageField(Serialisable):
    tagname: ClassVar[str]
    fld: Integer[Literal[False]]
    item: Integer[Literal[True]]
    hier: Integer[Literal[True]]
    name: String[Literal[True]]
    cap: String[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        fld: ConvertibleToInt,
        item: ConvertibleToInt | None = None,
        hier: ConvertibleToInt | None = None,
        name: str | None = None,
        cap: str | None = None,
        extLst: ExtensionList | None = None,
    ) -> None: ...

class RowColItem(Serialisable):
    tagname: ClassVar[str]
    t: Set[_ItemType]
    r: Integer[Literal[False]]
    i: Integer[Literal[False]]
    x: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, t: _ItemType = "data", r: ConvertibleToInt = 0, i: ConvertibleToInt = 0, x=()) -> None: ...

class RowColField(Serialisable):
    tagname: ClassVar[str]
    x: Integer[Literal[False]]
    def __init__(self, x: ConvertibleToInt) -> None: ...

class AutoSortScope(Serialisable):
    pivotArea: Typed[PivotArea, Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, pivotArea: PivotArea) -> None: ...

class FieldItem(Serialisable):
    tagname: ClassVar[str]
    n: String[Literal[True]]
    t: Set[_ItemType]
    h: Bool[Literal[True]]
    s: Bool[Literal[True]]
    sd: Bool[Literal[True]]
    f: Bool[Literal[True]]
    m: Bool[Literal[True]]
    c: Bool[Literal[True]]
    x: Integer[Literal[True]]
    d: Bool[Literal[True]]
    e: Bool[Literal[True]]
    def __init__(
        self,
        n: str | None = None,
        t: _ItemType = "data",
        h: _ConvertibleToBool | None = None,
        s: _ConvertibleToBool | None = None,
        sd: _ConvertibleToBool | None = True,
        f: _ConvertibleToBool | None = None,
        m: _ConvertibleToBool | None = None,
        c: _ConvertibleToBool | None = None,
        x: ConvertibleToInt | None = None,
        d: _ConvertibleToBool | None = None,
        e: _ConvertibleToBool | None = None,
    ) -> None: ...

class PivotField(Serialisable):
    tagname: ClassVar[str]
    items: Incomplete
    autoSortScope: Typed[AutoSortScope, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    name: String[Literal[True]]
    axis: NoneSet[_PivotAxis]
    dataField: Bool[Literal[True]]
    subtotalCaption: String[Literal[True]]
    showDropDowns: Bool[Literal[True]]
    hiddenLevel: Bool[Literal[True]]
    uniqueMemberProperty: String[Literal[True]]
    compact: Bool[Literal[True]]
    allDrilled: Bool[Literal[True]]
    numFmtId: Integer[Literal[True]]
    outline: Bool[Literal[True]]
    subtotalTop: Bool[Literal[True]]
    dragToRow: Bool[Literal[True]]
    dragToCol: Bool[Literal[True]]
    multipleItemSelectionAllowed: Bool[Literal[True]]
    dragToPage: Bool[Literal[True]]
    dragToData: Bool[Literal[True]]
    dragOff: Bool[Literal[True]]
    showAll: Bool[Literal[True]]
    insertBlankRow: Bool[Literal[True]]
    serverField: Bool[Literal[True]]
    insertPageBreak: Bool[Literal[True]]
    autoShow: Bool[Literal[True]]
    topAutoShow: Bool[Literal[True]]
    hideNewItems: Bool[Literal[True]]
    measureFilter: Bool[Literal[True]]
    includeNewItemsInFilter: Bool[Literal[True]]
    itemPageCount: Integer[Literal[True]]
    sortType: Set[_PivotFieldSortType]
    dataSourceSort: Bool[Literal[True]]
    nonAutoSortDefault: Bool[Literal[True]]
    rankBy: Integer[Literal[True]]
    defaultSubtotal: Bool[Literal[True]]
    sumSubtotal: Bool[Literal[True]]
    countASubtotal: Bool[Literal[True]]
    avgSubtotal: Bool[Literal[True]]
    maxSubtotal: Bool[Literal[True]]
    minSubtotal: Bool[Literal[True]]
    productSubtotal: Bool[Literal[True]]
    countSubtotal: Bool[Literal[True]]
    stdDevSubtotal: Bool[Literal[True]]
    stdDevPSubtotal: Bool[Literal[True]]
    varSubtotal: Bool[Literal[True]]
    varPSubtotal: Bool[Literal[True]]
    showPropCell: Bool[Literal[True]]
    showPropTip: Bool[Literal[True]]
    showPropAsCaption: Bool[Literal[True]]
    defaultAttributeDrillState: Bool[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        items=(),
        autoSortScope: AutoSortScope | None = None,
        name: str | None = None,
        axis: _PivotAxis | Literal["none"] | None = None,
        dataField: _ConvertibleToBool | None = None,
        subtotalCaption: str | None = None,
        showDropDowns: _ConvertibleToBool | None = True,
        hiddenLevel: _ConvertibleToBool | None = None,
        uniqueMemberProperty: str | None = None,
        compact: _ConvertibleToBool | None = True,
        allDrilled: _ConvertibleToBool | None = None,
        numFmtId: ConvertibleToInt | None = None,
        outline: _ConvertibleToBool | None = True,
        subtotalTop: _ConvertibleToBool | None = True,
        dragToRow: _ConvertibleToBool | None = True,
        dragToCol: _ConvertibleToBool | None = True,
        multipleItemSelectionAllowed: _ConvertibleToBool | None = None,
        dragToPage: _ConvertibleToBool | None = True,
        dragToData: _ConvertibleToBool | None = True,
        dragOff: _ConvertibleToBool | None = True,
        showAll: _ConvertibleToBool | None = True,
        insertBlankRow: _ConvertibleToBool | None = None,
        serverField: _ConvertibleToBool | None = None,
        insertPageBreak: _ConvertibleToBool | None = None,
        autoShow: _ConvertibleToBool | None = None,
        topAutoShow: _ConvertibleToBool | None = True,
        hideNewItems: _ConvertibleToBool | None = None,
        measureFilter: _ConvertibleToBool | None = None,
        includeNewItemsInFilter: _ConvertibleToBool | None = None,
        itemPageCount: ConvertibleToInt | None = 10,
        sortType: _PivotFieldSortType = "manual",
        dataSourceSort: _ConvertibleToBool | None = None,
        nonAutoSortDefault: _ConvertibleToBool | None = None,
        rankBy: ConvertibleToInt | None = None,
        defaultSubtotal: _ConvertibleToBool | None = True,
        sumSubtotal: _ConvertibleToBool | None = None,
        countASubtotal: _ConvertibleToBool | None = None,
        avgSubtotal: _ConvertibleToBool | None = None,
        maxSubtotal: _ConvertibleToBool | None = None,
        minSubtotal: _ConvertibleToBool | None = None,
        productSubtotal: _ConvertibleToBool | None = None,
        countSubtotal: _ConvertibleToBool | None = None,
        stdDevSubtotal: _ConvertibleToBool | None = None,
        stdDevPSubtotal: _ConvertibleToBool | None = None,
        varSubtotal: _ConvertibleToBool | None = None,
        varPSubtotal: _ConvertibleToBool | None = None,
        showPropCell: _ConvertibleToBool | None = None,
        showPropTip: _ConvertibleToBool | None = None,
        showPropAsCaption: _ConvertibleToBool | None = None,
        defaultAttributeDrillState: _ConvertibleToBool | None = None,
        extLst: Unused = None,
    ) -> None: ...

class Location(Serialisable):
    tagname: ClassVar[str]
    ref: String[Literal[False]]
    firstHeaderRow: Integer[Literal[False]]
    firstDataRow: Integer[Literal[False]]
    firstDataCol: Integer[Literal[False]]
    rowPageCount: Integer[Literal[True]]
    colPageCount: Integer[Literal[True]]
    def __init__(
        self,
        ref: str,
        firstHeaderRow: ConvertibleToInt,
        firstDataRow: ConvertibleToInt,
        firstDataCol: ConvertibleToInt,
        rowPageCount: ConvertibleToInt | None = None,
        colPageCount: ConvertibleToInt | None = None,
    ) -> None: ...

class TableDefinition(Serialisable):
    mime_type: str
    rel_type: str
    tagname: ClassVar[str]
    cache: Incomplete
    name: String[Literal[False]]
    cacheId: Integer[Literal[False]]
    dataOnRows: Bool[Literal[False]]
    dataPosition: Integer[Literal[True]]
    dataCaption: String[Literal[False]]
    grandTotalCaption: String[Literal[True]]
    errorCaption: String[Literal[True]]
    showError: Bool[Literal[False]]
    missingCaption: String[Literal[True]]
    showMissing: Bool[Literal[False]]
    pageStyle: String[Literal[True]]
    pivotTableStyle: String[Literal[True]]
    vacatedStyle: String[Literal[True]]
    tag: String[Literal[True]]
    updatedVersion: Integer[Literal[False]]
    minRefreshableVersion: Integer[Literal[False]]
    asteriskTotals: Bool[Literal[False]]
    showItems: Bool[Literal[False]]
    editData: Bool[Literal[False]]
    disableFieldList: Bool[Literal[False]]
    showCalcMbrs: Bool[Literal[False]]
    visualTotals: Bool[Literal[False]]
    showMultipleLabel: Bool[Literal[False]]
    showDataDropDown: Bool[Literal[False]]
    showDrill: Bool[Literal[False]]
    printDrill: Bool[Literal[False]]
    showMemberPropertyTips: Bool[Literal[False]]
    showDataTips: Bool[Literal[False]]
    enableWizard: Bool[Literal[False]]
    enableDrill: Bool[Literal[False]]
    enableFieldProperties: Bool[Literal[False]]
    preserveFormatting: Bool[Literal[False]]
    useAutoFormatting: Bool[Literal[False]]
    pageWrap: Integer[Literal[False]]
    pageOverThenDown: Bool[Literal[False]]
    subtotalHiddenItems: Bool[Literal[False]]
    rowGrandTotals: Bool[Literal[False]]
    colGrandTotals: Bool[Literal[False]]
    fieldPrintTitles: Bool[Literal[False]]
    itemPrintTitles: Bool[Literal[False]]
    mergeItem: Bool[Literal[False]]
    showDropZones: Bool[Literal[False]]
    createdVersion: Integer[Literal[False]]
    indent: Integer[Literal[False]]
    showEmptyRow: Bool[Literal[False]]
    showEmptyCol: Bool[Literal[False]]
    showHeaders: Bool[Literal[False]]
    compact: Bool[Literal[False]]
    outline: Bool[Literal[False]]
    outlineData: Bool[Literal[False]]
    compactData: Bool[Literal[False]]
    published: Bool[Literal[False]]
    gridDropZones: Bool[Literal[False]]
    immersive: Bool[Literal[False]]
    multipleFieldFilters: Bool[Literal[False]]
    chartFormat: Integer[Literal[False]]
    rowHeaderCaption: String[Literal[True]]
    colHeaderCaption: String[Literal[True]]
    fieldListSortAscending: Bool[Literal[False]]
    mdxSubqueries: Bool[Literal[False]]
    customListSort: Bool[Literal[True]]
    autoFormatId: Integer[Literal[True]]
    applyNumberFormats: Bool[Literal[False]]
    applyBorderFormats: Bool[Literal[False]]
    applyFontFormats: Bool[Literal[False]]
    applyPatternFormats: Bool[Literal[False]]
    applyAlignmentFormats: Bool[Literal[False]]
    applyWidthHeightFormats: Bool[Literal[False]]
    location: Typed[Location, Literal[False]]
    pivotFields: Incomplete
    rowFields: Incomplete
    rowItems: Incomplete
    colFields: Incomplete
    colItems: Incomplete
    pageFields: Incomplete
    dataFields: Incomplete
    formats: Incomplete
    conditionalFormats: Typed[ConditionalFormatList, Literal[True]]
    chartFormats: Incomplete
    pivotHierarchies: Incomplete
    pivotTableStyleInfo: Typed[PivotTableStyle, Literal[True]]
    filters: Incomplete
    rowHierarchiesUsage: Typed[RowHierarchiesUsage, Literal[True]]
    colHierarchiesUsage: Typed[ColHierarchiesUsage, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    id: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        name: str,
        cacheId: ConvertibleToInt,
        dataOnRows: _ConvertibleToBool = False,
        dataPosition: ConvertibleToInt | None = None,
        *,
        dataCaption: str,
        grandTotalCaption: str | None = None,
        errorCaption: str | None = None,
        showError: _ConvertibleToBool = False,
        missingCaption: str | None = None,
        showMissing: _ConvertibleToBool = True,
        pageStyle: str | None = None,
        pivotTableStyle: str | None = None,
        vacatedStyle: str | None = None,
        tag: str | None = None,
        updatedVersion: ConvertibleToInt = 0,
        minRefreshableVersion: ConvertibleToInt = 0,
        asteriskTotals: _ConvertibleToBool = False,
        showItems: _ConvertibleToBool = True,
        editData: _ConvertibleToBool = False,
        disableFieldList: _ConvertibleToBool = False,
        showCalcMbrs: _ConvertibleToBool = True,
        visualTotals: _ConvertibleToBool = True,
        showMultipleLabel: _ConvertibleToBool = True,
        showDataDropDown: _ConvertibleToBool = True,
        showDrill: _ConvertibleToBool = True,
        printDrill: _ConvertibleToBool = False,
        showMemberPropertyTips: _ConvertibleToBool = True,
        showDataTips: _ConvertibleToBool = True,
        enableWizard: _ConvertibleToBool = True,
        enableDrill: _ConvertibleToBool = True,
        enableFieldProperties: _ConvertibleToBool = True,
        preserveFormatting: _ConvertibleToBool = True,
        useAutoFormatting: _ConvertibleToBool = False,
        pageWrap: ConvertibleToInt = 0,
        pageOverThenDown: _ConvertibleToBool = False,
        subtotalHiddenItems: _ConvertibleToBool = False,
        rowGrandTotals: _ConvertibleToBool = True,
        colGrandTotals: _ConvertibleToBool = True,
        fieldPrintTitles: _ConvertibleToBool = False,
        itemPrintTitles: _ConvertibleToBool = False,
        mergeItem: _ConvertibleToBool = False,
        showDropZones: _ConvertibleToBool = True,
        createdVersion: ConvertibleToInt = 0,
        indent: ConvertibleToInt = 1,
        showEmptyRow: _ConvertibleToBool = False,
        showEmptyCol: _ConvertibleToBool = False,
        showHeaders: _ConvertibleToBool = True,
        compact: _ConvertibleToBool = True,
        outline: _ConvertibleToBool = False,
        outlineData: _ConvertibleToBool = False,
        compactData: _ConvertibleToBool = True,
        published: _ConvertibleToBool = False,
        gridDropZones: _ConvertibleToBool = False,
        immersive: _ConvertibleToBool = True,
        multipleFieldFilters: _ConvertibleToBool = None,
        chartFormat: ConvertibleToInt = 0,
        rowHeaderCaption: str | None = None,
        colHeaderCaption: str | None = None,
        fieldListSortAscending: _ConvertibleToBool = None,
        mdxSubqueries: _ConvertibleToBool = None,
        customListSort: _ConvertibleToBool | None = None,
        autoFormatId: ConvertibleToInt | None = None,
        applyNumberFormats: _ConvertibleToBool = False,
        applyBorderFormats: _ConvertibleToBool = False,
        applyFontFormats: _ConvertibleToBool = False,
        applyPatternFormats: _ConvertibleToBool = False,
        applyAlignmentFormats: _ConvertibleToBool = False,
        applyWidthHeightFormats: _ConvertibleToBool = False,
        location: Location,
        pivotFields=(),
        rowFields=(),
        rowItems=(),
        colFields=(),
        colItems=(),
        pageFields=(),
        dataFields=(),
        formats=(),
        conditionalFormats: ConditionalFormatList | None = None,
        chartFormats=(),
        pivotHierarchies=(),
        pivotTableStyleInfo: PivotTableStyle | None = None,
        filters=(),
        rowHierarchiesUsage: RowHierarchiesUsage | None = None,
        colHierarchiesUsage: ColHierarchiesUsage | None = None,
        extLst: ExtensionList | None = None,
        id=None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        name: str,
        cacheId: ConvertibleToInt,
        dataOnRows: _ConvertibleToBool,
        dataPosition: ConvertibleToInt | None,
        dataCaption: str,
        grandTotalCaption: str | None,
        errorCaption: str | None,
        showError: _ConvertibleToBool,
        missingCaption: str | None,
        showMissing: _ConvertibleToBool,
        pageStyle: str | None,
        pivotTableStyle: str | None,
        vacatedStyle: str | None,
        tag: str | None,
        updatedVersion: ConvertibleToInt,
        minRefreshableVersion: ConvertibleToInt,
        asteriskTotals: _ConvertibleToBool,
        showItems: _ConvertibleToBool,
        editData: _ConvertibleToBool,
        disableFieldList: _ConvertibleToBool,
        showCalcMbrs: _ConvertibleToBool,
        visualTotals: _ConvertibleToBool,
        showMultipleLabel: _ConvertibleToBool,
        showDataDropDown: _ConvertibleToBool,
        showDrill: _ConvertibleToBool,
        printDrill: _ConvertibleToBool,
        showMemberPropertyTips: _ConvertibleToBool,
        showDataTips: _ConvertibleToBool,
        enableWizard: _ConvertibleToBool,
        enableDrill: _ConvertibleToBool,
        enableFieldProperties: _ConvertibleToBool,
        preserveFormatting: _ConvertibleToBool,
        useAutoFormatting: _ConvertibleToBool,
        pageWrap: ConvertibleToInt,
        pageOverThenDown: _ConvertibleToBool,
        subtotalHiddenItems: _ConvertibleToBool,
        rowGrandTotals: _ConvertibleToBool,
        colGrandTotals: _ConvertibleToBool,
        fieldPrintTitles: _ConvertibleToBool,
        itemPrintTitles: _ConvertibleToBool,
        mergeItem: _ConvertibleToBool,
        showDropZones: _ConvertibleToBool,
        createdVersion: ConvertibleToInt,
        indent: ConvertibleToInt,
        showEmptyRow: _ConvertibleToBool,
        showEmptyCol: _ConvertibleToBool,
        showHeaders: _ConvertibleToBool,
        compact: _ConvertibleToBool,
        outline: _ConvertibleToBool,
        outlineData: _ConvertibleToBool,
        compactData: _ConvertibleToBool,
        published: _ConvertibleToBool,
        gridDropZones: _ConvertibleToBool,
        immersive: _ConvertibleToBool,
        multipleFieldFilters: _ConvertibleToBool,
        chartFormat: ConvertibleToInt,
        rowHeaderCaption: str | None,
        colHeaderCaption: str | None,
        fieldListSortAscending: _ConvertibleToBool,
        mdxSubqueries: _ConvertibleToBool,
        customListSort: _ConvertibleToBool | None,
        autoFormatId: ConvertibleToInt | None,
        applyNumberFormats: _ConvertibleToBool,
        applyBorderFormats: _ConvertibleToBool,
        applyFontFormats: _ConvertibleToBool,
        applyPatternFormats: _ConvertibleToBool,
        applyAlignmentFormats: _ConvertibleToBool,
        applyWidthHeightFormats: _ConvertibleToBool,
        location: Location,
        pivotFields=(),
        rowFields=(),
        rowItems=(),
        colFields=(),
        colItems=(),
        pageFields=(),
        dataFields=(),
        formats=(),
        conditionalFormats: ConditionalFormatList | None = None,
        chartFormats=(),
        pivotHierarchies=(),
        pivotTableStyleInfo: PivotTableStyle | None = None,
        filters=(),
        rowHierarchiesUsage: RowHierarchiesUsage | None = None,
        colHierarchiesUsage: ColHierarchiesUsage | None = None,
        extLst: ExtensionList | None = None,
        id=None,
    ) -> None: ...
    def to_tree(self) -> Element: ...  # type: ignore[override]
    @property
    def path(self) -> str: ...
    def formatted_fields(self) -> dict[Incomplete, list[Incomplete]]: ...
    @property
    def summary(self) -> str: ...
