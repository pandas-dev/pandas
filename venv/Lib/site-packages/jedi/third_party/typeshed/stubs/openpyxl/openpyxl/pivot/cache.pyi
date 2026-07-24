from _typeshed import ConvertibleToFloat, ConvertibleToInt, Incomplete, Unused
from datetime import datetime
from typing import ClassVar, Literal, overload
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Bool, DateTime, Float, Integer, NoneSet, Set, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedInteger
from openpyxl.descriptors.sequence import NestedSequence
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.pivot.fields import Error, Missing, Number, Text, TupleList
from openpyxl.pivot.table import PivotArea
from openpyxl.xml.functions import Element

from ..xml._functions_overloads import _HasTagAndGet

_RangePrGroupBy: TypeAlias = Literal["range", "seconds", "minutes", "hours", "days", "months", "quarters", "years"]
_CacheSourceType: TypeAlias = Literal["worksheet", "external", "consolidation", "scenario"]

class MeasureDimensionMap(Serialisable):
    tagname: ClassVar[str]
    measureGroup: Integer[Literal[True]]
    dimension: Integer[Literal[True]]
    def __init__(self, measureGroup: ConvertibleToInt | None = None, dimension: ConvertibleToInt | None = None) -> None: ...

class MeasureGroup(Serialisable):
    tagname: ClassVar[str]
    name: String[Literal[False]]
    caption: String[Literal[False]]
    def __init__(self, name: str, caption: str) -> None: ...

class PivotDimension(Serialisable):
    tagname: ClassVar[str]
    measure: Bool[Literal[False]]
    name: String[Literal[False]]
    uniqueName: String[Literal[False]]
    caption: String[Literal[False]]
    @overload
    def __init__(self, measure: _ConvertibleToBool = None, *, name: str, uniqueName: str, caption: str) -> None: ...
    @overload
    def __init__(self, measure: _ConvertibleToBool, name: str, uniqueName: str, caption: str) -> None: ...

class CalculatedMember(Serialisable):
    tagname: ClassVar[str]
    name: String[Literal[False]]
    mdx: String[Literal[False]]
    memberName: String[Literal[True]]
    hierarchy: String[Literal[True]]
    parent: String[Literal[True]]
    solveOrder: Integer[Literal[True]]
    set: Bool[Literal[False]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        name: str,
        mdx: str,
        memberName: str,
        hierarchy: str,
        parent: str,
        solveOrder: ConvertibleToInt,
        set: _ConvertibleToBool = None,
        extLst: Unused = None,
    ) -> None: ...

class CalculatedItem(Serialisable):
    tagname: ClassVar[str]
    field: Integer[Literal[True]]
    formula: String[Literal[False]]
    pivotArea: Typed[PivotArea, Literal[False]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(self, field: ConvertibleToInt | None = None, *, formula: str, pivotArea: PivotArea, extLst=None) -> None: ...
    @overload
    def __init__(self, field: ConvertibleToInt | None, formula: str, pivotArea: PivotArea, extLst=None) -> None: ...

class ServerFormat(Serialisable):
    tagname: ClassVar[str]
    culture: String[Literal[True]]
    format: String[Literal[True]]
    def __init__(self, culture: str | None = None, format: str | None = None) -> None: ...

class Query(Serialisable):
    tagname: ClassVar[str]
    mdx: String[Literal[False]]
    tpls: Typed[TupleList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, mdx: str, tpls: TupleList | None = None) -> None: ...

class OLAPSet(Serialisable):
    tagname: ClassVar[str]
    count: Integer[Literal[False]]
    maxRank: Integer[Literal[False]]
    setDefinition: String[Literal[False]]
    sortType: Incomplete
    queryFailed: Bool[Literal[False]]
    tpls: Typed[TupleList, Literal[True]]
    sortByTuple: Typed[TupleList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        count: ConvertibleToInt,
        maxRank: ConvertibleToInt,
        setDefinition: str,
        sortType=None,
        queryFailed: _ConvertibleToBool = None,
        tpls: TupleList | None = None,
        sortByTuple: TupleList | None = None,
    ) -> None: ...

class PCDSDTCEntries(Serialisable):
    tagname: ClassVar[str]
    count: Integer[Literal[True]]
    m: Typed[Missing, Literal[True]]
    n: Typed[Number, Literal[True]]
    e: Typed[Error, Literal[True]]
    s: Typed[Text, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, count: ConvertibleToInt, m: Missing, n: Number, e: Error, s: Text) -> None: ...

class TupleCache(Serialisable):
    tagname: ClassVar[str]
    entries: Typed[PCDSDTCEntries, Literal[True]]
    sets: NestedSequence[list[OLAPSet]]
    queryCache: NestedSequence[list[Query]]
    serverFormats: NestedSequence[list[ServerFormat]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        entries: PCDSDTCEntries | None = None,
        sets: list[OLAPSet] | tuple[OLAPSet, ...] = (),
        queryCache: list[Query] | tuple[Query, ...] = (),
        serverFormats: list[ServerFormat] | tuple[ServerFormat, ...] = (),
        extLst: ExtensionList | None = None,
    ) -> None: ...

class OLAPKPI(Serialisable):
    tagname: ClassVar[str]
    uniqueName: String[Literal[False]]
    caption: String[Literal[True]]
    displayFolder: String[Literal[True]]
    measureGroup: String[Literal[True]]
    parent: String[Literal[True]]
    value: String[Literal[False]]
    goal: String[Literal[True]]
    status: String[Literal[True]]
    trend: String[Literal[True]]
    weight: String[Literal[True]]
    time: String[Literal[True]]
    def __init__(
        self,
        uniqueName: str | None = None,
        caption: str | None = None,
        displayFolder: str | None = None,
        measureGroup: str | None = None,
        parent: str | None = None,
        value: str | None = None,
        goal: str | None = None,
        status: str | None = None,
        trend: str | None = None,
        weight: str | None = None,
        time: str | None = None,
    ) -> None: ...

class GroupMember(Serialisable):
    tagname: ClassVar[str]
    uniqueName: String[Literal[False]]
    group: Bool[Literal[False]]
    def __init__(self, uniqueName: str, group: _ConvertibleToBool = None) -> None: ...

class LevelGroup(Serialisable):
    tagname: ClassVar[str]
    name: String[Literal[False]]
    uniqueName: String[Literal[False]]
    caption: String[Literal[False]]
    uniqueParent: String[Literal[False]]
    id: Integer[Literal[False]]
    groupMembers: NestedSequence[list[GroupMember]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        name: str,
        uniqueName: str,
        caption: str,
        uniqueParent: str,
        id: ConvertibleToInt,
        groupMembers: list[GroupMember] | tuple[GroupMember, ...] = (),
    ) -> None: ...

class GroupLevel(Serialisable):
    tagname: ClassVar[str]
    uniqueName: String[Literal[False]]
    caption: String[Literal[False]]
    user: Bool[Literal[False]]
    customRollUp: Bool[Literal[False]]
    groups: NestedSequence[list[LevelGroup]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        uniqueName: str,
        caption: str,
        user: _ConvertibleToBool = None,
        customRollUp: _ConvertibleToBool = None,
        groups: list[LevelGroup] | tuple[LevelGroup, ...] = (),
        extLst: ExtensionList | None = None,
    ) -> None: ...

class FieldUsage(Serialisable):
    tagname: ClassVar[str]
    x: Integer[Literal[False]]
    def __init__(self, x: ConvertibleToInt) -> None: ...

class CacheHierarchy(Serialisable):
    tagname: ClassVar[str]
    uniqueName: String[Literal[False]]
    caption: String[Literal[True]]
    measure: Bool[Literal[False]]
    set: Bool[Literal[False]]
    parentSet: Integer[Literal[True]]
    iconSet: Integer[Literal[False]]
    attribute: Bool[Literal[False]]
    time: Bool[Literal[False]]
    keyAttribute: Bool[Literal[False]]
    defaultMemberUniqueName: String[Literal[True]]
    allUniqueName: String[Literal[True]]
    allCaption: String[Literal[True]]
    dimensionUniqueName: String[Literal[True]]
    displayFolder: String[Literal[True]]
    measureGroup: String[Literal[True]]
    measures: Bool[Literal[False]]
    count: Integer[Literal[False]]
    oneField: Bool[Literal[False]]
    memberValueDatatype: Integer[Literal[True]]
    unbalanced: Bool[Literal[True]]
    unbalancedGroup: Bool[Literal[True]]
    hidden: Bool[Literal[False]]
    fieldsUsage: NestedSequence[list[FieldUsage]]
    groupLevels: NestedSequence[list[GroupLevel]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        uniqueName: str = "",
        caption: str | None = None,
        measure: _ConvertibleToBool = None,
        set: _ConvertibleToBool = None,
        parentSet: ConvertibleToInt | None = None,
        iconSet: ConvertibleToInt = 0,
        attribute: _ConvertibleToBool = None,
        time: _ConvertibleToBool = None,
        keyAttribute: _ConvertibleToBool = None,
        defaultMemberUniqueName: str | None = None,
        allUniqueName: str | None = None,
        allCaption: str | None = None,
        dimensionUniqueName: str | None = None,
        displayFolder: str | None = None,
        measureGroup: str | None = None,
        measures: _ConvertibleToBool = None,
        *,
        count: ConvertibleToInt,
        oneField: _ConvertibleToBool = None,
        memberValueDatatype: ConvertibleToInt | None = None,
        unbalanced: _ConvertibleToBool | None = None,
        unbalancedGroup: _ConvertibleToBool | None = None,
        hidden: _ConvertibleToBool = None,
        fieldsUsage: list[FieldUsage] | tuple[FieldUsage, ...] = (),
        groupLevels: list[FieldUsage] | tuple[FieldUsage, ...] = (),
        extLst: ExtensionList | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        uniqueName: str,
        caption: str | None,
        measure: _ConvertibleToBool,
        set: _ConvertibleToBool,
        parentSet: ConvertibleToInt | None,
        iconSet: ConvertibleToInt,
        attribute: _ConvertibleToBool,
        time: _ConvertibleToBool,
        keyAttribute: _ConvertibleToBool,
        defaultMemberUniqueName: str | None,
        allUniqueName: str | None,
        allCaption: str | None,
        dimensionUniqueName: str | None,
        displayFolder: str | None,
        measureGroup: str | None,
        measures: _ConvertibleToBool,
        count: ConvertibleToInt,
        oneField: _ConvertibleToBool = None,
        memberValueDatatype: ConvertibleToInt | None = None,
        unbalanced: _ConvertibleToBool | None = None,
        unbalancedGroup: _ConvertibleToBool | None = None,
        hidden: _ConvertibleToBool = None,
        fieldsUsage: list[FieldUsage] | tuple[FieldUsage, ...] = (),
        groupLevels: list[FieldUsage] | tuple[FieldUsage, ...] = (),
        extLst: ExtensionList | None = None,
    ) -> None: ...

class GroupItems(Serialisable):
    tagname: ClassVar[str]
    m: Incomplete
    n: Incomplete
    b: Incomplete
    e: Incomplete
    s: Incomplete
    d: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    __attrs__: ClassVar[tuple[str, ...]]
    def __init__(self, count=None, m=(), n=(), b=(), e=(), s=(), d=()) -> None: ...
    @property
    def count(self) -> int: ...

class RangePr(Serialisable):
    tagname: ClassVar[str]
    autoStart: Bool[Literal[True]]
    autoEnd: Bool[Literal[True]]
    groupBy: NoneSet[_RangePrGroupBy]
    startNum: Float[Literal[True]]
    endNum: Float[Literal[True]]
    startDate: DateTime[Literal[True]]
    endDate: DateTime[Literal[True]]
    groupInterval: Float[Literal[True]]
    def __init__(
        self,
        autoStart: _ConvertibleToBool | None = True,
        autoEnd: _ConvertibleToBool | None = True,
        groupBy: _RangePrGroupBy = "range",
        startNum: ConvertibleToFloat | None = None,
        endNum: ConvertibleToFloat | None = None,
        startDate: datetime | str | None = None,
        endDate: datetime | str | None = None,
        groupInterval: ConvertibleToFloat | None = 1,
    ) -> None: ...

class FieldGroup(Serialisable):
    tagname: ClassVar[str]
    par: Integer[Literal[True]]
    base: Integer[Literal[True]]
    rangePr: Typed[RangePr, Literal[True]]
    discretePr: NestedSequence[list[NestedInteger[Literal[False]]]]
    groupItems: Typed[GroupItems, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        par: ConvertibleToInt | None = None,
        base: ConvertibleToInt | None = None,
        rangePr: RangePr | None = None,
        discretePr: list[NestedInteger[Literal[False]]] | tuple[NestedInteger[Literal[False]], ...] = (),
        groupItems: GroupItems | None = None,
    ) -> None: ...

class SharedItems(Serialisable):
    tagname: ClassVar[str]
    m: Incomplete
    n: Incomplete
    b: Incomplete
    e: Incomplete
    s: Incomplete
    d: Incomplete
    containsSemiMixedTypes: Bool[Literal[True]]
    containsNonDate: Bool[Literal[True]]
    containsDate: Bool[Literal[True]]
    containsString: Bool[Literal[True]]
    containsBlank: Bool[Literal[True]]
    containsMixedTypes: Bool[Literal[True]]
    containsNumber: Bool[Literal[True]]
    containsInteger: Bool[Literal[True]]
    minValue: Float[Literal[True]]
    maxValue: Float[Literal[True]]
    minDate: DateTime[Literal[True]]
    maxDate: DateTime[Literal[True]]
    longText: Bool[Literal[True]]
    __attrs__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        _fields=(),
        containsSemiMixedTypes: _ConvertibleToBool | None = None,
        containsNonDate: _ConvertibleToBool | None = None,
        containsDate: _ConvertibleToBool | None = None,
        containsString: _ConvertibleToBool | None = None,
        containsBlank: _ConvertibleToBool | None = None,
        containsMixedTypes: _ConvertibleToBool | None = None,
        containsNumber: _ConvertibleToBool | None = None,
        containsInteger: _ConvertibleToBool | None = None,
        minValue: ConvertibleToFloat | None = None,
        maxValue: ConvertibleToFloat | None = None,
        minDate: datetime | str | None = None,
        maxDate: datetime | str | None = None,
        count: Unused = None,
        longText: _ConvertibleToBool | None = None,
    ) -> None: ...
    @property
    def count(self) -> int: ...

class CacheField(Serialisable):
    tagname: ClassVar[str]
    sharedItems: Typed[SharedItems, Literal[True]]
    fieldGroup: Typed[FieldGroup, Literal[True]]
    mpMap: NestedInteger[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    name: String[Literal[False]]
    caption: String[Literal[True]]
    propertyName: String[Literal[True]]
    serverField: Bool[Literal[True]]
    uniqueList: Bool[Literal[True]]
    numFmtId: Integer[Literal[True]]
    formula: String[Literal[True]]
    sqlType: Integer[Literal[True]]
    hierarchy: Integer[Literal[True]]
    level: Integer[Literal[True]]
    databaseField: Bool[Literal[True]]
    mappingCount: Integer[Literal[True]]
    memberPropertyField: Bool[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        sharedItems: SharedItems | None = None,
        fieldGroup: FieldGroup | None = None,
        mpMap: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        extLst: ExtensionList | None = None,
        *,
        name: str,
        caption: str | None = None,
        propertyName: str | None = None,
        serverField: _ConvertibleToBool | None = None,
        uniqueList: _ConvertibleToBool | None = True,
        numFmtId: ConvertibleToInt | None = None,
        formula: str | None = None,
        sqlType: ConvertibleToInt | None = 0,
        hierarchy: ConvertibleToInt | None = 0,
        level: ConvertibleToInt | None = 0,
        databaseField: _ConvertibleToBool | None = True,
        mappingCount: ConvertibleToInt | None = None,
        memberPropertyField: _ConvertibleToBool | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        sharedItems: SharedItems | None,
        fieldGroup: FieldGroup | None,
        mpMap: Incomplete | None,
        extLst: ExtensionList | None,
        name: str,
        caption: str | None = None,
        propertyName: str | None = None,
        serverField: _ConvertibleToBool | None = None,
        uniqueList: _ConvertibleToBool | None = True,
        numFmtId: ConvertibleToInt | None = None,
        formula: str | None = None,
        sqlType: ConvertibleToInt | None = 0,
        hierarchy: ConvertibleToInt | None = 0,
        level: ConvertibleToInt | None = 0,
        databaseField: _ConvertibleToBool | None = True,
        mappingCount: ConvertibleToInt | None = None,
        memberPropertyField: _ConvertibleToBool | None = None,
    ) -> None: ...

class RangeSet(Serialisable):
    tagname: ClassVar[str]
    i1: Integer[Literal[True]]
    i2: Integer[Literal[True]]
    i3: Integer[Literal[True]]
    i4: Integer[Literal[True]]
    ref: String[Literal[False]]
    name: String[Literal[True]]
    sheet: String[Literal[True]]
    @overload
    def __init__(
        self,
        i1: ConvertibleToInt | None = None,
        i2: ConvertibleToInt | None = None,
        i3: ConvertibleToInt | None = None,
        i4: ConvertibleToInt | None = None,
        *,
        ref: str,
        name: str | None = None,
        sheet: str | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        i1: ConvertibleToInt | None,
        i2: ConvertibleToInt | None,
        i3: ConvertibleToInt | None,
        i4: ConvertibleToInt | None,
        ref: str,
        name: str | None = None,
        sheet: str | None = None,
    ) -> None: ...

class PageItem(Serialisable):
    tagname: ClassVar[str]
    name: String[Literal[False]]
    def __init__(self, name: str) -> None: ...

class Consolidation(Serialisable):
    tagname: ClassVar[str]
    autoPage: Bool[Literal[True]]
    pages: NestedSequence[list[PageItem]]
    rangeSets: NestedSequence[list[RangeSet]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        autoPage: _ConvertibleToBool | None = None,
        pages: list[PageItem] | tuple[PageItem, ...] = (),
        rangeSets: list[RangeSet] | tuple[RangeSet, ...] = (),
    ) -> None: ...

class WorksheetSource(Serialisable):
    tagname: ClassVar[str]
    ref: String[Literal[True]]
    name: String[Literal[True]]
    sheet: String[Literal[True]]
    def __init__(self, ref: str | None = None, name: str | None = None, sheet: str | None = None) -> None: ...

class CacheSource(Serialisable):
    tagname: ClassVar[str]
    type: Set[_CacheSourceType]
    connectionId: Integer[Literal[True]]
    worksheetSource: Typed[WorksheetSource, Literal[True]]
    consolidation: Typed[Consolidation, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        type: _CacheSourceType,
        connectionId: ConvertibleToInt | None = None,
        worksheetSource: WorksheetSource | None = None,
        consolidation: Consolidation | None = None,
        extLst: ExtensionList | None = None,
    ) -> None: ...

class CacheDefinition(Serialisable):
    mime_type: str
    rel_type: str
    records: Incomplete
    tagname: ClassVar[str]
    invalid: Bool[Literal[True]]
    saveData: Bool[Literal[True]]
    refreshOnLoad: Bool[Literal[True]]
    optimizeMemory: Bool[Literal[True]]
    enableRefresh: Bool[Literal[True]]
    refreshedBy: String[Literal[True]]
    refreshedDate: Float[Literal[True]]
    refreshedDateIso: DateTime[Literal[True]]
    backgroundQuery: Bool[Literal[True]]
    missingItemsLimit: Integer[Literal[True]]
    createdVersion: Integer[Literal[True]]
    refreshedVersion: Integer[Literal[True]]
    minRefreshableVersion: Integer[Literal[True]]
    recordCount: Integer[Literal[True]]
    upgradeOnRefresh: Bool[Literal[True]]
    supportSubquery: Bool[Literal[True]]
    supportAdvancedDrill: Bool[Literal[True]]
    cacheSource: Typed[CacheSource, Literal[True]]
    cacheFields: Incomplete
    cacheHierarchies: Incomplete
    kpis: NestedSequence[list[OLAPKPI]]
    tupleCache: Typed[TupleCache, Literal[True]]
    calculatedItems: Incomplete
    calculatedMembers: Incomplete
    dimensions: Incomplete
    measureGroups: Incomplete
    maps: Incomplete
    extLst: Typed[ExtensionList, Literal[True]]
    id: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        invalid: _ConvertibleToBool | None = None,
        saveData: _ConvertibleToBool | None = None,
        refreshOnLoad: _ConvertibleToBool | None = None,
        optimizeMemory: _ConvertibleToBool | None = None,
        enableRefresh: _ConvertibleToBool | None = None,
        refreshedBy: str | None = None,
        refreshedDate: ConvertibleToFloat | None = None,
        refreshedDateIso: datetime | str | None = None,
        backgroundQuery: _ConvertibleToBool | None = None,
        missingItemsLimit: ConvertibleToInt | None = None,
        createdVersion: ConvertibleToInt | None = None,
        refreshedVersion: ConvertibleToInt | None = None,
        minRefreshableVersion: ConvertibleToInt | None = None,
        recordCount: ConvertibleToInt | None = None,
        upgradeOnRefresh: _ConvertibleToBool | None = None,
        tupleCache: TupleCache | None = None,
        supportSubquery: _ConvertibleToBool | None = None,
        supportAdvancedDrill: _ConvertibleToBool | None = None,
        *,
        cacheSource: CacheSource,
        cacheFields=(),
        cacheHierarchies=(),
        kpis: list[OLAPKPI] | tuple[OLAPKPI, ...] = (),
        calculatedItems=(),
        calculatedMembers=(),
        dimensions=(),
        measureGroups=(),
        maps=(),
        extLst: ExtensionList | None = None,
        id=None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        invalid: _ConvertibleToBool | None,
        saveData: _ConvertibleToBool | None,
        refreshOnLoad: _ConvertibleToBool | None,
        optimizeMemory: _ConvertibleToBool | None,
        enableRefresh: _ConvertibleToBool | None,
        refreshedBy: str | None,
        refreshedDate: ConvertibleToFloat | None,
        refreshedDateIso: datetime | str | None,
        backgroundQuery: _ConvertibleToBool | None,
        missingItemsLimit: ConvertibleToInt | None,
        createdVersion: ConvertibleToInt | None,
        refreshedVersion: ConvertibleToInt | None,
        minRefreshableVersion: ConvertibleToInt | None,
        recordCount: ConvertibleToInt | None,
        upgradeOnRefresh: _ConvertibleToBool | None,
        tupleCache: TupleCache | None,
        supportSubquery: _ConvertibleToBool | None,
        supportAdvancedDrill: _ConvertibleToBool | None,
        cacheSource: CacheSource,
        cacheFields=(),
        cacheHierarchies=(),
        kpis=(),
        calculatedItems=(),
        calculatedMembers=(),
        dimensions=(),
        measureGroups=(),
        maps=(),
        extLst: ExtensionList | None = None,
        id=None,
    ) -> None: ...
    def to_tree(self) -> Element: ...  # type: ignore[override]
    @property
    def path(self) -> str: ...
