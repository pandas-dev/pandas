from _typeshed import ConvertibleToFloat, ConvertibleToInt, Incomplete
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Bool, Float, Integer, NoneSet, String, _ConvertibleToBool
from openpyxl.descriptors.serialisable import Serialisable

_WorkbookPropertiesShowObjects: TypeAlias = Literal["all", "placeholders"]
_WorkbookPropertiesUpdateLinks: TypeAlias = Literal["userSet", "never", "always"]
_CalcPropertiesCalcMode: TypeAlias = Literal["manual", "auto", "autoNoTable"]
_CalcPropertiesRefMode: TypeAlias = Literal["A1", "R1C1"]

class WorkbookProperties(Serialisable):
    tagname: ClassVar[str]
    date1904: Bool[Literal[True]]
    dateCompatibility: Bool[Literal[True]]
    showObjects: NoneSet[_WorkbookPropertiesShowObjects]
    showBorderUnselectedTables: Bool[Literal[True]]
    filterPrivacy: Bool[Literal[True]]
    promptedSolutions: Bool[Literal[True]]
    showInkAnnotation: Bool[Literal[True]]
    backupFile: Bool[Literal[True]]
    saveExternalLinkValues: Bool[Literal[True]]
    updateLinks: NoneSet[_WorkbookPropertiesUpdateLinks]
    codeName: String[Literal[True]]
    hidePivotFieldList: Bool[Literal[True]]
    showPivotChartFilter: Bool[Literal[True]]
    allowRefreshQuery: Bool[Literal[True]]
    publishItems: Bool[Literal[True]]
    checkCompatibility: Bool[Literal[True]]
    autoCompressPictures: Bool[Literal[True]]
    refreshAllConnections: Bool[Literal[True]]
    defaultThemeVersion: Integer[Literal[True]]
    def __init__(
        self,
        date1904: _ConvertibleToBool | None = None,
        dateCompatibility: _ConvertibleToBool | None = None,
        showObjects: _WorkbookPropertiesShowObjects | Literal["none"] | None = None,
        showBorderUnselectedTables: _ConvertibleToBool | None = None,
        filterPrivacy: _ConvertibleToBool | None = None,
        promptedSolutions: _ConvertibleToBool | None = None,
        showInkAnnotation: _ConvertibleToBool | None = None,
        backupFile: _ConvertibleToBool | None = None,
        saveExternalLinkValues: _ConvertibleToBool | None = None,
        updateLinks: _WorkbookPropertiesUpdateLinks | Literal["none"] | None = None,
        codeName: str | None = None,
        hidePivotFieldList: _ConvertibleToBool | None = None,
        showPivotChartFilter: _ConvertibleToBool | None = None,
        allowRefreshQuery: _ConvertibleToBool | None = None,
        publishItems: _ConvertibleToBool | None = None,
        checkCompatibility: _ConvertibleToBool | None = None,
        autoCompressPictures: _ConvertibleToBool | None = None,
        refreshAllConnections: _ConvertibleToBool | None = None,
        defaultThemeVersion: ConvertibleToInt | None = None,
    ) -> None: ...

class CalcProperties(Serialisable):
    tagname: ClassVar[str]
    calcId: Integer[Literal[False]]
    calcMode: NoneSet[_CalcPropertiesCalcMode]
    fullCalcOnLoad: Bool[Literal[True]]
    refMode: NoneSet[_CalcPropertiesRefMode]
    iterate: Bool[Literal[True]]
    iterateCount: Integer[Literal[True]]
    iterateDelta: Float[Literal[True]]
    fullPrecision: Bool[Literal[True]]
    calcCompleted: Bool[Literal[True]]
    calcOnSave: Bool[Literal[True]]
    concurrentCalc: Bool[Literal[True]]
    concurrentManualCount: Integer[Literal[True]]
    forceFullCalc: Bool[Literal[True]]
    def __init__(
        self,
        calcId: ConvertibleToInt = 124519,
        calcMode: _CalcPropertiesCalcMode | Literal["none"] | None = None,
        fullCalcOnLoad: _ConvertibleToBool | None = True,
        refMode: _CalcPropertiesRefMode | Literal["none"] | None = None,
        iterate: _ConvertibleToBool | None = None,
        iterateCount: ConvertibleToInt | None = None,
        iterateDelta: ConvertibleToFloat | None = None,
        fullPrecision: _ConvertibleToBool | None = None,
        calcCompleted: _ConvertibleToBool | None = None,
        calcOnSave: _ConvertibleToBool | None = None,
        concurrentCalc: _ConvertibleToBool | None = None,
        concurrentManualCount: ConvertibleToInt | None = None,
        forceFullCalc: _ConvertibleToBool | None = None,
    ) -> None: ...

class FileVersion(Serialisable):
    tagname: ClassVar[str]
    appName: String[Literal[True]]
    lastEdited: String[Literal[True]]
    lowestEdited: String[Literal[True]]
    rupBuild: String[Literal[True]]
    codeName: Incomplete
    def __init__(
        self,
        appName: str | None = None,
        lastEdited: str | None = None,
        lowestEdited: str | None = None,
        rupBuild: str | None = None,
        codeName=None,
    ) -> None: ...
