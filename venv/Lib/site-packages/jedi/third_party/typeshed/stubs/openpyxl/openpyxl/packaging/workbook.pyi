from _typeshed import ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl import _VisibilityType
from openpyxl.descriptors.base import Alias, Bool, Integer, NoneSet, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedString
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.workbook.defined_name import DefinedNameList
from openpyxl.workbook.function_group import FunctionGroupList
from openpyxl.workbook.properties import CalcProperties, FileVersion, WorkbookProperties
from openpyxl.workbook.protection import FileSharing, WorkbookProtection
from openpyxl.workbook.smart_tags import SmartTagList, SmartTagProperties
from openpyxl.workbook.web import WebPublishing, WebPublishObjectList
from openpyxl.xml.functions import Element

_WorkbookPackageConformance: TypeAlias = Literal["strict", "transitional"]

class FileRecoveryProperties(Serialisable):
    tagname: ClassVar[str]
    autoRecover: Bool[Literal[True]]
    crashSave: Bool[Literal[True]]
    dataExtractLoad: Bool[Literal[True]]
    repairLoad: Bool[Literal[True]]
    def __init__(
        self,
        autoRecover: _ConvertibleToBool | None = None,
        crashSave: _ConvertibleToBool | None = None,
        dataExtractLoad: _ConvertibleToBool | None = None,
        repairLoad: _ConvertibleToBool | None = None,
    ) -> None: ...

class ChildSheet(Serialisable):
    tagname: ClassVar[str]
    name: String[Literal[False]]
    sheetId: Integer[Literal[False]]
    state: NoneSet[_VisibilityType]
    id: Incomplete
    def __init__(
        self, name: str, sheetId: ConvertibleToInt, state: _VisibilityType | Literal["none"] | None = "visible", id=None
    ) -> None: ...

class PivotCache(Serialisable):
    tagname: ClassVar[str]
    cacheId: Integer[Literal[False]]
    id: Incomplete
    def __init__(self, cacheId: ConvertibleToInt, id=None) -> None: ...

class WorkbookPackage(Serialisable):
    tagname: ClassVar[str]
    conformance: NoneSet[_WorkbookPackageConformance]
    fileVersion: Typed[FileVersion, Literal[True]]
    fileSharing: Typed[FileSharing, Literal[True]]
    workbookPr: Typed[WorkbookProperties, Literal[True]]
    properties: Alias
    workbookProtection: Typed[WorkbookProtection, Literal[True]]
    bookViews: Incomplete
    sheets: Incomplete  # NestedSequence[ChildSheet]
    functionGroups: Typed[FunctionGroupList, Literal[True]]
    externalReferences: Incomplete
    definedNames: Typed[DefinedNameList, Literal[True]]
    calcPr: Typed[CalcProperties, Literal[True]]
    oleSize: NestedString[Literal[True]]
    customWorkbookViews: Incomplete
    pivotCaches: Incomplete
    smartTagPr: Typed[SmartTagProperties, Literal[True]]
    smartTagTypes: Typed[SmartTagList, Literal[True]]
    webPublishing: Typed[WebPublishing, Literal[True]]
    fileRecoveryPr: Typed[FileRecoveryProperties, Literal[True]]
    webPublishObjects: Typed[WebPublishObjectList, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    Ignorable: NestedString[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        conformance: _WorkbookPackageConformance | Literal["none"] | None = None,
        fileVersion: FileVersion | None = None,
        fileSharing: FileSharing | None = None,
        workbookPr: WorkbookProperties | None = None,
        workbookProtection: WorkbookProtection | None = None,
        bookViews=(),
        sheets=(),
        functionGroups: FunctionGroupList | None = None,
        externalReferences=(),
        definedNames: DefinedNameList | None = None,
        calcPr: CalcProperties | None = None,
        oleSize: object = None,
        customWorkbookViews=(),
        pivotCaches=(),
        smartTagPr: SmartTagProperties | None = None,
        smartTagTypes: SmartTagList | None = None,
        webPublishing: WebPublishing | None = None,
        fileRecoveryPr: FileRecoveryProperties | None = None,
        webPublishObjects: WebPublishObjectList | None = None,
        extLst: Unused = None,
        Ignorable: Unused = None,
    ) -> None: ...
    def to_tree(self) -> Element: ...  # type: ignore[override]
    @property
    def active(self) -> int: ...
