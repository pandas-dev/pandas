from typing import ClassVar, Literal

from openpyxl.descriptors.base import Bool, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.styles.colors import Color, ColorDescriptor

class Outline(Serialisable):
    tagname: ClassVar[str]
    applyStyles: Bool[Literal[True]]
    summaryBelow: Bool[Literal[True]]
    summaryRight: Bool[Literal[True]]
    showOutlineSymbols: Bool[Literal[True]]
    def __init__(
        self,
        applyStyles: _ConvertibleToBool | None = None,
        summaryBelow: _ConvertibleToBool | None = None,
        summaryRight: _ConvertibleToBool | None = None,
        showOutlineSymbols: _ConvertibleToBool | None = None,
    ) -> None: ...

class PageSetupProperties(Serialisable):
    tagname: ClassVar[str]
    autoPageBreaks: Bool[Literal[True]]
    fitToPage: Bool[Literal[True]]
    def __init__(self, autoPageBreaks: _ConvertibleToBool | None = None, fitToPage: _ConvertibleToBool | None = None) -> None: ...

class WorksheetProperties(Serialisable):
    tagname: ClassVar[str]
    codeName: String[Literal[True]]
    enableFormatConditionsCalculation: Bool[Literal[True]]
    filterMode: Bool[Literal[True]]
    published: Bool[Literal[True]]
    syncHorizontal: Bool[Literal[True]]
    syncRef: String[Literal[True]]
    syncVertical: Bool[Literal[True]]
    transitionEvaluation: Bool[Literal[True]]
    transitionEntry: Bool[Literal[True]]
    tabColor: ColorDescriptor[Literal[True]]
    outlinePr: Typed[Outline, Literal[True]]
    pageSetUpPr: Typed[PageSetupProperties, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        codeName: str | None = None,
        enableFormatConditionsCalculation: _ConvertibleToBool | None = None,
        filterMode: _ConvertibleToBool | None = None,
        published: _ConvertibleToBool | None = None,
        syncHorizontal: _ConvertibleToBool | None = None,
        syncRef: str | None = None,
        syncVertical: _ConvertibleToBool | None = None,
        transitionEvaluation: _ConvertibleToBool | None = None,
        transitionEntry: _ConvertibleToBool | None = None,
        tabColor: str | Color | None = None,
        outlinePr: Outline | None = None,
        pageSetUpPr: PageSetupProperties | None = None,
    ) -> None: ...
