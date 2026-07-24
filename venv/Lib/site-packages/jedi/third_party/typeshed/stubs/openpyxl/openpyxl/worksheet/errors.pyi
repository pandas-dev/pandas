from _typeshed import Incomplete
from typing import ClassVar, Literal

from openpyxl.descriptors.base import Bool, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.serialisable import Serialisable

class Extension(Serialisable):
    tagname: ClassVar[str]
    uri: String[Literal[True]]
    def __init__(self, uri: str | None = None) -> None: ...

class ExtensionList(Serialisable):
    tagname: ClassVar[str]
    ext: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, ext=()) -> None: ...

class IgnoredError(Serialisable):
    tagname: ClassVar[str]
    sqref: Incomplete
    evalError: Bool[Literal[True]]
    twoDigitTextYear: Bool[Literal[True]]
    numberStoredAsText: Bool[Literal[True]]
    formula: Bool[Literal[True]]
    formulaRange: Bool[Literal[True]]
    unlockedFormula: Bool[Literal[True]]
    emptyCellReference: Bool[Literal[True]]
    listDataValidation: Bool[Literal[True]]
    calculatedColumn: Bool[Literal[True]]
    def __init__(
        self,
        sqref=None,
        evalError: _ConvertibleToBool | None = False,
        twoDigitTextYear: _ConvertibleToBool | None = False,
        numberStoredAsText: _ConvertibleToBool | None = False,
        formula: _ConvertibleToBool | None = False,
        formulaRange: _ConvertibleToBool | None = False,
        unlockedFormula: _ConvertibleToBool | None = False,
        emptyCellReference: _ConvertibleToBool | None = False,
        listDataValidation: _ConvertibleToBool | None = False,
        calculatedColumn: _ConvertibleToBool | None = False,
    ) -> None: ...

class IgnoredErrors(Serialisable):
    tagname: ClassVar[str]
    ignoredError: Incomplete
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, ignoredError=(), extLst: ExtensionList | None = None) -> None: ...
