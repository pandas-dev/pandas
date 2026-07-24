from _typeshed import ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal, overload

from openpyxl.descriptors.base import Bool, Convertible, Integer, String, _ConvertibleToBool, _ConvertibleToMultiCellRange
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.worksheet.cell_range import MultiCellRange

class InputCells(Serialisable):
    tagname: ClassVar[str]
    r: String[Literal[False]]
    deleted: Bool[Literal[True]]
    undone: Bool[Literal[True]]
    val: String[Literal[False]]
    numFmtId: Integer[Literal[True]]

    @overload
    def __init__(
        self,
        r: str,
        deleted: _ConvertibleToBool | None = False,
        undone: _ConvertibleToBool | None = False,
        *,
        val: str,
        numFmtId: ConvertibleToInt | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        r: str,
        deleted: _ConvertibleToBool | None,
        undone: _ConvertibleToBool | None,
        val: str,
        numFmtId: ConvertibleToInt | None = None,
    ) -> None: ...

class Scenario(Serialisable):
    tagname: ClassVar[str]
    inputCells: Incomplete
    name: String[Literal[False]]
    locked: Bool[Literal[True]]
    hidden: Bool[Literal[True]]
    user: String[Literal[True]]
    comment: String[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    __attrs__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        inputCells=(),
        *,
        name: str,
        locked: _ConvertibleToBool | None = False,
        hidden: _ConvertibleToBool | None = False,
        count: Unused = None,
        user: str | None = None,
        comment: str | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        inputCells,
        name: str,
        locked: _ConvertibleToBool | None = False,
        hidden: _ConvertibleToBool | None = False,
        count: Unused = None,
        user: str | None = None,
        comment: str | None = None,
    ) -> None: ...
    @property
    def count(self) -> int: ...

class ScenarioList(Serialisable):
    tagname: ClassVar[str]
    scenario: Incomplete
    current: Integer[Literal[True]]
    show: Integer[Literal[True]]
    sqref: Convertible[MultiCellRange, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        scenario=(),
        current: ConvertibleToInt | None = None,
        show: ConvertibleToInt | None = None,
        sqref: _ConvertibleToMultiCellRange | None = None,
    ) -> None: ...
    def append(self, scenario) -> None: ...
    def __bool__(self) -> bool: ...
