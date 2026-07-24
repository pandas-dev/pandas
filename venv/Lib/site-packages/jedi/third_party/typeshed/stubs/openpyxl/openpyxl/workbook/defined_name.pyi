from _typeshed import ConvertibleToInt, Incomplete
from collections import defaultdict
from collections.abc import Generator, Iterator
from re import Pattern
from typing import ClassVar, Final, Literal

from openpyxl.descriptors import Sequence
from openpyxl.descriptors.base import Alias, Bool, Integer, String, _ConvertibleToBool
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.formula.tokenizer import _TokenOperandSubtypes, _TokenTypesNotOperand

RESERVED: Final[frozenset[str]]
RESERVED_REGEX: Final[Pattern[str]]

class DefinedName(Serialisable):
    tagname: ClassVar[str]
    name: String[Literal[False]]
    comment: String[Literal[True]]
    customMenu: String[Literal[True]]
    description: String[Literal[True]]
    help: String[Literal[True]]
    statusBar: String[Literal[True]]
    localSheetId: Integer[Literal[True]]
    hidden: Bool[Literal[True]]
    function: Bool[Literal[True]]
    vbProcedure: Bool[Literal[True]]
    xlm: Bool[Literal[True]]
    functionGroupId: Integer[Literal[True]]
    shortcutKey: String[Literal[True]]
    publishToServer: Bool[Literal[True]]
    workbookParameter: Bool[Literal[True]]
    attr_text: Incomplete
    value: Alias
    def __init__(
        self,
        name: str,
        comment: str | None = None,
        customMenu: str | None = None,
        description: str | None = None,
        help: str | None = None,
        statusBar: str | None = None,
        localSheetId: ConvertibleToInt | None = None,
        hidden: _ConvertibleToBool | None = None,
        function: _ConvertibleToBool | None = None,
        vbProcedure: _ConvertibleToBool | None = None,
        xlm: _ConvertibleToBool | None = None,
        functionGroupId: ConvertibleToInt | None = None,
        shortcutKey: str | None = None,
        publishToServer: _ConvertibleToBool | None = None,
        workbookParameter: _ConvertibleToBool | None = None,
        attr_text=None,
    ) -> None: ...
    @property
    def type(self) -> _TokenTypesNotOperand | _TokenOperandSubtypes: ...
    @property
    def destinations(self) -> Generator[tuple[str, str]]: ...
    @property
    def is_reserved(self) -> str | None: ...
    @property
    def is_external(self) -> bool: ...
    def __iter__(self) -> Iterator[tuple[str, str]]: ...

class DefinedNameDict(dict[str, DefinedName]):
    def add(self, value: DefinedName) -> None: ...

class DefinedNameList(Serialisable):
    tagname: ClassVar[str]
    definedName: Sequence[list[DefinedName]]
    def __init__(self, definedName: list[DefinedName] | tuple[DefinedName, ...] = ()) -> None: ...
    def by_sheet(self) -> defaultdict[int, DefinedNameDict]: ...
    def __len__(self) -> int: ...
