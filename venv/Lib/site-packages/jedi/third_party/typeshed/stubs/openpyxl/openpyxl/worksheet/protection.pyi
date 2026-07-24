from _typeshed import ConvertibleToInt, Incomplete
from typing import ClassVar, Literal, overload

from openpyxl.descriptors.base import Alias, Bool, Integer, String, _ConvertibleToBool
from openpyxl.descriptors.serialisable import Serialisable

class _Protected:
    @overload
    def set_password(self, value: str = "", already_hashed: Literal[False] = False) -> None: ...
    @overload
    def set_password(self, value: str | None, already_hashed: Literal[True]) -> None: ...
    @overload
    def set_password(self, value: str | None = "", *, already_hashed: Literal[True]) -> None: ...
    @property
    def password(self) -> str | None: ...
    @password.setter
    def password(self, value: str) -> None: ...

class SheetProtection(Serialisable, _Protected):
    tagname: ClassVar[str]
    sheet: Bool[Literal[False]]
    enabled: Alias
    objects: Bool[Literal[False]]
    scenarios: Bool[Literal[False]]
    formatCells: Bool[Literal[False]]
    formatColumns: Bool[Literal[False]]
    formatRows: Bool[Literal[False]]
    insertColumns: Bool[Literal[False]]
    insertRows: Bool[Literal[False]]
    insertHyperlinks: Bool[Literal[False]]
    deleteColumns: Bool[Literal[False]]
    deleteRows: Bool[Literal[False]]
    selectLockedCells: Bool[Literal[False]]
    selectUnlockedCells: Bool[Literal[False]]
    sort: Bool[Literal[False]]
    autoFilter: Bool[Literal[False]]
    pivotTables: Bool[Literal[False]]
    saltValue: Incomplete
    spinCount: Integer[Literal[True]]
    algorithmName: String[Literal[True]]
    hashValue: Incomplete
    __attrs__: ClassVar[tuple[str, ...]]
    password: Incomplete
    def __init__(
        self,
        sheet: _ConvertibleToBool = False,
        objects: _ConvertibleToBool = False,
        scenarios: _ConvertibleToBool = False,
        formatCells: _ConvertibleToBool = True,
        formatRows: _ConvertibleToBool = True,
        formatColumns: _ConvertibleToBool = True,
        insertColumns: _ConvertibleToBool = True,
        insertRows: _ConvertibleToBool = True,
        insertHyperlinks: _ConvertibleToBool = True,
        deleteColumns: _ConvertibleToBool = True,
        deleteRows: _ConvertibleToBool = True,
        selectLockedCells: _ConvertibleToBool = False,
        selectUnlockedCells: _ConvertibleToBool = False,
        sort: _ConvertibleToBool = True,
        autoFilter: _ConvertibleToBool = True,
        pivotTables: _ConvertibleToBool = True,
        password=None,
        algorithmName: str | None = None,
        saltValue=None,
        spinCount: ConvertibleToInt | None = None,
        hashValue=None,
    ) -> None: ...
    @overload
    def set_password(self, value: str = "", already_hashed: Literal[False] = False) -> None: ...
    @overload
    def set_password(self, value: str | None, already_hashed: Literal[True]) -> None: ...
    @overload
    def set_password(self, value: str | None = "", *, already_hashed: Literal[True]) -> None: ...
    def enable(self) -> None: ...
    def disable(self) -> None: ...
    def __bool__(self) -> bool: ...
