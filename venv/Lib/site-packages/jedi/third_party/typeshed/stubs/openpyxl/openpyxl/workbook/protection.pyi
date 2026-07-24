from _typeshed import ConvertibleToInt, Incomplete
from typing import ClassVar, Literal, overload
from typing_extensions import Self

from openpyxl.descriptors.base import Alias, Bool, Integer, String, _ConvertibleToBool
from openpyxl.descriptors.serialisable import Serialisable

from ..xml._functions_overloads import _SupportsIterAndAttribAndTextAndGet

class WorkbookProtection(Serialisable):
    tagname: ClassVar[str]
    workbook_password: Alias
    workbookPasswordCharacterSet: String[Literal[True]]
    revision_password: Alias
    revisionsPasswordCharacterSet: String[Literal[True]]
    lockStructure: Bool[Literal[True]]
    lock_structure: Alias
    lockWindows: Bool[Literal[True]]
    lock_windows: Alias
    lockRevision: Bool[Literal[True]]
    lock_revision: Alias
    revisionsAlgorithmName: String[Literal[True]]
    revisionsHashValue: Incomplete
    revisionsSaltValue: Incomplete
    revisionsSpinCount: Integer[Literal[True]]
    workbookAlgorithmName: String[Literal[True]]
    workbookHashValue: Incomplete
    workbookSaltValue: Incomplete
    workbookSpinCount: Integer[Literal[True]]
    __attrs__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        workbookPassword=None,
        workbookPasswordCharacterSet: str | None = None,
        revisionsPassword=None,
        revisionsPasswordCharacterSet: str | None = None,
        lockStructure: _ConvertibleToBool | None = None,
        lockWindows: _ConvertibleToBool | None = None,
        lockRevision: _ConvertibleToBool | None = None,
        revisionsAlgorithmName: str | None = None,
        revisionsHashValue=None,
        revisionsSaltValue=None,
        revisionsSpinCount: ConvertibleToInt | None = None,
        workbookAlgorithmName: str | None = None,
        workbookHashValue=None,
        workbookSaltValue=None,
        workbookSpinCount: ConvertibleToInt | None = None,
    ) -> None: ...
    @overload
    def set_workbook_password(self, value: str = "", already_hashed: Literal[False] = False) -> None: ...
    @overload
    def set_workbook_password(self, value: str | None, already_hashed: Literal[True]) -> None: ...
    @overload
    def set_workbook_password(self, value: str | None = "", *, already_hashed: Literal[True]) -> None: ...
    @property
    def workbookPassword(self) -> str | None: ...
    @workbookPassword.setter
    def workbookPassword(self, value: str) -> None: ...
    @overload
    def set_revisions_password(self, value: str = "", already_hashed: Literal[False] = False) -> None: ...
    @overload
    def set_revisions_password(self, value: str | None, already_hashed: Literal[True]) -> None: ...
    @overload
    def set_revisions_password(self, value: str | None = "", *, already_hashed: Literal[True]) -> None: ...
    @property
    def revisionsPassword(self) -> str | None: ...
    @revisionsPassword.setter
    def revisionsPassword(self, value: str) -> None: ...
    @classmethod
    def from_tree(cls, node: _SupportsIterAndAttribAndTextAndGet) -> Self: ...

DocumentSecurity = WorkbookProtection

class FileSharing(Serialisable):
    tagname: ClassVar[str]
    readOnlyRecommended: Bool[Literal[True]]
    userName: String[Literal[True]]
    reservationPassword: Incomplete
    algorithmName: String[Literal[True]]
    hashValue: Incomplete
    saltValue: Incomplete
    spinCount: Integer[Literal[True]]
    def __init__(
        self,
        readOnlyRecommended: _ConvertibleToBool | None = None,
        userName: str | None = None,
        reservationPassword=None,
        algorithmName: str | None = None,
        hashValue=None,
        saltValue=None,
        spinCount: ConvertibleToInt | None = None,
    ) -> None: ...
