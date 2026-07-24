from collections.abc import Callable, Iterable, Mapping
from typing import Final, Literal, Protocol, TypedDict, TypeVar, type_check_only
from typing_extensions import Unpack

from _win32typing import PyHKEY, PyIID

_T = TypeVar("_T", PyHKEY, int)

@type_check_only
class _RegisterClass(Protocol):
    _reg_clsid_: PyIID

@type_check_only
class _RegisterFlag(TypedDict, total=False):
    quiet: bool
    debug: bool
    finalize_register: Callable[[], None]

@type_check_only
class _UnregisterFlag(TypedDict, total=False):
    quiet: bool
    finalize_unregister: Callable[[], None]

@type_check_only
class _ElevatedFlag(TypedDict, total=False):
    quiet: bool
    unattended: bool
    hwnd: int

@type_check_only
class _CommandFlag(_RegisterFlag, _UnregisterFlag, _ElevatedFlag):  # type: ignore[misc]
    ...

CATID_PythonCOMServer: Final = "{B3EF80D0-68E2-11D0-A689-00C04FD658FF}"

def recurse_delete_key(path: str | None, base: PyHKEY | int = -2147483648) -> None: ...
def RegisterServer(
    clsid: PyIID,
    pythonInstString: str | None = None,
    desc: str | None = None,
    progID: str | None = None,
    verProgID: str | None = None,
    defIcon: str | None = None,
    threadingModel: Literal["apartment", "both", "free", "neutral"] = "both",
    policy: str | None = None,
    catids: list[PyIID] = [],
    other: Mapping[str, str] = {},
    addPyComCat: bool | None = None,
    dispatcher: str | None = None,
    clsctx: int | None = None,
    addnPath: str | None = None,
) -> None: ...
def GetUnregisterServerKeys(
    clsid: PyIID, progID: str | None = None, verProgID: str | None = None, customKeys: Iterable[tuple[str, _T]] | None = None
) -> list[tuple[str, _T | int]]: ...
def UnregisterServer(
    clsid: PyIID,
    progID: str | None = None,
    verProgID: str | None = None,
    customKeys: Iterable[tuple[str, PyHKEY | int]] | None = None,
) -> None: ...
def GetRegisteredServerOption(clsid: PyIID, optionName: str) -> str | None: ...
def RegisterClasses(*classes: type[_RegisterClass], **flags: Unpack[_RegisterFlag]) -> None: ...
def UnregisterClasses(*classes: type[_RegisterClass], **flags: Unpack[_UnregisterFlag]) -> None: ...
def UnregisterInfoClasses(*classes: type[_RegisterClass]) -> list[tuple[str, PyHKEY | int]]: ...
def ReExecuteElevated(flags: _ElevatedFlag) -> None: ...
def UseCommandLine(*classes: type[_RegisterClass], **flags: Unpack[_CommandFlag]) -> list[tuple[str, PyHKEY | int]] | None: ...
def RegisterPyComCategory() -> None: ...
