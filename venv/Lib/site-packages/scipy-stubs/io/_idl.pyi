from _typeshed import StrOrBytesPath
from typing import Any, Final, Generic, Literal, overload
from typing_extensions import TypeVar

__all__ = ["readsav"]

_VT = TypeVar("_VT", default=Any)

###

DTYPE_DICT: Final[dict[int, str]] = ...  # undocumented
RECTYPE_DICT: Final[dict[int, str]] = ...  # undocumented
STRUCT_DICT: Final[dict[str, dict[str, Any]]] = ...  # undocumented

class Pointer:  # undocumented
    index: int
    def __init__(self, /, index: int) -> None: ...

class ObjectPointer(Pointer): ...  # undocumented

class AttrDict(dict[str, _VT], Generic[_VT]):  # undocumented
    def __init__(self, /, init: dict[str, Any] | None = None) -> None: ...
    def __call__(self, /, name: str) -> Any: ...

@overload
def readsav(
    file_name: StrOrBytesPath,
    idict: dict[str, Any] | None = None,
    python_dict: Literal[False] = False,
    uncompressed_file_name: str | None = None,
    verbose: bool = False,
) -> AttrDict[Any]: ...
@overload
def readsav(
    file_name: StrOrBytesPath,
    idict: dict[str, Any] | None = None,
    *,
    python_dict: Literal[True],
    uncompressed_file_name: str | None = None,
    verbose: bool = False,
) -> dict[str, Any]: ...
@overload
def readsav(
    file_name: StrOrBytesPath,
    idict: dict[str, Any] | None,
    python_dict: Literal[True],
    uncompressed_file_name: str | None = None,
    verbose: bool = False,
) -> dict[str, Any]: ...
