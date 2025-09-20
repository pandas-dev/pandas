from os import PathLike
from typing import Final, Generic
from typing_extensions import TypeVar

__all__ = ["readsav"]

_VT = TypeVar("_VT", default=object)

###

DTYPE_DICT: Final[dict[int, str]] = ...
RECTYPE_DICT: Final[dict[int, str]] = ...
STRUCT_DICT: Final[dict[str, dict[str, object]]] = ...

class Pointer:
    index: int
    def __init__(self, /, index: int) -> None: ...

class ObjectPointer(Pointer): ...

class AttrDict(dict[str, _VT], Generic[_VT]):
    def __init__(self, /, init: dict[str, object] | None = None) -> None: ...
    def __call__(self, /, name: str) -> object: ...

def readsav(
    file_name: str | bytes | PathLike[str] | PathLike[bytes],
    idict: dict[str, object] | None = None,
    python_dict: bool = False,
    uncompressed_file_name: str | None = None,
    verbose: bool = False,
) -> AttrDict[object] | dict[str, object]: ...
