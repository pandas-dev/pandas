import enum
import types
from typing import Any, Final, Literal

type Array = Any  # evil; do not use
type ArrayLike = Any  # evil; do not use

SCIPY_ARRAY_API: Final[str | Literal[False]] = ...
SCIPY_DEVICE: Final[str] = ...

class _ArrayClsInfo(enum.Enum):
    skip = 0
    numpy = 1
    array_like = 2
    unknown = 3

def _validate_array_cls(cls: type, sparse_ok: bool = False) -> _ArrayClsInfo: ...
def array_namespace(*arrays: Any, sparse_ok: bool = False) -> types.ModuleType: ...
