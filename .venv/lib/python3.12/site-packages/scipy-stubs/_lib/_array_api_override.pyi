import enum
import sys
import types
from typing import Any, Final, Literal, TypeAlias

# mypy<=1.19.0 workaround, see https://github.com/python/mypy/pull/20392
if sys.version_info >= (3, 14):
    __conditional_annotations__: Final[set[int]] = ...

Array: TypeAlias = Any
ArrayLike: TypeAlias = Any

SCIPY_ARRAY_API: Final[str | Literal[False]] = ...
SCIPY_DEVICE: Final[str] = ...

class _ArrayClsInfo(enum.Enum):
    skip = 0
    numpy = 1
    array_like = 2
    unknown = 3

def _validate_array_cls(cls: type) -> _ArrayClsInfo: ...
def array_namespace(*arrays: Array) -> types.ModuleType: ...
