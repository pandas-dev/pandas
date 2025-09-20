# NOTE(scipy-stubs): This private module should not be used outside of scipy-stubs

from collections.abc import Sequence
from types import TracebackType
from typing import SupportsIndex, TypeAlias, final, type_check_only

__all__ = "AnyShape", "ExitMixin"

# helper mixins
@type_check_only
class ExitMixin:
    @final
    def __exit__(self, /, type: type[BaseException] | None, value: BaseException | None, tb: TracebackType | None) -> None: ...

# equivalent to `numpy._typing._shape._ShapeLike`
AnyShape: TypeAlias = SupportsIndex | Sequence[SupportsIndex]
