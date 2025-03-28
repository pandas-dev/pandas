from _typeshed import Incomplete
from typing import IO, Any, Literal, TypeVar

_T = TypeVar("_T")

__all__ = ["Require", "find_module", "get_module_constant", "extract_constant"]

def find_module(
    module, paths=None
) -> tuple[IO[Any], str | None, tuple[str, Literal["", "r", "rb"], Literal[7, 6, 1, 2, 3, -1]]]: ...

class Require:
    def __init__(
        self,
        name,
        requested_version,
        module,
        homepage: str = "",
        attribute: Incomplete | None = None,
        format: Incomplete | None = None,
    ) -> None: ...
    def full_name(self): ...
    def version_ok(self, version): ...
    def get_version(
        self, paths: Incomplete | None = None, default: _T | Literal["unknown"] = "unknown"
    ) -> _T | Literal["unknown"] | None | Any: ...
    def is_present(self, paths: Incomplete | None = None): ...
    def is_current(self, paths: Incomplete | None = None): ...

def get_module_constant(module, symbol, default: _T | int = -1, paths: Incomplete | None = None) -> _T | int | None | Any: ...
def extract_constant(code, symbol, default: _T | int = -1) -> _T | int | None | Any: ...
