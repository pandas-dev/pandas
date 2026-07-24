import types
from collections.abc import Callable, Collection
from typing import Any, NamedTuple

class _Package(NamedTuple):
    name: str
    version: str

def get_package_info(module: types.ModuleType) -> _Package: ...

class EnhancedModule(types.ModuleType):
    def __bool__(self) -> bool: ...
    def __getattribute__(self, attr: str) -> Any: ...

def passthrough_module(
    parent: types.ModuleType,
    child: str | types.ModuleType,
    allowed_attributes: Collection[str] = ...,
    *,
    callback: Callable[[str], object] = ...,
) -> types.ModuleType: ...
