from collections.abc import Callable, Iterable, MutableMapping, Sequence
from types import ModuleType
from typing import Any, Generic, TypeVar, overload
from typing_extensions import Self

_T = TypeVar("_T")

WRAPPER_ASSIGNMENTS: tuple[str, ...]
WRAPPER_UPDATES: tuple[str, ...]

def update_wrapper(
    wrapper: _T,
    wrapped: object,
    assigned: Sequence[str] = ("__module__", "__name__", "__qualname__", "__doc__", "__annotations__"),
    updated: Sequence[str] = ("__dict__",),
) -> _T: ...
def copy_globals(
    source: ModuleType,
    globs: MutableMapping[str, Any],
    only_names: Iterable[str] | None = None,
    ignore_missing_names: bool = False,
    names_to_ignore: Sequence[str] = (),
    dunder_names_to_keep: Sequence[str] = ("__implements__", "__all__", "__imports__"),
    cleanup_globs: bool = True,
) -> list[str]: ...
def import_c_accel(globs: MutableMapping[str, Any], cname: str) -> None: ...

class Lazy(Generic[_T]):
    data: _T
    def __init__(self, func: Callable[[Any], _T]) -> None: ...
    @overload
    def __get__(self, inst: None, class_: type[object]) -> Self: ...
    @overload
    def __get__(self, inst: object, class_: type[object]) -> _T: ...

class readproperty(Generic[_T]):
    func: Callable[[Any], _T]
    def __init__(
        self: readproperty[_T], func: Callable[[Any], _T]  # pyright: ignore[reportInvalidTypeVarUse]  #11780
    ) -> None: ...
    @overload
    def __get__(self, inst: None, class_: type[object]) -> Self: ...
    @overload
    def __get__(self, inst: object, class_: type[object]) -> _T: ...

class LazyOnClass(Generic[_T]):
    @classmethod
    def lazy(cls, cls_dict: MutableMapping[str, Any], func: Callable[[Any], _T]) -> None: ...
    name: str
    func: Callable[[Any], _T]
    def __init__(self, func: Callable[[Any], _T], name: str | None = None) -> None: ...
    @overload
    def __get__(self, inst: None, class_: type[object]) -> Self: ...
    @overload
    def __get__(self, inst: object, class_: type[object]) -> _T: ...

def gmctime() -> str: ...
def prereleaser_middle(data: MutableMapping[str, Any]) -> None: ...
def postreleaser_before(data: MutableMapping[str, Any]) -> None: ...
