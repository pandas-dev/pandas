from _typeshed import SupportsItems
from collections.abc import Callable
from typing import Any, Generic, Literal, overload
from typing_extensions import Self, TypeVar

_T = TypeVar("_T")
_DefaultT = TypeVar("_DefaultT", default=None)
_NoneLiteral = TypeVar("_NoneLiteral", default=None)
_ScopeT = TypeVar("_ScopeT", Literal["request"], Literal["response"], None, default=None)
_ScopeT2 = TypeVar("_ScopeT2", Literal["request"], Literal["response"], None)

class UpdateDict(dict[str, Any]):
    updated: Callable[..., Any] | None
    updated_args: tuple[Any, ...] | None

class exists_property(Generic[_ScopeT]):
    @overload
    def __init__(self: exists_property[None], prop: str) -> None: ...
    @overload
    def __init__(self, prop: str, type: _ScopeT) -> None: ...
    @overload
    def __get__(self, obj: None, type: type[CacheControl[Any]] | None = None) -> Self: ...
    @overload
    def __get__(self: exists_property[None], obj: CacheControl[Any], type: type[CacheControl[Any]] | None = None) -> bool: ...
    @overload
    def __get__(self, obj: CacheControl[_ScopeT], type: type[CacheControl[Any]] | None = None) -> bool: ...
    @overload
    def __set__(self: exists_property[None], obj: CacheControl[Any], value: bool | None) -> None: ...
    @overload
    def __set__(self, obj: CacheControl[_ScopeT], value: bool | None) -> None: ...
    @overload
    def __delete__(self, obj: CacheControl[Any]) -> None: ...
    @overload
    def __delete__(self, obj: CacheControl[_ScopeT]) -> None: ...

class value_property(Generic[_T, _DefaultT, _NoneLiteral, _ScopeT]):
    def __init__(self, prop: str, default: _DefaultT = None, none: _NoneLiteral = None, type: _ScopeT = None) -> None: ...  # type: ignore[assignment]
    @overload
    def __get__(self, obj: None, type: type[CacheControl[Any]] | None = None) -> Self: ...
    @overload
    def __get__(
        self: value_property[_T, _DefaultT, _NoneLiteral, None],
        obj: CacheControl[Any] | None,
        type: type[CacheControl[Any]] | None = None,
    ) -> _T | _DefaultT | _NoneLiteral: ...
    @overload
    def __get__(
        self, obj: CacheControl[_ScopeT] | None, type: type[CacheControl[Any]] | None = None
    ) -> _T | _DefaultT | _NoneLiteral: ...
    @overload
    def __set__(
        self: value_property[_T, _DefaultT, _NoneLiteral, None],
        obj: CacheControl[Any],
        value: _T | _DefaultT | Literal[True] | None,
    ) -> None: ...
    @overload
    def __set__(self, obj: CacheControl[_ScopeT], value: _T | _DefaultT | Literal[True] | None) -> None: ...
    @overload
    def __delete__(self, obj: CacheControl[Any]) -> None: ...
    @overload
    def __delete__(self, obj: CacheControl[_ScopeT]) -> None: ...

class CacheControl(Generic[_ScopeT]):
    header_value: str
    update_dict: type[UpdateDict]
    properties: dict[str, Any]
    type: _ScopeT
    def __init__(self, properties: dict[str, Any], type: _ScopeT) -> None: ...
    @overload
    @classmethod
    def parse(
        cls, header: str, updates_to: Callable[[dict[str, Any]], Any] | None = None, type: None = None
    ) -> CacheControl[None]: ...
    @overload
    @classmethod
    def parse(cls, header: str, updates_to: Callable[[dict[str, Any]], Any] | None, type: _ScopeT2) -> CacheControl[_ScopeT2]: ...
    @overload
    @classmethod
    def parse(
        cls, header: str, updates_to: Callable[[dict[str, Any]], Any] | None = None, *, type: _ScopeT2
    ) -> CacheControl[_ScopeT2]: ...
    max_stale: value_property[int, None, Literal["*"], Literal["request"]]
    min_fresh: value_property[int, None, None, Literal["request"]]
    only_if_cached: exists_property[Literal["request"]]
    public: exists_property[Literal["response"]]
    private: value_property[str, None, Literal["*"], Literal["response"]]
    no_cache: value_property[str, None, Literal["*"], None]
    no_store: exists_property[None]
    no_transform: exists_property[None]
    must_revalidate: exists_property[Literal["response"]]
    proxy_revalidate: exists_property[Literal["response"]]
    max_age: value_property[int, None, Literal[-1], None]
    s_maxage: value_property[int, None, None, Literal["response"]]
    s_max_age = s_maxage
    stale_while_revalidate: value_property[int, None, None, Literal["response"]]
    stale_if_error: value_property[int, None, None, Literal["response"]]
    def copy(self) -> Self: ...

def serialize_cache_control(properties: SupportsItems[str, Any] | CacheControl[Any]) -> str: ...
