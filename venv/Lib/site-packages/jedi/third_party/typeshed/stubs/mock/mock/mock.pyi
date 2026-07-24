from _typeshed import Incomplete
from collections.abc import Callable, Coroutine, Iterable, Mapping, Sequence
from contextlib import AbstractContextManager
from types import TracebackType
from typing import Any, ClassVar, Generic, Literal, TypeVar, overload, type_check_only
from typing_extensions import ParamSpec, Self

_F = TypeVar("_F", bound=Callable[..., Any])
_AF = TypeVar("_AF", bound=Callable[..., Coroutine[Any, Any, Any]])
_T = TypeVar("_T")
_TT = TypeVar("_TT", bound=type[Any])
_R = TypeVar("_R")
_P = ParamSpec("_P")

__all__ = (
    "Mock",
    "MagicMock",
    "patch",
    "sentinel",
    "DEFAULT",
    "ANY",
    "call",
    "create_autospec",
    "AsyncMock",
    "ThreadingMock",
    "FILTER_DIR",
    "NonCallableMock",
    "NonCallableMagicMock",
    "mock_open",
    "PropertyMock",
    "seal",
)

class InvalidSpecError(Exception): ...

FILTER_DIR: bool

class _SentinelObject:
    def __init__(self, name: str) -> None: ...
    name: str

class _Sentinel:
    def __getattr__(self, name: str) -> _SentinelObject: ...

sentinel: _Sentinel
DEFAULT: _SentinelObject

class _Call(tuple[Any, ...]):
    def __new__(
        cls, value: Any = (), name: Incomplete | None = "", parent=None, two: bool = False, from_kall: bool = True
    ) -> Self: ...
    name: Any
    parent: Any
    from_kall: Any
    def __init__(self, value: Any = (), name=None, parent=None, two: bool = False, from_kall: bool = True) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object, /) -> bool: ...
    def __call__(self, *args: Any, **kwargs: Any) -> _Call: ...
    def __getattr__(self, attr: str) -> Any: ...
    @property
    def args(self) -> tuple[Any, ...]: ...
    @property
    def kwargs(self) -> dict[str, Any]: ...
    def call_list(self) -> _CallList: ...

call: _Call

class _CallList(list[_Call]):
    def __contains__(self, value: Any) -> bool: ...

class Base:
    def __init__(self, *args: Any, **kwargs: Any) -> None: ...

# We subclass with "Any" because mocks are explicitly designed to stand in for other types,
# something that can't be expressed with our static type system.
class NonCallableMock(Base, Any):
    def __new__(
        cls,
        spec: list[str] | object | type[object] | None = None,
        wraps: Any | None = None,
        name: str | None = None,
        spec_set: list[str] | object | type[object] | None = None,
        parent: NonCallableMock | None = None,
        _spec_state=None,
        _new_name: str = "",
        _new_parent: NonCallableMock | None = None,
        _spec_as_instance: bool = False,
        _eat_self: bool | None = None,
        unsafe: bool = False,
        **kwargs: Any,
    ) -> Self: ...
    def __init__(
        self,
        spec: list[str] | object | type[object] | None = None,
        wraps: Any | None = None,
        name: str | None = None,
        spec_set: list[str] | object | type[object] | None = None,
        parent: NonCallableMock | None = None,
        _spec_state=None,
        _new_name: str = "",
        _new_parent: NonCallableMock | None = None,
        _spec_as_instance: bool = False,
        _eat_self: bool | None = None,
        unsafe: bool = False,
        **kwargs: Any,
    ) -> None: ...
    def __getattr__(self, name: str) -> Any: ...
    def _calls_repr(self) -> str: ...
    def assert_called_with(_mock_self, *args: Any, **kwargs: Any) -> None: ...
    def assert_not_called(_mock_self) -> None: ...
    def assert_called_once_with(_mock_self, *args: Any, **kwargs: Any) -> None: ...
    def _format_mock_failure_message(self, args: Any, kwargs: Any, action: str = "call") -> str: ...
    def assert_called(_mock_self) -> None: ...
    def assert_called_once(_mock_self) -> None: ...
    def reset_mock(self, visited: Any = None, *, return_value: bool = False, side_effect: bool = False) -> None: ...
    def _extract_mock_name(self) -> str: ...
    def assert_any_call(self, *args: Any, **kwargs: Any) -> None: ...
    def assert_has_calls(self, calls: Sequence[_Call], any_order: bool = False) -> None: ...
    def mock_add_spec(self, spec: Any, spec_set: bool = False) -> None: ...
    def _mock_add_spec(self, spec: Any, spec_set: bool, _spec_as_instance: bool = False, _eat_self: bool = False) -> None: ...
    def attach_mock(self, mock: NonCallableMock, attribute: str) -> None: ...
    def configure_mock(self, **kwargs: Any) -> None: ...
    return_value: Any
    side_effect: Any
    called: bool
    call_count: int
    call_args: Any
    call_args_list: _CallList
    mock_calls: _CallList
    def _format_mock_call_signature(self, args: Any, kwargs: Any) -> str: ...
    def _call_matcher(self, _call: tuple[_Call, ...]) -> _Call: ...
    def _get_child_mock(self, **kw: Any) -> NonCallableMock: ...

class CallableMixin(Base):
    side_effect: Any
    def __init__(
        self,
        spec=None,
        side_effect=None,
        return_value: Any = ...,
        wraps=None,
        name=None,
        spec_set=None,
        parent=None,
        _spec_state=None,
        _new_name: Any = "",
        _new_parent=None,
        **kwargs: Any,
    ) -> None: ...
    def __call__(_mock_self, *args: Any, **kwargs: Any) -> Any: ...

class Mock(CallableMixin, NonCallableMock): ...

class _patch(Generic[_T]):
    attribute_name: Any
    getter: Callable[[], Any]
    attribute: str
    new: _T
    new_callable: Any
    spec: Any
    create: bool
    has_local: Any
    spec_set: Any
    autospec: Any
    kwargs: Mapping[str, Any]
    additional_patchers: Any
    def __init__(
        self: _patch[_T],  # pyright: ignore[reportInvalidTypeVarUse]  #11780
        getter: Callable[[], Any],
        attribute: str,
        new: _T,
        spec: Incomplete | None,
        create: bool,
        spec_set: Incomplete | None,
        autospec: Incomplete | None,
        new_callable: Incomplete | None,
        kwargs: Mapping[str, Any],
        *,
        unsafe: bool = False,
    ) -> None: ...
    def copy(self) -> _patch[_T]: ...
    def __call__(self, func: Callable[_P, _R]) -> Callable[_P, _R]: ...
    def decorate_class(self, klass: _TT) -> _TT: ...
    def decorate_callable(self, func: _F) -> _F: ...
    def decorate_async_callable(self, func: _AF) -> _AF: ...
    def decoration_helper(
        self, patched: Any, args: tuple[Any, ...], keywargs: dict[str, Any]
    ) -> AbstractContextManager[tuple[tuple[Any, ...], dict[str, Any]]]: ...
    def get_original(self) -> tuple[Any, bool]: ...
    target: Any
    temp_original: Any
    is_local: bool
    def __enter__(self) -> _T: ...
    def __exit__(
        self, exc_type: type[BaseException] | None, exc_value: BaseException | None, traceback: TracebackType | None, /
    ) -> None: ...
    def start(self) -> _T: ...
    def stop(self) -> None: ...

class _patch_dict:
    in_dict: Any
    values: Any
    clear: Any
    def __init__(self, in_dict: Any, values: Any = (), clear: Any = False, **kwargs: Any) -> None: ...
    def __call__(self, f: Any) -> Any: ...
    def decorate_callable(self, f: _F) -> _F: ...
    def decorate_async_callable(self, f: _AF) -> _AF: ...
    def decorate_class(self, klass: Any) -> Any: ...
    def __enter__(self) -> Any: ...
    def __exit__(self, *args: object) -> Any: ...
    start: Any
    stop: Any

@type_check_only
class _patcher:
    TEST_PREFIX: str
    dict: type[_patch_dict]
    @overload
    def __call__(
        self,
        target: Any,
        *,
        spec: Incomplete | None = ...,
        create: bool = ...,
        spec_set: Incomplete | None = ...,
        autospec: Incomplete | None = ...,
        new_callable: Incomplete | None = ...,
        unsafe: bool = ...,
        **kwargs: Any,
    ) -> _patch[MagicMock | AsyncMock]: ...
    # This overload also covers the case, where new==DEFAULT. In this case, the return type is _patch[Any].
    # Ideally we'd be able to add an overload for it so that the return type is _patch[MagicMock],
    # but that's impossible with the current type system.
    @overload
    def __call__(
        self,
        target: Any,
        new: _T,
        spec: Incomplete | None = ...,
        create: bool = ...,
        spec_set: Incomplete | None = ...,
        autospec: Incomplete | None = ...,
        new_callable: Incomplete | None = ...,
        *,
        unsafe: bool = ...,
        **kwargs: Any,
    ) -> _patch[_T]: ...
    @overload
    def object(
        self,
        target: Any,
        attribute: str,
        *,
        spec: Incomplete | None = ...,
        create: bool = ...,
        spec_set: Incomplete | None = ...,
        autospec: Incomplete | None = ...,
        new_callable: Incomplete | None = ...,
        unsafe: bool = ...,
        **kwargs: Any,
    ) -> _patch[MagicMock | AsyncMock]: ...
    @overload
    def object(
        self,
        target: Any,
        attribute: str,
        new: _T,
        spec: Incomplete | None = ...,
        create: bool = ...,
        spec_set: Incomplete | None = ...,
        autospec: Incomplete | None = ...,
        new_callable: Incomplete | None = ...,
        *,
        unsafe: bool = ...,
        **kwargs: Any,
    ) -> _patch[_T]: ...
    def multiple(
        self,
        target: Any,
        spec: Incomplete | None = ...,
        create: bool = ...,
        spec_set: Incomplete | None = ...,
        autospec: Incomplete | None = ...,
        new_callable: Incomplete | None = ...,
        *,
        unsafe: bool = ...,
        **kwargs: _T,
    ) -> _patch[_T]: ...
    def stopall(self) -> None: ...

patch: _patcher

class MagicMixin:
    def __init__(self, *args: Any, **kw: Any) -> None: ...

class NonCallableMagicMock(MagicMixin, NonCallableMock):
    def mock_add_spec(self, spec: Any, spec_set: bool = False) -> None: ...

class MagicMock(MagicMixin, Mock):
    def mock_add_spec(self, spec: Any, spec_set: bool = False) -> None: ...

class AsyncMockMixin(Base):
    def __init__(self, *args: Any, **kwargs: Any) -> None: ...
    def assert_awaited(_mock_self) -> None: ...
    def assert_awaited_once(_mock_self) -> None: ...
    def assert_awaited_with(_mock_self, *args: Any, **kwargs: Any) -> None: ...
    def assert_awaited_once_with(_mock_self, *args: Any, **kwargs: Any) -> None: ...
    def assert_any_await(_mock_self, *args: Any, **kwargs: Any) -> None: ...
    def assert_has_awaits(_mock_self, calls: Iterable[_Call], any_order: bool = False) -> None: ...
    def assert_not_awaited(_mock_self) -> None: ...
    def reset_mock(self, *args: Any, **kwargs: Any) -> None: ...
    await_count: int
    await_args: _Call | None
    await_args_list: _CallList
    __name__: str
    __defaults__: tuple[Any, ...]
    __kwdefaults__: dict[str, Any]
    __annotations__: dict[str, Any] | None  # type: ignore[assignment]

class AsyncMagicMixin(MagicMixin): ...

class AsyncMock(AsyncMockMixin, AsyncMagicMixin, Mock):
    # Improving the `reset_mock` signature.
    # It is defined on `AsyncMockMixin` with `*args, **kwargs`, which is not ideal.
    # But, `NonCallableMock` super-class has the better version.
    def reset_mock(self, visited: Any = None, *, return_value: bool = False, side_effect: bool = False) -> None: ...

class MagicProxy(Base):
    name: str
    parent: Any
    def __init__(self, name: str, parent: Any) -> None: ...
    def create_mock(self) -> Any: ...
    def __get__(self, obj: Any, _type=None) -> Any: ...

class _ANY(Any):
    def __eq__(self, other: object) -> Literal[True]: ...
    def __ne__(self, other: object) -> Literal[False]: ...

ANY: _ANY

def create_autospec(
    spec: Any, spec_set: Any = False, instance: Any = False, _parent=None, _name=None, *, unsafe: bool = False, **kwargs: Any
) -> Any: ...

class _SpecState:
    spec: Any
    ids: Any
    spec_set: Any
    parent: Any
    instance: Any
    name: Any
    def __init__(self, spec: Any, spec_set: Any = False, parent=None, name=None, ids=None, instance: Any = False) -> None: ...

def mock_open(mock=None, read_data: Any = "") -> Any: ...

class PropertyMock(Mock):
    def __get__(self, obj: _T, obj_type: type[_T] | None = None) -> Self: ...
    def __set__(self, obj: Any, value: Any) -> None: ...

def seal(mock: Any) -> None: ...

class ThreadingMixin(Base):
    DEFAULT_TIMEOUT: ClassVar[float | None]

    def __init__(self, *args: Any, timeout: float | None = ..., **kwargs: Any) -> None: ...
    def reset_mock(self, *args: Any, **kwargs: Any) -> None: ...
    def wait_until_called(self, *, timeout: float | None = ...) -> None: ...
    def wait_until_any_call_with(self, *args: Any, **kwargs: Any) -> None: ...

class ThreadingMock(ThreadingMixin, MagicMixin, Mock):
    # Improving the `reset_mock` signature.
    # It is defined on `ThreadingMixin` with `*args, **kwargs`, which is not ideal.
    # But, `NonCallableMock` super-class has the better version.
    def reset_mock(self, visited: Any = None, *, return_value: bool = False, side_effect: bool = False) -> None: ...
