from collections.abc import Callable, Iterable
from contextlib import _GeneratorContextManager
from types import NotImplementedType
from typing import Any, ClassVar, Final, Generic, Literal, Protocol, TypedDict, Unpack, final, overload, type_check_only
from typing_extensions import ParamSpec, TypeVar

from scipy._typing import ExitMixin

__all__ = [
    "BackendNotImplementedError",
    "Dispatchable",
    "_BackendState",
    "_Function",
    "_SetBackendContext",
    "_SkipBackendContext",
    "all_of_type",
    "clear_backends",
    "create_multimethod",
    "determine_backend",
    "determine_backend_multi",
    "generate_multimethod",
    "get_state",
    "mark_as",
    "register_backend",
    "reset_state",
    "set_backend",
    "set_global_backend",
    "set_state",
    "skip_backend",
    "wrap_single_convertor",
    "wrap_single_convertor_instance",
]

_T_co = TypeVar("_T_co", covariant=True, default=Any)
_Tss = ParamSpec("_Tss", default=...)

type _DispatchType[T] = type[T] | str
type _DispatchTypeSlice[T] = (
    tuple[()]
    | tuple[T]
    | tuple[_DispatchType[T]]
    | tuple[_DispatchType[T], T]
    | tuple[T, _DispatchType[T]]
)  # fmt: skip

@type_check_only
class _Backend(Protocol[_T_co]):
    __ua_domain__: ClassVar[str] = ...
    @staticmethod
    def __ua_function__(method: Callable[..., Any], args: tuple[Any, ...], kwargs: dict[str, Any], /) -> _T_co: ...

@final
@type_check_only
class _DetermineBackendMultiKwargs[T](TypedDict, total=False):
    dispatch_type: type[T] | str

###

type ArgumentExtractorType = Callable[..., tuple[Dispatchable, ...]]
type ArgumentReplacerType = Callable[
    [tuple[Any, ...], dict[str, Any], tuple[Dispatchable, ...]],
    tuple[tuple[Any, ...], dict[str, Any]],
]  # fmt: skip

@final
class _BackendState: ...

@final
class _SetBackendContext(ExitMixin):
    def __init__(self, /, *args: object, **kwargs: object) -> None: ...
    def __enter__(self, /) -> None: ...

@final
class _SkipBackendContext(ExitMixin):
    def __init__(self, /, *args: object, **kwargs: object) -> None: ...
    def __enter__(self, /) -> None: ...

@final
class _Function(Generic[_Tss, _T_co]):
    @property
    def arg_extractor(self, /) -> ArgumentExtractorType: ...
    @property
    def arg_replacer(self, /) -> ArgumentReplacerType: ...
    @property
    def default(self, /) -> Callable[_Tss, _T_co] | None: ...
    @property
    def domain(self, /) -> str: ...
    def __init__(self, /, *args: object, **kwargs: object) -> None: ...
    def __call__(self, /, *args: _Tss.args, **kwargs: _Tss.kwargs) -> None: ...

class Dispatchable(Generic[_T_co]):
    value: _T_co
    type: _DispatchType[_T_co]
    coercible: Final[bool]

    def __init__(self, /, value: _T_co, dispatch_type: _DispatchType[_T_co], coercible: bool = True) -> None: ...
    @overload
    def __getitem__(self, index: Literal[1, -1], /) -> _T_co: ...
    @overload
    def __getitem__(self, index: Literal[0, -2], /) -> _DispatchType[_T_co]: ...
    @overload
    def __getitem__(self, index: slice, /) -> _DispatchTypeSlice[_T_co]: ...

class BackendNotImplementedError(NotImplementedError): ...

def get_state() -> _BackendState: ...
def reset_state() -> _GeneratorContextManager[None]: ...
def set_state(state: _BackendState) -> _GeneratorContextManager[None]: ...

#
def create_multimethod[**Tss, T](
    *args: ArgumentReplacerType | str | Callable[Tss, T], **kwargs: ArgumentReplacerType | str | Callable[Tss, T]
) -> Callable[[ArgumentExtractorType], _Function[Tss, T]]: ...
@overload
def generate_multimethod(
    argument_extractor: ArgumentExtractorType, argument_replacer: ArgumentReplacerType, domain: str, default: None = None
) -> _Function: ...
@overload
def generate_multimethod[**Tss, T](
    argument_extractor: ArgumentExtractorType, argument_replacer: ArgumentReplacerType, domain: str, default: Callable[Tss, T]
) -> _Function[Tss, T]: ...

#
def set_backend(backend: _Backend, coerce: bool = False, only: bool = False) -> _SetBackendContext: ...
def skip_backend(backend: _Backend) -> _SkipBackendContext: ...

#
def set_global_backend(backend: _Backend, coerce: bool = False, only: bool = False, *, try_last: bool = False) -> None: ...
def register_backend(backend: _Backend) -> None: ...
def clear_backends(domain: str | None, registered: bool = True, globals: bool = False) -> None: ...

#
def mark_as[T](dispatch_type: type[T] | str) -> Callable[[T], Dispatchable[T]]: ...

#
@type_check_only
class _AllOfTypeCallable[T](Protocol):
    def __call__[**Tss, T2](
        self, arg: Callable[Tss, Iterable[T | Dispatchable[T2]]], /
    ) -> Callable[Tss, tuple[Dispatchable[T | T2], ...]]: ...

#
def all_of_type[T](arg_type: type[T] | str) -> _AllOfTypeCallable[T]: ...

#
@overload
def wrap_single_convertor[V](
    convert_single: Callable[[V, _DispatchType[V], bool], NotImplementedType],
) -> Callable[[Iterable[Dispatchable[V]], bool], NotImplementedType]: ...
@overload
def wrap_single_convertor[V, C](
    convert_single: Callable[[V, _DispatchType[V], bool], C],
) -> Callable[[Iterable[Dispatchable[V]], bool], list[C]]: ...

#
@overload
def wrap_single_convertor_instance[S, V](
    convert_single: Callable[[S, V, _DispatchType[V], bool], NotImplementedType],
) -> Callable[[S, Iterable[Dispatchable[V]], bool], NotImplementedType]: ...
@overload
def wrap_single_convertor_instance[S, V, C](
    convert_single: Callable[[S, V, _DispatchType[V], bool], C],
) -> Callable[[S, Iterable[Dispatchable[V]], bool], list[C]]: ...

#
def determine_backend[V](
    value: V, dispatch_type: _DispatchType[V], *, domain: str, only: bool = True, coerce: bool = False
) -> _SetBackendContext: ...
def determine_backend_multi[V, T](
    dispatchables: Iterable[V | Dispatchable[V]],
    *,
    domain: str,
    only: bool = True,
    coerce: bool = False,
    **kwargs: Unpack[_DetermineBackendMultiKwargs[T]],
) -> _SetBackendContext: ...
