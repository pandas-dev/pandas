from collections.abc import Callable, Iterable
from contextlib import _GeneratorContextManager
from types import NotImplementedType
from typing import Any, ClassVar, Final, Generic, Literal, Protocol, TypeAlias, TypedDict, final, overload, type_check_only
from typing_extensions import ParamSpec, TypeVar, Unpack

import optype.numpy as onp

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

_V = TypeVar("_V")
_T = TypeVar("_T", default=object)
_T2 = TypeVar("_T2", default=object)
_S = TypeVar("_S")
_C = TypeVar("_C")
_T_co = TypeVar("_T_co", covariant=True, default=Any)
_Tss = ParamSpec("_Tss", default=...)

_DispatchType: TypeAlias = type[_T] | str
_DispatchTypeSlice: TypeAlias = (
    tuple[()]
    | tuple[_T]
    | tuple[_DispatchType[_T]]
    | tuple[_DispatchType[_T], _T]
    | tuple[_T, _DispatchType[_T]]
)  # fmt: skip

@type_check_only
class _Backend(Protocol[_T_co]):
    __ua_domain__: ClassVar[str] = ...
    @staticmethod
    def __ua_function__(method: Callable[..., Any], args: tuple[Any, ...], kwargs: dict[str, Any], /) -> _T_co: ...

@final
@type_check_only
class _DetermineBackendMultiKwargs(TypedDict, Generic[_T], total=False):
    dispatch_type: type[_T] | str

###

ArgumentExtractorType: TypeAlias = Callable[..., tuple[Dispatchable, ...]]
ArgumentReplacerType: TypeAlias = Callable[
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

    def __init__(self, /, value: _T_co, dispatch_type: _DispatchType[_T_co], coercible: onp.ToBool = True) -> None: ...
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
def create_multimethod(
    *args: ArgumentReplacerType | str | Callable[_Tss, _T], **kwargs: ArgumentReplacerType | str | Callable[_Tss, _T]
) -> Callable[[ArgumentExtractorType], _Function[_Tss, _T]]: ...
@overload
def generate_multimethod(
    argument_extractor: ArgumentExtractorType, argument_replacer: ArgumentReplacerType, domain: str, default: None = None
) -> _Function: ...
@overload
def generate_multimethod(
    argument_extractor: ArgumentExtractorType, argument_replacer: ArgumentReplacerType, domain: str, default: Callable[_Tss, _T]
) -> _Function[_Tss, _T]: ...

#
def set_backend(backend: _Backend, coerce: onp.ToBool = False, only: onp.ToBool = False) -> _SetBackendContext: ...
def skip_backend(backend: _Backend) -> _SkipBackendContext: ...

#
def set_global_backend(
    backend: _Backend, coerce: onp.ToBool = False, only: onp.ToBool = False, *, try_last: onp.ToBool = False
) -> None: ...
def register_backend(backend: _Backend) -> None: ...
def clear_backends(domain: str | None, registered: onp.ToBool = True, globals: onp.ToBool = False) -> None: ...

#
def mark_as(dispatch_type: type[_T] | str) -> Callable[[_T], Dispatchable[_T]]: ...
def all_of_type(
    arg_type: type[_T] | str,
) -> Callable[[Callable[_Tss, Iterable[_T | Dispatchable[_T2]]]], Callable[_Tss, tuple[Dispatchable[_T | _T2], ...]]]: ...

#
@overload
def wrap_single_convertor(
    convert_single: Callable[[_V, _DispatchType[_V], bool], _C],
) -> Callable[[Iterable[Dispatchable[_V]], onp.ToBool], list[_C]]: ...
@overload
def wrap_single_convertor(
    convert_single: Callable[[_V, _DispatchType[_V], bool], NotImplementedType],
) -> Callable[[Iterable[Dispatchable[_V]], onp.ToBool], NotImplementedType]: ...
@overload
def wrap_single_convertor(
    convert_single: Callable[[_V, _DispatchType[_V], bool], _C | NotImplementedType],
) -> Callable[[Iterable[Dispatchable[_V]], onp.ToBool], list[_C] | NotImplementedType]: ...
@overload
def wrap_single_convertor_instance(
    convert_single: Callable[[_S, _V, _DispatchType[_V], bool], _C],
) -> Callable[[_S, Iterable[Dispatchable[_V]], onp.ToBool], list[_C]]: ...
@overload
def wrap_single_convertor_instance(
    convert_single: Callable[[_S, _V, _DispatchType[_V], bool], NotImplementedType],
) -> Callable[[_S, Iterable[Dispatchable[_V]], onp.ToBool], NotImplementedType]: ...
@overload
def wrap_single_convertor_instance(
    convert_single: Callable[[_S, _V, _DispatchType[_V], bool], _C | NotImplementedType],
) -> Callable[[_S, Iterable[Dispatchable[_V]], onp.ToBool], list[_C] | NotImplementedType]: ...

#
def determine_backend(
    value: _V, dispatch_type: _DispatchType[_V], *, domain: str, only: onp.ToBool = True, coerce: onp.ToBool = False
) -> _SetBackendContext: ...
def determine_backend_multi(
    dispatchables: Iterable[_V | Dispatchable[_V]],
    *,
    domain: str,
    only: onp.ToBool = True,
    coerce: onp.ToBool = False,
    **kwargs: Unpack[_DetermineBackendMultiKwargs[_T]],
) -> _SetBackendContext: ...
